#!/usr/bin/env python3
"""MCM v3.3.1 (KPP) -> ASP translation, hybrid with the baseline lumped mechanism.

Every MCM rate is converted to a structural ASP rate flag, so GasChemicalMechanism.in is independent of temperature and
pressure (no rate is evaluated at a fixed T or P).  Run `python3 mcm2asp.py` for the self-check (numerical comparison of
every translated rate with the MCM Fortran expression over a T, P grid).
"""
import math, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcm_kpp import read_eqn, D
from mcm_rates import Term, parse, load_constants, fall_params

FN2, FO2 = 0.79, 0.21                # air composition used wherever MCM writes N2 / O2 explicitly (same as ASP's hard-wired flags 8, 9)
RENAME_MCM = {'SA': 'H2SO4', 'HCL': 'HCl'}      # MCM names that are hard-wired ASP names
# MCM J index -> (PhotoRates.in column of the baseline table, MCM J of the same process that column represents)
# (the rate is  J_asp = TUV column x  J_MCM(target)/J_MCM(analog)  evaluated at zenith 0; ratio 1 where the process is the same)
MCM_J = {   # l, m, n of the MCM parameterisation  l cos(chi)^m exp(-n sec(chi))
 'J_O3_O1D': (6.073E-05, 1.743, 0.474), 'J_O3_O3P': (4.775E-04, 0.298, 0.08), 'J_H2O2': (1.041E-05, 0.723, 0.279),
 'J_NO2': (1.165E-02, 0.244, 0.267), 'J_NO3_NO': (2.485E-02, 0.168, 0.108), 'J_NO3_NO2': (1.747E-01, 0.155, 0.125),
 'J_HONO': (2.644E-03, 0.261, 0.288), 'J_HNO3': (9.312E-07, 1.23, 0.307), 'J_HCHO_H': (4.642E-05, 0.762, 0.353),
 'J_HCHO_H2': (6.853E-05, 0.477, 0.323), 'J_CH3CHO': (7.344E-06, 1.202, 0.417), 'J_C2H5CHO': (2.879E-05, 1.067, 0.358),
 'J_C3H7CHO_HCO': (2.792E-05, 0.805, 0.338), 'J_C3H7CHO_C2H4': (1.675E-05, 0.805, 0.338),
 'J_IPRCHO': (7.914E-05, 0.764, 0.364), 'J_MACR_HCO': (1.482E-06, 0.396, 0.298), 'J_MACR_H': (1.482E-06, 0.396, 0.298),
 'J_C5HPALD1': (7.600E-04, 0.396, 0.298), 'J_CH3COCH3': (7.992E-07, 1.578, 0.271), 'J_MEK': (5.804E-06, 1.092, 0.377),
 'J_MVK_CO': (2.4246E-06, 0.395, 0.296), 'J_MVK_C2H3': (2.424E-06, 0.395, 0.296), 'J_GLYOX_H2': (6.845E-05, 0.13, 0.201),
 'J_GLYOX_HCHO': (1.032E-05, 0.13, 0.201), 'J_GLYOX_HCO': (3.802E-05, 0.644, 0.312), 'J_MGLYOX': (1.537E-04, 0.17, 0.208),
 'J_BIACET': (3.326E-04, 0.148, 0.215), 'J_CH3OOH': (7.649E-06, 0.682, 0.279), 'J_CH3NO3': (1.588E-06, 1.154, 0.318),
 'J_C2H5NO3': (1.907E-06, 1.244, 0.335), 'J_NC3H7NO3': (2.485E-06, 1.196, 0.328), 'J_IC3H7NO3': (4.095E-06, 1.111, 0.316),
 'J_TC4H9NO3': (1.135E-05, 0.974, 0.309), 'J_NOA': (4.365E-05, 1.089, 0.323)}
def j0(name):
    l, m, n = MCM_J[name]; return l * math.exp(-n)           # value at zenith angle 0
# target -> (column, [MCM J names whose sum the column represents])
J_MAP = {
 'J_O3_O1D': (2, ['J_O3_O1D']), 'J_O3_O3P': (3, ['J_O3_O3P']), 'J_H2O2': (5, ['J_H2O2']), 'J_NO2': (6, ['J_NO2']),
 'J_NO3_NO': (7, ['J_NO3_NO']), 'J_NO3_NO2': (8, ['J_NO3_NO2']), 'J_HONO': (12, ['J_HONO']), 'J_HNO3': (13, ['J_HNO3']),
 'J_HCHO_H': (17, ['J_HCHO_H']), 'J_HCHO_H2': (18, ['J_HCHO_H2']), 'J_CH3CHO': (19, ['J_CH3CHO']),
 'J_C2H5CHO': (22, ['J_C2H5CHO']), 'J_C3H7CHO_HCO': (22, ['J_C2H5CHO']), 'J_C3H7CHO_C2H4': (22, ['J_C2H5CHO']),
 'J_IPRCHO': (22, ['J_C2H5CHO']),
 'J_MACR_HCO': (25, ['J_MACR_HCO', 'J_MACR_H']), 'J_MACR_H': (25, ['J_MACR_HCO', 'J_MACR_H']),
 'J_C5HPALD1': (25, ['J_MACR_HCO', 'J_MACR_H']),
 'J_CH3COCH3': (26, ['J_CH3COCH3']), 'J_MEK': (28, ['J_MEK']),
 'J_MVK_CO': (27, ['J_MVK_CO', 'J_MVK_C2H3']), 'J_MVK_C2H3': (27, ['J_MVK_CO', 'J_MVK_C2H3']),
 'J_GLYOX_H2': (45, ['J_GLYOX_H2', 'J_GLYOX_HCHO']), 'J_GLYOX_HCHO': (45, ['J_GLYOX_H2', 'J_GLYOX_HCHO']),
 'J_GLYOX_HCO': (44, ['J_GLYOX_HCO']), 'J_MGLYOX': (46, ['J_MGLYOX']), 'J_BIACET': (47, ['J_BIACET']),
 'J_CH3OOH': (31, ['J_CH3OOH']), 'J_CH3NO3': (36, ['J_IC3H7NO3']), 'J_C2H5NO3': (36, ['J_IC3H7NO3']),
 'J_NC3H7NO3': (36, ['J_IC3H7NO3']), 'J_IC3H7NO3': (36, ['J_IC3H7NO3']), 'J_TC4H9NO3': (36, ['J_IC3H7NO3']),
 'J_NOA': (38, ['J_NOA'])}
def j_col_scale(name):
    col, ref = J_MAP[name]
    return col, j0(name) / sum(j0(r) for r in ref)

# baseline compound names -> MCM names (individual compounds the MCM treats explicitly)
BASE2MCM = {
 'HNO4': 'HO2NO2', 'MEOH': 'CH3OH', 'Formic_Acid': 'HCOOH', 'ETHE': 'C2H4', 'Ethane': 'C2H6', 'Acetaldehyde': 'CH3CHO',
 'ETOH': 'C2H5OH', 'Acetic_Acid': 'CH3CO2H', 'GLY': 'GLYOX', 'Acetylene': 'C2H2', 'Propene': 'C3H6', 'Propane': 'C3H8',
 'Propanal': 'C2H5CHO', 'Acetone': 'CH3COCH3', 'MGLY': 'MGLYOX', 'Butadiene': 'C4H6', '1-Butene': 'BUT1ENE',
 'cis-2-Butene': 'CBUT2ENE', 'trans-2-Butene': 'TBUT2ENE', 'n-Butane': 'NC4H10', 'i-Butane': 'IC4H10',
 'Butanal': 'C3H7CHO', 'MethylPropanal': 'IPRCHO', 'MCR': 'MACR', 'ISOP': 'C5H8', 'API': 'APINENE', 'LIM': 'LIMONENE',
 'BEN': 'BENZENE', 'BALD': 'BENZAL', 'PHEN': 'PHENOL', 'TOL': 'TOLUENE', 'XYM': 'MXYL', 'XYP': 'PXYL', 'XYO': 'OXYL',
 'HCl': 'HCl'}
# SOA-precursor counters: MCM precursor -> counter (product added to every reaction with OH, O3 or NO3, as in the baseline)
COUNTER = {'C5H8': 'ISOSOA', 'TOLUENE': 'TOLSOA', 'MXYL': 'XYLSOA', 'PXYL': 'XYLSOA', 'OXYL': 'XYLSOA',
           'APINENE': 'TERPSOA', 'LIMONENE': 'TERPSOA', 'BPINENE': 'TERPSOA'}
OXIDANTS = {'OH', 'O3', 'NO3'}
POOL = {'RO2_T', 'RCO3_T'}                    # ASP peroxy pools (not species of either mechanism)
NOT_COVERED_OK = {'M', 'H2O', 'hv', 'PROD'}

def mcm_species():
    sp = collections.OrderedDict()
    for l in open(D + 'mcm_3-3-1_species_complete.tsv'):
        if l.startswith('*') or l.startswith('Name') or not l.strip(): continue
        f = l.rstrip('\n').split('\t')
        sp[f[0]] = dict(smiles=f[1], formula=f[4], mass=float(f[5]) if f[5] else None, perox=(f[7] == 'true'))
    return sp

def ro2_list():
    return set(x.strip() for x in open('/home/claude/mcm_src/AtChem2/mcm/peroxy-radicals_v3.3.1') if x.strip())

def fmt(x):
    if abs(x) < 1e-300: return '0'
    return '%.8g' % x

def line(re_, pr, flag, params):
    def side(s, is_p):
        out = []
        for c, n in s:
            out.append(n if abs(c - 1) < 1e-9 else '%s %s' % (fmt(c), n))
        return ' + '.join(out)
    head = '%s => %s' % (side(re_, False), side(pr, True))
    ps = '; '.join(fmt(p) for p in params)
    fl = ('%s;' % flag) if flag is not None else ''
    return head + ' ; ' + fl + (' ' if fl and ps else '') + (ps + ';' if ps else '')

class MCMTranslator:
    def __init__(self):
        self.names, self.fall, self.assign = load_constants()
        self.allnames = dict(self.names)
        for k in self.fall: self.allnames[k] = [Term(sym=k)]
        for k in ('KMT05', 'KMT06', 'KMT11', 'KMT18'): self.allnames[k] = [Term(sym=k)]
        self.species = mcm_species()
        self.ro2 = ro2_list()
        self.notes = []
        self.skipped = []

    # ---- rate -> list of (extra_reactants, extra_products, flag, params, kind) given reactant count ---------------------
    def rate_items(self, r, nr):
        k = r['k']
        # forms with (1 - c exp(...)) that cannot be sums of positive terms
        m = re.fullmatch(r'3\.8E-13\*EXP\(780\./TEMP\)\*\(1\.-1\./\(1\.\+498\.\*EXP\(-1160\./TEMP\)\)\)', k)
        if m: return [([], [], 7, [3.8e-13, 780., 498., -1160.])]
        m = re.fullmatch(r'3\.8E-13\*EXP\(780\./TEMP\)\*\(1\./\(1\.\+498\.\*EXP\(-1160\./TEMP\)\)\)', k)
        if m: return [([], [], 8, [3.8e-13, 780., 498., -1160.])]
        if k == '2.*KCH3O2*RO2*0.5*(1.-7.18*EXP(-885./TEMP))':
            A = self.names['KCH3O2'][0]
            return [(['RO2_T'], ['RO2_T'], 9, [A.c, A.C, 7.18, -885.])]
        if k == '5.00E-12*O2*3.2*(1.-EXP(-550./TEMP))':
            return [(['O2'], [], 9, [1.6e-11, 0., 1., -550.])]
        if r['n'] == 20:       # HO2 + HO2 (+ H2O):  (2.2e-13 e^{600/T} + 1.9e-33 M e^{980/T}) (1 + 1.4e-21 [H2O] e^{2200/T})
            return [([], [], 4, [2.2e-13, 600., 1.9e-33, 980.]),
                    ([], [], 10, [2.2e-13 * 1.4e-21, 600. + 2200., 1.9e-33 * 1.4e-21, 980. + 2200.])]
        if r['n'] == 1:        # O + O2 + M: ASP hard-wired first-order flag 9 (same constants as the MCM expression)
            return [([], [], 9, [])]
        if r['n'] == 6:        # O1D + M  (CACM / JPL constants, identical to the MCM expression)
            return [([], [], 8, [])]
        terms = parse(k, self.allnames)
        key = collections.OrderedDict()
        for t in terms:                       # merge identical forms
            kk = (t.B, t.C, t.e3, t.o2, t.n2, t.h2o, t.m, t.ro2, t.J, t.sym)
            key[kk] = key.get(kk, 0.0) + t.c
        items = []
        for (B, C, e3, o2, n2, h2o, m_, ro2, J, sym), c in key.items():
            assert c > 0, (r['n'], k)
            extra = []
            if o2 > 1 or n2 or h2o > 1 or ro2 > 1: raise ValueError('unsupported powers %s' % k)
            if o2: extra.append('O2')
            if h2o: extra.append('H2O')
            if J:
                assert nr == 1 and not extra and not (m_ or ro2 or sym or e3), (r['n'], k)
                col, sc = j_col_scale(J)
                items.append(([], [], 0, [col] + ([] if abs(sc * c - 1) < 1e-9 else [sc * c]))); continue
            if sym:
                items.extend(self.fall_item(sym, c, nr, r)); continue
            if ro2:
                assert nr == 1 and not extra and not (m_ or e3), (r['n'], k)
                items.append(([], [], 12, [c, B, C, 1.])); continue
            if e3:
                assert nr == 1 and not extra and not m_, (r['n'], k)
                items.append(([], [], 14, [c, B, C, e3])); continue
            nb = nr + len(extra)
            if m_:
                assert m_ == 1 and nb == 2, (r['n'], k)
                items.append((extra + ['M'], ['M'], 'B3_2', [c, B, C, 1.])); continue
            if nb == 3:       # termolecular with O2 as third body (NO + NO + O2): O2 = FO2 * M
                assert nr == 2 and extra == ['O2'], (r['n'], k)
                items.append((['M'], ['M'], 'B3_2', [c * FO2, B, C, 1.])); continue
            assert nb in (1, 2), (r['n'], k)
            items.append((extra, [], 'AR1' if nb == 1 else 'AR2', [c, B, C, 1.]))
        return items

    def fall_item(self, sym, c, nr, r):
        if sym == 'KMT05':           # 1.44e-13 (1 + M/4.2e19)
            assert nr == 2; return [([], [], 4, [c * 1.44e-13, 0., c * 1.44e-13 / 4.2e19, 0.])]
        if sym == 'KMT11':
            assert nr == 2; return [([], [], 5, [c * 2.4e-14, 460., c * 2.7e-17, 2199., c * 6.5e-34, 1335.])]
        if sym == 'KMT18':           # 9.5e-39 [O2] e^{5270/T} / (1 + 7.5e-29 [O2] e^{5610/T}),  [O2] = FO2 M
            assert nr == 2; return [([], [], 5, [0., 0., c * 9.5e-39 / 7.5e-29, 5270. - 5610., c * 9.5e-39 * FO2, 5270.])]
        A0, n0, E0, Ai, ni, Ei, F, N = fall_params(self.fall[sym])
        Td = self.fall[sym]['Tdep']
        if Td:
            assert nr == 2
            return [(['M'], ['M'], 'B3_5', [c * A0, n0, E0, c * Ai, ni, Ei, Td[0], Td[1], Td[2]])]
        if nr == 1: return [([], [], 13, [c * A0, n0, E0, c * Ai, ni, Ei, F, N, 1., 0.])]
        assert nr == 2, (sym, r['n'])
        return [(['M'], ['M'], 'B3_4', [c * A0, n0, E0, c * Ai, ni, Ei, F, N])]

    # ---- one MCM reaction -> ASP lines ----------------------------------------------------------------------------------
    def translate(self, r):
        ren = lambda n: RENAME_MCM.get(n, n)
        reac = [n for c, n in r['re'] if n != 'hv' for _ in [0]]
        reac = [ren(n) for n in reac]
        if r['n'] == 15631:          # C522CO3H + OH + hv = C522CO3 : J(J_CH3OOH)  (J used as a bimolecular k: upstream error)
            self.skipped.append((r['n'], 'C522CO3H + OH + hv = C522CO3 : J(J_CH3OOH)',
                                 'photolysis rate constant used as a bimolecular rate constant in the MCM v3.3.1 file as published; omitted'))
            return []
        prods = collections.OrderedDict()
        for c, n in r['pr']:
            if n in ('PROD', 'hv'): continue
            n = ren(n); prods[n] = prods.get(n, 0.0) + c
        nr = len(reac)
        assert all(c == 1 for c, _ in r['re']) and nr in (1, 2), r
        cnt = COUNTER_FOR(reac)
        out = []
        for extra_re, extra_pr, flag, params in self.rate_items(r, nr):
            re_ = [(1.0, n) for n in reac + list(extra_re)]
            pr = [(c, n) for n, c in prods.items()] + [(1.0, n) for n in extra_pr]
            if not [1 for _, n in pr if n != 'M'] and not prods: pr = [(1.0, 'M')] if not extra_pr else pr
            if cnt: pr.append((1.0, cnt))
            if flag == 'AR1': out.append(line(re_, pr, 2, params))
            elif flag == 'AR2': out.append(line(re_, pr, None, params))
            elif flag == 'B3_2': out.append(line(re_, pr, 2, params))
            elif flag == 'B3_4': out.append(line(re_, pr, 4, params))
            elif flag == 'B3_5': out.append(line(re_, pr, 5, params))
            else: out.append(line(re_, pr, flag, params))
        return out

def COUNTER_FOR(reac):
    for n in reac:
        if n in COUNTER and any(o in OXIDANTS for o in reac): return COUNTER[n]
    return None
