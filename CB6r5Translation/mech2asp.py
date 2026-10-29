"""Generic translator: CMAQ mechanism definition (mech_*.def)  ->  ASP GasChemicalMechanism.in

Used for CRACMM2, CRACMM3, SAPRC07TC (AE6) and SAPRC07TIC (AE7i).  The carbon-bond translator
(cmaq2asp.py) keeps working unchanged; this module re-uses its rate helpers and ASP-rate port.

All rate forms are recognised *structurally* (not by reaction label):

  CMAQ form                                         ASP representation
  ------------------------------------------------  ---------------------------------------------------------
  A (T/300)^n exp(-E/T)  unimolecular               1-body flag 2
  A (T/300)^n exp(-E/T)  bimolecular                2-body (flag-less)
  photolysis  s/<LABEL>                             1-body, flag 0 + PhotoRates column
  %3, 2 terms  (k1 + k2 M)                          2-body flag 4
  %3, 2 terms with a H2O reactant                   2-body flag 10  (C_H2O (k1 + M k2))
  %3, 3 terms  (R051, NO + HO2 -> HNO3)             2-body flag 11
  %2, 3 terms  (HNO3 + OH)                          2-body flag 5
  fall-off, bimolecular (2 reactants + M)           3-body flag 3 (IUPAC N), flag-less (JPL, N=1) or flag 4 (any N, with E)
  fall-off, unimolecular (decomposition)            1-body flag 13  (Keq = 1)
  'A@E *E<ref>'  (reverse via equilibrium constant) 1-body flag 13  (with Keq)
  A*K<ref>                                          expanded statically
  O3P + O2 + M                                      1-body flag 9 (ASP hard-wired) - approximation, reported
  O1D + M  /  O1D + O2 and O1D + N2                 1-body flag 8 (ASP hard-wired) - approximation, reported
  NO + NO + O2                                      3-body flag 2 with [O2] = 0.2095 [M]
Flags 13 (1-body) and 4 (3-body) need common_asp >= 40c1760.
"""
import math
import re
from collections import OrderedDict

import cmaq2asp as c2a
from cmaq2asp import Mair, fmt, ncd, side_str, KB, O2_FRAC, C_H2O_TEST, TEST_POINTS
from cmaq_def import read_classified

MECHS = '/home/claude/usepa/cmaq/CCTM/src/MECHS/'
CMAQ_ENV = {'M', 'O2', 'N2'}


# ---------------------------------------------------------------------------------------------
# per-mechanism configuration
# ---------------------------------------------------------------------------------------------
CRACMM_CORE = r'^(R|TRP|ROCARO|RAM|RBAL|T|HC)\d+[a-z]*$'
SAPRC_CORE = r'^(\d+|BR\d+[a-z]*|R\d+[a-z]*|BP\d+[a-z]*|P\d+[a-z]*|PX\d+|PZ\d+|BE\d+[a-z]*|BT\d+[a-z]*|IS\d+[a-z]*|IA\d+[a-z]*|IC\d+[a-z]*|BL\d+[a-z]*)$'

# SOA first-step counters of the baseline (RACM2) mechanism: (parent, oxidant) -> counter.
# The ASP baseline appends the counter to the first oxidation step of these parents.
CRACMM_COUNTER_PAIRS = OrderedDict()
for _p, _c, _ox in [('HC5', 'HC5SOA', ('OH',)), ('HC10', 'HC8SOA', ('OH',)),
                    ('OLT', 'OLTSOA', ('OH', 'O3', 'NO3')), ('OLI', 'OLISOA', ('OH', 'O3', 'NO3')),
                    ('BDE13', 'OLISOA', ('OH', 'O3', 'NO3')), ('TOL', 'TOLSOA', ('OH',)),
                    ('XYL', 'XYLSOA', ('OH',)), ('CSL', 'XYLSOA', ('OH',)),
                    ('API', 'TERPSOA', ('OH', 'O3', 'NO3')), ('LIM', 'TERPSOA', ('OH', 'O3', 'NO3')),
                    ('ISO', 'ISOSOA', ('OH', 'O3', 'NO3'))]:
    for _o in _ox:
        CRACMM_COUNTER_PAIRS[(_p, _o)] = _c

SAPRC_COUNTER_MAP = {'TOLRO2': 'TOLSOA', 'XYLRO2': 'XYLSOA', 'ISOPRXN': 'ISOSOA', 'TRPRXN': 'TERPSOA'}
SAPRC_COUNTER_DROP = {'BENZRO2', 'PAHRO2', 'SULRXN', 'SESQRXN', 'TOLNRXN', 'TOLHRXN', 'XYLNRXN', 'XYLHRXN',
                      'BNZNRXN', 'BNZHRXN'}

CFG = OrderedDict()
CFG['CRACMM2'] = dict(
    label='CRACMM2', suffix='CRACMM2', def_=MECHS + 'cracmm2/mech_cracmm2.def', core=CRACMM_CORE,
    csqy=MECHS + 'cracmm2/CSQY_DATA_cracmm2', gc=MECHS + 'cracmm2/GC_cracmm2.nml', family='cracmm',
    rename={'HO': 'OH', 'O3P': 'O', 'SULF': 'H2SO4', 'ECH4': 'CH4', 'HCL': 'HCl'},
    skip={'R364'},                         # ECH4 + HO: the CH4 tracer reaction (CH4 + HO is R071)
    eliminate={'DELTA_C', 'DELTA_N', 'DELTA_SI'},
    drop_prod=r'^(V[A-Z0-9]+|A[A-Z0-9]+J|SULRXN)$',      # SVOC / aerosol bookkeeping products
    drop_react=r'^(V[A-Z0-9]+|A[A-Z0-9]+J)$',                        # reactions of SVOC bookkeeping species
    counter_map={}, counter_drop=set(), counter_pairs=CRACMM_COUNTER_PAIRS,
    ivoc_counters=['HC5SOA', 'HC8SOA', 'OLTSOA', 'OLISOA', 'TOLSOA', 'XYLSOA', 'ISOSOA', 'TERPSOA'],
    acetic='ORA2', o1d='split', ho2h2o='R046')
CFG['CRACMM3'] = dict(CFG['CRACMM2'], label='CRACMM3', suffix='CRACMM3', def_=MECHS + 'cracmm3/mech_cracmm3.def',
                      csqy=MECHS + 'cracmm3/CSQY_DATA_cracmm3', gc=MECHS + 'cracmm3/GC_cracmm3.nml')
CFG['SAPRC07TC'] = dict(
    label='SAPRC07TC', suffix='SAPRC07TC', def_=MECHS + 'saprc07tc_ae6_aq/mech_saprc07tc_ae6_aq.def', core=SAPRC_CORE,
    csqy=MECHS + 'saprc07tc_ae6_aq/CSQY_DATA_saprc07tc_ae6_aq',
    sp_table=MECHS + 'mechanism_information/saprc07tc_ae6_aq/saprc07tc_ae6_aq_species_table.md', family='saprc',
    rename={'O3P': 'O', 'HO2H': 'H2O2', 'SULF': 'H2SO4', 'HCL': 'HCl'},
    skip=set(), eliminate={'XC', 'XN'}, drop_prod=r'^$', drop_react=r'^$',
    counter_map=SAPRC_COUNTER_MAP, counter_drop=SAPRC_COUNTER_DROP,
    counter_pairs={('TERP', 'NO3'): 'TERPSOA'},
    ivoc_counters=['TOLSOA', 'XYLSOA', 'ISOSOA', 'TERPSOA'], acetic='CCOOH', o1d='M', ho2h2o='38')
CFG['SAPRC07TIC'] = dict(
    CFG['SAPRC07TC'], label='SAPRC07TIC', suffix='SAPRC07TIC',
    def_=MECHS + 'saprc07tic_ae7i_aq/mech_saprc07tic_ae7i_aq.def',
    csqy=MECHS + 'saprc07tic_ae7i_aq/CSQY_DATA_saprc07tic_ae7i_aq',
    sp_table=MECHS + 'mechanism_information/saprc07tic_ae7i_aq/saprc07tic_ae7i_aq_species_table.md')


# ---------------------------------------------------------------------------------------------
# CMAQ rate evaluation (extended with the forms found in these mechanisms)
# ---------------------------------------------------------------------------------------------
def _arr(t, T):
    A, n, E = t
    return A * (T / 300.0) ** n * math.exp(-E / T)


def cmaq_k(r, T, M, by):
    """CMAQ rate coefficient of reaction r (no pseudo-reactant concentration factors)"""
    if r.kind == 'arr':
        return _arr(r.terms[0], T)
    if r.kind == 'ref':
        return r.scale * cmaq_k(by[r.target], T, M, by)
    if r.kind == 'falloff':
        k0 = _arr(r.terms[0], T) * M
        ki = _arr(r.terms[1], T)
        F, N = r.terms[2][0], r.terms[3][0]
        return (k0 / (1.0 + k0 / ki)) * F ** (1.0 / (1.0 + (math.log10(k0 / ki) / N) ** 2))
    if r.kind == 'eref':
        return cmaq_k(by[r.target], T, M, by) / _arr(r.terms[0], T)
    if r.kind == 'pct':
        t = r.terms
        if r.ptype == '%3' and len(t) == 2:
            return _arr(t[0], T) + _arr(t[1], T) * M
        if r.ptype == '%3' and len(t) == 3:      # CMAQ FALLOFF_T11 (R051 of CRACMM)
            return _arr(t[0], T) + _arr(t[1], T) * M + _arr(t[2], T)
        if r.ptype == '%2':
            k0, k2, k3 = (_arr(x, T) for x in t)
            return k0 + k3 * M / (1.0 + k3 * M / k2)
    raise ValueError('cannot evaluate %s (%s %s)' % (r.label, r.kind, r.ptype))


# ---------------------------------------------------------------------------------------------
# port of the new ASP rate branches (common_asp 40c1760) + 2-body flag 11, on top of c2a.asp_rate
# ---------------------------------------------------------------------------------------------
def asp_rate(nb, params, T, M, C_H2O, P_mbar):
    RR = {i + 3: v for i, v in enumerate(params)}
    g = lambda i: RR.get(i, 0.0)
    f = int(g(3))
    if nb == 1 and f == 13:
        k0 = g(4) * (T / 300.) ** (-g(5)) * math.exp(g(6) / T) * M
        ki = g(7) * (T / 300.) ** (-g(8)) * math.exp(g(9) / T)
        kf = (k0 * ki) / (k0 + ki) * g(10) ** (1. / (1. + (math.log10(k0 / ki) / g(11)) ** 2))
        return kf / (g(12) * math.exp(g(13) / T))
    if nb == 3 and f == 4:
        k0 = g(4) * (T / 300.) ** (-g(5)) * math.exp(g(6) / T) * M
        ki = g(7) * (T / 300.) ** (-g(8)) * math.exp(g(9) / T)
        return (k0 * ki) / (k0 + ki) * g(10) ** (1. / (1. + (math.log10(k0 / ki) / g(11)) ** 2))
    if nb == 2 and f == 11:
        k1 = g(4) * math.exp(g(5) / T)
        k2 = g(6) / T + g(7) * P_mbar * 100.0 + g(8)
        return k1 * k2 / 100.0
    return c2a.asp_rate(nb, params, T, M, C_H2O, P_mbar)


# ---------------------------------------------------------------------------------------------
# photolysis columns
# ---------------------------------------------------------------------------------------------
class PhotoColumns:
    def __init__(self):
        self.cols = OrderedDict()                 # (label, scale) -> 1-based column

    def index(self, label, scale):
        key = (label, round(scale, 12))
        if key not in self.cols:
            self.cols[key] = len(self.cols) + 1
            if len(self.cols) > 51:
                raise ValueError('more than 51 photolysis columns needed')
        return self.cols[key]


# ---------------------------------------------------------------------------------------------
class Translation:
    def __init__(self, mech, T, P_mbar, fallback_species=None):
        self.cfg = cfg = CFG[mech]
        self.mech = mech
        self.T0, self.P0 = T, P_mbar
        self.all = read_classified(cfg['def_'])
        self.by = {r.label: r for r in self.all}
        core = re.compile(cfg['core'])
        self.core = [r for r in self.all if core.match(r.label)]
        self.photo = PhotoColumns()
        self.lines = []          # dict(label, text, line, nb, params, note)
        self.notes = []          # dict(label, level, text)
        self.species = OrderedDict()
        self.typos = [(r.label, r.typo) for r in self.core if r.typo]
        self._o1d_done = False
        self.n_dropped = 0
        for r in self.core:
            self._one(r)

    # ---- helpers -------------------------------------------------------------------------
    def note(self, r, level, text):
        self.notes.append(dict(label=r.label, level=level, text=text, cmaq=r.text))

    def _emit(self, r, reacts, prods, nb, params, add_M=False):
        rs = list(reacts) + (['M'] if add_M else [])
        ps = list(prods) if prods else [(1.0, 'M')]
        if add_M and prods:
            ps = ps + [(1.0, 'M')]
        line = '%s => %s ; %s' % (' + '.join(rs), side_str(ps), '; '.join(fmt(p) for p in params) + ';')
        self.lines.append(dict(label=r.label, text=r.text, line=line, nb=nb, params=list(params),
                               reacts=[x for x in reacts]))
        for n in list(reacts) + [n for _, n in ps]:
            if n not in ('M', 'H2O'):
                self.species.setdefault(n, None)

    def _name(self, n):
        return self.cfg['rename'].get(n, n)

    def _clean_products(self, r):
        cfg = self.cfg
        d = OrderedDict()
        dropped = []
        counters = []
        for c, n in r.products:
            if n in cfg['eliminate']:
                continue
            if n in cfg['counter_drop']:
                dropped.append(n)
                continue
            if re.match(cfg['drop_prod'], n):
                dropped.append(n)
                continue
            if n in cfg['counter_map']:
                counters.append((n, cfg['counter_map'][n]))
                n = cfg['counter_map'][n]
            n = self._name(n)
            d[n] = d.get(n, 0.0) + c
        return d, dropped, counters

    # ---- one reaction ---------------------------------------------------------------------
    def _one(self, r):
        cfg = self.cfg
        if r.label in cfg['skip']:
            self.note(r, 'skipped', 'not translated (CH4 tracer duplicate of the CH4 + OH reaction)')
            return
        reacts_raw = [n for _, n in r.reactants]
        if any(re.match(cfg['drop_react'], n) for n in reacts_raw):
            self.note(r, 'dropped', 'reactant is an SVOC/IVOC bookkeeping species of the CMAQ mechanism; ASP keeps its own '
                      'IVOC block (%s)' % ', '.join(n for n in reacts_raw if re.match(cfg['drop_react'], n)))
            self.n_dropped += 1
            return
        # follow *K<ref> chains
        base, fac, chain = r, 1.0, []
        while base.kind == 'ref':
            fac *= base.scale
            chain.append(base.target)
            base = self.by[base.target]
        if base.kind in ('hetero', 'op', 'unknown'):
            self.note(r, 'unsupported', 'rate type %s not translated' % base.kind)
            return
        if chain:
            self.note(r, 'expanded', "ASP cannot reference another reaction's rate; expanded exactly as %s x %s"
                      % (fmt(fac), '->'.join(chain)))
        reacts = [self._name(n) for n in reacts_raw]
        prods, dropped, counters = self._clean_products(r)
        if dropped:
            self.note(r, 'soa', 'CMAQ SOA/SVOC bookkeeping products not carried: ' + ', '.join(dropped))
        if counters:
            self.note(r, 'soa', 'CMAQ SOA counter re-keyed onto the ASP counter: ' +
                      ', '.join('%s->%s' % x for x in counters))
        # baseline-style counters on the first oxidation step of selected parents
        body = [x for x in reacts if x not in ('M', 'O2', 'N2', 'H2O')]
        if len(body) == 2:
            for (parent, ox), ctr in cfg['counter_pairs'].items():
                if sorted(body) == sorted([parent, ox]) and ctr not in prods:
                    prods[ctr] = prods.get(ctr, 0.0) + 1.0
                    self.note(r, 'soa', 'ASP SOA counter %s attached to %s + %s (as in the baseline mechanism)' %
                              (ctr, parent, ox))
        prods_l = [(c, n) for n, c in prods.items() if abs(c) > 1e-12]
        self._from_base(r, reacts, prods_l, base, fac)

    def _from_base(self, r, reacts, prods, base, fac):
        cfg = self.cfg
        body = [x for x in reacts if x not in ('M', 'O2', 'N2')]
        # ---- photolysis ---------------------------------------------------------------------
        if base.kind == 'photo':
            assert len(body) == 1, r.label
            col = self.photo.index(base.target, base.scale * fac)
            self._emit(r, body, prods, 1, [0.0, float(col)])
            return
        if base.kind == 'arr' and base.terms[0][0] * fac == 0.0:
            self.note(r, 'zero', 'rate constant is zero in CMAQ: reaction not emitted')
            return
        # ---- O3P + O2 + M ---------------------------------------------------------------------
        if reacts.count('M') == 1 and 'O2' in reacts and 'O' in reacts and len(reacts) == 3:
            self._emit(r, ['O'], prods, 1, [9.0])
            self.note(r, 'approx', 'ASP hard-wires O + O2 + M (flag 9): 5.6e-34 N2 / 6.0e-34 O2 (T/300)^-2.6 instead of the CMAQ '
                      'coefficient; O(3P) is >99%% converted to O3 either way')
            return
        # ---- O1D quenching ----------------------------------------------------------------------
        if reacts in (['O1D', 'M'], ['O1D', 'O2'], ['O1D', 'N2']):
            if not self._o1d_done:
                self._emit(r, ['O1D'], prods, 1, [8.0])
                self.note(r, 'approx', 'ASP hard-wires O1D quenching by air (flag 8, IUPAC N2/O2 form, one reaction); CMAQ has '
                          '%s' % ('O1D + M' if reacts[1] == 'M' else 'separate O1D + O2 and O1D + N2 reactions') +
                          ' (agrees within ~1%)')
                self._o1d_done = True
            else:
                self.note(r, 'merged', 'merged into the flag-8 O1D quenching reaction emitted for the first O1D quencher')
            return
        # ---- NO + NO + O2 -------------------------------------------------------------------------
        if sorted(reacts) == ['NO', 'NO', 'O2']:
            A, n, E = base.terms[0]
            self._emit(r, ['NO', 'NO'], prods, 3, [2.0, A * fac * O2_FRAC * 300.0 ** (-n), n, -E, 1.0], add_M=True)
            self.note(r, 'info', 'termolecular NO+NO+O2: 3-body rate flag 2 with [O2] = %.4f [M] (needs the common_asp flag-2 fix)'
                      % O2_FRAC)
            return
        # ---- H2O as third reactant ---------------------------------------------------------------------
        if 'H2O' in body and len(body) == 3:
            nonw = [x for x in body if x != 'H2O']
            assert base.kind == 'pct' and base.ptype == '%3' and len(base.terms) == 2, r.label
            (A1, n1, E1), (A2, n2, E2) = base.terms
            assert n1 == 0 and n2 == 0
            self._emit(r, nonw, prods, 2, [10.0, A1 * fac, -E1, A2 * fac, -E2])
            return
        nb = len(body)
        # ---- plain Arrhenius ------------------------------------------------------------------------
        if base.kind == 'arr':
            A, n, E = base.terms[0]
            Ap = A * fac * 300.0 ** (-n)
            if nb == 1:
                self._emit(r, body, prods, 1, [2.0, Ap, n, -E, 1.0])
            elif nb == 2:
                assert Ap < 1.0, (r.label, Ap)
                self._emit(r, body, prods, 2, [Ap, n, -E, 1.0])
            else:
                raise ValueError('%s: %d reactants' % (r.label, nb))
            return
        if base.kind == 'pct':
            t = base.terms
            if base.ptype == '%3' and len(t) == 2 and nb == 2:
                (A1, n1, E1), (A2, n2, E2) = t
                assert n1 == 0 and n2 == 0, r.label
                self._emit(r, body, prods, 2, [4.0, A1 * fac, -E1, A2 * fac, -E2])
                return
            if base.ptype == '%2' and nb == 2:
                assert all(x[1] == 0 for x in t), r.label
                (A0, n0, E0), (A2, n2, E2), (A3, n3, E3) = t
                self._emit(r, body, prods, 2, [5.0, A0 * fac, -E0, A2 * fac, -E2, A3 * fac, -E3])
                self.note(r, 'info', 'CMAQ %2 form k0 + k3[M]/(1+k3[M]/k2) = ASP 2-body flag 5 (needs common_asp >= 508ba8d)')
                return
            if base.ptype == '%3' and len(t) == 3 and nb == 2:
                (A1, n1, E1), (A2, n2, E2), (A3, n3, E3) = t
                assert n1 == -1.0 and n2 == 1.0 and n3 == 0.0 and E1 == E2 == E3, r.label
                # k = A1 (300/T) e^{-E/T} + A2 (T/300) M e^{-E/T} + A3 e^{-E/T},   M = P/(kT)
                #   = e^{-E/T} * (C/T + D P_Pa + F)/100  with A = 1 (ASP 2-body flag 11)
                C = 100.0 * A1 * 300.0 * fac
                D = 100.0 * A2 * fac / 300.0 / KB * 1e-6
                F = 100.0 * A3 * fac
                self._emit(r, body, prods, 2, [11.0, 1.0, -E1, C, D, F])
                self.note(r, 'info', 'CMAQ FALLOFF_T11 (k1 (T/300)^-1 + k2 (T/300) M + k3) exp(270/T) rewritten exactly as ASP '
                          '2-body flag 11 [T and P dependence of M cancels analytically]')
                return
        if base.kind == 'falloff':
            (A0, n0, E0), (Ai, ni, Ei) = base.terms[0], base.terms[1]
            F, N = base.terms[2][0], base.terms[3][0]
            if nb == 2:
                if E0 == 0.0 and Ei == 0.0 and abs(N - ncd(F)) < 1e-4:
                    self._emit(r, body, prods, 3, [3.0, A0 * fac, -n0, Ai * fac, -ni, F], add_M=True)
                elif E0 == 0.0 and Ei == 0.0 and abs(N - 1.0) < 1e-9:
                    self._emit(r, body, prods, 3, [A0 * fac, -n0, Ai * fac, -ni, F], add_M=True)
                else:
                    self._emit(r, body, prods, 3, [4.0, A0 * fac, -n0, -E0, Ai * fac, -ni, -Ei, F, N], add_M=True)
                    self.note(r, 'info', 'fall-off with F=%g, N=%g%s: ASP 3-body flag 4 (explicit N)' %
                              (F, N, ' and exponential terms' if (E0 or Ei) else ''))
                return
            if nb == 1:
                self._emit(r, body, prods, 1, [13.0, A0 * fac, -n0, -E0, Ai * fac, -ni, -Ei, F, N, 1.0, 0.0])
                self.note(r, 'info', 'unimolecular fall-off: ASP 1-body flag 13 (explicit N, Keq = 1)')
                return
        if base.kind == 'eref':
            fw = self.by[base.target]
            while fw.kind == 'ref':
                fw = self.by[fw.target]
            assert fw.kind == 'falloff' and nb == 1, r.label
            (A0, n0, E0), (Ai, ni, Ei) = fw.terms[0], fw.terms[1]
            F, N = fw.terms[2][0], fw.terms[3][0]
            Ke, ne, Ee = base.terms[0]
            assert ne == 0.0, r.label
            self._emit(r, body, prods, 1, [13.0, A0, -n0, -E0, Ai, -ni, -Ei, F, N, Ke, -Ee])
            self.note(r, 'info', 'reverse reaction via the equilibrium constant of %s: ASP 1-body flag 13 (fall-off / Keq)' %
                      base.target)
            return
        self.note(r, 'unsupported', 'rate form %s %s (%d reactants) not handled' % (base.kind, base.ptype, nb))


# ---------------------------------------------------------------------------------------------
def compare(tr):
    """max relative deviation of every translated non-photolysis reaction from the CMAQ expression"""
    out = []
    by = tr.by
    for ln in tr.lines:
        lab, nb, params = ln['label'], ln['nb'], ln['params']
        r = by[lab]
        if nb == 1 and int(params[0]) == 0 and len(params) == 2:
            continue
        reacts = [tr._name(n) for _, n in r.reactants]
        worst, wp = 0.0, None
        for T, P in TEST_POINTS:
            M = Mair(T, P)
            if reacts in (['O1D', 'M'], ['O1D', 'O2'], ['O1D', 'N2']):
                tot = 0.0
                for q in tr.core:
                    qr = [tr._name(n) for _, n in q.reactants]
                    if qr in (['O1D', 'M'], ['O1D', 'O2'], ['O1D', 'N2']):
                        k = cmaq_k(q, T, M, by)
                        tot += k * (M if qr[1] == 'M' else {'O2': O2_FRAC, 'N2': 0.7808}[qr[1]] * M)
                ref = tot
            else:
                base, fac = r, 1.0
                while base.kind == 'ref':
                    fac *= base.scale
                    base = by[base.target]
                if base.kind == 'eref':
                    ref = cmaq_k(r, T, M, by)
                else:
                    ref = cmaq_k(base, T, M, by) * fac
                if reacts.count('M') == 1 and 'O2' in reacts and 'O' in reacts:
                    ref *= M * O2_FRAC * M
                elif sorted(reacts) == ['NO', 'NO', 'O2']:
                    ref *= O2_FRAC * M
                elif 'H2O' in reacts and len(reacts) == 3:
                    ref *= C_H2O_TEST
            val = asp_rate(nb, params, T, M, C_H2O_TEST, P)
            rel = abs(val - ref) / abs(ref)
            if rel > worst:
                worst, wp = rel, (T, P)
        out.append(dict(label=lab, maxrel=worst, at=wp, line=ln['line']))
    return out


if __name__ == '__main__':
    import sys
    mech = sys.argv[1]
    tr = Translation(mech, 288.4, 880.0)
    print('%s: %d core statements -> %d ASP lines, %d dropped (SVOC reactants)' % (mech, len(tr.core), len(tr.lines),
                                                                                   tr.n_dropped))
    bad = [c for c in compare(tr) if c['maxrel'] > 1e-6]
    for c in sorted(bad, key=lambda c: -c['maxrel']):
        print('  %-8s %.3g  %s' % (c['label'], c['maxrel'], c['line'][:100]))
    from collections import Counter
    print(Counter(n['level'] for n in tr.notes))
    for n in tr.notes:
        if n['level'] in ('unsupported',):
            print('  UNSUPPORTED', n['label'], n['text'])
    print('species', len(tr.species), 'photo columns', len(tr.photo.cols))
