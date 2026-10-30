#!/usr/bin/env python3
"""GEOS-Chem fullchem v14.8.0 (KPP) -> ASP reactions.

Every gas-phase rate is expressed through ASP's structural rate flags (nothing is evaluated at a fixed T or P); the GEOS-Chem
user functions that ASP lacked were added to common_asp as new flags:
  1-body 17 (ARRPLUS/TUNPLUS);  2-body 12/13 (GC_NIT/GC_ALK), 14 (GC_EPO), 15 (DMS+OH), 16 (GCJPLAC), 17 (GLYX+NO3),
  18 (GLYC/HAC + OH branching), 19 (ARRPLUS), 20/21 (GC_ISO1/2, ISOP + OH).
Excluded: halogen (Cl, Br, I, F) chemistry, sea-salt, cloud and aqueous-sulfur reactions, and all aerosol-phase species.
Aerosol uptake: ASP's first-order heterogeneous flag 3 (fixed gamma; ASP's own surface area).
"""
import math, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aspfmt import line
from gc_kpp import read_eqn
import gc_species

R_ = 8.314462
NUM = r'[-+]?[0-9]*\.?[0-9]+(?:[eEdD][+-]?[0-9]+)?(?:_dp)?'
def num(s): return float(s.replace('_dp', '').replace('d', 'e').replace('D', 'e'))
def isnum(s): return re.fullmatch(NUM, s.strip()) is not None

# ------------------------------------------------------------------ species classification
AEROSOL = {'NIT', 'NITs', 'SO4s', 'O3A', 'O3C', 'AONITA', 'IONITA', 'MONITA', 'INDIOL', 'LVOCOA', 'SOAGX', 'SOAIE', 'HMS', 'OCS',
           'AERI', 'SALAAL', 'SALCAL', 'ISALA', 'ISALC', 'BrSALA', 'BrSALC', 'SALACL', 'SALCCL', 'SALASO2', 'SALCSO2', 'SALASO3', 'SALCSO3',
           'MDL'}
DUMMY_P = {'PH2SO4', 'PSO4AQ', 'PSO4MP', 'PHMSAQ', 'PHMSMP', 'LHMSAQ', 'LHMSMP'}
class Species:
    def __init__(self):
        self.T = gc_species.build()
        self.excluded = set(); self.dummy = set()
        for n, v in self.T.items():
            comp = v['comp'] or {}
            if n in AEROSOL or any(e in comp for e in ('Cl', 'Br', 'I', 'F')) or n in ('Cl', 'Br', 'I'): self.excluded.add(n)
            if v['com'].startswith('Dummy') or n in DUMMY_P: self.dummy.add(n)
        for n in ('Br', 'Br2', 'Cl', 'Cl2', 'I', 'I2', 'BrO', 'ClO', 'IO', 'HBr', 'HCl', 'HI'): self.excluded.add(n)
    def mw(self, n): return self.T[n]['mw']
SP = Species()
RENAME = {'SO4': 'H2SO4', 'HNO2': 'HONO'}          # ASP's enhancement-ratio output looks species up by the names H2SO4 / HONO
ACYL = {'MCO3', 'RCO3', 'ACO3', 'GCO3', 'BZCO3', 'MEKCO3', 'MACR1OO', 'AROMCO3', 'C96O2', 'PINO3', 'LIMO3'}
PEROX_EXTRA = {'XRO2', 'TRO2', 'NRO2', 'ZRO2', 'AROMRO2', 'TLFUO2', 'KO2', 'ATO2', 'PO2', 'PIO2', 'R4O2', 'R7O2', 'ALK4O2', 'ETO2', 'A3O2',
               'B3O2', 'MO2', 'BRO2', 'OLND', 'OLNN', 'OTHRO2', 'PRN1', 'R4N1', 'ALK4N1', 'R7N1', 'INO2B', 'INO2D', 'C4HVP1', 'C4HVP2'}
NOT_PEROX = {'NO2', 'HO2', 'SO2', 'CO2', 'O2', 'H2O2', 'CH2OO', 'CH3CHOO', 'OClO', 'ClOO'}
def peroxy_flag(n):
    if n in ACYL: return 2
    if n in NOT_PEROX or n in SP.excluded or n in SP.dummy: return 0
    if n in PEROX_EXTRA or re.search(r'(O2|OO|OO[12])$', n): return 1
    return 0

# ------------------------------------------------------------------ photolysis: (reactant, J index) -> [(TUV column, scale)]
PH = {
 ('O3', 2): [(3, 1)], ('O3', 3): [(2, 1)], ('O2', 1): [(1, 1)], ('NO2', 11): [(6, 1)], ('H2O2', 9): [(5, 1)], ('MP', 10): [(31, 1)],
 ('CH2O', 7): [(17, 1)], ('CH2O', 8): [(18, 1)], ('HNO3', 16): [(13, 1)], ('HNO2', 15): [(12, 1)],
 ('HNO4', 17): [(14, .33)], ('HNO4', 18): [(14, .67)], ('NO3', 12): [(8, 1)], ('NO3', 13): [(7, 1)], ('N2O5', 14): [(10, 1)],
 ('ALD2', 61): [(19, 1), (21, 1)], ('ALD2', 62): [(20, 1)], ('PAN', 59): [(40, 1), (41, 1)], ('APAN', 59): [(40, 1), (41, 1)],
 ('ACR', 66): [(24, 1)], ('AROMCHO', 70): [(25, 1)], ('RCHO', 70): [(22, 1)], ('ACET', 76): [(26, 1)], ('ACET', 77): [(26, .1)],
 ('MEK', 69): [(28, 1)], ('GLYC', 68): [(23, 1)], ('GLYX', 72): [(44, 1)], ('GLYX', 74): [(45, 1)], ('MGLY', 71): [(46, 1)],
 ('MVK', 63): [(27, 1 / 3.)], ('MVK', 64): [(27, 1 / 3.)], ('MVK', 65): [(27, 1 / 3.)], ('MACR', 66): [(25, 1)], ('HAC', 75): [(29, 1), (30, 1)],
 ('MAP', 99): [(33, 1)], ('MPN', 103): [(14, .33)], ('MPN', 104): [(14, .67)], ('ATOOH', 97): [(31, 1)], ('N2O', 36): [(11, 1)],
 ('PRPN', 79): [(31, 1)], ('ETP', 80): [(31, 1)], ('RA3P', 81): [(31, 1)], ('RB3P', 82): [(31, 1)], ('R4P', 83): [(31, 1)],
 ('R7P', 83): [(31, 1)], ('ALK4P', 83): [(31, 1)], ('PP', 84): [(31, 1)], ('RP', 85): [(31, 1)],
 ('R4N2', 98): [(36, 1)], ('R7N2', 98): [(36, 1)], ('RNO3', 98): [(36, 1)], ('ALK4N2', 98): [(36, 1)],
 ('PIP', 105): [(31, 1)], ('ETHLN', 107): [(38, 1)], ('MONITS', 111): [(36, 1)], ('MONITU', 112): [(36, 1)], ('HONIT', 113): [(37, 1)],
 ('MENO3', 134): [(34, 1)], ('ETNO3', 135): [(35, 1)], ('IPRNO3', 136): [(36, 1)], ('NPRNO3', 137): [(36, 1)],
 ('HMHP', 86): [(32, 1)], ('HPETHNL', 87): [(31, 1)], ('PYAC', 88): [(48, 1)], ('PROPNN', 89): [(38, 1)],
 ('MVKHC', 90): [(29, 1), (30, 1)], ('MVKHCB', 91): [(29, 1), (30, 1)], ('MVKHP', 92): [(31, 1)], ('MVKPC', 93): [(25, 1)],
 ('MCRENOL', 94): [(25, 1)], ('MCRHP', 95): [(31, 1)], ('MACR1OOH', 96): [(33, 1)], ('MVKN', 108): [(37, 1)], ('MCRHN', 109): [(37, 1)],
 ('MCRHNB', 110): [(37, 1)], ('RIPA', 138): [(31, 1)], ('RIPB', 139): [(31, 1)], ('RIPC', 140): [(31, 1)], ('RIPD', 141): [(31, 1)],
 ('HPALD1', 142): [(51, 1)], ('HPALD2', 143): [(51, 1)], ('HPALD3', 144): [(51, 1)], ('HPALD4', 145): [(51, 1)],
 ('IHN1', 146): [(37, 1)], ('IHN2', 147): [(37, 1)], ('IHN3', 148): [(37, 1)], ('IHN4', 149): [(37, 1)],
 ('INPB', 150): [(31, 1)], ('INPD', 151): [(31, 1)], ('INPD', 152): [(37, 1)], ('ICN', 106): [(38, 1)], ('IDN', 78): [(37, 1)],
 ('ICPDH', 153): [(31, 1)], ('ICPDH', 154): [(31, 1)], ('IDHDP', 155): [(31, 1)], ('IDHPE', 156): [(31, 1)], ('IDCHP', 157): [(31, 1)],
 ('ITHN', 158): [(31, 1)], ('ITHN', 159): [(37, 1)], ('ITCN', 160): [(38, 1)], ('ITCN', 161): [(37, 1)],
 ('ETHP', 162): [(31, 1)], ('BALD', 163): [(22, 1)], ('BZCO3H', 164): [(33, 1)], ('BENZP', 165): [(31, 1)], ('NPHEN', 166): [(12, 1)],
 ('PPN', 167): [(42, 1), (43, 1)], ('APINP', 168): [(31, 1)], ('PINAL', 169): [(22, 1)], ('PINO3H', 170): [(33, 1)], ('PINONIC', 171): [(26, 1)],
 ('C96O2H', 172): [(31, 1)], ('BPINP', 173): [(31, 1)], ('BPINOOH', 174): [(31, 1)], ('LIMO3H', 175): [(33, 1)], ('LIMO2H', 176): [(31, 1)],
 ('LIMAL', 177): [(22, 1)],
}
# heterogeneous uptake: function -> gamma (fixed; GEOS-Chem computes it from aerosol composition / RH / T; see README)
HET_GAMMA = {'HO2uptk1stOrd': 0.2, 'N2O5uptkByH2O': 0.02, 'NO2uptk1stOrdAndCloud': 5e-6, 'NO3uptk1stOrdAndCloud': 5e-3,
             'GLYXuptk1stOrd': 2.2e-3, 'MGLYuptk1stOrd': 3.6e-7, 'IEPOXuptk1stOrd': 0.01}
HMML_GAMMA = 3e-4
EXCLUDE_FUNCS = re.compile(r'K_MT|K_CLD|K_MDL|HMS|SO2uptk|N2O5uptkByCloud|N2O5uptkByStratHCl|Ibrkdn|Iuptk|hypsis|SALA|SALC|IONO2|ClNO|BrNO|HOCl|HOBr|HBr|OHuptk')

# ------------------------------------------------------------------ expression parsing
def split_top(s, ch):
    out, d, cur = [], 0, ''
    for i, c in enumerate(s):
        if c == '(': d += 1
        elif c == ')': d -= 1
        if c == ch and d == 0:
            prev = cur.rstrip()[-1:] if cur.strip() else ''
            if ch == '+' and (prev in 'dDeE' and prev != '' and re.search(r'[0-9.][dDeE]$', cur.rstrip())): cur += c; continue
            out.append(cur); cur = ''
        else: cur += c
    out.append(cur)
    return [x.strip() for x in out if x.strip()]

def parse_term(t):
    """-> dict(scale, numden, exps[list of c], funcs[list of (name, [args])])"""
    d = dict(scale=1.0, numden=False, exps=[], funcs=[])
    for f in split_top(t, '*'):
        f = f.strip()
        if f.startswith('(') and f.endswith(')') and isnum(f[1:-1]): f = f[1:-1]
        if isnum(f): d['scale'] *= num(f)
        elif f.upper() == 'NUMDEN': d['numden'] = True
        else:
            m = re.fullmatch(r'(?i:exp)\(\s*(' + NUM + r')\s*/\s*TEMP\s*\)', f)
            if m: d['exps'].append(num(m.group(1))); continue
            m = re.fullmatch(r'(\w+)\((.*)\)', f)
            if not m: raise ValueError('cannot parse factor %r in %r' % (f, t))
            args = [a.strip() for a in split_top(m.group(2), ',')]
            d['funcs'].append((m.group(1), args))
    return d

class GCTranslator:
    def __init__(self): self.log = collections.Counter()
    def translate(self, r, ren=lambda n: n, counter=None):
        """-> (lines, reason).  lines == [] with a reason when the reaction is not translated."""
        names = [n for c, n in r['re']]
        react0 = [n for n in names if n not in ('hv', 'M')]
        for n in react0:
            if n in SP.excluded: return [], 'excluded reactant %s' % n
        if EXCLUDE_FUNCS.search(r['k']): return [], 'excluded process (cloud / sea-salt / aqueous / halogen)'
        for c, n in r['pr']:
            if n in SP.excluded and n not in AEROSOL: return [], 'gas-phase halogen product %s' % n
        prods = collections.OrderedDict()
        for c, n in r['pr']:
            if n in ('M', 'H2O', 'O2', 'N2', 'hv') or n in SP.excluded or n in SP.dummy or n.startswith('L') and n in SP.dummy: continue
            n = ren(n); prods[n] = prods.get(n, 0.0) + c
        pr = [(c, n) for n, c in prods.items()]
        pr_sink = pr if pr else [(1.0, 'M')]
        react = [ren(n) for n in react0]
        cnt = counter(react) if counter else None
        out = []
        def emit(re_, pr_, flag, params):
            p = list(pr_)
            if cnt and cnt not in [n for _, n in p]: p.append((1.0, cnt))
            out.append(line([(1.0, x) for x in re_], p, flag, params))
        hv = 'hv' in names
        k = r['k']
        # ---- photolysis
        m = re.fullmatch(r'PHOTOL\((\d+)\)', k.replace(' ', ''))
        if m:
            key = (react0[0], int(m.group(1)))
            if key not in PH: return [], 'no TUV analog for photolysis %s' % (key,)
            for col, sc in PH[key]: emit(react, pr_sink, 0, [col] + ([] if abs(sc - 1) < 1e-12 else [sc]))
            return out, None
        # ---- heterogeneous uptake (ASP first-order heterogeneous flag 3)
        m = re.match(r'(\w+uptk1stOrd\w*|N2O5uptkByH2O)\(', k)
        if m:
            f = m.group(1); mm = re.search(r'SR_MW\(ind_(\w+)\)', k)
            sp = ren(mm.group(1)) if mm else ('HO2' if f.startswith('HO2') else 'NO2' if f.startswith('NO2') else 'NO3' if f.startswith('NO3') else 'N2O5')
            raw = mm.group(1) if mm else {'HO2': 'HO2', 'NO2': 'NO2', 'NO3': 'NO3', 'N2O5': 'N2O5'}[sp]
            if f.startswith('VOCuptk'): g = num(re.search(r'SR_MW\(ind_\w+\)\s*,\s*(' + NUM + r')', k).group(1))
            elif f.startswith('IEPOX') and raw == 'HMML': g = HMML_GAMMA
            else: g = HET_GAMMA[f]
            hp = [(c, ren(n)) for c, n in [(c, n) for c, n in r['pr']] if n not in SP.excluded and n not in ('H2O', 'M', 'hv') and n not in SP.dummy]
            emit([sp], hp if hp else [(1.0, 'M')], 3, [g, round(SP.mw(raw), 2), 1e-5]); return out, None
        # ---- sum of terms
        terms = split_top(k, '+')
        for t in terms:
            self.term(t, react, react0, pr_sink, emit, r)
        return out, None
    def term(self, t, react, react0, pr, emit, r):
        d = parse_term(t); sc = d['scale']; nb = len(react)
        if nb == 3 and react.count('H2O') == 2:       # X + 2 H2O: one H2O is the reactant, the other is carried by flag 10 (H2O (k1 + M k2))
            rr = list(react); rr.remove('H2O'); assert d['funcs'] and d['funcs'][0][0] == 'GCARR_ac'
            v = [num(x) for x in d['funcs'][0][1]]
            emit(rr, pr, 10, [sc * v[0], v[1], 0.0, 0.0]); return
        assert nb in (1, 2), (r, nb)
        flagless = lambda params: emit(react, pr, None if nb == 2 else 2, params)
        # constants and Arrhenius
        if not d['funcs']:
            assert not d['numden'], t
            if not d['exps']: flagless([sc, 0.0, 0.0, 1.0]); return
            assert len(d['exps']) == 1, t
            flagless([sc, 0.0, d['exps'][0], 1.0]); return
        assert len(d['funcs']) == 1, t
        fn, a = d['funcs'][0]; v = [num(x) for x in a if isnum(x)]
        if fn == 'GCARR_ac': flagless([sc * v[0], 0.0, v[1], 1.0]); return
        if fn == 'GCARR_abc': flagless([sc * v[0] * 300.0 ** v[1], -v[1], v[2], 1.0]); return
        if fn == 'GCARR_ab':
            A = sc * v[0] * 300.0 ** v[1]
            if d['numden']: emit(react + ['M'], pr + [(1.0, 'M')], 2, [A, -v[1], 0.0, 1.0])     # termolecular: M A T^B
            else: flagless([A, -v[1], 0.0, 1.0])
            return
        assert not d['numden'], t
        # fall-off (JPL) families -------------------------------------------------------
        def troe(k0, n0, e0, ki, ni, ei, F, N, keq=None):
            if nb == 2 and keq is None:
                emit(react + ['M'], pr + [(1.0, 'M')], 4, [sc * k0, n0, e0, sc * ki, ni, ei, F, N])
            else:
                assert nb == 1, t
                emit(react, pr, 13, [sc * k0, n0, e0, sc * ki, ni, ei, F, N] + (keq if keq else [1.0, 0.0]))
        if fn == 'GCJPLPR_aa': troe(v[0], 0, 0, v[1], 0, 0, v[2], 1.0); return
        if fn == 'GCJPLPR_aba': troe(v[0], v[1], 0, v[2], 0, 0, v[3], 1.0); return
        if fn == 'GCJPLPR_abab': troe(v[0], v[1], 0, v[2], v[3], 0, v[4], 1.0); return
        if fn == 'GCJPLPR_abcabc': troe(v[0], v[1], v[2], v[3], v[4], v[5], v[6], 1.0); return
        if fn == 'GCJPLEQ_acabab':
            assert nb == 1; troe(v[2], v[3], 0, v[4], v[5], 0, v[6], 1.0, keq=[v[0], v[1]]); return
        if fn == 'GCJPLAC_ababac':   # k3 (1 - kf/ki), Troe N = 1:  [A0; n0; Ai; ni; F; A3; C3]
            assert nb == 2; emit(react, pr, 16, [sc * v[0], v[1], v[2], v[3], v[6], sc * v[4], v[5]]); return
        if fn in ('GC_PAN_acac', 'GC_PAN_abab'):
            cf = v[4]; N = 0.75 - 1.27 * math.log10(cf)
            if fn == 'GC_PAN_acac': k0, n0, e0, ki, ni, ei = v[0], -v[1], 0.0, v[2], -v[3], 0.0
            else: k0, n0, e0, ki, ni, ei = v[0], 0.0, v[1], v[2], 0.0, v[3]
            if fn == 'GC_PAN_abab': k0 = k0; # exp(b/T) in the rate: energy term
            troe(k0, n0, e0, ki, ni, ei, cf, N); return
        # branching / pressure-dependent organic forms -----------------------------------
        if fn == 'GC_RO2NO_A1_ac': flagless([sc * v[0] * 3e-4, 0.0, v[1], 1.0]); return
        if fn == 'GC_RO2NO_B1_ac': flagless([sc * v[0] * (1 - 3e-4), 0.0, v[1], 1.0]); return
        if fn in ('GC_RO2NO_A2_aca', 'GC_RO2NO_B2_aca'):
            emit(react, pr, 1 if fn.endswith('A2_aca') else 2, [sc * v[0], 0.0, v[1], 1.0, v[2], 1.0, 1.0, 0.0]); return
        if fn in ('GC_NIT', 'GC_ALK'): emit(react, pr, 12 if fn == 'GC_NIT' else 13, [sc * v[0]] + v[1:]); return
        if fn == 'GC_TBRANCH_1_acac': emit(react, pr, 8, [sc * v[0], v[1], v[2], v[3]]); return
        if fn == 'GC_RO2HO2_aca': flagless([sc * v[0] * (1 - math.exp(-0.245 * v[2])), 0.0, v[1], 1.0]); return
        if fn == 'GC_DMSOH_acac': emit(react, pr, 15, [sc * v[0] * 0.2095, v[1], v[2] * 0.2095, v[3]]); return
        if fn == 'GC_GLYXNO3_ac': emit(react, pr, 17, [sc * v[0], v[1], 3.5e18, 0.2095]); return
        if fn == 'GC_GLYCOH_A_a': emit(react, pr, 18, [sc * v[0], 0.0, 11.0729, 73.0, 0]); return
        if fn == 'GC_GLYCOH_B_a': emit(react, pr, 18, [sc * v[0], 0.0, 11.0729, 73.0, 1]); return
        if fn == 'GC_HACOH_A_ac': emit(react, pr, 18, [sc * v[0], v[1], 23.7, 60.0, 0]); return
        if fn == 'GC_HACOH_B_ac': emit(react, pr, 18, [sc * v[0], v[1], 23.7, 60.0, 1]); return
        if fn == 'GC_EPO_a': emit(react, pr, 14, [sc * v[0], v[1], v[2]]); return
        if fn == 'GC_HO2HO2_acac':
            emit(react, pr, 4, [sc * v[0], v[1], sc * v[2], v[3]])
            emit(react, pr, 10, [sc * v[0] * 1.4e-21, v[1] + 2200.0, sc * v[2] * 1.4e-21, v[3] + 2200.0]); return
        if fn == 'GC_ISO1': emit(react, pr, 20, [sc * v[0]] + v[1:]); return
        if fn == 'GC_ISO2': emit(react, pr, 21, [sc * v[0]] + v[1:]); return
        # linear-in-T Arrhenius, floored at zero: ASP flags 17 (1-body) / 19 (2-body) [A; B; C; D; E]
        if fn == 'ARRPLUS_ade': emit(react, pr, 17 if nb == 1 else 19, [sc * v[0], 0.0, 0.0, v[1], v[2]]); return
        if fn == 'ARRPLUS_abde': emit(react, pr, 17 if nb == 1 else 19, [sc * v[0], -v[1], 0.0, v[2], v[3]]); return
        if fn == 'TUNPLUS_abcde':
            assert nb == 1; emit(react, pr, 17, [sc * v[0], v[1], v[2], v[3], v[4]]); return
        raise ValueError('unhandled rate function %s in %r' % (fn, r['k']))

def main():
    tr = GCTranslator(); R = read_eqn(); n_ok = 0; why = collections.Counter(); nl = 0
    for r in R:
        try: lines, reason = tr.translate(r)
        except Exception as e: print('ERR', r['n'], r['k'], e); continue
        if lines: n_ok += 1; nl += len(lines)
        else: why[reason.split(' ')[0] + ' ' + reason.split(' ')[1] if reason else '?'] += 1
    print(len(R), 'reactions:', n_ok, 'translated ->', nl, 'ASP lines;', sum(why.values()), 'excluded'); print(why)
if __name__ == '__main__': main()
