#!/usr/bin/env python3
"""MOZART-T1 (WRF-Chem KPP T1_MOZCART, tag MZ148_T1mozcart_20170323) -> ASP reactions.

Every rate is expressed through structural ASP rate flags (no rate evaluated at a fixed T or P).
Aerosol uptake (usr_*_aer) is expressed through ASP's first-order heterogeneous flag 3 (gamma, molecular weight),
which uses the model's own particle surface area.
"""
import math, os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aspfmt import line
from t1_kpp import read_eqn

# TUV columns of the baseline 51-column PhotoRates.in that carry each T1 photolysis label (T1 jmap -> ASP columns).
# Where T1 sums several TUV channels the reaction is split into one ASP reaction per column (rates add).
PJ = {'h2o2': [5], 'o2': [1], 'o33p': [3], 'o31d': [2], 'hno3': [13], 'hno4': [14], 'n2o': [11], 'n2o5': [10], 'no2': [6],
      'no3o': [7, 8], 'ch3o2h': [31], 'ch2or': [17], 'ch2om': [18], 'ch3cho': [19, 20, 21], 'ch3coch3': [26],
      'ch3cocho': [46], 'glyald': [23], 'gly': [44, 45], 'hyac': [29, 30], 'macr': [25], 'mvk': [27], 'pan': [40, 41]}
# heterogeneous uptake: function -> (gamma, c_coefficient a  [c = a sqrt(T) cm/s]  or ('M', mw)), products
R_GAS = 8.314462
def mw_from_a(a):                       # c = sqrt(8 R T / (pi M)) = a sqrt(T) (cm/s)
    return 8 * R_GAS * 1e3 / (math.pi * (a / 100.) ** 2)
HET = {'usr_HO2_aer': (0.2, 2.53e3), 'usr_N2O5_aer': (0.1, 1.40e3), 'usr_NO2_aer': (1e-4, 2.15e3),
       'usr_NO3_aer': (1e-3, 1.85e3), 'usr_HONITR_aer': (.005, 1.26e3), 'usr_ISOPNIT_aer': (.005, 1.20e3),
       'usr_ONITR_aer': (.005, 1.20e3), 'usr_TERPNIT_aer': (.01, .992e3), 'usr_NC4CH2OH_aer': (.005, 1.20e3),
       'usr_NC4CHO_aer': (.005, 1.21e3), 'usr_NTERPOOH_aer': (.01, .957e3), 'usr_GLYOXAL_aer': (2e-4, 1.455e4 / math.sqrt(58.))}
DROP_PROD = {'H2O', 'O2', 'N2', 'hv'}
NUM = r'([0-9]*\.?[0-9]+(?:[eEdD][+-]?[0-9]+)?)(?:_dp)?'
def num(s): return float(s.replace('_dp', '').replace('d', 'e').replace('D', 'e'))

def troee_flag13(A, B, k0, n, ki, m):   # TROEE(A,B,...)/c_m with M a reactant -> 1-body fall-off / Keq,  Keq = exp(B/T)/A
    return [k0, n, 0.0, ki, m, 0.0, 0.6, 1.0, 1.0 / A, B]

class T1Translator:
    def __init__(self):
        self.skipped = []
    def translate(self, r, ren=lambda n: n, counter=None):
        k = r['k']
        reac = [n for c, n in r['re'] if n != 'hv']
        hv = any(n == 'hv' for c, n in r['re'])
        has_M = 'M' in reac
        react = [ren(n) for n in reac if n != 'M']
        prods = collections.OrderedDict()
        for c, n in r['pr']:
            if n in DROP_PROD or n == 'M': continue
            n = ren(n); prods[n] = prods.get(n, 0.0) + c
        pr = [(c, n) for n, c in prods.items()]
        pr_sink = pr if pr else [(1.0, 'M')]
        cnt = counter(react) if counter else None
        def add(l, c=cnt): return l
        out = []
        def emit(re_, pr_, flag, params):
            p = list(pr_)
            if cnt and cnt not in [n for _, n in p]: p.append((1.0, cnt))
            out.append(line([(1.0, x) for x in re_], p, flag, params))
        # ---- photolysis
        m = re.fullmatch(r'(?:' + NUM + r'\*)?j\(Pj_(\w+)\)', k.replace(' ', ''))
        if m:
            sc = num(m.group(1)) if m.group(1) else 1.0
            for col in PJ[m.group(2)]:
                emit(react, pr_sink, 0, [col] + ([] if abs(sc - 1) < 1e-12 else [sc]))
            return out
        # ---- heterogeneous uptake (ASP first-order heterogeneous flag 3)
        fn = re.match(r'(usr_\w+)\(', k)
        if fn and fn.group(1) in HET:
            g, a = HET[fn.group(1)]
            emit(react, pr_sink, 3, [g, round(mw_from_a(a), 1), 1e-5]); return out
        # ---- constant
        m = re.fullmatch(NUM, k.replace(' ', ''))
        if m:
            c = num(m.group(1)); nb = len(react)
            emit(react, pr_sink, 2 if nb == 1 else None, [c, 0.0, 0.0, 1.0]); return out
        # ---- ARR2(A, C, TEMP) = A exp(-C/T)
        m = re.fullmatch(r'ARR2\(\s*' + NUM + r',\s*' + NUM.replace('(', '(-?', 1) + r',\s*TEMP\s*\)', k)
        if m:
            A, C = num(m.group(1)), num(m.group(2)); nb = len(react)
            assert nb in (1, 2), r
            emit(react, pr_sink, 2 if nb == 1 else None, [A, 0.0, -C, 1.0]); return out
        # ---- JPL_TROE
        m = re.fullmatch(r'JPL_TROE\(\s*' + '\\s*,\\s*'.join([NUM, NUM.replace('(', '(-?', 1), NUM, NUM.replace('(', '(-?', 1), NUM]) +
                         r'\s*,\s*TEMP\s*,\s*C_M\s*\)', k)
        if m:
            k0, n, ki, mm, base = [num(x) for x in m.groups()]
            assert len(react) == 2 and not has_M, r
            emit(react + ['M'], pr + [(1.0, 'M')], 4, [k0, n, 0.0, ki, mm, 0.0, base, 1.0]); return out
        # ---- named user functions
        f = fn.group(1) if fn else k
        if f == 'usr_O_O2':                  # 6.0e-34 (T/300)^-2.4 [O][O2][M]
            emit(react + ['M'], pr + [(1.0, 'M')], 2, [6.0e-34 * 300 ** 2.4, -2.4, 0.0, 1.0]); return out
        if f == 'usr_O_O_M':                 # 2.76e-34 exp(720/T) [O]^2 [M]
            emit(react + ['M'], pr + [(1.0, 'M')], 2, [2.76e-34, 0.0, 720.0, 1.0]); return out
        if f == 'usr_HO2_HO2':               # (ko + kinf)(1 + 1.4e-21 [H2O] exp(2200/T)) ; H2O in Fc, not a reactant
            emit(react, pr_sink, 4, [2.3e-13, 600., 1.7e-33, 1000.])
            emit(react, pr_sink, 10, [2.3e-13 * 1.4e-21, 600. + 2200., 1.7e-33 * 1.4e-21, 1000. + 2200.]); return out
        if f == 'usr_HNO3_OH':               # identical to MCM KMT11
            emit(react, pr_sink, 5, [2.4e-14, 460., 2.7e-17, 2199., 6.5e-34, 1335.]); return out
        if f == 'usr_HO2NO2_M': emit(react, pr_sink, 13, troee_flag13(4.76e26, 10900., 1.8e-31, 3.2, 4.7e-12, 1.4)); return out
        if f == 'usr_N2O5_M': emit(react, pr_sink, 13, troee_flag13(3.333e26, 10900., 2.2e-30, 4.4, 1.4e-12, .7)); return out
        if f == 'usr_PAN_M': emit(react, pr_sink, 13, troee_flag13(1.111e28, 14000., 8.5e-29, 6.5, 1.1e-11, 0.)); return out
        if f == 'usr_PBZNIT_M': emit(react, pr_sink, 13, troee_flag13(1.111e28, 14000., 9.7e-29, 5.6, 9.3e-12, 0.)); return out
        if f == 'usr_MPAN_M':                # 1.2221e17*300 exp(-14000/T)/T  (M reactant cancels c_m)
            emit(react, pr_sink, 2, [1.2221e17 * 300., -1.0, -14000., 1.0]); return out
        if f == 'usr_MCO3_NO2':              # 1.1e-11*300/(T c_m) [M] -> 3.3e-9 / T
            emit(react, pr_sink, None, [1.1e-11 * 300., -1.0, 0.0, 1.0]); return out
        if f == 'usr_CH3COCH3_OH':           # 3.82e-11 exp(-2000/T) + 1.33e-13
            emit(react, pr_sink, None, [3.82e-11, 0.0, -2000., 1.0])
            emit(react, pr_sink, None, [1.33e-13, 0.0, 0.0, 1.0]); return out
        if f == 'usr_DMS_OH':                # a M/(1+bM): a = .21*1.7e-42 e^{7810/T}, b = .21*5.5e-31 e^{7460/T}
            emit(react, pr_sink, 5, [0.0, 0.0, 1.7e-42 / 5.5e-31, 7810. - 7460., 0.21 * 1.7e-42, 7810.]); return out
        if f == 'usr_SO2_OH':                # fall-off, k0 = 3e-31 (300/T)^3.3 M, kinf = 1.5e-12, Fc = 0.6, N = 1
            emit(react + ['M'], pr + [(1.0, 'M')], 4, [3.0e-31, 3.3, 0.0, 1.5e-12, 0.0, 0.0, 0.6, 1.0]); return out
        if f == 'usr_CO_OH_a':               # 1.5e-13 (1 + 6e-7 k_B M T) = 1.5e-13 (1 + 6e-6 P[Pa]);  ASP flag 11: A [C/T + D P[mb] 100 + F]/100
            emit(react, pr_sink, 11, [1.5e-11, 0.0, 0.0, 6e-6, 1.0]); return out
        raise ValueError('unhandled rate: %s' % k)
