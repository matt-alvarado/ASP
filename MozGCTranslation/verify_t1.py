"""Independent check of every translated T1 rate against the WRF-Chem KPP Fortran expressions (transcribed from t1_mozcart.def and
WRFUserRateLaws.f90) over T = 220-320 K, P = 200-1013 mb, several H2O.  The ASP side mirrors GasPhaseChemistry.h."""
import math, re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MCMTranslation'))
import verify_rates as vr
from mozt1 import *

def ARR2(A, C, TEMP): return A * math.exp(-C / TEMP)
def JPL_TROE(k0, n, ki, m, base, T, cair):
    z = 300. / T; k0T = k0 * z ** n * cair; kiT = ki * z ** m; kr = k0T / kiT
    return k0T / (1 + kr) * base ** (1. / (1. + math.log10(kr) ** 2))
def TROEE(A, B, k0, n, ki, m, T, cair): return A * math.exp(-B / T) * JPL_TROE(k0, n, ki, m, 0.6, T, cair)
def make_env(T, cm, ch2o):
    e = dict(TEMP=T, temp=T, C_M=cm, c_m=cm, c_h2o=ch2o, ARR2=ARR2, JPL_TROE=JPL_TROE)
    def usr_O_O2(temp): return 6.00e-34 * (temp / 300.) ** (-2.4)
    def usr_O_O_M(temp): return 2.76e-34 * math.exp(720.0 / temp)
    def usr_HO2_HO2(temp, c_m, c_h2o):
        if c_h2o > 0: return (2.3e-13 * math.exp(600. / temp) + 1.7e-33 * c_m * math.exp(1000. / temp)) * (1. + 1.4e-21 * c_h2o * math.exp(2200. / temp))
        return 0.0
    def usr_HNO3_OH(temp, c_m):
        k0 = c_m * 6.5e-34 * math.exp(1335. / temp); k2 = math.exp(2199. / temp)
        k0 = k0 / (1.0 + k0 / (2.7e-17 * k2)); k2 = math.exp(460. / temp)
        return k0 + 2.4e-14 * k2
    def usr_HO2NO2_M(temp, c_m): return TROEE(4.76e26, 10900., 1.8e-31, 3.2, 4.7e-12, 1.4, temp, c_m) / c_m
    def usr_N2O5_M(temp, c_m): return TROEE(3.333e26, 10900., 2.2e-30, 4.4, 1.4e-12, .7, temp, c_m) / c_m
    def usr_PAN_M(temp, c_m): return TROEE(1.111e28, 14000., 8.5e-29, 6.5, 1.1e-11, 0., temp, c_m) / c_m
    def usr_PBZNIT_M(temp, c_m): return TROEE(1.111e28, 14000., 9.7e-29, 5.6, 9.3e-12, 0., temp, c_m) / c_m
    def usr_MPAN_M(temp, c_m): return 1.2221e17 * 300. * math.exp(-14000. / temp) / (temp * c_m)
    def usr_MCO3_NO2(temp, c_m): return 1.1e-11 * 300. / (temp * c_m)
    def usr_CH3COCH3_OH(temp): return 3.82e-11 * math.exp(-2000. / temp) + 1.33e-13
    def usr_DMS_OH(temp, c_m):
        wrk = .21 * c_m; ko = 1. + 5.5e-31 * math.exp(7460. / temp) * wrk
        return 1.7e-42 * math.exp(7810. / temp) * wrk / ko
    def usr_SO2_OH(temp, c_m):
        fc = 3.e-31 * (300. / temp) ** 3.3; wrk = fc * c_m; k0 = wrk / (1. + wrk / 1.5e-12)
        return k0 * .6 ** (1. / (1. + (math.log10(wrk / 1.5e-12)) ** 2))
    def usr_CO_OH_a(temp, c_m): return 1.5e-13 * (1. + 6.e-7 * 1.38044e-16 * c_m * temp)
    for f in (usr_O_O2, usr_O_O_M, usr_HO2_HO2, usr_HNO3_OH, usr_HO2NO2_M, usr_N2O5_M, usr_PAN_M, usr_PBZNIT_M, usr_MPAN_M,
              usr_MCO3_NO2, usr_CH3COCH3_OH, usr_DMS_OH, usr_SO2_OH, usr_CO_OH_a): e[f.__name__] = f
    return e

def asp_rate(l, T, M, h2o, P, colval):
    reac, nums = vr.parse_line(l)
    if len(reac) == 2 and 'M' not in reac and len(nums) == 6 and int(nums[0]) == 11:       # 2-body flag 11 (CO + OH form)
        A, B, C, D, F = nums[1:]
        return A * math.exp(B / T) * (C / T + D * P * 100.0 + F) / 100.0
    conc = dict(J=lambda col, sc=1.0: colval[col] * sc)
    return vr.asp_rate(l, T, M, h2o, 0.0, conc)

def main():
    tr = T1Translator(); R = read_eqn(); worst = []; n = 0
    colval = {c: 1e-3 * (c + 1) ** 1.3 for c in range(1, 52)}
    kB = 1.380649e-23
    for T in (220., 250., 280., 300., 320.):
        for P in (200., 600., 1013.):
            M = P * 100. / (kB * T) * 1e-6; h2o = 0.01 * M * 0.5
            env = make_env(T, M, h2o)
            fixed = dict(O2=0.21 * M, N2=0.79 * M, H2O=h2o, M=M)
            for r in R:
                if re.match(r'usr_\w+_aer', r['k']): continue
                lines = tr.translate(r)
                kk = r['k'].replace('_dp', '').replace('d', 'e') if False else r['k']
                pj = re.fullmatch(r'(?:(\S+?)\*)?j\(Pj_(\w+)\)', kk.replace(' ', ''))
                if pj:
                    sc = float(pj.group(1).replace('_dp', '')) if pj.group(1) else 1.0
                    ref = sc * sum(colval[c] for c in PJ[pj.group(2)])
                    # reference (k per reactant concentration); sum of ASP lines must equal it
                    tot = sum(asp_rate(l, T, M, h2o, P, colval) for l in lines)
                    if len(PJ[pj.group(2)]) > 1: pass
                    err = abs(tot - ref) / ref
                else:
                    ref = eval(kk.replace('_dp', ''), {}, env)
                    # KPP rate = k * product of reactant concentrations; fixed species (M, N2, O2, H2O) carry their concentration
                    f_ref = 1.0
                    for c, nme in r['re']:
                        if nme in fixed: f_ref *= fixed[nme]
                    tot = 0.0
                    for l in lines:
                        reac, _ = vr.parse_line(l)
                        k = asp_rate(l, T, M, h2o, P, colval)
                        # ASP lines: M is part of the flag formula (3-body) or explicit (never in 1/2-body lines)
                        f = 1.0
                        for x in reac:
                            if x in ('O2', 'N2', 'H2O'): f *= fixed[x]
                        if len(reac) == 3 and 'M' in reac: f_asp_M = 1.0
                        tot += k * f
                    refrate = ref * f_ref
                    # reference M factor: reactions with M among reactants and 3-body flags treat M inside k; first-order flag 13/2 absorb M
                    err = abs(tot - refrate) / max(abs(refrate), 1e-300)
                n += 1
                if err > 2e-3: worst.append((err, r['n'], r['k'], T, P, lines))
    worst.sort(key=lambda x: -x[0])
    print('checked', n, 'rate evaluations;', len(worst), 'outside 0.2%')
    seen = set()
    for w in worst:
        if w[1] in seen: continue
        seen.add(w[1]); print('%.3g' % w[0], w[1], w[2], w[3], w[4], w[5])
        if len(seen) > 20: break
    return worst
if __name__ == '__main__': main()
