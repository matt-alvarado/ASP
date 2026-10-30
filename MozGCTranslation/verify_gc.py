"""Check every translated GEOS-Chem rate against the KPP Fortran functions (fullchem_RateLawFuncs.F90 / rateLawUtilFuncs.F90, transcribed
below) over T = 220-320 K, P = 200-1013 mb.  The ASP side mirrors GasPhaseChemistry.h (common_asp flags).  Fixed reactants (O2, N2, H2O)
carry the same concentrations on both sides (ASP: O2 = 0.21 M, N2 = 0.79 M)."""
import math, re, sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gc_rates import *

def make_env(T, M, h2o):
    L10 = math.log10; ex = math.exp
    K300 = 300. / T; TK = T / 300.
    def GCARR_ab(a, b): return a * K300 ** b
    def GCARR_ac(a, c): return a * ex(c / T)
    def GCARR_abc(a, b, c): return a * ex(c / T) * K300 ** b
    def ARRPLUS_ade(a, d, e): return max(a * (d + T * e), 0.)
    def ARRPLUS_abde(a, b, d, e): return max(a * (d + T * e) * ex(-b / T), 0.)
    def TUNPLUS_abcde(a, b, c, d, e): return max(a * (d + T * e) * ex(b / T) * ex(c / T ** 3), 0.)
    def GC_ISO1(a, b, c, d, e, f, g):
        k0 = d * ex(e / T) * ex(1e8 / T ** 3); k1 = f * ex(g / T); return a * ex(b / T) * (1 - c * k0 / (k0 + k1))
    def GC_ISO2(a, b, c, d, e, f, g):
        k0 = d * ex(e / T) * ex(1e8 / T ** 3); k1 = f * ex(g / T); return a * ex(b / T) * c * k0 / (k0 + k1)
    def GC_EPO_a(a1, e1, m1): return a1 * ex(e1 / T) / (m1 * M + 1.)
    def pan(k0, k1, cf):
        kr = k0 / k1; nc = .75 - 1.27 * L10(cf); f = 10. ** (L10(cf) / (1. + (L10(kr) / nc) ** 2)); return k0 * k1 * f / (k0 + k1)
    def GC_PAN_abab(a0, b0, a1, b1, cf): return pan(a0 * ex(b0 / T) * M, a1 * ex(b1 / T), cf)
    def GC_PAN_acac(a0, c0, a1, c1, cf): return pan(a0 * TK ** c0 * M, a1 * TK ** c1, cf)
    def nitfrac(n):
        k0 = 2e-22 * ex(n) * M; k1 = .43 * (T / 298.) ** (-8); k1 = k0 / k1
        return (k0 / (1 + k1)) * .41 ** (1. / (1. + L10(k1) ** 2))
    def GC_NIT(a0, b0, c0, n, x0, y0):
        k2 = nitfrac(n); return max(a0 * (x0 - T * y0) * ex(b0 / T) * k2 / (k2 + c0), 0.)
    def GC_ALK(a0, b0, c0, n, x0, y0):
        k2 = nitfrac(n); return max(a0 * (x0 - T * y0) * ex(b0 / T) * c0 / (k2 + c0), 0.)
    def GC_HO2HO2_acac(a0, c0, a1, c1): return (a0 * ex(c0 / T) + a1 * ex(c1 / T) * M) * (1. + 1.4e-21 * h2o * ex(2200. / T))
    def GC_TBRANCH_1_acac(a0, c0, a1, c1): return a0 * ex(c0 / T) / (1. + a1 * ex(c1 / T))
    def GC_RO2HO2_aca(a0, c0, a1): return a0 * ex(c0 / T) * (1. - ex(-.245 * a1))
    def GC_DMSOH_acac(a0, c0, a1, c1): return a0 * ex(c0 / T) * M * .2095 / (1. + a1 * ex(c1 / T) * .2095)
    def GC_GLYXNO3_ac(a0, c0):
        O2 = M * .2095; return a0 * ex(c0 / T) * (O2 + 3.5e18) / (2. * O2 + 3.5e18)
    def frac(a, b): return max(1. - a * ex(-T / b), 0.)
    def GC_GLYCOH_A_a(a0): return a0 * frac(11.0729, 73.)
    def GC_GLYCOH_B_a(a0): return a0 * (1. - frac(11.0729, 73.))
    def GC_HACOH_A_ac(a0, c0): return a0 * ex(c0 / T) * frac(23.7, 60.)
    def GC_HACOH_B_ac(a0, c0): return a0 * ex(c0 / T) * (1. - frac(23.7, 60.))
    def GC_RO2NO_A1_ac(a0, c0): return a0 * ex(c0 / T) * 3e-4
    def GC_RO2NO_B1_ac(a0, c0): return a0 * ex(c0 / T) * (1. - 3e-4)
    def ro2no(a0, c0, a1):
        k0 = a0 * ex(c0 / T); xxyn = 1.94e-22 * ex(.97 * a1) * M; yyyn = .826 * (300. / T) ** 8.1
        aaa = L10(xxyn / yyyn); zzyn = 1. / (1. + aaa * aaa); rarb = (xxyn / (1. + xxyn / yyyn)) * .411 ** zzyn
        return k0, rarb / (1. + rarb)
    def GC_RO2NO_A2_aca(a0, c0, a1): k0, f = ro2no(a0, c0, a1); return k0 * f
    def GC_RO2NO_B2_aca(a0, c0, a1): k0, f = ro2no(a0, c0, a1); return k0 * (1. - f)
    def jpl(rlow, rhigh, fv):
        xy = rlow / rhigh; return rlow * fv ** (1. / (1. + L10(xy) ** 2)) / (1. + xy)
    def GCJPLPR_aa(a1, a2, fv): return jpl(a1 * M, a2, fv)
    def GCJPLPR_aba(a1, b1, a2, fv): return jpl(a1 * K300 ** b1 * M, a2, fv)
    def GCJPLPR_abab(a1, b1, a2, b2, fv): return jpl(a1 * K300 ** b1 * M, a2 * K300 ** b2, fv)
    def GCJPLPR_abcabc(a1, b1, c1, a2, b2, c2, fv): return jpl(a1 * K300 ** b1 * ex(c1 / T) * M, a2 * K300 ** b2 * ex(c2 / T), fv)
    def GCJPLEQ_acabab(a0, c0, a1, b1, a2, b2, fv): return GCJPLPR_abab(a1, b1, a2, b2, fv) / (a0 * ex(c0 / T))
    def GCJPLAC_ababac(a1, b1, a2, b2, a3, c3, fv):
        rlow = a1 * K300 ** b1 * M; rhigh = a2 * K300 ** b2; return a3 * ex(c3 / T) * (1. - jpl(rlow, rhigh, fv) / rhigh)
    e = {k: v for k, v in locals().items() if k[:2] in ('GC', 'AR', 'TU')}
    e.update(TEMP=T, NUMDEN=M, exp=ex, EXP=ex, H2O=h2o)
    return e

def asp_rate(l, T, M, h2o):
    ex = math.exp; L10 = math.log10
    head, rest = l.split(' ; ', 1)
    reac = [x.strip() for x in head.split('=>')[0].split('+')]
    nums = [float(t) for t in rest.strip().rstrip(';').split(';') if t.strip()]
    nbod = len(reac); has_M = 'M' in reac
    def falloff(p, keq=False):
        k0 = p[0] * (T / 300.) ** -p[1] * ex(p[2] / T) * M; ki = p[3] * (T / 300.) ** -p[4] * ex(p[5] / T)
        return (k0 * ki) / (k0 + ki) * p[6] ** (1. / (1. + (L10(k0 / ki) / p[7]) ** 2))
    if nbod == 1:
        f = int(nums[0]); p = nums[1:]
        if f == 2: return p[0] * T ** p[1] * ex(p[2] / T) ** p[3]
        if f == 13: return falloff(p) / (p[8] * ex(p[9] / T))
        if f == 17: return max(p[0] * (p[3] + T * p[4]) * ex(p[1] / T) * ex(p[2] / T ** 3), 0.)
        if f == 0: raise KeyError('photolysis')
        raise ValueError(l)
    if has_M:
        f = int(nums[0]); p = nums[1:]
        if f == 2: return M * p[0] * T ** p[1] * ex(p[2] / T) ** p[3]
        if f == 4: return falloff(p)
        raise ValueError(l)
    f = int(nums[0]); p = nums[1:]
    if f == 0: return nums[0] * T ** nums[1] * ex(nums[2] / T) ** nums[3]          # flagless: the parser reads INT(first number) = 0
    if f in (1, 2):
        k0 = p[0] * T ** p[1] * ex(p[2] / T) ** p[3]
        xxyn = 1.94e-22 * ex(.97 * p[4]) * M; yyyn = .826 * (T / 300.) ** -8.1
        z = 1. / (1. + L10(xxyn / yyyn) ** 2); ratio = p[5] * (xxyn / (1. + xxyn / yyyn)) * .411 ** z
        ratio = p[6] * ratio + p[7]; sc = ratio / (1. + ratio)
        return k0 * (sc if f == 1 else 1. - sc)
    if f in (4, 8, 10, 15): k1 = p[0] * ex(p[1] / T); k2 = p[2] * ex(p[3] / T)
    if f == 4: return k1 + k2 * M
    if f == 8: return k1 / (1. + k2)
    if f == 10: return h2o * (k1 + M * k2)
    if f in (12, 13):
        k0 = 2e-22 * ex(p[3]) * M; ki = .43 * (T / 298.) ** -8; kk = k0 / ki
        k2_ = (k0 / (1. + kk)) * .41 ** (1. / (1. + L10(kk) ** 2))
        k3 = k2_ / (k2_ + p[2]) if f == 12 else p[2] / (k2_ + p[2])
        return max(p[0] * (p[4] - T * p[5]) * ex(p[1] / T) * k3, 0.)
    if f == 14: return p[0] * ex(p[1] / T) / (p[2] * M + 1.)
    if f == 15: return M * k1 / (1. + k2)
    if f == 16:
        k0 = p[0] * (T / 300.) ** -p[1] * M; ki = p[2] * (T / 300.) ** -p[3]
        kf = (k0 * ki) / (k0 + ki) * p[4] ** (1. / (1. + L10(k0 / ki) ** 2)); return p[5] * ex(p[6] / T) * (1. - kf / ki)
    if f == 17: return p[0] * ex(p[1] / T) * (p[3] * M + p[2]) / (2. * p[3] * M + p[2])
    if f == 18:
        fr = max(1. - p[2] * ex(-T / p[3]), 0.)
        if p[4] > .5: fr = 1. - fr
        return p[0] * ex(p[1] / T) * fr
    if f == 19: return max(p[0] * (p[3] + T * p[4]) * ex(p[1] / T) * ex(p[2] / T ** 3), 0.)
    if f in (20, 21):
        k0 = p[3] * ex(p[4] / T) * ex(1e8 / T ** 3); k1 = p[5] * ex(p[6] / T); k2 = p[2] * k0 / (k0 + k1)
        return p[0] * ex(p[1] / T) * ((1. - k2) if f == 20 else k2)
    raise ValueError(l)

def main():
    tr = GCTranslator(); R = read_eqn(); worst = []; n = 0; seen_flags = collections.Counter(); nrx = 0
    kB = 1.380649e-23
    for T in (220., 250., 280., 300., 320.):
        for P in (200., 600., 1013.):
            M = P * 100. / (kB * T) * 1e-6; h2o = 0.01 * M * 0.5
            env = make_env(T, M, h2o)
            fixed = dict(O2=.21 * M, N2=.79 * M, H2O=h2o)
            for r in R:
                if re.search(r'PHOTOL|uptk|SR_MW', r['k']): continue
                lines, why = tr.translate(r)
                if not lines: continue
                k = re.sub(r'(\d)[dD]([+-]?\d)', r'\1e\2', r['k']).replace('_dp', '')
                ref = eval(k, {}, env)
                for c, nme in r['re']:
                    if nme in fixed: ref *= fixed[nme]
                tot = 0.
                for l in lines:
                    reac = [x.strip() for x in l.split(' ; ')[0].split('=>')[0].split('+')]
                    kk = asp_rate(l, T, M, h2o)
                    for x in reac:
                        if x in fixed: kk *= fixed[x]
                    tot += kk
                n += 1
                err = abs(tot - ref) / max(abs(ref), 1e-300) if ref != 0 or tot != 0 else 0.
                if err > 2e-3: worst.append((err, r['n'], r['k'], T, P, tot, ref, lines))
    worst.sort(key=lambda x: -x[0])
    print('checked', n, 'rate evaluations;', len(worst), 'outside 0.2%')
    s = set()
    for w in worst:
        if w[1] in s: continue
        s.add(w[1]); print('%.3g' % w[0], w[1], w[2][:70], w[3], w[4], w[5], w[6], w[7]); 
        if len(s) > 15: break
if __name__ == '__main__': main()
