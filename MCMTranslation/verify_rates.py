"""Check every translated ASP rate against the MCM Fortran expression over a T, P, H2O, RO2 grid (independent evaluation)."""
import math, re, sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mcm2asp import *

def ref_env(T, M, h2o, ro2):
    env = dict(TEMP=T, M=M, N2=FN2 * M, O2=FO2 * M, H2O=h2o, RO2=ro2, EXP=math.exp, LOG10=math.log10, cos=math.cos)
    txt = open(D + 'mcm_3-3-1_kpp_constants_mcm.f90').read()
    for ln in txt.split('SUBROUTINE define_constants_mcm()')[1].split('END SUBROUTINE')[0].split('\n'):
        m = re.match(r'^\s*([A-Z0-9_]+)\s*=\s*(.*?)\s*(?:!.*)?$', ln)
        if m and not m.group(1).startswith('J('):
            env[m.group(1)] = eval(fort_(m.group(2)), {}, env)
    for j in MCM_J: env[j] = j0(j)
    env['J'] = lambda v: v
    return env
def fort_(s): return s.replace('D+', 'E+')

def asp_k(flag, p, T, M, h2o, nrb):
    if flag is None: return p[0] * T ** p[1] * math.exp(p[2] / T) ** p[3] if nrb == 2 else None
    f = int(flag)
    return None

def parse_line(l):
    head, rest = l.split(' ; ', 1)
    reac = [x.strip() for x in head.split('=>')[0].split('+')]
    toks = [t.strip() for t in rest.strip().rstrip(';').split(';')]
    return reac, [float(t) for t in toks]

def asp_rate(l, T, M, h2o, ro2, conc):
    reac, nums = parse_line(l)
    nb = len([r for r in reac if r != 'M'])
    has_M = 'M' in reac
    nbod = len(reac)
    # decide flag: first number is the flag when the line has a flag.  Flagless 2-body: 4 numbers, first is A.
    # (flags are always integers followed by params; flagless lines have exactly 4 numbers)
    if nbod == 1:
        flag = int(nums[0]); p = nums[1:]
        if flag == 0: return conc['J'](int(p[0]), p[1] if len(p) > 1 else 1.0)
        if flag == 2: return p[0] * T ** p[1] * math.exp(p[2] / T) ** p[3]
        if flag == 12: return ro2 * p[0] * T ** p[1] * math.exp(p[2] / T) ** p[3]
        if flag == 14: return p[0] * T ** p[1] * math.exp(p[2] / T) * math.exp(p[3] / T ** 3)
        if flag == 8: return M * (0.21 * 3.2e-11 * math.exp(67. / T) + 0.79 * 2.0e-11 * math.exp(130. / T))
        if flag == 9:
            n2, o2 = .79 * M, .21 * M
            return (5.6e-34 * (T / 300.) ** -2.6 * n2 + 6.0e-34 * (T / 300.) ** -2.6 * o2) * o2
        if flag == 13:
            k0 = p[0] * (T / 300.) ** -p[1] * math.exp(p[2] / T) * M; ki = p[3] * (T / 300.) ** -p[4] * math.exp(p[5] / T)
            return (k0 * ki) / (k0 + ki) * p[6] ** (1. / (1. + (math.log10(k0 / ki) / p[7]) ** 2)) / (p[8] * math.exp(p[9] / T))
    if nbod == 3 or has_M:
        flag = int(nums[0]); p = nums[1:]
        if flag == 2: k = M * p[0] * T ** p[1] * math.exp(p[2] / T) ** p[3]
        elif flag == 4:
            k0 = p[0] * (T / 300.) ** -p[1] * math.exp(p[2] / T) * M; ki = p[3] * (T / 300.) ** -p[4] * math.exp(p[5] / T)
            k = (k0 * ki) / (k0 + ki) * p[6] ** (1. / (1. + (math.log10(k0 / ki) / p[7]) ** 2))
        elif flag == 5:
            k0 = p[0] * (T / 300.) ** -p[1] * math.exp(p[2] / T) * M; ki = p[3] * (T / 300.) ** -p[4] * math.exp(p[5] / T)
            Fc = p[6] * math.exp(-p[7] / T) + math.exp(-T / p[8]); N = 0.75 - 1.27 * math.log10(Fc)
            k = (k0 * ki) / (k0 + ki) * Fc ** (1. / (1. + (math.log10(k0 / ki) / N) ** 2))
        else: raise ValueError(l)
        return k
    # two-body
    if len(nums) == 4: return nums[0] * T ** nums[1] * math.exp(nums[2] / T) ** nums[3]
    flag = int(nums[0]); p = nums[1:]
    k1 = p[0] * math.exp(p[1] / T); k2 = p[2] * math.exp(p[3] / T)
    if flag == 4: return k1 + k2 * M
    if flag == 5:
        k3 = p[4] * math.exp(p[5] / T); return k1 + k3 * M / (1 + k3 * M / k2) if k2 else k1
    if flag == 7: return k1 * (1 - 1 / (1 + k2))
    if flag == 8: return k1 * (1 / (1 + k2))
    if flag == 9: return k1 * (1 - k2)
    if flag == 10: return h2o * (k1 + M * k2)
    raise ValueError(l)

def main():
    T_ = MCMTranslator(); R = read_eqn()
    worst = []; nchk = 0
    # J placeholders: TUV column value = sum of the MCM J's that column represents -> compare with j0(target)
    colval = {}
    for tgt, (col, ref) in J_MAP.items(): colval[col] = sum(j0(r) for r in ref)
    conc = dict(J=lambda col, sc=1.0: colval[col] * sc)
    for T in (230., 250., 270., 288., 300., 320.):
        for P in (250., 600., 1013.):
            M = 2.46e19 * (298. / T) * (P / 1013.); h2o = 0.01 * M * 0.5; ro2 = 3e9
            env = ref_env(T, M, h2o, ro2)
            for r in R:
                lines = T_.translate(r)
                if not lines: continue
                ref = eval(fort_(r['k']), {}, env)
                # concentrations: unit for MCM species; pools / background as set
                refrate = ref
                tot = 0.0
                for l in lines:
                    reac, nums = parse_line(l)
                    k = asp_rate(l, T, M, h2o, ro2, conc)
                    c = 1.0
                    for x in reac:
                        if x == 'O2': c *= FO2 * M
                        elif x == 'H2O': c *= h2o
                        elif x == 'RO2_T': c *= ro2
                    tot += k * c
                # reference: KPP rate is per unit reactant concentration product of MCM species only (O2/H2O/RO2 inside expr)
                # MCM reactions that list RO2-pool or M inside rate expressions are compared as-is
                err = abs(tot - refrate) / max(abs(refrate), 1e-300)
                nchk += 1
                if err > 2e-3: worst.append((err, r['n'], r['k'], T, P, lines))
    worst.sort(key=lambda x: -x[0])
    print('checked', nchk, 'rate evaluations;', len(worst), 'outside 0.2%')
    seen = set()
    for w in worst:
        if w[1] in seen: continue
        seen.add(w[1]); print('%.3g' % w[0], w[1], w[2], w[3], w[4], w[5])
        if len(seen) > 15: break
    return worst
if __name__ == '__main__': main()
