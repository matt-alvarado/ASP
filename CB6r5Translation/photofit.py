"""Photolysis rates for CMAQ photolysis labels that have no TUV column in ASP's PhotoRates.in.

CMAQ publishes, for every label, absorption cross section CS and effective quantum yield EQY
averaged over 7 wavelength bins (CSQY_DATA_<mech>), plus a reference photon flux per bin.
A label's rate is J = sum_b a_b * F_b * CS_b * EQY_b where a_b is the (unknown) attenuation of bin b
for the case at hand.  a_b(t) is obtained, for every time row of the case's TUV table, by non-negative
least squares against the TUV columns that have an unambiguous CMAQ analogue (relative residuals).
"""
import re
import numpy as np
from scipy.optimize import nnls

MECHS = '/home/claude/usepa/cmaq/CCTM/src/MECHS/'


def read_csqy(path):
    flux, temps, tab = None, [], {}
    lines = open(path, errors='ignore').read().replace('\r', '').split('\n')
    i = 0
    flux = []
    inband = False
    for ln in lines:
        s = ln.strip()
        if s.startswith('N_INLINE_BAND'):
            inband = True
            continue
        if inband:
            m = re.match(r'^(\d+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\dEe+.-]+),', s)
            if m:
                flux.append(float(m.group(5)))
                continue
            if s.startswith('!...') or s.startswith('!PHOTAB'):
                inband = False
        f = s.split()
        if len(f) == 10 and f[1] in ('CS', 'EQY'):
            lab, q, T = f[0], f[1], float(f[2])
            tab.setdefault(lab, {}).setdefault(q, {})[T] = np.array([float(x) for x in f[3:]])
    return np.array(flux), tab


def csqy_at(tab, lab, T):
    """per-bin CS*EQY at temperature T (linear interpolation in T)"""
    cs, eq = tab[lab]['CS'], tab[lab]['EQY']
    Ts = sorted(cs)
    out = []
    for q in (cs, eq):
        T_ = min(max(T, Ts[0]), Ts[-1])
        lo = max(t for t in Ts if t <= T_)
        hi = min(t for t in Ts if t >= T_)
        w = 0.0 if hi == lo else (T_ - lo) / (hi - lo)
        out.append((1 - w) * q[lo] + w * q[hi])
    return out[0] * out[1]


def fit_weights(Jtuv, A, rel=True):
    """Jtuv: (nlabel,), A: (nlabel,7) -> a (7,) >=0"""
    if Jtuv.max() <= 0:
        return np.zeros(A.shape[1])
    w = 1.0 / np.maximum(Jtuv, 1e-3 * Jtuv.max()) if rel else np.ones_like(Jtuv)
    a, _ = nnls(A * w[:, None], Jtuv * w)
    return a
