#!/usr/bin/env python3
"""dO3/dCO, dOA/dCO and dNOx/dCO vs time for the NEIVA MCE-dependent decks only.
usage: plot_mce_enhancement.py <runs_root> <out.png> [fire-type tag, default temperate_forest]"""
import glob, re, sys, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CB6r5Translation'))
root, out = sys.argv[1], sys.argv[2]; tag = sys.argv[3] if len(sys.argv) > 3 else 'temperate_forest'
T, P, BG = 288.4, 880.0, 129.213          # Williams deck


def table(path, strict=True):
    L = [re.sub(r'\*+', ' nan ', l).split() for l in open(path) if l.strip()]
    h = L[0]
    return h, np.array([[float(x) for x in r] for r in L[1:] if (len(r) == len(h) or not strict)])


def oa_series(run):
    """Time, TotPOA+TotSOA summed over the Aerosol*.txt bin files. Rows are read by fixed width (13 characters per field) because
    three-digit exponents (E-100) fuse neighbouring fields in later columns; Time/TotPOA/TotSOA come early in the row."""
    tot = {}
    for f in glob.glob(f'{run}/OutputFiles/Aerosol[0-9]*.txt'):
        L = [l.rstrip('\n') for l in open(f) if l.strip()]
        hd = L[0].split(); ip, isoa = hd.index('TotPOA'), hd.index('TotSOA')
        for l in L[1:]:
            try:
                t = float(l[0:13]); v = float(l[ip * 13:(ip + 1) * 13]) + float(l[isoa * 13:(isoa + 1) * 13])
            except ValueError:
                continue
            tot[t] = tot.get(t, 0.0) + v
    t = np.array(sorted(tot)); return t, np.array([tot[x] for x in t])


runs = sorted(glob.glob(f'{root}/NEIVA_{tag}_MCE*'))
mces = [int(re.search(r'MCE(\d+)', r).group(1)) / 1000 for r in runs]
norm = colors.Normalize(min(mces), max(mces)); cmap = cm.viridis
fig, AX = plt.subplots(2, 3, figsize=(15, 9), sharex=True); ax = AX[0]; bx = AX[1]
rows = []
for run, m in zip(runs, mces):
    c = cmap(norm(m))
    h, g = table(f'{run}/OutputFiles/GasOut1.txt'); tg, co = g[:, 0], g[:, h.index('CO')]
    he, e = table(f'{run}/OutputFiles/EnhanceRatios.txt'); te = e[:, 0]; keep = np.r_[True, np.diff(te) > 0]; te = te[keep]
    o3 = e[:, he.index('dO3/dCO')][keep]; nox = e[:, he.index('dNOx/dCO')][keep]
    ta, oa = oa_series(run); dco = np.interp(ta, tg, co) - BG
    ratio = oa * (1013.25 / P) * (T / 273.15) / (dco / 1000.0); k = ta >= 600
    ax[0].plot(te / 3600, o3, color=c, lw=1.8); ax[1].plot(ta[k] / 3600, ratio[k], color=c, lw=1.8); ax[2].plot(te / 3600, nox, color=c, lw=1.8)
    dcoe = np.interp(te, tg, co) - BG
    bx[0].plot(te / 3600, o3 * dcoe, color=c, lw=1.8); bx[1].plot(ta / 3600, oa, color=c, lw=1.8); bx[2].plot(te / 3600, nox * dcoe, color=c, lw=1.8)
    rows.append((m, o3[-1], ratio[-1], nox[-1], nox[np.searchsorted(te, 600)] if len(te) > 3 else np.nan))
for a, yl in zip(ax, [r'$\Delta$O$_3$/$\Delta$CO (ppb/ppb)', r'$\Delta$OA/$\Delta$CO ($\mu$g m$^{-3}$ STP / ppmv)', r'$\Delta$NO$_x$/$\Delta$CO (ppb/ppb)']):
    a.set_xlabel('Time (h)'); a.set_ylabel(yl); a.grid(alpha=0.3); a.set_xlim(0, 5.05)
for a, yl in zip(bx, [r'$\Delta$O$_3$ (ppb)', r'OA ($\mu$g m$^{-3}$, ambient)', r'$\Delta$NO$_x$ (ppb)']):
    a.set_xlabel('Time (h)'); a.set_ylabel(yl); a.grid(alpha=0.3)
ax[2].set_yscale('log'); bx[2].set_yscale('log'); bx[1].set_yscale('log')
fig.colorbar(cm.ScalarMappable(norm, cmap), ax=list(AX.ravel()), label='MCE', pad=0.01)
fig.suptitle('NEIVAv1.0 %s emission factors, MCE %.3f-%.3f, Williams baseline (RACM2-lumped) mechanism' % (tag.replace('_', ' '), min(mces), max(mces)))
fig.savefig(out, dpi=150, bbox_inches='tight')
print('MCE   dO3/dCO(5h)  dOA/dCO(5h)  dNOx/dCO(5h)')
for r in rows: print('%.3f  %.3f  %.1f  %.4f' % r[:4])
