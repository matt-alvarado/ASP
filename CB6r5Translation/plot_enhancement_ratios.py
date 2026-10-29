#!/usr/bin/env python3
"""Time plot of dO3/dCO (ASP EnhanceRatios.txt) and dOA/dCO (from Aerosol*.txt) for baseline vs CB6r5 runs.
usage: plot_enhancement_ratios.py <runs_root> <out.png>   (runs_root has base/<Case> and cb6/<Case>)"""
import glob, re, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root, out = sys.argv[1], sys.argv[2]
DECKS = '/home/claude/asp/InorganicInputDecks'
CASES = ['Williams', 'Otavi', 'Alaska', 'Timbavati']
TP = {'Williams': (288.4, 880.0), 'Otavi': (290.0, 730.0), 'Alaska': (280.6, 846.5), 'Timbavati': (288.0, 900.0)}


def table(path, strict=True):
    # Fortran overflow fields ('*****') are turned into NaN so the rows keep their column count
    L = [re.sub(r'\*+', ' nan ', l).split() for l in open(path) if l.strip()]
    h = L[0]
    return h, np.array([[float(x) for x in r] for r in L[1:] if (len(r) == len(h) or not strict)])


def co_bg(case):
    for ln in open(f'{DECKS}/{case}/GasPhaseChems.in', errors='ignore'):
        f = [x.strip() for x in ln.split('!')[0].split(';')]
        if f[0] == 'CO' and len(f) >= 7:
            return float(f[5])
    raise KeyError


def oa_series(run):
    tot = {}
    for f in glob.glob(f'{run}/OutputFiles/Aerosol[0-9]*.txt'):
        h, a = table(f, strict=False)
        for r in a:
            tot[r[0]] = tot.get(r[0], 0.0) + r[h.index('TotPOA')] + r[h.index('TotSOA')]
    t = np.array(sorted(tot))
    return t, np.array([tot[x] for x in t])


fig, ax = plt.subplots(2, 4, figsize=(15, 6.6), sharex='col', sharey='row')
sty = {'base': dict(color='#1f4e79', ls='-', label='Baseline mechanism'),
       'cb6': dict(color='#c4511a', ls='--', label='CB6r5')}
summary = []
for j, case in enumerate(CASES):
    T, P = TP[case]
    bg = co_bg(case)
    for tag in ('base', 'cb6'):
        run = f'{root}/{tag}/{case}'
        h, g = table(f'{run}/OutputFiles/GasOut1.txt')
        tg, co = g[:, 0], g[:, h.index('CO')]
        he, e = table(f'{run}/OutputFiles/EnhanceRatios.txt')
        te, o3r = e[:, 0], e[:, he.index('dO3/dCO')]
        keep = np.r_[True, np.diff(te) > 0]
        te, o3r = te[keep], o3r[keep]
        ax[0, j].plot(te / 3600, o3r, **sty[tag])
        ta, oa = oa_series(run)
        dco = np.interp(ta, tg, co) - bg                      # ppb
        oa_stp = oa * (1013.25 / P) * (T / 273.15)            # ug/m3 STP
        ratio = oa_stp / (dco / 1000.0)                       # ug/m3 STP per ppmv
        m = ta >= 600                                         # skip the first 10 min (dCO ~ 0 or OA ~ 0)
        ax[1, j].plot(ta[m] / 3600, ratio[m], **sty[tag])
        summary.append((case, tag, o3r[-1], ratio[-1], ta[-1] / 3600))
    ax[0, j].set_title(case)
    ax[1, j].set_xlabel('Time (h)')
    for r in (0, 1):
        ax[r, j].grid(alpha=0.3)
ax[0, 0].set_ylabel(r'$\Delta$O$_3$/$\Delta$CO (ppb/ppb)')
ax[1, 0].set_ylabel(r'$\Delta$OA/$\Delta$CO ($\mu$g m$^{-3}$ STP / ppmv)')
ax[0, 0].legend(frameon=False, fontsize=9)
fig.suptitle('Enhancement ratios: baseline vs CB6r5 gas mechanism (common_asp f3f1107)', y=0.995)
fig.tight_layout()
fig.savefig(out, dpi=150)
for s in summary:
    print('%-10s %-5s dO3/dCO %.3f  dOA/dCO %.1f ug/m3 STP/ppmv at %.1f h' % s)
