#!/usr/bin/env python3
"""Time plot of dO3/dCO (ASP EnhanceRatios.txt) and dOA/dCO (from Aerosol*.txt) for all gas mechanisms, with observations.
usage: plot_enhancement_ratios.py <runs_root> <out.png>
runs_root holds one run directory per deck, named like the deck: <Case>, <Case>_CB6r5, <Case>_CRACMM2, <Case>_CRACMM3,
<Case>_SAPRC07TC, <Case>_SAPRC07TIC  (each with OutputFiles/)."""
import glob, json, re, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root, out = sys.argv[1], sys.argv[2]
OBS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'observations', 'observed_raw.json')))
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


fig, ax = plt.subplots(3, 4, figsize=(16, 10.2), sharex='col', sharey='row')
# tag -> deck-name suffix, style
MECHS = [('base', '', dict(color='#0072B2', ls='-', lw=1.8, label='Baseline (RACM2-based)'), 'o'),
         ('cb6', '_CB6r5', dict(color='#D55E00', ls='--', lw=1.6, label='CB6r5'), 's'),
         ('cr2', '_CRACMM2', dict(color='#009E73', ls='-.', lw=1.6, label='CRACMM2'), '^'),
         ('cr3', '_CRACMM3', dict(color='#CC79A7', ls=(0, (5, 1, 1, 1, 1, 1)), lw=1.6, label='CRACMM3'), 'D'),
         ('stc', '_SAPRC07TC', dict(color='#222222', ls=':', lw=2.0, label='SAPRC07TC'), 'v'),
         ('sti', '_SAPRC07TIC', dict(color='#E0B000', ls=(0, (3, 1, 1, 1)), lw=1.8, label='SAPRC07TIC'), 'P'),
         ('mcm', '_MCM331', dict(color='#56B4E9', ls=(0, (1, 1)), lw=2.2, label='MCM v3.3.1 + lumped'), 'X'),
         ('mz1', '_MOZT1', dict(color='#7B3FA0', ls=(0, (4, 1.5)), lw=1.8, label='MOZART-T1'), '<'),
         ('gc', '_GCHEM', dict(color='#8C564B', ls=(0, (6, 1.5, 1, 1.5)), lw=1.8, label='GEOS-Chem 14.8'), 'h')]
MARK = {t: mk for t, _, _, mk in MECHS}
sty = {t: st for t, _, st, _ in MECHS}
summary = []
for j, case in enumerate(CASES):
    T, P = TP[case]
    bg = co_bg(case)
    for tag, suf, _, _ in MECHS:
        run = f'{root}/{case}{suf}'
        if not os.path.isdir(run + '/OutputFiles'):
            continue
        h, g = table(f'{run}/OutputFiles/GasOut1.txt')
        tg, co = g[:, 0], g[:, h.index('CO')]
        he, e = table(f'{run}/OutputFiles/EnhanceRatios.txt')
        te, o3r = e[:, 0], e[:, he.index('dO3/dCO')]
        keep = np.r_[True, np.diff(te) > 0]
        te, o3r = te[keep], o3r[keep]
        ax[0, j].plot(te / 3600, o3r, **sty[tag])
        xm = 5.0 - 0.14 * [t for t, _, _, _ in MECHS].index(tag)     # staggered so overlapping lines keep visible markers
        ax[0, j].plot(xm, np.interp(xm, te / 3600, o3r), MARK[tag], color=sty[tag]['color'], ms=9, mec='white', mew=0.8, zorder=6)
        nox = e[:, he.index('dNOx/dCO')][keep]
        ax[2, j].plot(te / 3600, nox, **sty[tag])
        ax[2, j].plot(xm, np.interp(xm, te / 3600, nox), MARK[tag], color=sty[tag]['color'], ms=9, mec='white', mew=0.8, zorder=6)
        ta, oa = oa_series(run)
        dco = np.interp(ta, tg, co) - bg                      # ppb
        oa_stp = oa * (1013.25 / P) * (T / 273.15)            # ug/m3 STP
        ratio = oa_stp / (dco / 1000.0)                       # ug/m3 STP per ppmv
        m = ta >= 600                                         # skip the first 10 min (dCO ~ 0 or OA ~ 0)
        ax[1, j].plot(ta[m] / 3600, ratio[m], **sty[tag])
        ax[1, j].plot(xm, np.interp(xm, ta / 3600, ratio), MARK[tag], color=sty[tag]['color'], ms=9, mec='white', mew=0.8, zorder=6)
        summary.append((case, tag, o3r[-1], ratio[-1], nox[-1], ta[-1] / 3600))
    def obs_plot(a, pts, first):
        for p in pts:
            xe = [[max(p['t_h'] - p['tlo_h'], 0)], [max(p['thi_h'] - p['t_h'], 0)]]
            ye = None
            if p.get('ylo') is not None:
                ye = [[max(p['y'] - p['ylo'], 0)], [max(p['yhi'] - p['y'], 0)]]
            a.errorbar([p['t_h']], [p['y']], xerr=xe, yerr=ye, fmt='o', ms=4, color='k', mfc='white', ecolor='0.35',
                       elinewidth=0.9, capsize=1.5, zorder=5, label='Observed (digitized from paper)' if first else None)
            first = False
    obs_plot(ax[0, j], OBS[case]['O3'], j == 0)
    obs_plot(ax[1, j], OBS[case]['OA'], False)
    obs_plot(ax[2, j], OBS[case].get('NOx', []), False)
    ax[0, j].set_title(case)
    ax[0, j].set_xlim(-0.15, 5.3)
    ax[2, j].set_xlabel('Time (h)')
    for r in (0, 1, 2):
        ax[r, j].grid(alpha=0.3)
ax[0, 0].set_ylabel(r'$\Delta$O$_3$/$\Delta$CO (ppb/ppb)')
ax[1, 0].set_ylabel(r'$\Delta$OA/$\Delta$CO ($\mu$g m$^{-3}$ STP / ppmv)')
ax[2, 0].set_ylabel(r'$\Delta$NO$_x$/$\Delta$CO (ppb/ppb)')
from matplotlib.lines import Line2D
h_ = [Line2D([], [], marker=MARK[t], mec='white', ms=8, **st) for t, _, st, _ in MECHS] + ax[0, 0].get_legend_handles_labels()[0][-1:]
l_ = [st['label'] for _, _, st, _ in MECHS] + ax[0, 0].get_legend_handles_labels()[1][-1:]
fig.legend(h_, l_, loc='lower center', ncol=5, frameon=False, fontsize=9)
fig.suptitle('Enhancement ratios, 5 h runs: nine gas-phase mechanisms (observations from the papers where available; none for NO$_x$ at Otavi/Timbavati)', y=0.995)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(out, dpi=150)
for s in summary:
    print('%-10s %-5s dO3/dCO %.3f  dOA/dCO %.1f ug/m3 STP/ppmv  dNOx/dCO %.4f at %.1f h' % s)
