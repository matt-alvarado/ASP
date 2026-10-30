#!/usr/bin/env python3
"""dOA/dCO and dO3/dCO of the CRACMM2/3 decks before (baseline IVOC ladder + counters) and after (CMAQ organic-aerosol chemistry).
usage: compare_cmaq_oa.py <runs_root> <out.png>   (runs_root holds <Case>_CRACMMn = new and old_<Case>_CRACMMn = old runs)"""
import glob, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
root, out = sys.argv[1], sys.argv[2]
MECHS = sys.argv[3].split(',') if len(sys.argv) > 3 else ['CRACMM2', 'CRACMM3']
COL = {'CRACMM2': '#009E73', 'CRACMM3': '#CC79A7', 'CB6r5': '#0072B2', 'SAPRC07TC': '#E69F00', 'SAPRC07TIC': '#D55E00'}
OAKIND = 'CMAQ CRACMM' if MECHS[0].startswith('CRACMM') else 'CMAQ AERO6/7'
OLDKIND = 'baseline IVOC-ladder organic aerosol'
usage = 'usage: compare_cmaq_oa.py <runs_root> <out.png> [MECH,MECH,...]'

HERE = os.path.dirname(os.path.abspath(__file__)); OBS = json.load(open(os.path.join(HERE, 'observations', 'observed_raw.json')))
CASES = ['Williams', 'Otavi', 'Alaska', 'Timbavati']
TP = {'Williams': (288.4, 880.0), 'Otavi': (290.0, 730.0), 'Alaska': (280.6, 846.5), 'Timbavati': (288.0, 900.0)}


def table(path):
    L = [re.sub(r'\*+', ' nan ', l).split() for l in open(path) if l.strip()]
    return L[0], np.array([[float(x) for x in r] for r in L[1:] if len(r) == len(L[0])])


def oa_series(run):
    """total organic aerosol (ug/m3, all bins) from SizeDist*min.txt TotOA (sums every organic species except BC; the Aerosol*.txt TotPOA/TotSOA
    columns hard-code the baseline organic-species ordering and are not meaningful for the CMAQ-OA decks)"""
    tot = {}
    for f in glob.glob(f'{run}/OutputFiles/SizeDist*min.txt'):
        L = [l.rstrip('\n') for l in open(f) if l.strip()]
        i = L[0].index('TotOA') // 13
        for l in L[1:]:
            try: t = 60.0 * float(l[0:13]); v = float(l[i * 13:(i + 1) * 13])
            except ValueError: continue
            tot[t] = tot.get(t, 0.0) + v
    t = np.array(sorted(tot)); return t, np.array([tot[x] for x in t])


def co_bg(case):
    for ln in open(f'/home/claude/asp/InorganicInputDecks/{case}/GasPhaseChems.in', errors='ignore'):
        f = [x.strip() for x in ln.split('!')[0].split(';')]
        if f[0] == 'CO' and len(f) >= 7: return float(f[5])


STY = {}
for _m, _c in COL.items():
    STY[(_m, 'old')] = dict(color=_c, ls='--', lw=1.3)
    STY[(_m, 'new')] = dict(color=_c, ls='-', lw=2.0)
fig, ax = plt.subplots(2, 4, figsize=(16, 7.2), sharex='col', sharey='row')
for j, case in enumerate(CASES):
    T, P = TP[case]; bg = co_bg(case)
    for mech in MECHS:
        for ver, pre in (('old', 'old_'), ('new', '')):
            run = f'{root}/{pre}{case}_{mech}'
            if not os.path.isdir(run + '/OutputFiles'): continue
            h, g = table(f'{run}/OutputFiles/GasOut1.txt'); tg, co = g[:, 0], g[:, h.index('CO')]
            he, e = table(f'{run}/OutputFiles/EnhanceRatios.txt'); te = e[:, 0]; k = np.r_[True, np.diff(te) > 0]
            ax[1, j].plot(te[k] / 3600, e[:, he.index('dO3/dCO')][k], label=f'{mech} {ver}', **STY[(mech, ver)])
            ta, oa = oa_series(run); dco = np.interp(ta, tg, co) - bg
            r = oa * (1013.25 / P) * (T / 273.15) / (dco / 1000.0); m = ta >= 600
            ax[0, j].plot(ta[m] / 3600, r[m], **STY[(mech, ver)])
    for row, key in ((1, 'O3'), (0, 'OA')):
        for p in OBS[case][key]:
            ye = None if p.get('ylo') is None else [[max(p['y'] - p['ylo'], 0)], [max(p['yhi'] - p['y'], 0)]]
            ax[row, j].errorbar([p['t_h']], [p['y']], xerr=[[max(p['t_h'] - p['tlo_h'], 0)], [max(p['thi_h'] - p['t_h'], 0)]], yerr=ye, fmt='o', ms=4,
                                color='k', mfc='white', ecolor='0.35', elinewidth=0.9, capsize=1.5, zorder=5, label='Observed' if (j == 0 and row == 0 and p is OBS[case][key][0]) else None)
    ax[0, j].set_title(case); ax[1, j].set_xlabel('Time (h)'); ax[0, j].set_xlim(-0.1, 5.1)
    for r_ in (0, 1): ax[r_, j].grid(alpha=0.3)
ax[0, 0].set_ylabel(r'$\Delta$OA/$\Delta$CO ($\mu$g m$^{-3}$ STP / ppmv)'); ax[1, 0].set_ylabel(r'$\Delta$O$_3$/$\Delta$CO (ppb/ppb)')
h_, l_ = ax[1, 0].get_legend_handles_labels(); h0, l0 = ax[0, 0].get_legend_handles_labels()
fig.legend(h_ + h0[-1:], l_ + l0[-1:], loc='lower center', ncol=len(MECHS) * 2 + 1, frameon=False)
fig.suptitle('%s: %s (dashed) vs %s organic-aerosol chemistry (solid)' % ('/'.join(MECHS), OLDKIND, OAKIND))
fig.tight_layout(rect=(0, 0.06, 1, 0.97)); fig.savefig(out, dpi=150)
print('case mech   dOA/dCO(5h) old new   dO3/dCO(5h) old new')
for case in CASES:
    T, P = TP[case]; bg = co_bg(case); row = []
    for mech in MECHS:
        vals = []
        for pre in ('old_', ''):
            run = f'{root}/{pre}{case}_{mech}'
            h, g = table(f'{run}/OutputFiles/GasOut1.txt'); tg, co = g[:, 0], g[:, h.index('CO')]
            he, e = table(f'{run}/OutputFiles/EnhanceRatios.txt'); ta, oa = oa_series(run)
            dco = np.interp(ta[-1], tg, co) - bg; vals.append((oa[-1] * (1013.25 / P) * (T / 273.15) / (dco / 1000.0), e[-1, he.index('dO3/dCO')]))
        print('%-10s %s  OA %.1f -> %.1f   O3 %.3f -> %.3f' % (case, mech, vals[0][0], vals[1][0], vals[0][1], vals[1][1]))
