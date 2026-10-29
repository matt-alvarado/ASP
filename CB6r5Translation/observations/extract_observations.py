#!/usr/bin/env python3
"""Digitize the observed O3 / CO (and OA) enhancement data from the papers' vector figures.

Sources (attachments of this project):
  Williams : Alvarado et al. (2015) ACP 15, 6667, Fig. 5a  (dO3/dCO vs smoke age, with age and measurement error bars)
             and Fig. 9 / text: dOA/dCO2 = (2.83 +/- 1.08)e-3 g/g at 4.0-4.5 h
  Otavi    : Alvarado & Prinn (2009) JGR 114, D09306, Fig. 2a,b  (plume-average CO and O3 of Jost et al., 2003)
  Alaska   : same paper, Fig. 3a,b (Mason et al., 2006)
  Timbavati: same paper, Fig. 4a,c (AFTIR CO and O3, Hobbs et al., 2003) and Fig. 5 (OC 34 ug/m3 at 47 min)
Output: observed_enhancement_ratios.json  (per case: list of points with x = time (h), y = ratio, and ranges)
"""
import glob, json, os
import pdfplumber
from pdfdigitize import Panel

HERE = os.path.dirname(os.path.abspath(__file__))
D = glob.glob('/mnt/user-data/attachments/7be0c5a0*')[0]
JGR = pdfplumber.open(glob.glob(D + '/f0200068*')[0])
ACP = pdfplumber.open(glob.glob(D + '/c53f7496*')[0])
# deck backgrounds (ppb) used by ASP's own dX/dCO definition  (X - X_bg)/(CO - CO_bg)
BG = {'Williams': dict(CO=129.213, O3=0.0), 'Otavi': dict(CO=250.0, O3=85.0),
      'Alaska': dict(CO=100.0, O3=32.5), 'Timbavati': dict(CO=182.0, O3=35.7)}


def pair(o3, co, case, tmax, tol):
    bg = BG[case]
    out = []
    for p in o3:
        if p['x'] > tmax:
            continue
        q = min(co, key=lambda c: abs(c['x'] - p['x']))
        if abs(q['x'] - p['x']) > tol:
            continue
        dco = q['y'] - bg['CO']
        if dco <= 0:
            continue
        num = p['y'] - bg['O3']
        lo = (p['ylo'] - bg['O3']) / ((q['yhi'] if q['yhi'] > q['y'] else q['y']) - bg['CO'])
        hi = (p['yhi'] - bg['O3']) / ((q['ylo'] if q['ylo'] < q['y'] else q['y']) - bg['CO'])
        out.append(dict(t_h=p['x'] / 60, tlo_h=p['xlo'] / 60, thi_h=p['xhi'] / 60, y=num / dco,
                        ylo=min(lo, hi), yhi=max(lo, hi), O3=p['y'], CO=q['y']))
    return sorted(out, key=lambda r: r['t_h'])


obs = {}

# ---- Williams: ACP Fig. 5a, x = smoke age (h), y = dO3/dCO ---------------------------------------
P = Panel(ACP.pages[10], (100.4, 290.6, 87.6, 228.1), ('x', 0, 5), ('y', 0.20, -0.05))
pts = []
for r in P.points(2.3):
    pts.append(dict(t_h=r['x'], tlo_h=r['xlo'], thi_h=r['xhi'], y=r['y'], ylo=r['ylo'], yhi=r['yhi']))
obs['Williams'] = dict(O3=sorted(pts, key=lambda r: r['t_h']), OA=[])
# OA: dOA/dCO2 = 2.83 +/- 1.08 e-3 g/g, 4.0-4.5 h  ->  ug m-3 (STP) per ppmv CO using the inert-tracer ratio of the deck
#   dCO2/dCO = (528440.21 - 390000) / (10000 - 129.213) mol/mol  (CO2, CO both conserved apart from dilution)
r_co2_co = (528440.21 - 390000.0) / (10000.0 - 129.213)
ug_co2_per_ppmv_co = r_co2_co * 44.01 / 22.414 * 1000.0       # ug/m3 (STP: 273.15 K, 1 atm) of CO2 per ppmv CO
obs['Williams']['OA'].append(dict(t_h=4.25, tlo_h=4.0, thi_h=4.5, y=2.83e-3 * ug_co2_per_ppmv_co,
                                  ylo=(2.83e-3 - 1.08e-3) * ug_co2_per_ppmv_co, yhi=(2.83e-3 + 1.08e-3) * ug_co2_per_ppmv_co,
                                  note='dOA/dCO2 = (2.83+/-1.08)e-3 g/g (ACP 2015 Fig. 9) converted with dCO2/dCO = %.2f mol/mol' % r_co2_co))

# ---- Otavi: JGR Fig. 2 (page rotated: time runs up the raw y axis, 0-150 min; concentration along raw x) ---
pg = JGR.pages[3]
o3 = Panel(pg, (131.7, 250.8, 312.6, 460.3), ('y', 150, 0), ('x', 110, 50)).points(2.8)
co = Panel(pg, (132.8, 251.8, 498.9, 646.6), ('y', 150, 0), ('x', 2000, 0)).points(2.8)
o3 = [p for p in o3 if p['y'] > 80]                  # drop the legend marker (69 ppb at 28 min)
co = [p for p in co if p['y'] < 1000]                # drop the legend marker (1367 ppb)
obs['Otavi'] = dict(O3=pair(o3, co, 'Otavi', 150, 1.5), OA=[])

# ---- Alaska: JGR Fig. 3 (0-180 min) -----------------------------------------------------------------
pg = JGR.pages[5]
o3 = Panel(pg, (67.3, 180.6, 221.8, 362.4), ('y', 180, 0), ('x', 200, 0)).points(2.7)
co = Panel(pg, (66.6, 179.9, 430.3, 570.8), ('y', 180, 0), ('x', 5000, 0)).points(2.7)
o3 = [p for p in o3 if p['y'] < 170]                 # legend marker at 180 ppb
co = [p for p in co if p['x'] > 10 and p['y'] < 4000]  # legend marker (4549) and the t=0 initial point (4498)
obs['Alaska'] = dict(O3=pair(o3, co, 'Alaska', 180, 3.0), OA=[])

# ---- Timbavati: JGR Fig. 4 (0-60 min), AFTIR circles ------------------------------------------------
pg = JGR.pages[7]
o3 = Panel(pg, (236.1, 403.7, 240.6, 376.2), ('x', 0, 60), ('y', 160, 0)).points(3.3)
co = Panel(pg, (122.6, 282.7, 73.8, 202.9), ('x', 0, 60), ('y', 7000, 0)).points(3.1)
o3 = [p for p in o3 if p['y'] < 140 or p['x'] > 10]   # legend marker at (5.6 min, 148 ppb)
o3 = [p for p in o3 if not (p['x'] < 6 and p['y'] > 100)]
co = [p for p in co if not (p['x'] > 20 and p['x'] < 21 and p['y'] > 5000)]   # legend marker
obs['Timbavati'] = dict(O3=pair(o3, co, 'Timbavati', 60, 1.5), OA=[])
# OA: Fig. 5 / text: observed OC 34 ug/m3 at 47 min downwind; no error bar shown.  Convert with the AFTIR CO at 46.8 min.
co47 = min(co, key=lambda c: abs(c['x'] - 46.8))
T, Pm = 288.0, 900.0
oa_stp = 34.0 * (1013.25 / Pm) * (T / 273.15)
obs['Timbavati']['OA'].append(dict(t_h=47 / 60, tlo_h=47 / 60, thi_h=47 / 60,
                                   y=oa_stp / ((co47['y'] - BG['Timbavati']['CO']) / 1000.0),
                                   ylo=None, yhi=None,
                                   note='OC 34 ug/m3 (JGR Fig. 5, ambient) -> STP; dCO from the AFTIR point at %.1f min (%.0f ppb)' % (co47['x'], co47['y'])))
json.dump(obs, open(os.path.join(HERE, 'observed_raw.json'), 'w'), indent=1)
print({k: (len(v['O3']) if 'O3' in v else None) for k, v in obs.items()})
