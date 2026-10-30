"""Build MCE-dependent ASP input decks from NEIVAv1.0 emission factors for one fire type, on top of the Williams baseline deck.

    python3 build_neiva_decks.py [--type "temperate forest"] [--mce 0.800 0.975 0.025] [--template Williams] [--no-unid-scale]

Per deck (MCE m):  dCO = template CO excess (10000 - 129.2 ppb), dC = dCO/(1-m) (ppb of emitted carbon as CO+CO2),
  CO2 excess = dCO m/(1-m);  every other species:  ppb = y_X(m) dC 12.011 / M_X   (y in g X per g C, see mce_fit.py)
Routes (neiva_map.py): explicit species | baseline lumps (carbon-conserving) | IVOC bins by C* (mass-conserving) | skipped.
Aerosol: OA (particle IVOC3-6, template volatility split), BC and the salts NH4Cl/NH4NO3/KCl/K2SO4 from the NEIVA ions; Nt scaled with mass.
Unidentified NMOG pool (template IVOC3-7 gas, Akagi-based) = r * identified NEIVA NMOG(m), r fixed from the template at its own MCE.
Species with no NEIVA data (H2, CH4, PAN, H2O2, HNO3, HCN, ...) keep the template value.
"""
import argparse, os, re, shutil, sys, csv
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import coverage as cv, mce_fit as mf, neiva_map as nm

DECKS = os.path.join(HERE, '..', 'InorganicInputDecks')
F = 0.0367 * 12.011          # ug C / m3 per ppb of carbon (ambient 288.4 K, 880 mb)
INORG = {'NO': 'NO', 'NO2': 'NO2', 'NH3': 'NH3', 'HNO2': 'HONO', 'HCl': 'HCl', 'SO2': 'SO2'}
SALT_MW = dict(NH4Cl=53.49, NH4NO3=80.04, KCl=74.55, K2SO4=174.26)
ION_MW = dict(sulfate=96.06, nitrate=62.0, chloride=35.45, ammonium=18.04, potassium=39.10)


def parse_gas(path):
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')
    rows = {}
    for i, l in enumerate(L):
        if l.startswith('!') or ';' not in l: continue
        f = l.split('!')[0].split(';')
        if len(f) >= 7:
            try: rows[f[0].strip()] = (i, f); float(f[4]); float(f[5])
            except ValueError: rows.pop(f[0].strip(), None)
    return L, rows


def setrow(L, rows, name, init):
    i, f = rows[name]
    L[i] = '%s\t;%s;%s;%s;%.6g;%s;%s!' % (f[0].strip(), f[1].strip(), f[2].strip(), f[3].strip(), init, f[5].strip(), f[6].strip())


def build(args):
    d, fit, ef, cols, mce, Cm, co, co2 = cv.classified(args.type)
    tdir = os.path.join(DECKS, args.template)
    L0, rows = parse_gas(os.path.join(tdir, 'GasPhaseChems.in'))
    g = lambda n: (float(rows[n][1][4]), float(rows[n][1][5]))
    dCO = g('CO')[0] - g('CO')[1]
    dCO2_t = g('CO2')[0] - g('CO2')[1]
    mce_ref = dCO2_t / (dCO2_t + dCO)
    ivoc_t = {k: g('IVOC%d' % k)[0] for k in range(3, 8)}
    pool_t = {k: ivoc_t[k] * nm.IVOC_MW[k] / F * 1.0 for k in ivoc_t}      # in "g per ppbC" bookkeeping: ppb * MW (mass, /0.0367 later)
    pool_mass_t = sum(ivoc_t[k] * nm.IVOC_MW[k] for k in ivoc_t) * 0.0367  # ug/m3
    allf = pd.concat([fit], ignore_index=True)

    def yv(r, m): return mf.evaluate(r, m)
    nmoc = d[d.route != 'skip']
    ident = lambda m: sum(yv(r, m) for _, r in nmoc.iterrows())
    dC_ref = dCO / (1 - mce_ref)
    r_unid = (pool_mass_t / (dC_ref * F)) / ident(mce_ref) if args.unid_scale else None
    # template aerosol
    A = open(os.path.join(tdir, 'AerosolModes.in'), errors='ignore').read().replace('\r', '')
    tot_t = float(sum(map(float, re.findall(r'^(?:NH4Cl|NH4NO3|KCl|K2SO4|BC|IVOC\d); *([0-9.]+)', A, re.M))))
    oa_split = {k: float(re.search(r'^IVOC%d; *([0-9.]+)' % k, A, re.M).group(1)) for k in range(3, 7)}
    s = sum(oa_split.values()); oa_split = {k: v / s for k, v in oa_split.items()}
    pool_split = {k: ivoc_t[k] * nm.IVOC_MW[k] for k in ivoc_t}; s = sum(pool_split.values()); pool_split = {k: v / s for k, v in pool_split.items()}
    bycmp = lambda c: fit[(fit.compound == c)].iloc[0]
    info = dict(mce_ref=mce_ref, r_unid=r_unid, dCO=dCO, pool_mass_t=pool_mass_t, n_study_cols=len(cols),
                mce_data=(float(mce.min()), float(mce.max())))
    out = []
    for m in np.arange(args.mce[0], args.mce[1] + 1e-9, args.mce[2]):
        m = round(float(m), 4)
        dC = dCO / (1 - m)
        L = list(L0); rr = dict(rows); mapping = []
        touched = {}
        def put(name, ppb, why):
            touched[name] = touched.get(name, 0.0) + ppb; mapping.append((name, why, ppb))
        # CO / CO2 from carbon balance
        L = L0.copy()
        touched['CO'] = dCO + g('CO')[1]; touched['CO2'] = dCO * m / (1 - m) + g('CO2')[1]
        ivb = {k: 0.0 for k in range(1, 10)}
        lump = {k: 0.0 for k in nm.LUMPS}; expl = {}; unid_ident = 0.0
        for _, r in d.iterrows():
            if r.route == 'skip': continue
            y = yv(r, m)
            if y <= 0: continue
            if r.route == 'explicit': expl[r.target] = expl.get(r.target, 0.0) + y * dC * 12.011 / r.mm
            elif r.route == 'lump':
                lump[r.target] += y * dC * 12.011 / r.mm * r.nC / nm.LUMPS[r.target]
            else: ivb[int(r.target)] += y * dC * 12.011 / nm.IVOC_MW[int(r.target)]
        for k, v in expl.items(): touched[k] = v + (g(k)[1] if k in rr else 0.0)
        for k, v in lump.items(): touched[k] = v + (g(k)[1] if k in rr else 0.0)
        # unidentified NMOG pool
        if r_unid is not None:
            pm = r_unid * ident(m) * dC * F      # ug/m3
            for k, fr in pool_split.items(): ivb[k] += pm * fr / (nm.IVOC_MW[k] * 0.0367)
        else:
            for k in ivoc_t: ivb[k] += ivoc_t[k] * dC / dC_ref
        for k, v in ivb.items(): touched['IVOC%d' % k] = v
        touched['XYP'] = 0.0          # template xylene lump replaced by NEIVA members (XYM / XYO)
        # inorganic gases
        for fm, sp in INORG.items():
            q = fit[(fit.formula == fm) & (fit.category == 'inorganic gas')]
            if len(q) and q.iloc[0].n:
                rrow = q.iloc[0]; M = float(ef.at[int(rrow.row), 'mm'])
                touched[sp] = yv(rrow, m) * dC * 12.011 / M + g(sp)[1]
        for n, v in touched.items():
            if n not in rr: raise SystemExit('template has no species ' + n)
            setrow(L, rr, n, v)
        # aerosol
        def yc(c, lg=None): return yv(bycmp(c), m)
        mass = lambda y: y * dC * F
        oa = mass(yc('Organic Aerosol (OA)')); bc = mass(yc('BC'))
        n_ = {k: mass(yc(c)) / ION_MW[k] for k, c in dict(sulfate='sulfate', nitrate='nitrate', chloride='chloride', ammonium='ammonium', potassium='potassium').items()}
        k2so4 = min(n_['sulfate'], n_['potassium'] / 2); kl = n_['potassium'] - 2 * k2so4
        kcl = min(kl, n_['chloride']); nh4no3 = min(n_['nitrate'], n_['ammonium']); nh4cl = min(n_['chloride'] - kcl, n_['ammonium'] - nh4no3)
        salts = dict(NH4Cl=nh4cl * SALT_MW['NH4Cl'], NH4NO3=nh4no3 * SALT_MW['NH4NO3'], KCl=kcl * SALT_MW['KCl'], K2SO4=k2so4 * SALT_MW['K2SO4'])
        modes = dict(salts, BC=bc); modes.update({'IVOC%d' % k: oa * f for k, f in oa_split.items()})
        tot = sum(modes.values())
        AA = A
        for k, v in modes.items(): AA = re.sub(r'^(%s); *[0-9.eE+-]+;' % k, lambda mo: '%s; %.5g;' % (mo.group(1), v), AA, flags=re.M)
        AA = re.sub(r'^234000 ;', '%d ;' % round(234000 * tot / tot_t), AA, count=1, flags=re.M)
        name = '%s_MCE%04d_baseline' % (args.type.replace(' ', '_').replace('temperate_forest', 'temperate_forest'), round(m * 1000))
        name = 'NEIVA_%s' % name if not name.startswith('NEIVA') else name
        nd_ = os.path.join(DECKS, name)
        if os.path.isdir(nd_): shutil.rmtree(nd_)
        shutil.copytree(tdir, nd_)
        L[0] = '!! GasPhaseChems.in - NEIVAv1.0 %s, MCE = %.3f, on the %s baseline deck' % (args.type, m, args.template)
        open(os.path.join(nd_, 'GasPhaseChems.in'), 'w').write('\n'.join(L) + '\n')
        open(os.path.join(nd_, 'AerosolModes.in'), 'w').write(AA)
        with open(os.path.join(nd_, 'NEIVA_SPECIES_MAPPING.csv'), 'w', newline='') as fh:
            w = csv.writer(fh); w.writerow(['deck_species', 'route', 'ppb_excess_or_value'])
            for a, b, c in mapping: w.writerow([a, b, c])
            for a in sorted(touched): w.writerow([a, 'final_initial_ppb', touched[a]])
        open(os.path.join(nd_, 'NEIVA_PROVENANCE.txt'), 'w').write(
            'NEIVAv1.0 (Shahid et al., GMD 17, 7679, 2024), fire type: %s, MCE = %.3f\nTemplate deck: %s (template MCE from its CO2/CO: %.3f)\n'
            'Study columns used: %d (MCE %.3f-%.3f)%s\nCO = %.1f ppb excess; carbon excess dC = %.0f ppb; CO2 excess = %.0f ppb (carbon balance)\n'
            'Unidentified-NMOG ratio r = %s\nAerosol (ug/m3): OA %.1f, BC %.2f, %s\n' % (
                args.type, m, args.template, mce_ref, len(cols), mce.min(), mce.max(),
                '  ** extrapolated beyond the data **' if m > mce.max() or m < mce.min() else '', dCO, dC, dCO * m / (1 - m),
                'off' if r_unid is None else '%.3f' % r_unid, oa, bc, ', '.join('%s %.1f' % kv for kv in salts.items())))
        out.append(dict(deck=name, mce=m, dC=dC, oa=oa, bc=bc, nox=touched['NO'] + touched['NO2'], unid=sum(ivb[k] for k in range(3, 8)), Nt=round(234000 * tot / tot_t)))
    return pd.DataFrame(out), info


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--type', default='temperate forest'); ap.add_argument('--template', default='Williams')
    ap.add_argument('--mce', nargs=3, type=float, default=[0.80, 0.975, 0.025]); ap.add_argument('--no-unid-scale', dest='unid_scale', action='store_false')
    a = ap.parse_args()
    df, info = build(a)
    print(info); print(df.to_string())
