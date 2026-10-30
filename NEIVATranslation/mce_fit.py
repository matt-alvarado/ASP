"""MCE model for one NEIVA fire type.

Per study column j (with MCE_j):   Cm_j = EF_CO2 12/44.01 + EF_CO 12/28.01   (g C emitted as CO + CO2 per kg fuel)
Per compound X:                    y_Xj = EF_Xj / Cm_j                       (g X per g C), fit versus MCE:
    'fit'   n >= NMIN, MCE span >= SPANMIN and slope p < PMAX :  y(MCE) = a + b MCE  (log-linear for OA, BC), clipped at zero
    'mean'  otherwise: y = mean of y_Xj (the fire type's average at its mean MCE), independent of MCE
CO and CO2 are not fitted: they follow from the carbon balance in the deck builder.
Lumped NEIVA rows (a measured sum of several compounds, id without an InChI) are apportioned to the member compounds of the same
formula, in the columns where no member is reported, with the members' mean fractions (cf. NEIVA's fractional-contribution step).
"""
import numpy as np, pandas as pd
from scipy import stats
import neiva_data as nd

NMIN, SPANMIN, PMAX = 10, 0.05, 0.10
LOGFIT = ('Organic Aerosol (OA)', 'BC')


def prepare(fire_type):
    ef, info = nd.efs(fire_type)
    cols = list(info.efcol); mce = info.set_index('efcol').MCE
    ef = ef.reset_index(drop=True)
    isinchi = ef.id.astype(str).str.startswith('InChI')
    nm = ef.pollutant_category == 'NMOC_g'
    # ---- apportion lumped rows
    ef['apportioned_from'] = ''
    for f, grp in ef[nm].groupby('formula'):
        mem = grp[isinchi[grp.index]]
        lump = grp[~isinchi[grp.index]]
        if mem.empty or lump.empty:
            continue
        mmean = mem[cols].mean(axis=1, skipna=True).fillna(0.0)
        if mmean.sum() <= 0:
            continue
        frac = mmean / mmean.sum()
        for li, lr in lump.iterrows():
            memb_any = mem[cols].notna().any(axis=0)          # columns in which a member is reported
            for c in cols:
                if pd.notna(lr[c]) and lr[c] > 0 and not memb_any[c]:
                    for mi in mem.index:
                        cur = ef.at[mi, c]
                        ef.at[mi, c] = (0.0 if pd.isna(cur) else cur) + lr[c] * frac[mi]
            ef.loc[li, cols] = np.nan                          # the lump itself is replaced by its members
    # ---- carbon basis
    co = ef[(ef.formula == 'CO') & (ef.pollutant_category == 'inorganic gas')][cols].iloc[0]
    co2 = ef[(ef.formula == 'CO2') & (ef.pollutant_category == 'inorganic gas')][cols].iloc[0]
    Cm = co2 * 12.011 / 44.01 + co * 12.011 / 28.01
    ok = [c for c in cols if pd.notna(Cm[c]) and Cm[c] > 0]
    Y = ef[ok].div(Cm[ok], axis=1)
    return ef, ok, mce[ok], Y, Cm[ok], co[ok], co2[ok]


def fit_all(fire_type):
    ef, cols, mce, Y, Cm, co, co2 = prepare(fire_type)
    rows = []
    x_all = mce.values
    for i in range(len(ef)):
        y = Y.iloc[i].values.astype(float)
        k = ~np.isnan(y)
        n = int(k.sum())
        rec = dict(row=i, formula=ef.at[i, 'formula'], compound=ef.at[i, 'compound'], category=ef.at[i, 'pollutant_category'],
                   id=ef.at[i, 'id'], n=n, method='none', a=np.nan, b=np.nan, r2=np.nan, p=np.nan, ymean=np.nan, mce_mean=np.nan,
                   span=np.nan, loo_rmse=np.nan, log=False)
        if n == 0:
            rows.append(rec); continue
        x = x_all[k]; yy = y[k]
        rec.update(ymean=float(yy.mean()), mce_mean=float(x.mean()), span=float(np.ptp(x)), method='mean')
        if n >= NMIN and np.ptp(x) >= SPANMIN:
            logf = ef.at[i, 'compound'] in LOGFIT and (yy > 0).all()
            yt = np.log(yy) if logf else yy
            r = stats.linregress(x, yt)
            if r.pvalue < PMAX:
                # leave-one-out check against the fire-type mean: keep the fit only when it predicts better
                e_fit, e_mean = [], []
                for j in range(n):
                    m = np.ones(n, bool); m[j] = False
                    rj = stats.linregress(x[m], yt[m])
                    e_fit.append(yt[j] - (rj.intercept + rj.slope * x[j])); e_mean.append(yt[j] - yt[m].mean())
                lf, lm = np.sqrt(np.mean(np.square(e_fit))), np.sqrt(np.mean(np.square(e_mean)))
                rec.update(loo_rmse=float(lf))
                if lf < lm:
                    rec.update(method='fit', a=float(r.intercept), b=float(r.slope), r2=float(r.rvalue ** 2), p=float(r.pvalue), log=bool(logf))
        rows.append(rec)
    return pd.DataFrame(rows), ef, cols, mce, Cm, co, co2


def evaluate(row, mce):
    """y (g X per g C) at the requested MCE"""
    if row['method'] == 'fit':
        v = row['a'] + row['b'] * mce
        return float(np.exp(v)) if row['log'] else max(float(v), 0.0)
    return float(row['ymean']) if row['n'] else 0.0
