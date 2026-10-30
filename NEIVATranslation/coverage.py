"""Coverage report: how every NEIVA temperate-forest gas row maps onto the baseline (Williams) deck. python3 coverage.py"""
import sys, collections
import pandas as pd, numpy as np
import mce_fit as mf, neiva_data as nd, neiva_map as nm

def classified(fire_type='temperate forest'):
    fit, ef, cols, mce, Cm, co, co2 = mf.fit_all(fire_type)
    prop = nd.load_property().drop_duplicates('id').set_index('id')
    rows = []
    for _, r in fit.iterrows():
        if r.category != 'NMOC_g' or r.n == 0:
            continue
        p = prop.loc[r['id']] if r['id'] in prop.index else None
        voc = p.voc_type if p is not None else ''
        smile = p.smile if p is not None else ''
        kOH = p.kOH if p is not None and p.kOH == p.kOH else None
        cs = p.cstar if p is not None and p.cstar == p.cstar else None
        mm = ef.at[int(r.row), 'mm']
        route, tgt, nC = nm.classify(r.formula, r.compound, voc, smile, kOH, cs)
        if route == 'ivoc': tgt = nm.ivoc_bin(cs)
        rows.append(dict(r, voc_type=voc, smile=smile, kOH=kOH, cstar=cs, mm=mm, route=route, target=tgt, nC=nC))
    return pd.DataFrame(rows), fit, ef, cols, mce, Cm, co, co2

if __name__ == '__main__':
    d = classified()[0]
    d['ym'] = d.ymean
    print(d.groupby('route').agg(n=('ym', 'size'), y=('ym', 'sum')).assign(share=lambda x: x.y / x.y.sum()))
    print(d[d.route == 'lump'].groupby('target').ym.agg(['size', 'sum']))
    print(d[d.route == 'ivoc'].groupby('target').ym.agg(['size', 'sum']))
    print(d[d.route == 'skip'].groupby('target').ym.agg(['size', 'sum']))
    print(d[d.route == 'explicit'][['compound', 'target', 'n', 'method', 'ym']].to_string())
    print(d[(d.route == 'ivoc') & (d.target == 9)].sort_values('ym', ascending=False)[['compound', 'formula', 'cstar', 'ym']].head(12).to_string())
    print(d[(d.route == 'ivoc')].sort_values('ym', ascending=False)[['compound', 'formula', 'voc_type', 'target', 'ym']].head(15).to_string())
