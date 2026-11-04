"""Collect the replay reports written by solver_diag.sh into WORK/solver.pkl (Solver sheet of the workbook)."""
import sys, os, re, pandas as pd; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from config import WORK, TAG_DOWN, TAG_UP
from cases import CASES, TEMPS
def parse(f):
    if not os.path.exists(f): return None
    s=open(f).read(); g=lambda pat: [int(x) for x in re.search(pat+r'\s+([\d\s]+)',s).group(1).split()]
    fb=g(r'coupled: fallbacks, singular, max iter, cycles')
    return dict(records=g('records')[0], fallbacks=fb[0], singular=fb[1], drops=g('coupled: DRH-ordered drops accepted')[0],
        resolved=g('calls with an unstable set resolved')[0], unstable=g(r'records with an unstable solid set: input, reference, replay')[2],
        violating=g(r'records violating by > 2 x tolerance: reference, replay')[1], cpu=float(re.search(r'solver CPU s\s+(\S+)',s).group(1)))
rows=[]
for tag,scan in ((TAG_DOWN,'down'),(TAG_UP,'up')):
    for c in CASES:
        for T in TEMPS:
            for solver,f in (('coupled','FEEReplay_cur.txt'),('worst-first','FEEReplay_worst.txt')):
                p=parse(f'{WORK}/runs/{tag}/{c}_{T:.2f}/{f}')
                if p: rows.append(dict(scan=scan,case=c,T=T,solver=solver,**p))
pd.DataFrame(rows).to_pickle(f'{WORK}/solver.pkl'); print(len(rows),'rows')
