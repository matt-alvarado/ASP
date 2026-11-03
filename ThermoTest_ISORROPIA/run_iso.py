"""Run ISORROPIA II (stable and metastable) for every case and temperature; RH 0.90 -> 0.20 in steps of 0.02."""
import sys, os, subprocess; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from config import WORK, ISO_EXE
from cases import CASES, TEMPS, totals
os.makedirs(f'{WORK}/isoruns', exist_ok=True)
for c in CASES:
    t=totals(c)
    for T in TEMPS:
        for m,name in ((0,'stable'),(1,'meta')):
            w=[t['Na'],t['SO4'],t['NH4'],t['NO3'],t['Cl'],0,t['K'],0]          # umol/m3 -> mol/m3 below
            inp=' '.join(f'{x*1e-6:.10e}' for x in w)+f' {T} {m}'
            out=subprocess.run([ISO_EXE],input=inp,capture_output=True,text=True,check=True).stdout
            open(f'{WORK}/isoruns/{c}_{T:.2f}_{name}.txt','w').write(out)
print('done')
