import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE
import numpy as np
from cases import *
MWI={'NH4':18.038,'NO3':62.004,'Cl':35.453,'SO4':96.06,'Na':22.99,'K':39.098}
GASK={'NH4':'NH3','NO3':'HNO3','Cl':'HCl','SO4':'H2SO4'}
def rdtab(f):
    L=[l for l in open(f).read().split('\n') if l.strip()]; h=L[0].split()
    return h, np.array([[float(x) for x in l.split()] for l in L[1:]])
def asp(tag,case,T):
    d=f'{WORK}/runs/{tag}/{case}_{T:.2f}/OutputFiles'
    ha,a=rdtab(f'{d}/AerosolOutThermo.txt'); hg,g=rdtab(f'{d}/GasOutThermo1.txt')
    out={'RH':a[:,0],'H2O':a[:,ha.index('TotH2O')]}
    ppb2umol=1000e2/(8.314462618*T)*1e-3
    for s in MWI:
        out['aer_'+s]=a[:,ha.index('Tot'+s)]/MWI[s]
        out['gas_'+s]=g[:,hg.index(GASK[s])]*ppb2umol if s in GASK else np.zeros(len(a))
    return out
def iso(case,T,mode):
    d=np.loadtxt(f'{WORK}/isoruns/{case}_{T:.2f}_{mode}.txt',usecols=range(38))
    t=totals(case); gas=d[:,1:4]*1e6; liq=d[:,4:19]*1e6
    out={'RH':d[:,0],'H2O':liq[:,7]*18.015,'sld':d[:,19:38]*1e6,'liq':liq}
    for k,s in enumerate(('NH4','NO3','Cl')): out['gas_'+s]=gas[:,k]; out['aer_'+s]=t[s]-gas[:,k]
    for s in ('SO4','Na','K'): out['aer_'+s]=np.full(len(d),t[s]); out['gas_'+s]=np.zeros(len(d))
    return out
