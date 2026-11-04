import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE, TAG_DOWN, TAG_UP
from load2 import *
import pandas as pd
MODELS=[('ASP_down','asp',TAG_DOWN),('ASP_up','asp',TAG_UP),('ASP_down_worstfirst','asp','worst'),('ASP_down_preClaude','asp','base'),
        ('ISO_stable','iso','stable'),('ISO_metastable','iso','meta')]
RHS=np.round(np.arange(0.20,0.9001,0.02),2)
def get(kind,tag,case,T):
    r=asp(tag,case,T) if kind=='asp' else iso(case,T,tag)
    rh=np.round(r['RH'],2)
    if len(rh)!=36: return None
    o=np.argsort(rh); return {k:(v[o] if isinstance(v,np.ndarray) and v.ndim==1 and len(v)==36 else v) for k,v in r.items()}
def build():
    frac=[]; water=[]
    for case in CASES:
        t=totals(case)
        for T in TEMPS:
            R={name:get(kind,tag,case,T) for name,kind,tag in MODELS}
            for k,rh in enumerate(RHS):
                water.append(dict(case=case,T=T,RH=rh,**{n:(R[n]['H2O'][k] if R[n] is not None else np.nan) for n,_,_ in MODELS}))
                for s in ION:
                    if t[s]==0: continue
                    frac.append(dict(case=case,T=T,RH=rh,species=s,total_umol=t[s],
                        **{n:(R[n]['aer_'+s][k]/t[s] if R[n] is not None else np.nan) for n,_,_ in MODELS}))
    return pd.DataFrame(frac), pd.DataFrame(water)
if __name__=='__main__':
    F,W=build(); print(F.shape, W.shape); print(F.isna().sum().to_dict())
    F.to_pickle(WORK+'/frac.pkl'); W.to_pickle(WORK+'/water.pkl')

# ---- dry state: solid composition at RH 0.20 (ASP species columns are mol cm-3 -> umol m-3 x 1e12)
SOLIDS=['NH4NO3','NH4Cl','(NH4)2SO4','NH4HSO4','(NH4)3H(SO4)2','KNO3','KCl','K2SO4','KHSO4','NaNO3','NaCl','Na2SO4','NaHSO4']
ISO_SLD={'NaNO3':0,'NH4NO3':1,'NaCl':2,'NH4Cl':3,'Na2SO4':4,'(NH4)2SO4':5,'NaHSO4':6,'NH4HSO4':7,'(NH4)3H(SO4)2':8,
         'K2SO4':12,'KHSO4':13,'KNO3':14,'KCl':15}
def asp_solids(tag,case,T,rh=0.20):
    f=f'{WORK}/runs/{tag}/{case}_{T:.2f}/OutputFiles/AerosolOutThermo.txt'
    L=[l for l in open(f).read().split('\n') if l.strip()]; h=L[0].split()
    for l in L[1:]:
        r=[float(x) for x in l.split()]
        if abs(r[0]-rh)<1e-6: return {s:r[h.index(s)]*1e12 for s in SOLIDS}, r[h.index('TotH2O')]
    return None, None
def build_solids():
    rows=[]
    for case in CASES:
        for T in TEMPS:
            i=iso(case,T,'stable'); k=list(np.round(i['RH'],2)).index(0.2)
            for s in SOLIDS:
                row=dict(case=case,T=T,solid=s,ISO_stable=float(i['sld'][k][ISO_SLD[s]]))
                for name,tag in (('ASP_down',TAG_DOWN),('ASP_up',TAG_UP)):
                    v,_=asp_solids(tag,case,T); row[name]=v[s] if v else np.nan
                rows.append(row)
    return pd.DataFrame(rows)
if __name__=='__main__':
    S=build_solids(); S.to_pickle(WORK+'/solids.pkl'); print('solids', S.shape)
