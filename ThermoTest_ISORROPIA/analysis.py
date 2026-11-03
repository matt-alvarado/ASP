import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE
from load2 import *
import pandas as pd
MODELS=[('ASP_down','asp','p34'),('ASP_up','asp','up'),('ASP_down_worstfirst','asp','worst'),('ASP_down_preClaude','asp','base'),
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
