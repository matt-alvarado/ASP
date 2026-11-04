import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE, ASP_LABEL
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from cases import *
F=pd.read_pickle(WORK+'/frac.pkl'); G=F[np.isclose(F.RH,0.9)]
LAB={'NH4':'Ammonium / NH$_3$','NO3':'Nitrate / HNO$_3$','Cl':'Chloride / HCl'}
STY=[('ISO_stable','ISORROPIA II, stable',dict(color='k',lw=2.2)),
     ('ASP_down','ASP, RH down',dict(color='#D55E00',lw=1.8,marker='o',ms=3)),
     ('ASP_up','ASP, RH up',dict(color='#009E73',lw=1.4,ls='--',marker='D',ms=3)),
     ('ASP_down_worstfirst','ASP worst-first loop, RH down',dict(color='#0072B2',lw=1.2,ls=':',marker='s',ms=3)),
     ('ASP_down_preClaude','ASP before Claude changes, RH down',dict(color='#777777',lw=1.2,ls='-.',marker='^',ms=3))]
fig,axs=plt.subplots(len(CASES),3,figsize=(12,2.25*len(CASES)),constrained_layout=True,sharex=True)
for r_,case in enumerate(CASES):
    t=totals(case)
    for c_,s in enumerate(['NH4','NO3','Cl']):
        ax=axs[r_,c_]; ax.grid(alpha=0.3); ax.set_ylim(-0.03,1.03)
        if r_==0: ax.set_title(LAB[s])
        if c_==0: ax.set_ylabel(case.replace('_',' ')+'\naerosol fraction',fontsize=9)
        if t[s]==0: ax.text(0.5,0.5,'not in this case',ha='center',va='center',transform=ax.transAxes,color='0.5'); continue
        d=G[(G.case==case)&(G.species==s)].sort_values('T')
        for m,name,kw in STY: ax.plot(d['T'],d[m],label=name,**kw)
for a in axs[-1]: a.set_xlabel('Temperature (K)')
h,l=axs[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=3,bbox_to_anchor=(0.5,-0.035),frameon=False)
fig.suptitle(f'Fraction in the aerosol at RH 0.90 (deliquesced) vs temperature ({ASP_LABEL})',fontsize=12)
fig.savefig(OUT+'/ThermoTest_summary_RH090_vs_T.png',dpi=120,bbox_inches='tight')
