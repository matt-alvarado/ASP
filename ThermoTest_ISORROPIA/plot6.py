import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE, ASP_LABEL
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from cases import *
TEMPS=[273.15,298.15,313.15]
MWI={'NH4':18.038,'NO3':62.004,'Cl':35.453,'SO4':96.06,'Na':22.99,'K':39.098}
F=pd.read_pickle(WORK+'/frac.pkl'); W=pd.read_pickle(WORK+'/water.pkl')
F=F[F['T'].isin([273.15,298.15,313.15])]; W=W[W['T'].isin([273.15,298.15,313.15])]
LAB={'NH4':'Ammonium (as NH$_4^+$)','NO3':'Nitrate (as NO$_3^-$)','Cl':'Chloride (as Cl$^-$)','SO4':'Sulfate (as SO$_4^{2-}$)','Na':'Sodium (as Na$^+$)','K':'Potassium (as K$^+$)'}
COL={273.15:'#0072B2',298.15:'#009E73',313.15:'#D55E00'}; col=lambda T: COL[T]
TITLE={'Urban':'Urban (ISORROPIA test case)','Marine':'Marine (ISORROPIA test case)','Remote':'Remote (ISORROPIA test case)',
 'Smoke_K':'Potassium: Timbavati sample 3, bulk (KCl, NH$_4$Cl, NH$_4$NO$_3$, (NH$_4$)$_2$SO$_4$)',
 'Recip_K_NH4':'Near-singular: K$_2$SO$_4$ + NH$_4$NO$_3$ (K/NH$_4$/NO$_3$/SO$_4$ reciprocal pair)',
 'Letovicite':'Near-singular: NH$_4$:SO$_4$ = 1.5 (letovicite composition)',
 'Recip_Na_NH4':'Near-singular: NaCl + NH$_4$NO$_3$ (Na/NH$_4$/Cl/NO$_3$ reciprocal pair)'}
STY=[('ASP_down',dict(lw=1.6)),('ASP_up',dict(lw=1.6,ls=':')),('ISO_stable',dict(lw=1.1,ls='--'))]
for case in CASES:
    t=totals(case); sp=[s for s in ION if t[s]>0]
    G=F[F.case==case].copy()
    for m,_ in STY: G[m+'_ug']=G[m]*G.total_umol*G.species.map(MWI)
    ymax=1.05*max(G[[m+'_ug' for m,_ in STY]].max().max(), max(t[s]*MWI[s] for s in sp))
    n=len(sp)+2; ncol=4; nrow=(n+ncol-1)//ncol
    fig,axs=plt.subplots(nrow,ncol,figsize=(15,3.8*nrow),constrained_layout=True); axs=axs.ravel()
    for k,s in enumerate(sp):
        ax=axs[k]; ax.set_title(LAB[s]); ax.grid(alpha=0.3); ax.set_xlim(0.18,0.92); ax.set_ylim(0,ymax)
        tot=t[s]*MWI[s]; ax.axhline(tot,color='0.6',lw=0.8); ax.text(0.19,tot,' total (gas + aerosol)',va='bottom',fontsize=7,color='0.4')
        for T in TEMPS:
            d=G[(G.species==s)&(G['T']==T)]
            for m,kw in STY: ax.plot(d.RH,d[m+'_ug'],color=col(T),**kw)
        ax.set_ylabel('aerosol, µg m$^{-3}$'); ax.set_xlabel('RH')
    ax=axs[len(sp)]; ax.set_title('Aerosol water (own scale)'); ax.grid(alpha=0.3); ax.set_xlim(0.18,0.92)
    for T in TEMPS:
        d=W[(W.case==case)&(W['T']==T)]
        for m,kw in STY: ax.plot(d.RH,d[m],color=col(T),**kw)
    ax.set_ylim(bottom=0); ax.set_ylabel('µg m$^{-3}$'); ax.set_xlabel('RH')
    lg=axs[len(sp)+1]; lg.axis('off')
    hs=[Line2D([],[],color='k',lw=1.6,label='ASP, RH scanned down (0.90→0.20)'),Line2D([],[],color='k',lw=1.6,ls=':',label='ASP, RH scanned up (0.20→0.90)'),
        Line2D([],[],color='k',lw=1.1,ls='--',label='ISORROPIA II, stable'),Line2D([],[],color='0.6',lw=0.8,label='total of the species (gas + aerosol)')]+[Line2D([],[],color=col(T),lw=3,label=f'{T:.2f} K') for T in TEMPS]
    lg.legend(handles=hs,loc='center',frameon=False,fontsize=9,ncol=2)
    for a in axs[len(sp)+2:]: a.axis('off')
    comp=', '.join(f'{s} {v:g}' for s,v in CASES[case])
    fig.suptitle(f'{TITLE[case]}: aerosol-phase amounts\ninput (µg m$^{{-3}}$): {comp}; {ASP_LABEL}; 273.15, 298.15 and 313.15 K; species panels share one y axis',fontsize=12)
    fig.savefig(f'{OUT}/ThermoTest_{case}_ugm3_3T.png',dpi=120); plt.close(fig)
print('ok')
