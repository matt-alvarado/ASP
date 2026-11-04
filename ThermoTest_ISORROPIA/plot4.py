import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE, ASP_LABEL
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from cases import *
F=pd.read_pickle(WORK+'/frac.pkl'); W=pd.read_pickle(WORK+'/water.pkl')
LAB={'NH4':'Ammonium / NH$_3$','NO3':'Nitrate / HNO$_3$','Cl':'Chloride / HCl','SO4':'Sulfate / H$_2$SO$_4$','Na':'Sodium','K':'Potassium'}
cmap=matplotlib.colormaps['plasma']; col=lambda T: cmap((T-273.15)/45*0.9)
TITLE={'Urban':'Urban (ISORROPIA test case)','Marine':'Marine (ISORROPIA test case)','Remote':'Remote (ISORROPIA test case)',
 'Smoke_K':'Potassium: Timbavati sample 3, bulk (KCl, NH$_4$Cl, NH$_4$NO$_3$, (NH$_4$)$_2$SO$_4$)',
 'Recip_K_NH4':'Near-singular: K$_2$SO$_4$ + NH$_4$NO$_3$ (K/NH$_4$/NO$_3$/SO$_4$ reciprocal pair)',
 'Letovicite':'Near-singular: NH$_4$:SO$_4$ = 1.5 (letovicite composition)',
 'Recip_Na_NH4':'Near-singular: NaCl + NH$_4$NO$_3$ (Na/NH$_4$/Cl/NO$_3$ reciprocal pair)'}
STY=[('ASP_down',dict(lw=1.6)),('ASP_up',dict(lw=1.6,ls=':')),('ISO_stable',dict(lw=1.1,ls='--'))]
for case in CASES:
    t=totals(case); sp=[s for s in ION if t[s]>0]
    n=len(sp)+2; ncol=4; nrow=(n+ncol-1)//ncol
    fig,axs=plt.subplots(nrow,ncol,figsize=(15,3.8*nrow),constrained_layout=True); axs=axs.ravel()
    for k,s in enumerate(sp):
        ax=axs[k]; ax.set_title(LAB[s]); ax.grid(alpha=0.3); ax.set_xlim(0.18,0.92)
        for T in TEMPS:
            d=F[(F.case==case)&(F.species==s)&(F['T']==T)]
            for m,kw in STY: ax.plot(d.RH,d[m],color=col(T),**kw)
        ax.set_ylim(-0.03,1.03); ax.set_ylabel('fraction in aerosol'); ax.set_xlabel('RH')
        sec=ax.secondary_yaxis('right',functions=(lambda x:1-x,lambda x:1-x)); sec.set_ylabel('fraction in gas')
    ax=axs[len(sp)]; ax.set_title('Aerosol water'); ax.grid(alpha=0.3); ax.set_xlim(0.18,0.92)
    for T in TEMPS:
        d=W[(W.case==case)&(W['T']==T)]
        for m,kw in STY: ax.plot(d.RH,d[m],color=col(T),**kw)
    ax.set_ylim(bottom=0); ax.set_ylabel('µg m$^{-3}$'); ax.set_xlabel('RH')
    lg=axs[len(sp)+1]; lg.axis('off')
    hs=[Line2D([],[],color='k',lw=1.6,label='ASP, RH scanned down (0.90→0.20)'),Line2D([],[],color='k',lw=1.6,ls=':',label='ASP, RH scanned up (0.20→0.90)'),
        Line2D([],[],color='k',lw=1.1,ls='--',label='ISORROPIA II, stable')]+[Line2D([],[],color=col(T),lw=3,label=f'{T:.2f} K') for T in TEMPS]
    lg.legend(handles=hs,loc='center',frameon=False,fontsize=9,ncol=2)
    for a in axs[len(sp)+2:]: a.axis('off')
    comp=', '.join(f'{s} {v:g}' for s,v in CASES[case])
    fig.suptitle(f'{TITLE[case]}\ninput (µg m$^{{-3}}$): {comp}; {ASP_LABEL}',fontsize=12)
    fig.savefig(f'{OUT}/ThermoTest_{case}.png',dpi=120); plt.close(fig)
# phase transitions: lowest RH with water (deliquescence on the up scan / efflorescence on the down scan)
thr=0.01
fig,axs=plt.subplots(2,4,figsize=(15,7),constrained_layout=True,sharey=True); axs=axs.ravel()
for k,case in enumerate(CASES):
    ax=axs[k]; ax.set_title(case.replace('_',' ')); ax.grid(alpha=0.3)
    for m,lab,kw in (('ASP_down','ASP down: lowest wet RH (efflorescence)',dict(color='#D55E00',marker='o')),
                     ('ASP_up','ASP up: lowest wet RH (deliquescence)',dict(color='#0072B2',marker='s',ls=':')),
                     ('ISO_stable','ISORROPIA stable: lowest wet RH',dict(color='k',marker='^',ls='--'))):
        ys=[]
        for T in TEMPS:
            d=W[(W.case==case)&(W['T']==T)&(W[m]>thr)]; ys.append(d.RH.min() if len(d) else np.nan)
        ax.plot(TEMPS,ys,label=lab,**kw)
    ax.set_xlabel('Temperature (K)'); ax.set_ylim(0.15,0.95)
    if k%4==0: ax.set_ylabel('RH')
axs[7].axis('off'); h,l=axs[0].get_legend_handles_labels(); axs[7].legend(h,l,loc='center',frameon=False)
fig.suptitle('Lowest RH at which the aerosol holds water (> 0.01 µg m$^{-3}$)',fontsize=12)
fig.savefig(OUT+'/ThermoTest_phase_transitions.png',dpi=120); plt.close(fig)
print('ok')
