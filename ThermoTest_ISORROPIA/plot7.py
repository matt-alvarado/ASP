import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from config import WORK, OUT, ASP_LABEL
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from cases import CASES
S=pd.read_pickle(WORK+'/solids.pkl')
fig,axs=plt.subplots(3,7,figsize=(20,8.5),constrained_layout=True,sharey='col')
for c_,case in enumerate(CASES):
    for r_,T in enumerate((273.15,298.15,313.15)):
        ax=axs[r_,c_]; d=S[(S.case==case)&(S['T']==T)]
        d=d[(d[['ISO_stable','ASP_down','ASP_up']].fillna(0).abs().sum(axis=1)>1e-6)]
        x=np.arange(len(d)); w=0.27
        ax.bar(x-w,d.ISO_stable,w,color='k',label='ISORROPIA II stable')
        ax.bar(x,d.ASP_down,w,color='#D55E00',label='ASP, RH down')
        ax.bar(x+w,d.ASP_up,w,color='#009E73',label='ASP, RH up')
        ax.set_xticks(x); ax.set_xticklabels(d.solid,rotation=60,ha='right',fontsize=8)
        if r_==0: ax.set_title(case.replace('_',' '))
        if c_==0: ax.set_ylabel(f'{T:.2f} K\nsolid, µmol m$^{{-3}}$')
        ax.grid(alpha=0.3,axis='y')
h,l=axs[0,0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=3,bbox_to_anchor=(0.5,-0.03),frameon=False)
fig.suptitle(f'Solids in the dry particle at RH 0.20 ({ASP_LABEL})',fontsize=13)
fig.savefig(OUT+'/ThermoTest_dry_solids_RH020.png',dpi=110,bbox_inches='tight'); print('ok')
