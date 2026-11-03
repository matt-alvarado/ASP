import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE
# Test cases for the ASP ThermoTest vs ISORROPIA II comparison: species amounts in ug/m3
MW={'HNO3':63.018,'NH3':17.034,'HCl':36.458,'H2SO4':98.076,'NH4Cl':53.492,'NH4NO3':80.052,'(NH4)2SO4':132.144,
    'KCl':74.5483,'K2SO4':174.2566,'KNO3':101.1083,'NaCl':58.44,'NaNO3':85.0,'Na2SO4':142.0395}
COMP={'HNO3':{'NO3':1},'NH3':{'NH4':1},'HCl':{'Cl':1},'H2SO4':{'SO4':1},'NH4Cl':{'NH4':1,'Cl':1},
      'NH4NO3':{'NH4':1,'NO3':1},'(NH4)2SO4':{'NH4':2,'SO4':1},'KCl':{'K':1,'Cl':1},'K2SO4':{'K':2,'SO4':1},
      'KNO3':{'K':1,'NO3':1},'NaCl':{'Na':1,'Cl':1},'NaNO3':{'Na':1,'NO3':1},'Na2SO4':{'Na':2,'SO4':1}}
CASES={'Urban':[('NH3',3.4),('H2SO4',9.143),('HNO3',1.953)],
       'Marine':[('NH3',0.02),('H2SO4',0.510),('HNO3',0.163),('NaCl',5.01),('HCl',0.002)],
       'Remote':[('NH3',4.25),('H2SO4',11.270),('HNO3',0.145)],
       # potassium case: the deck's own bulk fire case (Timbavati sample 3, bulk neutralized)
       'Smoke_K':[('KCl',18.12),('NH4Cl',31.45),('NH4NO3',9.28),('(NH4)2SO4',12.48)],
       # near-singular for the coupled Newton solve: dependent solid sets
       'Recip_K_NH4':[('K2SO4',4.0),('NH4NO3',4.0)],        # KNO3 + (NH4)2SO4 vs NH4NO3 + K2SO4
       'Letovicite':[('NH3',2.555),('H2SO4',9.808)],          # NH4:SO4 = 1.5 exactly: (NH4)3H(SO4)2 = (NH4)2SO4 + NH4HSO4
       'Recip_Na_NH4':[('NaCl',3.0),('NH4NO3',4.0)]}          # NaCl + NH4NO3 vs NaNO3 + NH4Cl
TEMPS=[273.15+5*i for i in range(9)]
ION=['NH4','NO3','Cl','SO4','Na','K']
def totals(case):
    t={k:0.0 for k in ION}
    for s,v in CASES[case]:
        for k,n in COMP[s].items(): t[k]+=n*v/MW[s]
    return t   # umol/m3
