"""Make ThermoTest run directories: python3 mkcase.py <tag> writes WORK/runs/<tag>/<case>_<T>/ from the default deck.
Tags ending in "up" start at RH 0.20 (for ASP_THERMO_SCAN=up)."""
import os, re, shutil, sys
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE
from cases import CASES, TEMPS
SRC=DECKS
def mk(case, d, T=298.15, rh0=None):
    deck=f'{d}/InorganicInputDecks'; os.makedirs(f'{d}/OutputFiles', exist_ok=True)
    os.makedirs(deck, exist_ok=True)
    for f in os.listdir(SRC):
        if f.endswith('.in') and os.path.isfile(f'{SRC}/{f}'): shutil.copy(f'{SRC}/{f}', deck)
    # AerosolModes: mass flag 1, one bulk mode with the case
    L=open(f'{deck}/AerosolModes.in','rb').read().decode('latin-1').split('\r\n')
    assert L[49].startswith('0 !!'); L[49]='1 !!'
    assert L[91].startswith('234000'); j=91
    while not L[j].startswith('!'): j+=1
    L[91:j]=['0 ; 0.2 ; 1.4;', f'{case} Case ISORROPIA;']+[f'{s}; {v};' for s,v in CASES[case]]
    open(f'{deck}/AerosolModes.in','wb').write('\r\n'.join(L).encode('latin-1'))
    # ASPInputDeck: ThermoTest on, parcel off, no coagulation, T
    L=open(f'{deck}/ASPInputDeck.in','rb').read().decode('latin-1').split('\r\n')
    def setline(n,v,check):
        assert check in L[n-1], (n,L[n-1]); L[n-1]=re.sub(r'^\S+', v, L[n-1])
    setline(35,'0','Include Coagulation'); setline(45,'0','Lagrangian Parcel'); setline(49,'1','Thermodynamics test')
    setline(63,f'{T}','Temperature')
    if rh0 is not None:
        setline(65,f'{rh0}','Initial Relative Humidity'); setline(66,f'{rh0}','Background Relative Humidity')
    open(f'{deck}/ASPInputDeck.in','wb').write('\r\n'.join(L).encode('latin-1'))
    # GasPhaseChems: no NH3/HNO3/HCl/H2SO4 in the gas initially (totals come from AerosolModes)
    s=open(f'{deck}/GasPhaseChems.in','rb').read().decode('latin-1'); out=[]
    for l in s.split('\n'):
        m=re.match(r'^(NH3|HNO3|HCl|H2SO4)(\s*;)',l)
        if m:
            p=l.split(';'); p[4]='0'; p[5]='0'; l=';'.join(p)
        out.append(l)
    open(f'{deck}/GasPhaseChems.in','wb').write('\n'.join(out).encode('latin-1'))
if __name__=='__main__':
    tag=sys.argv[1]
    for c in CASES:
        for T in TEMPS: mk(c, f'{WORK}/runs/{tag}/{c}_{T:.2f}', T, 0.2 if tag.endswith('up') else None)
