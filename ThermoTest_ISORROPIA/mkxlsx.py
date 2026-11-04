import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from config import WORK, OUT, DECKS, ISO_EXE, ASP_LABEL
import pandas as pd, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter as L
from openpyxl.drawing.image import Image
from openpyxl.comments import Comment
from cases import *
F=pd.read_pickle(WORK+'/frac.pkl'); W=pd.read_pickle(WORK+'/water.pkl'); S=pd.read_pickle(WORK+'/solver.pkl')
FONT='Arial'; f=lambda **k: Font(name=FONT, **k)
BLUE=f(color='0000FF'); BLACK=f(); BOLD=f(bold=True); HDR=f(bold=True,color='FFFFFF'); TITLE=f(bold=True,size=14)
HFILL=PatternFill('solid',fgColor='305496'); YEL=PatternFill('solid',fgColor='FFFF00'); GREY=PatternFill('solid',fgColor='F2F2F2')
wb=Workbook()
def header(ws,row,cols,widths=None):
    for j,c in enumerate(cols,1):
        x=ws.cell(row,j,c); x.font=HDR; x.fill=HFILL; x.alignment=Alignment(wrap_text=True,vertical='center',horizontal='center')
        if widths: ws.column_dimensions[L(j)].width=widths[j-1]
    ws.row_dimensions[row].height=45
def setfont(ws):
    for row in ws.iter_rows():
        for c in row:
            if c.font is None or c.font.name!=FONT: c.font=Font(name=FONT,bold=c.font.bold if c.font else False,color=c.font.color if c.font else None,size=c.font.size if c.font else 11)
SPN={'NH4':'Ammonium (NH4+ / NH3)','NO3':'Nitrate (NO3- / HNO3)','Cl':'Chloride (Cl- / HCl)','SO4':'Sulfate','Na':'Sodium','K':'Potassium'}
MOD=[('ASP_down','ASP, RH down (0.90→0.20)'),('ASP_up','ASP, RH up (0.20→0.90)'),('ASP_down_worstfirst','ASP down, worst-first loop'),
     ('ASP_down_preClaude','ASP down, before Claude changes'),('ISO_stable','ISORROPIA II stable'),('ISO_metastable','ISORROPIA II metastable')]

# ---------------- README
ws=wb.active; ws.title='README'
lines=[('ASP ThermoTest vs ISORROPIA II: gas–aerosol partitioning of the inorganic species',TITLE),('',None),
('What was run',BOLD),
('ASP: '+ASP_LABEL+' (common_asp; coupled Newton electrolyte solve), ThermoTest run mode, one bulk particle, T = 273.15–313.15 K every 5 K, P = 1000 mb.',None),
('  Down scan: initial equilibrium at RH 0.70, then RH 0.90 → 0.20 in steps of 0.02 (efflorescence branch). Up scan (ASP_THERMO_SCAN=up): initial equilibrium at RH 0.20, then 0.20 → 0.90 (deliquescence branch).',None),
('  Also on the down scan: the same code with the worst-first electrolyte loop (ASP_FEE_OUTER=worst), and the code before any Claude change (34b7ec7 + BinEdges build guard + IonReset fix + the letovicite output fix).',None),
('ISORROPIA II 2.3 (EPFL-LAPI, github.com/EPFL-LAPI/ISORROPIA-II), forward problem, default settings, stable (solids allowed) and metastable (liquid only), solved independently at each RH and T.',None),
('  Fountoukis, C. and Nenes, A.: ISORROPIA II, Atmos. Chem. Phys. 7, 4639–4659, 2007. Used under its academic, non-commercial licence.',None),
('',None),('Sheets',BOLD),
('Cases – case inputs (µg m⁻³), molar masses, totals per ion (µmol m⁻³, formulas).',None),
('RH0.90 – aerosol fraction of each volatile species at RH 0.90 (deliquesced), every model and temperature, with differences from ISORROPIA stable (formulas on Data).',None),
('MaxDiff – largest |difference in aerosol fraction| from ISORROPIA stable over RH, by phase state (formulas on Data).',None),
('PhaseRH – lowest RH at which each model holds water (formulas on Water).',None),
('Solver – coupled vs worst-first electrolyte solver diagnostics from the dump/replay harness.',None),
('Data – every case × T × RH × species point: aerosol fractions for all models (blue = model output), water, and formula columns.',None),
('Water – aerosol water (µg m⁻³) for all models at every point.',None),
('Plots – the PNG plots that go with this workbook.',None),
('',None),('Definitions',BOLD),
('Aerosol fraction = species in the aerosol (dissolved + solid) / total (gas + aerosol); gas fraction = 1 − aerosol fraction. Sulfate, sodium and potassium are entirely in the aerosol in every model.',None),
('A model state counts as wet when its aerosol water exceeds the threshold below; "both wet" / "both dry" / "one wet" compare ASP with ISORROPIA stable at the same point.',None),
('Blue text: model output or case input (hardcoded data). Black: formulas. Yellow: the setting you can change.',None)]
for i,(t,ft) in enumerate(lines,1):
    c=ws.cell(i,1,t); c.font=ft or BLACK; c.alignment=Alignment(wrap_text=True,vertical='top')
ws.column_dimensions['A'].width=150
r=len(lines)+2
ws.cell(r,1,'Wet threshold for aerosol water (µg m⁻³):').font=BOLD
ws.cell(r,2,0.01).font=BLUE; ws.cell(r,2).fill=YEL; ws.column_dimensions['B'].width=12
ws.cell(r,2).comment=Comment('Change to re-classify wet/dry points; MaxDiff and PhaseRH recalculate.','analysis')
THR=f"README!$B${r}"

# ---------------- Cases
wc=wb.create_sheet('Cases')
cols=['Case','Species','Input (µg m⁻³)','Molar mass (g mol⁻¹)','Input (µmol m⁻³)']+[f'{s} per mol' for s in ION]+[f'{s} (µmol m⁻³)' for s in ION]
header(wc,1,cols,[16,14,13,13,13]+[9]*6+[12]*6)
row=2; first=2
for case,lst in CASES.items():
    for sp,v in lst:
        wc.cell(row,1,case); wc.cell(row,2,sp); wc.cell(row,3,v).font=BLUE; wc.cell(row,4,MW[sp]).font=BLUE
        wc.cell(row,5,f'=C{row}/D{row}')
        for j,s in enumerate(ION):
            wc.cell(row,6+j,COMP[sp].get(s,0)).font=BLUE
            wc.cell(row,12+j,f'=E{row}*{L(6+j)}{row}')
        row+=1
last=row-1
wc.cell(row+1,1,'Molar masses: ASP AqPhaseChems.in (NaNO3 85.0, K2SO4 174.2566, KNO3 101.1083); inputs: default AerosolModes.in ISORROPIA cases, Timbavati sample 3, and the three near-singular test compositions.').font=f(italic=True)
r0=row+3; wc.cell(r0,1,'Totals per case (µmol m⁻³)').font=BOLD
header(wc,r0+1,['Case']+ION+['','',''])
CTOT={}
for i,case in enumerate(CASES):
    rr=r0+2+i; wc.cell(rr,1,case)
    for j,s in enumerate(ION):
        wc.cell(rr,2+j,f'=SUMIFS({L(12+j)}${first}:{L(12+j)}${last},$A${first}:$A${last},$A{rr})').number_format='0.0000'
    CTOT[case]=rr
for rr in range(2,last+1):
    for c in range(5,18): wc.cell(rr,c).number_format='0.0000' if c==5 or c>=12 else '0'
wc.freeze_panes='C2'

# ---------------- Data
wd=wb.create_sheet('Data')
D=F.merge(W[['case','T','RH','ASP_down','ASP_up','ISO_stable']].rename(columns={'ASP_down':'W_down','ASP_up':'W_up','ISO_stable':'W_iso'}),on=['case','T','RH'])
D=D.sort_values(['case','species','T','RH'],key=lambda c: c.map({k:i for i,k in enumerate(CASES)}) if c.name=='case' else c).reset_index(drop=True)
dcols=['Case','Species','T (K)','RH']+[m[1] for m in MOD]+['Water ASP down (µg m⁻³)','Water ASP up','Water ISORROPIA stable',
       '|Δ| ASP down vs ISO stable','|Δ| ASP up vs ISO stable','Phase: ASP down vs ISO','Phase: ASP up vs ISO']
header(wd,1,dcols,[13,8,8,6]+[11]*6+[11,10,11,11,11,11,11])
n=len(D)
for i,rw in enumerate(D.itertuples(index=False),start=2):
    wd.cell(i,1,rw.case); wd.cell(i,2,rw.species); wd.cell(i,3,rw.T); wd.cell(i,4,float(rw.RH))
    for j,(m,_) in enumerate(MOD):
        v=getattr(rw,m); c=wd.cell(i,5+j,None if pd.isna(v) else float(v)); c.font=BLUE; c.number_format='0.0000'
    for j,m in enumerate(('W_down','W_up','W_iso')):
        c=wd.cell(i,11+j,float(getattr(rw,m))); c.font=BLUE; c.number_format='0.000'
    wd.cell(i,14,f'=ABS(E{i}-I{i})').number_format='0.0000'
    wd.cell(i,15,f'=ABS(F{i}-I{i})').number_format='0.0000'
    wd.cell(i,16,f'=IF(AND(K{i}>{THR},M{i}>{THR}),"both wet",IF(AND(K{i}<={THR},M{i}<={THR}),"both dry","one wet"))')
    wd.cell(i,17,f'=IF(AND(L{i}>{THR},M{i}>{THR}),"both wet",IF(AND(L{i}<={THR},M{i}<={THR}),"both dry","one wet"))')
wd.freeze_panes='E2'; wd.auto_filter.ref=f'A1:Q{n+1}'
DR=lambda col: f"Data!${col}$2:${col}${n+1}"

# ---------------- Water
ww=wb.create_sheet('Water')
header(ww,1,['Case','T (K)','RH']+[m[1]+' (µg m⁻³)' for m in MOD],[13,8,6]+[13]*6)
Wd=W.sort_values(['case','T','RH'],key=lambda c: c.map({k:i for i,k in enumerate(CASES)}) if c.name=='case' else c).reset_index(drop=True)
for i,rw in enumerate(Wd.itertuples(index=False),start=2):
    ww.cell(i,1,rw.case); ww.cell(i,2,rw.T); ww.cell(i,3,float(rw.RH))
    for j,(m,_) in enumerate(MOD):
        v=getattr(rw,m); c=ww.cell(i,4+j,None if pd.isna(v) else float(v)); c.font=BLUE; c.number_format='0.000'
nw=len(Wd); ww.freeze_panes='D2'; ww.auto_filter.ref=f'A1:I{nw+1}'
WR=lambda col: f"Water!${col}$2:${col}${nw+1}"

# ---------------- RH0.90
wr=wb.create_sheet('RH0.90',1)
cols=['Case','Species','T (K)','ASP down','ASP up','ASP down worst-first','ASP down pre-Claude','ISO stable','ISO metastable','ASP down − ISO stable','ASP up − ISO stable']
header(wr,1,cols,[13,22,8,10,10,12,12,10,11,12,12])
row=2
for case in CASES:
    t=totals(case)
    for s in ('NH4','NO3','Cl'):
        if t[s]==0: continue
        for T in TEMPS:
            wr.cell(row,1,case); wr.cell(row,2,SPN[s]); wr.cell(row,3,T)
            sk=f'"{s}"'
            for j,col in enumerate(['E','F','G','H','I','J']):
                wr.cell(row,4+j,f'=IFERROR(AVERAGEIFS({DR(col)},{DR("A")},$A{row},{DR("B")},{sk},{DR("C")},$C{row},{DR("D")},0.9),"n/a")').number_format='0.000'
            wr.cell(row,10,f'=IFERROR(D{row}-H{row},"n/a")').number_format='+0.000;-0.000;0.000'
            wr.cell(row,11,f'=IFERROR(E{row}-H{row},"n/a")').number_format='+0.000;-0.000;0.000'
            row+=1
wr.cell(row+1,1,'Aerosol fraction at RH 0.90 (all models deliquesced). "n/a": the pre-Claude code stopped before finishing that scan (letovicite composition, negative molality).').font=f(italic=True)
wr.freeze_panes='D2'

# ---------------- MaxDiff
wm=wb.create_sheet('MaxDiff',2)
cols=['Case','Species','T (K)','ASP down: all RH','ASP down: both wet','ASP down: both dry','ASP down: one wet (# RH)','ASP up: all RH','ASP up: both wet','ASP up: both dry','ASP up: one wet (# RH)']
header(wm,1,cols,[13,22,8,11,11,11,12,11,11,11,12])
row=2
for case in CASES:
    t=totals(case)
    for s in ('NH4','NO3','Cl'):
        if t[s]==0: continue
        for T in TEMPS:
            wm.cell(row,1,case); wm.cell(row,2,SPN[s]); wm.cell(row,3,T)
            base=f'{DR("A")},$A{row},{DR("B")},"{s}",{DR("C")},$C{row}'
            for j,(dc,pc) in enumerate((('N','P'),('O','Q'))):
                c0=4+4*j
                wm.cell(row,c0,f'=_xlfn.MAXIFS({DR(dc)},{base})').number_format='0.000'
                for k,ph in enumerate(('both wet','both dry')):
                    wm.cell(row,c0+1+k,f'=IF(COUNTIFS({base},{DR(pc)},"{ph}")=0,"–",_xlfn.MAXIFS({DR(dc)},{base},{DR(pc)},"{ph}"))').number_format='0.000'
                wm.cell(row,c0+3,f'=COUNTIFS({base},{DR(pc)},"one wet")')
            row+=1
wm.cell(row+1,1,'Largest |aerosol fraction (ASP) − aerosol fraction (ISORROPIA stable)| over the 36 RH points, split by whether ASP and ISORROPIA hold water at that point; the last column of each block counts RH points where only one of them is wet (phase-transition offsets).').font=f(italic=True)
wm.freeze_panes='D2'

# ---------------- PhaseRH
wp=wb.create_sheet('PhaseRH',3)
cols=['Case','T (K)']+[f'{m[1]}: lowest wet RH' for m in MOD]+['ASP up − ISO stable','ASP down − ISO stable']
header(wp,1,cols,[13,8]+[14]*6+[12,12])
row=2
for case in CASES:
    for T in TEMPS:
        wp.cell(row,1,case); wp.cell(row,2,T)
        for j,col in enumerate(['D','E','F','G','H','I']):
            crit=f'{WR("A")},$A{row},{WR("B")},$B{row},{WR(col)},">"&{THR}'
            wp.cell(row,3+j,f'=IF(COUNTIFS({crit})=0,"never wet",_xlfn.MINIFS({WR("C")},{crit}))').number_format='0.00'
        wp.cell(row,9,f'=IFERROR(D{row}-G{row},"n/a")').number_format='+0.00;-0.00;0.00'
        wp.cell(row,10,f'=IFERROR(C{row}-G{row},"n/a")').number_format='+0.00;-0.00;0.00'
        row+=1
wp.cell(row+1,1,'Lowest RH (of 0.20–0.90) at which the model holds more water than the README threshold: the deliquescence RH for the ASP up scan, the efflorescence RH for the down scan.').font=f(italic=True)
wp.freeze_panes='C2'

# ---------------- Solver
wsv=wb.create_sheet('Solver',4)
cols=['RH scan','Case','T (K)','Electrolyte solver','Solves replayed','Fallbacks to worst-first','Singular Jacobians','DRH-ordered drops accepted','Calls with an unstable set resolved','Results with an unstable solid set','Results violating equilibrium (> 2× tol)','Solver CPU (s)','Violating fraction']
header(wsv,1,cols,[8,13,8,12,10,11,10,11,11,11,12,10,10])
Sd=S.sort_values(['scan','case','T','solver'],key=lambda c: c.map({k:i for i,k in enumerate(CASES)}) if c.name=='case' else c)
for i,rw in enumerate(Sd.itertuples(index=False),start=2):
    vals=[rw.scan,rw.case,rw.T,rw.solver,rw.records,rw.fallbacks,rw.singular,rw.drops,rw.resolved,rw.unstable,rw.violating,rw.cpu]
    for j,v in enumerate(vals,1):
        c=wsv.cell(i,j,v)
        if j>=5: c.font=BLUE
    wsv.cell(i,12).number_format='0.000'
    wsv.cell(i,13,f'=IF(E{i}=0,0,K{i}/E{i})').number_format='0.0%'
ns=len(Sd)+1
r=ns+2; wsv.cell(r,1,'Totals').font=BOLD
header(wsv,r+1,['RH scan','Solver','Solves','Fallbacks','Singular','Unstable sets','Violating','Violating fraction'])
for k,(sc,sv) in enumerate([(a,b) for a in ('down','up') for b in ('coupled','worst-first')]):
    rr=r+2+k; wsv.cell(rr,1,sc); wsv.cell(rr,2,sv)
    for j,col in enumerate(['E','F','G','J','K']):
        wsv.cell(rr,3+j,f'=SUMIFS(${col}$2:${col}${ns},$A$2:$A${ns},$A{rr},$D$2:$D${ns},$B{rr})')
    wsv.cell(rr,8,f'=IF(C{rr}=0,0,G{rr}/C{rr})').number_format='0.0%'
wsv.cell(r+7,1,'From the dump/replay harness (patch 27): every electrolyte solve of each ThermoTest scan dumped, then re-solved with the coupled solve and with the worst-first loop. A coupled fallback uses the worst-first result for that solve.').font=f(italic=True)
wsv.freeze_panes='E2'

# ---------------- Plots
wpl=wb.create_sheet('Plots')
wpl.cell(1,1,'Plots (also delivered as separate PNG files)').font=TITLE
imgs=['ThermoTest_summary_RH090_vs_T.png','ThermoTest_phase_transitions.png']+[f'ThermoTest_{c}.png' for c in CASES]
rowp=3
for im in imgs:
    wpl.cell(rowp,1,im).font=BOLD
    img=Image(OUT+'/'+im); sc=900/img.width; img.width*=sc; img.height*=sc
    wpl.add_image(img,f'A{rowp+1}'); rowp+=int(img.height/20)+4
for ws_ in wb.worksheets: setfont(ws_)
wb.save(OUT+'/ASP_ISORROPIA_ThermoTest_comparison.xlsx'); print('saved', n, nw)
