#!/usr/bin/env python3
"""Build InorganicInputDecks/<Case>_MOZT1/ : hybrid MOZART-T1 (WRF-Chem KPP T1_MOZCART) + baseline lumped chemistry."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hybrid, pure, lump_maps
from t1_kpp import read_eqn
from mozt1 import T1Translator, HET
from aspfmt import fmt

T1_RENAME = {'SO4': 'H2SO4', 'RO2': 'RO2_T1'}                       # T1 name -> ASP name
COUNTER = {'ISOP': 'ISOSOA', 'TOLUENE': 'TOLSOA', 'XYLENES': 'XYLSOA', 'APIN': 'TERPSOA', 'BPIN': 'TERPSOA', 'LIMON': 'TERPSOA',
           'MYRC': 'TERPSOA', 'BCARY': 'TERPSOA'}
PEROX1 = set('''CH3O2 C2H5O2 C3H7O2 PO2 EO2 ALKO2 ENEO2 HOCH2OO ISOPAO2 ISOPBO2 ISOPNO3 MACRO2 MBOO2 MBONO3O2 HMPROPO2 BENZO2 BZOO
C6H5O2 DICARBO2 MALO2 MDIALO2 PHENO2 TOLO2 XYLENO2 XYLOLO2 NTERPO2 TERPO2 TERP2O2 MEKO2 XO2 RO2_T1'''.split())
PEROX2 = {'CH3CO3', 'MCO3', 'ACBZO2'}
# approximate molecular weights (g/mol) for species without a baseline row (lumped MOZART species are approximate)
MW = dict(ALKNIT=119, ALKO2=117, ALKOOH=118, BCARY=204, BENZOOH=127, BEPOMUC=112, BIGALD=98, BIGALD1=84, BIGALD2=100, BIGALD3=98,
          BIGALD4=98, BIGALK=72, BIGENE=56, BPIN=136, BZOO=107, BZOOH=108, C3H7OOH=76, C6H5O2=109, C6H5OOH=110, CH3CN=41, DICARBO2=117,
          DMS=62, ENEO2=133, EO=61, EO2=77, EOOH=78, HMPROP=74, HMPROPO2=105, HONITR=119, HPALD=116, HYDRALD=100, ISOPAO2=117,
          ISOPNITA=147, ISOPNITB=147, ISOPNO3=147, ISOPNOOH=148, ISOPOOH=118, MACRO2=119, MACROOH=120, MALO2=115, MBO=86,
          MBONO3O2=164, MBOO2=117, MBOOOH=118, MCO3=101, MDIALO2=131, MEKO2=103, MEKOOH=104, MYRC=136, N2O=44, NC4CH2OH=149,
          NC4CHO=147, NTERPO2=230, NTERPOOH=231, ONITR=119, PBZNIT=183, PHENO=93, PHENO2=125, PHENOOH=126, PO2=91, POOH=92,
          RO2_T1=89, ROOH=90, TEPOMUC=140, TERP2O2=201, TERP2OOH=202, TERPNIT=215, TERPO2=185, TERPOOH=186, TERPROD1=168,
          TERPROD2=170, TOLO2=141, TOLOOH=142, XO2=133, XOOH=134, XYLENES=106, XYLENO2=155, XYLENOOH=156, XYLOL=122, XYLOLO2=153,
          XYLOLOOH=154, ACBZO2=149, ALKNIT_=0)
OMIT = {}

def counter_for(react):
    for n in react:
        if n in COUNTER and any(o in ('OH', 'O3', 'NO3') for o in react): return COUNTER[n]
    return None

def header(n_t1):
    return """!! GasChemicalMechanism.in - MOZART-T1 + SOA-precursor / IVOC chains
!! ASP
!! Translated by MozGCTranslation/mozt1.py from the WRF-Chem KPP mechanism T1_MOZCART (tag MZ148_T1mozcart_20170323; 142 species,
!! 344 reactions; the MOZART-4 / CAM-chem 'T1' tropospheric chemistry with the Emmons et al. (2010) and Knote et al. (2014) extensions).
!! Note: this is the 2017 WRF-Chem T1 (chem_opt 114), not the later CAM-chem MOZART-T1.1 / T1.2 updates.
!!
!! Contents
!!   1. T1 gas-phase reactions (%d ASP reactions; a T1 rate that is a sum of terms or of photolysis channels is split into one ASP
!!      reaction per term / TUV column).  T1 reaction numbers are in T1_reaction_index.csv.
!!   2. SOA-precursor counters (ISOSOA, TOLSOA, XYLSOA, TERPSOA) appended to the T1 OH / O3 / NO3 reactions of isoprene, toluene,
!!      xylenes and the terpenes, as in the baseline deck; their counter -> IVOC chains and the IVOC aging block are copied unchanged
!!      from the baseline deck (only Acetic_Acid -> CH3COOH).  T1's own VBS SOA yield function (vbs_yield) is not used.
!!   NO baseline RACM2-derived lumped gas reactions are kept: the baseline initial / background compounds are lumped onto T1 species
!!   (MOZT1_SPECIES_MAPPING.csv).
!!   Heterogeneous reactions: the T1 aerosol-uptake reactions (usr_*_aer: HO2, N2O5, NO2, NO3, HONITR, ISOPNITA/B, ONITR, TERPNIT,
!!      NC4CH2OH, NC4CHO, NTERPOOH, GLYOXAL) are expressed with ASP's first-order heterogeneous flag 3 (gamma, molecular weight), which
!!      uses ASP's own particle surface area and gas-phase diffusion (same kinetic-diffusion form as the T1 functions; T1 uses a fixed
!!      Dg = 0.1 cm2/s and its own prescribed sulfate / OC / BC surface areas).  gamma and molecular weight are those of T1
!!      (MW from the T1 mean-speed coefficient).
!!
!! TEMPERATURE / PRESSURE: nothing is evaluated at a fixed T or P.  Arrhenius, JPL fall-off (Troe, explicit M), equilibrium-
!! controlled decompositions (TROEE), the T1 user functions and photolysis (PhotoRates.in column x optional scale) are all expressed
!! through ASP rate flags.  verify_t1.py compares every translated rate with the KPP Fortran functions (T = 220-320 K, P = 200-1013 mb).
!!
!! Rate flags used (see MCM deck headers for the full list):
!!   1-body  0; col [; scale]  photolysis                  2; A; B; C; D   A T^B exp(C/T)^D
!!           3; gamma; MW; r   first-order heterogeneous    13; A0;n0;E0;Ai;ni;Ei;F;N;KeqA;KeqE  fall-off / equilibrium (TROEE)
!!   2-body  A; B; C; D        A T^B exp(C/T)^D            4 / 5 / 10 / 11  (k1 + k2 M ; HNO3 + OH form ; H2O term ; CO + OH pressure form)
!!   3-body  2; A; B; C; D     M A T^B exp(C/T)^D          4; A0;n0;E0;Ai;ni;Ei;F;N   Troe fall-off with explicit N (common_asp 40c1760)
!! REQUIRES common_asp >= e4888c4 (photolysis scale) and 40c1760 (3-body flag 4).
!!
""" % n_t1

def main():
    tr = T1Translator(); R = read_eqn()
    ren_t1 = lambda n: T1_RENAME.get(n, n)
    lines, idx = [], []
    for r in R:
        ls = tr.translate(r, ren=ren_t1, counter=counter_for)
        lines += ls; idx += [r['n']] * len(ls)
    species = collections.OrderedDict()
    for r in R:
        for c, n in r['re'] + r['pr']:
            n = ren_t1(n)
            if n not in hybrid.NEUTRAL and n != 'PROD': species[n] = 1
    info = {n: (MW.get(n, 100.0), 2 if n in PEROX2 else (1 if n in PEROX1 else 0)) for n in species}
    cfg = dict(suffix='_MOZT1', tag='MOZART-T1', lump=lump_maps.T1, species_info=info, counters=['TOLSOA', 'XYLSOA', 'ISOSOA', 'TERPSOA'],
               acetic='CH3COOH', mapping_csv='MOZT1_SPECIES_MAPPING.csv')
    with open(os.path.join(HERE, 'T1_reaction_index.csv'), 'w') as fh:
        fh.write('asp_line_in_block,t1_reaction\n')
        for i, n in enumerate(idx): fh.write('%d,%d\n' % (i + 1, n))
    print('T1: %d reactions -> %d ASP lines, %d species' % (len(R), len(lines), len(species)))
    for cs in [a for a in sys.argv[1:] if a in hybrid.CASES] or hybrid.CASES:
        res = pure.build_case(cs, cfg, lines, header)
        print('%-10s %d reactions + IVOC chains; %d species (%d peroxy)' % (cs, len(lines), res['n_species'], res['n_perox']))

if __name__ == '__main__': main()
