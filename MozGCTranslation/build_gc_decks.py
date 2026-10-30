#!/usr/bin/env python3
"""Build InorganicInputDecks/<Case>_GCHEM/ : hybrid GEOS-Chem fullchem v14.8.0 (non-halogen gas chemistry) + baseline lumped chemistry."""
import os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hybrid, gc_rates as g
from gc_kpp import read_eqn
from gc_rates import GCTranslator, SP, peroxy_flag, RENAME

BASE2GC = {'HCHO': 'CH2O', 'MEOH': 'MOH', 'Formic_Acid': 'HCOOH', 'ETHE': 'C2H4', 'Ethane': 'C2H6', 'Acetaldehyde': 'ALD2',
           'ETOH': 'EOH', 'Acetic_Acid': 'ACTA', 'CH3CO2H': 'ACTA', 'GLY': 'GLYX', 'Acetylene': 'C2H2', 'Propene': 'PRPE', 'Propane': 'C3H8',
           'Acetone': 'ACET', 'MCR': 'MACR', 'API': 'MTPA', 'LIM': 'LIMO', 'BEN': 'BENZ', 'TOL': 'TOLU', 'XYM': 'XYLE', 'XYP': 'XYLE',
           'XYO': 'XYLE', 'CH3OOH': 'MP', 'C2H5OOH': 'ETP', 'NOA': 'PROPNN', 'HOCH2CHO': 'GLYC', 'ACETOL': 'HAC', 'CH3CO3H': 'MAP',
           'CH3O2': 'MO2', 'CH3CO3': 'MCO3', 'n-Butane': 'ALK4', 'i-Butane': 'ALK4', 'Butadiene': 'C4H6', 'Propanal': 'RCHO',
           'C2H5NO3': 'ETNO3', 'IC3H7NO3': 'IPRNO3', 'NC3H7NO3': 'NPRNO3', 'CH3O2NO2': 'MPN',
           'C2H5O2': 'ETO2', 'NC3H7O2': 'A3O2', 'IC3H7O2': 'B3O2', 'CH3COCO2H': 'PYAC', 'HOCH2CO2H': 'HACTA', 'FURAN': 'FURA',
           'PROPACID': 'RCOOH', 'C2H5OOH': 'ETP', 'IC3H7OOH': 'RB3P', 'NC3H7OOH': 'RA3P', 'HOCH2OOH': 'HMHP'}
IC_SCALE = {'n-Butane': 0.9, 'i-Butane': 0.9}                       # ALK4 = C4-C5 alkanes (carbon-conserving, 4 C / 4.5 C)
COUNTER = {'ISOP': 'ISOSOA', 'TOLU': 'TOLSOA', 'XYLE': 'XYLSOA', 'MTPA': 'TERPSOA', 'MTPO': 'TERPSOA', 'LIMO': 'TERPSOA', 'PIN': 'TERPSOA'}
def counter_for(react):
    for n in react:
        if n in COUNTER and any(o in ('OH', 'O3', 'NO3') for o in react): return COUNTER[n]
    return None

def header(n_gc, n_keep, n_drop):
    return """!! GasChemicalMechanism.in - hybrid GEOS-Chem fullchem + baseline lumped chemistry
!! ASP
!! Translated by MozGCTranslation/gc_rates.py from the GEOS-Chem 'fullchem' KPP mechanism, version 14.8.0
!! (KPP/fullchem/fullchem.eqn, fullchem_RateLawFuncs.F90, rateLawUtilFuncs.F90; 1067 reactions, 356 species).
!!
!! Contents
!!   1. GEOS-Chem gas-phase reactions (%d ASP reactions).  A GEOS-Chem rate that is a sum of terms (CO + OH, HNO3 + OH, HO2 + HO2),
!!      or a photolysis with several channels, is split into one ASP reaction per term.
!!   2. Baseline (RACM2-derived, Williams et al. 2015) reactions that have at least one reactant that GEOS-Chem does not contain
!!      (%d kept, %d replaced by GEOS-Chem).  Baseline compounds GEOS-Chem treats explicitly are renamed to GEOS-Chem names
!!      (see GCHEM_SPECIES_MAPPING.csv).
!!   3. SOA-precursor counters (ISOSOA, TOLSOA, XYLSOA, TERPSOA) appended to the GEOS-Chem OH / O3 / NO3 reactions of isoprene,
!!      toluene, xylene, and the monoterpenes (MTPA, MTPO, LIMO, PIN), as in the baseline deck.  GEOS-Chem's own SOA schemes
!!      (aerosol-phase tracers SOAIE, SOAGX, LVOCOA, IONITA, ...) are not used.
!!
!! NOT INCLUDED from GEOS-Chem: all halogen chemistry (Cl, Br, I, F species and the CFC / HCFC / halon family), sea-salt
!!   (SALA*, SALC*, K_MT) and cloud / aqueous-sulfur reactions (K_CLD, K_MDL, SO2uptk*, HMS), the aerosol-phase species,
!!   and the dummy L* / P* production-loss tracers.  The gas-phase sulfate product SO4 is mapped to ASP's H2SO4.
!!   Photolysis: GEOS-Chem uses FAST-JX cross sections / quantum yields; ASP's PhotoRates.in holds 51 TUV columns, so each GEOS-Chem
!!   photolysis is mapped to the closest TUV column (analog compounds, e.g. all organic hydroperoxides use CH3OOH -> column 31,
!!   alkyl nitrates columns 34-37, HPALDs column 51); channel splits of HNO4 / MPN (33 / 67 %%) and MVK (three equal channels) follow GEOS-Chem.
!!
!! Heterogeneous reactions: GEOS-Chem computes uptake coefficients from aerosol composition, RH and T (HO2, NO2, NO3, N2O5, GLYX, MGLY,
!!   IEPOX, nitrates ...).  Here they are ASP first-order heterogeneous reactions (flag 3: fixed gamma, molecular weight; ASP's own
!!   surface area and gas diffusion): HO2 gamma 0.2, N2O5 0.02, NO2 5e-6, NO3 5e-3, GLYX 2.2e-3 (GEOS-Chem day 4.4e-3 / night 8e-6, mean),
!!   MGLY and PYAC 3.6e-7, IEPOX 0.01, HMML 3e-4, VOC nitrates as given in fullchem.eqn.  Aerosol-phase products are dropped.
!!
!! TEMPERATURE / PRESSURE: nothing is evaluated at a fixed T or P.  Every GEOS-Chem rate function (GCARR, GCJPL*, GC_NIT / GC_ALK, GC_RO2NO,
!!   GC_PAN, GC_ISO1/2, GC_EPO, ARRPLUS, TUNPLUS ...) is expressed through ASP rate flags; verify_gc.py compares every translated rate
!!   with the KPP Fortran functions (T = 220-320 K, P = 200-1013 mb), and probe_gc.py runs each new flag in ASP itself.
!!
!! Rate flags used (see the MCM deck headers for the earlier ones)
!!   1-body  0; col [; scale]  photolysis                2; A; B; C; D  A T^B exp(C/T)^D          3; gamma; MW; r   heterogeneous
!!           13; A0;n0;E0;Ai;ni;Ei;F;N;KeqA;KeqE  Troe fall-off / equilibrium                   17; A;B;C;D;E  max(A (D + T E) exp(B/T) exp(C/T^3), 0)
!!   2-body  A; B; C; D        A T^B exp(C/T)^D          1 / 2  RO2 + NO nitrate / alkoxy branches (Carter and Atkinson)
!!           4 / 8 / 10        k1 + k2 M ; k1/(1+k2) ; H2O (k1 + M k2)
!!           12 / 13  [a0;b0;c0;n;x0;y0] isoprene RO2 + NO nitrate / alkoxy (GC_NIT / GC_ALK)       14  [A;C;m]  A exp(C/T)/(1 + m M)
!!           15  [A1;C1;A2;C2]  M k1/(1+k2)  (DMS + OH)        16  [A0;n0;Ai;ni;F;A3;C3]  k3 (1 - kf/ki)  (JPL NO2 + O, HNO3 + OH, CO + OH)
!!           17  [A;C;K;fO2]  A exp(C/T)(X+K)/(2X+K), X = fO2 M  (GLYX + NO3)        18  [A;C;a;b;sel]  GLYC / HAC + OH branching
!!           19  [A;B;C;D;E]  as 1-body 17                      20 / 21  [a;b;c;d;e;f;g]  ISOP + OH  (GC_ISO1 / GC_ISO2)
!!   3-body  2; A; B; C; D     M A T^B exp(C/T)^D          4; A0;n0;E0;Ai;ni;Ei;F;N   Troe fall-off with explicit N
!! REQUIRES common_asp >= ba61060 (rate flags 12-21; photolysis scale e4888c4 and 3-body flag 4 40c1760 are earlier).
!!
""" % (n_gc, n_keep, n_drop)

def main():
    tr = GCTranslator(); R = read_eqn(); lines, idx, why = [], [], collections.Counter()
    ren = lambda n: RENAME.get(n, n)
    for r in R:
        ls, reason = tr.translate(r, ren=ren, counter=counter_for)
        lines += ls; idx += [r['n']] * len(ls)
        if not ls: why[reason.split(' ')[0] + ' ' + reason.split(' ')[1]] += 1
    species = collections.OrderedDict()
    for l in lines:
        for n in hybrid.side_names(l.split(' ; ')[0]):
            if n not in hybrid.NEUTRAL and n != 'PROD' and n not in ('ISOSOA', 'TOLSOA', 'XYLSOA', 'TERPSOA'): species[n] = 1
    inv = {v: k for k, v in RENAME.items()}
    info = {}
    for n in species:
        raw = inv.get(n, n); mw = SP.mw(raw)
        if mw is None: print('WARNING no MW for', n); mw = 100.0
        info[n] = (mw, peroxy_flag(raw))
    cfg = dict(suffix='_GCHEM', tag='GEOS-Chem', rename=BASE2GC, ic_scale=IC_SCALE, species_set=set(species), species_info=info,
               mapping_csv='GCHEM_SPECIES_MAPPING.csv', force_mw=set(species))
    with open(os.path.join(HERE, 'GC_reaction_index.csv'), 'w') as fh:
        fh.write('asp_line_in_block,geoschem_reaction\n')
        for i, n in enumerate(idx): fh.write('%d,%d\n' % (i + 1, n))
    print('GEOS-Chem: %d reactions -> %d ASP lines, %d species (%d peroxy); not translated: %s' %
          (len(R), len(lines), len(species), sum(1 for v in info.values() if v[1]), dict(why)))
    for cs in [a for a in sys.argv[1:] if a in hybrid.CASES] or hybrid.CASES:
        res = hybrid.build_case(cs, cfg, lines, header)
        print('%-10s kept %d / replaced %d baseline reactions; %d species (%d peroxy)' % (cs, res['n_keep'], res['n_drop'], res['n_species'], res['n_perox']))
if __name__ == '__main__': main()
