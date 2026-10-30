# MOZART-T1 and GEOS-Chem translations (ASP decks)

The decks contain the translated mechanism only, plus the SOA-counter -> IVOC chains and the IVOC aging block of the baseline deck. No baseline
RACM2-derived lumped gas reactions are kept (that policy is only for the MCM decks). The paper compounds of each case are lumped onto
mechanism species in `lump_maps.py` (mole-for-mole by OH reactivity / structure, as in the CRACMM and SAPRC decks; carbon is not conserved
for lumps smaller than the compound); `pure.py` builds the decks and writes `*_SPECIES_MAPPING.csv`.

* `t1_kpp.py`, `mozt1.py`, `build_t1_decks.py`, `verify_t1.py`: WRF-Chem KPP T1_MOZCART (tag MZ148_T1mozcart_20170323) -> ASP.
  `python3 build_t1_decks.py` writes `InorganicInputDecks/<Case>_MOZT1/`; `python3 verify_t1.py` checks every translated rate
  against the KPP Fortran functions (T 220-320 K, P 200-1013 mb).
* `pure.py` / `lump_maps.py`: deck builder and species lumping. `hybrid.py` (keep baseline reactions whose reactant is not in the mechanism) is no longer used by these two.  `aspfmt.py`: reaction-line formatting.
* This T1 is the 2017 WRF-Chem version (142 species, 344 reactions), not the later CAM-chem T1.1 / T1.2.
* Rates are structural (no T or P frozen).  T1 aerosol-uptake reactions (usr_*_aer) use ASP's first-order heterogeneous flag 3
  (gamma, MW) so they follow ASP's own particle surface area; T1's fixed Dg = 0.1 cm2/s and prescribed sulfate/OC/BC areas are not used.
* T1 species without a baseline row use approximate molecular weights (`build_t1_decks.MW`; lumped MOZART species are approximate).
* T1's VBS SOA yield function is not used; ASP's SOA counters are re-injected on first-generation VOC reactions.

## GEOS-Chem fullchem v14.8.0

* `gc_kpp.py` (fullchem.eqn parser), `gc_species.py` (formula -> molecular weight, exclusions), `gc_rates.py` (rate-function -> ASP flag
  translator), `build_gc_decks.py` (`InorganicInputDecks/<Case>_GCHEM/`), `verify_gc.py` (every translated rate vs. the KPP Fortran
  functions, T 220-320 K, P 200-1013 mb: 0 of 9495 evaluations outside 0.2 %), `probe_gc.py` + `probe_check.py` (run each new rate
  flag inside ASP and compare the measured decay rate with the Python mirror; agreement within 0.6 %).
* 1067 GEOS-Chem reactions -> 793 ASP lines.  Not translated: halogen (Cl, Br, I, F; incl. CFC/HCFC/halons), sea-salt, cloud, aqueous-sulfur
  and aerosol-phase species/reactions (286 reactions), three photolysis channels without a TUV analog (GLYX H2+2CO, SO4, NO).
* New common_asp rate flags (all T/P-structural): 2-body 12/13 (GC_NIT/GC_ALK), 14 (GC_EPO), 15 (DMS+OH), 16 (JPL activation, GCJPLAC),
  17 (GLYX+NO3), 18 (GLYC/HAC+OH branching), 19 (ARRPLUS), 20/21 (GC_ISO1/2, ISOP+OH); 1-body 17 (ARRPLUS/TUNPLUS, floored at zero).
* Names: SO4 -> H2SO4, HNO2 -> HONO (ASP's enhancement-ratio output looks these up by name; CH3NO3 stays as a placeholder species there).
  Paper compounds are lumped onto GEOS-Chem species by `lump_maps.GC` (butanes, pentanes, hexanes ... -> ALK4, >=C3 alkenes -> PRPE, ...).
* Photolysis: GEOS-Chem uses FAST-JX; each photolysis is mapped to the closest TUV column of PhotoRates.in (`gc_rates.PH`).  This is the
  largest departure from GEOS-Chem and is approximate for the many organic hydroperoxides / nitrates.
* Heterogeneous uptake: ASP flag 3 with fixed gamma (`gc_rates.HET_GAMMA`) instead of GEOS-Chem's composition/RH-dependent gammas.
* Run times (5 h, 2 cores, in pairs): Williams 11 s, Otavi 10 s, Alaska 36 s, Timbavati 40 s (MOZART-T1: 11 / 9 / 34 / 36 s). With the baseline reactions kept they were 18 / 16 / 41 / 45 s.
