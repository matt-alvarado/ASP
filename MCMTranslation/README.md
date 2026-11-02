# MCM v3.3.1 input decks for ASP

`InorganicInputDecks/<Case>_MCM331/` for Williams, Otavi, Alaska and Timbavati: the **Master Chemical Mechanism v3.3.1**
(York; 5,832 species, 16,698 reactions) combined with the baseline lumped chemistry for compounds the MCM does not contain.

Rebuild:  `cd MCMTranslation && PYTHONPATH=. python3 build_mcm_decks.py [Case ...]`  
Check:    `PYTHONPATH=. python3 verify_rates.py`  (every translated rate vs the MCM Fortran expression)

## What the hybrid contains

| part | content |
|---|---|
| MCM block | 16,696 MCM reactions -> 17,219 ASP reactions (a rate that is a sum of terms becomes one ASP reaction per term) |
| kept baseline | 305 baseline reactions that have at least one reactant the MCM does not contain: lumped RACM2 VOC/RO2 chemistry (HC5, HC8, OLT, OLI, DIEN, KET, ALD, UALD, ROH, FURAN, i-butene, 2-butenal ...), NH3, HCN, the SOA-counter -> IVOC chains and IVOC aging |
| replaced baseline | 1,302 baseline reactions whose reactants are all in the MCM (their chemistry now comes from the MCM) |
| counters | ISOSOA, TOLSOA, XYLSOA, TERPSOA added to every MCM OH / O3 / NO3 reaction of isoprene, toluene, the xylenes, and alpha-pinene / beta-pinene / limonene (one counter per reaction event, as in the baseline) so the IVOC -> OA chain still sees these precursors |

Species: 5,950 per deck (1,260 peroxy radicals, flag 1; acyl peroxy radicals flag 2). The RO2 pool is MCM's own list (AtChem2
`peroxy-radicals_v3.3.1`, 1,224 species).  Baseline compounds that the MCM treats explicitly are renamed to MCM names
(`MCM_SPECIES_MAPPING.csv` in each deck; table `BASE2MCM` in `mcm2asp.py`): e.g. ISOP -> C5H8, TOL -> TOLUENE, XYM/XYP/XYO ->
MXYL/PXYL/OXYL, API -> APINENE, LIM -> LIMONENE, Ethane -> C2H6, Propene -> C3H6, Acetone -> CH3COCH3, GLY -> GLYOX, MCR -> MACR,
HNO4 -> HO2NO2, ... Initial / background values are carried over unchanged. MCM `SA` -> `H2SO4` and `HCL` -> `HCl` (hard-wired ASP names).

## Temperature and pressure

Nothing is evaluated at a fixed T or P: `GasChemicalMechanism.in` is valid at any T and P (only `PhotoRates.in` and the initial conditions are
case specific). Every MCM rate is carried by an ASP rate flag:

* Arrhenius / power-law terms, `KRO2NO`, `KRO2HO2`, `KAPHO2`, `KDEC` ... : flagless 2-body, 1-body flag 2
* terms `k*RO2`: 1-body flag 12 (RO2_T pool); `KMT` Troe fall-off: 1-body flag 13 (decompositions) and 3-body flag 4 (`A + B + M`)
* `KMT05`: 2-body flag 4; `KMT11` (OH + HNO3) and `KMT18` (DMS + OH): 2-body flag 5; CH3O2 + HO2 and CH3O2 + RO2 branches: flags 7, 8, 9
* HO2 + HO2 (+H2O): flag 4 + flag 10; O + O2 + M and O(1D) + M: hard-wired flags 9, 8 (same constants as the MCM expressions)
* `NO + NO + O2`, `O + SO2`: 3-body flag 2 (O2 = 0.21 M); `O2`, `H2O` terms as explicit reactants
* photolysis: flag 0 with the optional scale (below)

New ASP rate forms added to common_asp for this (all optional, existing decks unchanged):

| flag | form | commit |
|---|---|---|
| photolysis `0; col; scale` | optional third term multiplies the PhotoRates.in column (PhotoRates is hard-wired to 51 columns, the MCM needs 66 distinct column x factor combinations) | e4888c4 |
| 1-body `14; A; B; C; E` | A T^B exp(C/T) exp(E/T^3): MCM isoprene H-shift isomerisations (`EXP(1.0E8/TEMP**3)`) | 718cca3 |
| 3-body `5; A0;n0;E0;Ai;ni;Ei;Fa;Fb;Fc` | Troe fall-off with F(T) = Fa exp(-Fb/T) + exp(-T/Fc): MCM C2H2 + OH (`KMT17`) | 710a126 |

(1-body flag 13 and 3-body flag 4, used here too, are from common_asp 40c1760.)

**Verification.** `verify_rates.py` evaluates all 17,219 translated ASP reactions with the formulas of `GasPhaseChemistry.h` and the
original MCM Fortran expressions (constants file executed as written) on T = 230-320 K x P = 250-1013 mb (x H2O, RO2): 300,546
comparisons, **all within 0.2 %**. Air composition where the MCM writes N2/O2: 0.79 M / 0.21 M (as in ASP flags 8, 9).

## Photolysis

The MCM's own J parameterisation (34 rates, l cos^m exp(-n sec) clear-sky) is **not** used; each MCM J is tied to the TUV 5.4 column
of the *same process* in the case's `PhotoRates.in` (unchanged from the baseline deck, so the retained lumped reactions keep working):

`J_ASP(j) = TUV column(j) x scale`, `scale = J_MCM(j)/J_MCM(analog)` (both at zenith 0; 1 where the TUV column is the same process).

Analog columns: J_O3_O1D 2, J_O3_O3P 3, J_H2O2 5, J_NO2 6, J_NO3_NO 7, J_NO3_NO2 8, J_HONO 12, J_HNO3 13, HCHO radical/molecular
17/18, CH3CHO 19, C2H5CHO 22 (also C3H7CHO, i-C3H7CHO channels with scale 0.99/0.59/2.73), CH3COCH3 26, MEK 28, MVK 27
(split 50/50), MACR 25 (split 50/50; HPALD J_C5HPALD1 = col 25 x 256, MCM ratio), glyoxal radical 44 / molecular 45 (split 0.87/0.13),
MGLYOX 46, biacetyl 47, CH3OOH 31, alkyl nitrates 36 (IC3H7NO3 analog; scales 0.39 CH3NO3, 0.46 C2H5NO3, 0.60 n-C3H7NO3, 2.79 t-C4H9NO3),
NOA 38. Midday J(ASP)/J(MCM clear-sky) ranges 0.3-3 (site-specific TUV vs the MCM default); the largest departures are the
glyoxal channel split (TUV/JPL: radical channel dominant; MCM parameterisation: molecular dominant), acetone, and the MACR/HPALD columns.

## Not translated / differences

* MCM 33, 34 (HNO3 -> NA, N2O5 -> 2 NA, crude first-order uptake to nitrate aerosol): omitted, ASP computes gas-aerosol exchange.
* MCM 15631 (`C522CO3H + OH + hv = C522CO3 : J(J_CH3OOH)`): omitted; the published file uses a photolysis rate as a bimolecular k (apparent upstream error).
* No CO2 is produced (the MCM carries no CO2 species; the baseline CO2 tracer is a passive species).
* The ISOP/TOL/aromatic SOA counters are an ASP/Williams device, not MCM chemistry; MCM SOA precursors other than those listed carry no counter.
* Baseline chemistry kept for species that exist in the MCM only if a lumped reactant is involved; orphaned lumped intermediates (e.g. ISOHC5, RIP,
  HPALD1) are formed only by retained reactions and stay near zero.
* Runtime is much larger than the lumped mechanisms (see the timing table in the run report).

## Results (5 h, 2026-11-02; common_asp c506969)

All four decks run (rc = 0, 349-385 s on one core). dOA/dCO at 5 h (ug/m3 STP per ppmv): Williams 107.6, Otavi 153.7, Alaska 37.9, Timbavati 46.4,
within 2 % of Baseline since both use the IVOC-ladder OA. dO3/dCO at 5 h: 0.144, 0.154, 0.070, 0.206, the highest or near-highest of the nine mechanisms in every case.
See `../CB6r5Translation/all_mechs_oa_o3_compare.png` and `../CB6r5Translation/NewMechanisms_README.md`.
