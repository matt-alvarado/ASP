# CRACMM2, CRACMM3, SAPRC07TC and SAPRC07TIC decks

Sixteen decks, `InorganicInputDecks/<Case>_<MECH>/`, translated from the CMAQ mechanism definitions
(`CCTM/src/MECHS/{cracmm2,cracmm3,saprc07tc_ae6_aq,saprc07tic_ae7i_aq}/mech_*.def`) with
`cmaq_def.py` (parser), `mech2asp.py` (rates and reactions), `mech_species.py` (paper-species mapping, molar masses),
`mech_photo.py` (photolysis columns) and `build_mech_decks.py` (deck writer).
"CRACCM" in the request was read as CRACMM. All runs are 300 min (baseline and CB6r5 decks changed too).

**Needs common_asp >= 40c1760** (new 1-body flag 13 and 3-body flag 4, explicit Troe N, optional equilibrium constant).

## Reactions
Core gas phase only: CRACMM groups R, TRP, ROCARO, RAM, RBAL, T, HC; SAPRC groups numeric, BR, R, BP, P, PX, PZ, BE, BT, BL, IS, IA, IC.
Not translated: chlorine, HET_*/HAL_*, aerosol/SVOC (AE), primary tracers, OLIG_*, POA/RPOA ageing, PCSOA, HYD_*, ROCALK/ROCOXY.
Reactions: CRACMM2 430, CRACMM3 441, SAPRC07TC 587, SAPRC07TIC 755 (see each deck header and "NOT TRANSLATED" list).

Every rate form was recognised structurally and checked numerically against the CMAQ expression at 7 (T,P) points
(`python3 mech2asp.py <MECH>`). All reactions agree to <1e-6 relative except the documented approximations:
* O3P+O2+M: ASP hard-wired flag 9 (up to 7% for CRACMM, 0.3% SAPRC; O(3P) goes to O3 either way).
* O1D quenching: ASP flag 8 (1-2%); CRACMM's O1D+O2 and O1D+N2 are merged into one reaction.
Exact mappings: fall-off with any F,N (flag 3/flag-less/flag 4), unimolecular fall-off and equilibrium reverse reactions
(`*E<ref>`: N2O5, HNO4, PAN, PPN) via flag 13, CRACMM R051 (T11 form) via 2-body flag 11, `%2` HNO3+OH via flag 5,
SAPRC `*K<ref>` operators expanded statically.
CMAQ data typo: SAPRC07TIC IS24 rate `3.94-15@1520` is read as 3.94e-15 (same for CB6 CL13, not translated).

## SVOC / SOA bookkeeping
* CRACMM reactions whose reactant is a V* (SVOC) or aerosol species (21) are dropped and V*/A*J products are not carried; ASP keeps its own IVOC block. This removes CRACMM's HOM/terpene-nitrate sink and its VROC chains.
* DELTA_C/N/SI and SAPRC XC/XN placeholders are dropped. CMAQ's ECH4+HO duplicate (R364) is skipped.
* Counters: CRACMM mirrors the baseline (HC5, HC10->HC8SOA, OLT, OLI, BDE13->OLISOA, TOL, XYL, CSL, API, LIM, ISO attach to the first oxidation step; all eight IVOC chains kept).
  SAPRC re-keys CMAQ's TOLRO2, XYLRO2, ISOPRXN, TRPRXN to ASP's TOLSOA, XYLSOA, ISOSOA, TERPSOA, adds TERPSOA to TERP+NO3 (TIC), and drops BENZRO2, PAHRO2, SULRXN, SESQRXN. HC5/HC8/OLT/OLI chains are not present (as in CB6r5).

## Species mapping (SPECIES_MAPPING.csv in each deck)
Mole-for-mole onto the lump CMAQ assigns by kOH class (SAPRC: ALK1..5, OLE1/2, ARO1/ARO2MN, PRD2) or structure (CRACMM: HC3, HC5, HC10, OLT, OLI, KET...).
Carbon is not conserved. Judgement calls: heptane -> ALK4/HC5 (boundary); furan -> OLE1 (SAPRC), FURAN (CRACMM); terpenes -> TERP (SAPRC), 0.5 API+0.5 LIM (CRACMM);
hydroxyacetone -> MEK (TC), HACET (TIC), HKET (CRACMM); glycolaldehyde -> HKET (CRACMM); biacetyl -> 1.33 MGLY (CRACMM);
RACM_HC8 -> ALK5/HC10; xylenes split equally into MXYL/OXYL/PXYL.

## Photolysis (PhotoRates.in, 32-35 columns)
Same-process TUV column x CMAQ scale where one exists. For labels without a TUV column the rate is J(reference TUV column) x band-weighted
CSQY ratio from CMAQ's CSQY_DATA, with band attenuations fitted to the analog labels at each time (`mech_photo.py`):
CH3COCH3B (vs A), UALD (vs acrolein), GLYH2 (vs GLYHX), BALD (vs MEK), AFG1 (vs NO2); HPALD uses column 51. These are proxies, uncertain by tens of percent.

## Results (5 h, dO3/dCO and dOA/dCO)
See `run_summary_5h.csv` and `enhancement_ratios_6mech.png`. Runs take 10-45 s each. Earlier run_summary.csv predates the flag-5 fix and the 5 h length.
Open: Timbavati OA about half of observed, SO2 background mismatch, Williams IVOC3-7 initial values, HET_/NTR2 not translated.
