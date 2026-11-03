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

## CRACMM2/3 organic aerosol now follows CMAQ (2026-10-30)

The CRACMM2 and CRACMM3 decks no longer use the baseline IVOC1-9 ladder, SOA counters or IVOC aging. They use CMAQ's V/A volatility species
(C*, MW, dHvap from `ORG_DEFN.F`), the CMAQ multigeneration gas chemistry of those species (ROCALK/ROCOXY/ROCARO, terpene HOM/nitrates) and
nonvolatile ASOAT (`cmaq_oa.py`, `mech2asp.py oa='cmaq'`).  Williams' IVOC distribution is mapped bin-by-bin to the lowest-O:C OXY species of the
same C*; other cases use CMAQ's GENERAL_POA profile.  Not translated: particle-phase hydrolysis/photolysis, IEPOX uptake, aqueous organic uptake.
Each deck has `CMAQ_OA_MAPPING.txt`.  Note: `TotPOA`/`TotSOA` in Aerosol*.txt assume the baseline species order; use `TotOA` in SizeDist*.txt.
`compare_cmaq_oa.py` makes `cmaq_oa_compare.png` (old vs new).  Rebuild: `python3 build_mech_decks.py CRACMM2 CRACMM3`.

## CB6r5, SAPRC07TC and SAPRC07TIC organic aerosol now follow CMAQ AERO6/7 (2026-10-30)

These 12 decks no longer use the baseline IVOC ladder, SOA counters or IVOC aging.  `cmaq_ae.py` reads the CMAQ `ORG_DEFN.F` (aero6 / aero7), AE and GC
namelists: SOA counters (ISOPRXN, TRPRXN, SESQRXN, ALKRXN, XYLNRXN ...) are replaced by their mass-yield products (alpha x MW(counter)/MW(particle)),
the AE-group yield reactions (SA*, AE5x, AALK), terpene nitrate, pcSOA and (CB6r5) POA ageing are translated, and every vapor/particle pair is one ASP
species with CMAQ's C*, dHvap and MW.  Williams IVOC3-6 map to VSVPO1-3 / VIVPO1 by C*; IVOC7-9 (C* >= 1e4) have no CMAQ species and are dropped.
Other cases use the CMAQ semivolatile POA profile.  Not translated: oligomerisation, nitrate hydrolysis, IEPOX/glyoxal uptake, MTNO3 photolysis (CB6r5).
`cmaq_ae_compare.png` (old vs new): `python3 compare_cmaq_oa.py <runs> cmaq_ae_compare.png CB6r5,SAPRC07TC,SAPRC07TIC`.
Rebuild: `python3 build_cb6r5_decks.py` (CB6r5) and `python3 build_mech_decks.py SAPRC07TC SAPRC07TIC`.  `OA_MODE = 'asp'` in build_cb6r5_decks.py
restores the old behaviour for CB6r5.

## Gas-to-particle uptake and particle-phase organic reactions (2026-10-30)

`cmaq_het.py` edits the 20 finished decks in place (run it after the deck builders; idempotent) and needs common_asp >= 5880438 (branch asp_v2.2.2).
* **Tier 1, uptake** (lines between `CMAQ_HET_BEGIN/END` in GasChemicalMechanism.in; product species are nonvolatile, like ASOAT):
  GLY, MGLY -> AGLY (AGLYOLIG in CRACMM3), flag 3, gamma 2.9e-3 and 0.09 x 2.9e-3 as CMAQ (not SAPRC07TC, which has no GLY uptake).
  IEPOX (EPOX in CB6r5), IMAE, IHMML, IPX (x2), INALD (x0.5, + HNO3) use the new **1-body rate flag 18**:
  `EPOX => AISO3 ; 18 ; scale ; MW ; kind ; channel ;`.  Per particle, the ring-opening rate sum k[acid][nucleophile] is built from H+, HSO4-, SO4--, NO3-
  and water of that particle (CMAQ Eddingsaas/Piletic constants, concentrations per litre of aqueous solution, CMAQ uses total particle volume), gamma
  follows CMAQ CALCISOPGAMMAS (alpha 0.02, D 1e-9 m2/s, H 3e7 / 1.2e5 M/atm) and the loss rate is summed over bins (flag-3 formula).
  kind 1 IEPOX with nitrate, 2 IMAE, 3 IEPOX without nitrate; channel 0 all, 1 non-sulfate product, 2 organosulfate.
  `cmaq_het.apply(dir, mech, acid=False)` gives constant gamma (flag 3) instead.
* **Tier 2, `OrgPhaseReactions.in`** (new common_asp routine `StepOrgPhaseReactions`, run after organic condensation each condensation step, exact exponential
  per particle): oligomerisation to AOLGA/AOLGB (9.48816e-6 /s), AMTNO3 / AISOPNN hydrolysis (9.259e-5 /s, HNO3 to the gas), CRACMM ATRPN/AHONIT -> AHOM + HNO3
  (9.26e-5 /s), CRACMM3 AHOM -> 0.8 AMTN1 + 0.2 ORA1 at 0.01 x jNO2 (`PHOTO <gas rxn #> <scale>`).  Format is in `Condensation/OrgPhaseReactions.h`.
* Not represented: NTR2 hydrolysis, NO3/HO2/N2O5 uptake, IEPOXP/sulfate consumption (epoxide goes directly to product), AERO7i dimers (folded into the
  tetrol / organosulfate channels).  Gas and particle mass of a species are one ASP species, so converted mass re-equilibrates with the gas phase.

### Effect of tiers 1 + 2 on OA (sum of SizeDist TotOA at 290 min, with vs without the `CMAQ_HET` block and OrgPhaseReactions.in)
Williams: CRACMM2 184.5 -> 206.2, CRACMM3 136.1 -> 136.6, CB6r5 62.7 -> 77.7, SAPRC07TC 24.5 -> 25.2, SAPRC07TIC 28.9 -> 52.9.
Otavi: changes 0-1.2 %.  Alaska: CRACMM2/3 +11 %, CB6r5 +15 %, SAPRC07TC +0.2 %, SAPRC07TIC +30 %.  Timbavati: CRACMM +7-8 %, CB6r5 +11 %, SAPRC07TC +0.1 %, SAPRC07TIC +7 %.
SAPRC07TC barely changes: it has only IEPOX uptake and oligomerisation of already-nonvolatile-ish species.  The runs were not split into uptake-only
and reactions-only, so the share of each tier is not separated.

### Update: CMAQ volume basis for flag 18, NTR2 hydrolysis, N2O5 / NO3 / HO2 uptake (common_asp >= asp_v2.2.2 head)
* Flag 18 now divides by the total wet-particle volume, as CMAQ.  Flag 19 (`NTR2 => HNO3 ; 19 ; 1.4 ;`, CB6r5) is CMAQ HETERO_NTR2 (4/day x particle
  fraction from OA x RH ramp 0.2-0.4).  Flag 20 (`N2O5 => 2 HNO3 ; 20 ; 1 ; 108.01 ;`) uses CMAQ's default N2O5 gamma (Davis et al. 2008).
  NO3 (gamma 1e-3 -> HNO3; SAPRC07TIC, CRACMM) and HO2 (gamma 0.2, no products; CRACMM) use flag 3.  N2O5 goes to 2 HNO3 in every mechanism (CMAQ's
  H2NO3P / NO2P tracer products are not represented).  N2O5 gamma is applied to all bins (CMAQ: I+J modes; coarse mode has its own value).
* OA at 290 min (none / previous uptake+reactions / now): Williams CRACMM2 184.5 / 206.2 / 193.1, CRACMM3 136.1 / 136.6 / 129.2, CB6r5 62.7 / 77.7 / 77.7,
  SAPRC07TC 24.5 / 25.2 / 25.2, SAPRC07TIC 28.9 / 52.9 / 52.9.  CB6r5 and SAPRC decks are unchanged by the new pathways and the volume basis (to
  0.01 ug/m3 in Williams), so their OA gain comes from oligomerisation, hydrolysis and glyoxal uptake, not the epoxide uptake.  Only CRACMM changes
  (HO2 / NO3 / N2O5 uptake, volume basis); the split among them was not run.  Alaska and Timbavati CRACMM change by 0.1-0.2 ug/m3.

### Update: NO2 uptake to HONO (flag 21; common_asp asp_v2.2.2 head)
* `NO2 => 0.5 HONO + 0.5 HNO3 ; 1 ; 21 ; 1e-4 ;` (CMAQ HETERO_NO2, Vogel et al. 2003): k = v x wet surface area per volume of air, v = 1e-4 m/s, summed over bins
  with r < 1.25 um (CMAQ I+J modes).  It is in all five mechanism decks (CMAQ defines HET_NO2 in the CB6r5, CRACMM2/3, SAPRC07TC and SAPRC07TIC mechanism files; checked in
  `CCTM/src/MECHS`).  It is a single removable line in the `CMAQ_HET` block.  Logged k is about 1.7e-6 /s in Williams at t=0.
* OA at 290 min is unchanged to <= 0.2 ug/m3 in every case (Williams CRACMM2 193.0 -> 193.0, CRACMM3 129.2 -> 129.4, CB6r5 77.7, SAPRC07TC 25.2, SAPRC07TIC 52.9),
  as expected for a gas-phase NOx/HONO pathway.  Comparison plots regenerated.

### All nine mechanisms, 5 h (2026-11-02; common_asp c506969)
`all_mechs_oa_o3_compare.png`: Baseline, the five CMAQ-derived decks, GCHEM, MOZT1 and MCM331 against observations (`all_mechs_compare.png`, the six-mechanism
version, is unchanged and still matches). Baseline and MCM331 use the IVOC-ladder OA. Notes and open questions: `Mechanism_notes_Part3_CRACMM_OA_chemistry.md`.

**dOA/dCO at 5 h (ug/m3 STP per ppmv)**

| Case | Baseline | CRACMM2 | CRACMM3 | CB6r5 | SAPRC07TC | SAPRC07TIC | GCHEM | MOZT1 | MCM331 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Williams | 106 | 419 | 285 | 172 | 55 | 118 | 377 | 127 | 108 |
| Otavi | 151 | 163 | 165 | 144 | 118 | 120 | 123 | 174 | 154 |
| Alaska | 38 | 41 | 41 | 26 | 18 | 23 | 47 | 42 | 38 |
| Timbavati | 46 | 47 | 48 | 39 | 29 | 31 | 47 | 44 | 46 |

**dO3/dCO at 5 h (ppb/ppb)**

| Case | Baseline | CRACMM2 | CRACMM3 | CB6r5 | SAPRC07TC | SAPRC07TIC | GCHEM | MOZT1 | MCM331 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Williams | 0.128 | 0.133 | 0.061 | 0.100 | 0.127 | 0.109 | 0.074 | 0.103 | 0.144 |
| Otavi | 0.142 | 0.098 | 0.095 | 0.141 | 0.143 | 0.153 | 0.105 | 0.113 | 0.154 |
| Alaska | 0.050 | 0.048 | 0.047 | 0.056 | 0.056 | 0.065 | 0.065 | 0.053 | 0.070 |
| Timbavati | 0.169 | 0.155 | 0.152 | 0.186 | 0.202 | 0.209 | 0.183 | 0.166 | 0.206 |

### Open items carried over from earlier steps
* Timbavati OA about half of observed; SO2 background mismatch; Williams IVOC3-7 initial values.
* NTR2 hydrolysis rate was only seen in the log at t=0 (0); not verified later in the run.
* Not represented: IEPOXP/sulfate consumption (epoxide goes directly to product), AERO7i dimers (folded into tetrol / organosulfate channels), CMAQ's
  H2NO3P / NO2P tracer products of N2O5 uptake; coarse-mode N2O5 gamma.
* SAPRC07TIC OA gain (+83 % Williams) is not from epoxide uptake or the volume basis; its cause is not isolated.  Share of uptake vs particle-phase reactions not split.
* Held by request: SurfaceTension simplification; thermodynamics NaN warnings stay uncommented.
* `src` submodule pin follows the `common_asp` `asp_v2.2.2` tip: moved whenever that branch gets a commit (7ffe10e as of patch 31); master untouched (Matt merges).
