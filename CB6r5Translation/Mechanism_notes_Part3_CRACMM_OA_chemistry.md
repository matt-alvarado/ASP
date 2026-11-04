# CRACMM OA chemistry: ASP baseline vs CMAQ

Oct 30, 2026 · @Matt

## Summary

The ASP baseline treats organic aerosol as one volatility ladder (IVOC1–9) fed by fixed-yield SOA counters and a symmetric OH aging step. CMAQ CRACMM instead carries a two-dimensional set of species (volatility × chemical class, with O:C), forms SOA through explicit multigeneration gas chemistry, keeps terpene HOM and nitrate species, and puts everything in one ideal absorbing phase. The CRACMM2/3 decks now use the CMAQ species, volatilities (C\*, molar mass, ΔHvap) and OH/NO/HO2 chemistry for all gas-phase OA formation. Particle-phase hydrolysis and photolysis, and isoprene IEPOX uptake, cannot be represented in ASP and are left out. Williams starts with the same 903 µg/m³ of OA as before, but it evaporates less: ΔOA/ΔCO at 5 h goes from about 110 to 393 (CRACMM2) and 294 (CRACMM3) µg/m³ per ppmv. The other three cases change by less than 8%.

## ASP baseline (and the CRACMM2/3 decks as they stand)

The baseline carries all organic aerosol on a single one-dimensional volatility ladder, IVOC1 to IVOC9, and the existing CRACMM2/3 decks inherit it unchanged.

- **Volatility:** nine bins with C\* = 10^(k-3) µg/m³ at 300 K (IVOC1 = 0.01, IVOC3 = 1, IVOC7 = 10,000, IVOC9 = 10^6). Each bin has a fixed molar mass (524 down to 163 g/mol) and a fixed enthalpy of vaporization (77 down to 46 kJ/mol). Oxidation state is not tracked.
- **Partitioning:** each IVOCk is one species that exists in both the gas and the hydrophobic organic particle phase. Gas and particle exchange as a condensation flux with mass accommodation 0.1, and the particle side obeys Raoult's law within the organic phase.
- **POA:** the Williams deck starts with 794 µg/m³ of particulate IVOC3–6 and 0.02–467 ppb of gas IVOC3–7 (the Akagi-style unidentified NMOG pool). The older POA1–8, LEVO, CBIO and CPD1–3 species exist in the files but are unused in these cases.
- **SOA formation:** the first oxidation step of HC5, HC8, OLT, OLI, TOL, XYL, CSL, ISO and API/LIM also produces a dummy counter (HC5SOA … TERPSOA). The counter reacts with NO, NO3 or HO2, or unimolecularly, and releases fixed yields into IVOC3–6. The NO-versus-HO2 branching sets the high-NOx and low-NOx yields.
- **Aging:** gas-phase IVOCk + OH (1.0e-11 cm³/s) gives 0.5375 IVOC(k-1) + 0.5375 IVOC(k+1) + 0.6 HO2 + 0.5 ORA2. The mass is conserved at 1.075 per OH hit, so material moves to both higher and lower volatility. Particle-phase material does not age.
- **Not represented:** nonvolatile SOA, highly oxygenated molecules, terpene nitrates and the O:C dependence of volatility. The existing CRACMM decks also drop CMAQ's V-species reactions and the ASOATJ products (21 reactions).

## CMAQ CRACMM2 and CRACMM3

CMAQ carries organics as a two-dimensional set of named species (volatility and chemical class), forms SOA by explicit gas-phase multigeneration chemistry, and partitions everything into one absorbing organic phase at equilibrium. The statements below come from the CMAQ source in `CCTM/src/aero/cracmm/ORG_DEFN.F`, the `AE_`/`GC_` namelists and `mech_cracmm2.def` / `mech_cracmm3.def`.

- **Species and volatility:** every V (vapor) species is paired with an A (particle) species. The name encodes the volatility, so VROCN2 is C\* = 10^-2 µg/m³ and VROCP6 is C\* = 10^6 µg/m³. Three classes exist: alkane-like (ALK, nine bins from 10^-2 to 10^6), aromatic (ARO, 10^5 and 10^6) and oxygenated (OXY, 15 species spanning C\* 10^-2 to 10^6 and O:C 0.1 to 0.8, with a number suffix for the O:C class). Only the bins up to C\* = 10^3 have a particle partner; the rest are gas-only. Each species has its own molar mass (for example 422.8 for ROCN2ALK, 158.2 for ROCP4OXY2) and enthalpy of vaporization (104 to 73 kJ/mol for ALK, 93 to 73 for OXY).
- **Terpene products:** VHOM (C\* = 6.3e-3, 130 kJ/mol), VELHOM (1e-5), VTRPN terpene nitrate (1,410), VHONIT (0.041) and, in CRACMM3 only, VMTN1. ASOAT, AISO3NOS, AISO3OS, AISO4 and AISO5 are nonvolatile (C\* = 1e-10).
- **Partitioning:** all organic mass, including nonvolatile species, is one ideal Raoult's-law solution with each species' own molar mass. C\* follows C\*(T) = C\*(298) (298/T) exp\[ΔH/R (1/298 − 1/T)\]. The solver finds the equilibrium total absorbing moles each step (Newton/bisection), assuming instantaneous equilibrium (the decks use accommodation 1 as the closest ASP equivalent) and no oligomerization (removed in 2018).
- **POA:** by default POC and PNCOM are emitted as nonvolatile APOC and APNCOM. The namelist offers optional semivolatile families (VLVPO1 to VIVPO1 with C\* 0.1 to 1,000) or the full-volatility profile (VROC\*/AROC\* species). CRACMM has no particle-phase POA aging reaction.
- **SOA from IVOCs:** ALK species react with OH (1.5 to 3.4e-11 cm³/s) to a peroxy radical. With NO this gives 72% of a second RO2 and 28% of an OXY product, and the second RO2 isomerizes at 0.188 /s (autoxidation). Aromatic IVOCs follow an MCM-style bicyclic RO2 route that also gives dicarbonyls.
- **Aging of OXY species:** each OXY species + OH (2 to 6e-11) goes to a lower-volatility OXY species, to lower-volatility OXY species through fragmentation products (MEK, KET, HCHO, ACD, ALD and RO2s), or back into the same bin. It is a multigeneration network with explicit fragmentation, not a volatility ladder.
- **SOA from VOCs:** BEN, TOL, XYL and EBZ peroxy radicals + NO or HO2 form OXY species at 0.1% to 28% yield. PHEN and CSL + OH give 15% to 20% nonvolatile ASOAT, and HC3/HC5 give 0.003% to 0.13%.
- **Terpenes:** API/LIM peroxy radicals form VTRPN, and autoxidation, RO2+RO2 and RO2+HO2 give VHOM. VTRPN + OH, O3 or NO3 gives VHONIT. In the particle, ATRPN and AHONIT hydrolyze to AHOM + HNO3 at 9.26e-5 /s, and in CRACMM3 AHOM photolyzes at 1% of jNO2 to AMTN1. Isoprene SOA (IEPOX/IPX/INALD uptake) goes to nonvolatile AISO species through reactive uptake that depends on aerosol surface area and acidity.

## Side-by-side

| Topic | ASP baseline | CMAQ CRACMM2/3 |
| --- | --- | --- |
| Volatility axis | IVOC1–9, C\* = 10^(k−3) µg/m³ at 300 K | V/A pairs, C\* 10^-2 to 10^6 at 298 K, plus HOM, nitrates and nonvolatile species |
| Other axes | None; MW and ΔHvap tied to bin | Class (ALK, ARO, OXY) and O:C, each species with own MW and ΔHvap |
| Gas/particle | One species per bin, accommodation 0.1 | One species per V/A pair, equilibrium (accommodation 1) |
| Absorbing phase | Hydrophobic organic phase, Raoult | Single ideal Raoult solution, nonvolatile species included |
| POA | Measured IVOC distribution (Williams); POA1–8 unused | Default nonvolatile APOC/APNCOM; optional semivolatile profiles |
| SOA formation | Dummy counters with fixed yields into IVOC3–6 | Explicit RO2+NO/HO2, autoxidation, ASOAT yields from phenols/cresols |
| Aging | IVOCk+OH → IVOC(k−1) + IVOC(k+1), symmetric | OXY multigeneration with fragmentation; ALK to OXY |
| Terpenes | Lumped into counter | VTRPN, VHOM, VHONIT (and VMTN1 in CRACMM3) |
| Particle-phase chemistry | None | ATRPN/AHONIT hydrolysis; AHOM photolysis in CRACMM3 |
| Isoprene SOA | None | Reactive IEPOX/IPX/INALD uptake to nonvolatile AISO species |

## Changes to the CRACMM2/3 decks

- **Chemistry:** all CMAQ gas-phase reactions of the V species, ROCALK/ROCOXY/ROCARO and terpene groups are now translated (515 reactions for CRACMM2, 525 for CRACMM3). The ASP counters, counter chains and IVOC aging reactions are removed. Particle-phase reactions (ATRPNJ, AHONITJ, AHOMJ, SULRXN) are dropped; ASOATJ is kept as nonvolatile ASOAT.
- **Species and volatility:** each partitioning V species uses its CMAQ C\*(298), ΔHvap and molar mass, converted to Psat(300 K) for ASP's organic dissolution (flag 1), mass accommodation 1. ASOAT has C\* = 1e-10 µg/m³.
- **Initial OA, Williams:** each IVOC bin (gas and particle) moves to the CMAQ oxygenated species with the same C\* and the lowest O:C, with mass conserved (for example particle IVOC3–6 → VROCP0OXY2, P1OXY1, P2OXY2, P3OXY2).
- **Initial OA, other cases:** the POA1 and CBIO particle mass is spread with CMAQ's documented GENERAL\_POA volatility profile, mass conserved.
- **Files:** OrgPhaseChems, OrganicDissolution, AerosolModes and the hydrophilic organic lists are regenerated. The hydrophilic lists carry the CMAQ species with a negligible Henry's constant because ASP's output routines index them in parallel with the organic list. A CMAQ\_OA\_MAPPING.txt in each deck records the mapping.
- **Results at 5 h, ΔOA/ΔCO (µg/m³ STP per ppmv), old → new:** Williams 110.6 → 393.2 (CRACMM2) and 109.4 → 293.7 (CRACMM3); Otavi 154 → 164 / 167; Alaska 38.9 → 36.4 / 36.7; Timbavati 46.8 → 44.0 / 44.6. ΔO3/ΔCO changes by at most 0.03 (Williams CRACMM3, 0.125 → 0.098).

## Assumptions and caveats

- Not translated: particle-phase hydrolysis and photolysis, IEPOX and related heterogeneous uptake (isoprene SOA), and aqueous uptake of organics. Isoprene and terpene-nitrate SOA is therefore underestimated relative to CMAQ.
- For non-Williams cases CMAQ's default would emit nonvolatile APOC; the decks use the semivolatile GENERAL\_POA profile instead. Switching is a small change in cmaq\_oa.py.
- The baseline IVOC distribution is mapped onto the lowest-O:C OXY species of the same C\*. A higher-O:C choice would change MW and ΔHvap slightly.
- ASP's TotPOA/TotSOA output columns assume the baseline species order, so for these decks use the SizeDist TotOA column for total OA.
- Williams OA evaporates less than in the baseline, partly because the CMAQ species have a higher ΔHvap and accommodation 1; the size of each effect was not isolated in separate runs.
- The earlier decks were overwritten in place; the previous versions are in git history.

## Part 2: CB6r5, SAPRC07TC and SAPRC07TIC (CMAQ AERO6/7)

These three mechanisms use CMAQ's older aerosol modules (AERO6 for SAPRC07TC, AERO7 for CB6r5, AERO7i for SAPRC07TIC). Source: `CCTM/src/aero/aero{6,7}/ORG_DEFN.F`, the `mech_*.def` files, the AE/GC namelists and `CMAQ_Chem_Control_*.nml`.

### How CMAQ differs from the baseline

- **Volatility:** semivolatile vapors have their own C\*(298 K), enthalpy of vaporization, molar mass and O:C, so there is no single ladder. Examples: SVAVB1–4 (aromatic/PAH/alkane products, C\* 0.01 to 100), SVMT1–6 (monoterpenes, 0.01 to 1,000), SVISO1/2, SVSQT, MTNO3 (12), ISOPNN (8.9), and POA species VLVPO1 to VIVPO1 (0.1 to 1,000).
- **SOA formation:** either explicit gas reactions (for example TOLRO2 + NO → 0.016 SVAVB2 + 0.051 SVAVB3 + 0.047 SVAVB4) or reaction counters (ISOPRXN, TRPRXN, SESQRXN, ALKRXN, XYLNRXN …). The aerosol module converts a counter to vapor with a fixed mass yield α, vapor mass = α × counter × MW(counter). NO and HO2 branching enters through separate high-NOx and low-NOx reactions. Some low-NOx aromatic products (AXYL3 and similar) are nonvolatile.
- **Partitioning:** one ideal Raoult organic phase at equilibrium; no particle-phase aging except slow oligomerization (about 29 h) and 3 h hydrolysis of particle nitrates.
- **POA:** the default is nonvolatile APOC/APNCOM. The control namelist documents a semivolatile alternative: 0.09 as ALVPO1 and 0.045 as ASVPO1 particle, and 0.045 VSVPO1, 0.14 VSVPO2, 0.18 VSVPO3, 0.50 VIVPO1 as gas. CB6r5 adds ten POA-ageing OH reactions into VLVOO\*/VSVOO\*. pcSOA (PCVOC + OH → low-volatility SOA) is 6.579 × PM in CMAQ but is zeroed for fire sources.

### Deck changes

- **Chemistry:** the baseline counters (TOLSOA, XYLSOA, ISOSOA, TERPSOA), their IVOC chains and the IVOC aging block are removed. Added from CMAQ: the SA/AE yield reactions, SOAALK, terpene nitrate (MTNO3), pcSOA, and CB6r5 POA ageing. Counters are replaced by their products, with coefficient α × MW(counter) / MW(particle species), so the same vapor mass forms as in CMAQ. CB6r5 now has 252 reactions, SAPRC07TC 597 and SAPRC07TIC 765.
- **Species:** each vapor/particle pair is one ASP species with CMAQ's C\*, ΔHvap and molar mass, converted to Psat(300 K) (flag 1), mass accommodation 1.
- **Williams initial OA:** the measured IVOC3–6 (gas and particle) move to VSVPO1, VSVPO2, VSVPO3 and VIVPO1 (same C\*), mass conserved. Gas IVOC7 (467 ppb, about 4,350 µg/m³, C\* 10^4) has no CMAQ species and is not carried; the option to send it to PCVOC exists but is off.
- **Other cases:** the POA1 and CBIO particle mass is spread with the CMAQ semivolatile profile above, so 86% of it starts as vapor and re-equilibrates.

### Results at 5 h, ΔOA/ΔCO (µg/m³ STP per ppmv), baseline → CMAQ-style

| Case | CB6r5 | SAPRC07TC | SAPRC07TIC |
| --- | --- | --- | --- |
| Williams | 109 → 137 | 110 → 54 | 109 → 63 |
| Otavi | 153 → 142 | 153 → 118 | 154 → 119 |
| Alaska | 38 → 22 | 38 → 18 | 38 → 18 |
| Timbavati | 46 → 35 | 46 → 29 | 46 → 29 |

ΔO3/ΔCO changes by at most 0.013 (Williams, CB6r5 0.138 → 0.125).

### Caveats

- Not translated: oligomerization, nitrate hydrolysis, IEPOX and glyoxal heterogeneous uptake, aqueous organic uptake, MTNO3 photolysis in CB6r5, and all chlorine, DMS, mercury and toxics chemistry.
- CMAQ's default nonvolatile POA is not used; the semivolatile profile is. For fires CMAQ's pcSOA is off, so the IVOC mass above C\* = 10^3 does not form SOA here.
- The SAPRC07TIC species table in the CMAQ docs lists shifted molar masses for SVAVB1–4; the AE namelist values are used.

## Part 3: CMAQ heterogeneous uptake and particle-phase reactions

All five mechanism decks now carry CMAQ's uptake and first-order particle-phase organic reactions, and the new pathways change OA beyond the earlier tier 1+2 step only in the CRACMM decks.

### What was added (common\_asp asp\_v2.2.2, one commit per patch)

- **Flag 3** (constant γ with molecular weight): glyoxal (γ 2.9e-3) and methylglyoxal (0.09 × γ) to nonvolatile AGLY; NO3 (γ 1e-3, SAPRC07TIC and CRACMM); HO2 (γ 0.2, CRACMM).
- **Flag 18**: acid-dependent epoxide uptake (CMAQ IEPOX/IMAE kinetics from H+, HSO4-, SO4--, NO3-, with the CMAQ 1/γ resistance form), on the total wet-particle volume basis as CMAQ.
- **Flag 19**: NTR2 hydrolysis (CB6r5), particle fraction from OA and an RH ramp of 0.2–0.4.
- **Flag 20**: N2O5 uptake with the Davis et al. 2008 γ (CMAQ default), products 2 HNO3 in every mechanism.
- **Flag 21**: NO2 uptake, k = 1e-4 m/s × wet surface area (bins below 1.25 µm), products 0.5 HONO + 0.5 HNO3, in all five mechanism decks (CMAQ defines HET\_NO2 in the CB6r5, CRACMM2/3, SAPRC07TC and SAPRC07TIC mechanism files). Logged k is about 1.7e-6 /s in Williams.
- **OrgPhaseReactions.in** (new routine `StepOrgPhaseReactions`): oligomerisation (9.48816e-6 /s), AMTNO3/AISOPNN hydrolysis (9.259e-5 /s), CRACMM ATRPN/AHONIT to AHOM + HNO3, and CRACMM3 AHOM photolysis-like loss at 0.01 × jNO2.

Decks are generated by `cmaq_het.py` (block between `CMAQ_HET_BEGIN` and `CMAQ_HET_END`, idempotent, flags work at any T and P).

### OA at 290 min (sum of TotOA, µg/m³)

| Case | Mechanism | No uptake | Tiers 1+2 | + volume basis, NTR2, N2O5, NO3, HO2 | + NO2 uptake |
| --- | --- | --- | --- | --- | --- |
| Williams | CRACMM2 | 184.5 | 206.2 | 193.1 | 193.0 |
| Williams | CRACMM3 | 136.1 | 136.6 | 129.2 | 129.4 |
| Williams | CB6r5 | 62.7 | 77.7 | 77.7 | 77.7 |
| Williams | SAPRC07TC | 24.5 | 25.2 | 25.2 | 25.2 |
| Williams | SAPRC07TIC | 28.9 | 52.9 | 52.9 | 52.9 |

Otavi, Alaska and Timbavati change by under 0.2 µg/m³ with the last two steps. The comparison plots (`cmaq_oa_compare.png`, `cmaq_ae_compare.png`, `all_mechs_compare.png`) are regenerated in the asp repo.

### Open items from earlier steps

- **Not verified:** the NTR2 rate was seen in the log only at t=0 (0); check later in the run.
- **Not isolated:** the SAPRC07TIC OA gain (+83 % in Williams) is not from epoxide uptake or the volume basis; the share of uptake vs particle-phase reactions was not split.
- **Not represented:** IEPOXP/sulfate consumption (epoxide goes straight to product), AERO7i dimers (folded into the tetrol and organosulfate channels), CMAQ's H2NO3P/NO2P tracer products of N2O5 uptake, and the coarse-mode N2O5 γ.
- **Unresolved from the base mechanism work:** Timbavati OA is about half of observed, SO2 background mismatch, Williams IVOC3–7 initial values.
- **Held by request:** SurfaceTension simplification; thermodynamics NaN warnings stay uncommented.
- **Standing constraints:** master untouched (you merge); `src` submodule pin moved to the `common_asp` `asp_v2.2.2` tip whenever that branch gets a commit (e3faf76 as of patch 37).

## Part 4: GEOS-Chem 14.8.0 and MOZART-T1 OA

The GEOS-Chem (`_GCHEM`) and MOZART-T1 (`_MOZT1`) decks for all four cases now carry each model's own OA scheme instead of the ASP IVOC ladder. The OA chemistry is added after the gas mechanism is built, by `MozGCTranslation/moz_oa.py` (asp commit 6841b37). Sources: GEOS-Chem 14.8.0 (geoschem/geos-chem), WRF-Chem KPP t1\_mozcart (tag MZ148\_T1mozcart\_20170323), and Tilmes et al. (2019) for the CESM2 VBS-SOA. MozGCTranslation/README.md documents the same treatment.

### GEOS-Chem: complex SOA with semivolatile POA

- **Species:** TSOA0–3 (terpene SOA), ASOA1–3 and ASOAN (aromatic and IVOC SOA), SVPOA1/2 (semivolatile POA) and OPOA1/2 (oxidized POA).
- **SOA yields:** aromatic and IVOC SOA follow GEOS-Chem's BETANO split between the RO2 + NO and RO2 + HO2 channels.
- **IVOC:** the Williams gas IVOC pool (IVOC7, 467 ppb) is carried as NAP, GEOS-Chem's IVOC surrogate, at 925 ppb after the ARO2CARB conversion, and forms ASOA with GEOS-Chem's yields.
- **Heterogeneous products:** SOAGX (glyoxal), SOAIE (IEPOX), LVOCOA, IONITA, MONITA and AONITA. IONITA and MONITA hydrolyze at 2.78e-4 /s (about a 1 h lifetime). Uptake uses fixed γ.
- **Molar mass:** every OA species uses MW 150, so ASP's molar Raoult's law gives the same partitioning as GEOS-Chem's mass-based absorption.

### MOZART-T1: CESM2 VBS-SOA

- **Species and volatility:** SOAG0–4 with C\* = 0.01 to 100 µg/m³ at 300 K.
- **Precursors:** SVOC = 0.6 × POM and IVOC = 0.2 × NMVOC, following the CESM2 scaling. Glyoxal goes to SOAG0.
- **Loss:** SOA photolysis at J(SOA) = 0.04 × J(NO2).
- **POA:** POM is nonvolatile.
- **Assumed constants:** MW 150 and ΔHvap 30 kJ/mol for SOAG, IVOC MW 184, SVOC/POM MW 250. These are not yet checked against CAM-chem `chem_mech.in`.

### Generation

- `python3 moz_oa.py <deck> <case> gc|t1` writes the OA block between `!! MOZ_OA_BEGIN` and `!! MOZ_OA_END`. `build_gc_decks.py` and `build_t1_decks.py` call it.
- Excluded at your request: all Cl-radical and Cl + VOC chemistry.

### Results

All eight decks run to completion (rc = 0, 9–41 s wall time on one core, common\_asp c506969, asp 6841b37 decks). In Williams, GEOS-Chem OA rises from 45 to 170 µg/m³ once the IVOC pool is carried as NAP, and is dominated by ASOA2, ASOA3, SVPOA2 and OPOA1. MOZART-T1 OA in Williams is dominated by POM (87 %) and SOAG0 (8 %).

| Case | OA at 290 min, GCHEM / MOZT1 (µg/m³) | ΔOA/ΔCO at 5 h, GCHEM / MOZT1 | ΔO3/ΔCO at 5 h, GCHEM / MOZT1 |
| --- | --- | --- | --- |
| Williams | 170.3 / 58.3 | 376.6 / 127.3 | 0.074 / 0.103 |
| Otavi | 7.9 / 11.2 | 123.1 / 173.7 | 0.105 / 0.113 |
| Alaska | 33.4 / 29.3 | 47.2 / 41.5 | 0.065 / 0.053 |
| Timbavati | 12.0 / 11.1 | 46.9 / 43.8 | 0.183 / 0.166 |

ΔOA/ΔCO is in µg/m³ STP per ppmv, computed as in Parts 1–3. GEOS-Chem in Williams grows to about five times the observed ΔOA/ΔCO (about 78 at 4.2 h), while MOZART-T1 is about 1.6 times observed and decreases slowly with time. In the other three cases the two mechanisms agree to within 15 % except Otavi, where MOZART-T1 is 40 % higher. The comparison plot is `MozGCTranslation/gchem_mozt1_oa_compare.png`.

### All nine mechanisms

`CB6r5Translation/all_mechs_oa_o3_compare.png` shows ΔOA/ΔCO and ΔO3/ΔCO against observations for Baseline, CRACMM2, CRACMM3, CB6r5, SAPRC07TC, SAPRC07TIC, GCHEM, MOZT1 and MCM331, all rerun with the same binary. Baseline and MCM331 use the baseline IVOC-ladder OA, the five CMAQ-derived decks use CMAQ OA (Parts 1–3), and GCHEM and MOZT1 use their own OA (above). The baseline and CMAQ-derived results match the committed `all_mechs_compare.png`.

**ΔOA/ΔCO at 5 h (µg/m³ STP per ppmv)**

| Case | Baseline | CRACMM2 | CRACMM3 | CB6r5 | SAPRC07TC | SAPRC07TIC | GCHEM | MOZT1 | MCM331 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Williams | 106 | 419 | 285 | 172 | 55 | 118 | 377 | 127 | 108 |
| Otavi | 151 | 163 | 165 | 144 | 118 | 120 | 123 | 174 | 154 |
| Alaska | 38 | 41 | 41 | 26 | 18 | 23 | 47 | 42 | 38 |
| Timbavati | 46 | 47 | 48 | 39 | 29 | 31 | 47 | 44 | 46 |

**ΔO3/ΔCO at 5 h (ppb/ppb)**

| Case | Baseline | CRACMM2 | CRACMM3 | CB6r5 | SAPRC07TC | SAPRC07TIC | GCHEM | MOZT1 | MCM331 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Williams | 0.128 | 0.133 | 0.061 | 0.100 | 0.127 | 0.109 | 0.074 | 0.103 | 0.144 |
| Otavi | 0.142 | 0.098 | 0.095 | 0.141 | 0.143 | 0.153 | 0.105 | 0.113 | 0.154 |
| Alaska | 0.050 | 0.048 | 0.047 | 0.056 | 0.056 | 0.065 | 0.065 | 0.053 | 0.070 |
| Timbavati | 0.169 | 0.155 | 0.152 | 0.186 | 0.202 | 0.209 | 0.183 | 0.166 | 0.206 |

- **Williams OA:** CRACMM2, GCHEM and CRACMM3 are far above the observed ΔOA/ΔCO of about 78. MOZT1 sits with Baseline, MCM331 and SAPRC07TIC; SAPRC07TC is closest to observed.
- **OA elsewhere:** GCHEM and MOZT1 fall within the spread of the other mechanisms; GCHEM is highest in Alaska. MCM331 follows Baseline closely, since both use the same OA scheme.
- **O3:** MCM331 gives the highest or near-highest ΔO3/ΔCO in every case. GCHEM and MOZT1 are mid-range.

### Open questions

1. **T1 constants:** the SOAG, IVOC and SVOC/POM molar masses and ΔHvap above are assumed and need CAM-chem `chem_mech.in` to confirm. A MW sensitivity run (150 vs 12.011) was not done.
2. **Williams GEOS-Chem OA:** the 45 → 170 µg/m³ jump comes from 925 ppb of NAP times the ASOA yields, and puts ΔOA/ΔCO at about five times observed. Its plausibility is not confirmed.
3. **Approximations:** NAP as the IVOC surrogate (ARO2CARB conversion) and the aromatic RO2 yield mapping.
4. **Aqueous sulfur:** should ASP's native S(IV) aqueous chemistry (disabled in all decks) be enabled for GEOS-Chem SO2 uptake and HMS?
5. **POA scheme:** confirm the complex SVPOA scheme rather than GEOS-Chem's simple SOA (SOAP/SOAS).
6. **Not replicated:** GEOS-Chem's composition- and RH-dependent γ (fixed γ is used), sea-salt reactions, Br/I heterogeneous chemistry (needs halogen gas chemistry and sea-salt tracers) and cloud/aqueous sulfur. T1 has no aqueous KPP set.
7. **Carried over:** CMAQ's H2NO3P/NO2P N2O5 tracer products remain flagged (Part 3), and earlier items are in `NewMechanisms_README.md` ("Open items carried over").
