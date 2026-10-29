# CB6r5 gas-phase mechanism for the four ASP case decks

Decks: `InorganicInputDecks/Williams_CB6r5`, `Otavi_CB6r5`, `Alaska_CB6r5`, `Timbavati_CB6r5`.
Only the **gas-phase chemistry** was changed. Aerosol, aqueous, thermodynamic, dissolution and organic-phase decks,
`ASPInputDeck.in` and all run settings are byte-for-byte copies of the original case decks.

Mechanism source: CMAQ `CB6r5-AE7-AQ` (`CCTM/src/MECHS/cb6r5_ae7_aq/mech_cb6r5_ae7_aq.def`, the species table `mechanism_information/cb6r5_ae7_aq/cb6r5_ae7_aq_species_table.md`, and `GC_cb6r5_ae7_aq.nml` for molecular weights), https://github.com/USEPA/CMAQ.

## 1. What is in each new directory

| File | Status |
|---|---|
| `GasChemicalMechanism.in` | **replaced**: 224 CB6r5 reactions (each preceded by a `!! R<n>: <CMAQ statement>` comment and any WARNING), then the TOLSOA/XYLSOA/ISOSOA/TERPSOA -> IVOC chains and the IVOC aging block of the original case deck |
| `GasPhaseChems.in` | **replaced**: CB6r5 species; initial/background values from the papers' compound tables mapped to CB6r5 (section 4) |
| `PhotoRates.in` | **replaced**: CB6r5 column set derived from the case's TUV 5.4 table (section 5) |
| `SPECIES_MAPPING.csv` | **new**: every paper compound, its CB6r5 representation, carbon check and basis |
| all other `.in` files | unchanged |

Rebuild: `python3 CB6r5Translation/build_cb6r5_decks.py && python3 CB6r5Translation/write_readme.py`.

## 2. WARNINGS - CB6r5 rate forms ASP cannot (exactly) handle

Of the 224 core `R<n>` reactions, 224 are translated to within 0.12% of the CMAQ rate over 260-310 K, 600-1013 hPa (most to <0.01%) except for those listed below.

### 2.1 Cannot be represented - a fixed value is used instead

* **R46, OH + HNO3** - CB6r5 uses the form `k0 + k3[M]/(1 + k3[M]/k2)` (CMAQ `%2`). ASP's 2-body rate flag 5 is meant for this form but is **coded wrongly** (see 2.4), so it cannot be used. R46 is therefore a fixed 2-body rate evaluated at the case T and P (Williams 1.72e-13, Timbavati 1.73e-13 cm3/molec/s). It is exact only at those conditions (34% off at 20 K / 300 hPa away). All four decks use constant T and P, so this is harmless for these cases.

### 2.2 Approximated using ASP's hard-wired forms

* **R2, O + O2 + M** - ASP flag 9 (5.6e-34 for N2, 6.0e-34 for O2) vs CB6r5 6.0e-34 (T/300)^-2.6 for all M (<=5% difference).
* **R10, O(1D) + M** - ASP flag 8 (IUPAC N2/O2 form) vs CB6r5 2.23e-11 exp(115/T) (<=0.2%).
* **R24, NO + NO + O2** (termolecular, O2 as reactant) - ASP has no working 3-body Arrhenius form (see 2.4), so k[O2] is fitted as k300 (T/300)^-n x M (flag 1) with value and slope matched at the case temperature (<=1.3% over 260-310 K). Negligible at ambient NO.
* **R37, N2O5 decomposition** - exact: ASP flag 7 is the same IUPAC expression as CB6r5.
* Fall-off reactions in IUPAC form use ASP's built-in N(F) = 0.75 - 1.27 log10(Fc); CMAQ rounds n to two decimals, so rates differ by ~0.1%.

### 2.3 Expressible only after rewriting

* **34 `A*K<Rn>` rate references** (a rate defined as a multiple of another reaction's rate): ASP has no such syntax, so each is expanded statically to the referenced reaction's rate x factor. Exact, but the rates are no longer linked: changing one rate in the deck does not change its dependents. Affected: R60, R62, R63, R65, R66, R67, R74, R77, R78, R79, R80, R81, R82, R83, R84, R85, R86, R115, R153, R154, R168, R169, R177, R179, R182, R184, R188, R189, R207, R208, R209, R210, R211, R212.
* **Pseudo-reactants.** `O2`, `M`, `N2`, `H2O` as CMAQ reactants (R2, R10, R20, R24, R41, R221) are handled as in ASP (`M`, flag 10 for water, etc.); R20, R41 and R221 use ASP's C_H2O flag 10.
* **Reactions with no products** (e.g. R7 `O + O3 =`) get the sink product `M`, which ASP's parser ignores as a species.

### 2.4 Bugs found in the ASP source (not patched - awaiting your decision)

1. **2-body rate flag 5 is wrong.** It evaluates `k1 + k3*M*(1 + k3*M/k2)`, while the documented/IUPAC form is `k1 + k3*M/(1 + k3*M/k2)`. It is used by the OH + HNO3 line of the **existing baseline decks** (`HNO3 + OH => NO3 + H2O ; 5.; ...`), where it gives 4.1e-11 instead of 1.7e-13 cm3/molec/s at 288 K / 880 hPa (**~240x too fast**), recycling HNO3 to NO3. A/B test (baseline deck with the correct constant rate, otherwise identical): Williams end-of-run O3 114.2 -> 120.3 ppb (+5%), HNO3 3.1 -> 6.8 ppb; Timbavati O3 60.3 -> 61.5 ppb (+2%), HNO3 1.46 -> 1.69 ppb. Not fixed in the source, because it changes all baseline results.
2. **3-body rate flag 2 ("thermal third order") is unusable**: it reads its parameters from the wrong positions of the rate array, so the flag value itself is taken as the first rate coefficient. No existing deck uses it; the CB6r5 decks do not.

### 2.5 CB6r5 content not translated (outside the gas-chemistry-only scope or without ASP equivalent)

* Reaction groups of the .def file that are **not included**: BL (4), BP (2), CL (30), CP (1), DMS (4), HAL (1), HET (28), HG (5), HYD (1), OLIG (6), PCSOA (1), POA (10), R (1), RPOAGELI (1), RPOAGELJ (1), RPOAGEPI (1), RPOAGEPJ (1), SA (12), T (21), TCL (6).
  These are chlorine (CL), DMS, mercury (HG), halogen/ozone (HAL), heterogeneous/hydrolysis (HET), source-apportionment (SA), biogenic/biomass tracer (BL, BP), POA aging and oligomerization (POA, RPOAGE, OLIG, PCSOA), and the **`T*` toxics diagnostics** (no feedback on the gas mechanism). `BUTADIENE13` and `ACROLEIN` have reactions only in `T*`, so they are not used as species (butadiene -> 2 OLE; acrolein -> ISPD).
* **R216 (ECH4 + OH)** - the tracer-CH4 duplicate; omitted (CH4 itself has R<n> CH4 + OH).
* **NTR2 heterogeneous hydrolysis** (`HET_NTR2`): NTR2 is only a NOx reservoir after this translation (CMAQ removes it by hydrolysis to HNO3 in the aerosol module).
* **Heterogeneous N2O5, NO2, HO2 uptake** (`HET_*`): not translated. The original decks' heterogeneous reactions are untouched (they were already commented out in Timbavati).
* **`CH3NO3`, `RO2_T`, `RCO3_T`** are hard-wired species names in ASP's output code; they exist as zero placeholders (CB6r5 has no CH3NO3: it is lumped in NTR1; `RO2_T`/`RCO3_T` are sums computed by ASP from the peroxy flags).

## 3. Reaction-by-reaction check

Every translated reaction was evaluated with a port of ASP's rate code and compared with CMAQ's own rate expression at 7 (T, P) points covering 260-310 K, 600-1013 hPa (`cmaq2asp.py --compare`). Largest deviations:

| Reaction | max rel. deviation | note |
|---|---|---|
| R46 | 0.34 | CB6r5 form k0 + k3[M]/(1+k3[M]/k2) cannot be represented: ASP rate flag 5 evalua |
| R2 | 0.05 | ASP hard-wires O + O2 + M (flag 9): 5.6e-34 N2 / 6.0e-34 O2 instead of CB6r5 6.0 |
| R24 | 0.013 | termolecular NO+NO+O2: ASP lacks a working 3-body Arrhenius form; k*[O2] fitted  |
| R10 | 0.0018 | ASP hard-wires O1D + M (flag 8, IUPAC N2/O2 form) instead of CB6r5 2.23e-11 exp( |
| R63 | 0.0011 |  |
| R209 | 0.0011 |  |
| R55 | 0.0011 |  |
| R62 | 0.0011 |  |
| R208 | 0.0011 |  |
| R54 | 0.0011 |  |

All other 185 reactions agree to < 0.1% (n rounding of IUPAC fall-off forms). Photolysis reactions are checked through the column map in section 5.

## 4. Species lumping and conversion

Procedure: the initial and background values come from the **papers' original compound tables** (Alvarado & Prinn 2009, JGR, auxiliary tables S1-S3 for Otavi, Alaska, Timbavati; the ACP 2015 supplement deck for Williams, whose nonzero values are identical to the old Williams deck), **not** from the old RACM2-style lumps. Each compound is expanded to CB6r5 species, and the contributions are summed per species (initial and background separately).

### 4.1 CB6 lumping rules used

Per the CB6 documentation (Ramboll CB6 mapping memo; Yarwood et al. 2010, Hildebrandt Ruiz & Yarwood 2013):

* `PAR` = 1 carbon (single-bonded); `OLE` = 2 C (terminal C=C); `IOLE` = 4 C (internal C=C; 2-butene = 1 IOLE); `KET` = 1 C (carbonyl); `ALDX` = 2 C (propanal = ALDX + PAR); `TOL` = 7 C; `XYLMN` = 8 C; `TERP` = 10 C.
* Explicit in CB6r5 / CMAQ: `ETHA, ETH, ETHY, ETOH, MEOH, FORM, ALD2, ACET, PRPA, ISOP, BENZENE, NAPH, APIN, TERP, AACD, FACD, GLY, GLYD, MGLY, CRES`.
* Biacetyl = 1.5 MGLY; methyl acetate = 2 PAR (+ unreactive carbon `UNR`, which is not a model species); ethylbenzene = TOL + PAR; furan ring = 2 OLE.
* `UNR` (unreactive carbon) is dropped, so carbon is **not** conserved where UNR applies (`carbon` columns below show it).

Basis codes: **direct** = explicit CB6r5 species; **doc** = stated in the CB6 documentation; **rule** = follows from the structural rules above; **assump** = my judgement (no explicit statement found) - these are the choices you may wish to change.

### 4.2 Compound -> CB6r5 map (all cases)

| Key | CB6r5 representation | C in compound | C in CB6r5 | Basis | Note |
|---|---|---|---|---|---|
| CO | CO |  | 1 | direct |  |
| CO2 | CO2 |  | 1 | direct | inert tracer (not in CB6r5) |
| NO | NO |  |  | direct |  |
| NO2 | NO2 |  |  | direct |  |
| O3 | O3 |  |  | direct |  |
| H2O2 | H2O2 |  |  | direct |  |
| HONO | HONO |  |  | direct |  |
| SO2 | SO2 |  |  | direct |  |
| HNO3 | HNO3 |  |  | direct |  |
| HCl | HCl |  |  | direct |  |
| NH3 | NH3 |  |  | direct |  |
| H2 | H2 |  |  | direct |  |
| CH4 | CH4 |  | 1 | direct |  |
| N2O5 | N2O5 |  |  | direct |  |
| PAN | PAN |  | 2 | direct |  |
| HCN | HCN |  | 1 | direct | inert tracer (not in CB6r5) |
| N2O | N2O |  |  | direct | inert tracer (not in CB6r5) |
| H2SO4 | H2SO4 |  |  | direct | CMAQ name SULF; ASP condensation/aqueous decks use H2SO4 |
| HCHO | FORM | 1 | 1 | direct |  |
| MEOH | MEOH | 1 | 1 | direct |  |
| FACD | FACD | 1 | 1 | direct |  |
| ETH | ETH | 2 | 2 | direct |  |
| ETHA | ETHA | 2 | 2 | direct |  |
| ETHY | ETHY | 2 | 2 | direct |  |
| ALD2 | ALD2 | 2 | 2 | direct |  |
| ETOH | ETOH | 2 | 2 | direct |  |
| AACD | AACD | 2 | 2 | direct |  |
| GLY | GLY | 2 | 2 | direct |  |
| GLYD | GLYD | 2 | 2 | direct | hydroxyacetaldehyde / glycolaldehyde (OHCH2CHO, HOCH2CHO) |
| MeFormate | PAR | 2 | 1 | rule | methyl formate: ester carbon is unreactive (UNR, not a CB6r5 model species) -> 1 PAR; cf. methyl acetate = 2 PAR + UNR |
| PRPA | PRPA | 3 | 3 | direct |  |
| Propene | OLE + PAR | 3 | 3 | doc | propene = OLE + PAR |
| Acrolein | 0.75 ISPD | 3 | 3 | assump | CMAQ BUTADIENE13/ACROLEIN are no-feedback toxics tracers, so acrolein is represented by ISPD (the CB6 unsaturated carbonyl lump, whose photolysis already uses acrolein cross sections), scaled to conserve carbon (3/4) |
| Propanal | ALDX + PAR | 3 | 3 | doc | propionaldehyde = ALDX + PAR |
| Acetone | ACET | 3 | 3 | direct |  |
| MGLY | MGLY | 3 | 3 | direct |  |
| Hydroxyacetone | KET + 2 PAR | 3 | 3 | assump | CH3C(O)CH2OH (ACETOL, C3H6O2): ketone group + 2 PAR; no explicit mapping found |
| nPropanol | 3 PAR | 3 | 3 | assump | C3+ alcohols other than methanol/ethanol: one PAR per carbon |
| C3H6O | OLE + PAR | 3 | 3 | assump | Timbavati Table S3 "C3H6O" is classified OLEL (lumped alkene) by the paper; mapped as propene |
| Butadiene | 2 OLE | 4 | 4 | assump | two terminal C=C bonds -> 2 OLE (CB treats each double bond separately; 4 C conserved) |
| 1-Butene | OLE + 2 PAR | 4 | 4 | rule |  |
| iButene | OLE + 2 PAR | 4 | 4 | rule | isobutene: terminal olefin + 2 PAR (branching ignored) |
| 2-Butene | IOLE | 4 | 4 | doc | cis/trans-2-butene = IOLE |
| C4H8 | OLE + 2 PAR | 4 | 4 | assump | Timbavati "C4H8" (unspeciated butenes): treated as 1-butene (as the baseline deck did) |
| Butane | 4 PAR | 4 | 4 | doc | butane = 4 PAR (n- and iso-) |
| Butanal | ALDX + 2 PAR | 4 | 4 | rule |  |
| MEK | KET + 3 PAR | 4 | 4 | doc | butanone = KET + 3 PAR |
| nButanol | 4 PAR | 4 | 4 | assump | C3+ alcohols: one PAR per carbon |
| Furan | 2 OLE | 4 | 4 | rule | furan ring = 2 C=C -> 2 OLE (cf. 2-butylfuran = 2 OLE + 4 PAR, pyrrole = 2 OLE in the Ramboll memo) |
| Biacetyl | 1.5 MGLY | 4 | 4.5 | doc | 2,3-butanedione = 1.5 MGLY (CB6 documentation) |
| ISOP | ISOP | 5 | 5 | direct |  |
| Pentane | 5 PAR | 5 | 5 | rule | n-, i-pentane, 2-methylbutane |
| 1-Pentene | OLE + 3 PAR | 5 | 5 | rule | 1-pentene, 3-methyl-1-butene, 2-methyl-1-butene |
| 2-Pentene | IOLE + PAR | 5 | 5 | doc | 2-pentene isomers = IOLE + PAR; also 2-methyl-2-butene, cyclopentene |
| Cyclopentanol | 5 PAR | 5 | 5 | assump | alcohol: one PAR per carbon |
| Pentanone | KET + 4 PAR | 5 | 5 | rule |  |
| Hexane | 6 PAR | 6 | 6 | rule | n-hexane, isohexanes |
| 1-Hexene | OLE + 4 PAR | 6 | 6 | rule | 1-hexene, 2-methyl-1-pentene, 4-methyl-1-pentene |
| Heptane | 7 PAR | 7 | 7 | rule |  |
| Heptanone | KET + 6 PAR | 7 | 7 | rule |  |
| Octene | OLE + 6 PAR | 8 | 8 | assump | octenes assumed terminal (1-octene) |
| Octanone | KET + 7 PAR | 8 | 8 | rule |  |
| MeAcetate | 2 PAR | 3 | 2 | doc | methyl acetate = 2 PAR + UNR |
| Benzene | BENZENE | 6 | 6 | direct |  |
| TOL | TOL | 7 | 7 | direct |  |
| Xylene | XYLMN | 8 | 8 | direct | xylene isomers |
| Ethylbenzene | TOL + PAR | 8 | 8 | doc | ethylbenzene = TOL + PAR |
| Styrene | TOL + OLE + -1 PAR | 8 | 8 | assump | Ramboll memo: styrenes = TOL, OLE and PAR (coefficients not given); TOL + OLE - PAR conserves 8 C |
| Phenol | CRES | 6 | 7 | assump | phenol -> CRES (cresol and higher phenols); per-mole, carbon not conserved |
| Benzaldehyde | TOL | 7 | 7 | assump | aromatic ring + one carbon: TOL (7 C); no explicit mapping found |
| PAH | NAPH |  | 10 | assump | PAHs represented as naphthalene (the only PAH in CB6r5) |
| APIN | APIN | 10 | 10 | direct | alpha-pinene (CMAQ explicit species) |
| TERP | TERP | 10 | 10 | direct | limonene and other monoterpenes; CB6 TERP = all monoterpenes |
| ALD_hi | ALDX + 3 PAR | 5 | 5 | assump | higher aldehydes (C3-C7): average C5 = ALDX + 3 PAR |
| RACM_ALD | ALDX + 3 PAR | 5 | 5 | assump | Williams "ALD" (acetaldehyde is listed separately, so C3+): same as C3-C7 aldehydes |
| RACM_HC5 | 5 PAR | 5 | 5 | assump | RACM2 HC5 (pentanes/hexanes) as 5 PAR |
| RACM_HC8 | 8 PAR | 8 | 8 | assump | RACM2 HC8 as 8 PAR |
| RACM_KET | KET + 3 PAR | 4 | 4 | assump | RACM2 KET (MEK and larger) as MEK |
| RACM_OLI | IOLE | 4 | 4 | assump | RACM2 OLI (internal olefins) as 2-butene |
| RACM_OLT | OLE + PAR | 3 | 3 | assump | RACM2 OLT (terminal olefins) as propene |
| RACM_ROH | 3 PAR | 3 | 3 | assump | RACM2 ROH (C3+ alcohols) as propanol |
| PPN | PANX | 3 | 3 | direct | peroxypropionyl nitrate -> PANX (higher PANs) |
| CH3NO3 | NTR1 | 1 | 1 | direct | methyl nitrate -> NTR1 (monofunctional organic nitrates) |

### 4.3 Renamed or placeholder species

* `SULF` -> `H2SO4`, `HCL` -> `HCl`, `ECH4` -> `CH4`: names matched to what ASP's dissolution/aqueous decks refer to.
* `CO2`, `HCN`, `N2O`, `Trace`: inert tracers not in CB6r5, carried over unchanged from the old deck.
* `NH3`, `HCl`, `H2SO4`, `HNO3`: interface species kept (zero if the paper gives no value).
* `CH3NO3`, `RO2_T`, `RCO3_T`: zero placeholders (section 2.5). `O2`/`N2` set to 2.095e8 / 7.808e8 ppb (CMAQ constants).
* **Peroxy flags**: `RO2` = 1; `C2O3`, `CXO3`, `OPO3` = 2; all other species 0. CB6r5 carries its peroxy radicals in the `RO2` pool and the `XO2*` operators, so only the pool is counted in `RO2_T` (avoids double counting; this controls the ASP autoxidation-type reactions of the SOA chains).
* Deposition velocities are 0 for every species (as in the old decks).
* POA1-8, LEVO, CBIO, CPD1-3 and IVOC1-9 rows are copied from the old deck unchanged because the organic dissolution decks refer to them.

### 4.4 Per-case values

Full per-compound tables are in each deck's `SPECIES_MAPPING.csv`. Notable per-case choices:

* **Williams** - the deck IVOC3-7 initial values (0.021, 0.255, 3.08, 37.8, 467 ppb) are carried unchanged; they differ from the paper supplement's Table S3 SVOC3-8 values (2.6e-3 ... 217 ppb). The old deck matches the supplement's own deck file, so this was left as is.
* **Timbavati** - `C4H8` (unspeciated butenes, 1.6 ppb) -> 1-butene; `C3H6O` (7.8 ppb, class OLEL) -> propene; O3 background 35.7 ppb (init 0) carried from the old deck (Table S3 gives no O3); CO2 380000 ppb carried; SO2 background 10 ppb as in S3 (the thesis says 1.45 ppb, an open item).
* **Otavi / Alaska** - values straight from Tables S1/S2.

## 5. Photolysis

`PhotoRates.in` keeps ASP's 51-column format. Column k of the new table is a linear combination of the case's existing TUV 5.4 columns (`PhotoRatesTUV/asp_photorates_map.py`); because ASP indexes photolysis by column number only, no source change is needed. The CB6r5 photolysis labels map as follows (`weights` apply to the old columns):

| New col | CMAQ label (x scale) | TUV species / channel | Old columns |
|---|---|---|---|
| 1 | NO2_IUPAC10 | NO2 -> NO + O(3P) | 1 x 6 |
| 2 | O3_O3P_IUPAC10 | O3 -> O2 + O(3P) | 1 x 3 |
| 3 | O3_O1D_IUPAC10 | O3 -> O2 + O(1D) | 1 x 2 |
| 4 | H2O2_IUPAC10 | H2O2 -> 2 OH | 1 x 5 |
| 5 | NO3NO2_06 | NO3 -> NO2 + O(3P) | 1 x 8 |
| 6 | NO3NO_06 | NO3 -> NO + O2 | 1 x 7 |
| 7 | N2O5_IUPAC10 | N2O5 -> NO3 + NO2 (+ NO3 + NO + O): sum of both TUV channels | 1 x 9 + 1 x 10 |
| 8 | HONO_IUPAC10 | HNO2 -> OH + NO | 1 x 12 |
| 9 | HNO3_IUPAC10 | HNO3 -> OH + NO2 | 1 x 13 |
| 10 | PNA_IUPAC10 | HNO4 -> HO2 + NO2 | 1 x 14 |
| 11 | PAN_IUPAC10 | PAN: sum of the NO2 and NO3 channels (CB6r5 splits 0.6/0.4 in the products) | 1 x 40 + 1 x 41 |
| 12 | MEPX_IUPAC10 | CH3OOH -> CH3O + OH | 1 x 31 |
| 13 | NTR_IUPAC10 | isopropyl nitrate -> products + NO2 | 1 x 36 |
| 14 | FORM_R_IUPAC13 | CH2O -> H + HCO | 1 x 17 |
| 15 | FORM_M_IUPAC13 | CH2O -> H2 + CO | 1 x 18 |
| 16 | ALD2_R_IUPAC13 | CH3CHO -> CH3 + HCO | 1 x 19 |
| 17 | ALDX_R_IUPAC13 | C2H5CHO -> C2H5 + HCO | 1 x 22 |
| 18 | GLYD_IUPAC13 | HOCH2CHO: sum of the three TUV channels | 1 x 23 |
| 19 | GLY_R_IUPAC13 | CHOCHO -> HCO + HCO | 1 x 44 |
| 20 | MGLY_IUPAC10 | CH3COCHO -> CH3CO + HCO | 1 x 46 |
| 21 | KET_IUPAC10 | MEK -> CH3CO + C2H5 | 1 x 28 |
| 22 | ACET_IUPAC10 | acetone -> CH3CO + CH3 | 1 x 26 |
| 23 | ISPD | acrolein -> products (CB6 ISPD uses acrolein cross sections) | 1 x 24 |
| 24 | HPALD | ASP convention: 100 x J(MACR), table column 51 | 1 x 51 |
| 25 | NO2_IUPAC10 x 0.05 | NO2 -> NO + O(3P) | 0.05 x 6 |
| 26 | NO2_IUPAC10 x 0.028 | NO2 -> NO + O(3P) | 0.028 x 6 |

Columns 27-51 are zero. Approximations: the CB6r5 `ISPD` label uses acrolein cross sections (old col 24); `HPALD` uses old col 51 (ASP convention, 100 x J(MACR)); NTR/NO3 photolysis use TUV's isopropyl nitrate and the two NO3 channels.

## 6. SOA coupling

The baseline mechanism counts every VOC oxidation event with a dummy species (TOLSOA, XYLSOA, ISOSOA, TERPSOA, ...), which then decays into the IVOC1-9 volatility bins (the SOA yields). The CB6r5 translation keeps this coupling:

* The CMAQ counters `TOLRO2, XYLRO2, ISOPRXN, TRPRXN, TERPNRO2` are re-keyed to `TOLSOA, XYLSOA, ISOSOA, TERPSOA, TERPSOA`; `BENZRO2, PAHRO2, SULRXN` are dropped (no ASP analogue).
* The yield chains for these four counters and the IVOC aging block (Grieshop for Otavi/Alaska/Timbavati; the ACP 2015 half-Ahmadov scheme for Williams, with `Acetic_Acid` renamed `AACD`) are copied verbatim from each case's old deck.

**Open decisions**: the HC5SOA/HC8SOA (alkane) and OLTSOA/OLISOA (olefin) chains are not re-attached (their old inventories were negligible or uncounted); ISOP + O3 is not counted as an SOA event in CB6r5 (the baseline counted it); explicit propene and the butenes carry no counter (as before).

## 7. Test runs (same ASP binary, same input decks except the three replaced files)

| Case | O3 end, baseline (ppb) | O3 end, CB6r5 (ppb) | NO2 end base / CB6r5 | SOA base / CB6r5 | run time (s) |
|---|---|---|---|---|---|
| Williams | 114.2 | 125.4 | 7.3 / 2.4 | 35.0 / 39.8 | 12.691 |
| Otavi | 96.5 | 98.2 | 1.9 / 1.0 | 11.0 / 11.0 | 6.185 |
| Alaska | 77.2 | 86.9 | 25.0 / 23.8 | 26.3 / 26.4 | 26.631 |
| Timbavati | 60.3 | 66.6 | 12.6 / 11.9 | 17.4 / 17.3 | 10.810 |

CB6r5 end-of-run O3 is 2-13% above the baseline mechanism and SOA within 1% (Williams +14%); all four decks load, integrate stably and finish.

## 8. Open items

1. Decide whether to fix ASP rate flag 5 (and flag 2 3-body) in `common_asp`; this changes all baseline results (section 2.4).
2. Review the **assump** rows of section 4.2 (alcohols, hydroxyacetone, benzaldehyde, phenol, styrene, octenes, Williams ALD/KET/OLI/OLT lumps).
3. Alkane/olefin SOA chains and ISOP+O3 counter (section 6).
4. Heterogeneous chemistry (HET_*) and the NTR2 hydrolysis sink are not translated.
5. Inherited open items from the baseline decks: Timbavati OA ~half the paper value; SO2 background 10 vs 1.45 ppb; Williams IVOC3-7 vs Table S3.

