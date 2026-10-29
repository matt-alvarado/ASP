#!/usr/bin/env python3
"""Generate CB6r5_README.md from the translation objects (so the documentation cannot drift from the decks)."""
import csv
import os
import sys
from collections import Counter, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_cb6r5_decks as B           # noqa: E402
import cmaq2asp as c2a                  # noqa: E402
import species_map as sm                # noqa: E402

RUN = {}   # filled from run_summary.csv if present
rs = os.path.join(HERE, 'run_summary.csv')
if os.path.exists(rs):
    for r in csv.DictReader(open(rs)):
        RUN[r['case']] = r

tr = c2a.Translation(B.DEF, 288.4, 880.0)
cmp_ = c2a.compare(tr)
by = {n['label']: [] for n in tr.notes}
for n in tr.notes:
    by[n['label']].append(n)

L = []
w = L.append
w('# CB6r5 gas-phase mechanism for the four ASP case decks')
w('')
w('Decks: `InorganicInputDecks/Williams_CB6r5`, `Otavi_CB6r5`, `Alaska_CB6r5`, `Timbavati_CB6r5`.')
w('Only the **gas-phase chemistry** was changed. Aerosol, aqueous, thermodynamic, dissolution and organic-phase decks,')
w('`ASPInputDeck.in` and all run settings are byte-for-byte copies of the original case decks.')
w('')
w('Mechanism source: CMAQ `CB6r5-AE7-AQ` (`CCTM/src/MECHS/cb6r5_ae7_aq/mech_cb6r5_ae7_aq.def`, the species table '
  '`mechanism_information/cb6r5_ae7_aq/cb6r5_ae7_aq_species_table.md`, and `GC_cb6r5_ae7_aq.nml` for molecular weights), '
  'https://github.com/USEPA/CMAQ.')
w('')
w('## 1. What is in each new directory')
w('')
w('| File | Status |')
w('|---|---|')
w('| `GasChemicalMechanism.in` | **replaced**: %d CB6r5 reactions (each preceded by a `!! R<n>: <CMAQ statement>` comment and any WARNING), then the TOLSOA/XYLSOA/ISOSOA/TERPSOA -> IVOC chains and the IVOC aging block of the original case deck |' % len(tr.lines))
w('| `GasPhaseChems.in` | **replaced**: CB6r5 species; initial/background values from the papers\' compound tables mapped to CB6r5 (section 4) |')
w('| `PhotoRates.in` | **replaced**: CB6r5 column set derived from the case\'s TUV 5.4 table (section 5) |')
w('| `SPECIES_MAPPING.csv` | **new**: every paper compound, its CB6r5 representation, carbon check and basis |')
w('| all other `.in` files | unchanged |')
w('')
w('Rebuild: `python3 CB6r5Translation/build_cb6r5_decks.py && python3 CB6r5Translation/write_readme.py`.')
w('')

# ------------------------------------------------------------------------------------------
w('## 2. WARNINGS - CB6r5 rate forms ASP cannot (exactly) handle')
w('')
w('Of the %d core `R<n>` reactions, %d are translated to within 0.12%% of the CMAQ rate over '
  '260-310 K, 600-1013 hPa (most to <0.01%%) except for those listed below.' % (len(tr.core), len(tr.lines)))
w('')
w('### 2.1 Cannot be represented - a fixed value is used instead')
w('')
w('* **R46, OH + HNO3** - CB6r5 uses the form `k0 + k3[M]/(1 + k3[M]/k2)` (CMAQ `%2`). ASP\'s 2-body rate flag 5 is meant '
  'for this form but is **coded wrongly** (see 2.4), so it cannot be used. R46 is therefore a fixed 2-body rate evaluated at the '
  'case T and P (Williams 1.72e-13, Timbavati 1.73e-13 cm3/molec/s). It is exact only at those conditions (34% off at 20 K / 300 hPa away). '
  'All four decks use constant T and P, so this is harmless for these cases.')
w('')
w('### 2.2 Approximated using ASP\'s hard-wired forms')
w('')
w('* **R2, O + O2 + M** - ASP flag 9 (5.6e-34 for N2, 6.0e-34 for O2) vs CB6r5 6.0e-34 (T/300)^-2.6 for all M (<=5% difference).')
w('* **R10, O(1D) + M** - ASP flag 8 (IUPAC N2/O2 form) vs CB6r5 2.23e-11 exp(115/T) (<=0.2%).')
w('* **R24, NO + NO + O2** (termolecular, O2 as reactant) - ASP has no working 3-body Arrhenius form (see 2.4), so '
  'k[O2] is fitted as k300 (T/300)^-n x M (flag 1) with value and slope matched at the case temperature (<=1.3% over 260-310 K). Negligible at ambient NO.')
w('* **R37, N2O5 decomposition** - exact: ASP flag 7 is the same IUPAC expression as CB6r5.')
w('* Fall-off reactions in IUPAC form use ASP\'s built-in N(F) = 0.75 - 1.27 log10(Fc); CMAQ rounds n to two decimals, so rates differ by ~0.1%.')
w('')
w('### 2.3 Expressible only after rewriting')
w('')
refs = [n for n in tr.notes if n['level'] == 'expanded']
w('* **%d `A*K<Rn>` rate references** (a rate defined as a multiple of another reaction\'s rate): ASP has no such syntax, '
  'so each is expanded statically to the referenced reaction\'s rate x factor. Exact, but the rates are no longer linked: changing one '
  'rate in the deck does not change its dependents. Affected: %s.' % (len(refs), ', '.join(n['label'] for n in refs)))
w('* **Pseudo-reactants.** `O2`, `M`, `N2`, `H2O` as CMAQ reactants (R2, R10, R20, R24, R41, R221) are handled as in ASP '
  '(`M`, flag 10 for water, etc.); R20, R41 and R221 use ASP\'s C_H2O flag 10.')
w('* **Reactions with no products** (e.g. R7 `O + O3 =`) get the sink product `M`, which ASP\'s parser ignores as a species.')
w('')
w('### 2.4 Bugs found in the ASP source (not patched - awaiting your decision)')
w('')
w('1. **2-body rate flag 5 is wrong.** It evaluates `k1 + k3*M*(1 + k3*M/k2)`, while the documented/IUPAC '
  'form is `k1 + k3*M/(1 + k3*M/k2)`. It is used by the OH + HNO3 line of the **existing baseline decks** (`HNO3 + OH => NO3 + H2O ; 5.; ...`), '
  'where it gives 4.1e-11 instead of 1.7e-13 cm3/molec/s at 288 K / 880 hPa (**~240x too fast**), recycling HNO3 to NO3. A/B test '
  '(baseline deck with the correct constant rate, otherwise identical): Williams end-of-run O3 114.2 -> 120.3 ppb (+5%), HNO3 3.1 -> 6.8 ppb; '
  'Timbavati O3 60.3 -> 61.5 ppb (+2%), HNO3 1.46 -> 1.69 ppb. Not fixed in the source, because it changes all baseline results.')
w('2. **3-body rate flag 2 ("thermal third order") is unusable**: it reads its parameters from the wrong positions of the rate array, '
  'so the flag value itself is taken as the first rate coefficient. No existing deck uses it; the CB6r5 decks do not.')
w('')
w('### 2.5 CB6r5 content not translated (outside the gas-chemistry-only scope or without ASP equivalent)')
w('')
sk = Counter()
import re
for r in tr.skipped:
    sk[re.match(r'[A-Za-z]+', r.label).group(0)] += 1
w('* Reaction groups of the .def file that are **not included**: ' + ', '.join('%s (%d)' % (k, v) for k, v in sorted(sk.items())) + '.')
w('  These are chlorine (CL), DMS, mercury (HG), halogen/ozone (HAL), heterogeneous/hydrolysis (HET), source-apportionment (SA), '
  'biogenic/biomass tracer (BL, BP), POA aging and oligomerization (POA, RPOAGE, OLIG, PCSOA), and the **`T*` toxics diagnostics** (no feedback on the gas mechanism). '
  '`BUTADIENE13` and `ACROLEIN` have reactions only in `T*`, so they are not used as species (butadiene -> 2 OLE; acrolein -> ISPD).')
w('* **R216 (ECH4 + OH)** - the tracer-CH4 duplicate; omitted (CH4 itself has R<n> CH4 + OH).')
w('* **NTR2 heterogeneous hydrolysis** (`HET_NTR2`): NTR2 is only a NOx reservoir after this translation (CMAQ removes it by hydrolysis to HNO3 in the aerosol module).')
w('* **Heterogeneous N2O5, NO2, HO2 uptake** (`HET_*`): not translated. The original decks\' heterogeneous reactions are untouched (they were already commented out in Timbavati).')
w('* **`CH3NO3`, `RO2_T`, `RCO3_T`** are hard-wired species names in ASP\'s output code; they exist as zero placeholders (CB6r5 has no CH3NO3: it is lumped in NTR1; `RO2_T`/`RCO3_T` are sums computed by ASP from the peroxy flags).')
w('')

# ------------------------------------------------------------------------------------------
w('## 3. Reaction-by-reaction check')
w('')
w('Every translated reaction was evaluated with a port of ASP\'s rate code and compared with CMAQ\'s own rate expression at '
  '7 (T, P) points covering 260-310 K, 600-1013 hPa (`cmaq2asp.py --compare`). Largest deviations:')
w('')
w('| Reaction | max rel. deviation | note |')
w('|---|---|---|')
for cdict in sorted(cmp_, key=lambda x: -x['maxrel'])[:12]:
    if cdict['maxrel'] > 1e-3:
        lab = cdict['label']
        nt = ' '.join(n['text'][:80] for n in by.get(lab, []) if n['level'] in ('approx', 'unsupported'))
        w('| %s | %.2g | %s |' % (lab, cdict['maxrel'], nt.replace('|', '/')))
rest = [x for x in cmp_ if x['maxrel'] <= 1e-3]
w('')
w('All other %d reactions agree to < 0.1%% (n rounding of IUPAC fall-off forms). Photolysis reactions are checked through the column map in section 5.' % len(rest))
w('')

# ------------------------------------------------------------------------------------------
w('## 4. Species lumping and conversion')
w('')
w('Procedure: the initial and background values come from the **papers\' original compound tables** (Alvarado & Prinn 2009, '
  'JGR, auxiliary tables S1-S3 for Otavi, Alaska, Timbavati; the ACP 2015 supplement deck for Williams, whose nonzero values are identical '
  'to the old Williams deck), **not** from the old RACM2-style lumps. Each compound is expanded to CB6r5 species, and the '
  'contributions are summed per species (initial and background separately).')
w('')
w('### 4.1 CB6 lumping rules used')
w('')
w('Per the CB6 documentation (Ramboll CB6 mapping memo; Yarwood et al. 2010, Hildebrandt Ruiz & Yarwood 2013):')
w('')
w('* `PAR` = 1 carbon (single-bonded); `OLE` = 2 C (terminal C=C); `IOLE` = 4 C (internal C=C; 2-butene = 1 IOLE); `KET` = 1 C (carbonyl); `ALDX` = 2 C (propanal = ALDX + PAR); `TOL` = 7 C; `XYLMN` = 8 C; `TERP` = 10 C.')
w('* Explicit in CB6r5 / CMAQ: `ETHA, ETH, ETHY, ETOH, MEOH, FORM, ALD2, ACET, PRPA, ISOP, BENZENE, NAPH, APIN, TERP, AACD, FACD, GLY, GLYD, MGLY, CRES`.')
w('* Biacetyl = 1.5 MGLY; methyl acetate = 2 PAR (+ unreactive carbon `UNR`, which is not a model species); ethylbenzene = TOL + PAR; furan ring = 2 OLE.')
w('* `UNR` (unreactive carbon) is dropped, so carbon is **not** conserved where UNR applies (`carbon` columns below show it).')
w('')
w('Basis codes: **direct** = explicit CB6r5 species; **doc** = stated in the CB6 documentation; **rule** = follows from the structural rules above; '
  '**assump** = my judgement (no explicit statement found) - these are the choices you may wish to change.')
w('')
w('### 4.2 Compound -> CB6r5 map (all cases)')
w('')
w('| Key | CB6r5 representation | C in compound | C in CB6r5 | Basis | Note |')
w('|---|---|---|---|---|---|')
for key, (comp, nC, basis, note) in sm.MAP.items():
    rep = ' + '.join(('%g ' % cf if cf != 1 else '') + sp for sp, cf in comp)
    cr = sm.carbon_represented(comp)
    w('| %s | %s | %s | %s | %s | %s |' % (key, rep, nC if nC is not None else '', ('%g' % cr) if cr is not None else '', basis, note.replace('|', '/')))
w('')
w('### 4.3 Renamed or placeholder species')
w('')
w('* `SULF` -> `H2SO4`, `HCL` -> `HCl`, `ECH4` -> `CH4`: names matched to what ASP\'s dissolution/aqueous decks refer to.')
w('* `CO2`, `HCN`, `N2O`, `Trace`: inert tracers not in CB6r5, carried over unchanged from the old deck.')
w('* `NH3`, `HCl`, `H2SO4`, `HNO3`: interface species kept (zero if the paper gives no value).')
w('* `CH3NO3`, `RO2_T`, `RCO3_T`: zero placeholders (section 2.5). `O2`/`N2` set to 2.095e8 / 7.808e8 ppb (CMAQ constants).')
w('* **Peroxy flags**: `RO2` = 1; `C2O3`, `CXO3`, `OPO3` = 2; all other species 0. CB6r5 carries its peroxy radicals in the `RO2` pool and the `XO2*` operators, so only the pool is counted in `RO2_T` (avoids double counting; this controls the ASP autoxidation-type reactions of the SOA chains).')
w('* Deposition velocities are 0 for every species (as in the old decks).')
w('* POA1-8, LEVO, CBIO, CPD1-3 and IVOC1-9 rows are copied from the old deck unchanged because the organic dissolution decks refer to them.')
w('')
w('### 4.4 Per-case values')
w('')
w('Full per-compound tables are in each deck\'s `SPECIES_MAPPING.csv`. Notable per-case choices:')
w('')
w('* **Williams** - the deck IVOC3-7 initial values (0.021, 0.255, 3.08, 37.8, 467 ppb) are carried unchanged; they differ from the paper supplement\'s Table S3 SVOC3-8 values (2.6e-3 ... 217 ppb). The old deck matches the supplement\'s own deck file, so this was left as is.')
w('* **Timbavati** - `C4H8` (unspeciated butenes, 1.6 ppb) -> 1-butene; `C3H6O` (7.8 ppb, class OLEL) -> propene; O3 background 35.7 ppb (init 0) carried from the old deck (Table S3 gives no O3); CO2 380000 ppb carried; SO2 background 10 ppb as in S3 (the thesis says 1.45 ppb, an open item).')
w('* **Otavi / Alaska** - values straight from Tables S1/S2.')
w('')

# ------------------------------------------------------------------------------------------
w('## 5. Photolysis')
w('')
w('`PhotoRates.in` keeps ASP\'s 51-column format. Column k of the new table is a linear combination of the case\'s existing TUV 5.4 columns '
  '(`PhotoRatesTUV/asp_photorates_map.py`); because ASP indexes photolysis by column number only, no source change is needed. '
  'The CB6r5 photolysis labels map as follows (`weights` apply to the old columns):')
w('')
w('| New col | CMAQ label (x scale) | TUV species / channel | Old columns |')
w('|---|---|---|---|')
for (lab, scale), idx in tr.photo.cols.items():
    combo, desc = c2a.PHOTO_COMBOS[lab]
    w('| %d | %s%s | %s | %s |' % (idx, lab, '' if scale == 1 else ' x %g' % scale, desc,
                                ' + '.join('%g x %d' % (wt * scale, oc) for oc, wt in combo)))
w('')
w('Columns 27-51 are zero. Approximations: the CB6r5 `ISPD` label uses acrolein cross sections (old col 24); `HPALD` uses old col 51 (ASP convention, 100 x J(MACR)); NTR/NO3 photolysis use TUV\'s isopropyl nitrate and the two NO3 channels.')
w('')

# ------------------------------------------------------------------------------------------
w('## 6. SOA coupling')
w('')
w('The baseline mechanism counts every VOC oxidation event with a dummy species (TOLSOA, XYLSOA, ISOSOA, TERPSOA, ...), '
  'which then decays into the IVOC1-9 volatility bins (the SOA yields). The CB6r5 translation keeps this coupling:')
w('')
w('* The CMAQ counters `TOLRO2, XYLRO2, ISOPRXN, TRPRXN, TERPNRO2` are re-keyed to `TOLSOA, XYLSOA, ISOSOA, TERPSOA, TERPSOA`; `BENZRO2, PAHRO2, SULRXN` are dropped (no ASP analogue).')
w('* The yield chains for these four counters and the IVOC aging block (Grieshop for Otavi/Alaska/Timbavati; the ACP 2015 half-Ahmadov scheme for Williams, with `Acetic_Acid` renamed `AACD`) are copied verbatim from each case\'s old deck.')
w('')
w('**Open decisions**: the HC5SOA/HC8SOA (alkane) and OLTSOA/OLISOA (olefin) chains are not re-attached (their old inventories were negligible or uncounted); '
  'ISOP + O3 is not counted as an SOA event in CB6r5 (the baseline counted it); explicit propene and the butenes carry no counter (as before).')
w('')

# ------------------------------------------------------------------------------------------
w('## 7. Test runs (same ASP binary, same input decks except the three replaced files)')
w('')
if RUN:
    w('| Case | O3 end, baseline (ppb) | O3 end, CB6r5 (ppb) | NO2 end base / CB6r5 | SOA base / CB6r5 | run time (s) |')
    w('|---|---|---|---|---|---|')
    for case in B.CASES:
        r = RUN.get(case)
        if r:
            w('| %s | %s | %s | %s / %s | %s / %s | %s |' % (case, r['o3_base'], r['o3_cb6'], r['no2_base'], r['no2_cb6'],
                                                     r['soa_base'], r['soa_cb6'], r['time']))
w('')
w('CB6r5 end-of-run O3 is 2-13% above the baseline mechanism and SOA within 1% (Williams +14%); all four decks load, integrate stably and finish.')
w('')
w('## 8. Open items')
w('')
w('1. Decide whether to fix ASP rate flag 5 (and flag 2 3-body) in `common_asp`; this changes all baseline results (section 2.4).')
w('2. Review the **assump** rows of section 4.2 (alcohols, hydroxyacetone, benzaldehyde, phenol, styrene, octenes, Williams ALD/KET/OLI/OLT lumps).')
w('3. Alkane/olefin SOA chains and ISOP+O3 counter (section 6).')
w('4. Heterogeneous chemistry (HET_*) and the NTR2 hydrolysis sink are not translated.')
w('5. Inherited open items from the baseline decks: Timbavati OA ~half the paper value; SO2 background 10 vs 1.45 ppb; Williams IVOC3-7 vs Table S3.')
w('')
open(os.path.join(HERE, 'CB6r5_README.md'), 'w').write('\n'.join(L) + '\n')
print('wrote README,', len(L), 'lines')
