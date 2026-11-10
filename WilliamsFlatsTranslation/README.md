# Williams Flats Fire (2019) input decks

`build_williamsflats_decks.py` writes nine ASP input-deck folders for the
**Williams Flats Fire** (Washington, 3 August 2019, sampled in FIREX-AQ):

| Folder | Mechanism |
|---|---|
| `InorganicInputDecks/WilliamsFlats/` | baseline (RACM2-derived, as `Williams/`) |
| `InorganicInputDecks/WilliamsFlats_<MECH>/` | CB6r5, CRACMM2, CRACMM3, GCHEM, MCM331, MOZT1, SAPRC07TC, SAPRC07TIC |

Rerun the script after changing anything below. `--check` only runs the
validation (rebuilding the `Williams_<MECH>` values from `Williams/`).

## Source

The SAM-ASP-TUV v1.0.0 input decks (Zenodo,
[10.5281/zenodo.13385240](https://doi.org/10.5281/zenodo.13385240);
`njune74/SAM-ASP-TUV`, `SRC_TUV/ASP/InorganicInputDecks`). Those decks were
used for the SAM-ASP-TUV Williams Flats runs (SAM case `SND1`, run
`SND1E1S1`, 47.94 N, 118.62 W). The case-specific files are kept unchanged in
`zenodo_samasp_decks/`.

The SAM-ASP-TUV repository keeps that deck set, with its original organic
aerosol scheme, as its own Williams Flats case.

## Baseline deck (`WilliamsFlats/`)

It is the `Williams/` baseline deck with these changes:

| File | Change |
|---|---|
| `GasPhaseChems.in` | The SAM deck gives emission factors (g/kg, flag 0). They are converted to initial mixing ratios with ASP's own flag-0 formula at an initial CO of 10,000 ppb, and written as totals (flag 2), like every ASP case deck. Backgrounds are the SAM deck's. |
| | Organic volatility bins: the SAM deck has two 4-bin sets (IVOC1-4 emitted, IVOC5-8 for SOA) whose C* and molecular weights equal the baseline bins IVOC1, 3, 5 and 7, so IVOC1/5 -> IVOC1, IVOC2/6 -> IVOC3, IVOC3/7 -> IVOC5, IVOC4/8 -> IVOC7. The baseline's OH aging of the 9 bins applies. |
| | `OXYARO` and `OTHERO`, two lumped SOA precursors of the SAM deck, are added. |
| `GasChemicalMechanism.in` | `OXYARO + OH` and `OTHERO + OH` from the SAM deck (OH is regenerated; products remapped to the 9 bins). |
| `AerosolModes.in` | The SAM deck's smoke mode ("Williams Flats Mod": 634,000 cm-3, 0.09 um, sigma 1.5; mass fractions NH4Cl 1, NH4NO3 2, KCl 0.5, K2SO4 1, BC 2.5, organics 93 in IVOC1/3/5/7). The baseline's 40-bin sectional set-up is kept (the SAM deck used 10 bins). |
| `ASPInputDeck.in` | T 297.0 K, P 937.3 hPa, RH 0.249: the SAM `SND1` sounding interpolated to its surface pressure. Photolysis start 19:36 UTC (SND1 `day0` 215.817). Mixing parameter 15 s (SAM deck). Everything else as `Williams/`, including the run length (300 min), the "hetero" particles and the fire CO emission line (40 kg/s, a Williams value). |
| `PhotoRates.in` | TUV table for the site and date. **Provisional**, see below. |

`EnvAerosolModes.in` is unchanged: the SAM deck's background aerosol has the
same composition as `Williams/`. The organic-aerosol property files
(`HydrophilicOrgChems.in` etc.) are the baseline's: the SAM deck gave the
organics a hygroscopicity kappa of 0.04; the baseline uses 0.

## Mechanism decks (`WilliamsFlats_<MECH>/`)

Each is its `Williams_<MECH>/` deck with:

* **Gas species**: the deck's own species-mapping CSV
  (`SPECIES_MAPPING.csv`, `GCHEM_`, `MCM_`, `MOZT1_`) applied to the
  `WilliamsFlats/` values. The CSV in each folder now lists the Williams
  Flats values.
* **Organic aerosol** (gas and particle), with the rules each mechanism's
  builder used for Williams:
  * CRACMM2/3: bin k -> the CMAQ species of the same C* (`cmaq_oa.BIN_MAP`), mass conserved.
    IVOC7 maps to VROCP4OXY2 (C* = 1e4), which is gas-only in these decks: its gas
    is kept, but its particle mass (40% of the smoke organics) is left out of the
    smoke mode (at plume loadings it would be almost entirely in the gas phase).
  * CB6r5, SAPRC07TC/TIC (AERO7): bins 1-2 -> VLVPO1, 3 -> VSVPO1, 5 -> VSVPO3
    (`cmaq_ae.BIN_MAP`); bins with C* >= 1e4 have no AERO7 species and are
    dropped, as for Williams. **This drops IVOC7: 40% of the smoke organic
    particle mass and 40% of the emitted organic gas mass in the bins.** The SAPRC07 decks gain `VLVPO1`
    (rows copied from `Williams_CB6r5`; partitioning only, like their other
    POA species).
  * GEOS-Chem: gas bins <= 6 and all particle organics -> SVPOA1/SVPOA2
    (0.49/0.51); gas bins >= 7 -> NAP (added to emitted naphthalene).
  * MOZART-T1: gas bins <= 6 -> SVOC, >= 7 -> IVOC; particle organics -> POM.
  * MCM331: the baseline bins, unchanged.

  `OA_MAPPING_WilliamsFlats.txt` in each folder lists the result.
* **OXYARO / OTHERO**: added, with their IVOC products mapped by the same rules.
* **`ASPInputDeck.in`**: the `WilliamsFlats/` one.
* **`PhotoRates.in`**: GEOS-Chem, MOZART-T1 and MCM331 use the baseline table.
  For CB6r5, CRACMM and SAPRC07, every column is rebuilt from the
  `WilliamsFlats/` TUV table as the `Williams_<MECH>` header describes it: a
  copy, a (weighted) sum of TUV columns, or a TUV column times a CMAQ
  cross-section x quantum-yield ratio. The ratio is taken from the Williams
  pair as a function of solar zenith angle; below the smallest Williams angle
  (53.9 deg; Williams Flats reaches 30.5 deg) it is held at its last value.
  This is an approximation for those few columns.

### Validation

Applying the same rules to `Williams/` rebuilds the gas values and first
aerosol mode of all eight `Williams_<MECH>` decks (within 0.02%), and the
copied and summed photolysis columns of the CMAQ-family tables (within 0.2%).
`--check` repeats this.

## Photolysis: provisional table

TUV 5.4 (used for every other ASP case table) was not available.
`PhotoRates_TUV531.in` was made with `PhotoRatesTUV/make_photorates.py` and the
TUV 5.3.1 copy in SAM-ASP-TUV (`SRC_TUV/TUV`; in a scratch copy, `params` kw
raised from 8 to 1000 and a debug-print block in `rtrans.f` commented out),
with the options of the Williams table:

```
make_photorates.py --tuv <TUV> --lat 47.9392 --lon -118.6183 --date 2019-08-03 --alt 0.66 \
  --pbl-height 1.0 --zout 0.5 --tau-aer 0.5 --wl-aer 330 --ssa-uv 0.9 --ssa-vis 0.9 \
  --o3-ppb 50 --no2-ppb 0.04 --so2-ppb 0.12
```

For the Williams inputs, this TUV 5.3.1 gives J-values about 30-60% higher than
the TUV 5.4 table in `Williams/PhotoRates.in` (J(NO2) +36%, J(O1D) +30% at
54 deg solar zenith angle); the cause was not found. **Regenerate the table
with TUV 5.4** (same command, `--tuv /path/to/V5.4 -o PhotoRates_TUV531.in`)
and rerun the script before using these decks for science.

## Test runs

ASP parcel model (common_ASP 111426c), all nine decks: see
`test_runs.txt`.
