# PhotoRatesTUV: make `PhotoRates.in` with NCAR TUV 5.4

`make_photorates.py` runs TUV 5.4 for one UT day at a site you specify and writes
an ASP `PhotoRates.in` (96 quarter-hour rows from UT 0.00, solar zenith angle,
51 J-values, end-of-file sentinel), with optional boundary-layer (PBL)
O3 / NO2 / SO2 / aerosol optical depths.

Needs: Python 3.8+, `gfortran`, `make`, `patch`, and the TUV 5.4 distribution
(the folder containing `TUV.f`, or the zip).  TUV is not redistributed here.

```
python3 make_photorates.py --tuv /path/to/V5.4 \
    --lat 35.2 --lon -112.2 --date 2010-06-21 --alt 2.0 \
    --pbl-height 1.0 --tau-aer 0.3 --tau-no2 0.02 --tau-so2 0.002 --tau-o3 0.05 \
    -o ../InorganicInputDecks/Williams/PhotoRates.in
```

First run compiles a patched copy of TUV (about 10 s) into `--workdir`
(default `./tuv_work`); later runs reuse it.  The TUV directory is never
modified.  A run takes about 10 s.

## Inputs

| Option | Meaning |
|---|---|
| `--lat --lon --date` | degrees N / E (west negative); **UT** date |
| `--alt` | ground altitude, km asl (TUV `zstart`) |
| `--zout` | height of the J-values above ground, km (default 0 = surface) |
| `--albedo`, `--o3-column`, `--psurf` | surface albedo (0.1), total O3 column in DU (300), surface pressure in mb (default: standard atmosphere at `--alt`) |
| `--pbl-height` | PBL depth above ground, km.  Rounded **down** to the 0.5 km grid.  0 = no PBL perturbation |
| `--tau-aer [--wl-aer]` | PBL aerosol optical depth at `--wl-aer` (default 330 nm), scaled to 330 nm as 1/wavelength (TUV's own PBL assumption) |
| `--ssa-uv --ssa-vis` | PBL aerosol single-scattering albedo below / above 400 nm (TUV's PBL defaults 0.6 / 0.9) |
| `--tau-o3 --tau-no2 --tau-so2` | PBL **gas absorption** optical depth at `--wl-o3` (310 nm), `--wl-no2` (400 nm), `--wl-so2` (300 nm) |
| `--o3-ppb --no2-ppb --so2-ppb` | alternative to the optical depths: PBL mixing ratio, ppb |

How the PBL inputs are used (TUV's own PBL mode, enabled by the patch):

* Inside the PBL the gas is uniformly mixed.  A gas optical depth is converted
  to a mixing ratio with TUV's surface-layer cross section in the wavelength
  bin containing the reference wavelength and the air column of the PBL; the
  derived ppb and bin are printed and written into the file header.
* Above the PBL: O3 follows the standard-atmosphere profile scaled to
  `--o3-column`; NO2 and SO2 are zero; aerosol follows Elterman (1968).
* A PBL gas that is not given is not inserted (O3 keeps its background value).
  Aerosol not given keeps the Elterman profile.
* Background aerosol is the Elterman profile scaled to a column optical depth
  above the **ground** of 0.235 at 550 nm (TUV default; `--bg-aer550`).  At a
  high-altitude site that is a heavy background, so lower it if appropriate.
* TUV uses the US Standard Atmosphere temperature and pressure profiles, not
  your parcel T and P.

## Output / mapping

Columns follow the header of the existing ASP table.  TUV 5.4 numbers its
reactions differently from TUV 5.0, so columns are mapped **by reaction name**
(`asp_photorates_map.py`).  Notes:

* Column 23 (glycolaldehyde) is the sum of TUV 5.4's three HOCH2CHO channels.
* Columns 49, 50, 51 are 0.33 x, 0.16 x column 22 and 100 x column 25, as in
  the existing table.
* Columns 15, 16 (nitrate in water) are TUV's snow/aqueous cross sections.
* Values below 1e-99 are written as 0.

## Patch (`tuv_asp_pbl.patch`)

TUV 5.4 already contains a PBL option, but it is hard-wired off and interactive.
The patch (3 files) makes TUV run non-interactively, reads PBL inputs from
`pbl.inp`, converts gas optical depths to mixing ratios, takes the aerosol SSA
as inputs, keeps the Elterman profile when no PBL aerosol is given, and keeps
the background O3 when no PBL O3 is given.

## Validation against the template table

The `PhotoRates.in` in `defaults/` is consistent with about 35N, 120E, 17 November
(fitted from its solar zenith angles, rms 0.3 deg), generated with TUV 5.0 at
z = 1.8 km.  Running this tool with `--lat 35 --lon 120 --date 2010-11-17 --alt 1.8
--o3-column 250` gives, at the solar-noon row, new / template: J(O1D) 0.96,
J(NO2) 0.98, J(O3P) 1.04, J(H2O2) 0.99, J(HNO3) 0.89, J(CH2O) 0.95,
J(acetone) 0.89, PAN 1.02.  Several spectra were revised between TUV 5.0 and
5.4 (for example HONO 0.76, MEK 0.64, glyoxal 1.2-1.3, HNO4 1.4, glycolaldehyde 2.1),
so those columns do not reproduce the old table.  The template's ozone column and surface
height are inferred from this comparison, not documented.
