#!/bin/bash
# Regenerate PhotoRates.in for the four case decks.  Usage: TUV_DIR=/path/to/V5.4 ./make_case_photorates.sh
# (Otavi 19.7 S 17.6 E; Alaska 63.46 N 158.65 W; Timbavati 24.38 S 31.25 E; Williams 34.696 N 120.206 W)
# Site, date, albedo and ozone column: Alvarado & Prinn (2009) Table 1 (Otavi, Alaska, Timbavati) and
# Alvarado et al. (2015) Sect. 4-5 (Williams).  J-values are taken at the altitude of the paper's TUV
# run (Williams 1.6 km, Otavi 3 km, Alaska 2 km, Timbavati 1 km asl), at the middle (0.5 km) of a 1 km
# PBL whose floor is therefore 0.5 km lower (Williams: plume 1.1-2.1 km).  PBL AOD(330 nm) = 0.5;
# PBL O3, NO2, SO2 = the background mixing ratios in each deck's GasPhaseChems.in; SSA 0.9.
set -e
D=$(cd "$(dirname "$0")" && pwd); R="$D/../InorganicInputDecks"; T=${TUV_DIR:?set TUV_DIR}
run() { c=$1; shift; python3 "$D/make_photorates.py" --tuv "$T" --pbl-height 1.0 --zout 0.5 \
        --tau-aer 0.5 --wl-aer 330 --ssa-uv 0.9 --ssa-vis 0.9 --workdir "${WORKDIR:-$D/tuv_work}" \
        -o "$R/$c/PhotoRates.in" "$@"; }
run Williams   --lat 34.6958 --lon -120.2064 --date 2009-11-17 --alt 1.1 --albedo 0.035 --o3-column 278 --o3-ppb 50   --no2-ppb 0.039  --so2-ppb 0.120
run Otavi      --lat -19.7   --lon 17.6      --date 2000-09-13 --alt 2.5 --albedo 0.15  --o3-column 349.13 --o3-ppb 85 --no2-ppb 1.3    --so2-ppb 0
run Alaska     --lat 63.46   --lon -158.65   --date 1997-06-27 --alt 1.5 --albedo 0.10  --o3-column 320    --o3-ppb 32.5 --no2-ppb 0.0166 --so2-ppb 0
run Timbavati  --lat -24.38  --lon 31.25     --date 2000-09-07 --alt 0.5 --albedo 0.15  --o3-column 288    --o3-ppb 0    --no2-ppb 0.46   --so2-ppb 10
