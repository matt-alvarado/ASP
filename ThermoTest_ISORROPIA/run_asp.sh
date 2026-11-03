#!/bin/bash
# run_asp.sh <tag> <ASP binary> [extra environment, e.g. ASP_FEE_OUTER=worst]
# Runs the ThermoTest in every run directory made by "python3 mkcase.py <tag>".
# Tags used by analysis.py: p34 (down scan), up (ASP_THERMO_SCAN=up; decks start at RH 0.20),
# worst (ASP_FEE_OUTER=worst), base (pre-Claude build).
HERE=$(cd $(dirname $0) && pwd); WORK=${THERMO_WORK:-$HERE/work}
tag=$1; bin=$2; shift 2
for d in $WORK/runs/$tag/*; do
  ( cd $d; env "$@" timeout 600 $bin > run.log 2> err.log; echo "$tag $(basename $d) rc=$?" >> $WORK/runs/log_$tag.txt )
done
