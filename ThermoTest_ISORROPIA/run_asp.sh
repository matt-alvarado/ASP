#!/bin/bash
# run_asp.sh <tag> <ASP binary>
# Runs the ThermoTest in every run directory made by "python3 mkcase.py <tag> [Name=value ...]".
# The RH scan direction and the solver choices are deck settings that mkcase.py writes into each run directory.
# Tags used by analysis.py (config.py TAG_DOWN, TAG_UP): down (RH scanned down), up (RH scanned up from 0.20),
# worst (made with AqOuterSolver=0, the worst-first electrolyte loop), base (pre-Claude build, made with no settings).
HERE=$(cd $(dirname $0) && pwd); WORK=${THERMO_WORK:-$HERE/work}
tag=$1; bin=$2
for d in $WORK/runs/$tag/*; do
  ( cd $d; timeout 600 $bin > run.log 2> err.log; echo "$tag $(basename $d) rc=$?" >> $WORK/runs/log_$tag.txt )
done
