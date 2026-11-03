#!/bin/bash
# solver_diag.sh <tag> <ASP binary> [up]: for every run directory of <tag>, dump every electrolyte solve of the
# ThermoTest scan (patch 27 harness), then replay the dump with the coupled solve and with the worst-first loop.
# Needs a replay directory $WORK/rp holding the default deck in Lagrangian mode (made here if missing).
HERE=$(cd $(dirname $0) && pwd); WORK=${THERMO_WORK:-$HERE/work}; DECKS=${ASP_DECKS:-$HERE/../InorganicInputDecks}
tag=$1; B=$2; scan=$3
if [ ! -d $WORK/rp ]; then mkdir -p $WORK/rp/OutputFiles $WORK/rp/InorganicInputDecks; cp $DECKS/*.in $WORK/rp/InorganicInputDecks/; fi
for d in $WORK/runs/$tag/*_*; do
  rm -f $d/fee.bin
  ( cd $d; if [ "$scan" = up ]; then export ASP_THERMO_SCAN=up; fi
    ASP_FEE_DUMP_ALL=1 ASP_FEE_DUMP_EVERY=1 ASP_FEE_DUMP=$PWD/fee.bin timeout 300 $B > /dev/null 2>&1 )
  ( cd $WORK/rp
    ASP_FEE_REPLAY=$d/fee.bin timeout 300 $B > run.log 2> err.log; cp FEEReplay.txt $d/FEEReplay_cur.txt 2>/dev/null
    ASP_FEE_OUTER=worst ASP_FEE_REPLAY=$d/fee.bin timeout 300 $B > run.log 2> err.log; cp FEEReplay.txt $d/FEEReplay_worst.txt 2>/dev/null )
  rm -f $d/fee.bin
done
