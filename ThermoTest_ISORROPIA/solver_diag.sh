#!/bin/bash
# solver_diag.sh <tag> <asp_feedump> <FEEReplay>
# For every run directory of <tag>, dumps every electrolyte solve of the ThermoTest scan with asp_feedump, then replays
# the dump with FEEReplay twice: with the coupled solve (FEEReplay_cur.txt) and with the worst-first loop
# (FEEReplay_worst.txt). Build asp_feedump and FEEReplay in src/Tests (see src/Tests/README.md).
# The scan direction comes from each run directory's deck (mkcase.py). The replays run in $WORK/rp (default deck) and
# $WORK/rp_worst (default deck plus AqOuterSolver = 0), made here if missing.
HERE=$(cd $(dirname $0) && pwd); WORK=${THERMO_WORK:-$HERE/work}; DECKS=${ASP_DECKS:-$HERE/../InorganicInputDecks}
tag=$1; DUMP=$2; REPLAY=$3

# make_replay_dir <directory> [deck setting]: a run directory with the default deck and, optionally, one setting
make_replay_dir() {
  [ -d $1 ] && return
  mkdir -p $1/OutputFiles $1/InorganicInputDecks; cp $DECKS/*.in $1/InorganicInputDecks/
  if [ -n "$2" ]; then printf '%s\r\n' "$2" >> $1/InorganicInputDecks/ASPInputDeck.in; fi
}
make_replay_dir $WORK/rp
make_replay_dir $WORK/rp_worst 'AqOuterSolver = 0'

for d in $WORK/runs/$tag/*_*; do
  rm -f $d/fee.bin
  ( cd $d; ASP_FEE_DUMP_ALL=1 ASP_FEE_DUMP_EVERY=1 ASP_FEE_DUMP=$PWD/fee.bin timeout 300 $DUMP > /dev/null 2>&1 )
  ( cd $WORK/rp;       timeout 300 $REPLAY $d/fee.bin > run.log 2> err.log; cp FEEReplay.txt $d/FEEReplay_cur.txt 2>/dev/null )
  ( cd $WORK/rp_worst; timeout 300 $REPLAY $d/fee.bin > run.log 2> err.log; cp FEEReplay.txt $d/FEEReplay_worst.txt 2>/dev/null )
  rm -f $d/fee.bin
done
