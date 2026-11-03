#!/bin/bash
# build_iso.sh <ISORROPIA-II checkout>: compile ISORROPIA II 2.3 with iso_driver.f into isobuild/iso
# ISORROPIA II: https://github.com/EPFL-LAPI/ISORROPIA-II (academic, non-commercial licence; not redistributed here)
set -e
HERE=$(cd $(dirname $0) && pwd); SRC=$1/iso_2_3Src
mkdir -p $HERE/isobuild; cd $HERE/isobuild
cp $SRC/isocom.for $SRC/isofwd.FOR $SRC/isorev.for $SRC/isrpia.inc .
gfortran -O2 -std=legacy -w -fno-automatic -o iso $HERE/iso_driver.f isocom.for isofwd.FOR isorev.for
echo built $HERE/isobuild/iso
