#!/bin/bash
if [ $1 ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
fi
source $HOME/fmippn-oper/config/set_common_config.sh
INPUT_PATH=${INPUT_PATH:-"/mnt/meru/data/prod/radman/latest/fmi/radar/composite/lowest"}
FILEPATTERN=${FILEPATTERN:-'*_fmi.radar.composite.lowest_FIN_RAVAKE.pgm'}
PGMNAME=${FILEPATTERN:1}

export ZR_AF=223.0
export ZR_BF=1.53
export ZR_AC=223.0
export ZR_BC=1.53
export ACCMINS=5

$RUNDIR/bin/acc1h $TIMESTAMP $INPUT_PATH $PGMNAME $DISTRIBDIR/ObsAcc1h_
