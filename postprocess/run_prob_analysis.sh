#!/bin/bash

# MEMBERS=15
# TIMESTEPS=12

if [ ! $PROBDIR ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

if [ ! -e $PPN_OUTPUT_FILE ]; then
   echo -e "Missing PPN output file $PPN_OUTPUT_FILE.\nRun the fmippn again, or copy the file from somewhere."
   exit 1
fi

MEMBERS=`h5dump -a /meta/configuration/ENSEMBLE_SIZE $PPN_OUTPUT_FILE | grep : | cut -d: -f2 | tr -d \"`
TIMESTEPS=`h5dump -a /meta/configuration/NUM_TIMESTEPS $PPN_OUTPUT_FILE | grep : | cut -d: -f2 | tr -d \"`
STORE_DETERM=`h5dump -a /meta/configuration/STORE_DETERMINISTIC $PPN_OUTPUT_FILE | grep : | cut -d: -f2 | tr -d \"' ' | tr '[:lower:]' '[:upper:]'`

echo "$STORE_DETERM $MEMBERS"

# If deterministic member is stored, it's always written to interpolated accumulation as member #0
if [ "$STORE_DETERM" == "TRUE" ]; then
   MEMBERS=`expr $MEMBERS + 1`
else
   # If deterministic member should be ignored in prob analysis, but data doesn't contain determ member
   export PROB_IGNORE_DETERM=FALSE
fi
cd $POSTPROCDIR

echo $TIMESTAMP $PROB_THRCFG $INTERPDIR $PROBDIR $DOMAIN $MEMBERS $TIMESTEPS
bin/prob_thresholding $TIMESTAMP $PROB_THRCFG $INTERPDIR $PROBDIR $DOMAIN $MEMBERS $TIMESTEPS $PROB_TIMESTEP
find $PROBDIR -name '*.*' -mmin +20 -exec rm -f {} \;
date
exit
