#!/bin/bash

if [ ! $DISTRIBDIR ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

if [ "$DOMAIN" == "test" ]; then
   ENSMEAN=`ls -1 $INTERPDIR/Ensmean_${TIMESTAMP}-*+060_${DOMAIN}.pgm`
   OUTTIF=$INTERPDIR/${TIMESTAMP}+acc1h_${DOMAIN}.tif
   $POSTPROCDIR/bin/reprojection_radardata --verbose --cfgfile=$CONFDIR/RAVAKE_3067.cfg --obstime=$TIMESTAMP $ENSMEAN $OUTTIF
fi
exit
