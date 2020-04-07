#!/bin/bash
# export DOMAIN=${DOMAIN}
if [ ! $ENSMEANDIR ]; then
   echo setting
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

cd $ENSMEANDIR
$POSTPROCDIR/bin/ensmean_determweighted $PPN_OUTPUT_FILE
exit
