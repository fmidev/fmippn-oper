#!/bin/bash
if [ ! $INTERPDIR ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

if [ ! -e $PPN_OUTPUT_FILE ]; then
   echo -e "Missing PPN output file $PPN_OUTPUT_FILE missing.\nRun the fmippn again, or copy the file from somewhere."
   exit 1
fi

cd $POSTPROCDIR
echo $TIMESTAMP $PPN_OUTPUT_FILE $LATEST_OBSFILE $INTERPDIR $INTERP_STEPS
find $INTERPDIR -name ${INTERP_NC_ACCPREF}'*' -exec rm -f {} \;
find $INTERPDIR -name '*.p?m' -exec rm -f {} \;
bin/thread_member_interp $TIMESTAMP $PPN_OUTPUT_FILE $LATEST_OBSFILE $INTERPDIR $DOMAIN $INTERP_STEPS

if [ $GENERATE_HOURMASKS ]; then
   MINS=${TIMESTAMP:10:2}
   if [ $MINS == 30 ] || [ $MINS == 00 ]; then
      LASTDBZ=`ls -1 $INTERPDIR/dBZ_M00_"$TIMESTAMP"*.pgm | tail -1`
      bin/hourmask $INTERPDIR `basename $LASTDBZ` mask1h
   fi
fi
