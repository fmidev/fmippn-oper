#!/bin/bash

if [ ! $DISTRIBDIR ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

cd $INTERPDIR

MINS=${TIMESTAMP:10:2}

if [ $(($MINS % 30)) == 0 ]; then
   for EMINS in 60 120 ; do
      EM=`expr $MINS + $EMINS`
      SM=`expr $EM \- 60`
      if [ $EM -lt 100 ]; then EM=0"$EM" ; fi
      if [ $SM -lt 100 ]; then SM=0"$SM" ; fi

      EPGM=`ls -1 EnsmeanAcc_"$TIMESTAMP"-*+"$EM"_"$DOMAIN".pgm 2> /dev/null`
      if [ "$EPGM" == "" ]; then continue ; fi

      TIMES=`echo $EPGM | cut -d_ -f2 | cut -d+ -f1`
      ENDTIME=${TIMES:13:12}
      ACCPGM="$TIMES"_acc1h_suomi.pgm
      ACCTIF=${DISTRIBDIR}/${ACCPGM%.*}.tif

      SPGM=`ls -1 EnsmeanAcc_"$TIMESTAMP"-*+"$SM"_"$DOMAIN".pgm 2> /dev/null`
      if [ "$SPGM" != "" ]; then
         pamarith -subtract $EPGM $SPGM > $ACCPGM
      else
         ACCPGM=$EPGM
      fi

      $POSTPROCDIR/bin/reprojection_radardata --cfgfile=$CONFDIR/RAVAKE_3067.cfg --obstime=$ENDTIME $ACCPGM $ACCTIF
   
      if [ -e $ACCTIF ]; then
         date -u
         echo "Copying $ACCTIF"
         cp $ACCTIF /mnt/meru/data/prod/radman/ppn/tif
      fi
   done
   cd $DISTRIBDIR
   find . -name '*.tif' -mmin +1500 -exec rm -f {} \;
fi

exit
