#!/bin/bash

if [ ! $DISTRIBDIR ]; then
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source ~/fmippn-oper/config/set_common_config.sh
fi

cd $INTERPDIR

MINS=${TIMESTAMP:10:2}
DISTRIB_NFS=/mnt/meru/data/prod/radman/ppn/tif

# if [ $(($MINS % 30)) == 0 ]; then
if [ $(($MINS % 60)) == 0 ]; then
   # observed accumulation of previous hour
   if [ "$MINS" == "00" ]; then 
      INPUTNAME=${FILEPATTERN:1}
      export ZR_AF=223.0
      export ZR_BF=1.53
      export ZR_AC=223.0
      export ZR_BC=1.43
      export ACCMINS=5
      OBSACC=$DISTRIBDIR/"$TIMESTAMP"_obsacc1h_suomi.pgm
      OBSACCTIF=${OBSACC%.*}.tif
      OBSACCH5=${OBSACC%.*}.h5
      echo "Creating observed accumulation for previous hour for $TIMESTAMP"
      $RUNDIR/bin/acc1h $TIMESTAMP $INPUT_PATH $INPUTNAME $OBSACC 
      $POSTPROCDIR/bin/reprojection_radardata --cfgfile=$CONFDIR/RAVAKE_3067.cfg --NoDataValue=65535 --obstime=$TIMESTAMP $OBSACC $OBSACCTIF $OBSACCH5
      echo "Copying $OBSACCTIF"
      cp $OBSACCTIF $DISTRIB_NFS
   fi

   # Nowcasted 1h accumulations of ensemble mean 
   for EMINS in 60 120 180 240 ; do
      EM=`expr $MINS + $EMINS`
      SM=`expr $EM \- 60`
      if [ $EM -lt 100 ]; then EM=0"$EM" ; fi
      if [ $SM -lt 100 ]; then SM=0"$SM" ; fi

      EPGM=`ls -1 EnsmeanAcc_"$TIMESTAMP"-*+"$EM"_"$DOMAIN".pgm 2> /dev/null`
      if [ "$EPGM" == "" ]; then continue ; fi

      TIMES=`echo $EPGM | cut -d_ -f2 | cut -d+ -f1`
      ENDTIME=${TIMES:13:12}
      ACCPGM="$TIMES"_acc1h_suomi.pgm
      ACCPGM_TMP="$TIMES"_acc1h_suomi_tmp.pgm
      ACCTIF=${DISTRIBDIR}/${ACCPGM%.*}.tif
      HSYACCTIF_ORIG=${DISTRIBDIR}/${TIMES}_PPN.tif
      HSYACCTIF=${DISTRIBDIR}/${TIMES#*-}_PPN.tif
      ACCH5=${DISTRIBDIR}/${ACCPGM%.*}.h5

      SPGM=`ls -1 EnsmeanAcc_"$TIMESTAMP"-*+"$SM"_"$DOMAIN".pgm 2> /dev/null`
      if [ "$SPGM" != "" ]; then
         pamarith -subtract $EPGM $SPGM > $ACCPGM
      else
         ACCPGM=$EPGM
      fi

      # Adding hourly radar coverage mask (see $POSTPROCDIR/run_interpolation)
      ACCPGM_MASKED="$TIMES"_acc1h_suomi_masked.pgm
      composite -compose plus mask1h_"$TIMES".pgm $ACCPGM $ACCPGM_TMP
      $RUNDIR/bin/replace_grey $ACCPGM_TMP $ACCPGM_MASKED 19270 65000 0

      $POSTPROCDIR/bin/reprojection_radardata --cfgfile=$CONFDIR/RAVAKE_3067.cfg --NoDataValue=65535 --obstime=$ENDTIME $ACCPGM_MASKED $ACCTIF $ACCH5

# HSYn pääkaupunkiseutu
      $POSTPROCDIR/bin/reprojection_radardata --cfgfile=$CONFDIR/HSY_3067.cfg --NoDataValue=65535 --obstime=$ENDTIME $ACCPGM_MASKED $HSYACCTIF_ORIG
#     scp $HSYACCTIF_ORIG fmi@dev.elmo.fmi.fi:"/smartmet/data/tutka/suomi/ppn/3067/tif/suomi_ppn_eureffin/${HSYACCTIF}" &

      if [ -e $ACCTIF ]; then
         echo "Copying $ACCTIF"
#         cp $ACCTIF $DISTRIB_NFS 
      fi
   done

   # Adds phase fields and copies corrected accumululations 
   ~/phase/bin/potprec_blur_ravake.tcsh $TIMESTAMP $DISTRIB_NFS

   cd $DISTRIBDIR
   find . -name '*.tif' -mmin +1500 -exec rm -f {} \;
   find . -name '*.pgm' -mmin +1500 -exec rm -f {} \;
   find . -name '*.h5' -mmin +1500 -exec rm -f {} \;
   cd $DISTRIB_NFS
   find . -name '*.tif' -mtime +30 -exec rm -f {} \;
fi
