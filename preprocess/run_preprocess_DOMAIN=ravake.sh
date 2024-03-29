#!/bin/bash

set_BeginTime
echo "$BeginStamp : BEGIN=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP obsfile=$LATEST_OBSFILE"

if [ 0 == 1 ]; then

# sleep 60
PROJCFG=config/ETRStoSTERE.cfg
ANDREDIR=/radar/storage/HDF5/latest/radar/rack/comp
ANDRESTOR=/dev/shm/andre
TPGM=$ANDRESTOR/t.pgm
if [ ! -e $ANDRESTOR ]; then mkdir $ANDRESTOR ; fi

ANDRESEEK=_radar.rack.comp_CONF=FMIPPN,ANDRE.h5
for ANDREH5 in `ls -1rt $ANDREDIR/*"$ANDRESEEK" | tail -4` ; do
   TIMESTAMP=`basename $ANDREH5 | cut -c-12`
   ANDREPGM=/dev/shm/andre/"$TIMESTAMP"_fmi.radar.composite.lowest_FIN_RAVAKE.pgm
   echo $ANDREPGM
   if [ ! -e $ANDREPGM ]; then
       $POSTPROCDIR/bin/reprojection_radardata --cfgfile="$PROJCFG" --obstime="$TIMESTAMP" $ANDREH5 $TPGM
       HDRLINES=`grep -m 1 -n ^255$ $TPGM | cut -d: -f1`
       head -"$HDRLINES" $TPGM > hdr
       HDRSIZE=`stat -c %s hdr`
       PGMSIZE=`stat -c %s $TPGM`
       DATASIZE=`expr $PGMSIZE \- $HDRSIZE`
       convert -fill '#FEFEFE' -floodfill +0+0 white -fill black -opaque white -fill white -opaque '#FEFEFE' $TPGM conv.pgm  
       tail --bytes="$DATASIZE" conv.pgm > data
       grep ^# hdr > commhdr
       grep -v ^# hdr > defhdr
       echo P5 > newhdr
       echo '# metersperpixel_x 1000.0' >> newhdr
       echo '# metersperpixel_y 1000.0' >> newhdr
       cat commhdr >> newhdr
       tail -2 defhdr >> newhdr
       cat newhdr data > $ANDREPGM
   fi
done
# rm $TPGM


# RUN_SODIS=1
# if [ $RUN_SODIS ]; then
#     scp $LATEST_OBSFILE hhohti@fmippn-1.fmi.fi:/data/input/orig_ravake & 
#     for MACHINE in 1 2 3 ; do
#        SODISLOG=${LOGDIR}/run_sodis_${MACHINE}.log
#        ssh hhohti@fmippn-${MACHINE}.fmi.fi /fmi/dev/ppn/fmippn-oper/run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
#        ssh hhohti@fmippn-${MACHINE}.fmi.fi ./run_test.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
#     done
# fi

endif
get_Runtime
echo "$EndStamp : END=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime"
exit 0
