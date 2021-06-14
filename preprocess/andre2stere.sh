#!/bin/bash
PROJCFG=config/ETRStoSTERE.cfg
ANDREDIR=/radar/storage/HDF5/latest/radar/rack/comp
ANDRESTOR=/dev/shm/andre
TPGM=$ANDRESTOR/t.pgm
if [ ! -e $ANDRESTOR ]; then mkdir $ANDRESTOR ; fi

ANDRESEEK=_radar.rack.comp_CONF=FMIPPN,ANDRE.h5
for ANDREH5 in `ls -1t $ANDREDIR/*"$ANDRESEEK" | tail -4` ; do
   TIMESTAMP=`basename $ANDREH5 | cut -c-12`
   ANDREPGM=/dev/shm/andre/"$TIMESTAMP"_fmi.radar.composite.lowest_FIN_RAVAKE.pgm
   echo $ANDREPGM
   if [ ! -e $ANDREPGM ]; then
       $POSTPROCDIR/bin/reprojection_radardata --cfgfile="$PROJCFG" --obstime="$TIMESTAMP" $ANDREH5 $TPGM
       convert -fill '#FEFEFE' -floodfill +0+0 white -fill black -opaque white -fill white -opaque '#FEFEFE' $TPGM $ANDREPGM 
   fi
done
rm $TPGM
 
