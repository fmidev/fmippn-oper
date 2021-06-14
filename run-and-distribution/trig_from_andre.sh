#!/bin/bash
export LD_LIBRARY_PATH=/var/opt/lib
cd ~/fmippn-oper/preprocess

# if [ "$1" == "" ]; then
# CRON=1
# fi
  
LOGDIR=/var/tmp/ppn/log
LOG=$LOGDIR/trig_andre.log
date -u >> $LOG
echo $* >> $LOG
echo $1
# exit

LATEST_TIMESTAMP=$1
MINS=`echo $LATEST_TIMESTAMP | cut -c11-12`
Q=`expr $MINS \% 15`

PROJCFG=config/ETRStoSTERE.cfg
ANDREDIR=/radar/storage/HDF5/latest/radar/rack/comp
ANDRESTOR=/dev/shm/andre
TRIGGER=$ANDRESTOR/TRIGGERED
TPGM=$ANDRESTOR/t.pgm
if [ ! -e $ANDRESTOR ]; then mkdir $ANDRESTOR ; fi

ANDRESEEK=_radar.rack.comp_CONF=FMIPPN,ANDRE.h5
for ANDREH5 in `ls -1rt $ANDREDIR/*"$ANDRESEEK" | tail -4` ; do
   TIMESTAMP=`basename $ANDREH5 | cut -c-12`
   ANDREPGM=/dev/shm/andre/"$TIMESTAMP"_fmi.radar.composite.lowest_FIN_RAVAKE.pgm
   echo $ANDREPGM >> $LOG 

# Tiedoston olemassaoloa ei tutkita, sillä backup-data tulee suoraan skriptiltä radman@radprod:/fmi/prod/pipe/bin/generate_lowest_composite.tcsh
#   if [ ! -e $ANDREPGM ]; then 
       echo $TIMESTAMP >> nodes.log
       h5dump -a /how/nodes $ANDREH5 | grep : >> nodes.log
       ../postprocess/bin/reprojection_radardata --cfgfile="$PROJCFG" --obstime="$TIMESTAMP" $ANDREH5 $TPGM >& /dev/null
       HDRLINES=`grep -m 1 -n ^255$ $TPGM | cut -d: -f1`
       head -"$HDRLINES" $TPGM > hdr
       HDRSIZE=`stat -c %s hdr`
       PGMSIZE=`stat -c %s $TPGM`
       DATASIZE=`expr $PGMSIZE \- $HDRSIZE`
       convert -fill '#FEFEFE' -floodfill +0+0 white -fill black -opaque white -fill white -opaque '#FEFEFE' $TPGM conv.pgm  
       tail --bytes="$DATASIZE" conv.pgm > data
       grep ^# hdr > commhdr
       COMMLINES=`cat commhdr | wc -l`
       TAILLINES=`expr $COMMLINES \- 2`
       head -2 commhdr > commhead
       tail -"$TAILLINES" commhdr > commtail
       echo P5 > newhdr
       cat commhead >> newhdr
       echo '# metersperpixel_x 1000.0' >> newhdr
       echo '# metersperpixel_y 1000.0' >> newhdr
       cat commtail >> newhdr
       tail -2 hdr >> newhdr
       cat newhdr data > "$ANDREPGM".tmp
       mv -f "$ANDREPGM".tmp $ANDREPGM
       convert -geometry 25% $ANDREPGM ~/test/small/"$MINS".png
#   fi
done

if [ $Q == 0 ]; then
  echo RUN $TIMESTAMP
  echo $1 > $TRIGGER
  date -u >> $TRIGGER
  ( ../run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> /var/tmp/ppn/log/trig_andre_err.log 2>&1 ; rm $TRIGGER ) &
fi
exit
