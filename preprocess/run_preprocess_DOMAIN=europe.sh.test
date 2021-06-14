#!/bin/bash

EURDIR_SRC=/mnt/meru/data/prod/radman/latest/fmi/radar/composite/Europe
DSTPATTERN=`echo "$FILEPATTERN" | tr -d \*`
echo filepattern:"$FILEPATTERN"
echo $DSTPATTERN

for SRC_PGM in `ls -1rt $EURDIR_SRC/*.pgm | tail -4` ; do
   TIMESTAMP=`basename $SRC_PGM | cut -c-12`
   DST_PGM=$OBSDIR/"$TIMESTAMP""$DSTPATTERN"
   if [ ! -e DST_PGM ]; then
      if [ ! $DATABYTES ]; then
         DIMSLINE=`pamfile $SRC_PGM | cut -d, -f2- | cut -d' ' -f-2,4`
         DATABYTES=`echo $DIMSLINE | tr ' ' \* | bc`
      fi
      echo $DST_PGM $DATABYTES
      cp $CONFDIR/europe_PGM_stere.hdr $DST_PGM
      echo -e $DIMSLINE"\n"255 >> $DST_PGM
      tail --bytes="$DATABYTES" $SRC_PGM >> $DST_PGM
   fi
done
exit $RECONFIG
