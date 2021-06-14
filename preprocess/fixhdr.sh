#!/bin/bash
grep ^# hdr > commhdr
COMMLINES=`cat commhdr | wc -l`
TAILLINES=`expr $COMMLINES \- 2`
head -2 commhdr > commhead
tail -"$TAILLINES" commhdr > commtail
grep -v ^# hdr > defhdr
echo P5 > newhdr
cat commhead >> newhdr
echo '# metersperpixel_x 1000.0' >> newhdr
echo '# metersperpixel_y 1000.0' >> newhdr
cat commtail >> newhdr
tail -2 defhdr >> newhdr
exit
