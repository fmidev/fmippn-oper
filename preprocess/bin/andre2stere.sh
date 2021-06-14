#!/bin/bash
export LD_LIBRARY_PATH=/var/opt/lib
export TZ=UTC

CFG=ETRStoSTERE.cfg
ANDREH5=/radar/storage/HDF5/latest/radar/rack/comp/202104140930_radar.rack.comp_CONF=FMIPPN,ANDRE.h5
TIMESTAMP=`basename $ANDREH5 | cut -c-12`
./reprojection_radardata --verbose --cfgfile="$CFG" --obstime="$TIMESTAMP" $ANDREH5 "$TIMESTAMP"_ANDRE.pgm
exit

 
