#!/bin/bash
# Extract motion field and archive

YY=${TIMESTAMP:0:4}
MM=${TIMESTAMP:4:2}
DD=${TIMESTAMP:6:2}
ARCHDIR=/arch/radar/raw/$YY/$MM/$DD/vector/
MOTIONFILE=${PPN_OUTPUT_DIR}/${TIMESTAMP}_${DOMAIN}_motion.h5
bin/extract_motion $PPN_OUTPUT_FILE $MOTIONFILE
scp $MOTIONFILE radman@esteri.fmi.fi:/radar/storage/HDF5/PPN/
rm $MOTIONFILE


