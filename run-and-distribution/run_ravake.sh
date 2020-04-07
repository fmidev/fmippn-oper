#!/bin/bash
export DOMAIN=ravake
source $HOME/ppn/config/set_common_config.sh

cd $PPNDIR
echo $OBSFILE
echo "PPN process for domain \"${DOMAIN}\" of $TIMESTAMP to output $PPN_OUTPUT_FILE started at"
date

conda info -e
conda deactivate
exit

$PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}  >> $PPNLOG
date

cd $POSTPROCDIR
./postprocess_ppndata.tcsh $PPN_OUTPUT_FILE >> $POSTPROCLOG 2>&1
echo "PPN process for $DOMAIN of $TIMESTAMP ended at"
date
echo
conda deactivate
exit
