#!/bin/bash
export DOMAIN=test
source $HOME/ppn/config/set_common_config.sh

cd $PPNDIR
echo "PPN process for domain \"${DOMAIN}\" at $TIMESTAMP to output $PPN_OUTPUT_FILE started at"
date

conda info -e

echo "PPN processing for domain \"${DOMAIN}\" of $TIMESTAMP started at"
date
# $PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}  # >> $PPNLOG 2>&1

cd $POSTPROCDIR
./run_postprocesses.sh # >> $POSTPROCLOG 2>&1
echo "PPN processes for domain \"${DOMAIN}\" of $TIMESTAMP ended at"
date
echo
conda deactivate
exit
