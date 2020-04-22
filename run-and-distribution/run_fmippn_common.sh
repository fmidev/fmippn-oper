#!/bin/bash
if [ $1 ]; then
   if [ `echo $1 | grep ^\-\-DOMAIN=` ]; then
      export DOMAIN=`echo $1 | cut -d= -f2`
   fi
fi

if [ ! $DOMAIN ]; then
      echo "No domain defined!"
      echo "Set DOMAIN environmet variable, or use option --DOMAIN=domain"
      exit 1
fi

source $HOME/fmippn-oper/config/set_common_config.sh

cd $PREPROCDIR
./run_preprocesses.sh  >> $PREPROCLOG 2>&1

cd $PPNDIR
echo "PPN processing for domain \"${DOMAIN}\" at $TIMESTAMP to output $PPN_OUTPUT_FILE started at"
date

conda info -e

echo "PPN processing for domain \"${DOMAIN}\" of $TIMESTAMP started at"
date
$PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}  # >> $PPNLOG 2>&1

cd $POSTPROCDIR
./run_postprocesses.sh  >> $POSTPROCLOG 2>&1
./run_distribution.sh >> $DISTRIBLOG 2>&1
echo "PPN processing for domain \"${DOMAIN}\" of $TIMESTAMP ended at"
date
echo
conda deactivate
exit
