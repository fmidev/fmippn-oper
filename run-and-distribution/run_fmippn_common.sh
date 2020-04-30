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

source $HOME/fmippn-oper/config/set_common_config.sh  #activate conda env in the end

cd $PREPROCDIR
./run_preprocesses.sh  >> $PREPROCLOG 2>&1

cd $PPNDIR
# conda info -e

echo "PPN processing for domain \"${DOMAIN}\" of $TIMESTAMP started at"
date
cmd="$PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}"
echo $cmd
$PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}  >> $PPNLOG 2>&1

cd $POSTPROCDIR
echo "Postprocessing started at"
date
./run_postprocesses.sh  >> $POSTPROCLOG 2>&1

cd $RUNDIR
DISTRIBUTE=distribute_DOMAIN=${DOMAIN}.sh
if [ -e $DISTRIBUTE ]; then
   echo "Distribution started at"
   date
   ./$DISTRIBUTE >> $DISTRIBLOG 2>&1 
fi

find $PPN_OUTPUT_DIR -name 'nc_????????????.h5' -mmin +15 -exec rm -f {} \;
echo "PPN processing for domain \"${DOMAIN}\" of $TIMESTAMP ended at"
date
echo
conda deactivate
exit
