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

#______________________________________________________________
# Configure 

export OPERDIR=${OPERDIR:-"$HOME/fmippn-oper"}
COMMONCONF=${OPERDIR}/config/set_common_config.sh
source $COMMONCONF  #activate conda env in the end

set_BeginTime
echo "$BeginStamp : BEGIN=fmippn domain=${DOMAIN} timestamp=$TIMESTAMP" >> $RUNLOG
fmippn_BeginTime=$BeginTime # store the begin time of whole fmippn process

#______________________________________________________________
# Preprocess

set_BeginTime
echo "$BeginStamp : BEGIN=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP" >> $RUNLOG
cd $PREPROCDIR
PREPROC=run_preprocess_DOMAIN=${DOMAIN}.sh
if [ -e $PREPROC ]; then
#   . ./$PREPROC  >> $PREPROCLOG 2>&1
   if [ $? == $RECONFIG ]; then # preprocess insisted reconfiguration
       source $COMMONCONF  
   fi
fi
get_Runtime
echo "$EndStamp : END=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime" >> $RUNLOG

#______________________________________________________________
# Run The PPN

set_BeginTime
echo "$BeginStamp : BEGIN=ppn domain=${DOMAIN} timestamp=$TIMESTAMP" >> $RUNLOG
cd $PPNDIR
# conda info -e
$PYTHON run_ppn.py --timestamp=${TIMESTAMP} --config=${DOMAIN}  >> $PPNLOG 2>&1
get_Runtime
echo "$EndStamp : END=ppn domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime" >> $RUNLOG

#______________________________________________________________
# Postprocess

set_BeginTime
echo "$BeginStamp : BEGIN=postprocess domain=${DOMAIN} timestamp=$TIMESTAMP" >> $RUNLOG
cd $POSTPROCDIR
./run_postprocesses.sh >> $RUNLOG
get_Runtime
echo "$EndStamp : END=postprocess domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime" >> $RUNLOG


#______________________________________________________________
# Distribution

cd $RUNDIR
DISTRIBUTE=distribute_DOMAIN=${DOMAIN}.sh
if [ -e $DISTRIBUTE ]; then
   set_BeginTime
   echo "$BeginStamp : BEGIN=distribution domain=${DOMAIN} timestamp=$TIMESTAMP" >> $RUNLOG
   ./$DISTRIBUTE >> $DISTRIBLOG 2>&1 
   get_Runtime
   echo "$EndStamp : END=distribution domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime" >> $RUNLOG
fi

#______________________________________________________________
# Cleanup

find $PPN_OUTPUT_DIR -name 'nc_????????????.h5' -mmin +15 -exec rm -f {} \;
find $OBSDIR -type f -mtime +1 -exec rm -f {} \;
conda deactivate
BeginTime=$fmippn_BeginTime
get_Runtime
echo "$EndStamp : END=fmippn domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime" >> $RUNLOG
exit
