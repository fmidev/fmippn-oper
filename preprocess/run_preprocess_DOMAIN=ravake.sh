#!/bin/bash

set_BeginTime
echo "$BeginStamp : BEGIN=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP obsfile=$LATEST_OBSFILE"
RUN_SODIS=1
if [ $RUN_SODIS ]; then
     scp $LATEST_OBSFILE hhohti@fmippn-1.fmi.fi:/data/input/orig_ravake & 
#     for MACHINE in 1 2 3 ; do
#        SODISLOG=${LOGDIR}/run_sodis_${MACHINE}.log
#        ssh hhohti@fmippn-${MACHINE}.fmi.fi /fmi/dev/ppn/fmippn-oper/run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
#        ssh hhohti@fmippn-${MACHINE}.fmi.fi ./run_test.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
#     done
fi
get_Runtime
echo "$EndStamp : END=preprocess domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime"
exit 0
