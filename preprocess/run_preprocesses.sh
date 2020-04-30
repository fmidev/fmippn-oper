#!/bin/bash

if [ "$DOMAIN" == "ravake" ]; then
   echo "preprocessing $LATEST_OBSFILE for domain \"$DOMAIN\""
   if [ $RUN_SODIS ]; then
#     scp $LATEST_OBSFILE hhohti@fmippn-1.fmi.fi:/mnt/data/input/$DOMAIN 
     for MACHINE in 1 2 3 ; do
        SODISLOG=${LOGDIR}/run_sodis_${MACHINE}.log
#        ssh hhohti@fmippn-${MACHINE}.fmi.fi /fmi/dev/ppn/fmippn-oper/run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
        ssh hhohti@fmippn-${MACHINE}.fmi.fi ./run_test.sh --DOMAIN=ravake >> $SODISLOG 2>&1 &     
     done
   fi
fi
exit