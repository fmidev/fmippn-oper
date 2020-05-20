#!/bin/bash
export DOMAIN=${DOMAIN:- 'test'}
if [ ! $CONFDIR ]; then 
   export PROCTIME=${PROCTIME:-"$1"}
   shift
   source $HOME/fmippn-oper/config/set_common_config.sh
fi

source $CONFDIR/common_functions.sh
LOGHOUR=${TIMESTAMP:8:2}
cd $POSTPROCDIR

for PROCESS in interpolation ensemble_mean prob_analysis visualization ; do
   RUNSCRIPT=run_${PROCESS}.sh
   DOMAINRUN=RUN_${PROCESS}
   if [ "${!DOMAINRUN}" == "TRUE" ]; then
      set_BeginTime
      echo "$BeginStamp : BEGIN=$PROCESS domain=${DOMAIN} timestamp=$TIMESTAMP"
      ./${RUNSCRIPT} >> $LOGDIR/${PROCESS}_DOMAIN=${DOMAIN}_${LOGHOUR}.log 2>&1
      get_Runtime
      echo "$EndStamp : END=$PROCESS domain=${DOMAIN} timestamp=$TIMESTAMP runtime=$Runtime"
#   else
#      echo "# No \"$PROCESS\" process defined for domain $DOMAIN"
   fi
done

