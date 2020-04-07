#!/bin/bash
export DOMAIN=${DOMAIN:- 'test'}
if [ ! $WORKDIR ]; then source $HOME/fmippn-oper/config/set_common_config.sh ; fi
cd $POSTPROCDIR

for PROCESS in interpolation ensemble_mean prob_analysis visualization ; do
   RUNSCRIPT=run_${PROCESS}.sh
   DOMAINRUN=RUN_${PROCESS}
   if [ "${!DOMAINRUN}" == "TRUE" ]; then
      echo "$PROCESS for domain \"$DOMAIN\" started at"
      date
      ./${RUNSCRIPT} >> $LOGDIR/${PROCESS}_DOMAIN=${DOMAIN}.log 2>&1
      echo "$PROCESS for domain \"$DOMAIN\" ready at"
      date
   else
      echo "No \"$PROCESS\" process defined for domain \"$DOMAIN\""
   fi
done
exit
