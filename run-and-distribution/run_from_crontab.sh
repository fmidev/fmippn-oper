#!/bin/bash
if [ -e /dev/shm/andre/TRIGGERED ]; then exit ; fi

date -u >> /var/tmp/ppn/log/CRONRUNS
~/fmippn-oper/run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> /var/tmp/ppn/log/run_from_cron_err.log 2>&1 &
exit

