#!/bin/bash
echo $* >> ~/fmippn-oper/log/trig.log
~/fmippn-oper/run-and-distribution/run_fmippn_common.sh --DOMAIN=ravake >> $HOME/fmippn-oper/log/trig_opertest.log 2>&1 &
exit
