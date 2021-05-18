#!/bin/bash

cmd="conda run -n ppn_pysteps14 python run_ppn_pot.py --config=pot --timestamp=202105171200"
echo $cmd
eval $cmd
