#!/bin/bash

cmd="conda run -n ppn_pysteps14 python run_ppn.py --config=ravake --timestamp=202011191050"
echo $cmd
eval $cmd
