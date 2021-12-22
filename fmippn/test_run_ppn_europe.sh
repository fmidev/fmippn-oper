#!/bin/bash

cmd="conda run -n ppn_pysteps14 python run_ppn.py --config=europe --timestamp=202108190815"
echo $cmd
eval $cmd
