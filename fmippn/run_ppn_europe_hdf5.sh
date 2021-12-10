#!/bin/bash

cmd="conda run -n ppn_pysteps14 python run_ppn.py --config=europe_hdf5 --timestamp=202108190815"
echo $cmd
eval $cmd
