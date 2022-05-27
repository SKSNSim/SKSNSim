#!/bin/csh -f 
cd /disk1/disk02/usr6/nakanisi/SKSNSim/submit
cd ..
source /usr/local/sklib_gcc8/skofl-trunk/env.csh
hostname
./main_snburst nakazato/intp2002.data 0 1 1 ./data/nakazato/intp2002.data/000000 1000
