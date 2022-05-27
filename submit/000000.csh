#!/bin/csh -f 
cd /disk1/disk02/usr6/nakanisi/SuperNova/VectGen-oxy/submit
cd ..
source /usr/local/sklib_gcc8/skofl-trunk/env.csh
hostname
./main_snburst nakazato/intp2002.data 0 1 0 1000
