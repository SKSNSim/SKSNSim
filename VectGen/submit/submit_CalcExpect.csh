#!/bin/csh -f

if ($1 == "") then
    echo "Usage : submit.csh filename"
    exit
endif


if(! -d ./script) then
    mkdir ./script;
endif

if(! -d ./out) then
    mkdir ./out;
endif

if(! -d ./err) then 
    mkdir ./err;
endif

set ofile = "script/"$1".csh"

echo $ofile
echo '#\!/bin/csh -f ' >! $ofile
echo 'cd '$PWD >> $ofile
echo 'cd ..' >> $ofile
echo 'source /usr/local/sklib_gcc8/skofl-trunk/env.csh' >> $ofile
echo 'hostname' >> $ofile
echo "./bin/CalcExpect "$1" 1" >> $ofile
chmod 755 $ofile

#pjsub -o out/$1 -e err/$1 $ofile
pjsub -L rscgrp=lowe -o out/$1 -e err/$1 $ofile
