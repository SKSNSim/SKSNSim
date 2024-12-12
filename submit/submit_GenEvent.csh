#!/bin/csh -f

if ($1 == "") then
    echo "Usage : submit_GenEvent.csh"
    echo " 1st: model name (e.g. nakazato/intp2001.data)"
    echo " 2nd: neutrino-oscillation (0:no, 1:normal, 2:inverted)"
    echo " 3rd: Distance (normalized to 10kpc)"
    echo " 4th: Generate event or not (0 : No just calculate expected, 1: Yes for detector simulation)"
    echo " 5th: Output directory, default: ./data/(SN model)"
    echo " 6th: Random seed (not supported yet)"
    exit
endif

set OUTDIR = "/disk2/disk03/lowe11/nakanisi/SKSNSim_241122/"
set LOGDIR = "/disk2/disk03/lowe11/nakanisi/SKSNSim_241122/"
#set num_simulation = 1000
set QUEUE = lowe
set JOBLIMIT = 250
set num_simulation = 10000


if(! -d /disk2/disk03/lowe11/nakanisi/SKSNSim_241122) then
    mkdir -p /disk2/disk03/lowe11/nakanisi/SKSNSim_241122;
    echo $OUTDIR;
endif

if(! -d $LOGDIR/script) then
    mkdir -p $LOGDIR/script;
    echo $LOGDIR/script;
endif

if(! -d $LOGDIR/script/$1) then
    mkdir -p $LOGDIR/script/$1;
    echo $LOGDIR/script/$1;
endif

if(! -d $OUTDIR/data/$1) then
    mkdir -p $OUTDIR/data/$1;
    echo $OUTDIR/data/$1;
endif

if(! -d $LOGDIR/out/$1) then
    mkdir -p $LOGDIR/out/$1;
    echo $LOGDIR/out/$1;
endif

if(! -d $LOGDIR/err/$1) then
    mkdir -p $LOGDIR/err/$1;
    echo $LOGDIR/err/$1;
endif

@ cur_num = 0
while ($cur_num < $num_simulation)

    if($cur_num < 10) then
       set fnum = "00000"$cur_num
    else if($cur_num < 100) then
       set fnum = "0000"$cur_num
    else if($cur_num < 1000) then
       set fnum = "000"$cur_num
    else if($cur_num < 10000) then
       set fnum = "00"$cur_num
    else if($cur_num < 100000) then
       set fnum = "0"$cur_num
    else if($cur_num < 1000000) then
       set fnum = $cur_num
    else
       echo "too much.."
       exit
    endif

    #set odir = $OUTDIR"/data/"$1"/784kpc_NO_"$fnum
    set odir = $OUTDIR"/data/"$1"/10kpc_WO_"$fnum
    if(! -d $odir/event) then
	mkdir -p $odir/event
    endif

    set ofile = $LOGDIR"/script/"$1"/"$fnum".csh"
    @ random = $cur_num + 1000

    echo $ofile
    echo '#\!/bin/csh -f ' >! $ofile
    #echo 'cd '$PWD >> $ofile
    #echo 'cd /disk2/disk03/usr8/nakanisi/SKSNSim/submit/' >> $ofile
    #echo 'cd /disk2/disk03/usr8/nakanisi/test_SKSNSim/submit/' >> $ofile
    echo 'cd /disk2/disk03/usr8/nakanisi/SKSNSim_241122/submit/' >> $ofile
    #echo 'cd /disk2/disk03/usr8/nakanisi/SN2023ixf/SKSNSim/submit/' >> $ofile
    echo 'cd ..' >> $ofile
    echo 'source /usr/local/sklib_gcc8/skofl-trunk/env.csh' >> $ofile
    echo 'source SKSNSimsource.csh' >> $ofile
    echo 'hostname' >> $ofile
    echo 'set SKSNSIMDATADIR = /disk2/disk03/usr8/nakanisi/SKSNSim_241122/data' >> $ofile
    echo "./bin/main_snburst "$1" "$2" "$3" "$4" "$odir" "$random "--time_min 0. --time_max 300. --time_nbins 300000 --energy_min 0.5 --energy_max 100. --energy_nbins 995">> $ofile
    #echo "./bin/main_snburst "$1" "$2" "$3" "$4" "$odir" "$random "--time_min 0. --time_max 20. --time_nbins 20000 --energy_min 0.5 --energy_max 100. --energy_nbins 995">> $ofile
    #echo "./bin/main_snburst_new "$1" "$2" "$3" "$4" "$odir" "$random "--time_min 0. --time_max 150. --time_nbins 150000 --energy_min 0. --energy_max 100. --energy_nbins 1000">> $ofile
    #echo "./bin/main_snburst "$1" "$2" "$3" "$4" "$odir" "$random >> $ofile

    chmod 755 $ofile

    set njobs = `pjstat -v | grep $QUEUE | wc -l` # cont number if running jobs
    echo "current jobs: $njobs"
    while ( $njobs > $JOBLIMIT ) # monitor number if running jobs
        sleep 10
        set njobs = `pjstat -v | grep $QUEUE | wc -l` # count number of running jobs
        echo "current jpbs: $njobs"
    end
    pjsub -z -o $LOGDIR/out/$1/$fnum.WO.Eth0.5 -e $LOGDIR/err/$1/$fnum.WO.Eth0.5 -L rscgrp=$QUEUE $ofile
    echo "pjsub -z -o $LOGDIR/out/$1/$fnum.WO.Eth0.5 -e $LOGDIR/err/$1/$fnum.WO.Eth0.5 -L rscgrp=$QUEUE $ofile"
    
    @ cur_num = $cur_num + 1
    
end

