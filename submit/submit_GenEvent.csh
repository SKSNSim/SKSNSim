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

#set num_simulation = 1000
set num_simulation = 1

if(! -d ./script/$1) then
    mkdir -p ./script/$1;
endif

if(! -d ./data) then
    mkdir -p ./data;
endif

if(! -d ./data/$1) then
    mkdir -p ./data/$1;
endif

if(! -d ./out/$1) then
    mkdir -p ./out/$1
endif

if(! -d ./err/$1) then
    mkdir -p ./err/$1
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

    set odir = "./data/"$1"/"$fnum
    if(! -d ../$odir/event) then
	mkdir -p ../$odir/event
    endif

    set ofile = "script/"$1"/"$fnum".csh"
    @ random = $cur_num + 1000

    echo $ofile
    echo '#\!/bin/csh -f ' >! $ofile
    echo 'cd '$PWD >> $ofile
    echo 'cd ..' >> $ofile
    echo 'source /usr/local/sklib_gcc8/skofl-trunk/env.csh' >> $ofile
    echo 'hostname' >> $ofile
    echo "./main "$1" "$2" "$3" "$4" "$odir" "$random >> $ofile

    chmod 755 $ofile

    pjsub -o out/$1/$fnum -e err/$1/$fnum $ofile
    
    @ cur_num = $cur_num + 1
    
end

