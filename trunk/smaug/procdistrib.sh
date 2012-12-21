#!/bin/bash


start=0
end=5  
step=1
proc="np0801"
current=$start

name="zeroOT"
#dir="/data/cs1mkg/smaug/out"
indir="out"
outdir="newout"
numfiles=99
numprocd=0
 

while [ "$numprocd" -le $numfiles ]
do
        numprocd=`expr $numprocd + 1`
        current=`expr $current + $step`

        echo "processing $numprocd step $current"
        fin="${indir}/${name}_${current}_${proc}.out"
        fout="${outdir}/${name}_${current}.out"
        echo "fileout is $fout"
        echo "file in is $fin"
        ./distribution $fin $fout
        echo

done
