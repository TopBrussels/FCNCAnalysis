#!/bin/bash 

mkdir Terminal_Output


for f in submit*.sh
do
#    echo $f
#    qsub $f -q express
    echo -n "qsub -l walltime=24:00:00 -o Terminal_Output/"$f".stdout -e Terminal_Output/"$f".stderr "$f >> massiveJobSumbission.txt
    echo " " >> massiveJobSumbission.txt
#    qsub $f -q express -l walltime=00:10:00 
done
