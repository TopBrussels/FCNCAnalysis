#!/bin/bash 

mkdir Terminal_Output


for f in submit*.sh
do
#    echo $f
#    qsub $f -q express
    qsub -o -l walltime=14:00:00 Terminal_Output/$f.stdout -e Terminal_Output/$f.stderr $f
#    qsub $f -q express -l walltime=00:05:00 
done
