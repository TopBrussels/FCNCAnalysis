#!/bin/bash 

cd SubmitScripts/
for f in *.sh
#for f in SubmitScripts/submit_Data1.sh SubmitScripts/submit_Data2.sh
do
#    echo $f
#    qsub $f -q express
    qsub -l walltime=00:55:00 -o ../Terminal_Output/$f.stdout -e ../Terminal_Output/$f.stderr $f
#    qsub $f -q express -l walltime=00:05:00 
done
