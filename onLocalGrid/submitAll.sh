#!/bin/bash 

for f in SubmitScripts/*.sh
do
    qsub -q localgrid@cream02  $f 
done
