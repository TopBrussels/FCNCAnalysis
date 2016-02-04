#!/bin/bash 

for f in SubmitScripts_DiMuon/*.sh
do

    #qsub -q localgrid@cream02  $f # used when our jobs need more than 1 hr to run over 
    #qsub -q short -l walltime=1:10:00  -o ../Terminal_Output/$f.stdout -e ../Terminal_Output/$f.stderr $f  # used when the jobs take less than one hr to run
    #qsub -l walltime=40:00:00 -o Terminal_Output/$f.stdout -e Terminal_Output/$f.stderr $f # used when the jobs take less than the walltime to run
    #qsub -q localgrid@cream02 -o Terminal_Output/$f.stdout -e Terminal_Output/$f.stderr $f # used when our jobs need more than 1 hr to run over
    qsub -l walltime=20:00:00 $f
done
