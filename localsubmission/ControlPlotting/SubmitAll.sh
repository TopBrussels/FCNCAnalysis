#!/bin/bash 

mkdir Terminal_Output


qsub -l walltime=10:00:00 -o Terminal_Output/b2j3.stdout -e Terminal_Output/b2j3.stderr submit_b2j3.sh
qsub -l walltime=27:00:00 -o Terminal_Output/b2j4.stdout -e Terminal_Output/b2j4.stderr submit_b2j4.sh
qsub -l walltime=10:30:00 -o Terminal_Output/b3j3.stdout -e Terminal_Output/b3j3.stderr submit_b3j3.sh
qsub -l walltime=10:00:00 -o Terminal_Output/b3j4.stdout -e Terminal_Output/b3j4.stderr submit_b3j4.sh
qsub -l walltime=10:30:00 -o Terminal_Output/b4j4.stdout -e Terminal_Output/b4j4.stderr submit_b4j4.sh
qsub -l walltime=15:30:00 -o Terminal_Output/Inclusive.stdout -e Terminal_Output/Inclusive.stderr Inclusive.sh

