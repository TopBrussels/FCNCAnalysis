#!/bin/bash 

mkdir Terminal_Output


qsub -l walltime=05:00:00 -o Terminal_Output/b2j3.stdout -e Terminal_Output/b2j3.stderr submit_b2j3.sh
qsub -l walltime=22:00:00 -o Terminal_Output/b2j4.stdout -e Terminal_Output/b2j4.stderr submit_b2j4.sh
qsub -l walltime=01:30:00 -o Terminal_Output/b3j3.stdout -e Terminal_Output/b3j3.stderr submit_b3j3.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b3j4.stdout -e Terminal_Output/b3j4.stderr submit_b3j4.sh
qsub -l walltime=01:30:00 -o Terminal_Output/b4j4.stdout -e Terminal_Output/b4j4.stderr submit_b4j4.sh

