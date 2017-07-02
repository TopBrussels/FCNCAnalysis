#!/bin/bash 

mkdir Terminal_Output


qsub -l walltime=05:00:00 -o Terminal_Output/b2j3hct.stdout -e Terminal_Output/b2j3hct.stderr submit_hct_b2j3.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b2j4hct.stdout -e Terminal_Output/b2j4hct.stderr submit_hct_b2j4.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b3j3hct.stdout -e Terminal_Output/b3j3hct.stderr submit_hct_b3j3.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b3j4hct.stdout -e Terminal_Output/b3j4hct.stderr submit_hct_b3j4.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b4j4hct.stdout -e Terminal_Output/b4j4hct.stderr submit_hct_b4j4.sh

qsub -l walltime=05:00:00 -o Terminal_Output/b2j3hut.stdout -e Terminal_Output/b2j3hut.stderr submit_hut_b2j3.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b2j4hut.stdout -e Terminal_Output/b2j4hut.stderr submit_hut_b2j4.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b3j3hut.stdout -e Terminal_Output/b3j3hut.stderr submit_hut_b3j3.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b3j4hut.stdout -e Terminal_Output/b3j4hut.stderr submit_hut_b3j4.sh
qsub -l walltime=05:00:00 -o Terminal_Output/b4j4hut.stdout -e Terminal_Output/b4j4hut.stderr submit_hut_b4j4.sh

