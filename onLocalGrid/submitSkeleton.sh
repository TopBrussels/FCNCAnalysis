
#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
cd /user/ivanpari/CMSSW_7_4_15/src/ 
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis

# want you really want to do!!
