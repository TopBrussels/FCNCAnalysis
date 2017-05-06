#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
source /user/kderoove/.bash_profile
cd /user/kderoove/FCNC/TopTreeFramework_Run2/CMSSW_8_0_26_patch1/src/TopBrussels/FCNCAnalysis
#cmsenv
eval `scramv1 runtime -sh`

# want you really want to do!!


#first make all the copies
./TreeProcessor_EventTraining_STandTTCombined  3 4 2D _All _18_3_2017 0 0  10  30 
