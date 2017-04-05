#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
source /user/kderoove/.bash_profile
cd /user/kderoove/FCNC/TopTreeFramework_Run2/CMSSW_8_0_24/src/TopBrussels/FCNCAnalysis
cmsenv
eval `scramv1 runtime -sh`

./TreeProcessor_FinalMVA 4 4 hut _All _4_2_2017 0 1 1 0