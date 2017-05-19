#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
source /user/kderoove/.bash_profile
cd /user/kderoove/FCNC/TopTreeFramework_Run2/CMSSW_8_0_26_patch1/src/TopBrussels/FCNCAnalysis
cmsenv
eval `scramv1 runtime -sh`

./TreeProcessor_FinalMVA 3 3 hct _All _12_5_2017 0 1 1 0 0 0 0
