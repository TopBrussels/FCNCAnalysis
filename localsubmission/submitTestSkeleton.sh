#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
source /user/kderoove/.bash_profile
cd /localgrid/kderoove/FCNC/TopTreeFramework/CMSSW_7_6_3/src/TopBrussels/FCNCAnalysis/
cmsenv
eval `scramv1 runtime -sh`

# want you really want to do!!


