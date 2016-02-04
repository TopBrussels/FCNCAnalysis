#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
source /user/sabuzeid/.bash_login
cd /localgrid/sabuzeid/FCNC_Study_LocalGrid/CMSSW_7_4_15/src/
#cmsenv
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis/
# want you really want to do!!
