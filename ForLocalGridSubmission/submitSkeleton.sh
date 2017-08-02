#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=24:00:00

# setting up your code and your env
source /user/sabuzeid/.bash_login
source $VO_CMS_SW_DIR/cmsset_default.sh
cd /user/sabuzeid/FCNC_Study/CMSSW_8_0_26_patch1/src/
#cmsenv
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis/
# want you really want to do!!
