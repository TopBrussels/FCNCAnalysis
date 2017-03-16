
#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=16:00:00

source /user/ivanpari/.bash_login
source $VO_CMS_SW_DIR/cmsset_default.sh
# setting up your code and your env
cd /user/ivanpari/CMSSW_8_0_24/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis

# want you really want to do!!
