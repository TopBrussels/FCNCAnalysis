
#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=03:00:00

source /user/ivanpari/.bash_login
source $VO_CMS_SW_DIR/cmsset_default.sh
# setting up your code and your env
cd /user/ivanpari/CMSSW_7_4_15/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis

# want you really want to do!!
./TreeMaker NP_overlay_FCNC FCNCttbartZqto3lNu 1 5 1 2 1 44786424.172 0.0220156 0.0 /pnfs/iihe/cms/store/user/ivanpari/20160114_SignalKiril_TopTree_CMSSW7_4_15_V9/TT_tZq_FCNC/TOPTREE.root   mumumu   1  0  10000
