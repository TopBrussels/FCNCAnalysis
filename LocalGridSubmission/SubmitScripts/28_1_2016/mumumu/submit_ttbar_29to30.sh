
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
./TreeMaker ttbar ttbar 1 633 1 2 1 114447.042416 831.76 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_160929/0000/TOPTREE_149.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_160929/0000/TOPTREE_59.root   mumumu   15  0  2000000
