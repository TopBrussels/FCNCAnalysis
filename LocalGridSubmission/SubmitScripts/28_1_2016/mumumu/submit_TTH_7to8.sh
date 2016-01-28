
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
./TreeMaker TTH t\bar{t}H 1 833 1 2 1 14793818.0002 0.2658816 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v9/ttHTobb_M125_13TeV_powheg_pythia8/crab_ttHTobbM12513TeVpowhegpythia8RunIISpring15DR74Asympt25nsMCRUN274V9v1CMSSW74Xv9MCRUN274V9/151215_171851/0000/TOPTREE_16.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v9/ttHTobb_M125_13TeV_powheg_pythia8/crab_ttHTobbM12513TeVpowhegpythia8RunIISpring15DR74Asympt25nsMCRUN274V9v1CMSSW74Xv9MCRUN274V9/151215_171851/0000/TOPTREE_13.root   mumumu   4  0  2000000
