
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
./TreeMaker TTJets t\bar{t}+jets_Madgraph_MLM 1 633 1 2 1 13632.8171588 831.76 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151020_160750/0000/TOPTREE_32.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151020_160750/0000/TOPTREE_30.root   mumumu   18  0  2000000
