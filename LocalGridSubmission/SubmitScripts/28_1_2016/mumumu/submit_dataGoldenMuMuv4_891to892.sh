
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
./TreeMaker dataGoldenMuMuv4 dataMuMu 1 1 1 2 1 1590.780213201 1 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8-GOLD/DoubleMuon/crab_DoubleMuon-Run2015D-PromptReco-v4-CMSSW_74X_v8-GOLD-74X_dataRun2_Prompt_v2/151126_112236/0000/TOPTREE_337.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8-GOLD/DoubleMuon/crab_DoubleMuon-Run2015D-PromptReco-v4-CMSSW_74X_v8-GOLD-74X_dataRun2_Prompt_v2/151126_112236/0000/TOPTREE_318.root   mumumu   446  0  2000000
