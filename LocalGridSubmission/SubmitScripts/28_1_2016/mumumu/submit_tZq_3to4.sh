
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
./TreeMaker tZq tZqto3lNu 1 5 1 2 1 38713601.5831 0.0758 0.0  dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151021_135224/0000/TOPTREE_14.root dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151021_135224/0000/TOPTREE_16.root   mumumu   2  0  2000000
