#!/bin/bash

#PBS -q express


# setting up your code and your env
source /user/sabuzeid/.bash_login
cd /localgrid/sabuzeid/FCNC_Study_LocalGrid/CMSSW_7_4_15/src/
#cmsenv
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis/
# want you really want to do!!
./Ntupler_WorkesOnLocalGrid DYJets DY+jets 1 kBlue-7 1 2 1 1502.46813384 6025.2 0.0 dcap://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_154341/0000/TOPTREE_17.root   28  0  2000000
