#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/:/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/tHFCNC/NtupleProducer/:/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/tHFCNC/NtupleProducer/obj:/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_25/lib/slc6_amd64_gcc530/

echo $LD_LIBRARY_PATH

./gen list.txt output
