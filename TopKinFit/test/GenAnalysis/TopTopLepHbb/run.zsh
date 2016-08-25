#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/:/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/NtupleProducer/:/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/NtupleProducer/obj

./gen list.txt output
