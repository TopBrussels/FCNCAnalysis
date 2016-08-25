#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/KinFit/:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/obj

./test list.txt output 5000
