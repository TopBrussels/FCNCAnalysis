#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/

nToys=1
isSig=1
applyMVA=0
nNonBJetMax=-1
coup="Hut"

./test list.txt output ${nToys} ${isSig} ${applyMVA} ${nNonBJetMax} ${coup}
