#!/bin/env zsh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/

nToys=20
isSig=1
applyMVA=1
nNonBJetMax=-1

./test list.txt output ${nToys} ${isSig} ${applyMVA} ${nNonBJetMax}
