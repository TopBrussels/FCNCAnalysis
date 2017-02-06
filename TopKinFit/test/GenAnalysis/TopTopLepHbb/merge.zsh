#!/bin/env zsh

dirName=${1}

#if [[ ${dirName} == "" ]]; then
#  echo "Please specify the run name"
#  exit
#fi

dpath="/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/test/GenAnalysis/TopTopLepHbb/run${dirName}/"
opath="/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/test/GenAnalysis/TopTopLepHbb/run${dirName}_MERGED/"

rm -rf ${opath}
mkdir ${opath}

flist=$(ls ${dpath})

echo $flist | while read dname
do
loci=${dpath}${dname}
loco=${opath}${dname}
mkdir ${loco}
rm -f ${loco}/data.root
flist=$(ls ${loci})
farr=()
idx=1
echo $flist | while read fname
do
farr[$idx]=${loci}/$fname
idx=$[$idx+1]
done
hadd -f ${loco}/data.root ${farr[*]}
done
