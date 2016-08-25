#!/bin/env zsh

#fpathMC="/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHbb/run_MERGED/"
fpathMC="/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopHLepbb/run/"

liMC=($(ls ${fpathMC}))

nFilesMC=1

outDir="lists/"

rm -rf ${outDir}
mkdir ${outDir}

rm -f /tmp/tempMC.txt
for line in $liMC
do
d1=$(echo $line)
liMC2=$(ls ${fpathMC}${d1})
echo $liMC2 | while read line2
do
f1=$(echo $line2)
file=$(echo ${fpathMC}${d1}/${f1})
echo "${file}" >> /tmp/tempMC.txt
done
split -a 5 -l ${nFilesMC} -d /tmp/tempMC.txt /tmp/${d1}_
lsfi=$(ls /tmp/${d1}_*)
jid=0
echo $lsfi | while read fil
do
mv ${fil} ${outDir}/${d1}_ID${jid}.txt
jid=$[$jid+1]
done
rm -f /tmp/tempMC.txt
done
