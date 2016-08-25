#!/bin/env zsh

cp /tmp/x509up_u7989 /home-pbs/kskovpen/proxy/.

que="dteam"

jName=${1}

export HOME=$(pwd)

dout="/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHad/"

fpath="${HOME}/lists/"

runName="run${jName}"
logName="log${jName}"

rm -rf ${runName}
mkdir ${runName}
rm -rf ${logName}
mkdir ${logName}

flist=$(ls ${fpath})

echo $flist | while read line
do
  jidx=0
  sample=$(echo $line | sed 's%.txt%%g')
  dataset=$(echo $sample | sed 's%_ID..*%%g')
  if [[ ! -d ${runName}/${dataset} ]]; then
    mkdir ${runName}/${dataset}
  fi

  fout=$(echo ${runName}/${dataset}/${line}_${jidx} | sed 's%.txt%%g')
  lout=$(echo ${line}_${jidx} | sed 's%.txt%%g')

  echo "${dataset}"
  qsub -N gen -q ${que} -o ${logName}/${sample}.log -j oe job.sh \
-v dout=${dout},line2=${fpath}${line},fout=${fout}
done
