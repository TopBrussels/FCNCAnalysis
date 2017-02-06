#!/bin/env zsh

cp /tmp/x509up_u20657 /user/kskovpen/proxy/.

que="localgrid"

jName=${1}

export HOME=$(pwd)

dout="/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/test/GenAnalysis/TopHLepbb/"

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
  while [[ $(qsub -N gen -q ${que} \
-l walltime=0:20:00 \
-o ${logName}/${sample}.log -j oe job.sh \
-v dout=${dout},line2=${fpath}${line},fout=${fout} 2>&1 | grep "Invalid credential") != "" ]];
  do echo "one more try";
  done
done
