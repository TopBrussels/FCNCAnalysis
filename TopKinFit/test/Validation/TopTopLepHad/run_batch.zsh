#!/bin/env zsh

cp /tmp/x509up_u20657 /user/kskovpen/proxy/.

#que="sbg_local"
#que="cms_local_short"
que="localgrid"
#que="cms_local_mdm"
#que="cms_local"

jName=${1}

export HOME=$(pwd)

#nToys=100
nToys=20
#nToys=5000
#nToys=100000
#nToys=30000

isSig=1
applyMVA=0
nNonBJetMax=-1

dout="/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/test/Validation/TopTopLepHad/"

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
  while [[ $(qsub -N TopKinFit -q ${que} \
  -l walltime=6:00:00 \
-o ${logName}/${sample}.log -j oe job.sh \
-v dout=${dout},line2=${fpath}${line},fout=${fout},nToys=${nToys},isSig=${isSig},applyMVA=${applyMVA},nNonBJetMax=${nNonBJetMax} 2>&1 | grep "Invalid credential") != "" ]];
  do echo "one more try";
  done
done
