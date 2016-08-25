#!/bin/env zsh

cp /tmp/x509up_u7989 /home-pbs/kskovpen/proxy/.

#que="sbg_local"
#que="cms_local_short"
que="dteam"
#que="cms"
#que="cms_local"

jName=${1}

export HOME=$(pwd)

#nToys=100
#nToys=1000
nToys=5000
#nToys=100000
#nToys=30000

dout="/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/KinFit/test/Validation/TopTopLepLep/"

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
  qsub -N KinFit -q ${que} -o ${logName}/${sample}.log -j oe job.sh \
-v dout=${dout},line2=${fpath}${line},fout=${fout},nToys=${nToys}
done
