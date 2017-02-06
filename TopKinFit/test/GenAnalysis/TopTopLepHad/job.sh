#!/bin/sh

export X509_USER_PROXY=/user/kskovpen/proxy/x509up_u20657

WDIR=$(pwd)

line2=${line2}
fout=${fout}
sample=${sample}
dout=${dout}
proc=${proc}

export ROOTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed4
ls $ROOTSYS/bin/thisroot.sh
source $ROOTSYS/bin/thisroot.sh
rootV=$(root-config --version)
echo "ROOT v${rootV} has been set up"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}:${dout}/../:\
/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/tHFCNC/NtupleProducer/:\
/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/tHFCNC/NtupleProducer/obj:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_25/lib/slc6_amd64_gcc530/

export LD_LIBRARY_PATH=\
/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib64:\
/usr/lib64:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/lib:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_25/external/slc6_amd64_gcc530/lib/:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_25/lib/slc6_amd64_gcc530/:\
$LD_LIBRARY_PATH

export PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gcc/5.3.0/bin:$PATH

echo "Executing ./test ${line2} ${dout}${fout}"
${dout}/./gen ${line2} ${dout}${fout}
