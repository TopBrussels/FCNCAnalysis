#!/bin/sh

export X509_USER_PROXY=/home-pbs/kskovpen/proxy/x509up_u7989

export LD_PRELOAD=/usr/lib64/libglobus_gssapi_gsi.so.4

WDIR=$(pwd)

line2=${line2}
fout=${fout}
sample=${sample}
dout=${dout}
nToys=${nToys}

export ROOTSYS=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/lcg/root/6.02.12-kpegke4
ls $ROOTSYS/bin/thisroot.sh
source $ROOTSYS/bin/thisroot.sh
rootV=$(root-config --version)
echo "ROOT v${rootV} has been set up"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${dout}:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/KinFit/:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/:/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/IPHCNtuple/NtupleProducer/obj

export LD_LIBRARY_PATH=\
/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib64:\
/usr/lib64:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/external/slc6_amd64_gcc493/lib/:\
/cvmfs/cms.cern.ch/slc6_amd64_gcc493/cms/cmssw/CMSSW_7_6_3/lib/slc6_amd64_gcc493/:\
$LD_LIBRARY_PATH

export PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/bin:$PATH

echo "Executing ./test ${line2} ${dout}${fout} ${nToys}"
#${dout}/./test ${line2} ${dout}${fout} ${nToys} > /dev/null
${dout}/./test ${line2} ${dout}${fout} ${nToys}
