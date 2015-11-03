# FCNCAnalysis

cmsrel CMSSW_7_6_0

cd CMSSW_7_6_0/src/

cmsenv


git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_76X #This is the developer's branch for TTP. Change to the tag according to the samples you are using

scram b clean

scram b

cd src

make

cd ../../



git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase

cd TopTreeAnalysisBase/

git checkout CMSSW_76X

make

cd ../



git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis

cd FCNCAnalysis/

git checkout CMSSW_76X_1L3B


