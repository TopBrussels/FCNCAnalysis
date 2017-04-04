# FCNCAnalysis
# Setup Enviroment for Analysis 
cmsrel CMSSW_8_0_26_patch2
cmsenv
# To clone TTP & TTAB
git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer
cd TopBrussels/TopTreeProducer/
git checkout CMSSW_80X_v7 #latest Version which work with reMINIAOD
scram b clean
scram b
cd src
make
cd ../../

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase
cd TopTreeAnalysisBase/
git checkout CMSSW_80X
make
cd ../
#setup FCNC 
git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis
cd FCNCAnalysis/
git checkout CMSSW_80X_2SSL
