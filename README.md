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



# Workflow

./CompileMacro

cd localsubmission

python createSubmitScript.py #creates the submission scripts. You can run 1 interactively to make sure there are no bugs in the code

source SubmitAll.sh

#check status with 'qstat -u kderoove'

cd ../../scripts #change the channel and date in Merger.py & MakePlots.C 

source PlotBomb.sh #This will merge your output controlPlots, as well as your ntuples into 1 file and store them in the directory ../Merged, with the appropriate tag. At the same time, the controlplots will be stacked into the directory ../Plots


