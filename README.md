# FCNCAnalysis

cmsrel CMSSW_8_0_21

cd CMSSW_8_0_21/src/

cmsenv


git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_80X #This is the developer's branch for TTP. Change to the tag according to the samples you are using

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



git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis

cd FCNCAnalysis/

git checkout CMSSW_80X_1L3B

cd TopKinFit

make clean

make


# Workflow

./CompileMacro

cd localsubmission

python createSubmitScriptWithCopy.py #creates the submission scripts.

cd SubmitScripts/DATE/CHANNEL #You can run 1 interactively to make sure there are no bugs in the code in the test directory

source SubmitAll.sh #submits the jobs. If you have a lot of jobs (>1500), you can run 'source MakeTXTforMassiveSubmission.sh; big-submission massiveJobSumbission.txt'

check status with 'qstat -u kderoove' #Or on the webpage http://mon.iihe.ac.be/jobview/overview.html

cd ../../scripts #change the channel and date in Merger.py

python Merge.py #merges the output from your Ntupler

./TreeProcessor 0 0 _CHANNEL _DATE 0 #Runs the code on the merged trees and makes plots


