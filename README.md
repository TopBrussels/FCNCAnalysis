# FCNCAnalysis
--------------------------------------------
--- FCNC analysis in the TopTreeFramework --
--------------------------------------------


1) Setup instructions:
----------------------

cmsrel CMSSW_7_4_2
cd CMSSW_7_4_2/src/
cmsenv
git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer
cd TopBrussels/TopTreeProducer/
git checkout -b CMSSW_70X phys14.v1   # this is where the BRANCH (CMSSW70X) and the TAG (phys14.v1) are defined!!
scram b clean
scram b
cd src
make
cd ../../

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase
cd TopTreeAnalysisBase/
git checkout master
make
cd ../

git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis
cd FCNCAnalysis/
git checkout CMSSW_74X



2) Compile and run:
-------------------

run the 'CompileMacro' to compile and launch the 'LaunchParallelMacro.py' as 'python LaunchParallelMacro.py'


