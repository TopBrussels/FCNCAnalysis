# FCNCAnalysis
--------------------------------------------
--- FCNC analysis in the TopTreeFramework --
--------------------------------------------

1. Getting started
---------------------


1) Set up CMSSW

--To get correct versions of ROOT, python, ... (You need ROOT 6)--

export SCRAM_ARCH=slc6_amd64_gcc491

cmsrel CMSSW_7_4_2

cd CMSSW_7_4_2/src

cmsenv

git cms-init

2) Get TopTreeProducer from git

--Make sure to add the 'TopBrussels' directory. Otherwise the compilation later on will fail.--

git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_74X

cd src

make

cd ../../..

3) Get TopTreeAnalysisBase from git

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

git checkout CMSSW_74X

make

cd ../..

4) Compile CMSSW and TopBrussels packages

scram b -j 16

5) Get private code directory from git

--Make a new local repository on the website and check the box to create a readme file. Now you can clone this repository in different ways: If you use the https link, you will need to enter your username and password each time you push, pull or fetch. If you use the ssh link, you need to enter the passphrase of your ssh key.--

git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis

git checkout CMSSW_74X

--Make a local branch--

git checkout -b CMSSW_74X_TriLepton





2) Compile and run:
-------------------
1) Way  of Kevin
run the 'CompileMacro' to compile and launch the 'LaunchParallelMacro.py' as 'python LaunchParallelMacro.py'


2) Way of Lies
source compile.sh
Run

./testAnalyser

make

cd ../


