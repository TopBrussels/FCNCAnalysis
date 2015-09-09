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

_In the next command the BRANCH (CMSSW74X) and the TAG (e.g. CMSSW_74X_v4 or CMSSW_74X_v5) are defined! Adapt this according to which TopTrees you want to analyze._

    git checkout -b CMSSW_74X CMSSW_74X_v4

_**Rule of thumb: for TopTree analysis one should indeed check out a TopTreeProducer tag (=fixed code), 
while for TopTreeProducer development one should check out only a (non-tag) branch (=code that can be updated -pulled and pushed- via git),
so for TopTreeProducer development the command should be just 'git checkout CMSSW_74X' without specifying a tag**_

    scram b clean
    scram b
    cd src
    make
    cd ../../

    git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase
    cd TopTreeAnalysisBase/
    git checkout CMSSW_74X
    make
    cd ../

    git clone https://github.com/TopBrussels/FCNCAnalysis FCNCAnalysis
    cd FCNCAnalysis/
    git checkout CMSSW_74X



2) Compile and run:
-------------------

_run the 'CompileMacro' to compile and launch the 'LaunchParallelMacro.py' as 'python LaunchParallelMacro.py'_


