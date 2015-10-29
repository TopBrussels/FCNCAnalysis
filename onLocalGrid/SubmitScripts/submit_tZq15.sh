
#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=02:00:00

# setting up your code and your env
cd /user/ivanpari/CMSSW_7_4_15/src/ 
eval `scramv1 runtime -sh`
cd TopBrussels/FCNCAnalysis

# want you really want to do!!
./Ntupler_localGrid tZq tZqto3lNu 1 5 1 2 1 30125151.4218 0.09741 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151021_135224/0000/TOPTREE_9.root  0  2000000   15
