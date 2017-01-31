import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
import shutil 


# getting the appropriate xml file & defining channel and production date of TreeMakerTrees
channel = "_All"
date = "_30_1_2017"

tree = ET.ElementTree(file='../config/FullMcBkgdSamples_Mu.xml')

inputdir_Ntuples_Mu = "../Merged/Ntuples_Mu/Ntuples_30_1_2017"
#inputdir_Ntuples_Mu_2 = "../Merged/Ntuples_Mu/Ntuples_28_1_2017"
inputdir_Ntuples_El = "../Merged/Ntuples_El/Ntuples_30_1_2017"
#inputdir_Ntuples_El_2 = "../Merged/Ntuples_El/Ntuples_28_1_2017"

mergedpath_1 = "../Merged/Ntuples_All"
if not os.path.exists(mergedpath_1):
  os.mkdir(str(mergedpath_1))
mergedpath = "../Merged/Ntuples_All/Ntuples"+str(date)
if not os.path.exists(mergedpath):
  os.mkdir(str(mergedpath))


root = tree.getroot()
datasets = root.find('datasets')
print "found  "  + str(len(datasets)) + " datasets"


datasetNames = []

# loop over the datasets to be added and fill the "topTrees" vector
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        datasetNames.append(str(d.attrib['name']))

for n in datasetNames:
    filenames_ntuples_Mu = inputdir_Ntuples_Mu + "/*" + n + ".root"
    filenames_ntuples_El = inputdir_Ntuples_El + "/*" + n + ".root"
#    filenames_ntuples_Mu_2 = inputdir_Ntuples_Mu_2 + "/*" + n + ".root"
#    filenames_ntuples_El_2 = inputdir_Ntuples_El_2 + "/*" + n + ".root"
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + ".root"
    commandString_ntuples = commandString_ntuples + " " + filenames_ntuples_Mu + " " + filenames_ntuples_El# + " " + filenames_ntuples_Mu_2 + " " + filenames_ntuples_El_2
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

