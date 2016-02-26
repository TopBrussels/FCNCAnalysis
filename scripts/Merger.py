import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
# libray to copy files
import shutil 


# getting the appropriate xml file
tree = ET.ElementTree(file='../config/FullMcBkgdSamplesV8_Mu.xml')
#defining channel and production date of TreeMakerTrees
channel = "_Mu"
date = "_17_2_2016"


#Making and setting all relevant paths
inputdir_ControlPlots = "../TreeMakerOutput/MACRO_histos"+str(channel)+"/MACRO_histos"+str(date)
inputdir_Ntuples = "../TreeMakerOutput/Ntuples"+str(channel)+"/Ntuples"+str(date)
mergedpath = "../Merged"
if not os.path.exists(mergedpath):
  os.mkdir(str(mergedpath))
pathdir_ctrplots_1 = mergedpath+"/ControlPlots"+ str(channel)
pathdir_ctrplots = pathdir_ctrplots_1+"/ControlPlots" + str(date)
if not os.path.exists(pathdir_ctrplots_1):
  os.mkdir(str(pathdir_ctrplots_1))
if not os.path.exists(pathdir_ctrplots):
  os.mkdir(str(pathdir_ctrplots))
pathdir_ntuples_1 = mergedpath+"/Ntuples"+ str(channel)
pathdir_ntuples = pathdir_ntuples_1+"/Ntuples" + str(date)
if not os.path.exists(pathdir_ntuples_1):
  os.mkdir(str(pathdir_ntuples_1))
if not os.path.exists(pathdir_ntuples):
  os.mkdir(str(pathdir_ntuples))

#reading xml file
root = tree.getroot()
datasets = root.find('datasets')
print "found  "  + str(len(datasets)) + " datasets"


# vector containing all the root file for a given dataset
topTrees = []
datasetNames = []

# loop over the datasets to be added and fill the "topTrees" vector
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        datasetNames.append(str(d.attrib['name']))

for n in datasetNames:
    filenames_ctrlplots = glob.glob(inputdir_ControlPlots + "/*" + n + "*.root")
    commandString_ctrplots = "hadd -f " + pathdir_ctrplots + "/FCNC_1L3B__Run2_TopTree_Study_" + n + ".root"
    for f in filenames_ctrlplots:
       commandString_ctrplots = commandString_ctrplots + " " + f
    filenames_ntuples = glob.glob(inputdir_Ntuples + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + pathdir_ntuples + "/FCNC_1L3B__Run2_TopTree_Study_" + n + ".root"
    for ff in filenames_ntuples:
       commandString_ntuples = commandString_ntuples + " " + ff
    print "Merging control plots for " + str(d.attrib['name'])
    os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)
