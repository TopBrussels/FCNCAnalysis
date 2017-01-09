import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
# libray to copy files
import shutil 


# getting the appropriate xml file & defining channel and production date of TreeMakerTrees
channel = "_Mu"
date = "_6_1_2017"


if channel == "_Mu":
    tree = ET.ElementTree(file='../config/FullMcBkgdSamples_Mu.xml')
elif channel == "_El":
    tree = ET.ElementTree(file='../config/FullMcBkgdSamples_El.xml')


#Making and setting all relevant paths
#inputdir_ControlPlots = "../NtuplerOutput/MACRO_histos"+str(channel)+"/MACRO_histos"+str(date)
inputdir_Ntuples = "../NtuplerOutput/Ntuples"+str(channel)+"/Ntuples"+str(date)
mergedpath = "../Merged"
if not os.path.exists(mergedpath):
  os.mkdir(str(mergedpath))
#pathdir_ctrplots_1 = mergedpath+"/ControlPlots"+ str(channel)
#pathdir_ctrplots = pathdir_ctrplots_1+"/ControlPlots" + str(date)
#if not os.path.exists(pathdir_ctrplots_1):
#  os.mkdir(str(pathdir_ctrplots_1))
#if not os.path.exists(pathdir_ctrplots):
#  os.mkdir(str(pathdir_ctrplots))
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
    filenames_ntuples = glob.glob(inputdir_Ntuples + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + pathdir_ntuples + "/FCNC_1L3B__Run2_TopTree_Study_" + n + ".root"
    for ff in filenames_ntuples:
       commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

