import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
import shutil 


# getting the appropriate xml file & defining channel and production date of TreeMakerTrees
channel = "_All"
date_output = "_18_3_2017"

tree = ET.ElementTree(file='../config/FullMcBkgdSamples_Mu.xml')

inputdir_Ntuples_Mu = "../NtuplerOutput/Ntuples_Mu/Ntuples_19_3_2017"
inputdir_Ntuples_Mu_2 = "../NtuplerOutput/Ntuples_Mu/Ntuples_20_3_2017"
inputdir_Ntuples_El = "../NtuplerOutput/Ntuples_El/Ntuples_19_3_2017"
inputdir_Ntuples_El_2 = "../NtuplerOutput/Ntuples_El/Ntuples_20_3_2017"

mergedpath_1 = "../Merged/Ntuples_All"
if not os.path.exists(mergedpath_1):
  os.mkdir(str(mergedpath_1))
mergedpath = "../Merged/Ntuples_All/Ntuples"+str(date_output)
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
    filenames_ntuples_Mu = glob.glob(inputdir_Ntuples_Mu + "/*" + n + "*.root")
    filenames_ntuples_El = glob.glob(inputdir_Ntuples_El + "/*" + n + "*.root")
    filenames_ntuples_Mu_2 = glob.glob(inputdir_Ntuples_Mu_2 + "/*" + n + "*.root")
    filenames_ntuples_El_2 = glob.glob(inputdir_Ntuples_El_2 + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + "JERMinus.root"
    for ff in filenames_ntuples_Mu:
       if 'JERMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El:
       if 'JERMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_Mu_2:
       if 'JERMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El_2:
       if 'JERMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

for n in datasetNames:
    filenames_ntuples_Mu = glob.glob(inputdir_Ntuples_Mu + "/*" + n + "*.root")
    filenames_ntuples_El = glob.glob(inputdir_Ntuples_El + "/*" + n + "*.root")
    filenames_ntuples_Mu_2 = glob.glob(inputdir_Ntuples_Mu_2 + "/*" + n + "*.root")
    filenames_ntuples_El_2 = glob.glob(inputdir_Ntuples_El_2 + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + "JERPlus.root"
    for ff in filenames_ntuples_Mu:
       if 'JERPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El:
       if 'JERPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_Mu_2:
       if 'JERPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El_2:
       if 'JERPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

for n in datasetNames:
    filenames_ntuples_Mu = glob.glob(inputdir_Ntuples_Mu + "/*" + n + "*.root")
    filenames_ntuples_El = glob.glob(inputdir_Ntuples_El + "/*" + n + "*.root")
    filenames_ntuples_Mu_2 = glob.glob(inputdir_Ntuples_Mu_2 + "/*" + n + "*.root")
    filenames_ntuples_El_2 = glob.glob(inputdir_Ntuples_El_2 + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + "JESMinus.root"
    for ff in filenames_ntuples_Mu:
       if 'JESMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El:
       if 'JESMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_Mu_2:
       if 'JESMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El_2:
       if 'JESMinus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

for n in datasetNames:
    filenames_ntuples_Mu = glob.glob(inputdir_Ntuples_Mu + "/*" + n + "*.root")
    filenames_ntuples_El = glob.glob(inputdir_Ntuples_El + "/*" + n + "*.root")
    filenames_ntuples_Mu_2 = glob.glob(inputdir_Ntuples_Mu_2 + "/*" + n + "*.root")
    filenames_ntuples_El_2 = glob.glob(inputdir_Ntuples_El_2 + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + "JESPlus.root"
    for ff in filenames_ntuples_Mu:
       if 'JESPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El:
       if 'JESPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_Mu_2:
       if 'JESPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El_2:
       if 'JESPlus' in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

for n in datasetNames:
    filenames_ntuples_Mu = glob.glob(inputdir_Ntuples_Mu + "/*" + n + "*.root")
    filenames_ntuples_El = glob.glob(inputdir_Ntuples_El + "/*" + n + "*.root")
    filenames_ntuples_Mu_2 = glob.glob(inputdir_Ntuples_Mu_2 + "/*" + n + "*.root")
    filenames_ntuples_El_2 = glob.glob(inputdir_Ntuples_El_2 + "/*" + n + "*.root")
    commandString_ntuples = "hadd -f " + mergedpath + "/FCNC_1L3B__Run2_TopTree_Study_" + n + ".root"
    for ff in filenames_ntuples_Mu:
       if 'JERPlus' not in ff and 'JERMinus' not in ff and 'JESPlus' not in ff and 'JESMinus' not in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El:
       if 'JERPlus' not in ff and 'JERMinus' not in ff and 'JESPlus' not in ff and 'JESMinus' not in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_Mu_2:
       if 'JERPlus' not in ff and 'JERMinus' not in ff and 'JESPlus' not in ff and 'JESMinus' not in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    for ff in filenames_ntuples_El_2:
       if 'JERPlus' not in ff and 'JERMinus' not in ff and 'JESPlus' not in ff and 'JESMinus' not in ff:
          commandString_ntuples = commandString_ntuples + " " + ff
    #print "Merging control plots for " + str(d.attrib['name'])
    #os.system(commandString_ctrplots)
    print "Merging ntuples for " + str(d.attrib['name'])
    os.system(commandString_ntuples)

