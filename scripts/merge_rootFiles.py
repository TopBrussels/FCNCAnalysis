from ROOT import TChain
import glob
import xml.etree.cElementTree as ET
import os

#Define path where ntuples are stored
channel = "MuMu" # to be adapted according to selected channel
runDate = "Test_OldJetCleaning_14Apr" # to be adapted to the date of analysis run
pathToDatedTrees = ""+runDate+"/"
pathToUnmergedrootFiles = "../Output_Ntuples/"
pathToUnmergedrootFiles += "Dilepton_"+channel+"_"+runDate+"/";
pathToMergedrootFiles = "../Merged_Ntuples/"
pathToMergedrootFiles += channel+"/"+pathToDatedTrees
#pathToMergedrootFiles = pathToUnmergedrootFiles+"merged_"+channel+"/"+pathToDatedTrees

if not os.path.exists(pathToMergedrootFiles):
    os.makedirs(pathToMergedrootFiles)



# get filenames from the xml!!!
tree = ET.ElementTree(file='../config/Run2SameSignDiLepton_76X_MuMu_Samples.xml')

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
    filenames = glob.glob(pathToUnmergedrootFiles + "/*" + n + "*.root")
    hadd = "hadd " + pathToMergedrootFiles + "" + n + ".root"
    for f in filenames:
        hadd = hadd + " " + f
    print "Merging ntuples for " + n
    os.system(hadd)