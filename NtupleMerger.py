
from ROOT import TChain
import glob
import xml.etree.cElementTree as ET
import os

#Define path where ntuples are stored
#pathDateTrees = "Trees_SelectionOutput_24_11_2015/" #needs to be changed for every different set of ntuples
pathNonMerged = "Ntuples_Bkg/"  #needs to be changed for different lepton channel
pathMerged = "Merged/" #+ pathDateTrees

if not os.path.exists(pathMerged):
    os.makedirs(pathMerged)


# get filenames from the xml!!!
tree = ET.ElementTree(file='config/Run2TriLepton_samples.xml')

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
    filenames = glob.glob(pathNonMerged + "/*" + n + "*.root")
    hadd = "hadd " + pathMerged + n + "tree.root"
    for f in filenames:
       hadd = hadd + " " + f
    print "Merging ntuples for " + n
    os.system(hadd)
