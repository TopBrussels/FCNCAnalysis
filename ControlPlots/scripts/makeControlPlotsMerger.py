import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
# libray to copy files
import shutil 


#mkdir ControlPlots/_ElEl_tZq/Merged/

#rm ControlPlots/_ElEl_tZq/Merged/merged_ElEl_tZq.root

#hadd -f ControlPlots/_ElEl_tZq/Merged/merged_ElEl_tZq.root ControlPlots/_ElEl_tZq/*.root


# getting the appropriate xml file
tree = ET.ElementTree(file='../../config/Run2TriLepton_samples.xml')
#tree = ET.ElementTree(file='../config/DataSamples.xml')


root = tree.getroot()
datasets = root.find('datasets')

print "found  "  + str(len(datasets)) + " datasets"

channel = "_ElEl_"
# create new dir if not already existing
#if not os.path.exists(""):
#    os.makedirs("SubmitScripts")


# vector containing all the root file for a given dataset
topTrees = []

filename="mergeControlplots"+str(channel)+".sh"
os.remove(str(filename))
oufile = open(filename,'a')
# loop over all the dataset with add="1"
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        pathdir = "../"+str(channel)+ "allSamples/"
        if not os.path.exists(pathdir):
           os.mkdir(str(pathdir))
        commandString = "hadd -f " + str(pathdir)+"/merged"+ str(channel) + str(d.attrib['name']) + ".root ../"+str(channel)+str(d.attrib['name'])+"/*.root" 
        print >> oufile, commandString
