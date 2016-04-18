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
tree = ET.ElementTree(file='../config/Run2SameSignDiLepton_samples.xml')
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

filename="mergedHistos"+str(channel)+".sh"
os.remove(str(filename))
oufile = open(filename,'a')
# loop over all the dataset with add="1"
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        pathdir = "../OutPutHistos/"+"merged"+str(channel)+"AllSamples"
        if not os.path.exists(pathdir):
           os.mkdir(str(pathdir))
        commandfirst = "rm " + str(pathdir)+"/merged"+ str(channel) + str(d.attrib['name']) + ".root"
        commandString = "hadd -f " + str(pathdir)+"/merged"+ str(channel) + str(d.attrib['name']) + ".root ../OutPutHistos/"+str(channel)+str(d.attrib['name'])+"/*.root" 
        print >> oufile, commandfirst
        print >> oufile, commandString
