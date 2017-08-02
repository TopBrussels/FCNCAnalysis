import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
# libray to copy files
import shutil 



# getting the appropriate xml file
#tree = ET.ElementTree(file='../config/Run2SameSignDiLepton_80X_ElEl_V10_Samples.xml')
#tree = ET.ElementTree(file='../config/Run2SameSignDiLepton_80X_MuMu_V10_Samples.xml')
tree = ET.ElementTree(file='../config/Run2SameSignDiLepton_80X_ElMu_V10_Samples.xml')


root = tree.getroot()
datasets = root.find('datasets')

print "found  "  + str(len(datasets)) + " datasets"

# create new dir if not already existing
#if not os.path.exists("SubmitScripts_DiMu"):
    #os.makedirs("SubmitScripts_DiMu")
#if not os.path.exists("SubmitScripts_DiElec"):
    #os.makedirs("SubmitScripts_DiElec")
if not os.path.exists("SubmitScripts_DiEMu"):
    os.makedirs("SubmitScripts_DiEMu")


# vector containing all the root file for a given dataset
topTrees = []


# loop over all the dataset with add="1"
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        commandString = "./Selection_80X_Ntupler "+str(d.attrib['name'])+" "+str(d.attrib['title'])+" "+str(d.attrib['add'])+" "+str(d.attrib['color'])+" "+str(d.attrib['ls'])+" "+str(d.attrib['lw'])+" "+str(d.attrib['normf'])+" "+str(d.attrib['EqLumi'])+" "+str(d.attrib['xsection'])+" "+str(d.attrib['PreselEff'])
#        print "commandString: %s" % commandString
        topTrees = glob.glob(d.attrib['filenames'])
	# setting the number of file per job depending whether it is data sample or not
            # this ca be tweaked
#        if "Data" in str(d.attrib['name']):
#                FilePerJob=20
#         else:
#                FilePerJob=2
		
#        print "filenames %s" % str(d.attrib['filenames'])  
#        print "glob.glob(): %s" % glob.glob(str(d.attrib['filenames']))
#        print "nb of files to be made: %d"  % len(topTrees)
        N_file = 1
        # loop over all the root files and make one job per root file
        for f in range(0,len(topTrees)):
            #filename="SubmitScripts_DiMu/submit_"+str(d.attrib['name'])+str(N_file)+".sh"
            #filename="SubmitScripts_DiElec/submit_"+str(d.attrib['name'])+str(N_file)+".sh"
	    filename="SubmitScripts_DiEMu/submit_"+str(d.attrib['name'])+str(N_file)+".sh"
	    # copy a skeleton file that set up the code environment, the wall time and the queue
            shutil.copyfile("submitSkeleton.sh", filename)
            # append to command to be run at the end of the skeleton
            outfile = open (filename, 'a')
            #print >> outfile, commandString, topTrees[f], " 0" , " 2000000" , " " , str(N_file)
	    print >> outfile, commandString, "dcap://maite.iihe.ac.be:"+topTrees[f], " " , str(N_file) , " 0" , " 2000000" 
            N_file=N_file+1
