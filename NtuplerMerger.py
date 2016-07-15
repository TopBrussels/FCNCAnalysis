from ROOT import TChain
import ROOT
import glob
import xml.etree.cElementTree as ET
import os
from datetime import datetime


# Define time variable
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
# pick one of the two above
date = "160705"
#date = "17_1_2016"

#channels = ["_MuMu","_ElEl"]
#channels = ["_ElEl"]
channels = ["_All"] #,"_ElElEl"]

for chan in channels:
    
    #Define path where ntuples are stored
    pathNonMerged = "NtupleMakerOutput/Ntuples"+chan+"/Ntuples_"+date+"/"
    if "MuMuMu" in chan:
        pathMerged = "NtupleMakerOutput/MergedTuples/mumumu/"+"160530"+"/"
    if "ElElEl" in chan:
        pathMerged = "NtupleMakerOutput/MergedTuples/eee/"+"160530"+"/"
    if "All" in chan:
        pathMerged = "NtupleMakerOutput/MergedTuples/all/"+date+"/"
    
    if not os.path.exists(pathMerged):
        os.makedirs(pathMerged)
    
    # get filenames from the xml!!!    
    if "MuMuMu" in chan:
        tree = ET.ElementTree(file='config/Run2TriLepton_samples_extra.xml')
    elif "ElEl" in chan:
        tree = ET.ElementTree(file='config/Run2TriLepton_samples_eee.xml')
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='config/FullSamplesElMuV0.xml')
    elif "All" in chan:
        tree = ET.ElementTree(file='config/Run2TriLepton_samples.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

    #tree = ET.ElementTree(file='config/FullMcBkgdSamplesV9.xml')
    #tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')
    #tree = ET.ElementTree(file='config/DataSamples.xml')

    # get the list of dataset
    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    datasetNames = []
    
    print ""
    # loop over the datasets to be added and fill the "topTrees" vector
    for d in datasets:
        if d.attrib['add'] == '1':
            print "found dataset to be added..." + str(d.attrib['name'])

            # select a subset of the existing root file
            if not "over" in str(d.attrib['name']) :
                datasetNames.append(str(d.attrib['name']))
                print str(d.attrib['name'])
    
    
    listOfZombie= []
    
    # loop over data set to search root files
    for n in datasetNames:
        filenames = glob.glob(pathNonMerged + "/*" + chan + "*" + n + "*.root")
        hadd = "hadd -f " + pathMerged + "/"+ n + ".root"

        if (len(filenames) == 0):
            print "no root files found in directory" , pathNonMerged ,  " for dataset " , n , " !!"
        else :
            # loop over root files
            for f in filenames:
                file=ROOT.TFile(f,"read")
                # check if the file is a zombie
                if (file.IsZombie()):
                    print "File" , f, "is a Zombie.... Skipping"
                    listOfZombie.append(f)
                else:
                    print f
                    hadd = hadd + " " + f
            print "Merging ntuples for " + n
            os.system(hadd)
        
    print "\n\n"
    
    # print the list of zombies
    print "The total number of zombie file is ", len(listOfZombie)
    if (len(listOfZombie) > 0):
        outfile = open (pathMerged+"/Zombie"+chan+".txt", 'a')
        print "And the list of the zombie is put in "+pathMerged+"/Zombie"+chan+".txt "
        for zombie in listOfZombie:
            print >> outfile, zombie
    
    mergeData=True
    
    if (mergeData):
    # combining all the Data in one
        dataList=glob.glob(pathMerged+"*Data*.root")
    
	cmd = "hadd " + pathMerged + "/"+ "Data.root"
        for data in dataList:
            cmd = cmd + " " + data
        os.system(cmd)
            
