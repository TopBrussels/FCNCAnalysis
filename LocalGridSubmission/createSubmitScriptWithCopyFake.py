import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
from shutil import copyfile
from datetime import datetime
# libray to copy files
import shutil 


# Define time variable                                                           
now = datetime.now()
dd = str(now.day)
mm = str(now.month)
yyyy = str(now.year)
hh = str(now.hour)
mn= str(now.minute)

# make a data string. Pick one of the two above                                                      
#date = dd+"_"+mm+"_"+yyyy+"_"+hh+"h"+mn+"min"
date = yyyy+mm+dd
#date = dd+"_"+mm+"_"+yyyy+"noTrig"

#channels = ["MuMu","ElEl"] 
#channels = ["mumumu","eee","all"] 
channels = ["all"]
doFakes = 1;
JES = 1; 
JER = 1; 
doJESJER = 0;

#if(doJESJERshift == 1) postfix = "_JESdown" ;
#    if(doJESJERshift == 2) postfix = "_JESup" ;
#    if(doJESJERshift == 3) postfix = "_JERdown" ;
#    if(doJESJERshift == 4) postfix = "_JERup" ;


# loop over channels
for chan in channels:
    tree = ET.ElementTree(file='../config/Run2TriLepton_samples.xml')
    print "\nSearching list of sample used for ", chan, " channel!"
    # getting the appropriate xml file
    if "mumumu" in chan:
        tree = ET.ElementTree(file='../config/Run2TriLepton_samples.xml')
#        tree = ET.ElementTree(file='../config/test.xml')
    elif "eee" in chan:
        tree = ET.ElementTree(file='../config/Run2TriLepton_samples.xml')
    elif "all" in chan:
        tree = ET.ElementTree(file='../config/Run2TriLepton_samples.xml')
    elif "ElMu" in chan:
        tree = ET.ElementTree(file='../config/FullSamplesElMuV9.xml')
    else:
        print "Channel '", chan , "' is not a correct channel name. No tree has been loaded!"
        #sys.exit()
    
    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"
    
    # create new dirs if not already existing
    if not os.path.exists("SubmitScripts/"+date):
        os.makedirs("SubmitScripts/"+date)
    if not os.path.exists("SubmitScripts/"+date+"/"+"fakes"):
        os.makedirs("SubmitScripts/"+date+"/"+"fakes")
    if not os.path.exists("SubmitScripts/"+date+"/"+"fakes"+"/output"):
        os.makedirs("SubmitScripts/"+date+"/"+"fakes"+"/output")
    #if not os.path.exists("SubmitScripts/"+date+"/"+chan+"/test"):
    #    os.makedirs("SubmitScripts/"+date+"/"+chan+"/test")
    #if not os.path.exists("SubmitScripts/"+date+"/fakes"):
#	os.makedirs("SubmitScripts/"+date+"/fakes")
#    if not os.path.exists("SubmitScripts/"+date+"/fakes/output"):
#        os.makedirs("SubmitScripts/"+date+"/fakes/output")
	
    # copy the submitAll macro
    copyfile("SubmitAll.sh","SubmitScripts/"+date+"/"+"fakes"+"/SubmitAll.sh")
    #copyfile("SubmitAll.sh","SubmitScripts/"+date+"/fakes/SubmitAll.sh")
    
    # list of variables 
    topTrees = []
    listOfFiles = []
    files_str=""
    FilePerJob=0
    addPrefix=True
    N_processed=0
    listOfScratchFiles = []
    listOfTmpDirFiles = []
    CopyCmdlistOfFiles = []
    scractFiles_str=""
    tmpdirFiles_str=""
  

    
    # loop over all the dataset with add="1"
    for d in datasets:
        if d.attrib['add'] == '1' and  not "80X" in str(d.attrib['name']):
            print "found dataset to be added..." + str(d.attrib['name'])
            commandString = "./Ntupler "+str(d.attrib['name'])+" "+str(d.attrib['title'])+" "+str(d.attrib['add'])+" "+str(d.attrib['color'])+" "+str(d.attrib['ls'])+" "+str(d.attrib['lw'])+" "+str(d.attrib['normf'])+" "+str(d.attrib['EqLumi'])+" "+str(d.attrib['xsection'])+" "+str(d.attrib['PreselEff'])
            topTrees = glob.glob(d.attrib['filenames'])

            # setting the number of file per job depending whether it is data sample or not
            # this ca be tweaked
            if "data" in str(d.attrib['name']):
                FilePerJob=80
            else:
                FilePerJob=50
            
            N_job = 0
            N_file = 1
            remainder= len(topTrees)%FilePerJob
#            print "remainder is", remainder
            
#            print "len(topTrees) is ", len(topTrees)
            # loop over all the root files 
            for f in range(0,len(topTrees)):
#                print "file number ", f , " is : ", topTrees[f]               

                # Combine multiple root files in a single job
                listOfFiles.append(topTrees[f])
                #temp copy 
                CopyCmdlistOfFiles.append("dccp dcap://maite.iihe.ac.be:"+topTrees[f]+" /$TMPDIR/TOPTREE_"+str(f)+".root")
                listOfScratchFiles.append(" /scratch/$PBS_JOBID/TOPTREE_"+str(f)+".root")
                listOfTmpDirFiles.append(" /$TMPDIR/TOPTREE_"+str(f)+".root")            


                # if the number of files is big enough, create one job with the list of files
                if (len(listOfFiles) == FilePerJob) or ((len(topTrees)- N_job * FilePerJob <= FilePerJob) and (len(listOfFiles) == remainder) ):
#                    print "len(listOfFiles) is ", len(listOfFiles) 
                    

                    # create a file for this job
                    filename="SubmitScripts/"+date+"/"+"fakes"+"/submit_"+str(d.attrib['name'])+"_"+str(N_job*FilePerJob+1)+"to"+str(N_job*FilePerJob+len(listOfFiles))+".sh"
		    #filename="SubmitScripts/"+date+"/fakes/submit_"+str(d.attrib['name'])+"_"+str(N_job*FilePerJob+1)+"to"+str(N_job*FilePerJob+len(listOfFiles))+".sh"
                    # copy a skeleton file that set up the code environment, the wall time and the queue
                    if "data" in str(d.attrib['name']):
	                shutil.copyfile("submitSkeletondata.sh", filename)
		    else: 
			shutil.copyfile("submitSkeleton.sh", filename)
                    # append to the file the actual command
                    outfile = open (filename, 'a')

                     # Loop over the files of the current job
                    for fpj in range (0,len(listOfFiles)):
#                        print listOfFiles[fpj]

                        # add prefix if need
                        if (addPrefix == True):
                            listOfFiles[fpj]="dcap://maite.iihe.ac.be"+listOfFiles[fpj]
                        # string contain the list of files separated by a space
                        files_str=files_str+ " " + listOfFiles[fpj]
                        scractFiles_str=scractFiles_str+ " " + listOfScratchFiles[fpj]
                        tmpdirFiles_str=tmpdirFiles_str+ " " + listOfTmpDirFiles [fpj]
                        N_processed=N_processed+1
                        # copy all the file
                        # print >> outfile , CopyCmdlistOfFiles[fpj]



                    print >> outfile, commandString, files_str, " ", JES, " " , JER, " " , doFakes, " ", doJESJER , " " , str(N_job+1) , " 0" , " 200000000" 

                    # cleaning
                    listOfFiles=[]
                    files_str=""
		    listOfScratchFiles=[]
                    CopyCmdlistOfFiles=[]
                    listOfTmpDirFiles =[]
		    scractFiles_str=""

                    N_job=N_job+1
#                    print N_job * FilePerJob
#                    print "Number of processed file is ", N_processed

                N_file=N_file+1



#                print lisfOflisOfFiles
                
# moving the newly created dir
#os.chdir("SubmitScripts/"+chan+"/"+date)
