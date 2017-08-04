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
date = dd+"_"+mm+"_"+yyyy
#date = dd+"_"+mm+"_"+yyyy+"noTrig"


# pick one of the following
channels = ["Mu"] 
#channels = ["MuMu"] 
#channels = ["ElEl"] 
#channels=["test"]


# loop over channels
for chan in channels:
    print "\nSearching list of sample used for ", chan, " channel!"
    # getting the appropriate xml file
    if "Mu" in chan:
        tree = ET.ElementTree(file='../config/tmp.xml')
    elif "El" in chan:
        tree = ET.ElementTree(file='../config/tmp.xml')
    else:
        print "Channel '", chan , "' is not a correct channel name. No tree has been loaded!"
        sys.exit()

    root = tree.getroot()
    datasets = root.find('datasets')
    print "found  "  + str(len(datasets)) + " datasets"

    
    # create new dirs if not already existing
    if not os.path.exists("SubmitScripts/"+date):
        os.makedirs("SubmitScripts/"+date)
    if not os.path.exists("SubmitScripts/"+date+"/"+chan):
        os.makedirs("SubmitScripts/"+date+"/"+chan)
    if not os.path.exists("SubmitScripts/"+date+"/"+chan+"/output"):
        os.makedirs("SubmitScripts/"+date+"/"+chan+"/output")
    if not os.path.exists("SubmitScripts/"+date+"/"+chan+"/test"):
        os.makedirs("SubmitScripts/"+date+"/"+chan+"/test")

    # copy the submitAll macro
    copyfile("SubmitAll.sh","SubmitScripts/"+date+"/"+chan+"/SubmitAll.sh")
    copyfile("MakeTXTforMassiveSubmission.sh","SubmitScripts/"+date+"/"+chan+"/MakeTXTforMassiveSubmission.sh")

    
    # list of variables 
    topTrees = []
    listOfFiles = []
    listOfScratchFiles = []
    listOfTmpDirFiles = []
    CopyCmdlistOfFiles = []
    files_str=""
    scractFiles_str=""
    tmpdirFiles_str=""
    FilePerJob=0 # add dccap...
    addPrefix=True
    N_processed=0
    Luminosity_=0.;
    

    for d in datasets:
        if "Data" in str(d.attrib['name']):
            Luminosity_=Luminosity_+float(d.attrib['EqLumi'])
    
    # loop over all the dataset with add="1"
    for d in datasets:
        if d.attrib['add'] == '1':
            print "found dataset to be added..." + str(d.attrib['name'])
            commandString = "./TreeMaker_noJEC "+str(d.attrib['name'])+" "+str(d.attrib['title'])+" "+str(d.attrib['add'])+" "+str(d.attrib['color'])+" "+str(d.attrib['ls'])+" "+str(d.attrib['lw'])+" "+str(d.attrib['normf'])+" "+str(d.attrib['EqLumi'])+" "+str(d.attrib['xsection'])+" "+str(d.attrib['PreselEff']) +" "+ str(Luminosity_)
            topTrees = glob.glob(d.attrib['filenames'])

            # setting the number of file per job depending whether it is data sample or not
            # this ca be tweaked
            if "Data" in str(d.attrib['name']):
                FilePerJob= 20
            else:
                FilePerJob=1

            # create a test job for each dataset
            # create a file for this job 
            filenameTest="SubmitScripts/"+date+"/"+chan+"/test"+"/submit_"+str(d.attrib['name'])+"_"+"Test"+".sh"
            # copy a skeleton file that set up the code environment, the wall time and the queue 
            shutil.copyfile("submitTestSkeleton.sh", filenameTest)
            # append to the file the actual command 
            outfileTest = open (filenameTest, 'a')
            if not len(topTrees) == 0:
                print >> outfileTest, commandString, "dcap://maite.iihe.ac.be:"+topTrees[0], " ", chan , " " , 1 , " 0" , " 10000"
                
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
                CopyCmdlistOfFiles.append("dccp dcap://maite.iihe.ac.be:"+topTrees[f]+" /$TMPDIR/TOPTREE_"+str(f)+".root")
                listOfScratchFiles.append(" /scratch/$PBS_JOBID/TOPTREE_"+str(f)+".root")
                listOfTmpDirFiles.append(" /$TMPDIR/TOPTREE_"+str(f)+".root")
                
#                print CopyCmdlistOfFiles[0]
                
                # if the number of files is big enough, create one job with the list of files
                if (len(listOfFiles) == FilePerJob) or ((len(topTrees)- N_job * FilePerJob <= FilePerJob) and (len(listOfFiles) == remainder) ):
#                    print "len(listOfFiles) is ", len(listOfFiles) 

                    # create a file for this job
                    filename="SubmitScripts/"+date+"/"+chan+"/submit_"+str(d.attrib['name'])+"_"+str(N_job*FilePerJob+1)+"to"+str(N_job*FilePerJob+len(listOfFiles))+".sh"
                    # copy a skeleton file that set up the code environment, the wall time and the queue
                    shutil.copyfile("submitSkeleton.sh", filename)

                    # add one copy cmd per file in that job
                    outfile = open (filename, 'a')
                    print >> outfile , "#first make all the copies"

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
                        print >> outfile , CopyCmdlistOfFiles[fpj]
                        
                    print >> outfile , "\n\n"

#                    print files_str


                    # run on the files
                    print >> outfile, "# now run on the file copied under /$TMPDIR/ "
                    print >> outfile, commandString, scractFiles_str , " ", chan , " " , str(N_job+1) , " 0" , " 10000000" 

                    # cleaning
                    listOfFiles=[]
                    listOfScratchFiles=[]
                    CopyCmdlistOfFiles=[]
                    listOfTmpDirFiles =[]
                    files_str=""
                    scractFiles_str=""

                    N_job=N_job+1
#                    print N_job * FilePerJob
#                    print "Number of processed file is ", N_processed

                N_file=N_file+1



#                print lisfOflisOfFiles
                
# moving the newly created dir
#os.chdir("SubmitScripts/"+chan+"/"+date)

