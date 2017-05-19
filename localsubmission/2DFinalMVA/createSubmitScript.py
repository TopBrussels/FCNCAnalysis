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

categories = [" 2 3 "," 2 4 "," 3 3 "," 3 4 "," 4 4 "]
catnames = ["b2j3","b2j4","b3j3","b3j4","b4j4"]
couplingVals = [" 0 "," 3 "," 6 "," 9 "," 12 "," 15 ", " 18 ", " 21 ", " 24 "]

for i in range(0,len(categories)):

  count = 0

  for khut in couplingVals:
    for khct in couplingVals:
    
      count = count + 1

      if not os.path.exists("SubmitScripts/"):
        os.makedirs("SubmitScripts/")
      if not os.path.exists("SubmitScripts/"+date):
        os.makedirs("SubmitScripts/"+date)


      filename="SubmitScripts/"+date+"/submit_"+catnames[i]+"_"+str(count)+".sh"
      # copy a skeleton file that set up the code environment, the wall time and the queue 
      shutil.copyfile("submitSkeleton.sh", filename)
            
      commandString = "./TreeProcessor_FinalMVA "+categories[i]+"2D _All _12_5_2017 0 1 1 0 0 "+khut+khct

      outfile = open (filename, 'a')
      print >> outfile , "#first make all the copies"
      print >> outfile , commandString

