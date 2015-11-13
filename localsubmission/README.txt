0) Work in /localgrid/<username>
1) Adapt your .cc code bby adding an extra argument, specifying the jobNum inside your code as well (see FCNC_EventSelection_TreeMaker_localsubmission.cc)
2) copy /user/<username>/lib to /localgrid/<username>
3) Adapt the submitSkeleton.sh & createSubmitScript.py pointing to your own directories and config files
4) make the scripts for submitting to local grid by executing 'python createSubmitScript.py'
5) Submit all the jobs at once by executing 'source SubmitAll.sh'
