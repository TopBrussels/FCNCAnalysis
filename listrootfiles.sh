
echo "Listing files from  pnfs directory"
for j in $( lcg-ls -D srmv2 srm://maite.iihe.ac.be:8443/srm/managerv2?SFN=/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/DoubleEG/crab_DoubleEG-Run2015D-05Oct2015-v1-CMSSW_74X_v8-74X_dataRun2_Prompt_v2/151020_153539/0000/) ; do

#	echo " in loop "
	fullfile=$j
#        echo $fullfile
	myfilename="${fullfile##*/}"
#	echo $myfilename
	filename="${fullfile##*_}"
#	echo $filename
	extension="${fullfile##*.}"
#	echo $extension
        if [   $extension == "root" ] && [ $filename != "matching.root" ]; then 
#		echo "It is a rootfile"
		echo -n "dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/DoubleEG/crab_DoubleEG-Run2015D-05Oct2015-v1-CMSSW_74X_v8-74X_dataRun2_Prompt_v2/151020_153539/0000/"$myfilename","
	fi
#	echo " out loop " 
done

echo "Done listing files from pnfs directory"

#
