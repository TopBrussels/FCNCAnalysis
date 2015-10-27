
echo "Listing files from  pnfs directory"
for j in $( lcg-ls -D srmv2 srm://maite.iihe.ac.be:8443/srm/managerv2?SFN=/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_161555/0000/) ; do

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
		echo -n "dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_161555/0000/"$myfilename","
	fi
#	echo " out loop " 
done

echo "Done listing files from pnfs directory"

#
