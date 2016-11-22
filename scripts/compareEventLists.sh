#outdir=~/www/30-10-2016/FCNC_sync/
outdir=~/www/14-11-2016/FCNC_sync/

rm -rf mine his
mkdir mine his

#FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/30-10-2016/FCNC_sync_2/EventInfo_el.txt
#FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_el.txt
#FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El.txt
FILE1=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El.txt
FILE2=http://kskovpen.web.cern.ch/kskovpen/tHFCNC/outputEl.txt

wget $FILE1 -P mine
wget $FILE2 -P his

FILE1=`basename $FILE1`
FILE2=`basename $FILE2`

rm prompt_el.txt
echo `cat mine/$FILE1 | wc -l` `cat his/$FILE2 | wc -l` > prompt_el.txt

echo "run,lumi,id,lepId,lepPt,lepEta,lepPhi,metpt,metphi,njets,nbjets" > mine/${FILE1/.txt/.csv}
cat mine/$FILE1 >> mine/${FILE1/.txt/.csv}
echo "run,lumi,id,lepId,lepPt,lepEta,lepPhi,metpt,metphi,njets,nbjets" > his/${FILE2/.txt/.csv}
cat his/$FILE2 >> his/${FILE2/.txt/.csv}

for i in mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv}; do
    for j in `seq 1 10`; do
        sed -i ':a;N;$!ba;s/\n /\n/g' $i
        sed -i ':a;N;$!ba;s/  / /g' $i
    done
    sed -i ':a;N;$!ba;s/ \n/\n/g' $i
    sed -i ':a;N;$!ba;s/ /,/g' $i
    sed -i '$ s/.$//' $i
done

python compare.py electrons mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv} -1 >> prompt_el.txt


#FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/30-10-2016/FCNC_sync_2/EventInfo_mu.txt
FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_mu.txt
#FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_mu.txt
FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_mu.txt
#FILE2=http://test-cirkovic.web.cern.ch/test-cirkovic/24-10-2016/FCNC_sync_2/EventInfo_mu.txt

wget $FILE1 -P mine
wget $FILE2 -P his

FILE1=`basename $FILE1`
FILE2=`basename $FILE2`

rm prompt_mu.txt
echo `cat mine/$FILE1 | wc -l` `cat his/$FILE2 | wc -l` > prompt_mu.txt

echo "run,lumi,id,lepId,lepPt,lepEta,lepPhi,metpt,metphi,njets,nbjets" > mine/${FILE1/.txt/.csv}
cat mine/$FILE1 >> mine/${FILE1/.txt/.csv}
echo "run,lumi,id,lepId,lepPt,lepEta,lepPhi,metpt,metphi,njets,nbjets" > his/${FILE2/.txt/.csv}
cat his/$FILE2 >> his/${FILE2/.txt/.csv}


for i in mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv}; do
    for j in `seq 1 10`; do
        sed -i ':a;N;$!ba;s/\n /\n/g' $i
        sed -i ':a;N;$!ba;s/  / /g' $i
    done
    sed -i ':a;N;$!ba;s/ \n/\n/g' $i
    sed -i ':a;N;$!ba;s/ /,/g' $i
    sed -i '$ s/.$//' $i
done

#sed -i '$ d' his/${FILE2/.txt/.csv} # last line incomplete

python compare.py muons mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv} -1 >> prompt_mu.txt



#FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/30-10-2016/FCNC_sync_2/EventInfo_mu.txt
FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_el_vars.txt
#FILE1=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El_extraElectronVariables.txt
#FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_mu.txt
FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El_extraElectronVariables.txt
#FILE2=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_el_vars.txt
#FILE2=http://test-cirkovic.web.cern.ch/test-cirkovic/24-10-2016/FCNC_sync_2/EventInfo_mu.txt

wget $FILE1 -P mine
wget $FILE2 -P his

FILE1=`basename $FILE1`
FILE2=`basename $FILE2`

rm prompt_el_vars.txt
echo `cat mine/$FILE1 | wc -l` `cat his/$FILE2 | wc -l` > prompt_el_vars.txt

echo "run,lumi,id,superClusterEta,deltaEtaIn,deltaPhiIn,sigmaIEtaIEta_full5x5,hadronicOverEm,ioEmIoP,pfElectronIso,missingHits" > mine/${FILE1/.txt/.csv}
cat mine/$FILE1 >> mine/${FILE1/.txt/.csv}
echo "run,lumi,id,superClusterEta,deltaEtaIn,deltaPhiIn,sigmaIEtaIEta_full5x5,hadronicOverEm,ioEmIoP,pfElectronIso,missingHits" > his/${FILE2/.txt/.csv}
cat his/$FILE2 >> his/${FILE2/.txt/.csv}


for i in mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv}; do
    for j in `seq 1 10`; do
        sed -i ':a;N;$!ba;s/\n /\n/g' $i
        sed -i ':a;N;$!ba;s/  / /g' $i
    done
    sed -i ':a;N;$!ba;s/ \n/\n/g' $i
    sed -i ':a;N;$!ba;s/ /,/g' $i
    sed -i '$ s/.$//' $i
done

#sed -i '$ d' his/${FILE2/.txt/.csv} # last line incomplete

python compare_el_vars.py electrons mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv} -1 >> prompt_el_vars.txt


#FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/30-10-2016/FCNC_sync_2/EventInfo_mu.txt
FILE1=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_el_vars_1.txt
#FILE1=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El_extraElectronVariables.txt
#FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_mu.txt
FILE2=http://mon.iihe.ac.be/~kderoove/FCNC_SynchExercise/EventInfo_El_extraElectronVariables_IsoVars.txt
#FILE2=http://test-cirkovic.web.cern.ch/test-cirkovic/14-11-2016/FCNC_sync_2/EventInfo_el_vars.txt
#FILE2=http://test-cirkovic.web.cern.ch/test-cirkovic/24-10-2016/FCNC_sync_2/EventInfo_mu.txt

wget $FILE1 -P mine
wget $FILE2 -P his

FILE1=`basename $FILE1`
FILE2=`basename $FILE2`

rm prompt_el_vars_1.txt
echo `cat mine/$FILE1 | wc -l` `cat his/$FILE2 | wc -l` > prompt_el_vars_1.txt

echo "run,lumi,id,superClusterEta,deltaEtaIn,deltaPhiIn,sigmaIEtaIEta_full5x5,hadronicOverEm,ioEmIoP,pfElectronIso,missingHits,sumChargedHadronPt,sumNeutralHadronEt,sumPhotonEt,rho,Aeff" > mine/${FILE1/.txt/.csv}
cat mine/$FILE1 >> mine/${FILE1/.txt/.csv}
echo "run,lumi,id,superClusterEta,deltaEtaIn,deltaPhiIn,sigmaIEtaIEta_full5x5,hadronicOverEm,ioEmIoP,pfElectronIso,missingHits,sumChargedHadronPt,sumNeutralHadronEt,sumPhotonEt,rho,Aeff" > his/${FILE2/.txt/.csv}
cat his/$FILE2 >> his/${FILE2/.txt/.csv}


for i in mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv}; do
    for j in `seq 1 10`; do
        sed -i ':a;N;$!ba;s/\n /\n/g' $i
        sed -i ':a;N;$!ba;s/  / /g' $i
    done
    sed -i ':a;N;$!ba;s/ \n/\n/g' $i
    sed -i ':a;N;$!ba;s/ /,/g' $i
    sed -i '$ s/.$//' $i
done

#sed -i '$ d' his/${FILE2/.txt/.csv} # last line incomplete

python compare_el_vars_1.py electrons mine/${FILE1/.txt/.csv} his/${FILE2/.txt/.csv} -1 >> prompt_el_vars_1.txt




cp prompt_el.txt prompt_mu.txt prompt_el_vars.txt prompt_el_vars_1.txt $outdir
cp -rf mine his $outdir


