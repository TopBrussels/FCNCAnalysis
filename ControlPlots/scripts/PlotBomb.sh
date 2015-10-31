#Need to set the channel in MakePlot and makeControlPlotsMerger
#Need to set the databool to true in MakePlot
echo "*** Make mergerscripts ***"
python makeControlPlotsMerger.py
echo "*** merge ***"
sh mergeControlplots_ElEl_.sh 
echo "*** make plots ***"
root -l MakePlot.C
