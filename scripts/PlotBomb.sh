#Need to set the channel in MakePlot and makeControlPlotsMerger
#Need to set the databool to true in MakePlot
echo "*** merge ***"
python Merger.py
echo "*** make plots ***"
root -l MakePlot.C
