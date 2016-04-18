rm ../OutPutHistos/merged_ElEl_AllSamples/merged_Data_Run2015D.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_Data_Run2015D.root ../OutPutHistos/_ElEl_Data_Run2015D/*.root

rm ../OutPutHistos/merged_ElEl_AllSamples/merged_DYJets50.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_DYJets50.root ../OutPutHistos/_ElEl_DYJets50/*.root
rm ../OutPutHistos/merged_ElEl_AllSamples/merged_DYJets10-50.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_DYJets10-50.root ../OutPutHistos/_ElEl_DYJets10-50/*.root
rm ../OutPutHistos/merged_ElEl_AllSamples/merged_WJets.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_WJets.root ../OutPutHistos/_ElEl_WJets/*.root
rm ../OutPutHistos/merged_ElEl_AllSamples/merged_TTJets.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_TTJets.root ../OutPutHistos/_ElEl_TTJets/*.root
rm ../OutPutHistos/merged_ElEl_AllSamples/merged_WZ.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_WZ.root ../OutPutHistos/_ElEl_WZ/*.root
rm ../OutPutHistos/merged_ElEl_AllSamples/merged_WW.root
hadd -f ../OutPutHistos/merged_ElEl_AllSamples/merged_WW.root ../OutPutHistos/_ElEl_WW/*.root
