mkdir ../_ElEl
hadd -f ../_ElEl/merged_ElEl_diboson.root ../_ElEl_allSamples//merged_ElEl_WW.root ../_ElEl_allSamples//merged_ElEl_ZZ.root ../_ElEl_allSamples//merged_ElEl_WZ.root
hadd -f ../_ElEl/merged_ElEl_ST.root ../_ElEl_allSamples//merged_ElEl_ST*.root
hadd -f ../_ElEl/merged_ElEl_Zjets.root ../_ElEl_allSamples//merged_ElEl_Zjets*.root
hadd -f ../_ElEl/merged_ElEl_data.root ../_ElEl_allSamples//merged_ElEl_data.root
hadd -f ../_ElEl/merged_ElEl_ttV.root ../_ElEl_allSamples//merged_ElEl_ttZ.root ../_ElEl_allSamples//merged_ElEl_ttW.root 
hadd -f ../_ElEl/merged_ElEl_ttbar.root ../_ElEl_allSamples//merged_ElEl_ttbar.root
hadd -f ../_ElEl/merged_ElEl_tZq.root ../_ElEl_allSamples//merged_ElEl_tZq.root
