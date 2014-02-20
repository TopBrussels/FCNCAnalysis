#include <iostream>
using namespace std;

void FCNC_photon(){

  AutoLibraryLoader::enable();

  gSystem->CompileMacro("TopAnalyzerLite.cc", "k");

  string treeName = "tree";
  string imageOutDir = "FCNC_photon_plot";

  TopAnalyzerLite* analyzer = new TopAnalyzerLite(treeName, imageOutDir);

  analyzer->addRealData("vallot_Run2012.root", 19700);

  // replace -1 with total number of events if you want to provide total number of MC events manually 
  analyzer->addMCSig("Hct", "Hct to #gamma #gamma", "vallot_tcHtoGG.root", 5, -1, kRed, false); //cross section = 1
  analyzer->addMCBkg("HToGG", "H to #gamma #gamma", "vallot_HtoGG.root", 19.5, -1, kYellow); //cross section = 19.5
  analyzer->addMCBkg("DiPhotons", "DiPhotons", "vallot_diphotons.root", 75.4, -1, kGreen); //cross section = 75.4
  analyzer->addMCBkg("TTJets", "TTJets", "vallot_TTJets.root", 250, -1, kRed+2); //cross section = 75.4
  analyzer->addMCBkg("DY", "ZJets", "vallot_ZJets.root", 3530, -1, kBlue); //cross section = 3.53*10^3

  analyzer->addMonitorPlot("njets", "njets", "Jet Multiplicity;Jet Multiplicity;Events", 15, 0, 15, 0.5, 2000, true);
  analyzer->addMonitorPlot("nbjets_CSVL", "nbjets_CSVL", "b-Jet Multiplicity;b-Jet Multiplicity (CSVL);Events", 5, 0, 5, 0.5, 2000, true);
  analyzer->addMonitorPlot("nbjets_CSVM", "nbjets_CSVM", "b-Jet Multiplicity;b-Jet Multiplicity (CSVM);Events", 5, 0, 5, 0.5, 2000, true);
  analyzer->addMonitorPlot("nbjets_CSVT", "nbjets_CSVT", "b-Jet Multiplicity;b-Jet Multiplicity (CSVT);Events", 5, 0, 5, 0.5, 2000, true);
  analyzer->addMonitorPlot("nphotons", "nphotons", "Photon Multiplicity;Photon Multiplicity;Events", 6, 0, 6, 0.5, 2000, true);
  analyzer->addMonitorPlot("photon1_pt", "photon1_pt", "Leading Photon P_{T};Photon P_{T} (GeV);Events", 20, 0, 200, 0.5, 3000, true);
  analyzer->addMonitorPlot("photon2_pt", "photon2_pt", "Second Leading Photon P_{T};Photon P_{T} (GeV);Events", 20, 0, 200, 0.5, 3000, true);
  analyzer->addMonitorPlot("photon1_eta", "photon1_eta", "Leading Photon #eta;Photon #eta;Events", 60, -3, 3, 0.5, 5000, true);
  analyzer->addMonitorPlot("photon2_eta", "photon2_eta", "Second Leading Photon #eta;Photon #eta;Events", 60, -3, 3, 0.5, 5000, true);
  analyzer->addMonitorPlot("diphoton_mass", "diphoton_mass", "Di-Photon Invariant Mass;Di-Photon Invariant Mass (GeV);Events", 80, 80, 160, 0.5, 2000, true);

  analyzer->addCutStep("nphotons >= 2 && njets >= 2 && photon1_pt > 40 && photon2_pt > 40 && abs(photon1_eta) < 2.5 && abs(photon2_eta) < 2.5 && photon1_relIso < 0.05 && photon2_relIso < 0.05", "nbjets_CSVL,nbjets_CSVM,nbjets_CSVT,njets,nphotons,photon1_pt,photon2_pt,photon1_eta,photon2_eta,diphoton_mass", 0.1);

  analyzer->addCutStep("nphotons >= 2 && njets >= 2 && photon1_pt > 40 && photon2_pt > 40 && abs(photon1_eta) < 2.5 && abs(photon2_eta) < 2.5 && nbjets_CSVT >= 1 && photon1_relIso < 0.1 && photon2_relIso < 0.1", "nbjets_CSVL,nbjets_CSVM,nbjets_CSVT,njets,nphotons,photon1_pt,photon2_pt,photon1_eta,photon2_eta,diphoton_mass", 0.1);

  analyzer->applyCutSteps();

  analyzer->saveHistograms();

}
