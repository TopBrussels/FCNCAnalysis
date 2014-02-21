#include <iostream>

using namespace std;

void FCNC_4L5L(){

  AutoLibraryLoader::enable();

  gSystem->CompileMacro("TopAnalyzerLite.cc", "k");
 
  string treeName = "tree";
  string imageOutDir = "../plots/FCNC_4L5L_plot";

  TopAnalyzerLite* analyzer = new TopAnalyzerLite(treeName, imageOutDir);
  
  analyzer->addRealData("45__GluGluHiggs4lep_tree.root", 19700);
  
  analyzer->addMCSig("HctL", "Hct", "45_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 5, -1, kRed, false);
  analyzer->addMCSig("HctR", "Hct", "45_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 5, -1, kRed, false);
  
  analyzer->addMCBkg("ggF", "ggF", "45__GluGluHiggs4lep_tree.root", 19.27, -1, kBlue, "");
  
  analyzer->addMonitorPlot("nJets", "nJets", "Jet Multiplicity;Jet Multiplicity;Events", 15, 0, 15, 0.5, 2000, true);
  
  analyzer->addCutStep("","nJets",0.1);
  
  analyzer->applyCutSteps();

  analyzer->saveHistograms();

}
