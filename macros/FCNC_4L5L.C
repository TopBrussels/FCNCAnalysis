#include <iostream>

using namespace std;

void FCNC_4L5L(){

  AutoLibraryLoader::enable();

  gSystem->CompileMacro("TopAnalyzerLite.cc", "k");
 
  //name of your tree
  string treeName = "tree";
  //where you store your information, plots.. 
  string imageOutDir = "../plots/FCNC_4L5L_plot";

  TopAnalyzerLite* analyzer = new TopAnalyzerLite(treeName, imageOutDir);
  
  analyzer->addRealData("45_GluGluHiggs4lep_tree.root", 19700);
  
  analyzer->addMCSig("HctL_HToZZ_ZToLL", "Hct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 5, -1, kRed, false);//cross section = 5 --> needs to be changed
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToNuNu", "Hct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL_tree.root", 5, -1, kRed, false);//cross section = 5 --> needs to be changed
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToBB", "Hct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL_tree.root", 5, -1, kRed, false);//cross section = 5 --> needs to be changed
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToUDC", "Hct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL_tree.root", 5, -1, kRed, false);//cross section = 5 --> needs to be changed
  analyzer->addMCSig("HctL_HToWW_WToLNu", "Hct", "../ntuples/45_TTJetsTocHbW_HToWW_WToLNuL_HctL_tree.root", 5, -1, kRed, false);//cross section = 5 --> needs to be changed
  analyzer->addMCSig("ZctL", "Zct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL_tree.root", 5, -1, kCyan, true);//cross section = 5 --> needs to be changed
  
  //analyzer->addMCSig("HctR", "Hct", "../ntuples/45_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 5, -1, kRed, false);
  
  analyzer->addMCBkg("ggF", "ggF", "../ntuples/45_GluGluHiggs4lep_tree.root", 19.27, -1, kBlue, ""); //cross section = 19.27
  analyzer->addMCBkg("WZ_To3LNu", "WZ", "../ntuples/45_WZ_To3LNu_tree.root", 1.087, -1, kGreen, ""); //cross section = 1.087
  //analyzer->addMCBkg("ZZ_To4L", "ZZ", "../ntuples/45_ZZ_To4L_tree.root", 0.18, -1, kViolet, "");  // //cross section = 0.18
  
  //analyzer->addMonitorPlot("name for cuts", "var in branch", "title, title x, title y", #bins, xmin, xmax, ymin, ymax,doLogy)
  analyzer->addMonitorPlot("nJets", "nJets", "Jet Multiplicity;Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("nLJets_CSVM", "nLJets", "Light Jet Multiplicity;Light Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("nBJets_CSVM", "nBJets", "B Jet Multiplicity;B Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("MET", "missingEt", "MET;missing Et;Events", 50, 0, 500, 0., 1000, true);
  analyzer->addMonitorPlot("OSSFpairs", "InvMass_4lept_Zdecay", "Reconstructed Z boson;Inv. Mass;Events", 50, 0, 300, 0., 1000, true);
  analyzer->addMonitorPlot("FCNC_top", "InvMass_FCNC_top_Zdecay", "Reconstructed top (FCNC);Inv. Mass;Events", 100, 0, 500, 0., 1000, true);
  analyzer->addMonitorPlot("SM_lb", "InvMass_SM_lb", "SM l+b;Inv. Mass;Events", 25, 0, 200, 0., 1000, true);
  analyzer->addMonitorPlot("SM_W", "InvMass_SM_W", "SM W;Inv. Mass;Events", 100, 0, 800, 0., 1000, true);
  analyzer->addMonitorPlot("SM_W_tr","TrMass_W", "SM W;Transverse Mass;Events", 100, 0, 1000, 0., 1000, true); 
  analyzer->addMonitorPlot("SM_top", "InvMass_SM_top", "SM top;Inv. Mass;Events", 25, 0, 850, 0., 1000, true);
  
  //addCutStep("cuts", "histograms, histograms",scaling") ! no spaces
  analyzer->addCutStep("","nJets,nLJets_CSVM,nBJets_CSVM,MET,OSSFpairs,FCNC_top,SM_lb,SM_W,SM_W_tr,SM_top",1);
  
  //writing the plots 
  analyzer->applyCutSteps();

  analyzer->saveHistograms();

}
