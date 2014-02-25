#include <iostream>

using namespace std;

void FCNC_3L(){

  AutoLibraryLoader::enable();

  gSystem->CompileMacro("TopAnalyzerLite.cc", "k");
 
   //scaling for signal
  double scaler = 5.; 
 
  //name of your tree
  string treeName = "tree";
  //where you store your information, plots.. 
  string imageOutDir = "../plots/FCNC_3L_plot";

  TopAnalyzerLite* analyzer = new TopAnalyzerLite(treeName, imageOutDir);
  
  analyzer->addRealData("../ntuples/3L_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 19700);
  
  analyzer->addMCSig("HctL_HToZZ_ZToLL", "Hct_ZZ", "../ntuples/3L_TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root", 0.00016516*scaler, -1, kRed, false);//cross section = 0.00016516
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToNuNu", "Hct_ZZ", "../ntuples/3L_TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL_tree.root", 0.00067929*scaler, -1, kRed, false);//cross section = 0.00067929
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToBB", "Hct_ZZ", "../ntuples/3L_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL_tree.root", 0.0005135*scaler, -1, kRed, false);//cross section = 0.0005135
  analyzer->addMCSig("HctL_HToZZ_ZToLL_ZToUDC", "Hct_ZZ", "../ntuples/3L_TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL_tree.root", 0.0018609*scaler, -1, kRed, false);//cross section = 0.0018609 --> needs to be changed
  analyzer->addMCSig("HctL_HToWW_WToLNu", "Hct_WW", "../ntuples/3L_TTJetsTocHbW_HToWW_WToLNuL_HctL_tree.root", 0.022659*scaler, -1, kRed+2, false);//cross section = 0.022659
  analyzer->addMCSig("HctL_HToWW_WToLNu_WToJets", "Hct_WW", "../ntuples/3L_TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL_tree.root", 0.090636*scaler, -1, kRed+2, false);//cross section = 0.090636
  analyzer->addMCSig("ZctL", "Zct", "../ntuples/3L_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL_tree.root", 0.1575*scaler, -1, kRed+4,false);//cross section = 0.1575
  
  analyzer->addMCBkg("HctL_HToWW_WToLNu", "Hct_WW", "../ntuples/3L_TTJetsTocHbW_HToWW_WToLNuL_HctL_tree.root", 0.022659*scaler, -1, kRed+2);//cross section = 0.022659
  
  /*
  analyzer->addMCBkg("ggF", "SM_H", "../ntuples/3L_GluGluHiggs4lep_tree.root", 19.27, -1, kBlue, ""); //cross section =  19.27pb
  analyzer->addMCBkg("VBF", "SM_H", "../ntuples/3L_VBHiggs4lep_tree.root", 0.00414341777, -1, kBlue, "");  // //cross section = 0.00414341777
  
  analyzer->addMCBkg("WW_To2L2Nu", "WW", "../ntuples/3L_WW_To2L2Nu_tree.root", 5.757, -1, kOrange, ""); //cross section =  5.757
  analyzer->addMCBkg("WZ_To2L2Q", "WZ", "../ntuples/3L_WZ_To2L2Q_tree.root", 2.267, -1, kOrange, ""); //cross section = 2.267
  analyzer->addMCBkg("WZ_To3LNu", "WZ", "../ntuples/3L_WZ_To3LNu_tree.root", 1.087, -1, kOrange, ""); //cross section = 1.087
  analyzer->addMCBkg("ZZ_To4L", "ZZ", "../ntuples/3L_ZZ_To4L_tree.root", 0.18, -1, kYellow, "");  // //cross section = 0.18
  analyzer->addMCBkg("ZZ_To2L2Nu", "ZZ", "../ntuples/3L_ZZ_To2L2Nu_tree.root", 0.713, -1, kYellow, "");  // //cross section = 0.713
  analyzer->addMCBkg("ZZ_To2L2Q", "ZZ", "../ntuples/3L_ZZ_To2L2Q_tree.root", 2.492, -1, kYellow, "");  // //cross section = 2.492
  
  analyzer->addMCBkg("ttH", "ttH", "../ntuples/3L_ttH_tree.root", 0.1293, -1, kViolet, "");  // //cross section = 0.1293
  analyzer->addMCBkg("ttZ", "ttZ", "../ntuples/3L_TTZ_tree.root", 0.172, -1, kPink, "");  // //cross section = 0.172
  analyzer->addMCBkg("ttW", "ttW", "../ntuples/3L_TTW_tree.root", 0.2148, -1, kPink+9, "");  // //cross section = 0.2148
  analyzer->addMCBkg("tbZ", "tbZ", "../ntuples/3L_TBZ_ToLL_4F_tree.root", 0.0114, -1, kPink+9, "");  // //cross section = 0.2148
  
  analyzer->addMCBkg("DY_1Jet", "DY", "../ntuples/3L_Z_1Jets_tree.root", 671.83, -1, kPink+9, "");  // //cross section = 671.83
  analyzer->addMCBkg("DY_2Jet", "DY", "../ntuples/3L_Z_2Jets_tree.root", 216.76, -1, kPink+9, "");  // //cross section = 216.76
  analyzer->addMCBkg("DY_3Jet", "DY", "../ntuples/3L_Z_3Jets_tree.root", 61.20, -1, kPink+9, "");  // //cross section = 61.20
  analyzer->addMCBkg("DY_4Jet", "DY", "../ntuples/3L_Z_4Jets_tree.root", 27.59, -1, kPink+9, "");  // //cross section = 27.59
  
  
  analyzer->addMCBkg("TT_FullLept", "tt", "../ntuples/3L_TT_FullLeptMGDecay_tree.root", 26.42, -1, kPink+9, "");  // //cross section = 26.42
  analyzer->addMCBkg("TT_SemiLept", "tt", "../ntuples/3L_TT_SemiLeptMGDecay_tree.root", 110.26, -1, kPink+9, "");  // //cross section = 110.26
  
  analyzer->addMCBkg("ST_T_s-ch", "ST", "../ntuples/3L_ST_T_s-ch_tree.root", 3.79, -1, kPink+9, "");  // //cross section = 3.79
  analyzer->addMCBkg("ST_Tbar_s-ch", "ST", "../ntuples/3L_STbar_T_s-ch_tree.root", 1.76, -1, kPink+9, "");  // //cross section = 1.76
  analyzer->addMCBkg("ST_T_tW-ch", "ST", "../ntuples/3L_ST_T_tW-ch_tree.root", 11.1, -1, kPink+9, "");  // //cross section = 11.1
  analyzer->addMCBkg("ST_Tbar_tW-ch", "ST", "../ntuples/3L_STbar_T_tW-ch_tree.root", 11.1, -1, kPink+9, "");  // //cross section = 11.1
  
  */
  
  //analyzer->addMonitorPlot("name for cuts", "var in branch", "title, title x, title y", #bins, xmin, xmax, ymin, ymax,doLogy)
  analyzer->addMonitorPlot("nJets", "nJets", "Jet Multiplicity;Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("nLJets_CSVM", "nLJets", "Light Jet Multiplicity;Light Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("nBJets_CSVM", "nBJets", "B Jet Multiplicity;B Jet Multiplicity;Events", 15, -0.5, 14.5, 0., 1000, true);
  analyzer->addMonitorPlot("MET", "missingEt", "MET;missing Et;Events", 50, 0, 500, 0., 1000, true);
  
  //addCutStep("cuts", "histograms, histograms",scaling") ! no spaces
  analyzer->addCutStep("","nJets,nLJets_CSVM,nBJets_CSVM,MET,OSSFpairs,FCNC_top,SM_lb,SM_W,SM_W_tr,SM_top",1);
  
  //writing the plots 
  analyzer->applyCutSteps();

  analyzer->saveHistograms();

}
