#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include "TStyle.h"
#include "TPaveText.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <TFile.h>
#include <TLeaf.h>
#include <utility>
#include "Style.C"
#include <map>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;


////////////////////////////////// mapping ///////////////////////////////
map<string,TFile*> tFileMap;
TFile *fin;
map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;
vector < Dataset* > datasets;


//////////////////////////// branches //////////////////////////////
// Declaration of leaf types
Float_t         MVA_Zboson_pt;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_charge_asym;
Float_t         MVA_bdiscCSVv2_jet_0;

Float_t         MVA_cdiscCvsB_jet_1;
Float_t         MVA_cdiscCvsB_jet_0;
Float_t         MVA_cdiscCvsL_jet_1;
Float_t         MVA_cdiscCvsL_jet_0;
Float_t         MVA_dRZc;
Float_t         MVA_dRWlepb;
Float_t         MVA_dRZWlep;
Float_t         MVA_mlb;
Float_t         MVA_FCNCtop_M;
Float_t           MVA_nJets_CharmL;
Float_t           MVA_NJets_CSVv2M;
Int_t         MVA_region;
Double_t       MVA_x1;
Double_t       MVA_x2;
Int_t       MVA_id1;
Int_t       MVA_id2;
Double_t       MVA_q;
//Float_t         MVA_weight;
Double_t       MVA_weight_nom;
Int_t         MVA_channel;
Double_t       MVA_BDT;
Double_t        MVA_EqLumi;
Double_t       MVA_weight_puSF_up;
Double_t       MVA_weight_puSF_down;
Double_t       MVA_weight_electronSF_up;
Double_t       MVA_weight_electronSF_down;
Double_t       MVA_weight_muonSF_up;
Double_t       MVA_weight_muonSF_down;
Double_t       MVA_weight_btagSF_cferr1_up;
Double_t       MVA_weight_btagSF_cferr1_down;
Double_t       MVA_weight_btagSF_cferr2_up;
Double_t       MVA_weight_btagSF_cferr2_down;
Double_t       MVA_weight_btagSF_hf_up;
Double_t       MVA_weight_btagSF_hf_down;
Double_t       MVA_weight_btagSF_hfstats1_up;
Double_t       MVA_weight_btagSF_hfstats1_down;
Double_t       MVA_weight_btagSF_hfstats2_up;
Double_t       MVA_weight_btagSF_hfstats2_down;
Double_t       MVA_weight_btagSF_lf_up;
Double_t       MVA_weight_btagSF_lf_down;
Double_t       MVA_weight_btagSF_lfstats1_up;
Double_t       MVA_weight_btagSF_lfstats1_down;
Double_t       MVA_weight_btagSF_lfstats2_up;
Double_t       MVA_weight_btagSF_lfstats2_down;

// List of branches
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_0;   //!

TBranch        *b_MVA_cdiscCvsB_jet_1;   //!
TBranch        *b_MVA_cdiscCvsB_jet_0;   //!
TBranch        *b_MVA_cdiscCvsL_jet_1;   //!
TBranch        *b_MVA_cdiscCvsL_jet_0;   //!
TBranch        *b_MVA_dRZc;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_mlb;   //!
TBranch        *b_MVA_FCNCtop_M;   //!
TBranch        *b_MVA_nJets_CharmL;   //!
TBranch        *b_MVA_NJets_CSVv2M;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_MVA_x1;   //!
TBranch        *b_MVA_x2;   //!
TBranch        *b_MVA_id1;   //!
TBranch        *b_MVA_id2;   //!
TBranch        *b_MVA_q;   //!
TBranch        *b_MVA_weight;   //!
TBranch        *b_MVA_channel;   //!
TBranch        *b_MVA_BDT;   //!
TBranch        *b_MVA_EqLumi;   //!
TBranch        *b_MVA_weight_puSF_up;   //!
TBranch        *b_MVA_weight_puSF_down;   //!
TBranch        *b_MVA_weight_electronSF_up;   //!
TBranch        *b_MVA_weight_electronSF_down;   //!
TBranch        *b_MVA_weight_muonSF_up;   //!
TBranch        *b_MVA_weight_muonSF_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_down;   //!
TBranch        *b_MVA_weight_btagSF_hf_up;   //!
TBranch        *b_MVA_weight_btagSF_hf_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_down;   //!
TBranch        *b_MVA_weight_btagSF_lf_up;   //!
TBranch        *b_MVA_weight_btagSF_lf_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_down;   //!



int main(int argc, char* argv[]){
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  //setTDRStyle(); // TO FIX
  //setMyStyle(); // TO FIX stat box title
  
  
  int testnr = -1;
  //////////// Settings of the analysis //////////////////
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      
      return 0;
    }
  }

