
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

#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/LHAPDF.h"
#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/Reweighting.h"

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH1F*> histo1DBDTvars;
map<string,TH1F*> histo1DPDF;
map<string,TH1F*> histo1DPDFMTW;
map<string,TH1F*> histo1DSys;
map<string,TH1F*> histo1DSysMTW;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotMTW;
map<string,MultiSamplePlot*> MSPlotTTZ;

double scaleFakes_uuu = 1.;
double scaleFakes_uue = 1.;
double scaleFakes_eeu = 1.;
double scaleFakes_eee = 1.;
bool FakeMuInW = false;
bool FakeEinW = false;

double nboffakesinWZ = 0.;
double nboffakesinTTSR = 0.;
double nboffakesinTTSRwith2jets = 0.;
double nboffakesinTTSRwith3jets = 0.;
double nboffakesinWZwith2jets = 0.;
double nboffakesinWZwith3jets = 0.;
////////////////////////////////// mapping ///////////////////////////////
map<string,TFile*> tFileMap;
TFile *fin;
map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;
vector < Dataset* > datasets;
vector < Dataset* > datasetsbf;
////////////////////////////////// functions ////////////////////////////////////////////
// bookkeeping
std::vector<std::string> split(const std::string &text,  char sep) ;
string ConvertIntToString(Int_t Number, Int_t pad);
string MakeTimeStamp();
string intToStr (Int_t number);
Double_t maximumValue(vector<double> array);
Double_t minimumValue(vector<double> array);
// initialisations

void InitMSPlotsMTW(string prefix, vector <int> decayChannels);
void InitMSPlotsTTZ(string prefix, vector <int> decayChannels);
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair, bool isZut);
void InitCalculatePDFWeightHisto(string dataSetName,bool doMTWtemplate,bool doTTZtemplate);
void InitMTWShapeHisto(string dataSetName, string systematic, Int_t isys,  vector <int> decayChannels);
void InitSystematicHisto(string dataSetName, string systematic, Int_t isys, bool doMTWtemplate);
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut);
void InitAnalyzerTree(TTree* tree);
void Init1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels);
// functions
vector<double> BDTCUT(string region, string coupling);
void CalculatePDFWeight(string dataSetName, Double_t BDT, Double_t MVA_weight_nom, Int_t MVA_channel,bool doMTWtemplate,bool doTTZtemplate);
void FillMTWPlots(Int_t d, string postfix, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm);
void FillTTZPlots(Int_t d, string postfix, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm);

void FillGeneralPlots(Int_t d, string prefix, vector <int> decayChannels, bool isZut , bool istoppair, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm,bool postfit);
void GetPDFEnvelope(string dataSetName,bool doMTWtemplate,bool doTTZtemplate);
void Fill1DHisto(string dataSetName,string systematic, bool istoppair, bool isZut, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel);
void FillMTWShapeHisto(string dataSetName, string systematic, Double_t weight_,Int_t isys, Int_t MVA_channel, vector <int> decayChannels);
void FillSystematicHisto(string dataSetName, string systematic, Double_t weight_, Int_t isys, bool doMTWtemplate);

//void Ini1DHistoBDTvariavles(bool istoppair, bool isZut, vector<int> decayChannels);


//////////////////////////////// settings ////////////////////////////////
bool makePlots = false;
bool doMTWtemplate = false;
bool doTTZtemplate = false;
bool singletoptemplate = false;
bool PlotSystematics = false;
bool PlotJeSystematics = false;
bool doPseudoData = false;
bool doSystematics = false;
Int_t channel = -999;
bool datafound = false;
bool testing = false;
bool doZut = false;
bool toppair = false;
bool addData = false;
bool DetermineCut = false;
bool CalculateSign = false;
bool doPDFunc  = false;
bool PlotMVAvars = false;
string placeNtup = "";
string tempstring = "";
string systematic = "";
string decaystring = "";
string coupling = "Zct";
string region = "singletop";
string placeOutputReading = "";
string combinetemplate_filename = "";
Double_t Luminosity = 36000.;
string dataSetName = "";
bool isData = false;
string tTreeName = "";
string postfix = "";
string output_histo_name = "";
string ntupleFileName ="";
Int_t nbin = 20; // 20 voor zct // 25 voor toppair zut en 20? voor st zut
//nbin = 20;
double BDT_begin = -1;
double BDT_end = 1;

bool doPostFit =false;
bool doErrorband = false;
Int_t nbinMTW = 20;
Double_t endMTW = 300.;

Int_t nbinTTZ = 1;
Double_t endTTZ= 140.;
Double_t beginTTZ= 50.;
Int_t nEntries = -1;
Int_t scaleNP = 10;

//////////////////////////// branches //////////////////////////////
// Declaration of leaf types
// Declaration of leaf types
Double_t        MVA_BDT;
Double_t        MVA_DeltaR_lep0Jet;
Double_t        MVA_DeltaR_lep1Jet;
Double_t        MVA_DeltaR_lep2Jet;
Double_t        MVA_DeltaR_MinLepJet;
Double_t        MVA_DeltaR_NonIsoLepJet;
Double_t        MVA_DeltaR_MinLepNonIsoJet;
Double_t        MVA_x1;
Double_t        MVA_x2;
Int_t           MVA_id1;
Int_t           MVA_id2;
Double_t        MVA_q;
Double_t        MVA_weight0;
Double_t        MVA_weight1;
Double_t        MVA_weight2;
Double_t        MVA_weight3;
Double_t        MVA_weight4;
Double_t        MVA_weight5;
Double_t        MVA_weight6;
Double_t        MVA_weight7;
Double_t        MVA_weight8;
Double_t        MVA_hdamp_up;
Double_t        MVA_hdamp_down;
Int_t           MVA_channel;
Float_t         MVA_weight;
Double_t        MVA_weight_nom;
Double_t        MVA_weight_puSF_up;
Double_t        MVA_weight_puSF_down;
Double_t        MVA_weight_electronSF_up;
Double_t        MVA_weight_electronSF_down;
Double_t        MVA_weight_puSF;
Double_t        MVA_weight_muonSF;
Double_t        MVA_weight_electronSF;
Double_t        MVA_weight_btagSF;
Double_t        MVA_weight_muonSF_up;
Double_t        MVA_weight_muonSF_down;
Double_t        MVA_weight_btagSF_cferr1_up;
Double_t        MVA_weight_btagSF_cferr1_down;
Double_t        MVA_weight_btagSF_cferr2_up;
Double_t        MVA_weight_btagSF_cferr2_down;
Double_t        MVA_weight_btagSF_hf_up;
Double_t        MVA_weight_btagSF_hf_down;
Double_t        MVA_weight_btagSF_hfstats1_up;
Double_t        MVA_weight_btagSF_hfstats1_down;
Double_t        MVA_weight_btagSF_hfstats2_up;
Double_t        MVA_weight_btagSF_hfstats2_down;
Double_t        MVA_weight_btagSF_lf_up;
Double_t        MVA_weight_btagSF_lf_down;
Double_t        MVA_weight_btagSF_lfstats1_up;
Double_t        MVA_weight_btagSF_lfstats1_down;
Double_t        MVA_weight_btagSF_lfstats2_up;
Double_t        MVA_weight_btagSF_lfstats2_down;
Double_t        MVA_weight_nloSF;
Float_t         MVA_region;
Double_t        MVA_EqLumi;
Double_t        MVA_Luminosity;
Float_t         MVA_lepton0_pt;
Float_t         MVA_lepton1_pt;
Float_t         MVA_lepton2_pt;
Float_t         MVA_Wlep_pt;
Float_t         MVA_Wboson_pt;
Float_t         MVA_SMbjet_pt;
Float_t         MVA_SMtop_pt;
Float_t         MVA_Zboson_pt;
Float_t         MVA_LightJet_pt;
Float_t         MVA_FCNCtop_pt;
Float_t         MVA_lepton0_eta;
Float_t         MVA_lepton1_eta;
Float_t         MVA_lepton2_eta;
Float_t         MVA_Wlep_eta;
Float_t         MVA_Wboson_eta;
Float_t         MVA_SMbjet_eta;
Float_t         MVA_SMtop_eta;
Float_t         MVA_Zboson_eta;
Float_t         MVA_LightJet_eta;
Float_t         MVA_FCNCtop_eta;
Float_t         MVA_lepton0_phi;
Float_t         MVA_lepton1_phi;
Float_t         MVA_lepton2_phi;
Float_t         MVA_Wlep_phi;
Float_t         MVA_Wboson_phi;
Float_t         MVA_SMbjet_phi;
Float_t         MVA_SMtop_phi;
Float_t         MVA_Zboson_phi;
Float_t         MVA_LightJet_phi;
Float_t         MVA_FCNCtop_phi;
Int_t           MVA_nElectrons;
Float_t         MVA_nJets;
Float_t         MVA_NJets_CSVv2L;
Float_t         MVA_NJets_CSVv2M;
Float_t         MVA_NJets_CSVv2T;
Int_t           MVA_nMuons;
Float_t         MVA_met;
Float_t         MVA_mWt;
Float_t         MVA_mWt2;
Float_t         MVA_SMtop_M;
Float_t         MVA_mlb;
Float_t         MVA_Wboson_M;
Float_t         MVA_SMtop_rap;
Float_t         MVA_dRWlepb;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_Wlep_Charge;
Float_t         MVA_charge_asym;
Float_t         MVA_TotalPt;
Float_t         MVA_TotalHt;
Float_t         MVA_TotalInvMass;
Float_t         MVA_TotalPt_jet;
Float_t         MVA_TotalHt_jet;
Float_t         MVA_TotalInvMass_jet;
Float_t         MVA_TotalPt_lep;
Float_t         MVA_TotalHt_lep;
Float_t         MVA_TotalInvMass_lep;
Float_t         MVA_bdiscCSVv2_jet_0;
Float_t         MVA_bdiscCSVv2_jet_1;
Float_t         MVA_CosTheta;
Float_t         MVA_CosTheta_alt;
Float_t         MVA_FCNCtop_rap;
Float_t         MVA_cdiscCvsB_cjet;
Float_t         MVA_cdiscCvsL_cjet;
Float_t         MVA_cdiscCvsB_highestdisjet;
Float_t         MVA_cdiscCvsL_highestdisjet;
Float_t         MVA_FCNCtop_M;
Float_t         MVA_Zboson_M;
Float_t         MVA_Bdis_Lightjet;
Float_t         MVA_Bdis_OtherJets;
Float_t         MVA_dRZc;
Float_t         MVA_dPhiZc;
Float_t         MVA_cdiscCvsB_jet_1;
Float_t         MVA_cdiscCvsL_jet_1;
Float_t         MVA_cdiscCvsB_jet_0;
Float_t         MVA_cdiscCvsL_jet_0;
Float_t         MVA_nJets_CharmL;
Float_t         MVA_nJets_CharmM;
Float_t         MVA_nJets_CharmT;
Float_t         MVA_dRSMFCNCtop;
Float_t         MVA_ptWQ;
Float_t         MVA_deltaRjj_max;
Float_t         MVA_deltaRjj_min;
Float_t         MVA_deltaRjj_sum;
Float_t         MVA_deltaRWlepJet_max;
Float_t         MVA_deltaRWlepJet_min;
Float_t         MVA_dRSMjetLightjet;
Float_t         MVA_dRZb;
Float_t         MVA_dRWlepc;
Float_t         MVA_dRZWlep;
Float_t         MVA_dRZSMtop;
Float_t         MVA_dPhiSMFCNCtop;
Float_t         MVA_dPhiZb;
Float_t         MVA_dPhiWlepc;
Float_t         MVA_dPhiZWlep;
Float_t         MVA_dPhiZSMtop;
Float_t         MVA_m3l;
// List of branches

TBranch       *b_MVA_DeltaR_lep0Jet;
TBranch       *b_MVA_DeltaR_lep1Jet;
TBranch       *b_MVA_DeltaR_lep2Jet;
TBranch       *b_MVA_DeltaR_MinLepJet;
TBranch       *b_MVA_DeltaR_NonIsoLepJet;
TBranch       *b_MVA_DeltaR_MinLepNonIsoJet;

TBranch         *b_MVA_BDT;
TBranch        *b_MVA_x1;   //!
TBranch        *b_MVA_x2;   //!
TBranch        *b_MVA_id1;   //!
TBranch        *b_MVA_id2;   //!
TBranch        *b_MVA_q;   //!
TBranch        *b_MVA_weight0;   //!
TBranch        *b_MVA_weight1;   //!
TBranch        *b_MVA_weight2;   //!
TBranch        *b_MVA_weight3;   //!
TBranch        *b_MVA_weight4;   //!
TBranch        *b_MVA_weight5;   //!
TBranch        *b_MVA_weight6;   //!
TBranch        *b_MVA_weight7;   //!
TBranch        *b_MVA_weight8;   //!
TBranch        *b_MVA_hdamp_up;   //!
TBranch        *b_MVA_hdamp_down;   //!
TBranch        *b_MVA_channel;   //!
TBranch        *b_MVA_weight;   //!
TBranch        *b_MVA_weight_nom;   //!
TBranch        *b_MVA_weight_puSF_up;   //!
TBranch        *b_MVA_weight_puSF_down;   //!
TBranch        *b_MVA_weight_electronSF_up;   //!
TBranch        *b_MVA_weight_electronSF_down;   //!
TBranch        *b_MVA_weight_puSF;   //!
TBranch        *b_MVA_weight_muonSF;   //!
TBranch        *b_MVA_weight_electronSF;   //!
TBranch        *b_MVA_weight_btagSF;   //!
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
TBranch        *b_MVA_weight_nloSF;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_MVA_EqLumi;   //!
TBranch        *b_MVA_Luminosity;   //!
TBranch        *b_MVA_lepton0_pt;   //!
TBranch        *b_MVA_lepton1_pt;   //!
TBranch        *b_MVA_lepton2_pt;   //!
TBranch        *b_MVA_Wlep_pt;   //!
TBranch        *b_MVA_Wboson_pt;   //!
TBranch        *b_MVA_SMbjet_pt;   //!
TBranch        *b_MVA_SMtop_pt;   //!
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_LightJet_pt;   //!
TBranch        *b_MVA_FCNCtop_pt;   //!
TBranch        *b_MVA_lepton0_eta;   //!
TBranch        *b_MVA_lepton1_eta;   //!
TBranch        *b_MVA_lepton2_eta;   //!
TBranch        *b_MVA_Wlep_eta;   //!
TBranch        *b_MVA_Wboson_eta;   //!
TBranch        *b_MVA_SMbjet_eta;   //!
TBranch        *b_MVA_SMtop_eta;   //!
TBranch        *b_MVA_Zboson_eta;   //!
TBranch        *b_MVA_LightJet_eta;   //!
TBranch        *b_MVA_FCNCtop_eta;   //!
TBranch        *b_MVA_lepton0_phi;   //!
TBranch        *b_MVA_lepton1_phi;   //!
TBranch        *b_MVA_lepton2_phi;   //!
TBranch        *b_MVA_Wlep_phi;   //!
TBranch        *b_MVA_Wboson_phi;   //!
TBranch        *b_MVA_SMbjet_phi;   //!
TBranch        *b_MVA_SMtop_phi;   //!
TBranch        *b_MVA_Zboson_phi;   //!
TBranch        *b_MVA_LightJet_phi;   //!
TBranch        *b_MVA_FCNCtop_phi;   //!
TBranch        *b_MVA_nElectrons;   //!
TBranch        *b_MVA_nJets;   //!
TBranch        *b_MVA_NJets_CSVv2L;   //!
TBranch        *b_MVA_NJets_CSVv2M;   //!
TBranch        *b_MVA_NJets_CSVv2T;   //!
TBranch        *b_MVA_nMuons;   //!
TBranch        *b_MVA_met;   //!
TBranch        *b_MVA_mWt;   //!
TBranch        *b_MVA_mWt2;   //!
TBranch        *b_MVA_SMtop_M;   //!
TBranch        *b_MVA_mlb;   //!
TBranch        *b_MVA_Wboson_M;   //!
TBranch        *b_MVA_SMtop_rap;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_Wlep_Charge;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_TotalPt;   //!
TBranch        *b_MVA_TotalHt;   //!
TBranch        *b_MVA_TotalInvMass;   //!
TBranch        *b_MVA_TotalPt_jet;   //!
TBranch        *b_MVA_TotalHt_jet;   //!
TBranch        *b_MVA_TotalInvMass_jet;   //!
TBranch        *b_MVA_TotalPt_lep;   //!
TBranch        *b_MVA_TotalHt_lep;   //!
TBranch        *b_MVA_TotalInvMass_lep;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_0;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_1;   //!
TBranch        *b_MVA_CosTheta;   //!
TBranch        *b_MVA_CosTheta_alt;   //!
TBranch        *b_MVA_FCNCtop_rap;   //!
TBranch        *b_MVA_cdiscCvsB_cjet;   //!
TBranch        *b_MVA_cdiscCvsL_cjet;   //!
TBranch        *b_MVA_cdiscCvsB_highestdisjet;   //!
TBranch        *b_MVA_cdiscCvsL_highestdisjet;   //!
TBranch        *b_MVA_FCNCtop_M;   //!
TBranch        *b_MVA_Zboson_M;   //!
TBranch        *b_MVA_Bdis_Lightjet;   //!
TBranch        *b_MVA_Bdis_OtherJets;   //!
TBranch        *b_MVA_dRZc;   //!
TBranch        *b_MVA_dPhiZc;   //!
TBranch        *b_MVA_cdiscCvsB_jet_1;   //!
TBranch        *b_MVA_cdiscCvsL_jet_1;   //!
TBranch        *b_MVA_cdiscCvsB_jet_0;   //!
TBranch        *b_MVA_cdiscCvsL_jet_0;   //!
TBranch        *b_MVA_nJets_CharmL;   //!
TBranch        *b_MVA_nJets_CharmM;   //!
TBranch        *b_MVA_nJets_CharmT;   //!
TBranch        *b_MVA_dRSMFCNCtop;   //!
TBranch        *b_MVA_ptWQ;   //!Branch
TBranch        *b_MVA_deltaRjj_max;   //!
TBranch        *b_MVA_deltaRjj_min;   //!
TBranch        *b_MVA_deltaRjj_sum;   //!
TBranch        *b_MVA_deltaRWlepJet_max;   //!
TBranch        *b_MVA_deltaRWlepJet_min;   //!
TBranch        *b_MVA_dRSMjetLightjet;   //!
TBranch        *b_MVA_dRZb;   //!
TBranch        *b_MVA_dRWlepc;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_dRZSMtop;   //!
TBranch        *b_MVA_dPhiSMFCNCtop;   //!
TBranch        *b_MVA_dPhiZb;   //!
TBranch        *b_MVA_dPhiWlepc;   //!
TBranch        *b_MVA_dPhiZWlep;   //!
TBranch        *b_MVA_dPhiZSMtop;   //!
TBranch        *b_MVA_m3l;   //!

Double_t MVA_weightWZcorr;
TBranch *b_MVA_weightWZcorr;



TH1F*  hist_BDT_tt_nonpromptinZ = new TH1F("hist_BDT_tt_nonpromptinZ","hist_BDT_tt_nonpromptinZ;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_nonpromptinW = new TH1F("hist_BDT_tt_nonpromptinW","hist_BDT_tt_nonpromptinW;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_uuu_nonpromptinZ = new TH1F("hist_BDT_tt_uuu_nonpromptinZ","hist_BDT_tt_uuu_nonpromptinZ;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_uuu_nonpromptinW = new TH1F("hist_BDT_tt_uuu_nonpromptinW","hist_BDT_tt_uuu_nonpromptinW;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_eee_nonpromptinZ = new TH1F("hist_BDT_tt_eee_nonpromptinZ","hist_BDT_tt_eee_nonpromptinZ;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_eee_nonpromptinW = new TH1F("hist_BDT_tt_eee_nonpromptinW","hist_BDT_tt_eee_nonpromptinW;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_eeu_nonpromptinZ = new TH1F("hist_BDT_tt_eeu_nonpromptinZ","hist_BDT_tt_eeu_nonpromptinZ;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_eeu_nonpromptinW = new TH1F("hist_BDT_tt_eeu_nonpromptinW","hist_BDT_tt_eeu_nonpromptinW;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_uue_nonpromptinZ = new TH1F("hist_BDT_tt_uue_nonpromptinZ","hist_BDT_tt_uue_nonpromptinZ;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_tt_uue_nonpromptinW = new TH1F("hist_BDT_tt_uue_nonpromptinW","hist_BDT_tt_uue_nonpromptinW;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_WZ_light = new TH1F("hist_BDT_WZ_light","hist_BDT_WZ_light;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_c = new TH1F("hist_BDT_WZ_c","hist_BDT_WZ_c;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uuu_light = new TH1F("hist_BDT_WZ_uuu_light","hist_BDT_WZ_uuu_light;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uuu_c = new TH1F("hist_BDT_WZ_uuu_c","hist_BDT_WZ_uuu_c;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eee_light = new TH1F("hist_BDT_WZ_eee_light","hist_BDT_WZ_eee_light;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eee_c = new TH1F("hist_BDT_WZ_eee_c","hist_BDT_WZ_eee_c;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eeu_light = new TH1F("hist_BDT_WZ_eeu_light","hist_BDT_WZ_eeu_light;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eeu_c = new TH1F("hist_BDT_WZ_eeu_c","hist_BDT_WZ_eeu_c;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uue_light = new TH1F("hist_BDT_WZ_uue_light","hist_BDT_WZ_uue_light;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uue_c = new TH1F("hist_BDT_WZ_uue_c","hist_BDT_WZ_uue_c;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_WZ_b = new TH1F("hist_BDT_WZ_b","hist_BDT_WZ_b;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uuu_b = new TH1F("hist_BDT_WZ_uuu_b","hist_BDT_WZ_uuu_b;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_uue_b = new TH1F("hist_BDT_WZ_uue_b","hist_BDT_WZ_uue_b;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eeu_b = new TH1F("hist_BDT_WZ_eeu_b","hist_BDT_WZ_eeu_b;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_WZ_eee_b = new TH1F("hist_BDT_WZ_eee_b","hist_BDT_WZ_eee_b;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);

map<string,TH1F*> histo1DMTW;
TH1F*  hist_BDT_JES_nom_sig = new TH1F("hist_BDT_JES_nom_sig","Effect of JES systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JES_nom_bkg = new TH1F("hist_BDT_JES_nom_bkg","Effect of JES systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JES_up_sig = new TH1F("hist_BDT_JES_up_sig","Effect of JES systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JES_up_bkg = new TH1F("hist_BDT_JES_up_bkg","Effect of JES systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JES_down_sig = new TH1F("hist_BDT_JES_down_sig","Effect of JES systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JES_down_bkg = new TH1F("hist_BDT_JES_down_bkg","Effect of JES systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_JER_nom_sig = new TH1F("hist_BDT_JER_nom_sig","Effect of JER systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JER_nom_bkg = new TH1F("hist_BDT_JER_nom_bkg","Effect of JER systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JER_up_sig = new TH1F("hist_BDT_JER_up_sig","Effect of JER systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JER_up_bkg = new TH1F("hist_BDT_JER_up_bkg","Effect of JER systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JER_down_sig = new TH1F("hist_BDT_JER_down_sig","Effect of JER systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_JER_down_bkg = new TH1F("hist_BDT_JER_down_bkg","Effect of JER systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);



TH1F*  hist_BDT_puSF_nom_sig = new TH1F("hist_BDT_puSF_nom_sig","Effect of pile up systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_puSF_nom_bkg = new TH1F("hist_BDT_puSF_nom_bkg","Effect of pile up systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_puSF_up_sig = new TH1F("hist_BDT_puSF_up_sig","Effect of pile up systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_puSF_up_bkg = new TH1F("hist_BDT_puSF_up_bkg","Effect of pile up systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_puSF_down_sig = new TH1F("hist_BDT_puSF_down_sig","Effect of pile up systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_puSF_down_bkg = new TH1F("hist_BDT_puSF_down_bkg","Effect of pile up systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_electronSF_nom_sig = new TH1F("hist_BDT_electronSF_nom_sig","Effect of electron SF systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_electronSF_nom_bkg = new TH1F("hist_BDT_electronSF_nom_bkg","Effect of electron SF systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_electronSF_up_sig = new TH1F("hist_BDT_electronSF_up_sig","Effect of electron SF systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_electronSF_up_bkg = new TH1F("hist_BDT_electronSF_up_bkg","Effect of electron SF systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_electronSF_down_sig = new TH1F("hist_BDT_electronSF_down_sig","Effect of electron SF systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_electronSF_down_bkg = new TH1F("hist_BDT_electronSF_down_bkg","Effect of electron SF systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_BDT_muonSF_nom_sig = new TH1F("hist_BDT_muonSF_nom_sig","Effect of muon SF systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_muonSF_nom_bkg = new TH1F("hist_BDT_muonSF_nom_bkg","Effect of muon SF systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_muonSF_up_sig = new TH1F("hist_BDT_muonSF_up_sig","Effect of muon SF systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_muonSF_up_bkg = new TH1F("hist_BDT_muonSF_up_bkg","Effect of muon SF systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_muonSF_down_sig = new TH1F("hist_BDT_muonSF_down_sig","Effect of muon SF systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_muonSF_down_bkg = new TH1F("hist_BDT_muonSF_down_bkg","Effect of muon SF systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_btagSF_cferr1_nom_sig = new TH1F("hist_BDT_btagSF_cferr1_nom_sig","Effect of btag SF cferr1 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr1_nom_bkg = new TH1F("hist_BDT_btagSF_cferr1_nom_bkg","Effect of btag SF cferr1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr1_up_sig = new TH1F("hist_BDT_btagSF_cferr1_up_sig","Effect of btag SF cferr1 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr1_up_bkg = new TH1F("hist_BDT_btagSF_cferr1_up_bkg","Effect of btag SF cferr1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr1_down_sig = new TH1F("hist_BDT_btagSF_cferr1_down_sig","Effect of btag SF cferr1 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr1_down_bkg = new TH1F("hist_BDT_btagSF_cferr1_down_bkg","Effect of btag SF cferr1 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_btagSF_cferr2_nom_sig = new TH1F("hist_BDT_btagSF_cferr2_nom_sig","Effect of btag SF cferr2 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr2_nom_bkg = new TH1F("hist_BDT_btagSF_cferr2_nom_bkg","Effect of btag SF cferr2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr2_up_sig = new TH1F("hist_BDT_btagSF_cferr2_up_sig","Effect of btag SF cferr2 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr2_up_bkg = new TH1F("hist_BDT_btagSF_cferr2_up_bkg","Effect of btag SF cferr2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr2_down_sig = new TH1F("hist_BDT_btagSF_cferr2_down_sig","Effect of btag SF cferr2 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_cferr2_down_bkg = new TH1F("hist_BDT_btagSF_cferr2_down_bkg","Effect of btag SF cferr2 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_BDT_btagSF_hf_nom_sig = new TH1F("hist_BDT_btagSF_hf_nom_sig","Effect of btag SF hf systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hf_nom_bkg = new TH1F("hist_BDT_btagSF_hf_nom_bkg","Effect of btag SF hf systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hf_up_sig = new TH1F("hist_BDT_btagSF_hf_up_sig","Effect of btag SF hf systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hf_up_bkg = new TH1F("hist_BDT_btagSF_hf_up_bkg","Effect of btag SF hf systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hf_down_sig = new TH1F("hist_BDT_btagSF_hf_down_sig","Effect of btag SF hf systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hf_down_bkg = new TH1F("hist_BDT_btagSF_hf_down_bkg","Effect of btag SF hf systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_btagSF_hfstats1_nom_sig = new TH1F("hist_BDT_btagSF_hfstats1_nom_sig","Effect of btag SF hfstats1 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats1_nom_bkg = new TH1F("hist_BDT_btagSF_hfstats1_nom_bkg","Effect of btag SF hfstats1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats1_up_sig = new TH1F("hist_BDT_btagSF_hfstats1_up_sig","Effect of btag SF hfstats1 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats1_up_bkg = new TH1F("hist_BDT_btagSF_hfstats1_up_bkg","Effect of btag SF hfstats1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats1_down_sig = new TH1F("hist_BDT_btagSF_hfstats1_down_sig","Effect of btag SF hfstats1 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats1_down_bkg = new TH1F("hist_BDT_btagSF_hfstats1_down_bkg","Effect of btag SF hfstats1 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_BDT_btagSF_hfstats2_nom_sig = new TH1F("hist_BDT_btagSF_hfstats2_nom_sig","Effect of btag SF hfstats2 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats2_nom_bkg = new TH1F("hist_BDT_btagSF_hfstats2_nom_bkg","Effect of btag SF hfstats2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats2_up_sig = new TH1F("hist_BDT_btagSF_hfstats2_up_sig","Effect of btag SF hfstats2 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats2_up_bkg = new TH1F("hist_BDT_btagSF_hfstats2_up_bkg","Effect of btag SF hfstats2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats2_down_sig = new TH1F("hist_BDT_btagSF_hfstats2_down_sig","Effect of btag SF hfstats2 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_hfstats2_down_bkg = new TH1F("hist_BDT_btagSF_hfstats2_down_bkg","Effect of btag SF hfstats2 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_BDT_btagSF_lf_nom_sig = new TH1F("hist_BDT_btagSF_lf_nom_sig","Effect of btag SF lf systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lf_nom_bkg = new TH1F("hist_BDT_btagSF_lf_nom_bkg","Effect of btag SF lf systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lf_up_sig = new TH1F("hist_BDT_btagSF_lf_up_sig","Effect of btag SF lf systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lf_up_bkg = new TH1F("hist_BDT_btagSF_lf_up_bkg","Effect of btag SF lf systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lf_down_sig = new TH1F("hist_BDT_btagSF_lf_down_sig","Effect of btag SF lf systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lf_down_bkg = new TH1F("hist_BDT_btagSF_lf_down_bkg","Effect of btag SF lf systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);

TH1F*  hist_BDT_btagSF_lfstats1_nom_sig = new TH1F("hist_BDT_btagSF_lfstats1_nom_sig","Effect of btag SF lfstats1 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats1_nom_bkg = new TH1F("hist_BDT_btagSF_lfstats1_nom_bkg","Effect of btag SF lfstats1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats1_up_sig = new TH1F("hist_BDT_btagSF_lfstats1_up_sig","Effect of btag SF lfstats1 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats1_up_bkg = new TH1F("hist_BDT_btagSF_lfstats1_up_bkg","Effect of btag SF lfstats1 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats1_down_sig = new TH1F("hist_BDT_btagSF_lfstats1_down_sig","Effect of btag SF lfstats1 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats1_down_bkg = new TH1F("hist_BDT_btagSF_lfstats1_down_bkg","Effect of btag SF lfstats1 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_BDT_btagSF_lfstats2_nom_sig = new TH1F("hist_BDT_btagSF_lfstats2_nom_sig","Effect of btag SF lfstats2 systematics on the BDT: Signal;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats2_nom_bkg = new TH1F("hist_BDT_btagSF_lfstats2_nom_bkg","Effect of btag SF lfstats2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats2_up_sig = new TH1F("hist_BDT_btagSF_lfstats2_up_sig","Effect of btag SF lfstats2 systematics on the BDT: Signal:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats2_up_bkg = new TH1F("hist_BDT_btagSF_lfstats2_up_bkg","Effect of btag SF lfstats2 systematics on the BDT: Background;BDT;Nb. of evts", nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats2_down_sig = new TH1F("hist_BDT_btagSF_lfstats2_down_sig","Effect of btag SF lfstats2 systematics on the BDT: Background:BDT:Nb. of evts" ,nbin,BDT_begin,BDT_end);
TH1F*  hist_BDT_btagSF_lfstats2_down_bkg = new TH1F("hist_BDT_btagSF_lfstats2_down_bkg","Effect of btag SF lfstats2 systematics on the BDT: Background:BDT:Nb. of evts", nbin,BDT_begin,BDT_end);


TH1F*  hist_mWt_JES_nom_sig = new TH1F("hist_mWt_JES_nom_sig","Effect of JES systematics on the mWt: Signal;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JES_nom_bkg = new TH1F("hist_mWt_JES_nom_bkg","Effect of JES systematics on the mWt: Background;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JES_up_sig = new TH1F("hist_mWt_JES_up_sig","Effect of JES systematics on the mWt: Signal:mWt:Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_JES_up_bkg = new TH1F("hist_mWt_JES_up_bkg","Effect of JES systematics on the mWt: Background;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JES_down_sig = new TH1F("hist_mWt_JES_down_sig","Effect of JES systematics on the mWt: Background:mWt:Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_JES_down_bkg = new TH1F("hist_mWt_JES_down_bkg","Effect of JES systematics on the mWt: Background:mWt:Nb. of evts", nbinMTW,0, endMTW);

TH1F*  hist_mWt_JER_nom_sig = new TH1F("hist_mWt_JER_nom_sig","Effect of JER systematics on the mWt: Signal;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JER_nom_bkg = new TH1F("hist_mWt_JER_nom_bkg","Effect of JER systematics on the mWt: Background;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JER_up_sig = new TH1F("hist_mWt_JER_up_sig","Effect of JER systematics on the mWt: Signal:mWt:Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_JER_up_bkg = new TH1F("hist_mWt_JER_up_bkg","Effect of JER systematics on the mWt: Background;mWt;Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_JER_down_sig = new TH1F("hist_mWt_JER_down_sig","Effect of JER systematics on the mWt: Background:mWt:Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_JER_down_bkg = new TH1F("hist_mWt_JER_down_bkg","Effect of JER systematics on the mWt: Background:mWt:Nb. of evts", nbinMTW,0, endMTW);




TH1F*  hist_mWt_puSF_nom_sig = new TH1F("hist_mWt_puSF_nom_sig","Effect of pile up systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_puSF_nom_bkg = new TH1F("hist_mWt_puSF_nom_bkg","Effect of pile up systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_puSF_up_sig = new TH1F("hist_mWt_puSF_up_sig","Effect of pile up systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_puSF_up_bkg = new TH1F("hist_mWt_puSF_up_bkg","Effect of pile up systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_puSF_down_sig = new TH1F("hist_mWt_puSF_down_sig","Effect of pile up systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_puSF_down_bkg = new TH1F("hist_mWt_puSF_down_bkg","Effect of pile up systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);

TH1F*  hist_mWt_electronSF_nom_sig = new TH1F("hist_mWt_electronSF_nom_sig","Effect of electron SF systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_electronSF_nom_bkg = new TH1F("hist_mWt_electronSF_nom_bkg","Effect of electron SF systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_electronSF_up_sig = new TH1F("hist_mWt_electronSF_up_sig","Effect of electron SF systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_electronSF_up_bkg = new TH1F("hist_mWt_electronSF_up_bkg","Effect of electron SF systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_electronSF_down_sig = new TH1F("hist_mWt_electronSF_down_sig","Effect of electron SF systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_electronSF_down_bkg = new TH1F("hist_mWt_electronSF_down_bkg","Effect of electron SF systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_muonSF_nom_sig = new TH1F("hist_mWt_muonSF_nom_sig","Effect of muon SF systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_muonSF_nom_bkg = new TH1F("hist_mWt_muonSF_nom_bkg","Effect of muon SF systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_muonSF_up_sig = new TH1F("hist_mWt_muonSF_up_sig","Effect of muon SF systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_muonSF_up_bkg = new TH1F("hist_mWt_muonSF_up_bkg","Effect of muon SF systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_muonSF_down_sig = new TH1F("hist_mWt_muonSF_down_sig","Effect of muon SF systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_muonSF_down_bkg = new TH1F("hist_mWt_muonSF_down_bkg","Effect of muon SF systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_btagSF_cferr1_nom_sig = new TH1F("hist_mWt_btagSF_cferr1_nom_sig","Effect of btag SF cferr1 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr1_nom_bkg = new TH1F("hist_mWt_btagSF_cferr1_nom_bkg","Effect of btag SF cferr1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr1_up_sig = new TH1F("hist_mWt_btagSF_cferr1_up_sig","Effect of btag SF cferr1 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr1_up_bkg = new TH1F("hist_mWt_btagSF_cferr1_up_bkg","Effect of btag SF cferr1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr1_down_sig = new TH1F("hist_mWt_btagSF_cferr1_down_sig","Effect of btag SF cferr1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr1_down_bkg = new TH1F("hist_mWt_btagSF_cferr1_down_bkg","Effect of btag SF cferr1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);

TH1F*  hist_mWt_btagSF_cferr2_nom_sig = new TH1F("hist_mWt_btagSF_cferr2_nom_sig","Effect of btag SF cferr2 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr2_nom_bkg = new TH1F("hist_mWt_btagSF_cferr2_nom_bkg","Effect of btag SF cferr2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr2_up_sig = new TH1F("hist_mWt_btagSF_cferr2_up_sig","Effect of btag SF cferr2 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr2_up_bkg = new TH1F("hist_mWt_btagSF_cferr2_up_bkg","Effect of btag SF cferr2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr2_down_sig = new TH1F("hist_mWt_btagSF_cferr2_down_sig","Effect of btag SF cferr2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_cferr2_down_bkg = new TH1F("hist_mWt_btagSF_cferr2_down_bkg","Effect of btag SF cferr2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_btagSF_hf_nom_sig = new TH1F("hist_mWt_btagSF_hf_nom_sig","Effect of btag SF hf systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hf_nom_bkg = new TH1F("hist_mWt_btagSF_hf_nom_bkg","Effect of btag SF hf systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hf_up_sig = new TH1F("hist_mWt_btagSF_hf_up_sig","Effect of btag SF hf systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hf_up_bkg = new TH1F("hist_mWt_btagSF_hf_up_bkg","Effect of btag SF hf systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hf_down_sig = new TH1F("hist_mWt_btagSF_hf_down_sig","Effect of btag SF hf systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hf_down_bkg = new TH1F("hist_mWt_btagSF_hf_down_bkg","Effect of btag SF hf systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);

TH1F*  hist_mWt_btagSF_hfstats1_nom_sig = new TH1F("hist_mWt_btagSF_hfstats1_nom_sig","Effect of btag SF hfstats1 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats1_nom_bkg = new TH1F("hist_mWt_btagSF_hfstats1_nom_bkg","Effect of btag SF hfstats1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats1_up_sig = new TH1F("hist_mWt_btagSF_hfstats1_up_sig","Effect of btag SF hfstats1 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats1_up_bkg = new TH1F("hist_mWt_btagSF_hfstats1_up_bkg","Effect of btag SF hfstats1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats1_down_sig = new TH1F("hist_mWt_btagSF_hfstats1_down_sig","Effect of btag SF hfstats1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats1_down_bkg = new TH1F("hist_mWt_btagSF_hfstats1_down_bkg","Effect of btag SF hfstats1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_btagSF_hfstats2_nom_sig = new TH1F("hist_mWt_btagSF_hfstats2_nom_sig","Effect of btag SF hfstats2 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats2_nom_bkg = new TH1F("hist_mWt_btagSF_hfstats2_nom_bkg","Effect of btag SF hfstats2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats2_up_sig = new TH1F("hist_mWt_btagSF_hfstats2_up_sig","Effect of btag SF hfstats2 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats2_up_bkg = new TH1F("hist_mWt_btagSF_hfstats2_up_bkg","Effect of btag SF hfstats2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats2_down_sig = new TH1F("hist_mWt_btagSF_hfstats2_down_sig","Effect of btag SF hfstats2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_hfstats2_down_bkg = new TH1F("hist_mWt_btagSF_hfstats2_down_bkg","Effect of btag SF hfstats2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_btagSF_lf_nom_sig = new TH1F("hist_mWt_btagSF_lf_nom_sig","Effect of btag SF lf systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lf_nom_bkg = new TH1F("hist_mWt_btagSF_lf_nom_bkg","Effect of btag SF lf systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lf_up_sig = new TH1F("hist_mWt_btagSF_lf_up_sig","Effect of btag SF lf systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lf_up_bkg = new TH1F("hist_mWt_btagSF_lf_up_bkg","Effect of btag SF lf systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lf_down_sig = new TH1F("hist_mWt_btagSF_lf_down_sig","Effect of btag SF lf systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lf_down_bkg = new TH1F("hist_mWt_btagSF_lf_down_bkg","Effect of btag SF lf systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);

TH1F*  hist_mWt_btagSF_lfstats1_nom_sig = new TH1F("hist_mWt_btagSF_lfstats1_nom_sig","Effect of btag SF lfstats1 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats1_nom_bkg = new TH1F("hist_mWt_btagSF_lfstats1_nom_bkg","Effect of btag SF lfstats1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats1_up_sig = new TH1F("hist_mWt_btagSF_lfstats1_up_sig","Effect of btag SF lfstats1 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats1_up_bkg = new TH1F("hist_mWt_btagSF_lfstats1_up_bkg","Effect of btag SF lfstats1 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats1_down_sig = new TH1F("hist_mWt_btagSF_lfstats1_down_sig","Effect of btag SF lfstats1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats1_down_bkg = new TH1F("hist_mWt_btagSF_lfstats1_down_bkg","Effect of btag SF lfstats1 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);


TH1F*  hist_mWt_btagSF_lfstats2_nom_sig = new TH1F("hist_mWt_btagSF_lfstats2_nom_sig","Effect of btag SF lfstats2 systematics on the m_{T}(W): DD non prompt lepton;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats2_nom_bkg = new TH1F("hist_mWt_btagSF_lfstats2_nom_bkg","Effect of btag SF lfstats2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats2_up_sig = new TH1F("hist_mWt_btagSF_lfstats2_up_sig","Effect of btag SF lfstats2 systematics on the m_{T}(W): DD non prompt lepton:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats2_up_bkg = new TH1F("hist_mWt_btagSF_lfstats2_up_bkg","Effect of btag SF lfstats2 systematics on the m_{T}(W): Background;m_{T}(W);Nb. of evts", nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats2_down_sig = new TH1F("hist_mWt_btagSF_lfstats2_down_sig","Effect of btag SF lfstats2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts" ,nbinMTW,0, endMTW);
TH1F*  hist_mWt_btagSF_lfstats2_down_bkg = new TH1F("hist_mWt_btagSF_lfstats2_down_bkg","Effect of btag SF lfstats2 systematics on the m_{T}(W): Background:m_{T}(W):Nb. of evts", nbinMTW,0, endMTW);




/*
 ///////// CUTFLOWS ///////////////
 TH1F*  CutflowTableHisto__WZ = new TH1F("CutflowTableHisto__WZ", "CutflowTableHisto__WZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__WZ = new TH1F("CutflowTableHisto_eee__WZ", "CutflowTableHisto_eee__WZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__WZ = new TH1F("CutflowTableHisto_eeu__WZ", "CutflowTableHisto_eeu__WZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__WZ = new TH1F("CutflowTableHisto_uue__WZ", "CutflowTableHisto_uue__WZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__WZ = new TH1F("CutflowTableHisto_uuu__WZ", "CutflowTableHisto_uuu__WZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 TH1F*  CutflowTableHisto__fake = new TH1F("CutflowTableHisto__fake", "CutflowTableHisto__fake" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__fake = new TH1F("CutflowTableHisto_eee__fake", "CutflowTableHisto_eee__fake" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__fake = new TH1F("CutflowTableHisto_eeu__fake", "CutflowTableHisto_eeu__fake" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__fake = new TH1F("CutflowTableHisto_uue__fake", "CutflowTableHisto_uue__fake" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__fake = new TH1F("CutflowTableHisto_uuu__fake", "CutflowTableHisto_uuu__fake" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 TH1F*  CutflowTableHisto__data = new TH1F("CutflowTableHisto__data", "CutflowTableHisto__data" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__data = new TH1F("CutflowTableHisto_eee__data", "CutflowTableHisto_eee__data" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__data = new TH1F("CutflowTableHisto_eeu__data", "CutflowTableHisto_eeu__data" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__data = new TH1F("CutflowTableHisto_uue__data", "CutflowTableHisto_uue__data" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__data = new TH1F("CutflowTableHisto_uuu__data", "CutflowTableHisto_uuu__data" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 
 TH1F*  CutflowTableHisto__tZq = new TH1F("CutflowTableHisto__tZq", "CutflowTableHisto__tZq" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__tZq = new TH1F("CutflowTableHisto_eee__tZq", "CutflowTableHisto_eee__tZq" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__tZq = new TH1F("CutflowTableHisto_eeu__tZq", "CutflowTableHisto_eeu__tZq" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__tZq = new TH1F("CutflowTableHisto_uue__tZq", "CutflowTableHisto_uue__tZq" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__tZq = new TH1F("CutflowTableHisto_uuu__tZq", "CutflowTableHisto_uuu__tZq" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 TH1F*  CutflowTableHisto__ttZ = new TH1F("CutflowTableHisto__ttZ", "CutflowTableHisto__ttZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__ttZ = new TH1F("CutflowTableHisto_eee__ttZ", "CutflowTableHisto_eee__ttZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__ttZ = new TH1F("CutflowTableHisto_eeu__ttZ", "CutflowTableHisto_eeu__ttZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__ttZ = new TH1F("CutflowTableHisto_uue__ttZ", "CutflowTableHisto_uue__ttZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__ttZ = new TH1F("CutflowTableHisto_uuu__ttZ", "CutflowTableHisto_uuu__ttZ" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 TH1F*  CutflowTableHisto__FCNCZut = new TH1F("CutflowTableHisto__FCNCZut", "CutflowTableHisto__FCNCZut" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__FCNCZut = new TH1F("CutflowTableHisto_eee__FCNCZut", "CutflowTableHisto_eee__FCNCZut" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__FCNCZut = new TH1F("CutflowTableHisto_eeu__FCNCZut", "CutflowTableHisto_eeu__FCNCZut" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__FCNCZut = new TH1F("CutflowTableHisto_uue__FCNCZut", "CutflowTableHisto_uue__FCNCZut" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__FCNCZut = new TH1F("CutflowTableHisto_uuu__FCNCZut", "CutflowTableHisto_uuu__FCNCZut" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 
 TH1F*  CutflowTableHisto__FCNCZct = new TH1F("CutflowTableHisto__FCNCZct", "CutflowTableHisto__FCNCZct" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__FCNCZct = new TH1F("CutflowTableHisto_eee__FCNCZct", "CutflowTableHisto_eee__FCNCZct" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__FCNCZct = new TH1F("CutflowTableHisto_eeu__FCNCZct", "CutflowTableHisto_eeu__FCNCZct" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__FCNCZct = new TH1F("CutflowTableHisto_uue__FCNCZct", "CutflowTableHisto_uue__FCNCZct" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__FCNCZct = new TH1F("CutflowTableHisto_uuu__FCNCZct", "CutflowTableHisto_uuu__FCNCZct" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 
 
 TH1F*  CutflowTableHisto__other = new TH1F("CutflowTableHisto__other", "CutflowTableHisto__other" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eee__other = new TH1F("CutflowTableHisto_eee__other", "CutflowTableHisto_eee__other" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_eeu__other = new TH1F("CutflowTableHisto_eeu__other", "CutflowTableHisto_eeu__other" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uue__other = new TH1F("CutflowTableHisto_uue__other", "CutflowTableHisto_uue__other" , thesystlistnames.size(), -0.5, thesystlistnames.size()-0.5);
 TH1F*  CutflowTableHisto_uuu__other = new TH1F("CutflowTableHisto_uuu__other", "CutflowTableHisto_uuu__other" , thesystlistnames.size(), -0.5, v_cutflow.size()-0.5);
 */






/*
 
 thesystlist.push_back("btagSF_cferr1Down");
 thesystlist.push_back("btagSF_cferr2Down");
 thesystlist.push_back("btagSF_hfDown");
 thesystlist.push_back("btagSF_hfstats1Down");
 thesystlist.push_back("btagSF_hfstats2Down");
 thesystlist.push_back("btagSF_lfDown");
 thesystlist.push_back("btagSF_lfstats1Down");
 thesystlist.push_back("btagSF_lfstats2Down");
 */



///////////////////////////////////// MAIN CODE /////////////////////////////////////////
Int_t main(Int_t argc, char* argv[]){
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
  bool applyiniweights = false;
  
  Int_t testnr = -1;
  //////////// Settings of the analysis //////////////////
  bool makeerrorbands =false;
  for(Int_t i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./AnalyzerBDT [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   Zct / Zut: make plots for this one" << endl;
      std::cout << "   toppair: make plots for this one" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. MVAoutput/placeNtup " << endl;
      std::cout << "   Test nb: loop over nb events" << endl;
      std::cout << "   Data: add CRcut" << endl;
      std::cout << "   PSdata: generate pseudo data" << endl;
      std::cout << "   doSystematics: loop over systematics" << endl;
      std::cout << "   doPDFunc: calculate PDF unc" << endl;
      std::cout << "   PlotSystematics: make sys plots " << endl;
      std::cout << "   PlotMVAvars: plot mva vars" << endl;
      std::cout << "  doMTWtemplate make mtw templates" << endl;
      std::cout << "  doTTZtemplate make mtw templates" << endl;
      std::cout << "  doIniWeight " << endl;
      std::cout << "  makeerrorbands " << endl;
      std::cout << "  doPostFit " << endl;
      return 0;
    }
    if(string(argv[i]).find("doIniWeight")!=std::string::npos){
      applyiniweights  = true;
      
    }
    if(string(argv[i]).find("doPostFit")!=std::string::npos){
      doPostFit = true;
      
    }
    if(string(argv[i]).find("doErrorband")!=std::string::npos){
      doErrorband = true;
      
    }
    if(string(argv[i]).find("ErrorBand")!=std::string::npos){
      makeerrorbands  = true;
      
    }
    if(string(argv[i]).find("doMTWtemplate")!=std::string::npos){
      doMTWtemplate = true;
      cout << "******* making MTW templates *********" << endl;
    }
    if(string(argv[i]).find("doTTZtemplate")!=std::string::npos){
      doTTZtemplate = true;
      i++;
      if(string(argv[i]).find("STtemplate")!=std::string::npos) singletoptemplate = true;
      cout << "******* making TTZ templates *********" << endl;
    }
    if(string(argv[i]).find("doPDFunc")!=std::string::npos){
      doPDFunc = true;
      cout << "******* calculating pdf uncertainties *********" << endl;
    }
    if(string(argv[i]).find("doSystematics")!=std::string::npos) {
      doSystematics= true;
    }
    if(string(argv[i]).find("PlotJeSystematics")!=std::string::npos) {
      PlotJeSystematics= true;
      //cout << "plotting sys" << endl;
    }
    if(string(argv[i]).find("PlotSystematics")!=std::string::npos) {
      PlotSystematics= true;
      //cout << "plotting sys" << endl;
    }
    if(string(argv[i]).find("PSdata")!=std::string::npos) {
      doPseudoData = true;
    }
    if(string(argv[i]).find("Data")!=std::string::npos) {
      addData = true;
    }
    if(string(argv[i]).find("toppair")!=std::string::npos) {
      toppair = true;
      region = "toppair";
    }
    if(string(argv[i]).find("singletop")!=std::string::npos) {
      toppair = false;
      region = "singletop";
    }
    if(string(argv[i]).find("Zut")!=std::string::npos) {
      doZut = true;
      coupling = "Zut";
    }
    if(string(argv[i]).find("Zct")!=std::string::npos) {
      doZut = false;
      coupling = "Zct";
    }
    if(string(argv[i]).find("Test")!=std::string::npos) {
      testing = true;
      testnr = strtol(argv[i+1], NULL, 10);
      i++;
    }
    if(string(argv[i]).find("Ntup")!=std::string::npos) {
      placeNtup = argv[i+1];
      i++;
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    if(string(argv[i]).find("DetermineCut")!=string::npos) {
      DetermineCut= true;
    }
    if(string(argv[i]).find("Significane")!=string::npos) {
      CalculateSign = true;
    }
    if(string(argv[i]).find("PlotMVAvars")!=string::npos) {
      PlotMVAvars= true;
    }
  }
  string xmlFileName = "";
  if(doPostFit) xmlFileName = "config/Run2TriLepton_samples_analycop.xml" ;
  else if(doPDFunc) xmlFileName = "config/Run2TriLepton_samples_analypdfall.xml" ;
  else xmlFileName = "config/Run2TriLepton_samples_analy.xml" ;
  const char* xmlFile = xmlFileName.c_str();
  cout << " - Using config file " << xmlFile << endl;
  cout << " - Using mvatrees of " << placeNtup << endl;
  placeOutputReading = "MVAoutput/";
  mkdir(placeOutputReading.c_str(), 0777);
  placeOutputReading += "outputtemplates";
  mkdir(placeOutputReading.c_str(), 0777);
  if(!doMTWtemplate && !doTTZtemplate) placeOutputReading += "/" + coupling + "_" + region;
  else if(doMTWtemplate) placeOutputReading += "/MTWtemplate" ;
  else if(doTTZtemplate && singletoptemplate) placeOutputReading += "/STTTZtemplate" ;
  else if(doTTZtemplate && !singletoptemplate) placeOutputReading += "/TTTTZtemplate" ;
  mkdir(placeOutputReading.c_str(), 0777);
  combinetemplate_filename = placeOutputReading+"/Reader_" + coupling +"_" + region + ".root";
  if(doMTWtemplate ) combinetemplate_filename = placeOutputReading+"/Reader_"+coupling+"_MTW.root";
  if(doTTZtemplate ) combinetemplate_filename = placeOutputReading+"/Reader_"+coupling+"_TTZ.root";
  cout <<" - Combine templates stored at " << combinetemplate_filename.c_str() << endl;
  
  
  string postfitSFfilename = "postfit/Zut_PostFitBinnedSF_both.root";
  TFile *postfitSFfile = 0;
  if(doPostFit) postfitSFfile = new TFile(postfitSFfilename.c_str(),"READ");
  
  /*
   if(!toppair && !doZut){ BDT_begin= -0.4; BDT_end = 0.8; nbin = 15;}
   else if(toppair&& !doZut){ BDT_begin= -0.9; BDT_end = 0.8; nbin = 20;}
   else if(!toppair && doZut){ BDT_begin= -0.4; BDT_end = 0.7; nbin = 14;}
   else if(toppair&& doZut){ BDT_begin= -0.9; BDT_end = 0.9; nbin = 22;}
   */
  std::vector < int>  decayChannels = {0,1,2,3,-9}; // uuu uue eeu eee all
  vector <string> thesystlist;
  vector <string> thesystlistnames;
  thesystlistnames.clear();
  thesystlist.clear();
  thesystlist.push_back(""); // nominal
  if(doSystematics){
    cout << "pushing back systematics" << endl;
    thesystlist.push_back("puSFDown");
    thesystlist.push_back("electronSFDown");
    thesystlist.push_back("muonSFDown");
    thesystlist.push_back("btagSFDown");
    thesystlist.push_back("btagSF_cferr1Down");
    thesystlist.push_back("btagSF_cferr2Down");
    thesystlist.push_back("btagSF_hfDown");
    thesystlist.push_back("btagSF_hfstats1Down");
    thesystlist.push_back("btagSF_hfstats2Down");
    thesystlist.push_back("btagSF_lfDown");
    thesystlist.push_back("btagSF_lfstats1Down");
    thesystlist.push_back("btagSF_lfstats2Down");
    thesystlist.push_back("JERUp"); // to plot je inmpact only tunr these one with do systematics and plotsystematics
    thesystlist.push_back("JESUp");
    
    thesystlist.push_back("puSFUp");
    thesystlist.push_back("electronSFUp");
    thesystlist.push_back("muonSFUp");
    thesystlist.push_back("btagSF_cferr1Up");
    thesystlist.push_back("btagSF_cferr2Up");
    thesystlist.push_back("btagSF_hfUp");
    thesystlist.push_back("btagSF_hfstats1Up");
    thesystlist.push_back("btagSF_hfstats2Up");
    thesystlist.push_back("btagSF_lfUp");
    thesystlist.push_back("btagSF_lfstats1Up");
    thesystlist.push_back("btagSF_lfstats2Up");
   
    
    thesystlist.push_back("JERDown");
    
    thesystlist.push_back("JESDown");
  }
  cout << "Number of systematics " << thesystlist.size() << endl;
  //for plotting
  thesystlistnames.push_back("puSF");
  thesystlistnames.push_back("electronSF");
  thesystlistnames.push_back("muonSF");
  thesystlistnames.push_back("btagSF_cferr1");
  thesystlistnames.push_back("btagSF_cferr2");
  thesystlistnames.push_back("btagSF_hf");
  thesystlistnames.push_back("btagSF_hfstats1");
  thesystlistnames.push_back("btagSF_hfstats2");
  thesystlistnames.push_back("btagSF_lf");
  thesystlistnames.push_back("btagSF_lfstats1");
  thesystlistnames.push_back("btagSF_lfstats2");
  thesystlistnames.push_back("JES");
  thesystlistnames.push_back("JER");
  
  
  ///////////////:  load datasets
  datasets.clear();
  datasetsbf.clear();
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasetsbf, xmlFile);
  cout << "datasets loaded " << datasetsbf.size() << " samples" <<  endl;
  for (Int_t d = 0; d < datasetsbf.size(); d++){   //Loop through datasets to get lumi setting
    
    dataSetName = datasetsbf[d]->Name();
    // cout << "sample " << dataSetName << endl;
    if( dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
      Luminosity = datasetsbf[d]->EquivalentLumi();
      cout << " - lumi set to " << Luminosity << endl;
    }
    
    if((dataSetName.find("Zct")!=std::string::npos || dataSetName.find("zct")!=std::string::npos) && doZut && !doMTWtemplate && !doTTZtemplate){
      cout << " - removing " << dataSetName << " from samples" << endl;
      continue;
    }
    else if ((dataSetName.find("Zut")!=std::string::npos || dataSetName.find("zut")!=std::string::npos) && !doZut && !doMTWtemplate && !doTTZtemplate){
      cout << " - removing " << dataSetName << " from samples" << endl;
      continue;
    }
    else{
      
      if( (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ) && doPseudoData ) {
        cout << " - removing " << dataSetName << " from samples" << endl;
      }
      else if(dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
        datafound = true;
        //cout << "pushing back " << dataSetName << endl;
        datasets.push_back(datasetsbf[d]);
      }
      else {
        //cout << "pushing back " << dataSetName << endl;
        datasets.push_back(datasetsbf[d]);
      }
    }
  }
  
  cout << datasets.size() << " samples will be used " << endl;
  ///////////////// Initialisation ////////////////////
  
  //if(toppair) nbin = 50;
  if(makePlots && !doMTWtemplate && !doTTZtemplate){
    for(Int_t isys = 0; isys < thesystlist.size() ; isys++){
      systematic = thesystlist[isys];
      tempstring = region + "_"+coupling;
      if(isys != 0 ) tempstring += "_"+ systematic;
      InitMSPlots(tempstring, decayChannels, toppair, doZut);
    }
  }
  
  if(makePlots && doMTWtemplate && !doTTZtemplate){
    for(Int_t isys = 0; isys < thesystlist.size() ; isys++){
      systematic = thesystlist[isys];
      if(isys != 0 ) tempstring = "_" + systematic;
      InitMSPlotsMTW(tempstring, decayChannels);
    }
  }
  
  if(makePlots && !doMTWtemplate && doTTZtemplate){
    for(Int_t isys = 0; isys < thesystlist.size() ; isys++){
      systematic = thesystlist[isys];
      if(isys != 0 ) tempstring = "_" + systematic;
      InitMSPlotsTTZ(tempstring, decayChannels);
    }
  }
  
  
  ///////////// General function //////////
  if(DetermineCut && !doMTWtemplate && !doTTZtemplate){
    vector<double> v_cut = BDTCUT(region, coupling);
    return 0;
  }
  
  //////////// START LOOPING ON SYS - DATASETS - EVENTS //////////
  ///*****************///
  ///   MAIN CODE  ///
  ///*****************///
  ntupleFileName = placeNtup;
  if(!doMTWtemplate && !doTTZtemplate) fin = new TFile((ntupleFileName).c_str(),"READ");
  bool onlynomforsys = false;
  Int_t WZregionEntries = 0;
  
  
  TH1::SetDefaultSumw2();
  TH1F* hist_WZ = new TH1F("MTW_WZ","trans. mass W boson in WZ region: WZ (GeV)",           nbinMTW, 0., endMTW);
  TH1F* hist_TT_FCNC = new TH1F("MTW_TT_FCNC","transv. mass W boson in TTSR (GeV)",           nbinMTW, 0., endMTW);
  TH1F* hist_fakes = new TH1F("MTW_fakes","transv. mass W boson WZ region: fakes (GeV)",           nbinMTW, 0., endMTW);
  
  
  
  
  histo1DMTW["MTW_WZ"] = hist_WZ ;
  histo1DMTW["MTW_TT_FCNC"] = hist_TT_FCNC ;
  histo1DMTW["MTW_fakes"] = hist_fakes ;
  
  for(Int_t isys = 0; isys < thesystlist.size() ; isys++){
    systematic = thesystlist[isys];
    cout << "looking at " << systematic << " systematics " << endl;
    
    for (Int_t d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      nboffakesinWZ = 0.;
      nboffakesinTTSR = 0.;
      nboffakesinTTSRwith2jets = 0.;
      nboffakesinTTSRwith3jets = 0.;
      nboffakesinWZwith2jets = 0.;
      nboffakesinWZwith3jets = 0.;
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
      // settings
      isData = false;
      onlynomforsys = false;
      dataSetName = datasets[d]->Name();
      if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
        isData = true;
      }
      
      if ((dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos || dataSetName.find("fake")!=std::string::npos   ) && isys != 0){
        onlynomforsys = true;
      }
      
      
      postfix = "";
      if(systematic.find("JESDown")!=std::string::npos) postfix = "_JESdown";
      else if(systematic.find("JESUp")!=std::string::npos) postfix = "_JESup";
      else if(systematic.find("JERDown")!=std::string::npos) postfix = "_JERdown";
      else if(systematic.find("JERUp")!=std::string::npos) postfix = "_JERup";
      else postfix = "";
      if(onlynomforsys) postfix = "";
      
      if(doMTWtemplate || doTTZtemplate) {
        ntupleFileName = placeNtup + "MVA_tree_" + dataSetName + postfix + ".root";
        tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
      }
      
      if(!doMTWtemplate && !doTTZtemplate) tTreeName = "Control_"+dataSetName + postfix;
      else tTreeName = "mvatree" + postfix;
      /// Get data
      cout << "   treename " << tTreeName << " from " << ntupleFileName <<  endl;
      if(doMTWtemplate || doTTZtemplate){
        tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      }
      else if(!doMTWtemplate && !doTTZtemplate) tTree[dataSetName.c_str()] = (TTree*)fin->Get(tTreeName.c_str());
      nEntries = -1;
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries: " << nEntries << endl;
      
      // Initialise tree
      if(!doMTWtemplate && !doTTZtemplate) InitTree(tTree[dataSetName.c_str()], isData, toppair, doZut);
      else if(doMTWtemplate || doTTZtemplate) InitAnalyzerTree(tTree[dataSetName.c_str()]);
      
      if(toppair && !doMTWtemplate && !doTTZtemplate) MVA_region = 1;
      
      // Initialise plots
      if(!doMTWtemplate && PlotMVAvars && isys == 0 && !doTTZtemplate){
        Init1DHisto(dataSetName, systematic, toppair, doZut, decayChannels);
      }
      // initialise combine output histograms
      TH1::SetDefaultSumw2();
      TH2::SetDefaultSumw2();
      //cout << "create template histo" << endl;
      TH1F *hist_uuu(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_uuu = new TH1F( (coupling + "_" + region+"_uuu").c_str(),           (coupling + "_" + region+"_uuu").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_uuu = new TH1F( (coupling + "_mTW_uuu").c_str(),           (coupling + "_mTW_uuu").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_uuu = new TH1F( (coupling + "_Zmass_uuu").c_str(),           (coupling + "_Zmass_uuu").c_str(),           nbinTTZ,beginTTZ, endTTZ);
      TH1F *hist_uue(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_uue = new TH1F( (coupling + "_" + region+"_uue").c_str(),           (coupling + "_" + region+"_uue").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_uue = new TH1F( (coupling + "_mTW_uue").c_str(),           (coupling + "_mTW_uue").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_uue = new TH1F( (coupling + "_Zmass_uue").c_str(),           (coupling + "_Zmass_uue").c_str(),            nbinTTZ,beginTTZ, endTTZ );
      TH1F *hist_eeu(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_eeu = new TH1F( (coupling + "_" + region+"_eeu").c_str(),           (coupling + "_" + region+"_eeu").c_str(),           nbin,BDT_begin,BDT_end );
      else  if(doMTWtemplate ) hist_eeu = new TH1F( (coupling + "_mTW_eeu").c_str(),           (coupling + "_mTW_eeu").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate)  hist_eeu = new TH1F( (coupling + "_Zmass_eeu").c_str(),           (coupling + "_Zmass_eeu").c_str(),           nbinTTZ,beginTTZ, endTTZ);
      TH1F *hist_eee(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_eee = new TH1F( (coupling + "_" + region+"_eee").c_str(),           (coupling + "_" + region+"_eee").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_eee = new TH1F( (coupling + "_mTW_eee").c_str(),           (coupling + "_mTW_eee").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_eee = new TH1F( (coupling + "_Zmass_eee").c_str(),           (coupling + "_Zmass_eee").c_str(),            nbinTTZ,beginTTZ, endTTZ );
      
      TH1F *hist_check_uuu(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_check_uuu = new TH1F( ("check_" +coupling+ "_" + region+"_uuu").c_str(),           ("check_" +coupling+ "_" + region+"_uuu").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_check_uuu = new TH1F( ("check_" +coupling+ "_mTW_uuu").c_str(),           ("check_" +coupling+ "_mTW_uuu").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_check_uuu = new TH1F( ("check_" +coupling+ "_Zmass_uuu").c_str(),           ("check_" +coupling+ "_Zmass_uuu").c_str(),           nbinTTZ,beginTTZ, endTTZ);
      TH1F *hist_check_uue(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_check_uue = new TH1F( ("check_" +coupling+ "_" + region+"_uue").c_str(),           ("check_" +coupling+ "_" + region+"_uue").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_check_uue = new TH1F( ("check_" +coupling+ "_mTW_uue").c_str(),           ("check_" +coupling+ "_mTW_uue").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_check_uue = new TH1F( ("check_" +coupling+ "_Zmass_uue").c_str(),           ("check_" +coupling+ "_Zmass_uue").c_str(),            nbinTTZ,beginTTZ, endTTZ );
      TH1F *hist_check_eeu(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_check_eeu = new TH1F( ("check_" +coupling+ "_" + region+"_eeu").c_str(),           ("check_" +coupling+ "_" + region+"_eeu").c_str(),           nbin,BDT_begin,BDT_end );
      else  if(doMTWtemplate ) hist_check_eeu = new TH1F( ("check_" +coupling+ "_mTW_eeu").c_str(),           ("check_" +coupling+ "_mTW_eeu").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate)  hist_check_eeu = new TH1F( ("check_" +coupling+ "_Zmass_eeu").c_str(),           ("check_" +coupling+ "_Zmass_eeu").c_str(),           nbinTTZ,beginTTZ, endTTZ);
      TH1F *hist_check_eee(0);
      if(!doMTWtemplate && !doTTZtemplate) hist_check_eee = new TH1F( ("check_" +coupling+ "_" + region+"_eee").c_str(),           ("check_" +coupling+ "_" + region+"_eee").c_str(),           nbin,BDT_begin,BDT_end );
      else if(doMTWtemplate) hist_check_eee = new TH1F( ("check_" +coupling+ "_mTW_eee").c_str(),           ("check_" +coupling+ "_mTW_eee").c_str(),           nbinMTW,0, endMTW );
      else if(doTTZtemplate) hist_check_eee = new TH1F( ("check_" +coupling+ "_Zmass_eee").c_str(),           ("check_" +coupling+ "_Zmass_eee").c_str(),            nbinTTZ,beginTTZ, endTTZ );
      
      
      
      /*
      TH2F *hist_lep2_pt(0);
      TH2F *hist_lep2_pt_uuu(0);
      TH2F *hist_lep2_pt_uue(0);
      TH2F *hist_lep2_pt_eeu(0);
      TH2F *hist_lep2_pt_eee(0);
      
      hist_lep2_pt = new TH2F(("hist_lep2_pt"+coupling+"_"+region+"_"+dataSetName).c_str(),("hist_lep2_pt"+coupling+"_"+region+"_"+dataSetName).c_str(), 25,0.,150.,nbin, BDT_begin, BDT_end);
      hist_lep2_pt_uuu = new TH2F(("hist_lep2_pt"+coupling+"_"+region+"_uuu_"+dataSetName).c_str(),("hist_lep2_pt"+coupling+"_"+region+"_uuu_"+dataSetName).c_str(), 25,0.,150.,nbin, BDT_begin, BDT_end);
      hist_lep2_pt_uue = new TH2F(("hist_lep2_pt"+coupling+"_"+region+"_uue_"+dataSetName).c_str(),("hist_lep2_pt"+coupling+"_"+region+"_uue_"+dataSetName).c_str(), 25,0.,150.,nbin, BDT_begin, BDT_end);
      hist_lep2_pt_eeu = new TH2F(("hist_lep2_pt"+coupling+"_"+region+"_eeu_"+dataSetName).c_str(),("hist_lep2_pt"+coupling+"_"+region+"_eeu_"+dataSetName).c_str(), 25,0.,150.,nbin, BDT_begin, BDT_end);
      hist_lep2_pt_eee = new TH2F(("hist_lep2_pt"+coupling+"_"+region+"_eee_"+dataSetName).c_str(),("hist_lep2_pt"+coupling+"_"+region+"_eee_"+dataSetName).c_str(), 25,0.,150.,nbin, BDT_begin, BDT_end);
      */
      //cout << "created template histo" << endl;
      /// Initialise WZ plots
      if((dataSetName.find("WZTo3LNu_3Jets_MLL50_80X")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("TTZ")!=std::string::npos || dataSetName.find("tZq")!=std::string::npos || dataSetName.find("ZZTo4")!=std::string::npos || dataSetName.find("NP")!=std::string::npos) && doPDFunc ){
        InitCalculatePDFWeightHisto(dataSetName, doMTWtemplate, doTTZtemplate);
      }
      
      
      
      
      if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && doMTWtemplate){
        // InitMTWShapeHisto(dataSetName, systematic, isys, decayChannels);
        
      }
      
      // safeties
      if(!doSystematics && isys != 0) continue;
      Int_t endEvent =nEntries;
      if(testing){
        if(endEvent > testnr) endEvent = testnr;
      }
      /// loop on events
      Double_t weight = 1.;
      WZregionEntries = 0;
      Int_t WZregionEntriesuuu = 0;
      for (Int_t ievt = 0; ievt < endEvent; ievt++)
      {
        if (ievt%100 == 0)
          std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
        
        
        /// Load event
        tTree[(dataSetName).c_str()]->GetEntry(ievt);
        /*if(datafound && MVA_BDT > -0.68 && !doMTWtemplate && !toppair && doZut){ continue;}
         else if(datafound && MVA_BDT > -0.6 && !doMTWtemplate && !toppair && !doZut){ continue;}
         else if(datafound && MVA_BDT > -0.12 && !doMTWtemplate && toppair && !doZut){ continue;}
         else if(datafound && MVA_BDT > -0.2 && !doMTWtemplate && toppair && doZut){ continue;}*/
        // cout << "region " <<MVA_region << endl;
        if(dataSetName.find("fake")!=std::string::npos || dataSetName.find("WZ")!=std::string::npos){
          if(MVA_region == 2){
            nboffakesinWZ++;
            if(MVA_nJets == 2)nboffakesinWZwith2jets++;
            if(MVA_nJets == 3)nboffakesinWZwith3jets++;
          }
        if(MVA_region == 1){
          nboffakesinTTSR++;
          if(MVA_nJets == 2) nboffakesinTTSRwith2jets++;
          if(MVA_nJets == 3) nboffakesinTTSRwith3jets++;
        }
        }
        
        if(doMTWtemplate &&MVA_region != 2 ){ continue ;} // only in WZ control region}
        else if(doMTWtemplate) { WZregionEntries++; }
        
        if(doTTZtemplate && MVA_region != 3 && MVA_region != 4  ){ continue ;} // only in TTZ control region}
        if(doTTZtemplate && singletoptemplate && MVA_region != 4 ){continue;}
        if(doTTZtemplate && !singletoptemplate && MVA_region != 3){continue;}
        
        if(doTTZtemplate && (MVA_Zboson_M > 76. && MVA_Zboson_M < 106.) ){continue;} // remove inner signal
        
        if(doMTWtemplate &&MVA_region == 2 && MVA_channel == 0) { WZregionEntriesuuu++;}
        // MVA_weight_nloSF = 1.;
        weight = 1.;
        if(!isData && dataSetName.find("fake")==std::string::npos)  weight = MVA_Luminosity /MVA_EqLumi;
        if(!isData && dataSetName.find("fake")==std::string::npos && !onlynomforsys && isys != 0 && (systematic.find("JES")==std::string::npos) && (systematic.find("JER")==std::string::npos) ){
          //cout <<  "getting weight for " << systematic << endl;
          if(systematic.find("puSFUp")) weight *= MVA_weight_puSF_up * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("puSFDown")) weight *= MVA_weight_puSF_down * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("electronSFUp")) weight *= MVA_weight_puSF * MVA_weight_electronSF_up * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("electronSFDown")) weight *= MVA_weight_puSF * MVA_weight_electronSF_down * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("muonSFUp")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF_up * MVA_weight_nloSF;
          else if(systematic.find("muonSFDown")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF_down * MVA_weight_nloSF;
          else if(systematic.find("btagSF_cferr1Up")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_cferr1_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_cferr1Down")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_cferr1_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_cferr2Up")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_cferr2_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_cferr2Down")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_cferr2_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfUp")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lf_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfDown")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lf_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfUp")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hf_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfDown")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hf_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfstats1Up")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfstats1Down")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfstats2Up")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_hfstats2Down")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_hfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfstats2Up")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfstats2Down")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfstats1Up")) weight *= MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF;
          else if(systematic.find("btagSF_lfstats1Down")) weight *=MVA_weight_puSF * MVA_weight_electronSF * MVA_weight_btagSF_lfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF;
          
        }
        else if(isData ) weight = 1.;
        else weight *= MVA_weight_nom;
        
        /*
        if(dataSetName.find("fake")!=std::string::npos && (MVA_channel == 0 || MVA_channel == 2)){ weight *= 0.545 ;}
        if(dataSetName.find("fake")!=std::string::npos && (MVA_channel == 1 || MVA_channel == 3)){ weight *= 0.590;}
        if(dataSetName.find("WZ")!=std::string::npos ){ weight *=0.841 ;}*/
        //if(dataSetName.find("fake")!=std::string::npos) weight *= 0.0001;
        
        if(applyiniweights){
                   if(doMTWtemplate){
          if(dataSetName.find("fake")!=std::string::npos  && (MVA_channel == 0 || MVA_channel == 1)){
            weight *= 1.825;
          }
          else if(dataSetName.find("fake")!=std::string::npos  && (MVA_channel == 2 || MVA_channel == 3)){
            weight *= 1.321;
          }
          else if(dataSetName.find("WZT")!=std::string::npos  ){
            weight *= 1.235;
          }
          }
          if(doTTZtemplate && ! singletoptemplate){
            if(dataSetName.find("fake")!=std::string::npos  && (MVA_channel == 0 || MVA_channel == 1)){
              weight *= 0.926; // fake el in Z
            }
            else if(dataSetName.find("fake")!=std::string::npos  && (MVA_channel == 2 || MVA_channel == 3)){
              weight *= 1.076;
            }
          }
         /* if(dataSetName.find("WZT")!=std::string::npos && MVA_channel == 2){ weight *= 1.776 ;}
          if(dataSetName.find("WZT")!=std::string::npos && MVA_channel == 1){
            weight *= 1.52;
          }
          if(dataSetName.find("WZT")!=std::string::npos && MVA_channel == 3){ weight *= 1.576;}*/
        }
        
        
        if(dataSetName.find("nonpromptwrong")!=std::string::npos) hist_BDT_tt_nonpromptinZ->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptcorrect")!=std::string::npos) hist_BDT_tt_nonpromptinW->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptwrong")!=std::string::npos && MVA_channel == 0) hist_BDT_tt_uuu_nonpromptinZ->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptcorrect")!=std::string::npos && MVA_channel == 0) hist_BDT_tt_uuu_nonpromptinW->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptwrong")!=std::string::npos && MVA_channel == 1) hist_BDT_tt_uue_nonpromptinZ->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptcorrect")!=std::string::npos && MVA_channel == 1) hist_BDT_tt_uue_nonpromptinW->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptwrong")!=std::string::npos && MVA_channel == 2) hist_BDT_tt_eeu_nonpromptinZ->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptcorrect")!=std::string::npos && MVA_channel == 2) hist_BDT_tt_eeu_nonpromptinW->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptwrong")!=std::string::npos && MVA_channel == 3) hist_BDT_tt_eee_nonpromptinZ->Fill(MVA_BDT, weight);
        if(dataSetName.find("nonpromptcorrect")!=std::string::npos && MVA_channel == 3) hist_BDT_tt_eee_nonpromptinW->Fill(MVA_BDT, weight);
        
        
        if(dataSetName.find("light")!=std::string::npos)  hist_BDT_WZ_light->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xcc")!=std::string::npos )   hist_BDT_WZ_c ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xbb")!=std::string::npos )   hist_BDT_WZ_b ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("light")!=std::string::npos && MVA_channel == 0)  hist_BDT_WZ_uuu_light->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xcc")!=std::string::npos && MVA_channel == 0)   hist_BDT_WZ_uuu_c ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xbb")!=std::string::npos && MVA_channel == 0)   hist_BDT_WZ_uuu_b ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        
        if(dataSetName.find("light")!=std::string::npos && MVA_channel == 3) hist_BDT_WZ_eee_light->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xcc")!=std::string::npos && MVA_channel == 3)   hist_BDT_WZ_eee_c ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xbb")!=std::string::npos && MVA_channel == 3)   hist_BDT_WZ_eee_b ->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("light")!=std::string::npos && MVA_channel == 2)  hist_BDT_WZ_eeu_light->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xcc")!=std::string::npos && MVA_channel == 2)   hist_BDT_WZ_eeu_c->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xbb")!=std::string::npos && MVA_channel == 2)   hist_BDT_WZ_eeu_b->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("light")!=std::string::npos && MVA_channel == 1)  hist_BDT_WZ_uue_light->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xcc")!=std::string::npos && MVA_channel == 1)   hist_BDT_WZ_uue_c->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        if(dataSetName.find("WZTo3LNu_amc_new_80Xbb")!=std::string::npos && MVA_channel == 1)   hist_BDT_WZ_uue_b->Fill(MVA_BDT, weight*MVA_weightWZcorr);
        
       
        
        
        if(!doMTWtemplate && !doTTZtemplate && (dataSetName.find("fake")==std::string::npos)){
          if(MVA_channel== 0) 		{hist_uuu->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 1) {hist_uue->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 2) {hist_eeu->Fill( MVA_BDT, weight);}
          else if(MVA_channel == 3) {hist_eee->Fill( MVA_BDT, weight);}
        }
        else if(doMTWtemplate ){
          if(MVA_channel== 0) 		{hist_uuu->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 1) {hist_uue->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 2) {hist_eeu->Fill( MVA_mWt, weight);}
          else if(MVA_channel == 3) {hist_eee->Fill( MVA_mWt, weight);}
        }
        else if(doTTZtemplate  && dataSetName.find("fake") ==std::string::npos){
          if(MVA_channel== 0) 		{hist_uuu->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel== 1) {hist_uue->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel== 2) {hist_eeu->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel == 3) {hist_eee->Fill( MVA_Zboson_M, weight);}
        }
        else if(!doMTWtemplate && !doTTZtemplate && dataSetName.find("fake") !=std::string::npos){
          hist_uuu->Fill( MVA_BDT, weight);
          hist_uue->Fill( MVA_BDT, weight);
          hist_eeu->Fill( MVA_BDT, weight);
          hist_eee->Fill( MVA_BDT, weight);
          
        }
        else if(doTTZtemplate  && dataSetName.find("fake") !=std::string::npos){
          hist_uuu->Fill( MVA_Zboson_M, weight);
          hist_uue->Fill( MVA_Zboson_M, weight);
          hist_eeu->Fill( MVA_Zboson_M, weight);
          hist_eee->Fill( MVA_Zboson_M, weight);
        }
        
        if(!doMTWtemplate && !doTTZtemplate){
          if(MVA_channel== 0) 		{hist_check_uuu->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 1) {hist_check_uue->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 2) {hist_check_eeu->Fill( MVA_BDT, weight);}
          else if(MVA_channel == 3) {hist_check_eee->Fill( MVA_BDT, weight);}
        }
        else if(doMTWtemplate ){
          if(MVA_channel== 0) 		{hist_check_uuu->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 1) {hist_check_uue->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 2) {hist_check_eeu->Fill( MVA_mWt, weight);}
          else if(MVA_channel == 3) {hist_check_eee->Fill( MVA_mWt, weight);}
        }
        else if(doTTZtemplate  ){
          if(MVA_channel== 0) 		{hist_check_uuu->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel== 1) {hist_check_uue->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel== 2) {hist_check_eeu->Fill( MVA_Zboson_M, weight);}
          else if(MVA_channel == 3) {hist_check_eee->Fill( MVA_Zboson_M, weight);}
        }
        
       /* hist_lep2_pt->Fill(MVA_lepton2_pt,MVA_BDT,weight);
        if(MVA_channel== 0) hist_lep2_pt_uuu->Fill(MVA_lepton2_pt,MVA_BDT,weight);
        if(MVA_channel== 1) hist_lep2_pt_uue->Fill(MVA_lepton2_pt,MVA_BDT,weight);
       if(MVA_channel== 2)  hist_lep2_pt_eeu->Fill(MVA_lepton2_pt,MVA_BDT,weight);
        if(MVA_channel== 3) hist_lep2_pt_eee->Fill(MVA_lepton2_pt,MVA_BDT,weight);
        */
        
        
        // for MS plots
        Double_t weightMSPlot = weight;//*MVA_weightWZcorr;
        string channelname = "";
        string regionname = "";
        string dataSetTitel = datasets[d]->Title();
        if(dataSetTitel.find("ZFake")!=std::string::npos && (MVA_channel == 0 || MVA_channel == 1)) dataSetTitel = "ZFakeMu";
        else if(dataSetTitel.find("ZFake")!=std::string::npos && (MVA_channel == 2 || MVA_channel == 3)) dataSetTitel = "ZFakeEl";
        else if(dataSetTitel.find("WFake")!=std::string::npos && (MVA_channel == 0 || MVA_channel == 2)) dataSetTitel = "WFakeMu";
        else if(dataSetTitel.find("WFake")!=std::string::npos && (MVA_channel == 1 || MVA_channel == 3)) dataSetTitel = "WFakeEl";
        
        if(doPostFit && !isData && dataSetName.find("NP")==std::string::npos){
          if(MVA_channel == 0) channelname = "_LepChan_3mu_";
          else if(MVA_channel == 1) channelname = "_LepChan_1e2mu_";
          else if(MVA_channel == 2) channelname = "_LepChan_2e1mu_";
          else if(MVA_channel == 3) channelname = "_LepChan_3e_";
          
          if(toppair) regionname = "TTSR";
          else regionname = "STSR";
          
          //cout << "getting " << ("SF_"+ dataSetTitel +channelname+regionname).c_str() << " from " << postfitSFfilename.c_str()  << endl;
          TH1F* h_postfitSF = (TH1F*)( postfitSFfile->Get(("SF_"+ dataSetTitel +channelname+regionname).c_str())->Clone("h_postfitSF"));
          if(h_postfitSF == 0) cout << "histo not found" << endl;
          // else cout << "got histo " << endl;
          double postfitSF = h_postfitSF->GetBinContent( h_postfitSF->GetXaxis()->FindBin(MVA_BDT));
          weightMSPlot *= postfitSF;
          
        }

        if(isData || dataSetName.find("fake")!=std::string::npos) weightMSPlot *= MVA_Luminosity;
        /// Fill plots
        if(doPDFunc){
          if(dataSetName.find("WZTo3LNu_3Jets_MLL50_80X")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("TTZ")!=std::string::npos || dataSetName.find("tZq")!=std::string::npos || dataSetName.find("ZZTo4")!=std::string::npos || dataSetName.find("NP")!=std::string::npos ) CalculatePDFWeight(dataSetName, MVA_BDT,MVA_weight_nom, MVA_channel,doMTWtemplate,doTTZtemplate);
        }
        if(PlotMVAvars  && isys == 0 && !doMTWtemplate && !doTTZtemplate){
          Fill1DHisto(dataSetName, systematic, toppair, doZut, decayChannels, weight, MVA_channel);
        }
        if(!isData && PlotJeSystematics && !doMTWtemplate && !doTTZtemplate){
          if(systematic.find("JESDown")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JES_down_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JESUp")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JES_up_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERDown")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JER_down_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERUp")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JER_up_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JES_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")!=std::string::npos) hist_BDT_JER_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          
          if(systematic.find("JESDown")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JES_down_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JESUp")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JES_up_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERDown")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JER_down_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERUp")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JER_up_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JES_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_BDT_JER_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
        }
        else if(!isData && PlotJeSystematics && doMTWtemplate ){
          if(systematic.find("JESDown")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JES_down_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JESUp")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JES_up_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERDown")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JER_down_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERUp")!=std::string::npos && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JER_up_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JES_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")!=std::string::npos) hist_mWt_JER_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          
          if(systematic.find("JESDown")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JES_down_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JESUp")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JES_up_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERDown")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JER_down_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(systematic.find("JERUp")!=std::string::npos && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JER_up_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JES_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
          if(isys ==0 && dataSetName.find("FCNC")==std::string::npos && dataSetName.find("fake")==std::string::npos) hist_mWt_JER_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
        }
        
        //cout << "booleans " <<  PlotSystematics << " "<< !isData << " " <<  !doMTWtemplate << endl;
        if(PlotSystematics && !isData && dataSetName.find("fake")==std::string::npos && !doMTWtemplate && !doSystematics && !doTTZtemplate){
          if(dataSetName.find("FCNC")!=std::string::npos){
            hist_BDT_puSF_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_puSF_up_sig->Fill(MVA_BDT, MVA_weight_puSF_up * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_puSF_down_sig->Fill(MVA_BDT, MVA_weight_puSF_down * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_BDT_electronSF_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_electronSF_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF_up * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_electronSF_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF_down * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_muonSF_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_muonSF_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_up * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_muonSF_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_down * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_BDT_btagSF_cferr1_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr1_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr1_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_cferr2_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr2_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr2_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hf_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hf_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hf_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hfstats1_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats1_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats1_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hfstats2_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats2_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats2_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lf_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lf_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lf_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lfstats1_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats1_up_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_up* MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats1_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lfstats2_nom_sig->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats2_up_sig->Fill(MVA_BDT, MVA_weight_nom*MVA_weight_btagSF_lfstats2_up/MVA_weight_btagSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats2_down_sig->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
          }
          else{
            
            hist_BDT_puSF_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_puSF_up_bkg->Fill(MVA_BDT, MVA_weight_puSF_up * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_puSF_down_bkg->Fill(MVA_BDT, MVA_weight_puSF_down * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_BDT_electronSF_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_electronSF_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF_up * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_electronSF_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF_down * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_muonSF_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_muonSF_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_up * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_muonSF_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_down * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_BDT_btagSF_cferr1_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr1_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr1_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_cferr2_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr2_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_cferr2_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hf_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hf_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hf_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hfstats1_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats1_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats1_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_hfstats2_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats2_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_hfstats2_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lf_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lf_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lf_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lfstats1_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats1_up_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_up* MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats1_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_BDT_btagSF_lfstats2_nom_bkg->Fill(MVA_BDT,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats2_up_bkg->Fill(MVA_BDT, MVA_weight_nom*MVA_weight_btagSF_lfstats2_up/MVA_weight_btagSF*MVA_Luminosity /MVA_EqLumi );
            hist_BDT_btagSF_lfstats2_down_bkg->Fill(MVA_BDT, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
          }
          
          
        }
        else if(PlotSystematics && !isData && dataSetName.find("FCNC")==std::string::npos && (doMTWtemplate )&& !doSystematics){
          
          if(dataSetName.find("fake")!=std::string::npos){
            hist_mWt_puSF_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_puSF_up_sig->Fill(MVA_mWt, MVA_weight_puSF_up * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_puSF_down_sig->Fill(MVA_mWt, MVA_weight_puSF_down * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_mWt_electronSF_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_electronSF_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF_up * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_electronSF_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF_down * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_muonSF_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_muonSF_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_up * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_muonSF_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_down * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_cferr1_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr1_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr1_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_cferr2_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr2_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr2_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hf_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hf_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hf_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hfstats1_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats1_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats1_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hfstats2_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats2_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats2_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lf_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lf_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lf_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lfstats1_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats1_up_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_up* MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats1_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lfstats2_nom_sig->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats2_up_sig->Fill(MVA_mWt, MVA_weight_nom*MVA_weight_btagSF_lfstats2_up/MVA_weight_btagSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats2_down_sig->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            // cout << "filling" <<endl;
          }
          else{
            hist_mWt_puSF_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_puSF_up_bkg->Fill(MVA_mWt,MVA_weight_puSF_up * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_puSF_down_bkg->Fill(MVA_mWt,MVA_weight_puSF_down * MVA_weight_electronSF * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_electronSF_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_electronSF_up_bkg->Fill(MVA_mWt,MVA_weight_puSF * MVA_weight_electronSF_up * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_electronSF_down_bkg->Fill(MVA_mWt,MVA_weight_puSF * MVA_weight_electronSF_down * MVA_weight_btagSF * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_muonSF_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_muonSF_up_bkg->Fill(MVA_mWt,MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_up * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_muonSF_down_bkg->Fill(MVA_mWt,MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF * MVA_weight_muonSF_down * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            
            hist_mWt_btagSF_cferr1_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr1_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr1_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_cferr2_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr2_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_cferr2_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_cferr2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hf_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hf_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hf_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hfstats1_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats1_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats1_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_hfstats2_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats2_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_hfstats2_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_hfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lf_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lf_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_up * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lf_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lf_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lfstats1_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats1_up_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_up* MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats1_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats1_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            hist_mWt_btagSF_lfstats2_nom_bkg->Fill(MVA_mWt,MVA_weight_nom*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats2_up_bkg->Fill(MVA_mWt, MVA_weight_nom*MVA_weight_btagSF_lfstats2_up/MVA_weight_btagSF*MVA_Luminosity /MVA_EqLumi );
            hist_mWt_btagSF_lfstats2_down_bkg->Fill(MVA_mWt, MVA_weight_puSF * MVA_weight_electronSF* MVA_weight_btagSF_lfstats2_down * MVA_weight_muonSF * MVA_weight_nloSF*MVA_Luminosity /MVA_EqLumi );
            
            // cout << "filling" <<endl;
          }
          
          
        }
        
        
        
        if (makePlots && !doMTWtemplate && !doTTZtemplate)
        {
          //cout << "ievt " << ievt << endl;
          tempstring = region + "_"+coupling;
          if(isys != 0) tempstring += "_"+ systematic;
          FillGeneralPlots(d, tempstring, decayChannels, doZut, toppair, weightMSPlot, MVA_channel, MVA_Zboson_M, MVA_NJets_CSVv2M,doPostFit);
        }
        if (makePlots &&( doMTWtemplate ))
        {
          if(isys != 0) tempstring = "_" + systematic;
          // if(isData) cout << "fill data " << endl;
          FillMTWPlots(d, tempstring, decayChannels, weightMSPlot, MVA_channel, MVA_Zboson_M, MVA_NJets_CSVv2M);
        }
        if (makePlots &&( doTTZtemplate ))
        {
          if(isys != 0) tempstring = "_" + systematic;
          // if(isData) cout << "fill data " << endl;
          FillTTZPlots(d, tempstring, decayChannels, weightMSPlot, MVA_channel, MVA_Zboson_M, MVA_NJets_CSVv2M);
        }
  
        
        if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos)&& isys == 0 && doMTWtemplate) histo1DMTW["MTW_WZ"]->Fill(MVA_mWt, weight);
        if(dataSetName.find("fake")!=std::string::npos && isys == 0 && doMTWtemplate) histo1DMTW["MTW_fakes"]->Fill(MVA_mWt, weight);
        if(dataSetName.find("TT_FCNC")!=std::string::npos && isys == 0 && doMTWtemplate) histo1DMTW["MTW_TT_FCNC"]->Fill(MVA_mWt, weight);
        
        
        
      } // events
      
      /*TCanvas* lept2 = new TCanvas();
      lept2->cd();
      hist_lep2_pt->Draw("colz");
      lept2->SaveAs(("hist_lep2_pt"+coupling+"_"+region+"_"+dataSetName+".png").c_str());
      hist_lep2_pt_uuu->Draw("colz");
      lept2->SaveAs(("hist_lep2_pt"+coupling+"_"+region+"_uuu_"+dataSetName+".png").c_str());
      hist_lep2_pt_uue->Draw("colz");
      lept2->SaveAs(("hist_lep2_pt"+coupling+"_"+region+"_uue_"+dataSetName+".png").c_str());
      hist_lep2_pt_eeu->Draw("colz");
      lept2->SaveAs(("hist_lep2_pt"+coupling+"_"+region+"_eeu_"+dataSetName+".png").c_str());
      hist_lep2_pt_eee->Draw("colz");
      lept2->SaveAs(("hist_lep2_pt"+coupling+"_"+region+"_eee_"+dataSetName+".png").c_str());
      
      delete lept2;
      */
      
      cout << "jet numbers" << endl;
      cout << "#WZ " << nboffakesinWZ  << " #2j " << nboffakesinWZwith2jets << " #3j " << nboffakesinWZwith3jets << endl;
      cout << "#TTSR " << nboffakesinTTSR << " #2j " << nboffakesinTTSRwith2jets << " #3j " << nboffakesinTTSRwith3jets << endl;
      
      
      cout << endl;
      if(doMTWtemplate) cout << "                WZ entries " << WZregionEntries << " uuu " << WZregionEntriesuuu << endl;
      /// Write combine histograms
      // --- Write histograms
      //cout << "DATASET " << dataSetName << " ISYS " << isys << endl;
      if((dataSetName.find("data")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && isys != 0) {
        delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
        delete  hist_check_eee; delete hist_check_uuu; delete hist_check_uue; delete hist_check_eeu;
        continue;
      }
      if(doMTWtemplate && WZregionEntries == 0) {
        delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
        //delete  hist_check_eee; delete hist_check_uuu; delete hist_check_uue; delete hist_check_eeu;
        continue;
      }
      //cout << "making template" << endl;
      TFile* combinetemplate_file(0);
      
      if(d == 0 && isys == 0) combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "RECREATE" );
      else combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "UPDATE" );
      combinetemplate_file->cd();
      
      //cout << "opened " << combinetemplate_filename.c_str() << endl;
      //NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)] FIX ME
      output_histo_name = "";
      if(!doMTWtemplate && !doTTZtemplate){
        if (dataSetName.find("fake")!=std::string::npos   ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          scaleFakes_uuu = hist_check_uuu->Integral()/hist_uuu->Integral();
          scaleFakes_uue = hist_check_uue->Integral()/hist_uue->Integral();
          scaleFakes_eeu = hist_check_eeu->Integral()/hist_eeu->Integral();
          scaleFakes_eee = hist_check_eee->Integral()/hist_eee->Integral();
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uuu_FakeMu_80X"  ;
          
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_uuu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_uuu_ZFakeMu_80X"  ;
            
          }
          
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Scale(hist_check_uuu->Integral()/hist_uuu->Integral());
          /*if(!toppair && !doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_uuu->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uue_FakeEl_80X"  ;
          
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_WFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_uue_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_uue_ZFakeMu_80X"  ;
            
          }
          
          hist_uue->SetTitle(output_histo_name.c_str());
          /*if(!toppair && !doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_uue->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_uue->Scale(hist_check_uue->Integral()/hist_uue->Integral());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eeu_FakeMu_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eeu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eeu_ZFakeEl_80X"  ;
            
          }
          /*if(!toppair && !doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_eeu->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Scale(hist_check_eeu->Integral()/hist_eeu->Integral());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eee_FakeEl_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_WFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eee_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eee_ZFakeEl_80X"  ;
            
          }
          /* if(!toppair && !doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_eee->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Scale(hist_check_eee->Integral()/hist_eee->Integral());
          hist_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uuu_"  + dataSetName ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          /* if(!toppair && !doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_uuu->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_uuu->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uue_"  + dataSetName ;
          /*if(!toppair && !doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_uue->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_uue->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eeu_"  + dataSetName ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eeu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eeu_ZFakeEl_80X"  ;
            
          }
          hist_eeu->SetTitle(output_histo_name.c_str());
          /* if(!toppair && !doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_eeu->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_eeu->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eee_"  + dataSetName ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_WFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_ee_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = coupling + "_BDT_" + region+"_eee_ZFakeEl_80X"  ;
            
          }
          hist_eee->SetTitle(output_histo_name.c_str());
          /* if(!toppair && !doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.4,0.8); }
           else if(toppair&& !doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.9,0.8);}
           else if(!toppair && doZut){hist_eee->GetXaxis()->SetRangeUser(-0.4,0.7);}
           else if(toppair&& doZut){ hist_eee->GetXaxis()->SetRangeUser(-0.9,0.9);}*/
          hist_eee->Write(output_histo_name.c_str());
        }
      }
      else if(doMTWtemplate){
       if (dataSetName.find("fake")!=std::string::npos  ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = "MTW_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "MTW_uuu_FakeMu_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_uuu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = "MTW_uuu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_uuu_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = "MTW_uuu_ZFakeMu_80X"  ;
            
          }
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = "MTW_uue_FakeEl_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_uue_WFakeEl_80X_"  + systematic ;
            else output_histo_name = "MTW_uue_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_uue_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = "MTW_uue_ZFakeMu_80X"  ;
            
          }
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "MTW_eeu_FakeMu_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_eeu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = "MTW_eeu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_eeu_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = "MTW_eeu_ZFakeEl_80X"  ;
            
          }
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "MTW_eee_FakeEl_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_eee_WFakeEl_80X_"  + systematic ;
            else output_histo_name = "MTW_eee_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "MTW_eee_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = "MTW_eee_ZFakeEl_80X"  ;
            
          }
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "MTW_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_uuu_"  + dataSetName ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_uue_"  + dataSetName ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_eeu_"  + dataSetName ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_eee_"  + dataSetName ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
      }
      else if(doTTZtemplate){
       if (dataSetName.find("fake")!=std::string::npos) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          scaleFakes_uuu = hist_check_uuu->Integral()/hist_uuu->Integral();
          scaleFakes_uue = hist_check_uue->Integral()/hist_uue->Integral();
          scaleFakes_eeu = hist_check_eeu->Integral()/hist_eeu->Integral();
          scaleFakes_eee = hist_check_eee->Integral()/hist_eee->Integral();
          
          
          if(isys!=0) output_histo_name = "Zmass_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "Zmass_uuu_FakeMu_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_uuu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = "Zmass_uuu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_uuu_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = "Zmass_uuu_ZFakeMu_80X"  ;
            
          }
          hist_uuu->Scale(hist_check_uuu->Integral()/hist_uuu->Integral());
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_uue_FakeMu_80X_"  + systematic ;
          else output_histo_name = "Zmass_uue_FakeMu_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_uue_WFakeEl_80X_"  + systematic ;
            else output_histo_name = "Zmass_uue_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_uue_ZFakeMu_80X_"  + systematic ;
            else output_histo_name = "Zmass_uue_ZFakeMu_80X"  ;
            
          }
          hist_uue->Scale(hist_check_uue->Integral()/hist_uue->Integral());
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_eeu_FakeEl_80X_"  + systematic ;
          else output_histo_name = "Zmass_eeu_FakeEl_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_eeu_WFakeMu_80X_"  + systematic ;
            else output_histo_name = "Zmass_eeu_WFakeMu_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_eeu_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = "Zmass_eeu_ZFakeEl_80X"  ;
            
          }
          hist_eeu->Scale(hist_check_eeu->Integral()/hist_eeu->Integral());
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "Zmass_eee_FakeEl_80X"  ;
          if(dataSetName.find("nonpromptinW")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_eee_WFakeEl_80X_"  + systematic ;
            else output_histo_name = "Zmass_eee_WFakeEl_80X"  ;
            
          }
          else if(dataSetName.find("nonpromptinZ")!=std::string::npos){
            if(isys!=0) output_histo_name = "Zmass_eee_ZFakeEl_80X_"  + systematic ;
            else output_histo_name = "Zmass_eee_ZFakeEl_80X"  ;
            
          }
          hist_eee->Scale(hist_check_eee->Integral()/hist_eee->Integral());
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "Zmass_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "Zmass_uuu_"  + dataSetName ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "Zmass_uue_"  + dataSetName ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "Zmass_eeu_"  + dataSetName ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "Zmass_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "Zmass_eee_"  + dataSetName ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
      }

      if(!doMTWtemplate && !doTTZtemplate){
        if (dataSetName.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_uuu_FakeMu_80X"  ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          // hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_uue_FakeEl_80X"  ;
          hist_check_uue->SetTitle(output_histo_name.c_str());
          //hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_eeu_FakeMu_80X"  ;
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          // hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_eee_FakeEl_80X"  ;
          hist_check_eee->SetTitle(output_histo_name.c_str());
          //hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_uuu_"  + dataSetName ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          // hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_uue_"  + dataSetName ;
          hist_check_uue->SetTitle(output_histo_name.c_str());
          // hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_eeu_"  + dataSetName ;
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          // hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "check" + coupling + "_BDT_" + region+"_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "check" + coupling + "_BDT_" + region+"_eee_"  + dataSetName ;
          hist_check_eee->SetTitle(output_histo_name.c_str());
          //hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->Write(output_histo_name.c_str());
        }
      }
      else if(doMTWtemplate){
        if (dataSetName.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = "checkMTW_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "checkMTW_uuu_FakeMu_80X"  ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          //hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = "checkMTW_uue_FakeEl_80X"  ;
          //hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->SetTitle(output_histo_name.c_str());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "checkMTW_eeu_FakeMu_80X"  ;
          //hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "checkMTW_eee_FakeEl_80X"  ;
          // hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->SetTitle(output_histo_name.c_str());
          hist_check_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "checkMTW_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkMTW_uuu_"  + dataSetName ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          //hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkMTW_uue_"  + dataSetName ;
          // hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->SetTitle(output_histo_name.c_str());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkMTW_eeu_"  + dataSetName ;
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          //hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkMTW_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkMTW_eee_"  + dataSetName ;
          hist_check_eee->SetTitle(output_histo_name.c_str());
          // hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->Write(output_histo_name.c_str());
        }
      }
      else if(doTTZtemplate){
        if (dataSetName.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = "checkZmass_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "checkZmass_uuu_FakeMu_80X"  ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          //hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = "checkZmass_uue_FakeEl_80X"  ;
          //hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->SetTitle(output_histo_name.c_str());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "checkZmass_eeu_FakeMu_80X"  ;
          //hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "checkZmass_eee_FakeEl_80X"  ;
          // hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->SetTitle(output_histo_name.c_str());
          hist_check_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "checkZmass_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkZmass_uuu_"  + dataSetName ;
          hist_check_uuu->SetTitle(output_histo_name.c_str());
          //hist_check_uuu->Scale(1/hist_check_uuu->Integral());
          hist_check_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkZmass_uue_"  + dataSetName ;
          // hist_check_uue->Scale(1/hist_check_uue->Integral());
          hist_check_uue->SetTitle(output_histo_name.c_str());
          hist_check_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkZmass_eeu_"  + dataSetName ;
          hist_check_eeu->SetTitle(output_histo_name.c_str());
          //hist_check_eeu->Scale(1/hist_check_eeu->Integral());
          hist_check_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "checkZmass_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "checkZmass_eee_"  + dataSetName ;
          hist_check_eee->SetTitle(output_histo_name.c_str());
          // hist_check_eee->Scale(1/hist_check_eee->Integral());
          hist_check_eee->Write(output_histo_name.c_str());
        }
      }
      
      /*
      if(isys == 0 && dataSetName.find("fake")!=std::string::npos){
        
        cout << endl;
        cout << "************************* KOLMOGOROV TESTING : fakes ************************ " << endl;
        cout << hist_uuu->KolmogorovTest(hist_check_uuu,"D")<< endl;
        cout << hist_uue->KolmogorovTest(hist_check_uue,"D")<< endl;
        cout << hist_eee->KolmogorovTest(hist_check_eee,"D")<< endl;
        cout << hist_eeu->KolmogorovTest(hist_check_eeu,"D")<< endl;
        cout << "************************* KOLMOGOROV TESTING : fakes ************************ " << endl;
        cout << endl;
        
      }
      */
      
      
      /*if(dataSetName.find("fake")!=std::string::npos){
       TCanvas* c1 = new TCanvas();
       gStyle->SetOptStat(0);
       c1->cd();
       hist_check_uuu->SetLineColor(kRed);
       hist_check_uuu->Draw("hist");
       hist_uuu->Draw("hist sames");
       c1->Modified();
       c1->SaveAs("checkMTW_uuu_FakeMu_80X.png");
       
       c1 = new TCanvas();
       c1->cd();
       hist_check_uue->SetLineColor(kRed);
       hist_check_uue->Draw("hist");
       hist_uue->Draw("hist sames");
       c1->Modified();
       c1->SaveAs("checkMTW_uue_FakeEl_80X.png");
       
       c1 = new TCanvas();
       c1->cd();
       hist_check_eeu->SetLineColor(kRed);
       hist_check_eeu->Draw("hist");
       hist_eeu->Draw("hist sames");
       c1->Modified();
       c1->SaveAs("checkMTW_eeu_FakeMu_80X.png");
       
       c1 = new TCanvas();
       c1->cd();
       hist_check_eee->SetLineColor(kRed);
       hist_check_eee->Draw("hist");
       hist_eee->Draw("hist sames");
       c1->Modified();
       c1->SaveAs("checkMTW_eee_FakeEl_80X.png");
       
       }*/
      combinetemplate_file->Close();
      
      //cout << "closed " << combinetemplate_filename.c_str() << endl;
      delete combinetemplate_file;
      delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
      delete  hist_check_eee; delete hist_check_uuu; delete hist_check_uue; delete hist_check_eeu;
     // if(doMTWtemplate || doTTZtemplate) tFileMap[dataSetName.c_str()]->Close();
    }// datasets
    
    /*if(isys == 0){
      cout << endl;
      cout << "************************* KOLMOGOROV TESTING : Non prompts ************************ " << endl;
      cout << hist_BDT_tt_nonpromptinZ->KolmogorovTest(hist_BDT_tt_nonpromptinW,"D") << endl;
      cout << hist_BDT_tt_uuu_nonpromptinZ->KolmogorovTest(hist_BDT_tt_uuu_nonpromptinW,"D")<< endl;
      cout << hist_BDT_tt_uue_nonpromptinZ->KolmogorovTest(hist_BDT_tt_uue_nonpromptinW,"D")<< endl;
      cout << hist_BDT_tt_eee_nonpromptinZ->KolmogorovTest(hist_BDT_tt_eee_nonpromptinW,"D")<< endl;
      cout << hist_BDT_tt_eeu_nonpromptinZ->KolmogorovTest(hist_BDT_tt_eeu_nonpromptinW,"D")<< endl;
      cout << "************************* KOLMOGOROV TESTING : Non prompts ************************ " << endl;
      cout << endl;
      
      cout << endl;
      cout << "************************* KOLMOGOROV TESTING : light vs c ************************ " << endl;
      cout << hist_BDT_WZ_light->KolmogorovTest(hist_BDT_WZ_c,"D") << endl;
      cout << hist_BDT_WZ_uuu_light->KolmogorovTest(hist_BDT_WZ_uuu_c,"D")<< endl;
      cout << hist_BDT_WZ_uue_light->KolmogorovTest(hist_BDT_WZ_uue_c,"D")<< endl;
      cout << hist_BDT_WZ_eee_light->KolmogorovTest(hist_BDT_WZ_eee_c,"D")<< endl;
      cout << hist_BDT_WZ_eeu_light->KolmogorovTest(hist_BDT_WZ_eeu_c,"D")<< endl;
      cout << "************************* KOLMOGOROV TESTING : light vs c ************************ " << endl;
      cout << endl;
     
      
      cout << endl;
      cout << "************************* KOLMOGOROV TESTING : light vs b ************************ " << endl;
      cout << hist_BDT_WZ_light->KolmogorovTest(hist_BDT_WZ_c,"D") << endl;
      cout << hist_BDT_WZ_uuu_light->KolmogorovTest(hist_BDT_WZ_uuu_b,"D")<< endl;
      cout << hist_BDT_WZ_uue_light->KolmogorovTest(hist_BDT_WZ_uue_b,"D")<< endl;
      cout << hist_BDT_WZ_eee_light->KolmogorovTest(hist_BDT_WZ_eee_b,"D")<< endl;
      cout << hist_BDT_WZ_eeu_light->KolmogorovTest(hist_BDT_WZ_eeu_b,"D")<< endl;
      cout << "************************* KOLMOGOROV TESTING : light vs b ************************ " << endl;
      cout << endl;
      
      cout << endl;
      cout << "************************* KOLMOGOROV TESTING : b vs c ************************ " << endl;
      cout << hist_BDT_WZ_b->KolmogorovTest(hist_BDT_WZ_c,"D") << endl;
      cout << hist_BDT_WZ_uuu_b->KolmogorovTest(hist_BDT_WZ_uuu_c,"D")<< endl;
      cout << hist_BDT_WZ_uue_b->KolmogorovTest(hist_BDT_WZ_uue_c,"D")<< endl;
      cout << hist_BDT_WZ_eee_b->KolmogorovTest(hist_BDT_WZ_eee_c,"D")<< endl;
      cout << hist_BDT_WZ_eeu_b->KolmogorovTest(hist_BDT_WZ_eeu_c,"D")<< endl;
      cout << "************************* KOLMOGOROV TESTING : b vs c ************************ " << endl;
      cout << endl;
    }*/
    
    
    
    if(isys != 0) cout<<"Done with "<< systematic <<" systematic"<<endl;
    else cout<<"Done with nominal sample"<<endl;
    systematic = "";
  } // systematics
  
  //cout << "ENTRIES " << hist_WZ->GetEntries() << endl;
  
  cout << "ENTRIES " << histo1DMTW["MTW_WZ"]->GetEntries() << endl;
  if(!doMTWtemplate && !doTTZtemplate && !doPostFit){
    
   // delete fin;
  }
  //if(doTTZtemplate) fin->Close();
  ///*****************///
  ///   PDF envelope   ///
  ///*****************///
  
  
  if(doPDFunc ) {
    GetPDFEnvelope("WZTo3LNu_amc_80X", doMTWtemplate,doTTZtemplate);
    GetPDFEnvelope("tZq_amc_80X", doMTWtemplate, doTTZtemplate);
    GetPDFEnvelope("TTZToLLNuNu_amc_80X", doMTWtemplate,doTTZtemplate);
    GetPDFEnvelope("ZZTo4L_80X",doMTWtemplate, doTTZtemplate);
    if(!doMTWtemplate && !doTTZtemplate){
   if(doZut ){
     GetPDFEnvelope("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X",doMTWtemplate, doTTZtemplate);
     GetPDFEnvelope("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X",doMTWtemplate, doTTZtemplate);
     GetPDFEnvelope("NP_overlay_ST_FCNC_zut_80X",doMTWtemplate, doTTZtemplate);
     
   }
   else {GetPDFEnvelope("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X",doMTWtemplate, doTTZtemplate);
     GetPDFEnvelope("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X",doMTWtemplate, doTTZtemplate);
     GetPDFEnvelope("NP_overlay_ST_FCNC_zct_80X",doMTWtemplate, doTTZtemplate);
   }}
    else {
      GetPDFEnvelope("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X",doMTWtemplate, doTTZtemplate);
      GetPDFEnvelope("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X",doMTWtemplate, doTTZtemplate);
      GetPDFEnvelope("NP_overlay_ST_FCNC_zut_80X",doMTWtemplate, doTTZtemplate);
      GetPDFEnvelope("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X",doMTWtemplate, doTTZtemplate);
      GetPDFEnvelope("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X",doMTWtemplate, doTTZtemplate);
      GetPDFEnvelope("NP_overlay_ST_FCNC_zct_80X",doMTWtemplate, doTTZtemplate);


    }
  }
 
  
  ///*****************///
  ///   Pseudodata   ///
  ///*****************///
  if(doPseudoData && !addData){
    cout << "generating pseudo data" << endl;
    TRandom3 therand(0); //Randomization
    
    string pseudodata_input_name = placeOutputReading+"/Reader_" + coupling + "_" + region + ".root";
    
    if(doMTWtemplate ) pseudodata_input_name = placeOutputReading+"/Reader_"+coupling+"_MTW.root";
    TFile* pseudodata_file = TFile::Open( pseudodata_input_name.c_str(), "UPDATE" );
    
    cout << "Generating pseudo data from " << pseudodata_input_name << endl;
    
    
    
    TH1F *h_sum = 0, *h_tmp = 0;
    string histo_name = "";
    string template_fake_name = "";
    
    
    std::vector<string> channel_list;
    channel_list.push_back("eee");
    channel_list.push_back("uue");
    channel_list.push_back("eeu");
    channel_list.push_back("uuu");
    for(Int_t ichan=0; ichan<channel_list.size(); ichan++)
    {
      h_sum = 0;
      histo_name = "";
      if(channel_list[ichan] == "uuu" || channel_list[ichan] == "eeu") {template_fake_name = "FakeMu_80X";}
      else {template_fake_name = "FakeEl_80X";}
      
      
      
      
      for(Int_t isample = 0; isample < datasets.size(); isample++)
      {
        string dataSetName = datasets[isample]->Name();
        // cout << dataSetName << endl;
        if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {continue; } // no signal in data
        else if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // safety
        else if(datasets[isample]->Name().find("fake")==std::string::npos){
          // cout << "  -- sample " << datasets[isample]->Name() << endl;
          h_tmp = 0;
          if(!doMTWtemplate) histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
          else histo_name = "MTW_" + channel_list[ichan] + "_" + datasets[isample]->Name();
          //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
          //  cout << "  --- histo " << histo_name << endl;
          if(!pseudodata_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
          h_tmp = (TH1F*) pseudodata_file->Get(histo_name.c_str());
          if(h_sum == 0) {h_sum = (TH1F*) h_tmp->Clone();}
          else {h_sum->Add(h_tmp);}
        }
        else{
          
          if(!doMTWtemplate) histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + template_fake_name;
          else if(doMTWtemplate) histo_name = "MTW_" + channel_list[ichan] + "_" + template_fake_name;
          // cout << "  --- histo " << histo_name << endl;
          if(!pseudodata_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<histo_name<<" : not found"<<endl;}
          else
          {
            h_tmp = (TH1F*) pseudodata_file->Get(histo_name.c_str());
            if(h_sum == 0) {h_sum = (TH1F*) h_tmp->Clone();}
            else {h_sum->Add(h_tmp);}
          }
          
        }
      }
      
      if(h_sum == 0) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; }
      Int_t nofbins = h_sum->GetNbinsX();
      
      for(Int_t i=0; i<nofbins; i++)
      {
        Double_t bin_content = h_sum->GetBinContent(i+1); //cout<<"bin "<<i+1<<endl; cout<<"initial content = "<<bin_content<<endl;
        Double_t new_bin_content = therand.Poisson(bin_content);// cout<<"new content = "<<new_bin_content<<endl;
        h_sum->SetBinContent(i+1, new_bin_content);
        h_sum->SetBinError(i+1, sqrt(new_bin_content)); //Poissonian error
      }
      
      pseudodata_file->cd();
      string output_histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_data";
      if(doMTWtemplate) output_histo_name = "MTW_" + channel_list[ichan] + "_data";
      h_sum->SetTitle(output_histo_name.c_str());
      h_sum->SetName(output_histo_name.c_str());
      h_sum->Write(output_histo_name.c_str(), TObject::kOverwrite);
      
    } // chan
    
    pseudodata_file->Close();
    
    cout<<"--- Done with generation of pseudo-data"<<endl<<endl;
    
    delete pseudodata_file;
    
  }
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  if((makePlots || doPDFunc || PlotMVAvars || PlotSystematics || PlotJeSystematics) ){
    string pathOutput = "OutputPlots/";
    mkdir(pathOutput.c_str(),0777);
    string pathOutputdate = pathOutput + dateString + "/"  ;
    mkdir(pathOutputdate.c_str(),0777);
    string place =pathOutputdate+"MSPlot/";
    if(doMTWtemplate) place = pathOutputdate + "MSPlotMTW/";
    if(doTTZtemplate) place = pathOutputdate + "MSPlotTTZ/";
    string placeTH1F = pathOutputdate+"TH1F/";
    string placeTH2F = pathOutputdate+"TH2F/";
    vector <string> vlabel_chan = {"3#mu", "1e2#mu", "2e1#mu", "3e"};
    mkdir(place.c_str(),0777);
    mkdir(placeTH1F.c_str(),0777);
    mkdir(placeTH2F.c_str(),0777);
    string rootFileName ="NtuplePlotsMVA.root";
    if(doMTWtemplate) rootFileName = "NtuplePlotsMTW.root";
    if(doTTZtemplate) rootFileName = "NtuplePlotsTTZ.root";
    cout << " - Recreate output file ..." << endl;
    
    TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutputdate+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    
    if(makePlots && (doMTWtemplate )){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlotMTW.begin(); it != MSPlotMTW.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        if(!datafound){
          cout << "no data found, setting lumi as " << Luminosity << endl;
          temp->setDataLumi(Luminosity);
        }
        //temp->SetPreliminary(false);
        double scalefakes =  1.;
        if(!datafound) temp->setDataLumi(Luminosity);
        if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
        if(name.find("eee")!=std::string::npos){ temp->setChannel(true, "3e"); scalefakes = scaleFakes_eee; cout << "scalefakes_eee " << scaleFakes_eee << endl; }
        if(name.find("eeu")!=std::string::npos){ temp->setChannel(true, "2e1#mu"); scalefakes = scaleFakes_eeu;  cout << "scalefakes_eeu " << scaleFakes_eeu << endl;  }
        if(name.find("uue")!=std::string::npos){ temp->setChannel(true, "1e2#mu"); scalefakes = scaleFakes_uue;  cout << "scalefakes_uue " << scaleFakes_uue << endl; }
        if(name.find("uuu")!=std::string::npos){ temp->setChannel(true, "3#mu"); scalefakes = scaleFakes_uuu;  cout << "scalefakes_uuu " << scaleFakes_uuu<< endl;  }
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        //temp->SetPreliminary(false);
        temp->Draw(name, 1, false, false, false, 10,scalefakes);  // string label, unsigned Int_t RatioType, bool addRatioErrorBand, bool addErrorBand, bool    cout << "writing to " << pathOutputdate+"MSPlotMTW" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        if(doMTWtemplate) temp->Write(fout, name, true, (pathOutputdate+"MSPlotMTW").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
        
        
      }
      
    }
    if(makePlots && doTTZtemplate){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlotTTZ.begin(); it != MSPlotTTZ.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        if(!datafound){
          cout << "no data found, setting lumi as " << Luminosity << endl;
          temp->setDataLumi(Luminosity);
        }
        //temp->SetPreliminary(false);
        double scalefakes =  1.;
        if(!datafound) temp->setDataLumi(Luminosity);
        if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
        if(name.find("eee")!=std::string::npos){ temp->setChannel(true, "3e"); scalefakes = scaleFakes_eee; }
        if(name.find("eeu")!=std::string::npos){ temp->setChannel(true, "2e1#mu"); scalefakes = scaleFakes_eeu; }
        if(name.find("uue")!=std::string::npos){ temp->setChannel(true, "1e2#mu"); scalefakes = scaleFakes_uue; }
        if(name.find("uuu")!=std::string::npos){ temp->setChannel(true, "3#mu"); scalefakes = scaleFakes_uuu; }
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        //temp->SetPreliminary(false);
        temp->Draw(name, 1, false, false, false, 10,scalefakes);  // string label, unsigned Int_t RatioType, bool addRatioErrorBand, bool addErrorBand, bool    cout << "writing to " << pathOutputdate+"MSPlotMTW" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        if(doTTZtemplate) temp->Write(fout, name, true, (pathOutputdate+"MSPlotTTZ").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
        
        
      }
      
    }
    if(makePlots && !doMTWtemplate && !doTTZtemplate){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
        cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        double scalefakes = 1.;
        if(!datafound) temp->setDataLumi(Luminosity);
        if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
          if(name.find("eee")!=std::string::npos){ temp->setChannel(true, "3e"); scalefakes = scaleFakes_eee; cout << "scalefakes_eee " << scaleFakes_eee << endl;  }
          if(name.find("eeu")!=std::string::npos){ temp->setChannel(true, "2e1#mu"); scalefakes = scaleFakes_eeu; cout << "scalefakes_eeu " << scaleFakes_eeu << endl;  }
          if(name.find("uue")!=std::string::npos){ temp->setChannel(true, "1e2#mu"); scalefakes = scaleFakes_uue;cout << "scalefakes_uue " << scaleFakes_uue << endl;  }
          if(name.find("uuu")!=std::string::npos){ temp->setChannel(true, "3#mu"); scalefakes = scaleFakes_uuu; cout << "scalefakes_uuu " << scaleFakes_uuu << endl;  }
        if(doPostFit){
        if(name.find("eee")!=std::string::npos){ temp->setChannel(true, "3e"); scalefakes =  1.; cout << "scalefakes_eee " << scaleFakes_eee << endl;  }
        if(name.find("eeu")!=std::string::npos){ temp->setChannel(true, "2e1#mu"); scalefakes = 1.; cout << "scalefakes_eeu " << scaleFakes_eeu << endl;  }
        if(name.find("uue")!=std::string::npos){ temp->setChannel(true, "1e2#mu"); scalefakes = 1.;cout << "scalefakes_uue " << scaleFakes_uue << endl;  }
        if(name.find("uuu")!=std::string::npos){ temp->setChannel(true, "3#mu"); scalefakes =  1.; cout << "scalefakes_uuu " << scaleFakes_uuu << endl;  }
        
        }
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        //temp->SetPreliminary(false);
        if(doPostFit && doErrorband){
          //temp->showNumberEntries(false);
          if(doZut ) temp->setErrorBandFile("Zut_PostFitErrorBand.root");
          else if(!doZut ) temp->setErrorBandFile("Zct_PostFitErrorBand.root");
          temp->Draw(name, 1, true, true, true, 10,scalefakes);
        }
        else temp->Draw(name, 1, false, false, false, 10,scalefakes);  // string label, unsigned Int_t RatioType, bool addRatioErrorBand, bool addErrorBand, bool    cout << "writing to " << cout << "writing to " << pathOutputdate+"MSPlot" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        temp->Write(fout, name, true, (pathOutputdate+"MSPlot").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
    }
    if(doPDFunc ){
      string pdfrootFileName ="PDFhistograms_BDT_toppair_Zutttz.root";
      if(!doMTWtemplate && !doTTZtemplate){
       if(toppair && doZut) pdfrootFileName ="PDFhistograms_BDT_toppair_Zutttz.root";
       else if(!toppair && doZut) pdfrootFileName ="PDFhistograms_BDT_singletop_Zutttz.root";
       else if(toppair && !doZut) pdfrootFileName ="PDFhistograms_BDT_toppair_Zctttz.root";
       else if(!toppair && !doZut) pdfrootFileName ="PDFhistograms_BDT_singletop_Zctttz.root";
      }
      else if(!doMTWtemplate ){
        if(toppair && doZut) pdfrootFileName ="PDFhistograms_MTW_toppair_Zutttz.root";
        else if(!toppair && doZut) pdfrootFileName ="PDFhistograms_MTW_singletop_Zutttz.root";
        else if(toppair && !doZut) pdfrootFileName ="PDFhistograms_MTW_toppair_Zctttz.root";
        else if(!toppair && !doZut) pdfrootFileName ="PDFhistograms_MTW_singletop_Zcttz.root";
      }
      else if( doTTZtemplate){
        if(toppair && doZut) pdfrootFileName ="PDFhistograms_Zmass_toppair_Zutttz.root";
        else if(!toppair && doZut) pdfrootFileName ="PDFhistograms_Zmass_singletop_Zutttz.root";
        else if(toppair && !doZut) pdfrootFileName ="PDFhistograms_Zmass_toppair_Zctttz.root";
        else if(!toppair && !doZut) pdfrootFileName ="PDFhistograms_Zmass_singletop_Zctttz.root";
      }
      TFile *pdffout = new TFile ((pdfrootFileName).c_str(), "RECREATE");
      cout << "   Output file is " << pdfrootFileName << endl;
      pdffout->cd();
      TDirectory* th1dir = pdffout->mkdir("1D_PDF_histograms");
      th1dir->cd();
      gStyle->SetOptStat(1110);
      if(!doMTWtemplate ){
      for (std::map<std::string,TH1F*>::const_iterator it = histo1DPDF.begin(); it != histo1DPDF.end(); it++)
      {
        TH1F *temp = it->second;
        /*Int_t N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...*/
        temp->Write();
        if(it->first.find("PDFEnvelope")==std::string::npos) continue; 
        if(it->first.find("nominal")!=std::string::npos || it->first.find("Up")!=std::string::npos || it->first.find("Down")!=std::string::npos) {
          TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
          tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
        }
      }
      }
      else if(doMTWtemplate ){
        for (std::map<std::string,TH1F*>::const_iterator it = histo1DPDFMTW.begin(); it != histo1DPDFMTW.end(); it++)
        {
          TH1F *temp = it->second;
          /*Int_t N = temp->GetNbinsX();
           temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
           temp->SetBinContent(N+1,0);
           temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...*/
          temp->Write();
          if(it->first.find("PDFEnvelope")==std::string::npos) continue;
          if(it->first.find("nominal")!=std::string::npos || it->first.find("Up")!=std::string::npos || it->first.find("Down")!=std::string::npos) {
            TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
            tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
          }
        }
      }
      pdffout->Write();
      pdffout->Close();
      delete pdffout;
    }
    if(PlotJeSystematics && !doMTWtemplate && !doTTZtemplate){
      double maximum_sig = hist_BDT_JES_nom_sig->GetMaximum()*1.5;
      double maximum_bkg = hist_BDT_JES_nom_bkg->GetMaximum()*1.5;
      
      hist_BDT_JES_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_JES_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_JER_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_JER_nom_sig->SetMaximum(maximum_sig);
      
      hist_BDT_JER_nom_sig->SetLineColor(kRed);
      hist_BDT_JER_up_sig->SetLineColor(kBlue);
      hist_BDT_JER_down_sig->SetLineColor(kGreen+2);
      hist_BDT_JER_nom_bkg->SetLineColor(kRed);
      hist_BDT_JER_up_bkg->SetLineColor(kBlue);
      hist_BDT_JER_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_JER_nom_sig->SetLineWidth(2);
      hist_BDT_JER_up_sig->SetLineWidth(2);
      hist_BDT_JER_down_sig->SetLineWidth(2);
      hist_BDT_JER_nom_bkg->SetLineWidth(2);
      hist_BDT_JER_up_bkg->SetLineWidth(2);
      hist_BDT_JER_down_bkg->SetLineWidth(2);
      
      hist_BDT_JES_nom_sig->SetLineColor(kRed);
      hist_BDT_JES_up_sig->SetLineColor(kBlue);
      hist_BDT_JES_down_sig->SetLineColor(kGreen+2);
      hist_BDT_JES_nom_bkg->SetLineColor(kRed);
      hist_BDT_JES_up_bkg->SetLineColor(kBlue);
      hist_BDT_JES_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_JES_nom_sig->SetLineWidth(2);
      hist_BDT_JES_up_sig->SetLineWidth(2);
      hist_BDT_JES_down_sig->SetLineWidth(2);
      hist_BDT_JES_nom_bkg->SetLineWidth(2);
      hist_BDT_JES_up_bkg->SetLineWidth(2);
      hist_BDT_JES_down_bkg->SetLineWidth(2);
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *legsig = new TLegend(xl1,yl1,xl2,yl2);
      legsig->AddEntry(hist_BDT_JER_nom_sig,"nominal","L");   // h1 and h2 are histogram pointers
      legsig->AddEntry(hist_BDT_JER_up_sig,"unc +","L");
      legsig->AddEntry(hist_BDT_JER_down_sig,"unc -","L");
      
      gStyle->SetOptStat(0);
      
      //JER
      TCanvas* tempCanvas = new TCanvas("", "");
      TPad *histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_JER_nom_sig->Draw("e histo");
      hist_BDT_JER_up_sig->Draw("e same histo");
      hist_BDT_JER_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      TPad *ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      TH1F* ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_JER_nom_sig->GetNbinsX(),hist_BDT_JER_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_JER_nom_sig,hist_BDT_JER_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      TH1F* ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_JER_nom_sig->GetNbinsX(),hist_BDT_JER_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_JER_nom_sig,hist_BDT_JER_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      TLine *line = new TLine(hist_BDT_JER_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_JER_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_JER_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_JER_nom_bkg->Draw("e histo");
      hist_BDT_JER_up_bkg->Draw("e same histo");
      hist_BDT_JER_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_JER_nom_bkg->GetNbinsX(),hist_BDT_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_JER_nom_bkg,hist_BDT_JER_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_JER_nom_bkg->GetNbinsX(),hist_BDT_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_JER_nom_bkg,hist_BDT_JER_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JER_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_JER_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_JER_nom_bkg_LogY.png");
      
      
      
      //JES
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_JES_nom_sig->Draw("e histo");
      hist_BDT_JES_up_sig->Draw("e same histo");
      hist_BDT_JES_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_JES_nom_sig->GetNbinsX(),hist_BDT_JES_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_JES_nom_sig,hist_BDT_JES_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_JES_nom_sig->GetNbinsX(),hist_BDT_JES_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_JES_nom_sig,hist_BDT_JES_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_JES_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_JES_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_JES_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_JES_nom_bkg->Draw("e histo");
      hist_BDT_JES_up_bkg->Draw("e same histo");
      hist_BDT_JES_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_JES_nom_bkg->GetNbinsX(),hist_BDT_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_JES_nom_bkg,hist_BDT_JES_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_JES_nom_bkg->GetNbinsX(),hist_BDT_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_JES_nom_bkg,hist_BDT_JES_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_JES_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_JES_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_JES_nom_bkg_LogY.png");
      
      ////// CUTFLOWS
      string sysratiosfilename = "sysratio_" ;
      if(doZut && toppair) sysratiosfilename += "TTZut";
      else if(!doZut && toppair) sysratiosfilename += "TTZct";
      else if(doZut && !toppair) sysratiosfilename += "STZut";
      else if(!doZut && !toppair) sysratiosfilename += "STZct";
      TFile* sysratiosfile = TFile::Open((sysratiosfilename+"_jerjes.root").c_str(),"RECREATE");
      sysratiosfile->cd();
      cout << "opening" << (sysratiosfilename+"_jerjes.root").c_str() << endl;
      hist_BDT_JER_nom_sig->Write();
      hist_BDT_JER_up_sig->Write();
      hist_BDT_JER_down_sig->Write();
      hist_BDT_JER_nom_bkg->Write();
      hist_BDT_JER_up_bkg->Write();
      hist_BDT_JER_down_bkg->Write();
      hist_BDT_JES_nom_sig->Write();
      hist_BDT_JES_up_sig->Write();
      hist_BDT_JES_down_sig->Write();
      hist_BDT_JES_nom_bkg->Write();
      hist_BDT_JES_up_bkg->Write();
      hist_BDT_JES_down_bkg->Write();
      
      sysratiosfile->Close();
      
      
    }
    if(PlotJeSystematics && (doMTWtemplate)){
      double maximum_sig = hist_mWt_JES_nom_sig->GetMaximum()*1.5;
      double maximum_bkg = hist_mWt_JES_nom_bkg->GetMaximum()*1.5;
      
      hist_mWt_JES_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_JES_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_JER_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_JER_nom_sig->SetMaximum(maximum_sig);
      
      hist_mWt_JER_nom_sig->SetLineColor(kRed);
      hist_mWt_JER_up_sig->SetLineColor(kBlue);
      hist_mWt_JER_down_sig->SetLineColor(kGreen+2);
      hist_mWt_JER_nom_bkg->SetLineColor(kRed);
      hist_mWt_JER_up_bkg->SetLineColor(kBlue);
      hist_mWt_JER_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_JER_nom_sig->SetLineWidth(2);
      hist_mWt_JER_up_sig->SetLineWidth(2);
      hist_mWt_JER_down_sig->SetLineWidth(2);
      hist_mWt_JER_nom_bkg->SetLineWidth(2);
      hist_mWt_JER_up_bkg->SetLineWidth(2);
      hist_mWt_JER_down_bkg->SetLineWidth(2);
      
      hist_mWt_JES_nom_sig->SetLineColor(kRed);
      hist_mWt_JES_up_sig->SetLineColor(kBlue);
      hist_mWt_JES_down_sig->SetLineColor(kGreen+2);
      hist_mWt_JES_nom_bkg->SetLineColor(kRed);
      hist_mWt_JES_up_bkg->SetLineColor(kBlue);
      hist_mWt_JES_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_JES_nom_sig->SetLineWidth(2);
      hist_mWt_JES_up_sig->SetLineWidth(2);
      hist_mWt_JES_down_sig->SetLineWidth(2);
      hist_mWt_JES_nom_bkg->SetLineWidth(2);
      hist_mWt_JES_up_bkg->SetLineWidth(2);
      hist_mWt_JES_down_bkg->SetLineWidth(2);
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *legsig = new TLegend(xl1,yl1,xl2,yl2);
      legsig->AddEntry(hist_mWt_JER_nom_sig,"nominal","L");   // h1 and h2 are histogram pointers
      legsig->AddEntry(hist_mWt_JER_up_sig,"unc +","L");
      legsig->AddEntry(hist_mWt_JER_down_sig,"unc -","L");
      
      gStyle->SetOptStat(0);
      
      //JER
      TCanvas* tempCanvas = new TCanvas("", "");
      TPad *histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_JER_nom_sig->Draw("e histo");
      hist_mWt_JER_up_sig->Draw("e same histo");
      hist_mWt_JER_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      TPad *ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      TH1F* ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_JER_nom_sig->GetNbinsX(),hist_mWt_JER_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_JER_nom_sig,hist_mWt_JER_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      TH1F* ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_JER_nom_sig->GetNbinsX(),hist_mWt_JER_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_JER_nom_sig,hist_mWt_JER_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      TLine *line = new TLine(hist_mWt_JER_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_JER_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_JER_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_JER_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_JER_nom_bkg->Draw("e histo");
      hist_mWt_JER_up_bkg->Draw("e same histo");
      hist_mWt_JER_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_JER_nom_bkg->GetNbinsX(),hist_mWt_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_JER_nom_bkg,hist_mWt_JER_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_JER_nom_bkg->GetNbinsX(),hist_mWt_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_JER_nom_bkg,hist_mWt_JER_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_JER_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_JER_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JER_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_JER_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_JER_nom_bkg_LogY.png");
      
      
      
      //JES
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_JES_nom_sig->Draw("e histo");
      hist_mWt_JES_up_sig->Draw("e same histo");
      hist_mWt_JES_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_JES_nom_sig->GetNbinsX(),hist_mWt_JES_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_JES_nom_sig,hist_mWt_JES_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_JES_nom_sig->GetNbinsX(),hist_mWt_JES_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_JES_nom_sig,hist_mWt_JES_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_JES_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_JES_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_JES_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_JES_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_JES_nom_bkg->Draw("e histo");
      hist_mWt_JES_up_bkg->Draw("e same histo");
      hist_mWt_JES_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_JES_nom_bkg->GetNbinsX(),hist_mWt_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_JES_nom_bkg,hist_mWt_JES_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_JES_nom_bkg->GetNbinsX(),hist_mWt_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_JES_nom_bkg,hist_mWt_JES_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_JES_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_JES_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_JES_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_JES_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_JES_nom_bkg_LogY.png");
      
      
    }
    
    if(PlotSystematics && !doMTWtemplate && !doTTZtemplate){
      double maximum_sig = hist_BDT_puSF_nom_sig->GetMaximum()*1.5;
      double maximum_bkg = hist_BDT_puSF_nom_bkg->GetMaximum()*1.5;
      
      hist_BDT_puSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_puSF_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_electronSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_electronSF_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_muonSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_muonSF_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_cferr1_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_cferr1_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_cferr2_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_cferr2_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_hf_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_hf_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_hfstats1_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_hfstats1_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_hfstats2_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_hfstats2_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_lf_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_lf_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_lfstats1_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_lfstats1_nom_sig->SetMaximum(maximum_sig);
      hist_BDT_btagSF_lfstats2_nom_bkg->SetMaximum(maximum_bkg);
      hist_BDT_btagSF_lfstats2_nom_sig->SetMaximum(maximum_sig);
      
      hist_BDT_puSF_nom_sig->SetLineColor(kRed);
      hist_BDT_puSF_up_sig->SetLineColor(kBlue);
      hist_BDT_puSF_down_sig->SetLineColor(kGreen+2);
      hist_BDT_puSF_nom_bkg->SetLineColor(kRed);
      hist_BDT_puSF_up_bkg->SetLineColor(kBlue);
      hist_BDT_puSF_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_puSF_nom_sig->SetLineWidth(2);
      hist_BDT_puSF_up_sig->SetLineWidth(2);
      hist_BDT_puSF_down_sig->SetLineWidth(2);
      hist_BDT_puSF_nom_bkg->SetLineWidth(2);
      hist_BDT_puSF_up_bkg->SetLineWidth(2);
      hist_BDT_puSF_down_bkg->SetLineWidth(2);
      
      
      hist_BDT_electronSF_nom_sig->SetLineColor(kRed);
      hist_BDT_electronSF_up_sig->SetLineColor(kBlue);
      hist_BDT_electronSF_down_sig->SetLineColor(kGreen+2);
      hist_BDT_electronSF_nom_bkg->SetLineColor(kRed);
      hist_BDT_electronSF_up_bkg->SetLineColor(kBlue);
      hist_BDT_electronSF_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_electronSF_nom_sig->SetLineWidth(2);
      hist_BDT_electronSF_up_sig->SetLineWidth(2);
      hist_BDT_electronSF_down_sig->SetLineWidth(2);
      hist_BDT_electronSF_nom_bkg->SetLineWidth(2);
      hist_BDT_electronSF_up_bkg->SetLineWidth(2);
      hist_BDT_electronSF_down_bkg->SetLineWidth(2);
      
      hist_BDT_muonSF_nom_sig->SetLineColor(kRed);
      hist_BDT_muonSF_up_sig->SetLineColor(kBlue);
      hist_BDT_muonSF_down_sig->SetLineColor(kGreen+2);
      hist_BDT_muonSF_nom_bkg->SetLineColor(kRed);
      hist_BDT_muonSF_up_bkg->SetLineColor(kBlue);
      hist_BDT_muonSF_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_muonSF_nom_sig->SetLineWidth(2);
      hist_BDT_muonSF_up_sig->SetLineWidth(2);
      hist_BDT_muonSF_down_sig->SetLineWidth(2);
      hist_BDT_muonSF_nom_bkg->SetLineWidth(2);
      hist_BDT_muonSF_up_bkg->SetLineWidth(2);
      hist_BDT_muonSF_down_bkg->SetLineWidth(2);
      
      
      
      hist_BDT_btagSF_cferr1_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_cferr1_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_cferr1_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_cferr1_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_cferr1_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_cferr1_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_cferr1_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr1_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr1_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr1_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_cferr1_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_cferr1_down_bkg->SetLineWidth(2);
      
      
      hist_BDT_btagSF_cferr2_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_cferr2_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_cferr2_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_cferr2_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_cferr2_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_cferr2_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_cferr2_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr2_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr2_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_cferr2_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_cferr2_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_cferr2_down_bkg->SetLineWidth(2);
      
      hist_BDT_btagSF_hf_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_hf_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_hf_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hf_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_hf_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_hf_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hf_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_hf_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_hf_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_hf_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hf_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hf_down_bkg->SetLineWidth(2);
      
      hist_BDT_btagSF_hfstats1_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_hfstats1_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_hfstats1_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hfstats1_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_hfstats1_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_hfstats1_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hfstats1_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats1_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats1_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats1_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hfstats1_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hfstats1_down_bkg->SetLineWidth(2);
      
      
      hist_BDT_btagSF_hfstats2_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_hfstats2_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_hfstats2_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hfstats2_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_hfstats2_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_hfstats2_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_hfstats2_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats2_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats2_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_hfstats2_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hfstats2_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_hfstats2_down_bkg->SetLineWidth(2);
      
      hist_BDT_btagSF_lf_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_lf_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_lf_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lf_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_lf_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_lf_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lf_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_lf_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_lf_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_lf_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lf_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lf_down_bkg->SetLineWidth(2);
      
      hist_BDT_btagSF_lfstats1_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_lfstats1_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_lfstats1_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lfstats1_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_lfstats1_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_lfstats1_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lfstats1_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats1_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats1_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats1_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lfstats1_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lfstats1_down_bkg->SetLineWidth(2);
      
      
      hist_BDT_btagSF_lfstats2_nom_sig->SetLineColor(kRed);
      hist_BDT_btagSF_lfstats2_up_sig->SetLineColor(kBlue);
      hist_BDT_btagSF_lfstats2_down_sig->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lfstats2_nom_bkg->SetLineColor(kRed);
      hist_BDT_btagSF_lfstats2_up_bkg->SetLineColor(kBlue);
      hist_BDT_btagSF_lfstats2_down_bkg->SetLineColor(kGreen+2);
      hist_BDT_btagSF_lfstats2_nom_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats2_up_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats2_down_sig->SetLineWidth(2);
      hist_BDT_btagSF_lfstats2_nom_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lfstats2_up_bkg->SetLineWidth(2);
      hist_BDT_btagSF_lfstats2_down_bkg->SetLineWidth(2);
      
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *legsig = new TLegend(xl1,yl1,xl2,yl2);
      legsig->AddEntry(hist_BDT_puSF_nom_sig,"nominal","L");   // h1 and h2 are histogram pointers
      legsig->AddEntry(hist_BDT_puSF_up_sig,"unc +","L");
      legsig->AddEntry(hist_BDT_puSF_down_sig,"unc -","L");
      
      gStyle->SetOptStat(0);
      //PILE UP
      TCanvas* tempCanvas = new TCanvas("", "");
      TPad *histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_puSF_nom_sig->Draw("e histo");
      hist_BDT_puSF_up_sig->Draw("e same histo");
      hist_BDT_puSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      TPad *ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      TH1F* ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_puSF_nom_sig->GetNbinsX(),hist_BDT_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_puSF_nom_sig,hist_BDT_puSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      TH1F* ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_puSF_nom_sig->GetNbinsX(),hist_BDT_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_puSF_nom_sig,hist_BDT_puSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      TLine *line = new TLine(hist_BDT_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_puSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_puSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_puSF_nom_bkg->Draw("e histo");
      hist_BDT_puSF_up_bkg->Draw("e same histo");
      hist_BDT_puSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_puSF_nom_bkg->GetNbinsX(),hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_puSF_nom_bkg,hist_BDT_puSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_puSF_nom_bkg->GetNbinsX(),hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_puSF_nom_bkg,hist_BDT_puSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_puSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_puSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_puSF_nom_bkg_LogY.png");
      
      
      //ELECTRON
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_electronSF_nom_sig->Draw("e histo");
      hist_BDT_electronSF_up_sig->Draw("e same histo");
      hist_BDT_electronSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_electronSF_nom_sig->GetNbinsX(),hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_electronSF_nom_sig,hist_BDT_electronSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_electronSF_nom_sig->GetNbinsX(),hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_electronSF_nom_sig,hist_BDT_electronSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_electronSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_electronSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_electronSF_nom_bkg->Draw("e histo");
      hist_BDT_electronSF_up_bkg->Draw("e same histo");
      hist_BDT_electronSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_electronSF_nom_bkg->GetNbinsX(),hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_electronSF_nom_bkg,hist_BDT_electronSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_electronSF_nom_bkg->GetNbinsX(),hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_electronSF_nom_bkg,hist_BDT_electronSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_electronSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_electronSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_electronSF_nom_bkg_LogY.png");
      
      
      
      //MUON
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_muonSF_nom_sig->Draw("e histo");
      hist_BDT_muonSF_up_sig->Draw("e same histo");
      hist_BDT_muonSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_muonSF_nom_sig->GetNbinsX(),hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_muonSF_nom_sig,hist_BDT_muonSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_muonSF_nom_sig->GetNbinsX(),hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_muonSF_nom_sig,hist_BDT_muonSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_muonSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_muonSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_muonSF_nom_bkg->Draw("e histo");
      hist_BDT_muonSF_up_bkg->Draw("e same histo");
      hist_BDT_muonSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_muonSF_nom_bkg->GetNbinsX(),hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_muonSF_nom_bkg,hist_BDT_muonSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_muonSF_nom_bkg->GetNbinsX(),hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_muonSF_nom_bkg,hist_BDT_muonSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_muonSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_muonSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_muonSF_nom_bkg_LogY.png");
      
      //BTAG CFERR1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_cferr1_nom_sig->Draw("e histo");
      hist_BDT_btagSF_cferr1_up_sig->Draw("e same histo");
      hist_BDT_btagSF_cferr1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_cferr1_nom_sig->GetNbinsX(),hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_cferr1_nom_sig,hist_BDT_btagSF_cferr1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_cferr1_nom_sig->GetNbinsX(),hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_cferr1_nom_sig,hist_BDT_btagSF_cferr1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_cferr1_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_cferr1_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_cferr1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_cferr1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_cferr1_nom_bkg,hist_BDT_btagSF_cferr1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_cferr1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_cferr1_nom_bkg,hist_BDT_btagSF_cferr1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr1_nom_bkg_LogY.png");
      
      //BTAG cferr2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_cferr2_nom_sig->Draw("e histo");
      hist_BDT_btagSF_cferr2_up_sig->Draw("e same histo");
      hist_BDT_btagSF_cferr2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_cferr2_nom_sig->GetNbinsX(),hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_cferr2_nom_sig,hist_BDT_btagSF_cferr2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_cferr2_nom_sig->GetNbinsX(),hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_cferr2_nom_sig,hist_BDT_btagSF_cferr2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_cferr2_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_cferr2_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_cferr2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_cferr2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_cferr2_nom_bkg,hist_BDT_btagSF_cferr2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_cferr2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_cferr2_nom_bkg,hist_BDT_btagSF_cferr2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_cferr2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_cferr2_nom_bkg_LogY.png");
      
      //BTAG lf
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lf_nom_sig->Draw("e histo");
      hist_BDT_btagSF_lf_up_sig->Draw("e same histo");
      hist_BDT_btagSF_lf_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lf_nom_sig->GetNbinsX(),hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lf_nom_sig,hist_BDT_btagSF_lf_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lf_nom_sig->GetNbinsX(),hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lf_nom_sig,hist_BDT_btagSF_lf_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lf_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lf_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lf_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_lf_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_lf_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lf_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lf_nom_bkg,hist_BDT_btagSF_lf_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lf_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lf_nom_bkg,hist_BDT_btagSF_lf_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lf_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lf_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lf_nom_bkg_LogY.png");
      
      //BTAG lfstats1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lfstats1_nom_sig->Draw("e histo");
      hist_BDT_btagSF_lfstats1_up_sig->Draw("e same histo");
      hist_BDT_btagSF_lfstats1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lfstats1_nom_sig->GetNbinsX(),hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lfstats1_nom_sig,hist_BDT_btagSF_lfstats1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lfstats1_nom_sig->GetNbinsX(),hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lfstats1_nom_sig,hist_BDT_btagSF_lfstats1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lfstats1_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_lfstats1_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_lfstats1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lfstats1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lfstats1_nom_bkg,hist_BDT_btagSF_lfstats1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lfstats1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lfstats1_nom_bkg,hist_BDT_btagSF_lfstats1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats1_nom_bkg_LogY.png");
      
      
      
      //BTAG lfstats2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lfstats2_nom_sig->Draw("e histo");
      hist_BDT_btagSF_lfstats2_up_sig->Draw("e same histo");
      hist_BDT_btagSF_lfstats2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lfstats2_nom_sig->GetNbinsX(),hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lfstats2_nom_sig,hist_BDT_btagSF_lfstats2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lfstats2_nom_sig->GetNbinsX(),hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lfstats2_nom_sig,hist_BDT_btagSF_lfstats2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_lfstats2_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_lfstats2_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_lfstats2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_lfstats2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_lfstats2_nom_bkg,hist_BDT_btagSF_lfstats2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_lfstats2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_lfstats2_nom_bkg,hist_BDT_btagSF_lfstats2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_lfstats2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_lfstats2_nom_bkg_LogY.png");
      
      
      //BTAG hf
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hf_nom_sig->Draw("e histo");
      hist_BDT_btagSF_hf_up_sig->Draw("e same histo");
      hist_BDT_btagSF_hf_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hf_nom_sig->GetNbinsX(),hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hf_nom_sig,hist_BDT_btagSF_hf_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hf_nom_sig->GetNbinsX(),hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hf_nom_sig,hist_BDT_btagSF_hf_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hf_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hf_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hf_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_hf_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_hf_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hf_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hf_nom_bkg,hist_BDT_btagSF_hf_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hf_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hf_nom_bkg,hist_BDT_btagSF_hf_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hf_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hf_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hf_nom_bkg_LogY.png");
      
      //BTAG hfstats1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hfstats1_nom_sig->Draw("e histo");
      hist_BDT_btagSF_hfstats1_up_sig->Draw("e same histo");
      hist_BDT_btagSF_hfstats1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hfstats1_nom_sig->GetNbinsX(),hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hfstats1_nom_sig,hist_BDT_btagSF_hfstats1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hfstats1_nom_sig->GetNbinsX(),hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hfstats1_nom_sig,hist_BDT_btagSF_hfstats1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hfstats1_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_hfstats1_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_hfstats1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hfstats1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hfstats1_nom_bkg,hist_BDT_btagSF_hfstats1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hfstats1_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hfstats1_nom_bkg,hist_BDT_btagSF_hfstats1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats1_nom_bkg_LogY.png");
      
      
      
      //BTAG hfstats2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hfstats2_nom_sig->Draw("e histo");
      hist_BDT_btagSF_hfstats2_up_sig->Draw("e same histo");
      hist_BDT_btagSF_hfstats2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hfstats2_nom_sig->GetNbinsX(),hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hfstats2_nom_sig,hist_BDT_btagSF_hfstats2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hfstats2_nom_sig->GetNbinsX(),hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hfstats2_nom_sig,hist_BDT_btagSF_hfstats2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_BDT_btagSF_hfstats2_nom_bkg->Draw("e histo");
      hist_BDT_btagSF_hfstats2_up_bkg->Draw("e same histo");
      hist_BDT_btagSF_hfstats2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_BDT_btagSF_hfstats2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_BDT_btagSF_hfstats2_nom_bkg,hist_BDT_btagSF_hfstats2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_BDT_btagSF_hfstats2_nom_bkg->GetNbinsX(),hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_BDT_btagSF_hfstats2_nom_bkg,hist_BDT_btagSF_hfstats2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_BDT_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_BDT_btagSF_hfstats2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_BDT_btagSF_hfstats2_nom_bkg_LogY.png");
      
      
      
      ////// CUTFLOWS
      string sysratiosfilename = "sysratio_" ;
      if(doZut && toppair) sysratiosfilename += "TTZut";
      else if(!doZut && toppair) sysratiosfilename += "TTZct";
      else if(doZut && !toppair) sysratiosfilename += "STZut";
      else if(!doZut && !toppair) sysratiosfilename += "STZct";
      TFile* sysratiosfile = TFile::Open((sysratiosfilename+".root").c_str(),"RECREATE");
      sysratiosfile->cd();
      cout << "opening" << (sysratiosfilename+".root").c_str() << endl;
      hist_BDT_puSF_nom_bkg->Write();
      hist_BDT_puSF_nom_sig->Write();
      hist_BDT_electronSF_nom_bkg->Write();
      hist_BDT_electronSF_nom_sig->Write();
      hist_BDT_muonSF_nom_bkg->Write();
      hist_BDT_muonSF_nom_sig->Write();
      hist_BDT_btagSF_cferr1_nom_bkg->Write();
      hist_BDT_btagSF_cferr1_nom_sig->Write();
      hist_BDT_btagSF_cferr2_nom_bkg->Write();
      hist_BDT_btagSF_cferr2_nom_sig->Write();
      hist_BDT_btagSF_hf_nom_bkg->Write();
      hist_BDT_btagSF_hf_nom_sig->Write();
      hist_BDT_btagSF_hfstats1_nom_bkg->Write();
      hist_BDT_btagSF_hfstats1_nom_sig->Write();
      hist_BDT_btagSF_hfstats2_nom_bkg->Write();
      hist_BDT_btagSF_hfstats2_nom_sig->Write();
      hist_BDT_btagSF_lf_nom_bkg->Write();
      hist_BDT_btagSF_lf_nom_sig->Write();
      hist_BDT_btagSF_lfstats1_nom_bkg->Write();
      hist_BDT_btagSF_lfstats1_nom_sig->Write();
      hist_BDT_btagSF_lfstats2_nom_bkg->Write();
      hist_BDT_btagSF_lfstats2_nom_sig->Write();
      
      hist_BDT_puSF_up_bkg->Write();
      hist_BDT_puSF_up_sig->Write();
      hist_BDT_electronSF_up_bkg->Write();
      hist_BDT_electronSF_up_sig->Write();
      hist_BDT_muonSF_up_bkg->Write();
      hist_BDT_muonSF_up_sig->Write();
      hist_BDT_btagSF_cferr1_up_bkg->Write();
      hist_BDT_btagSF_cferr1_up_sig->Write();
      hist_BDT_btagSF_cferr2_up_bkg->Write();
      hist_BDT_btagSF_cferr2_up_sig->Write();
      hist_BDT_btagSF_hf_up_bkg->Write();
      hist_BDT_btagSF_hf_up_sig->Write();
      hist_BDT_btagSF_hfstats1_up_bkg->Write();
      hist_BDT_btagSF_hfstats1_up_sig->Write();
      hist_BDT_btagSF_hfstats2_up_bkg->Write();
      hist_BDT_btagSF_hfstats2_up_sig->Write();
      hist_BDT_btagSF_lf_up_bkg->Write();
      hist_BDT_btagSF_lf_up_sig->Write();
      hist_BDT_btagSF_lfstats1_up_bkg->Write();
      hist_BDT_btagSF_lfstats1_up_sig->Write();
      hist_BDT_btagSF_lfstats2_up_bkg->Write();
      hist_BDT_btagSF_lfstats2_up_sig->Write();
      
      
      hist_BDT_puSF_down_bkg->Write();
      hist_BDT_puSF_down_sig->Write();
      hist_BDT_electronSF_down_bkg->Write();
      hist_BDT_electronSF_down_sig->Write();
      hist_BDT_muonSF_down_bkg->Write();
      hist_BDT_muonSF_down_sig->Write();
      hist_BDT_btagSF_cferr1_down_bkg->Write();
      hist_BDT_btagSF_cferr1_down_sig->Write();
      hist_BDT_btagSF_cferr2_down_bkg->Write();
      hist_BDT_btagSF_cferr2_down_sig->Write();
      hist_BDT_btagSF_hf_down_bkg->Write();
      hist_BDT_btagSF_hf_down_sig->Write();
      hist_BDT_btagSF_hfstats1_down_bkg->Write();
      hist_BDT_btagSF_hfstats1_down_sig->Write();
      hist_BDT_btagSF_hfstats2_down_bkg->Write();
      hist_BDT_btagSF_hfstats2_down_sig->Write();
      hist_BDT_btagSF_lf_down_bkg->Write();
      hist_BDT_btagSF_lf_down_sig->Write();
      hist_BDT_btagSF_lfstats1_down_bkg->Write();
      hist_BDT_btagSF_lfstats1_down_sig->Write();
      hist_BDT_btagSF_lfstats2_down_bkg->Write();
      hist_BDT_btagSF_lfstats2_down_sig->Write();
      
      sysratiosfile->Close();
      
      
    }
    else if(PlotSystematics && (doMTWtemplate )){
      double maximum_sig = hist_mWt_puSF_nom_sig->GetMaximum()*1.5;
      double maximum_bkg = hist_mWt_puSF_nom_bkg->GetMaximum()*1.5;
      
      hist_mWt_puSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_puSF_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_electronSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_electronSF_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_muonSF_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_muonSF_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_cferr1_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_cferr1_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_cferr2_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_cferr2_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_hf_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_hf_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_hfstats1_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_hfstats1_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_hfstats2_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_hfstats2_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_lf_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_lf_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_lfstats1_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_lfstats1_nom_sig->SetMaximum(maximum_sig);
      hist_mWt_btagSF_lfstats2_nom_bkg->SetMaximum(maximum_bkg);
      hist_mWt_btagSF_lfstats2_nom_sig->SetMaximum(maximum_sig);
      
      
      
      hist_mWt_puSF_nom_sig->SetLineColor(kRed);
      hist_mWt_puSF_up_sig->SetLineColor(kBlue);
      hist_mWt_puSF_down_sig->SetLineColor(kGreen+2);
      hist_mWt_puSF_nom_bkg->SetLineColor(kRed);
      hist_mWt_puSF_up_bkg->SetLineColor(kBlue);
      hist_mWt_puSF_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_puSF_nom_sig->SetLineWidth(2);
      hist_mWt_puSF_up_sig->SetLineWidth(2);
      hist_mWt_puSF_down_sig->SetLineWidth(2);
      hist_mWt_puSF_nom_bkg->SetLineWidth(2);
      hist_mWt_puSF_up_bkg->SetLineWidth(2);
      hist_mWt_puSF_down_bkg->SetLineWidth(2);
      
      
      hist_mWt_electronSF_nom_sig->SetLineColor(kRed);
      hist_mWt_electronSF_up_sig->SetLineColor(kBlue);
      hist_mWt_electronSF_down_sig->SetLineColor(kGreen+2);
      hist_mWt_electronSF_nom_bkg->SetLineColor(kRed);
      hist_mWt_electronSF_up_bkg->SetLineColor(kBlue);
      hist_mWt_electronSF_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_electronSF_nom_sig->SetLineWidth(2);
      hist_mWt_electronSF_up_sig->SetLineWidth(2);
      hist_mWt_electronSF_down_sig->SetLineWidth(2);
      hist_mWt_electronSF_nom_bkg->SetLineWidth(2);
      hist_mWt_electronSF_up_bkg->SetLineWidth(2);
      hist_mWt_electronSF_down_bkg->SetLineWidth(2);
      
      hist_mWt_muonSF_nom_sig->SetLineColor(kRed);
      hist_mWt_muonSF_up_sig->SetLineColor(kBlue);
      hist_mWt_muonSF_down_sig->SetLineColor(kGreen+2);
      hist_mWt_muonSF_nom_bkg->SetLineColor(kRed);
      hist_mWt_muonSF_up_bkg->SetLineColor(kBlue);
      hist_mWt_muonSF_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_muonSF_nom_sig->SetLineWidth(2);
      hist_mWt_muonSF_up_sig->SetLineWidth(2);
      hist_mWt_muonSF_down_sig->SetLineWidth(2);
      hist_mWt_muonSF_nom_bkg->SetLineWidth(2);
      hist_mWt_muonSF_up_bkg->SetLineWidth(2);
      hist_mWt_muonSF_down_bkg->SetLineWidth(2);
      
      
      
      hist_mWt_btagSF_cferr1_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_cferr1_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_cferr1_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_cferr1_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_cferr1_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_cferr1_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_cferr1_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr1_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr1_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr1_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_cferr1_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_cferr1_down_bkg->SetLineWidth(2);
      
      
      hist_mWt_btagSF_cferr2_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_cferr2_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_cferr2_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_cferr2_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_cferr2_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_cferr2_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_cferr2_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr2_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr2_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_cferr2_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_cferr2_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_cferr2_down_bkg->SetLineWidth(2);
      
      hist_mWt_btagSF_hf_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_hf_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_hf_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hf_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_hf_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_hf_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hf_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_hf_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_hf_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_hf_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hf_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hf_down_bkg->SetLineWidth(2);
      
      hist_mWt_btagSF_hfstats1_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_hfstats1_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_hfstats1_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hfstats1_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_hfstats1_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_hfstats1_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hfstats1_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats1_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats1_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats1_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hfstats1_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hfstats1_down_bkg->SetLineWidth(2);
      
      
      hist_mWt_btagSF_hfstats2_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_hfstats2_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_hfstats2_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hfstats2_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_hfstats2_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_hfstats2_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_hfstats2_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats2_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats2_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_hfstats2_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hfstats2_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_hfstats2_down_bkg->SetLineWidth(2);
      
      hist_mWt_btagSF_lf_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_lf_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_lf_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lf_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_lf_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_lf_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lf_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_lf_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_lf_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_lf_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lf_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lf_down_bkg->SetLineWidth(2);
      
      hist_mWt_btagSF_lfstats1_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_lfstats1_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_lfstats1_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lfstats1_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_lfstats1_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_lfstats1_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lfstats1_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats1_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats1_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats1_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lfstats1_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lfstats1_down_bkg->SetLineWidth(2);
      
      
      hist_mWt_btagSF_lfstats2_nom_sig->SetLineColor(kRed);
      hist_mWt_btagSF_lfstats2_up_sig->SetLineColor(kBlue);
      hist_mWt_btagSF_lfstats2_down_sig->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lfstats2_nom_bkg->SetLineColor(kRed);
      hist_mWt_btagSF_lfstats2_up_bkg->SetLineColor(kBlue);
      hist_mWt_btagSF_lfstats2_down_bkg->SetLineColor(kGreen+2);
      hist_mWt_btagSF_lfstats2_nom_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats2_up_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats2_down_sig->SetLineWidth(2);
      hist_mWt_btagSF_lfstats2_nom_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lfstats2_up_bkg->SetLineWidth(2);
      hist_mWt_btagSF_lfstats2_down_bkg->SetLineWidth(2);
      
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *legsig = new TLegend(xl1,yl1,xl2,yl2);
      legsig->AddEntry(hist_mWt_puSF_nom_sig,"nominal","L");   // h1 and h2 are histogram pointers
      legsig->AddEntry(hist_mWt_puSF_up_sig,"unc +","L");
      legsig->AddEntry(hist_mWt_puSF_down_sig,"unc -","L");
      
      gStyle->SetOptStat(0);
      
      
      gStyle->SetOptStat(0);
      //PILE UP
      TCanvas* tempCanvas = new TCanvas("", "");
      TPad *histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_puSF_nom_sig->Draw("e histo");
      hist_mWt_puSF_up_sig->Draw("e same histo");
      hist_mWt_puSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      TPad *ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      TH1F* ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_puSF_nom_sig->GetNbinsX(),hist_mWt_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_puSF_nom_sig,hist_mWt_puSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      TH1F* ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_puSF_nom_sig->GetNbinsX(),hist_mWt_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_puSF_nom_sig,hist_mWt_puSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      TLine *line = new TLine(hist_mWt_puSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_puSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_puSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_puSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_puSF_nom_bkg->Draw("e histo");
      hist_mWt_puSF_up_bkg->Draw("e same histo");
      hist_mWt_puSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_puSF_nom_bkg->GetNbinsX(),hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_puSF_nom_bkg,hist_mWt_puSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_puSF_nom_bkg->GetNbinsX(),hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_puSF_nom_bkg,hist_mWt_puSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_puSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_puSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_puSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_puSF_nom_bkg_LogY.png");
      
      
      //ELECTRON
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_electronSF_nom_sig->Draw("e histo");
      hist_mWt_electronSF_up_sig->Draw("e same histo");
      hist_mWt_electronSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_electronSF_nom_sig->GetNbinsX(),hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_electronSF_nom_sig,hist_mWt_electronSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_electronSF_nom_sig->GetNbinsX(),hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_electronSF_nom_sig,hist_mWt_electronSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_electronSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_electronSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_electronSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_electronSF_nom_bkg->Draw("e histo");
      hist_mWt_electronSF_up_bkg->Draw("e same histo");
      hist_mWt_electronSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_electronSF_nom_bkg->GetNbinsX(),hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_electronSF_nom_bkg,hist_mWt_electronSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_electronSF_nom_bkg->GetNbinsX(),hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_electronSF_nom_bkg,hist_mWt_electronSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_electronSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_electronSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_electronSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_electronSF_nom_bkg_LogY.png");
      
      
      
      //MUON
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_muonSF_nom_sig->Draw("e histo");
      hist_mWt_muonSF_up_sig->Draw("e same histo");
      hist_mWt_muonSF_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_muonSF_nom_sig->GetNbinsX(),hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_muonSF_nom_sig,hist_mWt_muonSF_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_muonSF_nom_sig->GetNbinsX(),hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_muonSF_nom_sig,hist_mWt_muonSF_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_muonSF_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_muonSF_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_muonSF_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_muonSF_nom_bkg->Draw("e histo");
      hist_mWt_muonSF_up_bkg->Draw("e same histo");
      hist_mWt_muonSF_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_muonSF_nom_bkg->GetNbinsX(),hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_muonSF_nom_bkg,hist_mWt_muonSF_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_muonSF_nom_bkg->GetNbinsX(),hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_muonSF_nom_bkg,hist_mWt_muonSF_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_muonSF_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_muonSF_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_muonSF_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_muonSF_nom_bkg_LogY.png");
      
      //BTAG CFERR1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_cferr1_nom_sig->Draw("e histo");
      hist_mWt_btagSF_cferr1_up_sig->Draw("e same histo");
      hist_mWt_btagSF_cferr1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_cferr1_nom_sig->GetNbinsX(),hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_cferr1_nom_sig,hist_mWt_btagSF_cferr1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_cferr1_nom_sig->GetNbinsX(),hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_cferr1_nom_sig,hist_mWt_btagSF_cferr1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_cferr1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_cferr1_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_cferr1_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_cferr1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_cferr1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_cferr1_nom_bkg,hist_mWt_btagSF_cferr1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_cferr1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_cferr1_nom_bkg,hist_mWt_btagSF_cferr1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_cferr1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr1_nom_bkg_LogY.png");
      
      //BTAG cferr2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_cferr2_nom_sig->Draw("e histo");
      hist_mWt_btagSF_cferr2_up_sig->Draw("e same histo");
      hist_mWt_btagSF_cferr2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_cferr2_nom_sig->GetNbinsX(),hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_cferr2_nom_sig,hist_mWt_btagSF_cferr2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_cferr2_nom_sig->GetNbinsX(),hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_cferr2_nom_sig,hist_mWt_btagSF_cferr2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_cferr2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_cferr2_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_cferr2_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_cferr2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_cferr2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_cferr2_nom_bkg,hist_mWt_btagSF_cferr2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_cferr2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_cferr2_nom_bkg,hist_mWt_btagSF_cferr2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_cferr2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_cferr2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_cferr2_nom_bkg_LogY.png");
      
      //BTAG lf
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lf_nom_sig->Draw("e histo");
      hist_mWt_btagSF_lf_up_sig->Draw("e same histo");
      hist_mWt_btagSF_lf_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lf_nom_sig->GetNbinsX(),hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lf_nom_sig,hist_mWt_btagSF_lf_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lf_nom_sig->GetNbinsX(),hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lf_nom_sig,hist_mWt_btagSF_lf_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lf_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lf_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lf_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_lf_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_lf_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lf_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lf_nom_bkg,hist_mWt_btagSF_lf_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lf_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lf_nom_bkg,hist_mWt_btagSF_lf_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lf_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lf_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lf_nom_bkg_LogY.png");
      
      //BTAG lfstats1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lfstats1_nom_sig->Draw("e histo");
      hist_mWt_btagSF_lfstats1_up_sig->Draw("e same histo");
      hist_mWt_btagSF_lfstats1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lfstats1_nom_sig->GetNbinsX(),hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lfstats1_nom_sig,hist_mWt_btagSF_lfstats1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lfstats1_nom_sig->GetNbinsX(),hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lfstats1_nom_sig,hist_mWt_btagSF_lfstats1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lfstats1_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_lfstats1_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_lfstats1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lfstats1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lfstats1_nom_bkg,hist_mWt_btagSF_lfstats1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lfstats1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lfstats1_nom_bkg,hist_mWt_btagSF_lfstats1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats1_nom_bkg_LogY.png");
      
      
      
      //BTAG lfstats2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lfstats2_nom_sig->Draw("e histo");
      hist_mWt_btagSF_lfstats2_up_sig->Draw("e same histo");
      hist_mWt_btagSF_lfstats2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lfstats2_nom_sig->GetNbinsX(),hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lfstats2_nom_sig,hist_mWt_btagSF_lfstats2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lfstats2_nom_sig->GetNbinsX(),hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lfstats2_nom_sig,hist_mWt_btagSF_lfstats2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_lfstats2_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_lfstats2_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_lfstats2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_lfstats2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_lfstats2_nom_bkg,hist_mWt_btagSF_lfstats2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_lfstats2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_lfstats2_nom_bkg,hist_mWt_btagSF_lfstats2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_lfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_lfstats2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_lfstats2_nom_bkg_LogY.png");
      
      
      //BTAG hf
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hf_nom_sig->Draw("e histo");
      hist_mWt_btagSF_hf_up_sig->Draw("e same histo");
      hist_mWt_btagSF_hf_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hf_nom_sig->GetNbinsX(),hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hf_nom_sig,hist_mWt_btagSF_hf_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hf_nom_sig->GetNbinsX(),hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hf_nom_sig,hist_mWt_btagSF_hf_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hf_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hf_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hf_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hf_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_hf_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_hf_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hf_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hf_nom_bkg,hist_mWt_btagSF_hf_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hf_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hf_nom_bkg,hist_mWt_btagSF_hf_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hf_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hf_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hf_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hf_nom_bkg_LogY.png");
      
      //BTAG hfstats1
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hfstats1_nom_sig->Draw("e histo");
      hist_mWt_btagSF_hfstats1_up_sig->Draw("e same histo");
      hist_mWt_btagSF_hfstats1_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hfstats1_nom_sig->GetNbinsX(),hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hfstats1_nom_sig,hist_mWt_btagSF_hfstats1_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hfstats1_nom_sig->GetNbinsX(),hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hfstats1_nom_sig,hist_mWt_btagSF_hfstats1_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hfstats1_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats1_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats1_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hfstats1_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_hfstats1_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_hfstats1_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hfstats1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hfstats1_nom_bkg,hist_mWt_btagSF_hfstats1_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hfstats1_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hfstats1_nom_bkg,hist_mWt_btagSF_hfstats1_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hfstats1_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats1_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats1_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats1_nom_bkg_LogY.png");
      
      
      
      //BTAG hfstats2
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hfstats2_nom_sig->Draw("e histo");
      hist_mWt_btagSF_hfstats2_up_sig->Draw("e same histo");
      hist_mWt_btagSF_hfstats2_down_sig->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hfstats2_nom_sig->GetNbinsX(),hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_sig->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hfstats2_nom_sig,hist_mWt_btagSF_hfstats2_up_sig);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hfstats2_nom_sig->GetNbinsX(),hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_sig->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hfstats2_nom_sig,hist_mWt_btagSF_hfstats2_down_sig);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hfstats2_nom_sig->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_sig->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats2_nom_sig.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats2_nom_sig_LogY.png");
      
      
      tempCanvas = new TCanvas("", "");
      histoPad = new TPad("histoPad","histoPad",0,0.3,1,1);
      // pad1->SetLogy(1);
      histoPad->Draw();
      histoPad->SetBottomMargin(0.);
      histoPad->cd();
      
      hist_mWt_btagSF_hfstats2_nom_bkg->Draw("e histo");
      hist_mWt_btagSF_hfstats2_up_bkg->Draw("e same histo");
      hist_mWt_btagSF_hfstats2_down_bkg->Draw("e same histo");
      legsig->Draw("same");
      histoPad->Update();
      tempCanvas->cd();
      ratioPad = new TPad("ratioPad","ratioPad",0,0,1,0.3);
      ratioPad->SetTopMargin(0);
      ratioPad->SetBottomMargin(0.35);
      //   pad2->SetLogy();
      ratioPad->Draw();
      ratioPad->cd();
      ratioPad->SetGridy();
      ratioUp = new TH1F("ratioUp","ratioUp",hist_mWt_btagSF_hfstats2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_bkg->GetNbinsX()));
      ratioUp->Divide(hist_mWt_btagSF_hfstats2_nom_bkg,hist_mWt_btagSF_hfstats2_up_bkg);
      ratioUp->SetMarkerColor(kBlue);
      ratioUp->SetMarkerStyle(9);
      ratioUp->SetStats(0);
      ratioUp->SetTitle("");
      ratioUp->GetYaxis()->SetRangeUser(0.8, 1.2);
      ratioUp->GetYaxis()->SetTitle("#frac{nominal}{uncertainty}");
      ratioUp->GetYaxis()->CenterTitle(kTRUE);
      ratioUp->GetYaxis()->SetTitleSize(0.11);
      ratioUp->GetYaxis()->SetLabelSize(0.1);
      ratioUp->GetYaxis()->SetTitleOffset(0.4);
      ratioUp->GetXaxis()->SetTitle("D");
      ratioUp->GetXaxis()->SetTitleSize(0.12);
      ratioUp->GetXaxis()->SetLabelSize(0.1);
      ratioUp->Draw("P");
      
      ratioDown = new TH1F("ratioDown","ratioDown",hist_mWt_btagSF_hfstats2_nom_bkg->GetNbinsX(),hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0), hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_bkg->GetNbinsX()));
      ratioDown->Divide(hist_mWt_btagSF_hfstats2_nom_bkg,hist_mWt_btagSF_hfstats2_down_bkg);
      ratioDown->SetMarkerColor(kGreen+2);
      ratioDown->SetMarkerStyle(9);
      ratioDown->SetLineColor(kGreen+2);
      ratioDown->Draw("P SAME");
      
      
      line = new TLine(hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinLowEdge(0),1,hist_mWt_btagSF_hfstats2_nom_bkg->GetXaxis()->GetBinUpEdge(hist_mWt_btagSF_hfstats2_nom_bkg->GetNbinsX()),1);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->Draw();
      tempCanvas->cd();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats2_nom_bkg.png");
      histoPad->SetLogy();
      ratioPad->SetLogy();
      tempCanvas->SetLogy();
      tempCanvas->Update();
      tempCanvas->SaveAs("hist_mWt_btagSF_hfstats2_nom_bkg_LogY.png");
      
    }
    
    
    if(PlotMVAvars && !doMTWtemplate && !doTTZtemplate){
      // cout << "plot mva vars " << endl;
      TDirectory* th1dirmva = fout->mkdir("1D_MVA_histograms");
      th1dirmva->cd();
      gStyle->SetOptStat(0);
      
      std::vector<string> channellist;
      channellist.push_back("all");
      channellist.push_back("eee");
      channellist.push_back("uue");
      channellist.push_back("eeu");
      channellist.push_back("uuu");
      // find all variables
      std::string splitname = "";
      char seperator = '_';
      std::vector < std::string > variables;
      for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
      {
        
        if(toppair && it->first.find("TT_FCNC")!=std::string::npos && it->first.find("all")!=std::string::npos){}
        else if(!toppair && it->first.find("ST_FCNC")!=std::string::npos && it->first.find("all")!=std::string::npos){
          //cout << it->first << endl;
          
        }
        else continue;
        
        splitname = (split(it->first, seperator))[0];
        
        variables.push_back(splitname);
      }
      for(Int_t i = 0 ; i < variables.size(); i++){
        splitname = variables[i];
        // cout << "name " << splitname << endl;
        for(Int_t iChan = 0; iChan < channellist.size(); iChan++){
          TH1F *tempBKG(0);
          TH1F *tempSignal(0);
          TH1F *tempfakes = new TH1F();
          string title = "";
          for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
          {
            if(it->first.find(channellist[iChan].c_str())==std::string::npos) continue;
            if(it->first.find(splitname.c_str())==std::string::npos) continue;
            if(it->first.find("data")!=std::string::npos) continue;
            cout << "looking at " << it->first << endl;
            TH1F *temp = it->second;
            
            /*  if(it->first.find("fake")!=std::string::npos) {
             // cout << "filling fake " << tempfake << " " << temp << endl;
             tempfakes = (TH1F*) temp->Clone("");
             
             }
             else */if(it->first.find("T_FCNC")!=std::string::npos) {
               if(tempSignal == 0){ tempSignal = (TH1F*) temp->Clone(); title =  temp->GetTitle();}
               else tempSignal->Add(temp);
             }
             else if(it->first.find("T_FCNC")==std::string::npos){
               if(tempBKG == 0) tempBKG = (TH1F*) temp->Clone();
               else tempBKG->Add(temp);
             }
            //delete temp;
          }
          //cout << "filled histos " << endl;
          if(tempBKG == 0) cout << "ERROR tempBKG is null" << endl ;
          // if(tempfakes == 0) cout << "ERROR tempfake is null" << endl ;
          if(tempSignal == 0) cout << "ERROR tempSignal is null" << endl ;
          
          // tempfake->SetLineColor(kGreen);
          tempBKG->SetName(splitname.c_str());
          if(doZut && !toppair) tempBKG->SetTitle(("Shape comparison ST tZu - " + channellist[iChan]).c_str());
          if(!doZut && !toppair) tempBKG->SetTitle(("Shape comparison ST tZc - " + channellist[iChan]).c_str());
          if(doZut && toppair) tempBKG->SetTitle(("Shape comparison TT tZu - " + channellist[iChan]).c_str());
          if(!doZut && toppair) tempBKG->SetTitle(("Shape comparison TT tZc - " + channellist[iChan]).c_str());
          
          
          
          
          Double_t scaleBKG = 1./tempBKG->Integral();
          tempBKG->Scale(scaleBKG);
          Double_t scaleSig = 1./tempSignal->Integral();
          tempSignal->Scale(scaleSig);
          // Double_t scalefake = 1./tempfake->Integral();
          // tempfake->Scale(scalefake);
          Double_t max = TMath::Max(tempSignal->GetMaximum(), tempBKG->GetMaximum());
          // Double_t max = TMath::Max(max0, tempfake->GetMaximum());
          
          
          
          tempBKG->SetLineColor(kBlue);
          tempSignal->SetLineColor(kRed);
          tempBKG->SetMarkerColor(kBlue);
          tempSignal->SetMarkerColor(kRed);
          tempBKG->SetLineWidth(2);
          tempBKG->SetMarkerStyle(2);
          tempSignal->SetLineWidth(2);
          tempSignal->SetMarkerStyle(2);
          // cout << "title " << title.c_str() << endl;
          tempBKG->GetXaxis()->SetTitle(title.c_str());
          tempBKG->GetYaxis()->SetTitle("Nb. Norm. Events");
          tempBKG->SetMaximum(max*1.7);
          Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
          TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
          leg->AddEntry(tempSignal,"FCNC: ST + TT","P");   // h1 and h2 are histogram pointers
          leg->AddEntry(tempBKG,"SM bkg","P");
          //leg->AddEntry(tempfake,"DD non prompt","L");
          
          TCanvas* tempCanvas = TCanvasCreator(tempBKG,"Normalised MVA distribution" );
          tempBKG->Draw("e hist");
          tempSignal->Draw("e hist same ");
          //tempfake->Draw("L,Sames");
          tempCanvas->Update();
          leg->Draw("Same");
          tempCanvas->Update();
          tempCanvas->SaveAs( (placeTH1F+splitname+"_"+channellist[iChan]+".png").c_str() );
          
          tempBKG->Write();
          tempSignal->Write();
          // tempfake->Write();
          delete tempCanvas;
          delete tempSignal;
          //delete tempfake;
          delete tempBKG;
        } // channellist
      }
    }
    
    
    // if(false){ // TO FIX
    if(makePlots && doMTWtemplate && !doTTZtemplate){
      hist_WZ->SetLineColor(kBlue);
      hist_fakes->SetLineColor(kGreen);
      hist_TT_FCNC->SetLineColor(kRed-2);
      
      
      
      hist_WZ->SetName("M_T(W)");
      hist_WZ->SetTitle("Shape comparison");
      
      Double_t scaleBKG_nom = 1./hist_WZ->Integral();
      hist_WZ->Scale(scaleBKG_nom);
      Double_t scaleSigTT_nom= 1./hist_TT_FCNC->Integral();
      hist_TT_FCNC->Scale(scaleSigTT_nom);
      Double_t scalefake_nom = 1./hist_fakes->Integral();
      cout << "scale fakes " << scalefake_nom << endl;
      hist_fakes->Scale(scalefake_nom);
      Double_t max0 = TMath::Max(hist_fakes->GetMaximum(), hist_TT_FCNC->GetMaximum());
      Double_t max1 = TMath::Max(hist_TT_FCNC->GetMaximum(), hist_WZ->GetMaximum());
      Double_t maximum = TMath::Max(max0, max1);
      hist_WZ->SetMaximum(maximum*1.2);
      hist_WZ->GetXaxis()->SetTitle("M_T(W)");
      hist_WZ->GetYaxis()->SetTitle("Nb. Norm. Events");
      
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
      leg->AddEntry(hist_TT_FCNC,"Signal TT","L");
      leg->AddEntry(hist_WZ,"WZ background","L");
      leg->AddEntry(hist_fakes,"DD non prompt","L");
      
      TCanvas* tempCanvas = TCanvasCreator(hist_WZ,"Normalised distribution" );
      hist_WZ->Draw("h");
      hist_TT_FCNC->Draw("h,sames");
      hist_fakes->Draw("h,sames");
      leg->Draw("Same");
      tempCanvas->SaveAs( (placeTH1F+"MWT.png").c_str() );
      
      hist_WZ->Write();
      hist_TT_FCNC->Write();
      hist_fakes->Write();
      
      
      /*
       cout << "plotting mtW shapes " << endl;
       std::vector<string> channellist;
       channellist.push_back("all");
       channellist.push_back("eee");
       channellist.push_back("uue");
       channellist.push_back("eeu");
       channellist.push_back("uuu");
       
       for(Int_t iChan = 0; iChan < channellist.size(); iChan++){
       TH1F *tempBKG_nom(0);
       TH1F *tempSignalST_nom(0);
       TH1F *tempfake_nom(0);
       TH1F *tempBKG_up(0);
       TH1F *tempSignalST_up(0);
       TH1F *tempBKG_down(0);
       TH1F *tempSignalST_down(0);
       TH1F *tempSignalTT_down(0);
       TH1F *tempSignalTT_up(0);
       TH1F *tempSignalTT_nom(0);
       
       //cout << "histo mtw size " << histo1DMTW.size() << endl;
       for (map<string,TH1F*>::const_iterator it = histo1DMTW.begin(); it != histo1DMTW.end(); it++)
       {
       // cout << "looking at " << it->first << " and the channel to keep " << channellist[iChan].c_str() << endl;
       
       if(it->first.find(channellist[iChan].c_str())==std::string::npos){
       //cout << "continuing " << endl;
       continue;
       }
       
       TH1F* temp = it->second;
       if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
       if(tempBKG_nom == 0) tempBKG_nom = (TH1F*) temp->Clone();
       else tempBKG_nom->Add(temp);
       }
       else if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
       if(tempBKG_up == 0) tempBKG_up = (TH1F*) temp->Clone();
       else tempBKG_up->Add(temp);
       }
       else if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
       if(tempBKG_down == 0) tempBKG_down= (TH1F*) temp->Clone();
       else tempBKG_down->Add(temp);
       }
       if(it->first.find("fake")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
       if(tempfake_nom == 0) tempfake_nom = (TH1F*) temp->Clone();
       else tempfake_nom->Add(temp);
       }
       if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
       cout << "st nom" << endl;
       if(tempSignalST_nom == 0) tempSignalST_nom = (TH1F*) temp->Clone();
       else tempSignalST_nom->Add(temp);
       }
       else if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
       if(tempSignalST_up == 0) tempSignalST_up = (TH1F*) temp->Clone();
       else tempSignalST_up->Add(temp);
       }
       else if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
       if(tempSignalST_down == 0) tempSignalST_down= (TH1F*) temp->Clone();
       else tempSignalST_down->Add(temp);
       }
       if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
       cout << "tt nom" << endl;
       if(tempSignalTT_nom == 0) tempSignalTT_nom = (TH1F*) temp->Clone();
       else tempSignalTT_nom->Add(temp);
       }
       else if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
       if(tempSignalTT_up == 0) tempSignalTT_up = (TH1F*) temp->Clone();
       else tempSignalTT_up->Add(temp);
       }
       else if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
       if(tempSignalTT_down == 0) tempSignalTT_down= (TH1F*) temp->Clone();
       else tempSignalTT_down->Add(temp);
       }
       delete temp;
       }
       //cout << "filled histos " << endl;
       if(tempBKG_nom == 0) cout << "ERROR tempBKG is null" << endl ;
       if(tempfake_nom == 0) cout << "ERROR tempfake is null" << endl ;
       if(tempSignalST_nom == 0) cout << "ERROR tempSignalST is null" << endl ;
       if(tempSignalTT_nom == 0) cout << "ERROR tempSignalTT is null" << endl ;
       
       tempBKG_nom->SetLineColor(kBlue);
       tempBKG_up->SetLineColor(kBlue);
       tempBKG_down->SetLineColor(kBlue);
       tempBKG_up->SetLineStyle(8);
       tempBKG_down->SetLineStyle(3);
       tempSignalST_nom->SetLineColor(kRed);
       tempSignalST_up->SetLineColor(kRed);
       tempSignalST_down->SetLineColor(kRed);
       tempSignalST_up->SetLineStyle(8);
       tempSignalST_down->SetLineStyle(3);
       tempfake_nom->SetLineColor(kGreen);
       tempSignalTT_nom->SetLineColor(kRed-2);
       tempSignalTT_up->SetLineColor(kRed-2);
       tempSignalTT_down->SetLineColor(kRed-2);
       tempSignalTT_up->SetLineStyle(8);
       tempSignalTT_down->SetLineStyle(3);
       
       
       tempBKG_nom->SetName("M_T(W)");
       tempBKG_nom->SetTitle(("Shape comparison - " + channellist[iChan]).c_str());
       
       Double_t scaleBKG_nom = 1./tempBKG_nom->Integral();
       tempBKG_nom->Scale(scaleBKG_nom);
       Double_t scaleBKG_up= 1./tempBKG_up->Integral();
       tempBKG_up->Scale(scaleBKG_up);
       Double_t scaleBKG_down = 1./tempBKG_down->Integral();
       tempBKG_down->Scale(scaleBKG_down);
       Double_t scaleSigST_nom = 1./tempSignalST_nom->Integral();
       tempSignalST_nom->Scale(scaleSigST_nom);
       Double_t scaleSigST_down = 1./tempSignalST_down->Integral();
       tempSignalST_down->Scale(scaleSigST_down);
       Double_t scaleSigST_up = 1./tempSignalST_up->Integral();
       tempSignalST_up->Scale(scaleSigST_up);
       Double_t scaleSigTT_nom = 1./tempSignalTT_nom->Integral();
       tempSignalST_nom->Scale(scaleSigST_nom);
       Double_t scaleSigTT_down = 1./tempSignalTT_down->Integral();
       tempSignalTT_down->Scale(scaleSigTT_down);
       Double_t scaleSigTT_up = 1./tempSignalTT_up->Integral();
       tempSignalTT_up->Scale(scaleSigTT_up);
       Double_t scalefake_nom = 1./tempfake_nom->Integral();
       tempfake_nom->Scale(scalefake_nom);
       Double_t max0 = TMath::Max(tempSignalST_nom->GetMaximum(), tempBKG_nom->GetMaximum());
       Double_t max1 = TMath::Max(tempSignalST_nom->GetMaximum(), tempfake_nom->GetMaximum());
       Double_t max2 = TMath::Max(tempSignalST_nom->GetMaximum(), tempSignalTT_nom->GetMaximum());
       Double_t max00 = TMath::Max(tempSignalST_up->GetMaximum(), tempBKG_up->GetMaximum());
       Double_t max10 = TMath::Max(tempSignalST_up->GetMaximum(), tempfake_nom->GetMaximum());
       Double_t max20 = TMath::Max(tempSignalST_up->GetMaximum(), tempSignalTT_up->GetMaximum());
       Double_t max01 = TMath::Max(tempSignalST_down->GetMaximum(), tempBKG_down->GetMaximum());
       Double_t max11 = TMath::Max(tempSignalST_down->GetMaximum(), tempfake_nom->GetMaximum());
       Double_t max21 = TMath::Max(tempSignalST_down->GetMaximum(), tempSignalTT_down->GetMaximum());
       Double_t maxA = TMath::Max(max0, max1);
       Double_t maxB = TMath::Max(max2,max00);
       Double_t maxC = TMath::Max(max10,max20);
       Double_t maxD = TMath::Max(max01,max11);
       Double_t maxX = TMath::Max(max21,maxA);
       Double_t maxY = TMath::Max(maxB,maxC);
       Double_t maxi = TMath::Max(maxD,maxX);
       Double_t maximum = TMath::Max(maxi, maxY);
       tempBKG_nom->SetMaximum(maximum*1.2);
       tempBKG_nom->GetXaxis()->SetTitle("M_T(W)");
       tempBKG_nom->GetYaxis()->SetTitle("Nb. Events");
       
       
       Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
       TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
       leg->AddEntry(tempSignalST_nom,"Signal ST","L");   // h1 and h2 are histogram pointers
       leg->AddEntry(tempSignalTT_nom,"Signal TT","L");
       leg->AddEntry(tempBKG_nom,"WZ background","L");
       leg->AddEntry(tempfake_nom,"DD non prompt","L");
       
       TCanvas* tempCanvas = TCanvasCreator(tempBKG_nom,"Normalised distribution" );
       tempBKG_nom->Draw("h");
       tempBKG_up->Draw("h,sames");
       tempBKG_down->Draw("h,sames");
       tempSignalST_nom->Draw("h,Sames");
       tempSignalST_up->Draw("h,Sames");
       tempSignalST_down->Draw("h,Sames");
       tempSignalTT_nom->Draw("h,Sames");
       tempSignalTT_up->Draw("h,Sames");
       tempSignalTT_down->Draw("h,Sames");
       tempfake_nom->Draw("L,Sames");
       leg->Draw("Same");
       tempCanvas->SaveAs( (placeTH1F+"MWT_"+channellist[iChan]+".png").c_str() );
       
       tempBKG_nom->Write();
       tempSignalST_nom->Write();
       tempSignalTT_nom->Write();
       tempfake_nom->Write();
       tempBKG_up->Write();
       tempSignalST_up->Write();
       tempSignalTT_up->Write();
       tempBKG_down->Write();
       tempSignalST_down->Write();
       tempSignalTT_down->Write();
       
       delete tempCanvas;
       delete tempSignalST_nom;
       delete tempSignalTT_nom;
       delete tempfake_nom;
       delete tempBKG_nom;
       delete tempBKG_up;
       delete tempSignalST_up;
       delete tempSignalTT_up;
       delete tempBKG_down;
       delete tempSignalST_down;
       delete tempSignalTT_down;
       
       } // channellist
       */
    }
    fout->Write();
    fout->Close();
    
    delete fout;
    
  }
  
  /* delete hist_WZ;
   delete hist_TT_FCNC;
   delete hist_fakes;*/
  ///************************************///
  ///   ADD PDF UNC TO COMBINE TEMPLATE  ///
  ///************************************///
 /* if(doPDFunc && doMTWtemplate && !doTTZtemplate){
    TFile* combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "UPDATE" );
    combinetemplate_file->cd();
    
    TH1::SetDefaultSumw2();
    TH1F *hist_uuu     = new TH1F( (coupling + "_" + region+"_uuu").c_str(),           (coupling + "_" + region+"_uuu").c_str(),           nbin,BDT_begin,BDT_end );
    TH1F *hist_uue     = new TH1F( (coupling + "_" + region+"_uue").c_str(),           (coupling + "_" + region+"_uue").c_str(),           nbin,BDT_begin,BDT_end );
    TH1F *hist_eeu     = new TH1F( (coupling + "_" + region+"_eeu").c_str(),           (coupling + "_" + region+"_eeu").c_str(),           nbin,BDT_begin,BDT_end );
    TH1F *hist_eee     = new TH1F( (coupling + "_" + region+"_eee").c_str(),           (coupling + "_" + region+"_eee").c_str(),           nbin,BDT_begin,BDT_end );
    
    //NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)] FIX ME
    output_histo_name = "";
    vector<string> v_sys = {"PDFup", "PDFdown"};
    dataSetName  = "WZTo3LNu_3Jets_MLL50_80X";
    for(Int_t isys = 0; isys < v_sys.size() ; isys++)
    {
      systematic = v_sys[isys];
      output_histo_name = coupling + "_" + region+"_uuu_"  + dataSetName + "_" + systematic ;
      hist_uuu->SetTitle(output_histo_name.c_str());
      hist_uuu->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_uue_"  + dataSetName + "_" + systematic ;
      hist_uue->SetTitle(output_histo_name.c_str());
      hist_uue->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_eeu_"  + dataSetName + "_" + systematic ;
      hist_eeu->SetTitle(output_histo_name.c_str());
      hist_eeu->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_eee_"  + dataSetName + "_" + systematic ;
      hist_eee->SetTitle(output_histo_name.c_str());
      hist_eee->Write(output_histo_name.c_str());
    }
    combinetemplate_file->Close();
    delete combinetemplate_file;
    //delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
  }
  
  */
  
  
  
  Double_t time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    Int_t mins = time/60;
    float secs = time - mins*60;
    
    if (mins >= 60 )
    {
      Int_t hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
  
  //fin->Close();
   delete fin;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  return 0;
}



///// BOOK KEEPING FUNCTIONS
std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    if (end != start) {
      tokens.push_back(text.substr(start, end - start));
    }
    start = end + 1;
  }
  if (end != start) {
    tokens.push_back(text.substr(start));
  }
  return tokens;
}
string ConvertIntToString(Int_t Number, Int_t pad){
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}
string MakeTimeStamp(){
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  Int_t year = now->tm_year - 100;  /// + 1900 to get current year
  Int_t month = now->tm_mon + 1;
  Int_t day = now->tm_mday;
  Int_t hour = now->tm_hour;
  Int_t min = now->tm_min;
  //Int_t sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}
string intToStr (Int_t number){
  ostringstream buff;
  buff<<number;
  return buff.str();
}
Double_t maximumValue(vector<double> array){
  Int_t length = array.size();  // establish size of array
  Double_t max = array[0];       // start with max = first element
  
  for(Int_t i = 1; i<length; i++)
  {
    if(array[i] > max)
      max = array[i];
  }
  return max;                // return highest value in array
}
Double_t minimumValue(vector<double> array){
  Int_t length = array.size();  // establish size of array
  Double_t max = array[0];       // start with max = first element
  
  for(Int_t i = 1; i<length; i++)
  {
    if(array[i] < max)
      max = array[i];
  }
  return max;                // return highest value in array
}

//// INITIALISATIONS
void InitMSPlotsMTW(string prefix, vector <int> decayChannels){
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    decaystring += prefix;
    
    //cout << "init msplots " << endl;
    MSPlotMTW[("MTW_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("MTW_"+decaystring).c_str(), nbinMTW, 0,endMTW, "transv. mass W boson","GeV");
    
/*
    MSPlotMTW[("DeltaR_NonIsoLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_NonIsoLepJet_"+decaystring).c_str(), 20,0.,8., "#Delta R (non iso, jet)","units");
    MSPlotMTW[("DeltaR_lep0Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep0Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep0, jet)","units");
    MSPlotMTW[("DeltaR_lep1Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep1Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep1, jet)","units");
    MSPlotMTW[("DeltaR_lep2Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep2Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep2, jet)","units");
    //MSPlotMTW[("DeltaR_MinLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_MinLepJet_"+decaystring).c_str(), 20,0.,8., "min. #Delta R (lep, jet)","units");
  //  MSPlotMTW[("DeltaR_MinLepNonIsoJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_MinLepNonIsoJet_"+decaystring).c_str(), 20,0.,8., "min #Delta R (non iso, jet)","units");
   
    MSPlotMTW[("ZBOSON_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("ZBOSONM_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
    MSPlotMTW[("ZBOSON_Mcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("ZBOSONMcut_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
    //MSPlotMTW[("NJETSCSVM_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("NJETSCSVM_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");
    //MSPlotMTW[("NJETSCSVMcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("NJETSCSVMcut_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");*/
  }
  decaystring = "";
}
void InitMSPlotsTTZ(string prefix, vector <int> decayChannels){
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    decaystring += prefix;
    
    //cout << "init msplots " << endl;
    MSPlotTTZ[("TTZ_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("TTZ_"+decaystring).c_str(), nbinTTZ, beginTTZ,endTTZ, "","units");
    
    /*
     MSPlotMTW[("DeltaR_NonIsoLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_NonIsoLepJet_"+decaystring).c_str(), 20,0.,8., "#Delta R (non iso, jet)","units");
     MSPlotMTW[("DeltaR_lep0Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep0Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep0, jet)","units");
     MSPlotMTW[("DeltaR_lep1Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep1Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep1, jet)","units");
     MSPlotMTW[("DeltaR_lep2Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_lep2Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep2, jet)","units");
     //MSPlotMTW[("DeltaR_MinLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_MinLepJet_"+decaystring).c_str(), 20,0.,8., "min. #Delta R (lep, jet)","units");
     //  MSPlotMTW[("DeltaR_MinLepNonIsoJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("DeltaR_MinLepNonIsoJet_"+decaystring).c_str(), 20,0.,8., "min #Delta R (non iso, jet)","units");
     
     MSPlotMTW[("ZBOSON_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("ZBOSONM_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
     MSPlotMTW[("ZBOSON_Mcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("ZBOSONMcut_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
     //MSPlotMTW[("NJETSCSVM_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("NJETSCSVM_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");
     //MSPlotMTW[("NJETSCSVMcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("NJETSCSVMcut_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");*/
  }
  decaystring = "";
}
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair, bool isZut){
  clock_t start_sub = clock();
  
  prefix = prefix + "_";
  
  // control plots
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    //cout << "init " << (prefix+"_BDT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), nbin,BDT_begin,BDT_end, "D");
    MSPlot[ (prefix+"channel_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"channel_"+decaystring).c_str(), 5,-0.5, 4.5, "decay");
  //  MSPlot[ (prefix+"weight_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"weight_"+decaystring).c_str(), 100,0, 0.3, "eventweight");
    /*
    MSPlot[(prefix+"DeltaR_NonIsoLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_NonIsoLepJet_"+decaystring).c_str(), 20,0.,8., "#Delta R (non iso, jet)","units");
   
    MSPlot[(prefix+"DeltaR_lep0Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_lep0Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep0, jet)","units");
    MSPlot[(prefix+"DeltaR_lep1Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_lep1Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep1, jet)","units");
    MSPlot[(prefix+"DeltaR_lep2Jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_lep2Jet_"+decaystring).c_str(), 20,0.,8., "#Delta R (lep2, jet)","units");
   // MSPlot[(prefix+"DeltaR_MinLepJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_MinLepJet_"+decaystring).c_str(), 20,0.,8., "min. #Delta R (lep, jet)","units");
   // MSPlot[(prefix+"DeltaR_MinLepNonIsoJet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"DeltaR_MinLepNonIsoJet_"+decaystring).c_str(), 20,0.,8., "min #Delta R (non iso, jet)","units");
    
    //cout << (prefix+"ZBOSON_M_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"ZBOSON_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"ZBOSONM_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
    MSPlot[(prefix+"ZBOSON_Mcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"ZBOSONMcut_"+decaystring).c_str(), 20,70,110, "mass Zboson","GeV");
    cout << (prefix+"NJETSCSVM_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"NJETSCSVM_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"NJETSCSVM_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");
    MSPlot[(prefix+"NJETSCSVMcut_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"NJETSCSVMcut_"+decaystring).c_str(), 6,-0.5,5.5, "nb CSVM jets","units");
   */
  /*  if(!istoppair){
      MSPlot[(prefix+"mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"mlb_"+decaystring).c_str(),10, 0, 500, "inv. mass l_{W}b (GeV)","GeV");
      MSPlot[(prefix+"dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRWlepb_"+decaystring).c_str(),10,0, 5, "#Delta R(l_{W},b)");
      MSPlot[(prefix+"dPhiWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dPhiWlepb_"+decaystring).c_str(),10,-4, 4, "#Delta #Phi (l_{W},b)");
      MSPlot[(prefix+"Zboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"Zboson_pt_"+decaystring).c_str(), 20,0, 500, "Z boson p_{T} (GeV)", "GeV");
      MSPlot[(prefix+"dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZWlep_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,l_{W})");
      MSPlot[(prefix+"bdiscCSVv2_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"bdiscCSVv2_jet_0_"+decaystring).c_str(),20, 0.5, 1, "leading jet CSVv2");
      MSPlot[(prefix+"cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_0_"+decaystring).c_str(),10, 0, 0.8, "leading jet charm vs b disc.");
      
      if(isZut){
        MSPlot[(prefix+"charge_asym_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"charge_asym_"+decaystring).c_str(),10, -4, 4, "charge l_{W} x|W boson #eta|");
        
      }
      else{
        MSPlot[(prefix+"cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_0_"+decaystring).c_str(),25, 0, 1, "leading jet charm vs light disc.");
        
      }
      
    }
    else if(istoppair ){
      MSPlot[(prefix+"mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"mlb_"+decaystring).c_str(),25, 0, 500, "inv mass l_{W}b (GeV)","GeV");
      MSPlot[(prefix+"FCNCtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"FCNCtop_M_"+decaystring).c_str(),20, 100, 500, "inv. mass Zq (GeV)","GeV");
      MSPlot[(prefix+"dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRWlepb_"+decaystring).c_str(),20,0, 6, "#Delta R(b,l_{W})");
      MSPlot[(prefix+"nJets_CSVv2M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"nJets_CSVv2M_"+decaystring).c_str(),10,-0.5, 9.5, "# CSVv2M jets");
      MSPlot[(prefix+"nJets_CharmL_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"nJets_CharmL_"+decaystring).c_str(),10,-0.5, 9.5, "# charm loose jets");
      MSPlot[(prefix+"dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZWlep_"+decaystring).c_str(),20,0, 6, "#Delta R(Z,l_{W})");
      MSPlot[(prefix+"dRZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZc_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,q)");
      
      if(isZut){
        MSPlot[(prefix+"cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_0_"+decaystring).c_str(),20, 0, 0.8, "leading jet vharm vs b disc.");
        MSPlot[(prefix+"cdiscCvsB_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_1_"+decaystring).c_str(),20, 0, 0.85, "2nd leading charm vs b disc. ");
        
      }
      else{
        MSPlot[(prefix+"cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_0_"+decaystring).c_str(),20, 0, 1, "leading jet charm vs light disc.");
        MSPlot[(prefix+"cdiscCvsL_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_1_"+decaystring).c_str(),20, 0, 1, "2nd leading jet charm vs light disc.");
        
        
      }
    }*/
    
  }
  decaystring = "";
  
  
  
  
}
void InitCalculatePDFWeightHisto(string dataSetName,bool doMTWtemplate,bool doTTZtemplate){
  TH1::SetDefaultSumw2();
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  string template_name = "BDT";
  if(doMTWtemplate) template_name = "MTW";
  else if (doTTZtemplate) template_name = "Zmass";
  if(!doMTWtemplate ){
    for(Int_t iChan = 0; iChan < channel_list.size(); iChan++){
      channel = channel_list[iChan];
      for ( Int_t i=0; i<101; i++)
      {
        
        output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_"+intToStr(i);
      //  cout << output_histo_name << endl;
        histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,BDT_begin,BDT_end);
      }
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_nominal";
     // cout << output_histo_name << endl;
      histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,BDT_begin,BDT_end);
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_PDFEnvelopeUp";
      histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,BDT_begin,BDT_end);
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_PDFEnvelopeDown";
      histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,BDT_begin,BDT_end);
      
      output_histo_name = "";
    }
  }
  else   if(doMTWtemplate ){
    for(Int_t iChan = 0; iChan < channel_list.size(); iChan++){
      channel = channel_list[iChan];
      for ( Int_t i=0; i<101; i++)
      {
        
        output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_"+intToStr(i);
     //   cout << output_histo_name << endl;
        histo1DPDFMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0., endMTW);
      }
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_nominal";
      //cout << output_histo_name << endl;
      histo1DPDFMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0.,endMTW);
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_PDFEnvelopeUp";
      histo1DPDFMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0.,endMTW);
      output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_PDFEnvelopeDown";
      histo1DPDFMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0.,endMTW);
      
      output_histo_name = "";
    }
  }
}
void InitMTWShapeHisto(string dataSetName, string systematic, Int_t isys,  vector <int> decayChannels){
  TH1::SetDefaultSumw2();
  //histo1DMTW.clear();
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    
    if(isys == 0) output_histo_name = dataSetName+"_MTW_nominal_"+decaystring;
    else output_histo_name = dataSetName+"_MTW_"+systematic + "_" + decaystring;
    //cout << "init " << output_histo_name.c_str() << endl;
    histo1DMTW[output_histo_name.c_str()] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0.,endMTW);
    
    output_histo_name = "";
  }
  decaystring = "";
}
void InitSystematicHisto(string dataSetName, string systematic, Int_t isys, bool doMTWtemplate){
  TH1::SetDefaultSumw2();
  
  if(isys == 0 && !doMTWtemplate) output_histo_name = dataSetName+"_BDT_nominal";
  else if(isys == 0 && doMTWtemplate) output_histo_name = dataSetName+"_MTW_nominal";
  else if(!doMTWtemplate) output_histo_name = dataSetName+"_BDT_"+systematic;
  else if(doMTWtemplate) output_histo_name = dataSetName+"_MTW_"+systematic;
  
  if(!doMTWtemplate) histo1DSys[output_histo_name] = new TH1F(output_histo_name.c_str(), dataSetName.c_str(), nbin,BDT_begin,BDT_end);
  else if(doMTWtemplate) histo1DSysMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), dataSetName.c_str(), nbinMTW,0.,endMTW);
  output_histo_name = "";
  
}
void InitAnalyzerTree(TTree* tree){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress( "MVA_weight_nom", &MVA_weight_nom, &b_MVA_weight_nom);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress("MVA_EqLumi", &MVA_EqLumi, &b_MVA_EqLumi);
  tree->SetBranchAddress("MVA_Luminosity", &MVA_Luminosity, &b_MVA_Luminosity);
  tree->SetBranchAddress("MVA_mWt2", &MVA_mWt2, &b_MVA_mWt2);
  tree->SetBranchAddress("MVA_mWt", &MVA_mWt, &b_MVA_mWt);
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  
  tree->SetBranchAddress("MVA_x1", &MVA_x1, &b_MVA_x1);
  tree->SetBranchAddress("MVA_x2", &MVA_x2, &b_MVA_x2);
  tree->SetBranchAddress("MVA_id1", &MVA_id1, &b_MVA_id1);
  tree->SetBranchAddress("MVA_id2", &MVA_id2, &b_MVA_id2);
  tree->SetBranchAddress("MVA_q", &MVA_q, &b_MVA_q);
  
  tree->SetBranchAddress("MVA_DeltaR_NonIsoLepJet", &MVA_DeltaR_NonIsoLepJet, &b_MVA_DeltaR_NonIsoLepJet);
  tree->SetBranchAddress("MVA_DeltaR_lep0Jet", &MVA_DeltaR_lep0Jet, &b_MVA_DeltaR_lep0Jet);
  tree->SetBranchAddress("MVA_DeltaR_lep1Jet", &MVA_DeltaR_lep1Jet, &b_MVA_DeltaR_lep1Jet);
  tree->SetBranchAddress("MVA_DeltaR_lep2Jet", &MVA_DeltaR_lep2Jet, &b_MVA_DeltaR_lep2Jet);
  tree->SetBranchAddress("MVA_DeltaR_MinLepJet", &MVA_DeltaR_MinLepJet, &b_MVA_DeltaR_MinLepJet);
  tree->SetBranchAddress("MVA_DeltaR_MinLepNonIsoJet", &MVA_DeltaR_MinLepNonIsoJet, &b_MVA_DeltaR_MinLepNonIsoJet);
  tree->SetBranchAddress("MVA_Zboson_M", &MVA_Zboson_M, &b_MVA_Zboson_M);
  tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
  tree->SetBranchAddress( "MVA_nJets",&MVA_nJets, &b_MVA_nJets);
  
  tree->SetBranchAddress("MVA_weightWZcorr", &MVA_weightWZcorr, &b_MVA_weightWZcorr);
  tree->SetBranchAddress( "MVA_weight_nloSF", &MVA_weight_nloSF, &b_MVA_weight_nloSF);
  tree->SetBranchAddress( "MVA_weight_puSF_up", &MVA_weight_puSF_up, &b_MVA_weight_puSF_up);
  tree->SetBranchAddress( "MVA_weight_puSF_down", &MVA_weight_puSF_down, &b_MVA_weight_puSF_down);
  tree->SetBranchAddress( "MVA_weight_electronSF_up", &MVA_weight_electronSF_up, &b_MVA_weight_electronSF_up);
  tree->SetBranchAddress( "MVA_weight_electronSF_down", &MVA_weight_electronSF_down, &b_MVA_weight_electronSF_down);
  tree->SetBranchAddress( "MVA_weight_muonSF_up", &MVA_weight_muonSF_up, &b_MVA_weight_muonSF_up);
  tree->SetBranchAddress( "MVA_weight_muonSF_down", &MVA_weight_muonSF_down, &b_MVA_weight_muonSF_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_cferr1_up", &MVA_weight_btagSF_cferr1_up, &b_MVA_weight_btagSF_cferr1_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_cferr1_down", &MVA_weight_btagSF_cferr1_down, &b_MVA_weight_btagSF_cferr1_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_cferr2_up", &MVA_weight_btagSF_cferr2_up, &b_MVA_weight_btagSF_cferr2_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_cferr2_down", &MVA_weight_btagSF_cferr2_down, &b_MVA_weight_btagSF_cferr2_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_hf_up", &MVA_weight_btagSF_hf_up, &b_MVA_weight_btagSF_hf_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_hf_down", &MVA_weight_btagSF_hf_down, &b_MVA_weight_btagSF_hf_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_hfstats1_up", &MVA_weight_btagSF_hfstats1_up, &b_MVA_weight_btagSF_hfstats1_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_hfstats1_down", &MVA_weight_btagSF_hfstats1_down, &b_MVA_weight_btagSF_hfstats1_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_hfstats2_up", &MVA_weight_btagSF_hfstats2_up, &b_MVA_weight_btagSF_hfstats2_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_hfstats2_down", &MVA_weight_btagSF_hfstats2_down, &b_MVA_weight_btagSF_hfstats2_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_lf_up", &MVA_weight_btagSF_lf_up, &b_MVA_weight_btagSF_lf_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_lf_down", &MVA_weight_btagSF_lf_down, &b_MVA_weight_btagSF_lf_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_lfstats1_up", &MVA_weight_btagSF_lfstats1_up, &b_MVA_weight_btagSF_lfstats1_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_lfstats1_down", &MVA_weight_btagSF_lfstats1_down, &b_MVA_weight_btagSF_lfstats1_down);
  tree->SetBranchAddress( "MVA_weight_btagSF_lfstats2_up", &MVA_weight_btagSF_lfstats2_up, &b_MVA_weight_btagSF_lfstats2_up);
  tree->SetBranchAddress( "MVA_weight_btagSF_lfstats2_down", &MVA_weight_btagSF_lfstats2_down, &b_MVA_weight_btagSF_lfstats2_down);
  
  
  tree->SetBranchAddress( "MVA_weight_puSF",&MVA_weight_puSF, &b_MVA_weight_puSF);
  tree->SetBranchAddress( "MVA_weight_muonSF",&MVA_weight_muonSF, &b_MVA_weight_muonSF);
  tree->SetBranchAddress( "MVA_weight_electronSF",&MVA_weight_electronSF, &b_MVA_weight_electronSF);
  tree->SetBranchAddress( "MVA_weight_btagSF",&MVA_weight_btagSF, &b_MVA_weight_btagSF);
  
  
  
}
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  InitAnalyzerTree(tree);
  
  tree->SetBranchAddress("MVA_lepton2_pt", &MVA_lepton2_pt, &b_MVA_lepton2_pt);
  
  tree->SetBranchAddress("MVA_BDT", &MVA_BDT, &b_MVA_BDT);
  if(!istoppair && !doZut){
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_SMtop_eta", &MVA_SMtop_eta, &b_MVA_SMtop_eta);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_SMtop_rap", &MVA_SMtop_rap, &b_MVA_SMtop_rap);
    tree->SetBranchAddress("MVA_TotalInvMass_lep", &MVA_TotalInvMass_lep, &b_MVA_TotalInvMass_lep);
  }
  else if(!istoppair && doZut){
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_SMtop_rap", &MVA_SMtop_rap, &b_MVA_SMtop_rap);
    tree->SetBranchAddress("MVA_TotalInvMass_lep", &MVA_TotalInvMass_lep, &b_MVA_TotalInvMass_lep);
    tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
    tree->SetBranchAddress("MVA_ptWQ", &MVA_ptWQ, &b_MVA_ptWQ);
  }
  else if(istoppair && doZut){
    tree->SetBranchAddress("MVA_dPhiZWlep", &MVA_dPhiZWlep, &b_MVA_dPhiZWlep);
    tree->SetBranchAddress("MVA_SMtop_rap", &MVA_SMtop_rap, &b_MVA_SMtop_rap);
    //tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
    tree->SetBranchAddress("MVA_SMtop_M", &MVA_SMtop_M, &b_MVA_SMtop_M);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_TotalHt_lep", &MVA_TotalHt_lep, &b_MVA_TotalHt_lep);
    //tree->SetBranchAddress("MVA_Zboson_M", &MVA_Zboson_M, &b_MVA_Zboson_M);
    tree->SetBranchAddress("MVA_dPhiZb", &MVA_dPhiZb, &b_MVA_dPhiZb);
    tree->SetBranchAddress("MVA_dRZb", &MVA_dRZb, &b_MVA_dRZb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_FCNCtop_rap", &MVA_FCNCtop_rap, &b_MVA_FCNCtop_rap);
    tree->SetBranchAddress("MVA_FCNCtop_M", &MVA_FCNCtop_M, &b_MVA_FCNCtop_M);
  }
  else if(istoppair && !doZut){
    tree->SetBranchAddress("MVA_FCNCtop_rap", &MVA_FCNCtop_rap, &b_MVA_FCNCtop_rap);
    tree->SetBranchAddress("MVA_SMtop_rap", &MVA_SMtop_rap, &b_MVA_SMtop_rap);
   // tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_FCNCtop_M", &MVA_FCNCtop_M, &b_MVA_FCNCtop_M);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_TotalInvMass", &MVA_TotalInvMass, &b_MVA_TotalInvMass);
    tree->SetBranchAddress("MVA_Bdis_Lightjet", &MVA_Bdis_Lightjet, &b_MVA_Bdis_Lightjet);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_dPhiZWlep", &MVA_dPhiZWlep, &b_MVA_dPhiZWlep);
   // tree->SetBranchAddress("MVA_Zboson_M", &MVA_Zboson_M, &b_MVA_Zboson_M);
    tree->SetBranchAddress("MVA_dRZb", &MVA_dRZb, &b_MVA_dRZb);
    
  }
  
  
  
  
  
}
void Init1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels){
  TH1::SetDefaultSumw2();
  cout << "initialising MVA var histo" << endl;
  // control plots
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    output_histo_name = "BDT_"+dataSetName +"_"+decaystring+"_"+systematic;
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,BDT_begin,BDT_end);
    output_histo_name = "SMtoprap_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "SM top rapidity",10,-6,6);
    output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Inv. Mass l_{W}b (GeV)",10,0,400);
    output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta R(l_{W},Z)",15,0,6);
    
    
    if(!istoppair){
      output_histo_name = "bdiscCSVv2jet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "CSVv2 disc. highest p_{T} jet",20,0.5,1);
      output_histo_name = "dPhiWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta #phi(l_{W},b)",10,-4,4);
      output_histo_name = "TotalInvMasslep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Inv. Mass of the leptons (GeV)",20,100,500);
      
      if(isZut){
        output_histo_name = "chargeAsym_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Q(l_{W}|#eta of l_{W}|",6,-3,3);
        output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta R(l_{W},b)",10,0,5);
        output_histo_name = "ptWQ_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "p_{T}(l_{W})Q(l_{W}",10,-200,200);
        
        
      }
      if(!isZut){
        output_histo_name = "SMtopeta_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "SM top #eta",10,-6,6);
      }
      
      
      
    }
    
    else if(istoppair ){
      output_histo_name = "FCNCtoprap_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "FCNC top rapidity",10,-3,3);
      output_histo_name = "nJetsCSVv2M_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Nb. of CSVv2 M jets",7,-0.5,6.5);
      output_histo_name = "FCNCtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Inv. Mass FCNC top (GeV)",10,0,400);
      output_histo_name = "ZbosonM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Inv. Mass Zboson (GeV)",10,60,120);
      output_histo_name = "dRZb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta R(Z,b)",10,0,5);
      output_histo_name = "dPhiZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta #phi(l_{W},Z)",10,-4,4);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta R(l_{W},b)",10,0,5);
      
      
      if(doZut){
        output_histo_name = "dPhiZb_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "#Delta #phi(b,Z)",10,-4,4);
        output_histo_name = "SMtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Inv. Mass SM top (GeV)",10,0,400);
        output_histo_name = "TotalHtlep_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Total H_{T} of the leptons (GeV)",10,0,400);
        
      }
      if(!doZut){
        output_histo_name = "TotalInvMass_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "Total Inv. Mass (GeV)",10,200,1500);
        output_histo_name = "bdisclightjet_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), "CSVv2 disc. FCNC jet",20,0.5,1);
      }
    }
    
    
    
  }
  
  output_histo_name = "";
  decaystring = "";
}
////////// FUNCTIONS
vector<double> BDTCUT(string region, string coupling){
  
  cout << "Determine BDT cut " << endl;
  
  
  string bdtinput_name = placeOutputReading+"/Reader_" + coupling + "_" + region + ".root";
  TFile* bdt_file = TFile::Open( bdtinput_name.c_str(), "READ" );
  
  cout << bdtinput_name << endl;
  
  
  
  
  string histo_name = "";
  string template_fake_name = "";
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  
  Double_t cut_eee, cut_eeu, cut_uuu, cut_uue;
  cut_eee = cut_uue = cut_eeu = cut_uuu = -1.;
  for(Int_t ichan=0; ichan<channel_list.size(); ichan++)
  {
    string channel = channel_list[ichan];
    TH1F *h_sum_bkg(0), *h_sum_sig(0), *h_tmp(0);
    histo_name = "";
    if(channel.find("uuu") != std::string::npos || channel.find("eeu") != std::string::npos) {template_fake_name = "FakeMu";}
    else {template_fake_name = "FakeEl";}
    
    
    for(Int_t isample = 0; isample < datasets.size(); isample++)
    {
      string dataSetName = datasets[isample]->Name();
      //cout << dataSetName << endl;
      if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // no signal in data
      if(datasets[isample]->Name().find("fake")==std::string::npos && datasets[isample]->Name().find("FCNC")==std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << endl;
        h_tmp = 0;
        histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
        if(h_tmp->GetEntries() != 0){
          if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
          else {h_sum_bkg->Add(h_tmp);}}
      }
      else if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << " is signal "<< endl;
        h_tmp = 0;
        histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        if(h_tmp->GetEntries() != 0){
          if(h_sum_sig == 0) {h_sum_sig = (TH1F*) h_tmp->Clone();}
          else {h_sum_sig->Add(h_tmp);}
        }
      }
      else{
        h_tmp = 0;
        histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + template_fake_name;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<histo_name<<" : not found"<<endl;}
        else
        {
          h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
          //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
          if(h_tmp->GetEntries() != 0){
            if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
            else {h_sum_bkg->Add(h_tmp);}
          }
        }
        
      }
    }
    //cout << "h_sum_bkg " << h_sum_bkg << " h_sum_sig "<< h_sum_sig << endl;
    if(h_sum_bkg == 0 || h_sum_sig == 0) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; }
    //S+B histogram
    TH1F* h_total = (TH1F*) h_sum_bkg->Clone();
    h_total->Add(h_sum_sig);
    
    //Normalization
    h_total->Scale(1/h_total->Integral()); //Integral() : Return integral of bin contents in range [binx1,binx2] (inclusive !)
    h_sum_bkg->Scale(1/h_sum_bkg->Integral());
    h_sum_sig->Scale(1/h_sum_sig->Integral());
    
    Double_t sig_over_total = 100; //initialize to unreasonable value
    Int_t bin_cut = -1; //initialize to false value
    Int_t nofbins = h_total->GetNbinsX();
    
    for(Int_t ibin=nofbins; ibin>0; ibin--)
    {
      //Search the bin w/ lowest sig/total, while keeping enough bkg events (criterion needs to be optimized/tuned)
      if( (h_sum_sig->Integral(1, ibin) / h_total->Integral(1, ibin)) < sig_over_total && (h_sum_bkg->Integral(1, ibin) / h_sum_bkg->Integral()) >= 0.6 )
      {
        bin_cut = ibin;
        sig_over_total = h_sum_sig->Integral(1, bin_cut) / h_total->Integral(1,bin_cut);
      }
    }
    
    Double_t cut = h_total->GetBinLowEdge(bin_cut+1); //Get the BDT cut value to apply to create a BDT CR control tree
    Double_t min = TMath::Min(h_sum_bkg->GetMinimum(), h_sum_sig->GetMinimum());
    Double_t max = TMath::Max(h_sum_bkg->GetMaximum(), h_sum_sig->GetMaximum());
    //Create plot to represent the cut on BDT
    TCanvas* c = new TCanvas("c", "Signal vs Background");
    gStyle->SetOptStat(0);
    h_sum_bkg->GetYaxis()->SetRange(min,1.049*max);
    h_sum_bkg->GetXaxis()->SetTitle("D");
    h_sum_bkg->SetTitle("Signal vs Background");
    h_sum_bkg->SetLineColor(kBlue);
    h_sum_sig->SetLineColor(kRed-3);
    h_sum_bkg->SetLineWidth(3);
    h_sum_sig->SetLineWidth(3);
    //h_sum_sig->Scale(100);
    h_sum_bkg->Draw("HIST");
    h_sum_sig->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.7,0.75,0.88,0.85);
    leg->SetHeader("");
    leg->AddEntry(h_sum_sig,"Signal","L");
    leg->AddEntry(h_sum_bkg,"Background","L");
    leg->Draw();
    //Draw vertical line at cut value
    
    TLine* l = new TLine(cut,min, cut, 1.049*max);
    l->SetLineWidth(3);
    l->SetLineColor(921);
    l->Draw("");
    TString outputsaving = "CutPlots";
    mkdir(outputsaving,0777);
    outputsaving += "/"+region + "_" + coupling + "_";
    outputsaving += "Signal_Background_BDT_"+channel+".png";
    c->SaveAs(outputsaving.Data());
    
    //Cout some results
    cout<<"---------------------------------------"<<endl;
    cout<<"* Cut Value = "<<cut<<endl;
    cout<<"-> BDT_CR defined w/ all events inside bins [1 ; "<<bin_cut<<"] of the BDT distribution!"<<endl<<endl;
    cout<<"* Signal integral = "<<h_sum_sig->Integral(1, bin_cut)<<" / Total integral "<<h_total->Integral(1, bin_cut)<<endl;
    cout<<"Signal contamination in CR --> Sig/Total = "<<sig_over_total<<endl;
    cout<<"Bkg(CR) / Bkg(Total) = "<<h_sum_bkg->Integral(1,bin_cut) / h_sum_bkg->Integral()<<endl;
    cout<<"---------------------------------------"<<endl<<endl;
    
    //for(Int_t i=0; i<h_sig->GetNbinsX(); i++) {cout<<"bin content "<<i+1<<" = "<<h_sig->GetBinContent(i+1)<<endl;} //If want to verify that the signal is computed correctly
    delete c; delete leg; delete l;
    if(channel.find("uuu") != std::string::npos) cut_uuu = cut;
    else if(channel.find("eeu") != std::string::npos) cut_eeu = cut;
    else if(channel.find("eee") != std::string::npos) cut_eee = cut;
    else if(channel.find("uue") != std::string::npos) cut_uue = cut;
    
    delete h_sum_sig;
    delete h_sum_bkg;
    delete h_tmp;
  }
  bdt_file->Close();
  vector<double> v_return;
  v_return.push_back(cut_uuu);
  v_return.push_back(cut_uue);
  v_return.push_back(cut_eeu);
  v_return.push_back(cut_eee);
  return v_return;
  
}
void CalculatePDFWeight(string dataSetName, Double_t BDT, Double_t MVA_weight_nom, Int_t MVA_channel,bool doMTWtemplate,bool doTTZtemplate){
 // cout << "calculate pdf" << endl;
  //std::vector<double> pdfweights;
  //cout << "MVA channel " << MVA_channel << endl;
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  
  //PDF weights calculation
  LHAPDF::setVerbosity(0);
  string sbase = "";
  if(dataSetName.find("WZTo3LNu_amc")!=std::string::npos || dataSetName.find("TTZ")!=std::string::npos) sbase = "NNPDF30_nlo_nf_5_pdfas";
  else if(dataSetName.find("tZq")!=std::string::npos) sbase = "NNPDF30_nlo_nf_5_pdfas";
  else if(dataSetName.find("ZZTo4L")!=std::string::npos) sbase ="NNPDF30_nlo_as_0118";
  else if(dataSetName.find("NP")!=std::string::npos) sbase ="NNPDF23_nlo_as_0119"; 
  //cout << "base set " << sbase << endl;
   LHAPDF::PDFSet basepdfSet(sbase.c_str());
  //LHAPDF::PDFSet basepdfSet("NNPDF30_lo_as_0130"); // base from main MC // WZ
  //LHAPDF::PDFSet basepdfSet("NNPDF30_nlo_nf_5_pdfas"); // WZ amc at nlo
  LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_100"); // give the correct name see https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#PDF_uncertainties
  //LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_30");
  const LHAPDF::PDF* basepdf = basepdfSet.mkPDF(0);
  channel = channel_list[MVA_channel];
  //cout << "MVA_x1 "<< MVA_x1 <<" MVA_x2 "<< MVA_x2 <<" MVA_id1 "<< MVA_id1 <<" MVA_id2 "<< MVA_id2 <<" MVA_q "<< MVA_q << endl;
  
  string template_name = "BDT";
  if(doMTWtemplate) template_name = "MTW";
  else if (doTTZtemplate) template_name = "Zmass";
  
  for ( size_t i=0; i<newpdfSet.size(); i++)
  {
    const LHAPDF::PDF* newpdf = newpdfSet.mkPDF(i);
    
    Double_t weightpdf = LHAPDF::weightxxQ(MVA_id1, MVA_id2, MVA_x1, MVA_x2, MVA_q, *basepdf, *newpdf);
    output_histo_name = dataSetName+"_"+template_name +"_"+channel+"_"+intToStr(i);
    
    if(!doTTZtemplate && !doMTWtemplate) histo1DPDF[output_histo_name]->Fill(MVA_BDT,MVA_weight_nom*weightpdf);
    else if(doMTWtemplate) histo1DPDFMTW[output_histo_name]->Fill(MVA_mWt,MVA_weight_nom*weightpdf); 
    else if (doTTZtemplate) histo1DPDF[output_histo_name]->Fill(MVA_Zboson_M,MVA_weight_nom*weightpdf);

   /* 
    if(!doMTWtemplate && !doTTZtemplate) cout << "fill " << output_histo_name.c_str() << " with " << MVA_BDT << " " <<  weightpdf << endl;
    else  if(doMTWtemplate && !doTTZtemplate) cout << "fill " << output_histo_name.c_str() << " with " << MVA_mWt << " " <<  weightpdf << endl;
    else  if(doTTZtemplate) cout << "fill " << output_histo_name.c_str() << " with " << MVA_Zboson_M << " " <<  weightpdf << endl;
    //cout << "pdf weight " << weight << endl;
    //pdfweights.push_back(weight);
    */
    delete newpdf;
    
  }
  output_histo_name = dataSetName+"_"+template_name+"_"+channel+"_nominal";
  if(!doTTZtemplate && !doMTWtemplate) histo1DPDF[output_histo_name]->Fill(MVA_BDT,MVA_weight_nom);
  else if(doMTWtemplate) histo1DPDFMTW[output_histo_name]->Fill(MVA_mWt,MVA_weight_nom);
  else if (doTTZtemplate) histo1DPDF[output_histo_name]->Fill(MVA_Zboson_M,MVA_weight_nom);
  output_histo_name = "";
  delete basepdf;
  
  
}
void FillMTWPlots(Int_t d, string postfix, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm){
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  
  
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    decaystring += postfix;
    
    //if(datasets[d]->Name().find("data")!=std::string::npos) cout << "filling " << ("MTW_"+decaystring).c_str() << " with " << MVA_mWt2 << " " << weight_ << endl;
    MSPlotMTW[("MTW_"+decaystring).c_str()]->Fill(MVA_mWt , datasets[d], true, weight_);
    /*
    MSPlotMTW[("DeltaR_NonIsoLepJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_NonIsoLepJet, datasets[d], true, weight_);
    MSPlotMTW[("DeltaR_lep0Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep0Jet, datasets[d], true, weight_);
    MSPlotMTW[("DeltaR_lep1Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep1Jet, datasets[d], true, weight_);
    MSPlotMTW[("DeltaR_lep2Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep2Jet, datasets[d], true, weight_);
    //MSPlotMTW[("DeltaR_MinLepJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_MinLepJet, datasets[d], true, weight_);
   // MSPlotMTW[("DeltaR_MinLepNonIsoJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_MinLepNonIsoJet, datasets[d], true, weight_);*/
    
   // MSPlotMTW[("ZBOSON_M_"+decaystring).c_str()]->Fill(MVA_Zboson_M, datasets[d], true, weight_);
    //if(  83.5 < MVA_Zboson_M && MVA_Zboson_M <  98.5) MSPlotMTW[("ZBOSON_Mcut_"+decaystring).c_str()]->Fill(MVA_Zboson_M, datasets[d], true, weight_);
   }
  decaystring = "";
}
void FillTTZPlots(Int_t d, string postfix, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm){
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  
  
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    decaystring += postfix;
    
    //if(datasets[d]->Name().find("data")!=std::string::npos) cout << "filling " << ("MTW_"+decaystring).c_str() << " with " << MVA_mWt2 << " " << weight_ << endl;
    MSPlotTTZ[("TTZ_"+decaystring).c_str()]->Fill(MVA_Zboson_M , datasets[d], true, weight_);
      }
  decaystring = "";
}
void FillGeneralPlots(Int_t d, string prefix, vector <int> decayChannels, bool isZut , bool istoppair, Double_t weight_, Int_t MVA_channel, Float_t zbosonmass, Float_t njetscsvm, bool postfit){
  
  //cout << "fill plots" << endl;
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  

  double originalweight = weight_;
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)&& ((datasets[d]->Name()).find("fake")==std::string::npos)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    
    if(decayChannels[iChan] != -9){
    if( ((datasets[d]->Name()).find("fake")!=std::string::npos) && postfit && !istoppair && isZut){
      if(decayChannels[iChan]  == 3) weight_ = originalweight *0.0714286;
      else if(decayChannels[iChan]  == 2) weight_ = originalweight * 0.0714286;
      else if(decayChannels[iChan]  == 1) weight_ = originalweight *0.357143;
      else if(decayChannels[iChan]  == 0) weight_= originalweight *0.5;
    }
    else if( ((datasets[d]->Name()).find("fake")!=std::string::npos) && postfit && istoppair && isZut){
      if(decayChannels[iChan]  == 3) weight_ = originalweight *0.166667;
      else if(decayChannels[iChan]  == 2) weight_ = originalweight * 0.0833333;
      else if(decayChannels[iChan]  == 1) weight_ = originalweight *0.458333;
      else if(decayChannels[iChan]  == 0) weight_= originalweight *0.291667;
    }
    }
    
    //cout << "bdt " << MVA_BDT << " in " << (prefix+"_BDT_"+decaystring).c_str()<< endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()]->Fill(MVA_BDT , datasets[d], true, weight_);
    MSPlot[ (prefix+"_channel_"+decaystring).c_str()]->Fill(MVA_channel, datasets[d], true, weight_);
   // MSPlot[ (prefix+"_weight_"+decaystring).c_str()]->Fill(weight_, datasets[d], true, 1.);
    
    /*
    MSPlot[(prefix+"_DeltaR_NonIsoLepJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_NonIsoLepJet, datasets[d], true, weight_);
    //cout << (prefix+"_DeltaR_lep0Jet_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_DeltaR_lep0Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep0Jet, datasets[d], true, weight_);
    MSPlot[(prefix+"_DeltaR_lep1Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep1Jet, datasets[d], true, weight_);
    MSPlot[(prefix+"_DeltaR_lep2Jet_"+decaystring).c_str()]->Fill(MVA_DeltaR_lep2Jet, datasets[d], true, weight_);
   // MSPlot[(prefix+"_DeltaR_MinLepJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_MinLepJet, datasets[d], true, weight_);
   // MSPlot[(prefix+"_DeltaR_MinLepNonIsoJet_"+decaystring).c_str()]->Fill(MVA_DeltaR_MinLepNonIsoJet, datasets[d], true, weight_);
    
   // cout << (prefix+"_ZBOSON_M_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_ZBOSON_M_"+decaystring).c_str()]->Fill(MVA_Zboson_M, datasets[d], true, weight_);
    if( 83.5 < MVA_Zboson_M  &&  MVA_Zboson_M < 98.5) MSPlot[(prefix+"_ZBOSON_Mcut_"+decaystring).c_str()]->Fill(MVA_Zboson_M, datasets[d], true, weight_);
   
   //cout << (prefix+"_NJETSCSVM_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_NJETSCSVM_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M, datasets[d], true, weight_);
    if(MVA_NJets_CSVv2M > 0 ) MSPlot[(prefix+"_NJETSCSVMcut_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M, datasets[d], true, weight_);
*/
    /*
    if(!istoppair){
      MSPlot[(prefix+"_mlb_"+decaystring).c_str()] ->Fill(MVA_mlb, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRWlepb_"+decaystring).c_str()] ->Fill(MVA_dRWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_dPhiWlepb_"+decaystring).c_str()] ->Fill(MVA_dPhiWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZWlep_"+decaystring).c_str()] ->Fill(MVA_dRZWlep, datasets[d], true, weight_);
      MSPlot[(prefix+"_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0, datasets[d], true, weight_);
      MSPlot[(prefix+"_cdiscCvsB_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_0, datasets[d], true, weight_);
      
      if(isZut){
        MSPlot[(prefix+"_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym, datasets[d], true, weight_);
        
      }
      else{
        MSPlot[(prefix+"_cdiscCvsL_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_0, datasets[d], true, weight_);
        
      }
      
    }
    else if(istoppair ){
      MSPlot[(prefix+"_mlb_"+decaystring).c_str()] ->Fill(MVA_mlb, datasets[d], true, weight_);
      MSPlot[(prefix+"_FCNCtop_M_"+decaystring).c_str()] ->Fill(MVA_FCNCtop_M, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRWlepb_"+decaystring).c_str()] ->Fill(MVA_dRWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_nJets_CSVv2M_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M,  datasets[d], true, weight_);
      MSPlot[(prefix+"_nJets_CharmL_"+decaystring).c_str()] ->Fill(MVA_nJets_CharmL, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZWlep_"+decaystring).c_str()] ->Fill(MVA_dRZWlep, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZc_"+decaystring).c_str()] ->Fill(MVA_dRZc, datasets[d], true, weight_);
      
      if(isZut){
        MSPlot[(prefix+"_cdiscCvsB_jet_0_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsB_jet_0, datasets[d], true, weight_);
        MSPlot[(prefix+"_cdiscCvsB_jet_1_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsB_jet_1, datasets[d], true, weight_);
        
      }
      else{
        MSPlot[(prefix+"_cdiscCvsL_jet_0_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsL_jet_0, datasets[d], true, weight_);
        MSPlot[(prefix+"_cdiscCvsL_jet_1_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsL_jet_1, datasets[d], true, weight_);
        
        
      }
    }*/
  }
  decaystring = "";
}

void GetPDFEnvelope(string dataSetName,bool doMTWtemplate,bool doTTZtemplate){
  Double_t binContentMax = -1.;
  Double_t binContentMin = 1000000000000000.;
  
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  
  string template_name = "BDT";
  if(doMTWtemplate) template_name = "MTW";
  else if (doTTZtemplate) template_name = "Zmass";

  // loop over channels
  for(Int_t iChan = 0; iChan < channel_list.size(); iChan++){
    channel = channel_list[iChan];
    vector<double> bincontents;
    output_histo_name = dataSetName+"_"+ template_name+"_" + channel + "_nominal";
    //cout <<  output_histo_name << endl;
    //cout << histo1DPDF[output_histo_name]->GetTitle() << endl;
    // get nominal th1F
    TH1F* histo_nom;
    if(!doMTWtemplate) histo_nom = (TH1F*) histo1DPDF[output_histo_name]->Clone();
    if(doMTWtemplate) histo_nom = (TH1F*) histo1DPDFMTW[output_histo_name]->Clone();
    
    // loop over bins
    for( Int_t ibin = 1; ibin <histo_nom->GetNbinsX(); ibin++)
    {
      binContentMax = -1.;
      binContentMin = 1000000000000000.;
      bincontents.clear();
     // cout << "ibin " << ibin << endl;
      // get bincontents of each histo for this bin
      for(Int_t iCount = 0; iCount < 101; iCount++)
      {
        output_histo_name = dataSetName+"_"+template_name+"_"+channel + "_" +intToStr(iCount);
        
         if(doMTWtemplate ) cout << histo1DPDFMTW.size() << endl;
     //   cout <<" pushing "<< endl;
        if(!doMTWtemplate)bincontents.push_back(histo1DPDF[output_histo_name]->GetBinContent(ibin));
        if(doMTWtemplate ) bincontents.push_back(histo1DPDFMTW[output_histo_name]->GetBinContent(ibin));
       // cout << output_histo_name << endl;
      }
      if(binContentMin > minimumValue(bincontents)) binContentMin = minimumValue(bincontents);
      else binContentMin = histo_nom->GetBinContent(ibin);
      if(binContentMax < maximumValue(bincontents)) binContentMax = maximumValue(bincontents);
      else binContentMax = histo_nom->GetBinContent(ibin);
      
      output_histo_name = dataSetName+"_"+template_name+"_"+channel + "_PDFEnvelopeUp";
    // cout << output_histo_name << endl;
      if(!doMTWtemplate && !doTTZtemplate) histo1DPDF[output_histo_name]->SetBinContent(ibin, binContentMax);
      else if(doMTWtemplate ) histo1DPDFMTW[output_histo_name]->SetBinContent(ibin, binContentMax);
      output_histo_name = dataSetName+"_"+template_name+"_"+channel + "_PDFEnvelopeDown";
      //cout << output_histo_name << endl;
     if(!doMTWtemplate && !doTTZtemplate) histo1DPDF[output_histo_name]->SetBinContent(ibin, binContentMin);
      else if(doMTWtemplate ) histo1DPDFMTW[output_histo_name]->SetBinContent(ibin, binContentMin);
      output_histo_name = "";
    }// bins
  }//channels
  
}
void Fill1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels, Double_t weight_, Int_t MVA_channel){
  
  
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    output_histo_name = "BDT_"+dataSetName +"_"+decaystring+"_"+systematic;
    histo1D[output_histo_name] ->Fill(  MVA_BDT        ,weight_);
    
    output_histo_name = "SMtoprap_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name]->Fill(MVA_SMtop_rap,weight_);
    output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name]->Fill(MVA_mlb, weight_);
    output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
    histo1D[output_histo_name]->Fill(MVA_dRZWlep, weight_);
    
    
    if(!istoppair){
      output_histo_name = "bdiscCSVv2jet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_bdiscCSVv2_jet_0,weight_);
      output_histo_name = "dPhiWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_dPhiWlepb, weight_);
      output_histo_name = "TotalInvMasslep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_TotalInvMass_lep,weight_);
      
      if(isZut){
        output_histo_name = "chargeAsym_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_charge_asym,weight_);
        output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_dRWlepb,weight_);
        output_histo_name = "ptWQ_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(MVA_ptWQ,weight_);
        
        
      }
      if(!isZut){
        output_histo_name = "SMtopeta_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(MVA_SMtop_eta,weight_);
      }
      
      
      
    }
    
    else if(istoppair ){
      output_histo_name = "FCNCtoprap_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(MVA_FCNCtop_rap,weight_);
      output_histo_name = "nJetsCSVv2M_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(MVA_NJets_CSVv2M,weight_);
      output_histo_name = "FCNCtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_FCNCtop_M,weight_);
      output_histo_name = "ZbosonM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_Zboson_M,weight_);
      output_histo_name = "dRZb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(MVA_dRZb,weight_);
      output_histo_name = "dPhiZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name]->Fill(MVA_dPhiZWlep,weight_);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(MVA_dRWlepb,weight_);
      
      
      if(doZut){
        output_histo_name = "dPhiZb_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_dPhiZb,weight_);
        output_histo_name = "SMtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(MVA_SMtop_M,weight_);
        output_histo_name = "TotalHtlep_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_TotalHt_lep,weight_);
        
      }
      if(!doZut){
        output_histo_name = "TotalInvMass_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_TotalInvMass,weight_);
        output_histo_name = "bdisclightjet_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name]->Fill(MVA_Bdis_Lightjet,weight_);
      }
    }
    
    
    
    
    
  }
  
  output_histo_name = "";
  decaystring = "";
  
}
void FillMTWShapeHisto(string dataSetName, string systematic, Double_t weight_,Int_t isys, Int_t MVA_channel, vector <int> decayChannels){
  
  for(Int_t iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    if(isys == 0) output_histo_name = dataSetName+"_MTW_nominal_"+decaystring;
    else output_histo_name = dataSetName+"_MTW_"+systematic + "_" + decaystring;
    
    // cout << "fill " << output_histo_name << " " << MVA_mWt2 << " " << weight_ <<  endl;
    histo1DMTW[output_histo_name.c_str()]->Fill(MVA_mWt2, 1.);
    // histo1DMTW[output_histo_name.c_str()]->Fill(1, 1);
    output_histo_name = "";
  }
  decaystring = "";
}

void FillSystematicHisto(string dataSetName, string systematic, Double_t weight_, Int_t isys, bool doMTWtemplate ){
  
  if(isys == 0 && !doMTWtemplate) output_histo_name = dataSetName+"_BDT_nominal";
  else if(isys == 0 && doMTWtemplate) output_histo_name = dataSetName+"_MTW_nominal";
  else if(!doMTWtemplate) output_histo_name = dataSetName+"_BDT_"+systematic;
  else if(doMTWtemplate) output_histo_name = dataSetName+"_MTW_"+systematic;
  
  
  if(!doMTWtemplate) histo1DSys[output_histo_name]->Fill(MVA_BDT, weight_);
  else histo1DSysMTW[output_histo_name]->Fill(MVA_mWt2, weight_);
  
  
  output_histo_name = "";
  
}


void Ini1DHistoBDTvariavles(bool istoppair, bool isZut, vector<int> decayChannels){
  
}



