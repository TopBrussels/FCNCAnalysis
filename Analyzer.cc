#define TreeAnalyzer_cxx
//#include "TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include "TStyle.h"
#include "TPaveText.h"

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

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"


struct HighestPt
{
  bool operator()(TLorentzVector j1, TLorentzVector j2) const
  {
    return j1.Pt() > j2.Pt();
  }
  
  
};


using namespace std;
using namespace TopTree;



// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCutfl;


map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;


vector < Dataset* > datasets;

std::vector < int>  decayChannels = {0,1,2,3,-9}; // uuu uue eeu eee all
//std::vector < int>  decayChannels = {-9};
bool firstevent = false;

// Decleration of MVA variables
TFile * tupfile = 0;
TTree* mvatree = 0;
Float_t MVA_channel = -999.;
Float_t MVA_weight = 1.;
Float_t MVA_region = -999.;


Float_t MVA_lepton0_pt = -999.;
Float_t MVA_lepton1_pt= -999.;
Float_t MVA_lepton2_pt= -999.;
Float_t MVA_lepton0_eta= -999.;
Float_t MVA_lepton1_eta= -999.;
Float_t MVA_lepton2_eta= -999.;
Float_t MVA_lepton0_phi= -999.;
Float_t MVA_lepton1_phi= -999.;
Float_t MVA_lepton2_phi= -999.;

Float_t MVA_jet0_pt = -999.;
Float_t MVA_jet0_eta= -999.;
Float_t MVA_jet0_phi= -999.;
Float_t MVA_jet1_pt = -999.;
Float_t MVA_jet1_eta = -999.;
Float_t MVA_jet1_phi = -999.;



// SM side
Float_t MVA_Wlep_pt= -999.;
Float_t MVA_Wlep_eta= -999.;
Float_t MVA_Wlep_phi= -999.;
Float_t MVA_SMbjet_pt= -999.;
Float_t MVA_SMbjet_eta= -999.;
Float_t MVA_SMbjet_phi= -999.;
Float_t MVA_Wboson_pt= -999.;
Float_t MVA_Wboson_eta= -999.;
Float_t MVA_Wboson_phi= -999.;
Float_t MVA_met= -999.;
Float_t MVA_SMtop_pt= -999.;
Float_t MVA_SMtop_eta= -999.;
Float_t MVA_SMtop_phi= -999.;


// FCNC side
Float_t MVA_Zboson_pt= -999.;
Float_t MVA_Zboson_eta= -999.;
Float_t MVA_Zboson_phi= -999.;


Float_t MVA_LightJet_pt= -999.;
Float_t MVA_LightJet_eta= -999.;
Float_t MVA_LightJet_phi= -999.;

Float_t MVA_FCNCtop_pt= -999.;
Float_t MVA_FCNCtop_eta= -999.;
Float_t MVA_FCNCtop_phi= -999.;


// nbrs
Float_t MVA_nMuons = -999.;
Float_t MVA_NJets_CSVv2T = -999.;
Float_t MVA_NJets_CSVv2M = -999.;
Float_t MVA_NJets_CSVv2L = -999.;
Float_t MVA_nJets = -999.;
Float_t MVA_nElectrons = -999.;
Float_t MVA_nJets_CharmL = -999.;
Float_t MVA_nJets_CharmM = -999.;
Float_t MVA_nJets_CharmT = -999.;


//SM kinematics
Float_t MVA_mWt = -999.;
Float_t MVA_SMtop_M = -999.;
Float_t MVA_mlb = -999.;
Float_t MVA_Wboson_M = -999.;

Float_t MVA_dRWlepb = -999.;

Float_t MVA_dPhiWlepb = -999.;

Float_t MVA_Wlep_Charge = -999.;
Float_t MVA_charge_asym = -999.;
Float_t MVA_TotalPt = -999.;
Float_t MVA_TotalHt = -999.;
Float_t MVA_TotalInvMass = -999.;
Float_t MVA_TotalPt_lep = -999.;
Float_t MVA_TotalHt_lep = -999.;
Float_t MVA_TotalInvMass_lep = -999.;
Float_t MVA_TotalPt_jet = -999.;
Float_t MVA_TotalHt_jet = -999.;
Float_t MVA_TotalInvMass_jet = -999.;
Float_t MVA_bdiscCSVv2_jet_0 = -999.;
Float_t MVA_bdiscCSVv2_jet_1 = -999.;
Float_t MVA_CosTheta = -999.;
Float_t MVA_CosTheta_alt = -999.;



// FCNC kinematics
Float_t MVA_FCNCtop_M = -999.;
Float_t MVA_Zboson_M = -999.;

Float_t MVA_dRZc = -999.;
Float_t MVA_dPhiZc = -999.;

Float_t MVA_cdiscCvsB_jet_1 = -999.;
Float_t MVA_cdiscCvsL_jet_1 = -999.;
Float_t MVA_cdiscCvsB_jet_0 = -999.;
Float_t MVA_cdiscCvsL_jet_0 = -999.;

// interplay
Float_t MVA_dRSMFCNCtop = -999.;
Float_t MVA_dRZb = -999.;
Float_t MVA_dRWlepc = -999.;
Float_t MVA_dRZWlep = -999.;
Float_t MVA_dRZSMtop = -999.;

Float_t MVA_dPhiSMFCNCtop = -999.;
Float_t MVA_dPhiZb = -999.;
Float_t MVA_dPhiWlepc = -999.;
Float_t MVA_dPhiZWlep = -999.;
Float_t MVA_dPhiZMET = -999.;
Float_t MVA_dPhiZSMtop = -999.;

Float_t MVA_m3l = -999.;

// Declaration of leaf types

Int_t nEv;
Int_t nofNegWeights;
Int_t nofPosWeights;
Int_t sumW;


TBranch        *b_nEv;
TBranch *b_nofNegWeights;
TBranch *b_nofPosWeights;
TBranch *b_sumW;
int TotalEvents;
double Luminosity = 32000.;
Double_t WPb_L;
Double_t WPb_M;
Double_t WPb_T;
Double_t WPc_CvsB_Loose;
Double_t WPc_CvsB_Medium;
Double_t WPc_CvsB_Tight;
Double_t WPc_CvsL_Loose;
Double_t WPc_CvsL_Medium;
Double_t WPc_CvsL_Tight;
bool datafound = false;

TBranch *b_WPb_L;
TBranch *b_WPb_M;
TBranch *b_WPb_T;
TBranch *b_WPc_CvsB_Loose;
TBranch *b_WPc_CvsB_Medium;
TBranch *b_WPc_CvsB_Tight;
TBranch *b_WPc_CvsL_Loose;
TBranch *b_WPc_CvsL_Medium;
TBranch *b_WPc_CvsL_Tight;




Int_t nbEv_3lep_all;
Int_t nbEv_AtLeast1jet_all;
Int_t nbEv_Initial_all;
Int_t nbEv_LooseLepVeto_all;
Int_t nbEv_Trigged_all;
Int_t nbEv_ZbosonWindow_all;

TBranch *b_nbEv_3lep_all;
TBranch *b_nbEv_AtLeast1jet_all;
TBranch *b_nbEv_Initial_all;
TBranch *b_nbEv_LooseLepVeto_all;
TBranch *b_nbEv_Trigged_all;
TBranch *b_nbEv_ZbosonWindow_all;

Int_t nbEv_3lep_uuu;
Int_t nbEv_AtLeast1jet_uuu;
Int_t nbEv_Initial_uuu;
Int_t nbEv_LooseLepVeto_uuu;
Int_t nbEv_Trigged_uuu;
Int_t nbEv_ZbosonWindow_uuu;

TBranch *b_nbEv_3lep_uuu;
TBranch *b_nbEv_AtLeast1jet_uuu;
TBranch *b_nbEv_Initial_uuu;
TBranch *b_nbEv_LooseLepVeto_uuu;
TBranch *b_nbEv_Trigged_uuu;
TBranch *b_nbEv_ZbosonWindow_uuu;


Int_t nbEv_3lep_uue;
Int_t nbEv_AtLeast1jet_uue;
Int_t nbEv_Initial_uue;
Int_t nbEv_LooseLepVeto_uue;
Int_t nbEv_Trigged_uue;
Int_t nbEv_ZbosonWindow_uue;

TBranch *b_nbEv_3lep_eeu;
TBranch *b_nbEv_AtLeast1jet_eeu;
TBranch *b_nbEv_Initial_eeu;
TBranch *b_nbEv_LooseLepVeto_eeu;
TBranch *b_nbEv_Trigged_eeu;
TBranch *b_nbEv_ZbosonWindow_eeu;


Int_t nbEv_3lep_eeu;
Int_t nbEv_AtLeast1jet_eeu;
Int_t nbEv_Initial_eeu;
Int_t nbEv_LooseLepVeto_eeu;
Int_t nbEv_Trigged_eeu;
Int_t nbEv_ZbosonWindow_eeu;

TBranch *b_nbEv_3lep_uue;
TBranch *b_nbEv_AtLeast1jet_uue;
TBranch *b_nbEv_Initial_uue;
TBranch *b_nbEv_LooseLepVeto_uue;
TBranch *b_nbEv_Trigged_uue;
TBranch *b_nbEv_ZbosonWindow_uue;

Int_t nbEv_3lep_eee;
Int_t nbEv_AtLeast1jet_eee;
Int_t nbEv_Initial_eee;
Int_t nbEv_LooseLepVeto_eee;
Int_t nbEv_Trigged_eee;
Int_t nbEv_ZbosonWindow_eee;

TBranch *b_nbEv_3lep_eee;
TBranch *b_nbEv_AtLeast1jet_eee;
TBranch *b_nbEv_Initial_eee;
TBranch *b_nbEv_LooseLepVeto_eee;
TBranch *b_nbEv_Trigged_eee;
TBranch *b_nbEv_ZbosonWindow_eee;



// Declaration of leaf types
Int_t           channelInt;
Double_t        nloWeight;
Int_t           run_num;
Long64_t        evt_num;
Int_t           lumi_num;
Int_t           nvtx;
Int_t           npu;
Double_t        puSF;
Double_t        btagSF;
Int_t           PassedMETFilter;
Int_t           PassedGoodPV;
Int_t           nElectrons;
Double_t        ElectronSF[3];   //[nElectrons]
Float_t         pt_electron[3];   //[nElectrons]
Float_t         phi_electron[3];   //[nElectrons]
Float_t         eta_electron[3];   //[nElectrons]
Float_t         eta_superCluster_electron[3];   //[nElectrons]
Float_t         E_electron[3];   //[nElectrons]
Float_t         chargedHadronIso_electron[3];   //[nElectrons]
Float_t         neutralHadronIso_electron[3];   //[nElectrons]
Float_t         photonIso_electron[3];   //[nElectrons]
Float_t         pfIso_electron[3];   //[nElectrons]
Int_t           charge_electron[3];   //[nElectrons]
Float_t         d0_electron[3];   //[nElectrons]
Float_t         d0BeamSpot_electron[3];   //[nElectrons]
Float_t         sigmaIEtaIEta_electron[3];   //[nElectrons]
Float_t         deltaEtaIn_electron[3];   //[nElectrons]
Float_t         deltaPhiIn_electron[3];   //[nElectrons]
Float_t         hadronicOverEm_electron[3];   //[nElectrons]
Int_t           missingHits_electron[3];   //[nElectrons]
Bool_t          passConversion_electron[3];   //[nElectrons]
Bool_t          isId_electron[3];   //[nElectrons]
Bool_t          isIso_electron[3];   //[nElectrons]
Bool_t          isEBEEGap[3];   //[nElectrons]
Int_t           nMuons;
Double_t        MuonIDSF[3];   //[nMuons]
Double_t        MuonIsoSF[3];   //[nMuons]
Double_t        MuonTrigSFv2[3];   //[nMuons]
Double_t        MuonTrigSFv3[3];   //[nMuons]
Float_t         pt_muon[3];   //[nMuons]
Float_t         phi_muon[3];   //[nMuons]
Float_t         eta_muon[3];   //[nMuons]
Float_t         E_muon[3];   //[nMuons]
Float_t         chargedHadronIso_muon[3];   //[nMuons]
Float_t         neutralHadronIso_muon[3];   //[nMuons]
Float_t         photonIso_muon[3];   //[nMuons]
Bool_t          isId_muon[3];   //[nMuons]
Bool_t          isIso_muon[3];   //[nMuons]
Float_t         pfIso_muon[3];   //[nMuons]
Int_t           charge_muon[3];   //[nMuons]
Float_t         d0_muon[3];   //[nMuons]
Float_t         d0BeamSpot_muon[3];   //[nMuons]
Int_t           nJets;
Float_t         pt_jet[7];   //[nJets]
Float_t         px_jet[7];   //[nJets]
Float_t         py_jet[7];   //[nJets]
Float_t         pz_jet[7];   //[nJets]
Float_t         phi_jet[7];   //[nJets]
Float_t         eta_jet[7];   //[nJets]
Float_t         E_jet[7];   //[nJets]
Int_t           charge_jet[7];   //[nJets]
Float_t         bdisc_jet[7];   //[nJets]
Float_t         jet_Pt_before_JER[7];   //[nJets]
Float_t         jet_Pt_before_JES[7];   //[nJets]
Float_t         jet_Pt_after_JER[7];   //[nJets]
Float_t         jet_Pt_after_JES[7];   //[nJets]
Float_t         cdiscCvsL_jet[7];   //[nJets]
Float_t         cdiscCvsB_jet[7];   //[nJets]
Float_t         met_Pt;
Float_t         met_Eta;
Float_t         met_Phi;
Float_t         met_Px;
Float_t         met_Py;
Float_t         met_before_JES;
Float_t         met_after_JES;

// List of branches
TBranch        *b_channelInt;   //!
TBranch        *b_nloWeight;   //!
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nvtx;   //!
TBranch        *b_npu;   //!
TBranch        *b_puSF;   //!
TBranch        *b_btagSF;   //!
TBranch        *b_PassedMETFilter;   //!
TBranch        *b_PassedGoodPV;   //!
TBranch        *b_nElectrons;   //!
TBranch        *b_ElectronSF;   //!
TBranch        *b_pt_electron;   //!
TBranch        *b_phi_electron;   //!
TBranch        *b_eta_electron;   //!
TBranch        *b_eta_superCluster_electron;   //!
TBranch        *b_E_electron;   //!
TBranch        *b_chargedHadronIso_electron;   //!
TBranch        *b_neutralHadronIso_electron;   //!
TBranch        *b_photonIso_electron;   //!
TBranch        *b_pfIso_electron;   //!
TBranch        *b_charge_electron;   //!
TBranch        *b_d0_electron;   //!
TBranch        *b_d0BeamSpot_electron;   //!
TBranch        *b_sigmaIEtaIEta_electron;   //!
TBranch        *b_deltaEtaIn_electron;   //!
TBranch        *b_deltaPhiIn_electron;   //!
TBranch        *b_hadronicOverEm_electron;   //!
TBranch        *b_missingHits_electron;   //!
TBranch        *b_passConversion_electron;   //!
TBranch        *b_isId_electron;   //!
TBranch        *b_isIso_electron;   //!
TBranch        *b_isEBEEGap;   //!
TBranch        *b_nMuons;   //!
TBranch        *b_MuonIDSF;   //!
TBranch        *b_MuonIsoSF;   //!
TBranch        *b_MuonTrigSFv2;   //!
TBranch        *b_MuonTrigSFv3;   //!
TBranch        *b_pt_muon;   //!
TBranch        *b_phi_muon;   //!
TBranch        *b_eta_muon;   //!
TBranch        *b_E_muon;   //!
TBranch        *b_chargedHadronIso_muon;   //!
TBranch        *b_neutralHadronIso_muon;   //!
TBranch        *b_photonIso_muon;   //!
TBranch        *b_isId_muon;   //!
TBranch        *b_isIso_muon;   //!
TBranch        *b_pfIso_muon;   //!
TBranch        *b_charge_muon;   //!
TBranch        *b_d0_muon;   //!
TBranch        *b_d0BeamSpot_muon;   //!
TBranch        *b_nJets;   //!
TBranch        *b_pt_jet;   //!
TBranch        *b_px_jet;   //!
TBranch        *b_py_jet;   //!
TBranch        *b_pz_jet;   //!
TBranch        *b_phi_jet;   //!
TBranch        *b_eta_jet;   //!
TBranch        *b_E_jet;   //!
TBranch        *b_charge_jet;   //!
TBranch        *b_bdisc_jet;   //!
TBranch        *b_jet_Pt_before_JER;   //!
TBranch        *b_jet_Pt_before_JES;   //!
TBranch        *b_jet_Pt_after_JER;   //!
TBranch        *b_jet_Pt_after_JES;   //!
TBranch        *b_cdiscCvsL_jet;   //!
TBranch        *b_cdiscCvsB_jet;   //!
TBranch        *b_met_Pt;   //!
TBranch        *b_met_Eta;   //!
TBranch        *b_met_Phi;   //!
TBranch        *b_met_Px;   //!
TBranch        *b_met_Py;   //!
TBranch        *b_met_before_JES;   //!
TBranch        *b_met_after_JES;   //!





int verbose = 2;


string MakeTimeStamp();
void InitMSPlots(string prefix, vector<int> decayChannels);
void InitMVAMSPlotsSingletop(string prefix,vector<int> decayChannels);
void InitMVAMSPlotsTopPair(string prefix, vector<int> decayChannels);
void InitMVAMSPlotsWZ(string prefix, vector <int> decayChannels);
void Init1DPlots();
void Init2DPlots();
void InitTree(TTree* tree, bool isData);
// data from global tree
void ClearMetaData();
void GetMetaData(TTree* tree, bool isData,int Entries, bool isAMC);
// put everything to default values
void ClearObjects();
void ClearVars();
void ClearMVAVars();
void ClearTLVs();
void FillGeneralPlots(int d, string prefix, vector<int>decayChannels, bool isData);
void FillMVAPlots(int d, string dataSetName, int Region, string prefix, vector<int>decayChannels);
string ConvertIntToString(int nb, bool pad);
void ReconstructObjects(vector<TLorentzVector> Muons,vector<TLorentzVector> selectedElectrons, vector<TLorentzVector> selectedJets,int Region);
void MakeMVAvars(int Region, double scaleFactor);
void createMVAtree(string dataSetName);
void writeMVAtree();
int SMjetCalculator(vector<TLorentzVector> Jets,int verb);
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met);
int FCNCjetCalculator(vector < TLorentzVector>  Jets, TLorentzVector recoZ ,int index, int verb);
int FCNCjetCalculatorCvsBTagger(vector < TLorentzVector>  Jets, int index, int verb);
int FCNCjetCalculatorCvsLTagger(vector < TLorentzVector>  Jets, int index, int verb);
int FCNCjetCalculatorCwp(vector < TLorentzVector>  Jets, std::vector <int> cjetindex, int index, int verb);
std::pair <Float_t,Float_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo);
void LeptonAssigner(vector<TLorentzVector> electrons, vector<TLorentzVector> muons, std::vector<int> electronsCharge,std::vector<int> muonsCharge);

vector<TLorentzVector> selectedMuons;
vector<TLorentzVector> selectedElectrons;
vector<TLorentzVector> selectedLeptons;
vector<TLorentzVector> selectedJets;
vector<int> selectedElectronsCharge;
vector<int> selectedMuonsCharge;
bool Assigned = false;
TLorentzVector temp0;
TLorentzVector temp1;
TLorentzVector temp2;
TLorentzVector electron;
TLorentzVector jet;
TLorentzVector muon;
TLorentzVector Wboson;
TLorentzVector metTLV;
TLorentzVector metTLVbf;
TLorentzVector Zboson;
TLorentzVector SMtop;
TLorentzVector FCNCtop;
TLorentzVector SMbjet;
TLorentzVector LightJet;
TLorentzVector Wlep;
vector <int> selectednonCSVLJetID;
vector <int>  selectednonCSVMJetID;
vector <int>  selectednonCSVTJetID;
vector <int>  selectednonCvsBLJetID;
vector <int>  selectednonCvsBMJetID;
vector <int>  selectednonCvsBTJetID;
vector <int> selectednonCvsLLJetID;
vector <int> selectednonCvsLMJetID;
vector <int> selectednonCvsLTJetID;
vector <int> selectedCSVLJetID;
vector <int> selectedCSVMJetID;
vector <int> selectedCSVTJetID;
vector <int> selectedCvsBLJetID;
vector <int> selectedCvsBMJetID;
vector <int> selectedCvsBTJetID;
vector <int> selectedCvsLLJetID;
vector <int> selectedCvsLMJetID;
vector <int> selectedCvsLTJetID;
vector <int> selectedCharmLJetsindex;
vector <int> selectedCharmMJetsindex;
vector <int> selectedCharmTJetsindex;
double scaleFactor;
double muonSFtemp;
double electronSFtemp;
double EquilumiSF = 1.;
double nloSF = 1.;
double Xsect = 1.;
float mWT;
Float_t mWT2;
int globalnEntries;
int nEntries;
int SMjetIndex;
int Region;
int cjetindex;
int cjetindex_CvsLtagger;
int cjetindex_CvsBtagger;
int cjetindex_Cloose ;
int cjetindex_Cmedium;
int cjetindex_Ctight ;

Int_t           WmuIndiceF;
Int_t           WelecIndiceF;
Int_t           ZelecIndiceF_0;
Int_t           ZelecIndiceF_1;
Int_t           ZmuIndiceF_0;
Int_t           ZmuIndiceF_1;

int total_nbEv_3lep_all= 0;
int total_nbEv_AtLeast1jet_all= 0;
int total_nbEv_Initial_all= 0;
int total_nbEv_LooseLepVeto_all= 0;
int total_nbEv_Trigged_all= 0;
int total_nbEv_ZbosonWindow_all= 0;

int total_nbEv_3lep_uuu= 0;
int total_nbEv_AtLeast1jet_uuu= 0;
int total_nbEv_Initial_uuu= 0;
int total_nbEv_LooseLepVeto_uuu= 0;
int total_nbEv_Trigged_uuu= 0;
int total_nbEv_ZbosonWindow_uuu= 0;
int total_nbEv_3lep_uue= 0;
int total_nbEv_AtLeast1jet_uue= 0;
int total_nbEv_Initial_uue= 0;
int total_nbEv_LooseLepVeto_uue= 0;
int total_nbEv_Trigged_uue= 0;
int total_nbEv_ZbosonWindow_uue= 0;
int total_nbEv_3lep_eeu= 0;
int total_nbEv_AtLeast1jet_eeu= 0;
int total_nbEv_Initial_eeu= 0;
int total_nbEv_LooseLepVeto_eeu= 0;
int total_nbEv_Trigged_eeu= 0;
int total_nbEv_ZbosonWindow_eeu= 0;
int total_nbEv_3lep_eee= 0;
int total_nbEv_AtLeast1jet_eee= 0;
int total_nbEv_Initial_eee= 0;
int total_nbEv_LooseLepVeto_eee= 0;
int total_nbEv_Trigged_eee= 0;
int total_nbEv_ZbosonWindow_eee= 0;



Float_t tempPx;
Float_t tempPy;
Float_t tempHt;
Float_t tempPx_jet;
Float_t tempPy_jet;
Float_t tempHt_jet;

TLorentzVector tempInvMassObj;
TLorentzVector tempInvMassObj_jet;

string pathOutputdate = "";

int main(int argc, char* argv[])
{
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  setTDRStyle();
  
  string xmlFileName = "";
  xmlFileName = "config/Run2TriLepton_samples_analy.xml" ;
  const char* xmlFile = xmlFileName.c_str();
  cout << " - Using config file " << xmlFile << endl;
  
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  pathOutputdate = pathOutput + dateString + "/"  ;
  mkdir(pathOutputdate.c_str(),0777);
  
  
  
  //  load datasets
  datasets.clear();

  
  
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasets, xmlFile);
  cout << "datasets loaded" <<endl;
  bool makePlots = false;
  bool makeMVAtree = false;
  bool applyMuonSF = false;
  bool applyElectronSF = false;
  bool applyBTagSF = false; // To fix
  bool applyPUSF = false;
  bool applyNloSF = false;
  bool applyMETfilter = false;
  string placeNtup = "singletop/170214";
  int channel = -999;
  datafound = false;
  
  
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./Analyzer [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   decay uuu : add decay channel uuu, default value is all decays: uuu / uue / eeu / eee / all = 0 / 1 / 2 / 3 / -9" << endl;
      std::cout << "   debug: set verbose to high value" << endl;
      std::cout << "   applySF: applyPUSF, applyElectronSF, applyMuonsSF, applyMET, applyBtagSF, applyAMC. These can also be set seperatly" << endl;
      std::cout << "   MakeMVATree: make trees for the MVA" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. NtupleMakerOutput/MergedTuples/placeNtup " << endl;
      return 0;
    }
    if(string(argv[i]).find("Ntup")!=std::string::npos) {
      placeNtup = argv[i+1];
      i++;
      
    }
    if(string(argv[i]).find("decay")!=std::string::npos) {
      channel = strtol(argv[i+1], NULL, 10);
      i++;
      decayChannels.clear();
      decayChannels.push_back(channel);
    }
    if(string(argv[i]).find("debug")!=std::string::npos) {
      verbose = 4;
    }
    if(string(argv[i]).find("applyPUSF")!=std::string::npos) {
      applyPUSF =true;
    }
    if(string(argv[i]).find("applySF")!=string::npos) {
      applyElectronSF =true;
      applyNloSF = true;
      applyMETfilter = true;
      applyMuonSF = true;
      applyPUSF = true;
      applyBTagSF = true;
      cout << "                applying scalefactors " << endl;
    }
    if(string(argv[i]).find("applyElectronSF")!=string::npos) {
      applyElectronSF =true;
    }
    if(string(argv[i]).find("applyMuonSF")!=string::npos) {
      applyMuonSF =true;
    }
    if(string(argv[i]).find("applyMET")!=string::npos) {
      applyMETfilter =true;
    }
    if(string(argv[i]).find("applyBtagSF")!=string::npos) {
      applyBTagSF =true;
    }
    if(string(argv[i]).find("applyAMC")!=string::npos) {
      applyNloSF =true;
    }
    if(string(argv[i]).find("MakeMVATree")!=string::npos) {
      makeMVAtree = true;
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    
  }
  
  if(makePlots){
    firstevent = true;
    InitMSPlots("control_afterAtLeast1Jet", decayChannels);
    InitMSPlots("control_afterAtLeast1Jet_afterZWindow", decayChannels);
    InitMSPlots("control_afterAtLeast1Jet_afterZWindow_afterAtLeast1BJet", decayChannels);
    Init1DPlots();
    Init2DPlots();
    if(makeMVAtree) {
      InitMVAMSPlotsSingletop("singletop", decayChannels);
      //InitMVAMSPlotsTopPair("toppair", decayChannels);
      //InitMVAMSPlotsWZ("wzcontrol", decayChannels);
    }
  }
  
  
  

  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      Luminosity = datasets[d]->EquivalentLumi();
      datafound = true;
    }

    
  }
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  bool isData = false;
  bool isAMC = false;
  int nSelectedEntriesST = 0;
  int nSelectedEntriesTT = 0;
  int nSelectedEntriesWZ = 0;
  Double_t nSelectedEntriesSTweighted = 0.;
  double nSelectedEntriesTTweighted = 0.;
  double nSelectedEntriesWZweighted = 0.;
  /// Loop over datasets
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    
    firstevent = true;
    ClearMetaData();
    //cout << "meta data cleared" << endl;
    dataSetName = datasets[d]->Name();
    Xsect = datasets[d]->Xsection();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false;
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data")!= std::string::npos || dataSetName.find("DATA")!= std::string::npos )
    {
      isData = true;
    }
    isAMC= false;
    if ( dataSetName.find("amc")!= std::string::npos || dataSetName.find("AMC") != std::string::npos  )
    {
      isAMC = true;
      //cout << "amc at nlo sample" <<endl;
    }
    
    
    string ntupleFileName = "NtupleMakerOutput/MergedTuples/"+placeNtup+"/"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "globaltree";
    
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    globalnEntries = (int) tStatsTree[dataSetName.c_str()]->GetEntries();
    //cout << "getting meta data " << endl;
    GetMetaData(tStatsTree[dataSetName.c_str()], isData, globalnEntries, isAMC);
    //cout << "meta data gotten" << endl;
    
    
    /// Get data
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries: " << nEntries << endl;
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], isData);
    
    if(makeMVAtree){
      
      TString output_file_name = pathOutputdate+"/MVAtrees/";
      mkdir(output_file_name, 0777);
      output_file_name = pathOutputdate+"/MVAtrees/MVA_tree_" + dataSetName+ ".root";
      cout << "                making MVA tree file with name  " << output_file_name << endl;
      tupfile = new TFile(output_file_name,"RECREATE");
      mvatree = new TTree("mvatree", "mvatree");
      createMVAtree(dataSetName);
      
    }
    
    int endEvent = nEntries;
    int istartevt = 0;
    nSelectedEntriesST = 0;
    nSelectedEntriesTT = 0;
    nSelectedEntriesWZ = 0;
    nSelectedEntriesSTweighted = 0.;
    nSelectedEntriesTTweighted = 0.;
    nSelectedEntriesWZweighted = 0.;
    

    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects(); // put everything to default values
      
      if(ievt == 0 ){ firstevent = true;}
      else{ firstevent = false;}
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      
      if(applyMETfilter && !PassedMETFilter) continue;
      

      // fill objects,
      clock_t startFill = clock();
      muon.Clear();
      tempPx = 0.;
      tempPy = 0.;
      tempHt = 0.;
      tempInvMassObj.SetPtEtaPhiE(0.,0., 0.,0.);
      for(unsigned int iMu = 0; iMu < nMuons ; iMu++){
        muon.SetPtEtaPhiE(pt_muon[iMu], eta_muon[iMu], phi_muon[iMu], E_muon[iMu]);
        selectedMuons.push_back(muon);
        selectedMuonsCharge.push_back(charge_muon[iMu]);
        selectedLeptons.push_back(muon);
        tempPx = tempPx + muon.Px();
        tempPy = tempPy + muon.Py();
        tempHt = tempHt + muon.Pt();
        tempInvMassObj = tempInvMassObj + muon;
      }
      electron.Clear();
      for(unsigned int iEl = 0; iEl < nElectrons ; iEl++){
        electron.SetPtEtaPhiE(pt_electron[iEl], eta_electron[iEl], phi_electron[iEl], E_electron[iEl]);
        selectedElectrons.push_back(electron);
        selectedElectronsCharge.push_back(charge_electron[iEl]);
        selectedLeptons.push_back(electron);
        tempPx = tempPx + electron.Px();
        tempPy = tempPy + electron.Py();
        tempHt = tempHt + electron.Pt();
        tempInvMassObj = tempInvMassObj + electron;
      }
      sort(selectedLeptons.begin(), selectedLeptons.end(), HighestPt());
      MVA_TotalHt_lep = tempHt,
      MVA_TotalPt_lep = sqrt(tempPx*tempPx + tempPy*tempPx);
      MVA_TotalInvMass_lep = tempInvMassObj.M();
      jet.Clear();
      tempHt_jet = 0.;
      tempPy_jet = 0.;
      tempPx_jet = 0.;
      tempInvMassObj_jet.Clear();
      for(unsigned int iJet = 0; iJet < nJets ; iJet++){
        jet.SetPtEtaPhiE(pt_jet[iJet], eta_jet[iJet], phi_jet[iJet], E_jet[iJet]);
        selectedJets.push_back(jet);
        tempPx = tempPx + jet.Px();
        tempPy = tempPy + jet.Py();
        tempHt = tempHt + jet.Pt();
        tempInvMassObj = tempInvMassObj + jet;
        
        tempPx_jet = tempPx_jet + jet.Px();
        tempPy_jet = tempPy_jet + jet.Py();
        tempHt_jet = tempHt_jet + jet.Pt();
        tempInvMassObj_jet = tempInvMassObj_jet + jet;
        
        if(bdisc_jet[iJet] >  WPb_L) selectedCSVLJetID.push_back(iJet);
        else selectednonCSVLJetID.push_back(iJet);
        if(bdisc_jet[iJet] >  WPb_M) selectedCSVMJetID.push_back(iJet);
        else selectednonCSVMJetID.push_back(iJet);
        if(bdisc_jet[iJet] >  WPb_T) selectedCSVTJetID.push_back(iJet);
        else selectednonCSVTJetID.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Loose) selectedCvsBLJetID.push_back(iJet);
        else selectednonCvsBLJetID.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Medium) selectedCvsBMJetID.push_back(iJet);
        else selectednonCvsBMJetID.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Tight) selectedCvsBTJetID.push_back(iJet);
        else selectednonCvsBTJetID.push_back(iJet);
        if(cdiscCvsL_jet[iJet] > WPc_CvsL_Loose) selectedCvsLLJetID.push_back(iJet);
        else selectednonCvsLLJetID.push_back(iJet);
        if(cdiscCvsL_jet[iJet] > WPc_CvsL_Medium) selectedCvsLMJetID.push_back(iJet);
        else selectednonCvsLMJetID.push_back(iJet);
        if(cdiscCvsL_jet[iJet] > WPc_CvsL_Tight) selectedCvsLTJetID.push_back(iJet);
        else selectednonCvsLTJetID.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Loose && cdiscCvsL_jet[iJet] > WPc_CvsL_Loose) selectedCharmLJetsindex.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Medium && cdiscCvsL_jet[iJet] > WPc_CvsL_Medium) selectedCharmMJetsindex.push_back(iJet);
        if(cdiscCvsB_jet[iJet] > WPc_CvsB_Tight && cdiscCvsL_jet[iJet] > WPc_CvsL_Tight) selectedCharmTJetsindex.push_back(iJet);
        
        
      }
      
      MVA_TotalHt = tempHt + met_Pt;
      MVA_TotalPt = sqrt(tempPx*tempPx + tempPy*tempPx);
      MVA_TotalInvMass = tempInvMassObj.M();
      
      MVA_TotalHt_jet = tempHt_jet + met_Pt;
      MVA_TotalPt_jet = sqrt(tempPx_jet*tempPx_jet + tempPy_jet*tempPx_jet);
      MVA_TotalInvMass_jet = tempInvMassObj_jet.M();
      
      double time_fill = ((double)clock() - startFill) / CLOCKS_PER_SEC;
      if(firstevent && verbose > 3){
        cout << "It took us " << time_fill << " s to fill the vectors" << endl;
        if ( time_fill >= 60 )
        {
          int mins = time_fill/60;
          float secs = time_fill - mins*60;
          
          if (mins >= 60 )
          {
            int hours = mins/60;
            mins = mins - hours*60;
            cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
          }
          else
            cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
        }
      }
      
      //cout << "in assigner" <<endl;
      
      //cout << "WmuIndiceF " << WmuIndiceF <<" WelecIndiceF "<< WelecIndiceF <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
      // to be sure
      if((nElectrons + nMuons) != 3) continue;
      if(selectedJets.size() == 0) continue;
      
      LeptonAssigner(selectedElectrons, selectedMuons,selectedElectronsCharge ,selectedMuonsCharge);
      //cout << "WmuIndiceF " << WmuIndiceF <<" WelecIndiceF "<< WelecIndiceF <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
      if(!Assigned) continue;
      //cout << "in reco" << endl;
      ReconstructObjects(selectedMuons, selectedElectrons, selectedJets, Region);
      
      // apply SF
      scaleFactor = 1.;
      if (! isData)
      {
        if (applyMuonSF) {
          if(ievt == 2)cout << "                - applying muon factors " << endl;
          muonSFtemp = 1.;
          for(unsigned int iMu = 0; iMu < nMuons ; iMu++){
            scaleFactor *= MuonIDSF[iMu] * MuonIsoSF[iMu] ; //* (fracHLT[0]*muonTrigSFv2[0] + fracHLT[1]*muonTrigSFv3[0])
            muonSFtemp *= MuonIDSF[iMu] * MuonIsoSF[iMu] ;
            //cout << "iMu " << iMu << " MuonIDSF[iMu] " << MuonIDSF[0] << " MuonIsoSF[iMu] " << MuonIsoSF[0] << endl;
          }
          //cout << "muon  SF " << muonSFtemp << endl;
        }
        else muonSFtemp = 1.;
        
        if (applyElectronSF) {
          if(ievt == 2)cout << "                - applying electron factors " << endl;
          electronSFtemp = 1.;
          for(unsigned int iEl = 0; iEl < nElectrons ; iEl++){
            scaleFactor *= ElectronSF[iEl] ;
            electronSFtemp *= ElectronSF[iEl];
          }
        }
        else electronSFtemp = 1.;
        
        if (applyPUSF) {
          scaleFactor *= puSF;
          if(ievt == 2)cout << "                - applying pu factors " << endl;
        }
        else puSF  =1.;
        
        if (applyBTagSF) {
          scaleFactor *= btagSF;
          if(ievt == 2)cout << "                - applying btag factors " << endl;
        }
        else btagSF = 1.;
        
        if (applyNloSF && isAMC) {
          scaleFactor *= nloWeight * nloSF;
          if(ievt == 2) cout << "                - applying nlo factors " << endl;
        }  // additional SF due to number of events with neg weight!!
        
      }
      else if(isData) scaleFactor = 1.;
      if(ievt == 2) cout << "                ==> scaleFactor " << scaleFactor << endl;
      
      if (makePlots)
      {
        //cout << "ievt " << ievt << endl;
        FillGeneralPlots(d, "control_afterAtLeast1Jet", decayChannels, isData);
        
      }
      //cout << "zmass" << endl;
      if(Zboson.M() < 76 || Zboson.M() > 106) continue;

      if (makePlots)
      {
        FillGeneralPlots(d, "control_afterAtLeast1Jet_afterZWindow", decayChannels,isData);

      }
      
      if(selectednonCSVLJetID.size()>0 && makePlots){
        FillGeneralPlots(d, "control_afterAtLeast1Jet_afterZWindow_afterAtLeast1BJet", decayChannels,isData);
      }
    
     
      
      if(selectedJets.size() == 1 && selectedCSVLJetID.size() > 0){ Region = 0; nSelectedEntriesST++; } // ST region
      else if(selectedJets.size() > 1 && selectedCSVLJetID.size() > 0){ Region = 1; nSelectedEntriesTT++;} // ttbar region
      else if(selectedJets.size() >0 && selectedCSVLJetID.size() == 0){ Region = 2; nSelectedEntriesWZ++;}// WZ control region
      else {continue; }
      
      

      
      
      if(Region == 0 ) nSelectedEntriesSTweighted += scaleFactor*Luminosity/EquilumiSF;
      if(Region == 1 ) nSelectedEntriesTTweighted += scaleFactor*Luminosity/EquilumiSF;
      if(Region == 2 ) nSelectedEntriesWZweighted += scaleFactor*Luminosity/EquilumiSF;
      
      if(makeMVAtree){
        //cout << "ievt " << ievt << endl;
        MakeMVAvars(Region, scaleFactor);
      }
      
      /// Make plots
      if (makePlots )
      {
        if(makeMVAtree && Region == 0) FillMVAPlots(d,dataSetName, Region, "singletop", decayChannels);
        //if(makeMVAtree && Region == 1) FillMVAPlots(d,dataSetName, Region, "toppair", decayChannels);
       // if(makeMVAtree && Region == 2) FillMVAPlots(d,dataSetName, Region, "wzcontrol", decayChannels);
      }
      
    } // events
    
    if(makeMVAtree){
      firstevent = true;
      writeMVAtree();
    }
    if(isData) {
      nSelectedEntriesSTweighted = nSelectedEntriesST;
      nSelectedEntriesTTweighted = nSelectedEntriesTT;
      nSelectedEntriesWZweighted = nSelectedEntriesWZ;
    }
    cout << "                nSelectedEntries ST region: " << nSelectedEntriesST << " weighted " << nSelectedEntriesSTweighted << endl;
    cout << "                nSelectedEntries TT region: " << nSelectedEntriesTT << " weighted " << nSelectedEntriesTTweighted << endl;
    cout << "                nSelectedEntries WZ region: " << nSelectedEntriesWZ  << " weighted " << nSelectedEntriesWZweighted << endl;
    cout << endl;
    
  } // data
  
  
  
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string rootFileName ="NtuplePlots.root";
  string place =pathOutputdate+"MSPlot/";
  vector <string> vlabel_chan = {"uuu", "uue", "eeu", "eee"};
  mkdir(place.c_str(),0777);
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
  cout << "   Output file is " << pathOutputdate+rootFileName << endl;
  
  ///Write histograms
  fout->cd();
  
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    //cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    if(!datafound) temp->setDataLumi(Luminosity);
    if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
    if(name.find("eee")!=std::string::npos) temp->setChannel(true, "eee");
    if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "eeu");
    if(name.find("uue")!=std::string::npos) temp->setChannel(true, "uue");
    if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "uuu");
    if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
    temp->Draw(name, 1, false, false, false, 5);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathOutputdate+"MSPlot", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  
  // 1D
  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  gStyle->SetOptStat(1111);
  for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (pathOutputdate+it->first+".png").c_str() );
  }
  
  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms");
  th2dir->cd();
  gStyle->SetPalette(55);
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
    tempCanvas->SaveAs( (pathOutputdate+it->first+".png").c_str() );
  }
  
  fout->Close();
  
  delete fout;
  
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    float secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  
  return 0;
  
}  // end main



void LeptonAssigner(vector<TLorentzVector> electrons, vector<TLorentzVector> muons, std::vector<int> electronsCharge,std::vector<int> muonsCharge){
  clock_t start_sub =  clock();
  
  //  cout << " in assigner " << endl;
 
  Assigned = false;
  if(electronsCharge.size() != electrons.size() || muons.size() != muonsCharge.size()) cout << "ERROR vectors not filled properly" << endl;
  if(electronsCharge.size() + muonsCharge.size() != 3){
    cout << " WARNING: not 3 leptons " << endl;
    cout << "muons " << muonsCharge.size() << " electrons " << electronsCharge.size() << endl;
  }

  ZelecIndiceF_0 = -999;
  ZelecIndiceF_1 = -999;
  ZmuIndiceF_0 = -999;
  ZmuIndiceF_1 = -999;
  WelecIndiceF = -999;
  WmuIndiceF = -999;
  //cout << " in 3 lep " << endl;
  
  
  bool can01 = false;
  bool can02= false;
  bool can12 = false;
  
  double mass01 = 9999.;
  double mass02 = 9999.;
  double mass12 = 9999.;
  
  if(electronsCharge.size() == 2){
    //cout << "2 electr " << electronsCharge[0] << " " << electronsCharge[1] << endl;
    if(electronsCharge[0] != electronsCharge[1]){
      Assigned = true;
      ZelecIndiceF_0 = 0;
      ZelecIndiceF_1 = 1;
      WmuIndiceF = 0;
    }
    else Assigned = false;
  }
  else if(muonsCharge.size() == 2){
   // cout << "2 muonsCharge" << endl;
    if(muonsCharge[0] != muonsCharge[1]){
      Assigned = true;
      
      ZmuIndiceF_0 = 0;
      ZmuIndiceF_1 = 1;
      WelecIndiceF = 0;
    }
    else Assigned = false;
  }
  else if(electronsCharge.size() ==3){
   // cout << " 3 electronsCharge " << endl;
    can01 = false;
    can02= false;
    can12 = false;
    if(electronsCharge[0] != electronsCharge[1]) can01 = true;
    if(electronsCharge[0] != electronsCharge[2]) can02 = true;
    if(electronsCharge[2] != electronsCharge[1]) can12 = true;
    
    mass01 = 9999.;
    mass02 = 9999.;
    mass12 = 9999.;
    
    
    temp0.SetPxPyPzE(electrons[0].Px(), electrons[0].Py(),electrons[0].Pz(),electrons[0].Energy());
    temp1.SetPxPyPzE(electrons[1].Px(), electrons[1].Py(),electrons[1].Pz(),electrons[1].Energy());
    temp2.SetPxPyPzE(electrons[2].Px(), electrons[2].Py(),electrons[2].Pz(),electrons[2].Energy());
    
    if(can01) mass01 = fabs(91.1-(temp1+temp0).M());
    if(can02) mass02 = fabs(91.1-(temp2+temp0).M());
    if(can12) mass12 = fabs(91.1-(temp1+temp2).M());
    
    if(mass01 <= mass02 && mass01 <= mass12 && can01){
       Assigned = true;
      ZelecIndiceF_0 = 0;ZelecIndiceF_1 = 1;WelecIndiceF = 2;
    }
    
    else if(mass02 <= mass12 && mass02 < mass01 && can02){
        Assigned = true;
      ZelecIndiceF_0 = 0; ZelecIndiceF_1=2;WelecIndiceF = 1;
    }
    else if(mass12 < mass01 && mass12 < mass02 && can12){
      Assigned = true;
      ZelecIndiceF_0 = 1; ZelecIndiceF_1=2;WelecIndiceF = 0;
    }
    else Assigned = false;
  }
  else if(muonsCharge.size() == 3){
    can01 = false;
    can02= false;
    can12 = false;
    
    if(muonsCharge[0] != muonsCharge[1]) can01 = true;
    if(muonsCharge[0] != muonsCharge[2]) can02 = true;
    if(muonsCharge[2] != muonsCharge[1]) can12 = true;
    
    mass01 = 9999.;
    mass02 = 9999.;
    mass12 = 9999.;

    temp0.SetPxPyPzE(muons[0].Px(), muons[0].Py(),muons[0].Pz(),muons[0].Energy());
    temp1.SetPxPyPzE(muons[1].Px(), muons[1].Py(),muons[1].Pz(),muons[1].Energy());
    temp2.SetPxPyPzE(muons[2].Px(), muons[2].Py(),muons[2].Pz(),muons[2].Energy());
    
    if(can01) mass01 = fabs(91.1-(temp1+temp0).M());
    if(can02) mass02 = fabs(91.1-(temp2+temp0).M());
    if(can12) mass12 = fabs(91.1-(temp1+temp2).M());
    
    if(mass01 <= mass02 && mass01 <= mass12 && can01){
       Assigned = true;
      ZmuIndiceF_0 = 0; ZmuIndiceF_1=1; WmuIndiceF = 2;
    }
    else if(mass02 <= mass12 && mass02 < mass01 && can02){
       Assigned = true;
      ZmuIndiceF_0 = 0; ZmuIndiceF_1=2;WmuIndiceF = 1;
    }
    else if(mass12 < mass01 && mass12 < mass02 && can12){
       Assigned = true;
      ZmuIndiceF_0 = 1; ZmuIndiceF_1=2;WmuIndiceF = 0;
    }
    else Assigned = false;
  }

  
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the reconstruction" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }

  
};





void ReconstructObjects(vector<TLorentzVector> Muons,vector<TLorentzVector> selectedElectrons, vector<TLorentzVector> selectedJets,int Region){
  clock_t start_sub =  clock();
  
  
  
  
  //cout << "SMjetIndex "<< SMjetIndex << " selectedJets.size() " << selectedJets.size() << endl;
  SMjetIndex = SMjetCalculator(selectedJets,verbose);
  if(SMjetIndex > selectedJets.size())  cout << "WARNING SMjetIndex "<< SMjetIndex << endl;
  //SMbjet.SetPtEtaPhiE(selectedJets[SMjetIndex].Pt(),selectedJets[SMjetIndex].Eta(), selectedJets[SMjetIndex].Phi(), selectedJets[SMjetIndex].Energy());
  SMbjet.SetPxPyPzE(selectedJets[SMjetIndex].Px(),selectedJets[SMjetIndex].Py(), selectedJets[SMjetIndex].Pz(), selectedJets[SMjetIndex].Energy());
 // cout << "SMjetIndex new"<< SMjetIndex << endl;
  
  
  
  if(WmuIndiceF != -999 && WelecIndiceF != -999) cout << "ERROR: 2 W leptons found" << endl;
  else if(WmuIndiceF!=-999) Wlep.SetPtEtaPhiE(selectedMuons[WmuIndiceF].Pt(),selectedMuons[WmuIndiceF].Eta(), selectedMuons[WmuIndiceF].Phi(), selectedMuons[WmuIndiceF].Energy());
  else if(WelecIndiceF!=-999)Wlep.SetPtEtaPhiE(selectedElectrons[WelecIndiceF].Pt(), selectedElectrons[WelecIndiceF].Eta(), selectedElectrons[WelecIndiceF].Phi(), selectedElectrons[WelecIndiceF].Energy());
  else cout << " ERROR: Wlep not found" << endl;
  
  metTLVbf.SetPtEtaPhiE(met_Pt, met_Eta, met_Phi, TMath::Sqrt(met_Px*met_Px+met_Py*met_Py));
  metTLV = MetzCalculator(Wlep, metTLVbf);
  
  //cout << "met" << endl;
  Wboson = metTLV + Wlep;
  SMtop = metTLV + Wlep  + SMbjet;
  
  //cout << "Wboson mass " << Wboson.M() << " SMtop " << SMtop.M() << " mlb " << (Wlep+SMbjet).M() << endl;
  mWT = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_Px)*(Wlep.Px() + met_Px) - (Wlep.Py() + met_Py)* (Wlep.Py() +met_Py));
  double phis = Wlep.Phi() - met_Phi;
  double cosphis = TMath::Cos(phis);
  mWT2 = TMath::Sqrt(2*Wlep.Pt() * met_Pt*(1-cosphis));


  if(ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999 && ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999 ) cout << "ERROR: 2 Zbosons found " << endl;
  else if(ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999 ) Zboson = selectedMuons[ZmuIndiceF_0] + selectedMuons[ZmuIndiceF_1];
  else if(ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999 ) Zboson = selectedElectrons[ZelecIndiceF_0] + selectedElectrons[ZelecIndiceF_1];
  else cout << "ERROR: Zboson not found " << endl;
  
  if(selectedJets.size()>1 )
  {
    //cout << "cjet" << endl;
    cjetindex = FCNCjetCalculator(selectedJets,Zboson ,SMjetIndex, 3);
    cjetindex_CvsLtagger = FCNCjetCalculatorCvsLTagger(selectedJets,SMjetIndex, 3); //TO FIX
    cjetindex_CvsBtagger = FCNCjetCalculatorCvsBTagger(selectedJets,SMjetIndex, 3); //TO FIX
    cjetindex_Cloose = FCNCjetCalculatorCwp(selectedJets, selectedCharmLJetsindex,  SMjetIndex, 3);
    cjetindex_Cmedium = FCNCjetCalculatorCwp(selectedJets, selectedCharmMJetsindex,  SMjetIndex, 3);
    cjetindex_Ctight = FCNCjetCalculatorCwp(selectedJets, selectedCharmTJetsindex,  SMjetIndex, 3);
    
    if(cjetindex > selectedJets.size() || cjetindex == SMjetIndex) cout << "Warning wrong cjet index " << endl;
    LightJet.SetPxPyPzE(selectedJets[cjetindex].Px(),selectedJets[cjetindex].Py(),selectedJets[cjetindex].Pz(),selectedJets[cjetindex].Energy());
    
    //cout << "fcnctop" << endl;
    FCNCtop = LightJet + Zboson;
    
  }
  //cout << "end reco" << endl;
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the reconstruction" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
};

void MakeMVAvars(int Region, double scaleFactor){
  clock_t start_sub = clock();
  if(true){
    MVA_channel = (float) channelInt;
    MVA_weight = (float) scaleFactor;
    MVA_region = (float) Region;
    
    MVA_lepton0_pt = (float) selectedLeptons[0].Pt();
    MVA_lepton1_pt= (float) selectedLeptons[1].Pt();
    MVA_lepton2_pt= (float) selectedLeptons[2].Pt();
    
    MVA_lepton0_eta= (float) selectedLeptons[0].Eta();
    MVA_lepton1_eta= (float) selectedLeptons[1].Eta();
    MVA_lepton2_eta = (float) selectedLeptons[2].Eta();
    MVA_lepton0_phi = (float) selectedLeptons[0].Phi();
    MVA_lepton1_phi = (float) selectedLeptons[1].Phi();
    MVA_lepton2_phi = (float) selectedLeptons[2].Phi();
    
    MVA_nMuons = (float) selectedMuons.size();
    MVA_nElectrons = (float) selectedElectrons.size();
    // cout << "leptons done" << endl;
    if(selectedJets.size()>0){
      MVA_jet0_pt = (float) selectedJets[0].Pt();
      MVA_jet0_eta = (float) selectedJets[0].Eta();
      MVA_jet0_phi = (float) selectedJets[0].Phi();
    }
    if(selectedJets.size()>1){
      MVA_jet1_pt = (float) selectedJets[1].Pt();
      MVA_jet1_eta = (float) selectedJets[1].Eta();
      MVA_jet1_phi = (float) selectedJets[1].Phi();
    }
    
    MVA_nJets = (float) selectedJets.size();
    MVA_NJets_CSVv2L = (float) selectedCSVLJetID.size();
    MVA_NJets_CSVv2M= (float) selectedCSVMJetID.size();
    MVA_NJets_CSVv2T= (float) selectedCSVTJetID.size();
    //cout << "jets done" << endl;
    
    // SM side
    MVA_Wlep_pt = (float) Wlep.Pt();
    MVA_Wlep_eta = (float) Wlep.Eta();
    MVA_Wlep_phi = (float) Wlep.Phi();
    MVA_SMbjet_pt = (float) selectedJets[SMjetIndex].Pt();
    MVA_SMbjet_eta = (float) selectedJets[SMjetIndex].Eta();
    MVA_SMbjet_phi = (float) selectedJets[SMjetIndex].Phi();
    
    MVA_Wboson_pt = (float) Wboson.Pt();
    MVA_Wboson_eta = (float) Wboson.Eta();
    MVA_Wboson_phi = (float) Wboson.Phi();
    MVA_met = (float) met_Pt;
    MVA_SMtop_pt = (float) SMtop.Pt();
    MVA_SMtop_eta = (float) SMtop.Eta();
    MVA_SMtop_phi = (float) SMtop.Phi();
    
    //cout << "SM side done" << endl;
    
    // FCNC side
    MVA_Zboson_pt = (float) Zboson.Pt();
    MVA_Zboson_eta = (float) Zboson.Eta();
    MVA_Zboson_phi = (float) Zboson.Phi();
    
    if(Region == 1) // ttbar region
    {
      MVA_LightJet_pt = (float) LightJet.Pt();
      MVA_LightJet_eta = (float) LightJet.Eta();
      MVA_LightJet_phi = (float) LightJet.Phi();
      
      MVA_FCNCtop_pt = (float) FCNCtop.Pt();
      MVA_FCNCtop_eta = (float) FCNCtop.Eta();
      MVA_FCNCtop_phi =(float) FCNCtop.Phi();
      
      MVA_nJets_CharmL = (float) selectedCharmLJetsindex.size();
      MVA_nJets_CharmM = (float) selectedCharmMJetsindex.size();
      MVA_nJets_CharmT = (float) selectedCharmTJetsindex.size();
    }
    
    // cout << "FCNC side " << endl;
    
    //SM kinematics
    
    MVA_SMtop_M = SMtop.M();
    MVA_mlb =(SMbjet+Wlep).M();
    MVA_Wboson_M = Wboson.M();
    
    MVA_dRWlepb =Wlep.DeltaR(SMbjet);
    
    MVA_dPhiWlepb =Wlep.DeltaPhi(SMbjet);
    
    if(WmuIndiceF != -999) MVA_Wlep_Charge = selectedMuonsCharge[WmuIndiceF];
    else if(WelecIndiceF != -999) MVA_Wlep_Charge = selectedElectronsCharge[WelecIndiceF];
    MVA_charge_asym = MVA_Wlep_Charge*fabs(Wlep.Eta());
    if(selectedJets.size()>1) MVA_bdiscCSVv2_jet_1 = bdisc_jet[1];
    if(selectedJets.size()>0) MVA_bdiscCSVv2_jet_0 =bdisc_jet[0];
    MVA_CosTheta = (CosThetaCalculation(Wlep, metTLV, SMbjet, false)).first;;
    MVA_CosTheta_alt = (CosThetaCalculation(Wlep, metTLV, SMbjet,false)).second;
    
    //cout << "sm kin" << endl;
    
    // FCNC kinematics
    if(Region == 1)   MVA_FCNCtop_M = FCNCtop.M();
    MVA_Zboson_M = Zboson.M(); // TO FIX
    //cout << "Zboson mass " << Zboson_M << endl;
    
    if(Region == 1) MVA_dRZc = Zboson.DeltaR(LightJet);
    if(Region == 1) MVA_dPhiZc = Zboson.DeltaPhi(LightJet);
    
    if(selectedJets.size()>0) MVA_cdiscCvsB_jet_0 = cdiscCvsB_jet[0];
    if(selectedJets.size()>0) MVA_cdiscCvsL_jet_0 = cdiscCvsL_jet[0];
    if(selectedJets.size()>1) MVA_cdiscCvsB_jet_1 = cdiscCvsB_jet[1];
    if(selectedJets.size()>1)  MVA_cdiscCvsL_jet_1 = cdiscCvsL_jet[1];
    
    //cout << "fcnc kin " << endl;
    
    // interplay
    if(Region == 1) MVA_dRSMFCNCtop = SMtop.DeltaR(FCNCtop);
    MVA_dRZb = Zboson.DeltaR(SMbjet);
    if(Region == 1) MVA_dRWlepc = Wlep.DeltaR(LightJet);
    MVA_dRZWlep = Zboson.DeltaR(Wlep);
    MVA_dRZSMtop = Zboson.DeltaR(SMtop);
    MVA_dPhiZMET = Zboson.DeltaPhi(metTLV);
    
    if(Region == 1) MVA_dPhiSMFCNCtop = SMtop.DeltaPhi(FCNCtop);
    MVA_dPhiZb = Zboson.DeltaPhi(SMbjet);
    if(Region == 1) MVA_dPhiWlepc = Wlep.DeltaPhi(LightJet);
    MVA_dPhiZWlep = Zboson.DeltaPhi(Wlep);
    MVA_dPhiZSMtop = Zboson.DeltaPhi(SMtop);
    MVA_m3l = (selectedLeptons[0]+selectedLeptons[1]+selectedLeptons[2]).M();

    
    
    
  //cout << "interplay done " << endl;
  }
  MVA_mWt = mWT2;
  
  mvatree->Fill();
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the Mva var filler" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
}


void createMVAtree(string dataSetName)
{
  clock_t start_sub = clock();
  tupfile->cd();
  //cout << "does tree exist? " << !tupfile->GetListOfKeys()->Contains("mvatree") <<endl;
  
  
  // event
  
  mvatree->Branch("MVA_channel", &MVA_channel , "MVA_channel/F");
  mvatree->Branch("MVA_weight", &MVA_weight, "MVA_weight/F");
  mvatree->Branch("MVA_region", &MVA_region, "MVA_region/F");
  
  
  // pt
  mvatree->Branch("MVA_lepton0_pt", &MVA_lepton0_pt,"MVA_lepton0_pt/F");
  mvatree->Branch("MVA_lepton1_pt", &MVA_lepton1_pt,"MVA_lepton1_pt/F");
  mvatree->Branch("MVA_lepton2_pt", &MVA_lepton2_pt,"MVA_lepton2_pt/F");
  mvatree->Branch("MVA_Wlep_pt", &MVA_Wlep_pt,"MVA_Wlep_pt/F");
  mvatree->Branch("MVA_Wboson_pt", &MVA_Wboson_pt,"MVA_Wboson_pt/F");
  mvatree->Branch("MVA_SMbjet_pt", &MVA_SMbjet_pt,"MVA_SMbjet_pt/F");
  mvatree->Branch("MVA_SMtop_pt", &MVA_SMtop_pt,"MVA_SMtop_pt/F");
  mvatree->Branch("MVA_Zboson_pt", &MVA_Zboson_pt,"MVA_Zboson_pt/F");
  mvatree->Branch("MVA_LightJet_pt", &MVA_LightJet_pt,"MVA_LightJet_pt/F");
  mvatree->Branch("MVA_FCNCtop_pt", &MVA_FCNCtop_pt,"MVA_FCNCtop_pt/F");
  
  
  // eta
  mvatree->Branch("MVA_lepton0_eta", &MVA_lepton0_eta,"MVA_lepton0_eta/F");
  mvatree->Branch("MVA_lepton1_eta", &MVA_lepton1_eta,"MVA_lepton1_eta/F");
  mvatree->Branch("MVA_lepton2_eta", &MVA_lepton2_eta,"MVA_lepton2_eta/F");
  mvatree->Branch("MVA_Wlep_eta", &MVA_Wlep_eta,"MVA_Wlep_eta/F");
  mvatree->Branch("MVA_Wboson_eta", &MVA_Wboson_eta,"MVA_Wboson_eta/F");
  mvatree->Branch("MVA_SMbjet_eta", &MVA_SMbjet_eta,"MVA_SMbjet_eta/F");
  mvatree->Branch("MVA_SMtop_eta", &MVA_SMtop_eta,"MVA_SMtop_eta/F");
  mvatree->Branch("MVA_Zboson_eta", &MVA_Zboson_eta,"MVA_Zboson_eta/F");
  
  mvatree->Branch("MVA_LightJet_eta", &MVA_LightJet_eta,"MVA_LightJet_eta/F");
  mvatree->Branch("MVA_FCNCtop_eta", &MVA_FCNCtop_eta,"MVA_FCNCtop_eta/F");
  
  // phi
  mvatree->Branch("MVA_lepton0_phi", &MVA_lepton0_phi,"MVA_lepton0_phi/F");
  mvatree->Branch("MVA_lepton1_phi", &MVA_lepton1_phi,"MVA_lepton1_phi/F");
  mvatree->Branch("MVA_lepton2_phi", &MVA_lepton2_phi,"MVA_lepton2_phi/F");
  mvatree->Branch("MVA_Wlep_phi", &MVA_Wlep_phi,"MVA_Wlep_phi/F");
  mvatree->Branch("MVA_Wboson_phi", &MVA_Wboson_phi,"MVA_Wboson_phi/F");
  mvatree->Branch("MVA_SMbjet_phi", &MVA_SMbjet_phi,"MVA_SMbjet_phi/F");
  mvatree->Branch("MVA_SMtop_phi", &MVA_SMtop_phi,"MVA_SMtop_phi/F");
  mvatree->Branch("MVA_Zboson_phi", &MVA_Zboson_phi,"MVA_Zboson_phi/F");
  
  mvatree->Branch("MVA_LightJet_phi", &MVA_LightJet_phi,"MVA_LightJet_phi/F");
  mvatree->Branch("MVA_FCNCtop_phi", &MVA_FCNCtop_phi,"MVA_FCNCtop_phi/F");
  
  
  
  //nbrs
  mvatree->Branch("MVA_nElectrons", &MVA_nElectrons, "MVA_nElectrons/F");
  mvatree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
  mvatree->Branch("MVA_NJets_CSVv2L", &MVA_NJets_CSVv2L, "MVA_NJets_CSVv2L/F");
  mvatree->Branch("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, "MVA_NJets_CSVv2M/F");
  mvatree->Branch("MVA_NJets_CSVv2T", &MVA_NJets_CSVv2T, "MVA_NJets_CSVv2T/F");
  mvatree->Branch("MVA_nMuons", &MVA_nMuons, "MVA_nMuons/F");
  
  
  mvatree->Branch("MVA_met", &MVA_met, "MVA_met/F");
  
  
  //SM kinematics
  mvatree->Branch("MVA_mWt", &MVA_mWt,"MVA_mWt/F");
  mvatree->Branch("MVA_SMtop_M", &MVA_SMtop_M, "MVA_SMtop_M/F");
  mvatree->Branch("MVA_mlb", &MVA_mlb,"MVA_mlb/F");
  mvatree->Branch("MVA_Wboson_M", &MVA_Wboson_M, "MVA_Wboson_M/F");
  
  mvatree->Branch("MVA_dRWlepb", &MVA_dRWlepb,"MVA_dRWlepb/F");
  
  mvatree->Branch("MVA_dPhiWlepb", &MVA_dPhiWlepb, "MVA_dPhiWlepb/F");
  
  mvatree->Branch("MVA_Wlep_Charge", &MVA_Wlep_Charge,"MVA_Wlep_Charge/F");
  mvatree->Branch("MVA_charge_asym", &MVA_charge_asym, "MVA_charge_asym/F");
  mvatree->Branch("MVA_TotalPt", &MVA_TotalPt, "MVA_TotalPt/F");
  mvatree->Branch("MVA_TotalHt", &MVA_TotalHt,"MVA_TotalHt/F");
  mvatree->Branch("MVA_TotalInvMass", &MVA_TotalInvMass,"MVA_TotalInvMass/F");
  
  mvatree->Branch("MVA_TotalPt_jet", &MVA_TotalPt_jet, "MVA_TotalPt_jet/F");
  mvatree->Branch("MVA_TotalHt_jet", &MVA_TotalHt_jet,"MVA_TotalHt_jet/F");
  mvatree->Branch("MVA_TotalInvMass_jet", &MVA_TotalInvMass_jet,"MVA_TotalInvMass_jet/F");
  
  mvatree->Branch("MVA_TotalPt_lep", &MVA_TotalPt_lep, "MVA_TotalPt_lep/F");
  mvatree->Branch("MVA_TotalHt_lep", &MVA_TotalHt_lep,"MVA_TotalHt_lep/F");
  mvatree->Branch("MVA_TotalInvMass_lep", &MVA_TotalInvMass_lep,"MVA_TotalInvMass_lep/F");
  
  mvatree->Branch("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0,"MVA_bdiscCSVv2_jet_0/F");
  mvatree->Branch("MVA_bdiscCSVv2_jet_1", &MVA_bdiscCSVv2_jet_0,"MVA_bdiscCSVv2_jet_1/F");
  mvatree->Branch("MVA_CosTheta", &MVA_CosTheta,"MVA_CosTheta/F");
  mvatree->Branch("MVA_CosTheta_alt", &MVA_CosTheta_alt,"MVA_CosTheta_alt/F");
  
  
  
  // FCNC kinematics
  
  mvatree->Branch("MVA_FCNCtop_M", &MVA_FCNCtop_M, "MVA_FCNCtop_M/F");
  mvatree->Branch("MVA_Zboson_M", &MVA_Zboson_M,"MVA_Zboson_M/F");
  
  
  mvatree->Branch("MVA_dRZc", &MVA_dRZc,"MVA_dRZc/F");
  mvatree->Branch("MVA_dPhiZc", &MVA_dPhiZc,"MVA_dPhiZc/F");
  
  mvatree->Branch("MVA_cdiscCvsB_jet_1", &MVA_cdiscCvsB_jet_1,"MVA_cdiscCvsB_jet_1/F");
  mvatree->Branch("MVA_cdiscCvsL_jet_1", &MVA_cdiscCvsL_jet_1,"MVA_cdiscCvsL_jet_1/F");
  mvatree->Branch("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0,"MVA_cdiscCvsB_jet_0/F");
  mvatree->Branch("MVA_cdiscCvsL_jet_0", &MVA_cdiscCvsL_jet_0,"MVA_cdiscCvsL_jet_0/F");
  
  mvatree->Branch("MVA_nJets_CharmL", &MVA_nJets_CharmL, "MVA_nJets_CharmL/F");
  mvatree->Branch("MVA_nJets_CharmM", &MVA_nJets_CharmM, "MVA_nJets_CharmM/F");
  mvatree->Branch("MVA_nJets_CharmT", &MVA_nJets_CharmT, "MVA_nJets_CharmT/F");
  
  // interplay
  mvatree->Branch("MVA_dRSMFCNCtop", &MVA_dRSMFCNCtop,"MVA_dRSMFCNCtop/F");
  
  mvatree->Branch("MVA_dRZb", &MVA_dRZb,"MVA_dRZb/F");
  mvatree->Branch("MVA_dRWlepc", &MVA_dRWlepc,"MVA_dRWlepc/F");
  mvatree->Branch("MVA_dRZWlep", &MVA_dRZWlep,"MVA_dRZWlep/F");
  mvatree->Branch("MVA_dRZSMtop", &MVA_dRZSMtop,"MVA_dRZSMtop/F");
  
  
  mvatree->Branch("MVA_dPhiSMFCNCtop", &MVA_dPhiSMFCNCtop,"MVA_dPhiSMFCNCtop/F");
  
  mvatree->Branch("MVA_dPhiZb", &MVA_dPhiZb,"MVA_dPhiZb/F");
  mvatree->Branch("MVA_dPhiWlepc", &MVA_dPhiWlepc,"MVA_dPhiWlepc/F");
  mvatree->Branch("MVA_dPhiZWlep", &MVA_dPhiZWlep,"MVA_dPhiZWlep/F");
  mvatree->Branch("MVA_dPhiZSMtop", &MVA_dPhiZSMtop,"MVA_dPhiZSMtop/F");
  
  mvatree->Branch("MVA_m3l", &MVA_m3l, "MVA_m3l/F");
  
  
  // clear values in the vars
  ClearMVAVars();
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the mva tree creator" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
}
void writeMVAtree()
{
  clock_t start_sub = clock();
  tupfile->cd();
  mvatree->Write();
  tupfile->Close();
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run MVA tree writer" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
}

std::pair <Float_t,Float_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo){
  
  clock_t start_sub = clock();
  // see https://github.com/TopBrussels/TopTreeAnalysis/blob/CMSSW_53X/WHelicities/src/BTagCosThetaCalculation.cc
  // ttbar http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2012/157
  // single top : https://cds.cern.ch/record/1601800
  
  
  float CosTheta = -999.;
  
  //----------------------------------------------
  //  Calculating cos theta value
  //----------------------------------------------
  
  TRootMCParticle WLeptonic = (Neutrino+lepton);
  TRootMCParticle TopLeptonic = (Neutrino+lepton+leptonicBJet);
  TLorentzVector TopWRF = (Neutrino+lepton+leptonicBJet);
  TLorentzVector leptWRF = lepton;
  TLorentzVector WTRF;
  WTRF.SetPxPyPzE(WLeptonic.Px(), WLeptonic.Py(), WLeptonic.Pz(), WLeptonic.Energy());
  // TLorentzVector leptTRF = lepton;
  
  //Angle between Top in WRF and lepton in WRF
  TopWRF.Boost(-WLeptonic.BoostVector());
  leptWRF.Boost(-WLeptonic.BoostVector());
  
  // boost to top RF
  WTRF.Boost(-TopLeptonic.BoostVector());
  //leptTRF.Boost(-TopLeptonic.BoostVector());
  
  //Calculating cos:
  float ThetaTevatron = ROOT::Math::VectorUtil::Angle( TopWRF, leptWRF );
  CosTheta = -(TMath::Cos(ThetaTevatron));
  //Cos theta is defined as the angle between the lepton and the reversed direction of the top quark, both boosted to the W-boson rest frame.
  //Still reversed direction doesn't need to be calculated since angle between lepton and top and between lepton and reversed top is proportional to theta and Pi-theta.
  //For these the angles the following relation holds: cos(theta) = - cos(Pi-theta)
  // --> Hence the need of the minus sign in the CosTheta definition!!
  
  float ThetaSTAR = ROOT::Math::VectorUtil::Angle( WTRF, leptWRF );
  float CosThetaSTAR= TMath::Cos(ThetaSTAR);
  
  
  
  if(WLeptonic.E() < 0.){
    cout << " Event with negative WLept energy!!! (BTagCosThetaCalculation class) Cos theta = " << CosTheta << endl;
  }
  
  //cout << "cos " << CosThetaSTAR << " cos " << CosTheta <<  endl;
  /*
   
   As the W-boson helicity fractions are very sensitive to the Wtb vertex, their measurements can be used to investigate contribution from non-standard models.
   In this study, the t ̄t signal events are reconstructed where both W-bosons decay leptonically (electron or muon).
   In this case, the relevant angle θ∗ is defined in top quark rest frame, as the angle between the 3-momentum of the charged lepton in the rest
   frame of the decaying W-boson and the momentum of the W-boson in the rest frame of the decaying top quark.
   */
  
  
  std::pair <Float_t,Float_t> CosThetaPair(CosTheta, CosThetaSTAR);
  
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the Cos Theta Calculator" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
  
  return CosThetaPair;
}



int FCNCjetCalculator(vector < TLorentzVector>  Jets, TLorentzVector recoZ ,int index, int verb){
  clock_t start_sub = clock();
  // cout << "calculating fcnc jet" << endl;
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -5;
  TLorentzVector Jet;
  if(Jets.size() > 0){
    
    for( int iJ = 0; iJ < Jets.size(); iJ++)
    {
      if(iJ == index) continue;
      Jet.SetPxPyPzE(Jets[iJ].Px(),Jets[iJ].Py(),Jets[iJ].Pz(),Jets[iJ].Energy());
      //cout << iJ << " tempMinM " << TempMinMass << " newmass " << (recoZ+Jet).M() ;
      if(fabs((recoZ+Jet).M() - TopMass) < TempMinMass)
      {
        TempMinMass = fabs((recoZ+Jet).M() - TopMass);
        Jetcandidate.SetPxPyPzE(Jet.Px(), Jet.Py(), Jet.Pz(), Jet.E());
        NbInColl = iJ;
        
      }
      //cout << " NbInColl is " << iJ << endl;
    }
  }
  else{
    NbInColl = -999,
    cout << "no cjets available" << endl;
  }
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the FCNC jet reconstruction" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
 // cout << "end finding light jet" <<endl;
  return NbInColl;
};

int FCNCjetCalculatorCvsBTagger(vector < TLorentzVector>  Jets,int index, int verb){
  
  double tempcdis = -9999.;
  int NbInColl = -999;

  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      
      if(cdiscCvsB_jet[iJ]>tempcdis){
        NbInColl = iJ;
        tempcdis = cdiscCvsB_jet[iJ];
      }
      
    }
  }
  else if(Jets.size() == 2){
    if(index == 0) NbInColl = 1;
    if(index == 1) NbInColl = 0;
  }
  else{
    NbInColl = -999,
    cout << "no cjets available" << endl;
  }
  return NbInColl;
};
int FCNCjetCalculatorCvsLTagger(vector < TLorentzVector>  Jets,int index, int verb){
  double tempcdis = -9999.;
  int NbInColl = -999;
  
  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      
      if(cdiscCvsL_jet[iJ]>tempcdis){
        NbInColl = iJ;
        tempcdis = cdiscCvsL_jet[iJ];
      }
      
    }
  }
  else if(Jets.size() == 2){
    if(index == 0) NbInColl = 1;
    if(index == 1) NbInColl = 0;
  }
  else{
    NbInColl = -999,
    cout << "no cjets available" << endl;
  }
  return NbInColl;
};
int FCNCjetCalculatorCwp(vector < TLorentzVector>  Jets, std::vector <int> cjetindex, int index, int verb){
  clock_t start_sub = clock();
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -5;
  bool isCjet = false;
  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      for(int iC = 0; iC < cjetindex.size() ; iC++){
        if(iJ == cjetindex[iJ]) isCjet = true;
      }
      if(!isCjet) continue;
      for(int kJ = 1; kJ < Jets.size(); kJ++){
        if(kJ == index) continue;
        for(int iC = 0; iC < cjetindex.size() ; iC++){
          if(kJ == cjetindex[kJ]) isCjet = true;
        }
        if(!isCjet) continue;
        if(Jets[iJ].Pt()>=Jets[kJ].Pt()) NbInColl = iJ;
        else NbInColl = kJ;
      }
    }
  }
  else if(Jets.size() == 2 && cjetindex.size() > 0){
    if(index == 0) NbInColl = 1;
    if(index == 1) NbInColl = 0;
  }
  else{
    NbInColl = -999;
    // cout << "no cjets available" << endl;
  }
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the FCNC jet reconstruction Cwp" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
  return NbInColl;
};


TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met){
  
  clock_t start_sub = clock();
 // cout << "calculating metz" << endl;
  
  double term1 = leptW.Pz() * ( leptW.Px()* v_met.Px() + leptW.Py()*v_met.Py() + pow(80.399, 2)/2.);
  
  double det = pow(leptW.Px() * v_met.Px() + leptW.Py() * v_met.Py() + pow(80.399, 2)/2., 2) - v_met.Pt()*v_met.Pt() * (leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz() );
  
  if(det<0) det=0;
  
  double term2 = leptW.E() * pow(det, 0.5);
  double denom = leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz();
  double sol1 = (term1 - term2) / denom;
  //double sol2 = (term1 + term2) / denom;
  double nu_E = 0;
  
  TLorentzVector neutrino;
  
  nu_E = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(sol1,2), 0.5);//neglecting neutrino mass
  neutrino.SetPxPyPzE( v_met.Px(), v_met.Py(), sol1, nu_E);
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the metz reconstruction" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
  
  // cout << "end calculating metz" << endl;
  
  return neutrino;
  
  
};



void FillMVAPlots(int d, string dataSetName, int Region, string prefix, vector<int> decayChannels )
{
  clock_t start_sub = clock();
  
  string sregion;
  sregion = prefix+"_";
  string decaystring = "";
  //cout <<  "region " << sregion << endl;
  for(int iChan = 0;iChan < decayChannels.size() ;iChan++){
    if(decayChannels[iChan]== -9) continue;
    if((decayChannels[iChan]) != MVA_channel) continue;
    decaystring= "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
   
    if(Region != 2){
      MSPlot[(sregion+"MVA_channel_"+decaystring).c_str()]->Fill(MVA_channel, datasets[d], true, 1);
      MSPlot[(sregion+"MVA_weight_"+decaystring).c_str()]->Fill(MVA_weight, datasets[d], true, 1);
      
      MSPlot[(sregion+"MVA_lepton0_pt_"+decaystring).c_str()]->Fill(MVA_lepton0_pt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_pt_"+decaystring).c_str()]->Fill(MVA_lepton1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_pt_"+decaystring).c_str()]->Fill(MVA_lepton2_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton0_eta_"+decaystring).c_str()]->Fill(MVA_lepton0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_eta_"+decaystring).c_str()]->Fill(MVA_lepton1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_eta_"+decaystring).c_str()]->Fill(MVA_lepton2_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton0_phi_"+decaystring).c_str()]->Fill(MVA_lepton0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_phi_"+decaystring).c_str()]->Fill(MVA_lepton1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_phi_"+decaystring).c_str()]->Fill(MVA_lepton2_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_jet0_pt_"+decaystring).c_str()]->Fill(MVA_jet0_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_jet0_eta_"+decaystring).c_str()]->Fill(MVA_jet0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_jet0_phi_"+decaystring).c_str()]->Fill(MVA_jet0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_pt_"+decaystring).c_str()]->Fill(MVA_jet1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_eta_"+decaystring).c_str()]->Fill(MVA_jet1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_phi_"+decaystring).c_str()]->Fill(MVA_jet1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      // SM side
      MSPlot[(sregion+"MVA_Wlep_pt_"+decaystring).c_str()]->Fill(MVA_Wlep_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wlep_eta_"+decaystring).c_str()]->Fill(MVA_Wlep_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wlep_phi_"+decaystring).c_str()]->Fill(MVA_Wlep_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_pt_"+decaystring).c_str()]->Fill(MVA_SMbjet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_eta_"+decaystring).c_str()]->Fill(MVA_SMbjet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_phi_"+decaystring).c_str()]->Fill(MVA_SMbjet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_pt_"+decaystring).c_str()]->Fill(MVA_Wboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_eta_"+decaystring).c_str()]->Fill(MVA_Wboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_phi_"+decaystring).c_str()]->Fill(MVA_Wboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_met_"+decaystring).c_str()]->Fill(MVA_met, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_pt_"+decaystring).c_str()]->Fill(MVA_SMtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_eta_"+decaystring).c_str()]->Fill(MVA_SMtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_phi_"+decaystring).c_str()]->Fill(MVA_SMtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // FCNC side
      MSPlot[(sregion+"MVA_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_eta_"+decaystring).c_str()]->Fill(MVA_Zboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_phi_"+decaystring).c_str()]->Fill(MVA_Zboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // nbrs
      MSPlot[(sregion+"MVA_nMuons_"+decaystring).c_str()]->Fill(MVA_nMuons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2T_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2T, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2M_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2L_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2L, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_"+decaystring).c_str()]->Fill(MVA_nJets, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nElectrons_"+decaystring).c_str()]->Fill(MVA_nElectrons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      
      //SM kinematics
      MSPlot[(sregion+"MVA_SMtop_M_"+decaystring).c_str()]->Fill(MVA_SMtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_mlb_"+decaystring).c_str()]->Fill(MVA_mlb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_M_"+decaystring).c_str()]->Fill(MVA_Wboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dRWlepb_"+decaystring).c_str()]->Fill(MVA_dRWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dPhiWlepb_"+decaystring).c_str()]->Fill(MVA_dPhiWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_Wlep_Charge_"+decaystring).c_str()]->Fill(MVA_Wlep_Charge , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalPt_"+decaystring).c_str()]->Fill(MVA_TotalPt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalHt_"+decaystring).c_str()]->Fill(MVA_TotalHt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalInvMass_"+decaystring).c_str()]->Fill( MVA_TotalInvMass, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_TotalPt_jet_"+decaystring).c_str()]->Fill(MVA_TotalPt_jet , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalHt_jet_"+decaystring).c_str()]->Fill(MVA_TotalHt_jet , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalInvMass_jet_"+decaystring).c_str()]->Fill( MVA_TotalInvMass_jet, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_TotalPt_lep_"+decaystring).c_str()]->Fill(MVA_TotalPt_lep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalHt_lep_"+decaystring).c_str()]->Fill(MVA_TotalHt_lep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalInvMass_lep_"+decaystring).c_str()]->Fill( MVA_TotalInvMass_lep, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_CosTheta_"+decaystring).c_str()]->Fill(MVA_CosTheta , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_CosTheta_alt_"+decaystring).c_str()]->Fill(MVA_CosTheta_alt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // FCNC kinematics
      
      MSPlot[(sregion+"MVA_dRZb_"+decaystring).c_str()]->Fill(MVA_dRZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRZWlep_"+decaystring).c_str()]->Fill(MVA_dRZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRZSMtop_"+decaystring).c_str()]->Fill(MVA_dRZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      MSPlot[(sregion+"MVA_dPhiZb_"+decaystring).c_str()]->Fill(MVA_dPhiZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZWlep_"+decaystring).c_str()]->Fill(MVA_dPhiZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZMET_"+decaystring).c_str()]->Fill(MVA_dPhiZMET , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZSMtop_"+decaystring).c_str()]->Fill(MVA_dPhiZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_m3l_"+decaystring).c_str()]->Fill(MVA_m3l , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    }
    
    MSPlot[(sregion+"MVA_mWt_"+decaystring).c_str()]->Fill(MVA_mWt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    if(Region == 1){
      MSPlot[(sregion+"MVA_LightJet_pt_"+decaystring).c_str()]->Fill(MVA_LightJet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_LightJet_eta_"+decaystring).c_str()]->Fill(MVA_LightJet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_LightJet_phi_"+decaystring).c_str()]->Fill(MVA_LightJet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_FCNCtop_pt_"+decaystring).c_str()]->Fill(MVA_FCNCtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_FCNCtop_eta_"+decaystring).c_str()]->Fill(MVA_FCNCtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_FCNCtop_phi_"+decaystring).c_str()]->Fill(MVA_FCNCtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_nJets_CharmL_"+decaystring).c_str()]->Fill(MVA_nJets_CharmL, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_CharmM_"+decaystring).c_str()]->Fill(MVA_nJets_CharmM, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_CharmT_"+decaystring).c_str()]->Fill(MVA_nJets_CharmT, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_FCNCtop_M_"+decaystring).c_str()]->Fill(MVA_FCNCtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_M_"+decaystring).c_str()]->Fill(MVA_Zboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dRZc_"+decaystring).c_str()]->Fill(MVA_dRZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZc_"+decaystring).c_str()]->Fill(MVA_dPhiZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_cdiscCvsB_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_cdiscCvsL_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      // interplay
      MSPlot[(sregion+"MVA_dRSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dRSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dPhiSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_bdiscCSVv2_jet_1_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    }
    
  }
  
  for(int iChan = 0;iChan < decayChannels.size() ;iChan++){
    if(decayChannels[iChan]!= -9) continue;
    if((decayChannels[iChan]) != MVA_channel) continue;
    decaystring= "all";

    
    if(Region != 2){
      MSPlot[(sregion+"MVA_channel_"+decaystring).c_str()]->Fill(MVA_channel, datasets[d], true, 1);
      MSPlot[(sregion+"MVA_weight_"+decaystring).c_str()]->Fill(MVA_weight, datasets[d], true, 1);
      
      MSPlot[(sregion+"MVA_lepton0_pt_"+decaystring).c_str()]->Fill(MVA_lepton0_pt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_pt_"+decaystring).c_str()]->Fill(MVA_lepton1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_pt_"+decaystring).c_str()]->Fill(MVA_lepton2_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton0_eta_"+decaystring).c_str()]->Fill(MVA_lepton0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_eta_"+decaystring).c_str()]->Fill(MVA_lepton1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_eta_"+decaystring).c_str()]->Fill(MVA_lepton2_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton0_phi_"+decaystring).c_str()]->Fill(MVA_lepton0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton1_phi_"+decaystring).c_str()]->Fill(MVA_lepton1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_lepton2_phi_"+decaystring).c_str()]->Fill(MVA_lepton2_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_jet0_pt_"+decaystring).c_str()]->Fill(MVA_jet0_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_jet0_eta_"+decaystring).c_str()]->Fill(MVA_jet0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_jet0_phi_"+decaystring).c_str()]->Fill(MVA_jet0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_pt_"+decaystring).c_str()]->Fill(MVA_jet1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_eta_"+decaystring).c_str()]->Fill(MVA_jet1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      if(selectedJets.size()>1) MSPlot[(sregion+"MVA_jet1_phi_"+decaystring).c_str()]->Fill(MVA_jet1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      // SM side
      MSPlot[(sregion+"MVA_Wlep_pt_"+decaystring).c_str()]->Fill(MVA_Wlep_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wlep_eta_"+decaystring).c_str()]->Fill(MVA_Wlep_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wlep_phi_"+decaystring).c_str()]->Fill(MVA_Wlep_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_pt_"+decaystring).c_str()]->Fill(MVA_SMbjet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_eta_"+decaystring).c_str()]->Fill(MVA_SMbjet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMbjet_phi_"+decaystring).c_str()]->Fill(MVA_SMbjet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_pt_"+decaystring).c_str()]->Fill(MVA_Wboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_eta_"+decaystring).c_str()]->Fill(MVA_Wboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_phi_"+decaystring).c_str()]->Fill(MVA_Wboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_met_"+decaystring).c_str()]->Fill(MVA_met, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_pt_"+decaystring).c_str()]->Fill(MVA_SMtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_eta_"+decaystring).c_str()]->Fill(MVA_SMtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_SMtop_phi_"+decaystring).c_str()]->Fill(MVA_SMtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // FCNC side
      MSPlot[(sregion+"MVA_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_eta_"+decaystring).c_str()]->Fill(MVA_Zboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_phi_"+decaystring).c_str()]->Fill(MVA_Zboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // nbrs
      MSPlot[(sregion+"MVA_nMuons_"+decaystring).c_str()]->Fill(MVA_nMuons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2T_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2T, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2M_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_NJets_CSVv2L_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2L, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_"+decaystring).c_str()]->Fill(MVA_nJets, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nElectrons_"+decaystring).c_str()]->Fill(MVA_nElectrons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      
      //SM kinematics
      MSPlot[(sregion+"MVA_SMtop_M_"+decaystring).c_str()]->Fill(MVA_SMtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_mlb_"+decaystring).c_str()]->Fill(MVA_mlb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Wboson_M_"+decaystring).c_str()]->Fill(MVA_Wboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dRWlepb_"+decaystring).c_str()]->Fill(MVA_dRWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dPhiWlepb_"+decaystring).c_str()]->Fill(MVA_dPhiWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_Wlep_Charge_"+decaystring).c_str()]->Fill(MVA_Wlep_Charge , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalPt_"+decaystring).c_str()]->Fill(MVA_TotalPt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalHt_"+decaystring).c_str()]->Fill(MVA_TotalHt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_TotalInvMass_"+decaystring).c_str()]->Fill( MVA_TotalInvMass, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_CosTheta_"+decaystring).c_str()]->Fill(MVA_CosTheta , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_CosTheta_alt_"+decaystring).c_str()]->Fill(MVA_CosTheta_alt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      // FCNC kinematics
      
      MSPlot[(sregion+"MVA_dRZb_"+decaystring).c_str()]->Fill(MVA_dRZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRZWlep_"+decaystring).c_str()]->Fill(MVA_dRZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRZSMtop_"+decaystring).c_str()]->Fill(MVA_dRZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      MSPlot[(sregion+"MVA_dPhiZb_"+decaystring).c_str()]->Fill(MVA_dPhiZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZWlep_"+decaystring).c_str()]->Fill(MVA_dPhiZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZMET_"+decaystring).c_str()]->Fill(MVA_dPhiZMET , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZSMtop_"+decaystring).c_str()]->Fill(MVA_dPhiZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_m3l_"+decaystring).c_str()]->Fill(MVA_m3l , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    }
    
    MSPlot[(sregion+"MVA_mWt_"+decaystring).c_str()]->Fill(MVA_mWt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    if(Region == 1){
      MSPlot[(sregion+"MVA_LightJet_pt_"+decaystring).c_str()]->Fill(MVA_LightJet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_LightJet_eta_"+decaystring).c_str()]->Fill(MVA_LightJet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_LightJet_phi_"+decaystring).c_str()]->Fill(MVA_LightJet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_FCNCtop_pt_"+decaystring).c_str()]->Fill(MVA_FCNCtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_FCNCtop_eta_"+decaystring).c_str()]->Fill(MVA_FCNCtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_FCNCtop_phi_"+decaystring).c_str()]->Fill(MVA_FCNCtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_nJets_CharmL_"+decaystring).c_str()]->Fill(MVA_nJets_CharmL, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_CharmM_"+decaystring).c_str()]->Fill(MVA_nJets_CharmM, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_nJets_CharmT_"+decaystring).c_str()]->Fill(MVA_nJets_CharmT, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_FCNCtop_M_"+decaystring).c_str()]->Fill(MVA_FCNCtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_Zboson_M_"+decaystring).c_str()]->Fill(MVA_Zboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_dRZc_"+decaystring).c_str()]->Fill(MVA_dRZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiZc_"+decaystring).c_str()]->Fill(MVA_dPhiZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_cdiscCvsB_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_cdiscCvsL_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      // interplay
      MSPlot[(sregion+"MVA_dRSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dRSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dRWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dPhiSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"MVA_dPhiWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"MVA_bdiscCSVv2_jet_1_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    }
    
  }
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the fill mva plots" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
}





/////////////////SM B
int SMjetCalculator(vector<TLorentzVector> Jets,int verb){
  clock_t start_sub = clock();
  // cout << "calculating bjet" << endl;
  int index_ = -999 ;
  double tempbdis = -999.;
  
  for(int iJ = 0; iJ < Jets.size() ; iJ++){
    
    if(bdisc_jet[iJ] > tempbdis)
    {
      index_ = iJ;
      tempbdis = bdisc_jet[iJ];
      //cout << "index is " << iJ << endl ;
    }
  }
  
  
  
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the SM bjet reconstruction" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  // cout << "end calculating bjet" << endl;
  //cout << "SMjetIndex "<< SMjetIndex << endl;
  return index_;
};




void GetMetaData(TTree* tree, bool isData,int Entries, bool isAMC)
{
  clock_t start_sub = clock();
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("nEv", &nEv, &b_nEv);
  
  tree->SetBranchAddress("WPb_L", &WPb_L, &b_WPb_L);
  tree->SetBranchAddress("WPb_M", &WPb_M, &b_WPb_M);
  tree->SetBranchAddress("WPb_T", &WPb_T, &b_WPb_T);
  tree->SetBranchAddress("WPc_CvsB_Loose",&WPc_CvsB_Loose,&b_WPc_CvsB_Loose);
   tree->SetBranchAddress("WPc_CvsB_Medium",&WPc_CvsB_Medium,&b_WPc_CvsB_Medium);
   tree->SetBranchAddress("WPc_CvsB_Tight",&WPc_CvsB_Tight,&b_WPc_CvsB_Tight);
   tree->SetBranchAddress("WPc_CvsL_Loose",&WPc_CvsL_Loose,&b_WPc_CvsL_Loose);
   tree->SetBranchAddress("WPc_CvsL_Medium",&WPc_CvsL_Medium,&b_WPc_CvsL_Medium);
   tree->SetBranchAddress("WPc_CvsL_Tight",&WPc_CvsL_Tight,&b_WPc_CvsL_Tight);
  tree->SetBranchAddress("nofPosWeights",&nofPosWeights, &b_nofPosWeights);
  tree->SetBranchAddress("nofNegWeights",&nofNegWeights, &b_nofNegWeights);
  tree->SetBranchAddress("sumW",&sumW, &b_sumW);
  
  /* TO FIX
  tree->SetBranchAddress("nbEv_Initial_all",&nbEv_Initial_all,&b_nbEv_Initial_all);
  tree->SetBranchAddress("nbEv_Trigged_all",&nbEv_Trigged_all,&b_nbEv_Trigged_all);
  tree->SetBranchAddress("nbEv_3lep_all",&nbEv_3lep_all,&b_nbEv_3lep_all);
  tree->SetBranchAddress("nbEv_LooseLepVeto_all",&nbEv_LooseLepVeto_all,&b_nbEv_LooseLepVeto_all);
  tree->SetBranchAddress("nbEv_ZbosonWindow_all",&nbEv_ZbosonWindow_all,&b_nbEv_ZbosonWindow_all);
  tree->SetBranchAddress("nbEv_AtLeast1jet_all",&nbEv_AtLeast1jet_all,&b_nbEv_AtLeast1jet_all);
  
  tree->SetBranchAddress("nbEv_Initial_uuu",&nbEv_Initial_uuu,&b_nbEv_Initial_uuu);
  tree->SetBranchAddress("nbEv_Trigged_uuu",&nbEv_Trigged_uuu,&b_nbEv_Trigged_uuu);
  tree->SetBranchAddress("nbEv_3lep_uuu",&nbEv_3lep_uuu,&b_nbEv_3lep_uuu);
  tree->SetBranchAddress("nbEv_LooseLepVeto_uuu",&nbEv_LooseLepVeto_uuu,&b_nbEv_LooseLepVeto_uuu);
  tree->SetBranchAddress("nbEv_ZbosonWindow_uuu",&nbEv_ZbosonWindow_uuu,&b_nbEv_ZbosonWindow_uuu);
  tree->SetBranchAddress("nbEv_AtLeast1jet_uuu",&nbEv_AtLeast1jet_uuu,&b_nbEv_AtLeast1jet_uuu);
  
  tree->SetBranchAddress("nbEv_Initial_uue",&nbEv_Initial_uue,&b_nbEv_Initial_uue);
  tree->SetBranchAddress("nbEv_Trigged_uue",&nbEv_Trigged_uue,&b_nbEv_Trigged_uue);
  tree->SetBranchAddress("nbEv_3lep_uue",&nbEv_3lep_uue,&b_nbEv_3lep_uue);
  tree->SetBranchAddress("nbEv_LooseLepVeto_uue",&nbEv_LooseLepVeto_uue,&b_nbEv_LooseLepVeto_uue);
  tree->SetBranchAddress("nbEv_ZbosonWindow_uue",&nbEv_ZbosonWindow_uue,&b_nbEv_ZbosonWindow_uue);
  tree->SetBranchAddress("nbEv_AtLeast1jet_uue",&nbEv_AtLeast1jet_uue,&b_nbEv_AtLeast1jet_uue);
  
  tree->SetBranchAddress("nbEv_Initial_eeu",&nbEv_Initial_eeu,&b_nbEv_Initial_eeu);
  tree->SetBranchAddress("nbEv_Trigged_eeu",&nbEv_Trigged_eeu,&b_nbEv_Trigged_eeu);
  tree->SetBranchAddress("nbEv_3lep_eeu",&nbEv_3lep_eeu,&b_nbEv_3lep_eeu);
  tree->SetBranchAddress("nbEv_LooseLepVeto_eeu",&nbEv_LooseLepVeto_eeu,&b_nbEv_LooseLepVeto_eeu);
  tree->SetBranchAddress("nbEv_ZbosonWindow_eeu",&nbEv_ZbosonWindow_eeu,&b_nbEv_ZbosonWindow_eeu);
  tree->SetBranchAddress("nbEv_AtLeast1jet_eeu",&nbEv_AtLeast1jet_eeu,&b_nbEv_AtLeast1jet_eeu);
  
  tree->SetBranchAddress("nbEv_Initial_eee",&nbEv_Initial_eee,&b_nbEv_Initial_eee);
  tree->SetBranchAddress("nbEv_Trigged_eee",&nbEv_Trigged_eee,&b_nbEv_Trigged_eee);
  tree->SetBranchAddress("nbEv_3lep_eee",&nbEv_3lep_eee,&b_nbEv_3lep_eee);
  tree->SetBranchAddress("nbEv_LooseLepVeto_eee",&nbEv_LooseLepVeto_eee,&b_nbEv_LooseLepVeto_eee);
  tree->SetBranchAddress("nbEv_ZbosonWindow_eee",&nbEv_ZbosonWindow_eee,&b_nbEv_ZbosonWindow_eee);
  tree->SetBranchAddress("nbEv_AtLeast1jet_eee",&nbEv_AtLeast1jet_eee,&b_nbEv_AtLeast1jet_eee);
  
  */
  
  
  
  
  TotalEvents = 0;
  
  int nPos = 0;
  int nNeg = 0;
  int Weights = 0;
  
  
  for (int k = 0; k< Entries; k++)
  {
    tree->GetEntry(k);
    TotalEvents += nEv;
    /* TO FIX
    total_nbEv_3lep_all += nbEv_3lep_all;
    total_nbEv_AtLeast1jet_all += nbEv_AtLeast1jet_all;
    total_nbEv_Initial_all += nbEv_Initial_all;
    total_nbEv_LooseLepVeto_all += nbEv_LooseLepVeto_all;
    total_nbEv_Trigged_all += nbEv_Trigged_all;
    total_nbEv_ZbosonWindow_all += nbEv_ZbosonWindow_all;
    
    total_nbEv_3lep_uuu += nbEv_3lep_uuu;
    total_nbEv_AtLeast1jet_uuu += nbEv_AtLeast1jet_uuu;
    total_nbEv_Initial_uuu += nbEv_Initial_uuu;
    total_nbEv_LooseLepVeto_uuu += nbEv_LooseLepVeto_uuu;
    total_nbEv_Trigged_uuu += nbEv_Trigged_uuu;
    total_nbEv_ZbosonWindow_uuu += nbEv_ZbosonWindow_uuu;
    
    total_nbEv_3lep_uue += nbEv_3lep_uue;
    total_nbEv_AtLeast1jet_uue += nbEv_AtLeast1jet_uue;
    total_nbEv_Initial_uue += nbEv_Initial_uue;
    total_nbEv_LooseLepVeto_uue += nbEv_LooseLepVeto_uue;
    total_nbEv_Trigged_uue += nbEv_Trigged_uue;
    total_nbEv_ZbosonWindow_uue += nbEv_ZbosonWindow_uue;
    
    total_nbEv_3lep_eeu += nbEv_3lep_eeu;
    total_nbEv_AtLeast1jet_eeu += nbEv_AtLeast1jet_eeu;
    total_nbEv_Initial_eeu += nbEv_Initial_eeu;
    total_nbEv_LooseLepVeto_eeu += nbEv_LooseLepVeto_eeu;
    total_nbEv_Trigged_eeu += nbEv_Trigged_eeu;
    total_nbEv_ZbosonWindow_eeu += nbEv_ZbosonWindow_eeu;
    
    total_nbEv_3lep_eee += nbEv_3lep_eee;
    total_nbEv_AtLeast1jet_eee += nbEv_AtLeast1jet_eee;
    total_nbEv_Initial_eee += nbEv_Initial_eee;
    total_nbEv_LooseLepVeto_eee += nbEv_LooseLepVeto_eee;
    total_nbEv_Trigged_eee += nbEv_Trigged_eee;
    total_nbEv_ZbosonWindow_eee += nbEv_ZbosonWindow_eee;
    */
    if( isAMC && !isData) nPos += nofPosWeights;
    if(isAMC && !isData) nNeg += nofNegWeights;
    if(isAMC && !isData) Weights += sumW;
    
  }
  
  if(!isData){
    EquilumiSF = TotalEvents / Xsect;
    cout << "                equilumi = " <<  TotalEvents <<" / " << Xsect <<" = " << EquilumiSF << endl;
  }
  else EquilumiSF = 1.;
  if(isAMC && !isData) nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
  if(isAMC && !isData) cout << "                 nloSF: " << nloSF << endl;
  
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the GetMetaData" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
  
  
}

string ConvertIntToString(int Number, int pad)
{
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}


string MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}
void ClearObjects()
{
  //  ClearLeaves();
  ClearTLVs();
  ClearVars();
  
}
void ClearMVAVars()
{
  MVA_channel = -999.;
  MVA_weight = 1.;
  MVA_region = -999.;
  
  MVA_lepton0_pt = -999.;
  MVA_lepton1_pt= -999.;
  MVA_lepton2_pt= -999.;
  MVA_lepton0_eta= -999.;
  MVA_lepton1_eta= -999.;
  MVA_lepton2_eta= -999.;
  MVA_lepton0_phi= -999.;
  MVA_lepton1_phi= -999.;
  MVA_lepton2_phi= -999.;
  
  MVA_jet0_pt = -999.;
  MVA_jet0_eta= -999.;
  MVA_jet0_phi= -999.;
  MVA_jet1_pt = -999.;
  MVA_jet1_eta = -999.;
  MVA_jet1_phi = -999.;
  
  
  
  // SM side
  MVA_Wlep_pt= -999.;
  MVA_Wlep_eta= -999.;
  MVA_Wlep_phi= -999.;
  MVA_SMbjet_pt= -999.;
  MVA_SMbjet_eta= -999.;
  MVA_SMbjet_phi= -999.;
  MVA_Wboson_pt= -999.;
  MVA_Wboson_eta= -999.;
  MVA_Wboson_phi= -999.;
  MVA_met= -999.;
  MVA_SMtop_pt= -999.;
  MVA_SMtop_eta= -999.;
  MVA_SMtop_phi= -999.;
  
  
  // FCNC side
  MVA_Zboson_pt= -999.;
  MVA_Zboson_eta= -999.;
  MVA_Zboson_phi= -999.;
  
  
  MVA_LightJet_pt= -999.;
  MVA_LightJet_eta= -999.;
  MVA_LightJet_phi= -999.;
  
  MVA_FCNCtop_pt= -999.;
  MVA_FCNCtop_eta= -999.;
  MVA_FCNCtop_phi= -999.;
  
  
  //nbrs
  MVA_nMuons = -999.;
  MVA_NJets_CSVv2T = -999.;
  MVA_NJets_CSVv2M = -999.;
  MVA_NJets_CSVv2L = -999.;
  MVA_nJets = -999.;
  MVA_nJets_CharmL = -999.;
  MVA_nJets_CharmM = -999.;
  MVA_nJets_CharmT = -999.;
  MVA_nElectrons = -999.;
  
  //SM kinematics
  MVA_mWt = -999.;
  MVA_SMtop_M = -999.;
  MVA_mlb = -999.;
  MVA_Wboson_M = -999.;
  
  MVA_dRWlepb = -999.;
  
  MVA_dPhiWlepb = -999.;
  
  MVA_Wlep_Charge = -999.;
  MVA_charge_asym = -999.;
  MVA_TotalPt = -999.;
  MVA_TotalHt = -999.;
  MVA_TotalInvMass = -999.;
  MVA_TotalPt_jet = -999.;
  MVA_TotalHt_jet = -999.;
  MVA_TotalInvMass_jet = -999.;
  MVA_TotalPt_lep = -999.;
  MVA_TotalHt_lep = -999.;
  MVA_TotalInvMass_lep = -999.;
  MVA_bdiscCSVv2_jet_0 = -999.;
  MVA_bdiscCSVv2_jet_1 = -999.;
  MVA_CosTheta = -999.;
  MVA_CosTheta_alt = -999.;
  
  
  
  // FCNC kinematics
  MVA_FCNCtop_M = -999.;
  MVA_Zboson_M = -999.;
  
  MVA_dRZc = -999.;
  MVA_dPhiZc = -999.;
  
  MVA_cdiscCvsB_jet_1 = -999.;
  MVA_cdiscCvsL_jet_1 = -999.;
  MVA_cdiscCvsB_jet_0 = -999.;
  MVA_cdiscCvsL_jet_0 = -999.;
  
  // interplay
  MVA_dRSMFCNCtop = -999.;
  MVA_dRZb = -999.;
  MVA_dRWlepc = -999.;
  MVA_dRZWlep = -999.;
  MVA_dRZSMtop = -999.;
  
  MVA_dPhiSMFCNCtop = -999.;
  MVA_dPhiZb = -999.;
  MVA_dPhiWlepc = -999.;
  MVA_dPhiZWlep = -999.;
  MVA_dPhiZMET = -999.;
  MVA_dPhiZSMtop = -999.;
  
  
  
  MVA_m3l = -999.;
  
}
void ClearTLVs(){
  selectedMuons.clear();
  selectedElectrons.clear();
  selectedLeptons.clear();
  selectedJets.clear();
  muon.Clear();
  electron.Clear();
  jet.Clear();
  metTLV.Clear();
  metTLVbf.Clear();
  Wboson.Clear();
  Zboson.Clear();
  SMtop.Clear();
  FCNCtop.Clear();
  Wlep.Clear();
  LightJet.Clear();
  tempInvMassObj.Clear();
  tempInvMassObj_jet.Clear();
  temp0.Clear();
  temp1.Clear();
  temp2.Clear();

  
  
}
void ClearVars(){
  Assigned = false;
  scaleFactor = 1.;
  muonSFtemp = 1.;
  mWT = -999.;
  mWT2 = -999.;
  electronSFtemp = 1.;
  selectedElectronsCharge.clear();
  selectedMuonsCharge.clear();
  selectednonCvsLTJetID.clear();
  selectednonCvsLMJetID.clear();
  selectednonCvsLLJetID.clear();
  selectednonCvsBTJetID.clear();
  selectednonCvsBMJetID.clear();
  selectednonCvsBLJetID.clear();
  selectednonCSVTJetID.clear();
  selectednonCSVMJetID.clear();
  selectednonCSVLJetID.clear();
  selectedCvsLTJetID.clear();
  selectedCvsLMJetID.clear();
  selectedCvsLLJetID.clear();
  selectedCvsBTJetID.clear();
  selectedCvsBMJetID.clear();
  selectedCvsBLJetID.clear();
  selectedCSVTJetID.clear();
  selectedCSVMJetID.clear();
  selectedCSVLJetID.clear();
  selectedCharmLJetsindex.clear();
  selectedCharmMJetsindex.clear();
  selectedCharmTJetsindex.clear();
  
  SMjetIndex = -999;
  cjetindex = -999;
  cjetindex_CvsLtagger = -999;
  cjetindex_CvsBtagger = -999;
  cjetindex_Cloose = -999;
  cjetindex_Cmedium = -999;
  cjetindex_Ctight = -999;
  WmuIndiceF = -999;
  WelecIndiceF = -999;
  ZelecIndiceF_0= -999;
  ZelecIndiceF_1 = -999;
  ZmuIndiceF_1 = -999;
  ZmuIndiceF_0 = -999;
  Region = -9;
  
  tempPx = -999.;
  tempPy = -999.;
  tempHt = -999.;
  
  //cout << "WmuIndiceF " << WmuIndiceF <<" WelecIndiceF "<< WelecIndiceF <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
  
}
void ClearMetaData()
{
  nEv = 0;
  nloSF = 1.;
  nofPosWeights = 0;
  nofNegWeights = 0;
  sumW = 0;
  TotalEvents = 0;
  EquilumiSF = 1.;
  Xsect = 1.;
  nEntries = 0;
  globalnEntries = 0;
  WPb_L= 0.;
  WPb_M= 0.;
  WPb_T= 0.;
  WPc_CvsB_Loose= 0.;
  WPc_CvsB_Medium= 0.;
  WPc_CvsB_Tight= 0.;
  WPc_CvsL_Loose= 0.;
  WPc_CvsL_Medium= 0.;
  Double_t WPc_CvsL_Tight= 0.;
  
  // TO FIX
  total_nbEv_3lep_all= 0;
  total_nbEv_AtLeast1jet_all= 0;
  total_nbEv_Initial_all= 0;
  total_nbEv_LooseLepVeto_all= 0;
  total_nbEv_Trigged_all= 0;
  total_nbEv_ZbosonWindow_all= 0;
  
  total_nbEv_3lep_uuu= 0;
  total_nbEv_AtLeast1jet_uuu= 0;
  total_nbEv_Initial_uuu= 0;
  total_nbEv_LooseLepVeto_uuu= 0;
  total_nbEv_Trigged_uuu= 0;
  total_nbEv_ZbosonWindow_uuu= 0;
  total_nbEv_3lep_uue= 0;
  total_nbEv_AtLeast1jet_uue= 0;
  total_nbEv_Initial_uue= 0;
  total_nbEv_LooseLepVeto_uue= 0;
  total_nbEv_Trigged_uue= 0;
  total_nbEv_ZbosonWindow_uue= 0;
  total_nbEv_3lep_eeu= 0;
  total_nbEv_AtLeast1jet_eeu= 0;
  total_nbEv_Initial_eeu= 0;
  total_nbEv_LooseLepVeto_eeu= 0;
  total_nbEv_Trigged_eeu= 0;
  total_nbEv_ZbosonWindow_eeu= 0;
  total_nbEv_3lep_eee= 0;
  total_nbEv_AtLeast1jet_eee= 0;
  total_nbEv_Initial_eee= 0;
  total_nbEv_LooseLepVeto_eee= 0;
  total_nbEv_Trigged_eee= 0;
  total_nbEv_ZbosonWindow_eee= 0;
}


void Init1DPlots()
{
  TH1::SetDefaultSumw2();
  
}





void InitMSPlots(string prefix, vector <int> decayChannels)
{
  clock_t start_sub = clock();
      string decaystring = "";
  // control plots
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //MSPlot[plotname.c_str()]->setChannel(true, decayChan);
    
    
    MSPlot[(prefix+"_NbOfVertices_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_NbOfVertices_"+decaystring).c_str(), 60, 0, 60, "Nb Of vertices");
    //MSPlot[(prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_NbOfVertices_bfPU_"+decaystring).c_str(), 60, 0, 60, "Nb Of vertices before PU reweighing"); TO FIX
    MSPlot[(prefix+"_puSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_puSF_"+decaystring).c_str(), 200, 0, 2, "Pile up SF");
    
    //cout << "init " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_bdisc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_bdisc_"+decaystring).c_str(), 23, 0., 1, "CSVv2 discriminant");
    MSPlot[(prefix+"_bdisc_bfBT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_bdisc_bfBT_"+decaystring).c_str(),30, 0., 1, "CSVv2 discriminant before Btag SF");
    MSPlot[(prefix+"_btagSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_btagSF_"+decaystring).c_str(), 80, 0.9, 1.3, "Btag SF");
    
    
    MSPlot[(prefix+"_nMu_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nMu_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of Muons");
    MSPlot[(prefix+"_nMu_bfMuSF_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nMu_bfMuSF_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of Muons before Muon SF");
    MSPlot[(prefix+"_muSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_muSF_"+decaystring).c_str(), 60, 0.75, 1.05, "Muon SF");
    
    MSPlot[(prefix+"_nEl_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nEl_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of Electrons");
    MSPlot[(prefix+"_nEl_bfElSF_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nEl_bfElSF_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of Electrons before Electron SF");
    MSPlot[(prefix+"_elSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_elSF_"+decaystring).c_str(), 60, 0.8, 1.05, "Electron SF");
    
     MSPlot[(prefix+"_JetPt_bfJER_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_JetPt_bfJER_"+decaystring).c_str(), 100, 0, 1, "Jet Pt before JER");
     MSPlot[(prefix+"_JetPt_bfJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_JetPt_bfJES_"+decaystring).c_str(), 100, 0, 1, "Jet Pt before JES, after JER");
     MSPlot[(prefix+"_JetPt_afJER_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_JetPt_afJER_"+decaystring).c_str(), 100, 0, 1, "Jet Pt after JER");
     MSPlot[(prefix+"_JetPt_afJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_JetPt_afJES_"+decaystring).c_str(), 100, 0, 1, "Jet Pt after JES, after JER");
     MSPlot[(prefix+"_met_bfJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_met_bfJES_"+decaystring).c_str(), 100, 0, 1, "MET before JES");
     MSPlot[(prefix+"_met_afJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_met_afJES_"+decaystring).c_str(), 100, 0, 1, "MET after JES");
    
    
    
    // vars
    MSPlot[(prefix+"_ZbosonMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_ZbosonMass_"+decaystring).c_str(), 70, 60, 130, "Inv Mass Zboson");
    MSPlot[(prefix+"_WbosonMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_WbosonMass_"+decaystring).c_str(), 200, 0, 200, "Inv Mass Wboson");
    MSPlot[(prefix+"_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_mlb_"+decaystring).c_str(), 100, 0, 500, "Inv Mass (l_{W},SMbjet)");
    MSPlot[(prefix+"_SMTopMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_SMTopMass_"+decaystring).c_str(), 100, 0, 500, "Inv Mass SMTop");
    MSPlot[(prefix+"_mWT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_mWT_"+decaystring).c_str(), 100, 0, 400, "Transv. Mass Wboson");
    MSPlot[(prefix+"_mWT2_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_mWT2_"+decaystring).c_str(), 50, 0, 300, "Transv. Mass Wboson");
    
    MSPlot[(prefix+"_LeadingJetPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_LeadingJetPt_"+decaystring).c_str(), 125, 0, 500, "Leading Jet Pt [GeV]");
    MSPlot[(prefix+"_LeadingLepPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_LeadingLepPt_"+decaystring).c_str(), 125, 0, 500, "Leading Lepton Pt [GeV]");
    MSPlot[(prefix+"_2ndLeadingJetPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_2ndLeadingJetPt_"+decaystring).c_str(), 50, 30, 250, "2ndLeading Jet Pt [GeV]");
    MSPlot[(prefix+"_2ndLeadingLepPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_2ndLeadingLepPt_"+decaystring).c_str(), 40, 20, 250, "2nd Leading Lepton Pt [GeV]");
    
    MSPlot[(prefix+"_nJets_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJets_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of Jets");
    MSPlot[(prefix+"_nJetsCSVL_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJetsCSVL_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of CSVL");
    MSPlot[(prefix+"_nJetsCSVM_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJetsCSVM_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of CSVM");
    MSPlot[(prefix+"_nJetsCSVT_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJetsCSVT_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of CSVT");
    MSPlot[(prefix+"_nJetsCvsL_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJetsCvsL_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of CvsL Jets");
    MSPlot[(prefix+"_nJetsCvsB_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefix+"_nJetsCvsB_"+decaystring).c_str(), 10, -0.5, 9.5, "Nb of CvsB Jets");
    
  }
  
  
  
   MSPlot[(prefix+"_Decay").c_str()]= new MultiSamplePlot(datasets, (prefix+"_Decay").c_str(), 10, -0.5, 9.5, "Nb. Events with decay");
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the InitMSPlots" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
}

void InitMVAMSPlotsWZ(string prefix, vector <int> decayChannels){
  
  clock_t start_sub = clock();
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";

    MSPlot[(prefix+"_MVA_mWt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mWt_"+decaystring).c_str(),100, 0, 100, "M_{T}(W)");
    
  }
  
}

void InitMVAMSPlotsSingletop(string prefix, vector <int> decayChannels)
{
  clock_t start_sub = clock();
  
  InitMVAMSPlotsWZ(prefix, decayChannels);
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    MSPlot[ (prefix+"_MVA_channel_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_channel_"+decaystring).c_str(), 4,-0.5, 4.5, "decay");
    MSPlot[ (prefix+"_MVA_weight_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_weight_"+decaystring).c_str(), 100,-0.5, 9.5, "eventweight");
    
    MSPlot[ (prefix+"_MVA_lepton0_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(lep0) [GeV]");
    MSPlot[ (prefix+"_MVA_lepton1_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(lep1) [GeV]");
    MSPlot[ (prefix+"_MVA_lepton2_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(lep2) [GeV]");
    MSPlot[ (prefix+"_MVA_lepton0_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_eta_"+decaystring).c_str(),60,-6, 6, "#eta (lep0)");
    MSPlot[ (prefix+"_MVA_lepton1_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_eta_"+decaystring).c_str(),60,-6, 6, "#eta (lep1)");
    MSPlot[ (prefix+"_MVA_lepton2_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_eta_"+decaystring).c_str(),60,-6, 6, "#eta (lep2)");
    MSPlot[ (prefix+"_MVA_lepton0_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_phi_"+decaystring).c_str(),40,-4, 4, "#phi (lep0)");
    MSPlot[ (prefix+"_MVA_lepton1_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_phi_"+decaystring).c_str(),40,-4, 4, "#phi (lep1)");
    MSPlot[ (prefix+"_MVA_lepton2_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_phi_"+decaystring).c_str(),40,-4, 4, "#phi (lep2)");
    
    MSPlot[ (prefix+"_MVA_jet0_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (jet0) [GeV]");
    MSPlot[ (prefix+"_MVA_jet0_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_eta_"+decaystring).c_str(),60,-6, 6, "#eta (jet0)");
    MSPlot[ (prefix+"_MVA_jet0_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_phi_"+decaystring).c_str(),40,-4, 4, "#phi (jet0)");

    // SM side
    MSPlot[ (prefix+"_MVA_Wlep_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(l_{W})[GeV]");
    MSPlot[ (prefix+"_MVA_Wlep_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_eta_"+decaystring).c_str(),60,-6, 6, "#eta (l_{W})");
    MSPlot[ (prefix+"_MVA_Wlep_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_phi_"+decaystring).c_str(),40,-4, 4, "#phi (l_{W})");
    MSPlot[ (prefix+"_MVA_SMbjet_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (Bjet) [GeV]");
    MSPlot[ (prefix+"_MVA_SMbjet_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_eta_"+decaystring).c_str(),60,-6, 6, "#eta (Bjet)");
    MSPlot[ (prefix+"_MVA_SMbjet_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_phi_"+decaystring).c_str(),40,-4, 4, "#phi (Bjet)");
    MSPlot[ (prefix+"_MVA_Wboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(W) [GeV]");
    MSPlot[ (prefix+"_MVA_Wboson_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_eta_"+decaystring).c_str(),60,-6, 6, "#eta (W)");
    MSPlot[ (prefix+"_MVA_Wboson_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_phi_"+decaystring).c_str(),40,-4, 4, "#phi (W)");
    MSPlot[ (prefix+"_MVA_met_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_met_"+decaystring).c_str(), 500,0, 500, "met");
    MSPlot[ (prefix+"_MVA_SMtop_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_pt_"+decaystring).c_str(), 500,0, 500, "p_{T}(SMtop) [GeV]");
    MSPlot[ (prefix+"_MVA_SMtop_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_eta_"+decaystring).c_str(),60,-6, 6, "#eta (SMtop)");
    MSPlot[ (prefix+"_MVA_SMtop_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_phi_"+decaystring).c_str(),40,-4, 4, "#phi (SMtop)");
    
    // FCNC side
    MSPlot[ (prefix+"_MVA_Zboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (Z)[GeV]");
    MSPlot[ (prefix+"_MVA_Zboson_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_eta_"+decaystring).c_str(),60,-6, 6, "#eta (Z)");
    MSPlot[ (prefix+"_MVA_Zboson_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_phi_"+decaystring).c_str(),40,-4, 4, "#phi (Z)");
    
    // nbrs
    MSPlot[ (prefix+"_MVA_nMuons_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nMuons_"+decaystring).c_str(), 10,-0.5, 9.5, "nb Muons");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2T_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2T_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CSVv2T");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2M_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2M_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CSVv2M");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2L_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2L_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CSVv2L");
    MSPlot[ (prefix+"_MVA_nJets_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_"+decaystring).c_str(), 10,-0.5, 9.5, "nb Jets");
    MSPlot[ (prefix+"_MVA_nElectrons_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nElectrons_"+decaystring).c_str(), 10,-0.5, 9.5, "nb Electrons");
    
    //SM kinematics
    MSPlot[(prefix+"_MVA_SMtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_M_"+decaystring).c_str(),300, 0,300, "M(SMtop)");
    MSPlot[(prefix+"_MVA_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mlb_"+decaystring).c_str(),100, 0, 100, "M(l_{W}b)");
    MSPlot[(prefix+"_MVA_Wboson_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_M_"+decaystring).c_str(),100, 0, 100, "M(W)");
    
    MSPlot[(prefix+"_MVA_dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRWlepb_"+decaystring).c_str(),100,-10, 10, "dR(l_{W}b)");
    
    MSPlot[(prefix+"_MVA_dPhiWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiWlepb_"+decaystring).c_str(),40,-4, 4, "d#phi(l_{W}b)");
    
    MSPlot[(prefix+"_MVA_Wlep_Charge_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_Charge_"+decaystring).c_str(),4, -2, 2, "Q(l_{W})");
    MSPlot[(prefix+"_MVA_charge_asym_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_charge_asym_"+decaystring).c_str(),100, -10, 10, "Q(l_{W})|#eta(W)|");
    MSPlot[(prefix+"_MVA_TotalPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_"+decaystring).c_str(), 500, 0, 5000, "total P_{T}");
    MSPlot[(prefix+"_MVA_TotalHt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_"+decaystring).c_str(),500, 0, 5000, "total H_{T}");
    MSPlot[(prefix+"_MVA_TotalInvMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_"+decaystring).c_str(),500, 0, 5000, "total Inv Mass");
    MSPlot[(prefix+"_MVA_TotalPt_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_jet_"+decaystring).c_str(), 500, 0, 5000, "total jet and met P_{T}");
    MSPlot[(prefix+"_MVA_TotalHt_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_jet_"+decaystring).c_str(),500, 0, 5000, "total jet and met H_{T}");
    MSPlot[(prefix+"_MVA_TotalInvMass_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_jet_"+decaystring).c_str(),500, 0, 5000, "total jet and met Inv Mass");
    MSPlot[(prefix+"_MVA_TotalPt_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_lep_"+decaystring).c_str(), 500, 0, 5000, "total lepton P_{T}");
    MSPlot[(prefix+"_MVA_TotalHt_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_lep_"+decaystring).c_str(),500, 0, 5000, "total lepton H_{T}");
    MSPlot[(prefix+"_MVA_TotalInvMass_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_lep_"+decaystring).c_str(),500, 0, 5000, "total lepton Inv Mass");
    
    
    MSPlot[(prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str(),100, 0, 1, "CSVv2 Highest pt jet");
    MSPlot[(prefix+"_MVA_CosTheta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_CosTheta_"+decaystring).c_str(),25, -1, 1, "Cos(#theta *)");
    MSPlot[(prefix+"_MVA_CosTheta_alt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_CosTheta_alt_"+decaystring).c_str(),25, -1, 1, "Cos(#theta *)");
    
    // FCNC kinematics
    MSPlot[(prefix+"_MVA_Zboson_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_M_"+decaystring).c_str(),300, 0,300, "M(Z)");
    
    // interplay
    MSPlot[(prefix+"_MVA_dRZb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZb_"+decaystring).c_str(),100,-10, 10, "dR(Z,Bjet)");
    MSPlot[(prefix+"_MVA_dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZWlep_"+decaystring).c_str(),100,-10, 10, "dR(Z,l_{W})");
    MSPlot[(prefix+"_MVA_dRZSMtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZSMtop_"+decaystring).c_str(),100,-10, 10, "dR(Z,SMtop)");
    
    MSPlot[(prefix+"_MVA_dPhiZb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZb_"+decaystring).c_str(),40,-4, 4, "d#phi (Z,Bjet)");
    MSPlot[(prefix+"_MVA_dPhiZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZWlep_"+decaystring).c_str(),40,-4, 4, "d#phi (Z, l_{W})");
    MSPlot[(prefix+"_MVA_dPhiZMET_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZMET_"+decaystring).c_str(),40,-4, 4, "d#phi (Z,met)");
    MSPlot[(prefix+"_MVA_dPhiZSMtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZSMtop_"+decaystring).c_str(),40,-4, 4, "d#phi (Z,SMtop)");
    
    MSPlot[(prefix+"_MVA_m3l_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_m3l_"+decaystring).c_str(),200,0, 200, "Inv Mass Leptons");
  }
  
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the IntitMVAMSplots ST" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
}

void InitMVAMSPlotsTopPair(string prefix, vector <int> decayChannels)
{
  clock_t start_sub = clock();
  
  InitMVAMSPlotsSingletop(prefix, decayChannels);
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    MSPlot[(prefix+"_MVA_bdiscCSVv2_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_bdiscCSVv2_jet_1_"+decaystring).c_str(),100, 0, 1, "CSVv2 2nd Highest pt jet");
    MSPlot[(prefix+"_MVA_cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsB_jet_0_"+decaystring).c_str(),100, 0, 1, "CvsB Highest pt jet");
    MSPlot[(prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str(),100, 0, 1, "CvsL Highest pt jet");
    
    
    MSPlot[(prefix+"_MVA_jet1_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (jet1) [GeV]");
    MSPlot[(prefix+"_MVA_jet1_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_eta_"+decaystring).c_str(),60,-6, 6, "#eta (jet1)");
    MSPlot[(prefix+"_MVA_jet1_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_phi_"+decaystring).c_str(),40,-4, 4, "#phi (jet1)");
    
    
    MSPlot[(prefix+"_MVA_LightJet_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (Ljet) [GeV]");
    MSPlot[(prefix+"_MVA_LightJet_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_eta_"+decaystring).c_str(),60,-6, 6, "#eta (Ljet)");
    MSPlot[(prefix+"_MVA_LightJet_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_phi_"+decaystring).c_str(),40,-4, 4, "#phi (Ljet)");
    
    MSPlot[(prefix+"_MVA_FCNCtop_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_pt_"+decaystring).c_str(), 500,0, 500, "p_{T} (FCNCtop) [GeV]");
    MSPlot[(prefix+"_MVA_FCNCtop_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_eta_"+decaystring).c_str(),60,-6, 6, "#eta (FCNCtop)");
    MSPlot[(prefix+"_MVA_FCNCtop_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_phi_"+decaystring).c_str(),40,-4, 4, "#phi (FCNCtop)");
    
    MSPlot[(prefix+"_MVA_nJets_CharmL_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmL_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CharmL jets");
    MSPlot[(prefix+"_MVA_nJets_CharmM_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmM_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CharmM jets");
    MSPlot[(prefix+"_MVA_nJets_CharmT_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmT_"+decaystring).c_str(), 10,-0.5, 9.5, "nb CharmT jets");
    
    
    // FCNC kinematics
    MSPlot[(prefix+"_MVA_FCNCtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_M_"+decaystring).c_str(),300, 0,300, "M(FCNCtop)");
    
    MSPlot[(prefix+"_MVA_dRZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZc_"+decaystring).c_str(),100,-10, 10, "dR(Z,Ljet)");
    MSPlot[(prefix+"_MVA_dPhiZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZc_"+decaystring).c_str(),40,-4, 4, "d#phi (Z,Ljet)");
    
    MSPlot[(prefix+"_MVA_cdiscCvsB_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsB_jet_1_"+decaystring).c_str(),100, 0, 1, "CvsB 2nd Highest pt jet");
    MSPlot[(prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str(),100, 0, 1, "CvsL 2nd Highest pt jet");
    
    // interplay
    MSPlot[(prefix+"_MVA_dRSMFCNCtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRSMFCNCtop_"+decaystring).c_str(),100,-10, 10, "dR(SMtop,FCNCtop)");
    MSPlot[(prefix+"_MVA_dRWlepc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRWlepc_"+decaystring).c_str(),100,-10, 10, "dR(l_{W},Ljet)");
    
    MSPlot[(prefix+"_MVA_dPhiSMFCNCtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiSMFCNCtop_"+decaystring).c_str(),40,-4, 4, "d#phi (SMtop,FCNCtop)");
    MSPlot[(prefix+"_MVA_dPhiWlepc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiWlepc_"+decaystring).c_str(),40,-4, 4, "d#phi (l_{W},c)");
  }
  double time_sub = ((double)clock() - start_sub) / CLOCKS_PER_SEC;
  if(firstevent && verbose > 3){
    cout << "It took us " << time_sub << " s to run the IntitMVAMSplots ST" << endl;
    if ( time_sub >= 60 )
    {
      int mins = time_sub/60;
      float secs = time_sub - mins*60;
      
      if (mins >= 60 )
      {
        int hours = mins/60;
        mins = mins - hours*60;
        cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
      }
      else
        cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
    }
  }
  
  
}

void Init2DPlots()
{
  TH2::SetDefaultSumw2();
  //histo2D["CosTheta"]= new TH2F("CosTheta", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);
}

void InitTree(TTree* tree, bool isData)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("channelInt", &channelInt, &b_channelInt);
  tree->SetBranchAddress("nloWeight", &nloWeight, &b_nloWeight);
  tree->SetBranchAddress("run_num", &run_num, &b_run_num);
  tree->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
  tree->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
  tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree->SetBranchAddress("npu", &npu, &b_npu);
  tree->SetBranchAddress("puSF", &puSF, &b_puSF);
  tree->SetBranchAddress("btagSF", &btagSF, &b_btagSF);
  tree->SetBranchAddress("PassedMETFilter", &PassedMETFilter, &b_PassedMETFilter);
  tree->SetBranchAddress("PassedGoodPV", &PassedGoodPV, &b_PassedGoodPV);
  tree->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
  tree->SetBranchAddress("ElectronSF", ElectronSF, &b_ElectronSF);
  tree->SetBranchAddress("pt_electron", pt_electron, &b_pt_electron);
  tree->SetBranchAddress("phi_electron", phi_electron, &b_phi_electron);
  tree->SetBranchAddress("eta_electron", eta_electron, &b_eta_electron);
  tree->SetBranchAddress("eta_superCluster_electron", eta_superCluster_electron, &b_eta_superCluster_electron);
  tree->SetBranchAddress("E_electron", E_electron, &b_E_electron);
  tree->SetBranchAddress("chargedHadronIso_electron", chargedHadronIso_electron, &b_chargedHadronIso_electron);
  tree->SetBranchAddress("neutralHadronIso_electron", neutralHadronIso_electron, &b_neutralHadronIso_electron);
  tree->SetBranchAddress("photonIso_electron", photonIso_electron, &b_photonIso_electron);
  tree->SetBranchAddress("pfIso_electron", pfIso_electron, &b_pfIso_electron);
  tree->SetBranchAddress("charge_electron", charge_electron, &b_charge_electron);
  tree->SetBranchAddress("d0_electron", d0_electron, &b_d0_electron);
  tree->SetBranchAddress("d0BeamSpot_electron", d0BeamSpot_electron, &b_d0BeamSpot_electron);
  tree->SetBranchAddress("sigmaIEtaIEta_electron", sigmaIEtaIEta_electron, &b_sigmaIEtaIEta_electron);
  tree->SetBranchAddress("deltaEtaIn_electron", deltaEtaIn_electron, &b_deltaEtaIn_electron);
  tree->SetBranchAddress("deltaPhiIn_electron", deltaPhiIn_electron, &b_deltaPhiIn_electron);
  tree->SetBranchAddress("hadronicOverEm_electron", hadronicOverEm_electron, &b_hadronicOverEm_electron);
  tree->SetBranchAddress("missingHits_electron", missingHits_electron, &b_missingHits_electron);
  tree->SetBranchAddress("passConversion_electron", passConversion_electron, &b_passConversion_electron);
  tree->SetBranchAddress("isId_electron", isId_electron, &b_isId_electron);
  tree->SetBranchAddress("isIso_electron", isIso_electron, &b_isIso_electron);
  tree->SetBranchAddress("isEBEEGap", isEBEEGap, &b_isEBEEGap);
  tree->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
  tree->SetBranchAddress("MuonIDSF", MuonIDSF, &b_MuonIDSF);
  tree->SetBranchAddress("MuonIsoSF", MuonIsoSF, &b_MuonIsoSF);
  tree->SetBranchAddress("MuonTrigSFv2", MuonTrigSFv2, &b_MuonTrigSFv2);
  tree->SetBranchAddress("MuonTrigSFv3", MuonTrigSFv3, &b_MuonTrigSFv3);
  tree->SetBranchAddress("pt_muon", pt_muon, &b_pt_muon);
  tree->SetBranchAddress("phi_muon", phi_muon, &b_phi_muon);
  tree->SetBranchAddress("eta_muon", eta_muon, &b_eta_muon);
  tree->SetBranchAddress("E_muon", E_muon, &b_E_muon);
  tree->SetBranchAddress("chargedHadronIso_muon", chargedHadronIso_muon, &b_chargedHadronIso_muon);
  tree->SetBranchAddress("neutralHadronIso_muon", neutralHadronIso_muon, &b_neutralHadronIso_muon);
  tree->SetBranchAddress("photonIso_muon", photonIso_muon, &b_photonIso_muon);
  tree->SetBranchAddress("isId_muon", isId_muon, &b_isId_muon);
  tree->SetBranchAddress("isIso_muon", isIso_muon, &b_isIso_muon);
  tree->SetBranchAddress("pfIso_muon", pfIso_muon, &b_pfIso_muon);
  tree->SetBranchAddress("charge_muon", charge_muon, &b_charge_muon);
  tree->SetBranchAddress("d0_muon", d0_muon, &b_d0_muon);
  tree->SetBranchAddress("d0BeamSpot_muon", d0BeamSpot_muon, &b_d0BeamSpot_muon);
  tree->SetBranchAddress("nJets", &nJets, &b_nJets);
  tree->SetBranchAddress("pt_jet", pt_jet, &b_pt_jet);
  tree->SetBranchAddress("px_jet", px_jet, &b_px_jet);
  tree->SetBranchAddress("py_jet", py_jet, &b_py_jet);
  tree->SetBranchAddress("pz_jet", pz_jet, &b_pz_jet);
  tree->SetBranchAddress("phi_jet", phi_jet, &b_phi_jet);
  tree->SetBranchAddress("eta_jet", eta_jet, &b_eta_jet);
  tree->SetBranchAddress("E_jet", E_jet, &b_E_jet);
  tree->SetBranchAddress("charge_jet", charge_jet, &b_charge_jet);
  tree->SetBranchAddress("bdisc_jet", bdisc_jet, &b_bdisc_jet);
  tree->SetBranchAddress("jet_Pt_before_JER", jet_Pt_before_JER, &b_jet_Pt_before_JER);
  tree->SetBranchAddress("jet_Pt_before_JES", jet_Pt_before_JES, &b_jet_Pt_before_JES);
  tree->SetBranchAddress("jet_Pt_after_JER", jet_Pt_after_JER, &b_jet_Pt_after_JER);
  tree->SetBranchAddress("jet_Pt_after_JES", jet_Pt_after_JES, &b_jet_Pt_after_JES);
  tree->SetBranchAddress("cdiscCvsL_jet", cdiscCvsL_jet, &b_cdiscCvsL_jet);
  tree->SetBranchAddress("cdiscCvsB_jet", cdiscCvsB_jet, &b_cdiscCvsB_jet);
  tree->SetBranchAddress("met_Pt", &met_Pt, &b_met_Pt);
  tree->SetBranchAddress("met_Eta", &met_Eta, &b_met_Eta);
  tree->SetBranchAddress("met_Phi", &met_Phi, &b_met_Phi);
  tree->SetBranchAddress("met_Px", &met_Px, &b_met_Px);
  tree->SetBranchAddress("met_Py", &met_Py, &b_met_Py);
  tree->SetBranchAddress("met_before_JES", &met_before_JES, &b_met_before_JES);
  tree->SetBranchAddress("met_after_JES", &met_after_JES, &b_met_after_JES);
}


void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData)
{
  //cout << "fill plots" << endl;
      string decaystring = "";
  double eventW = 1.;
  if(!isData) eventW = Luminosity*scaleFactor/EquilumiSF;
  else eventW = Luminosity*scaleFactor/EquilumiSF;
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
   // cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    if(decayChannels[iChan] == -9) continue;;
   // cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    if(decayChannels[iChan] != channelInt) continue;
   // cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    
    
    //cout << "nvtx "  << nvtx << endl;
    MSPlot[(prefix+"_NbOfVertices_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true,eventW);
    // MSPlot[(prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true, Luminosity*(scaleFactor/puSF)/(EquilumiSF)); TO FIX
    MSPlot[(prefix+"_puSF_"+decaystring).c_str()]->Fill(puSF , datasets[d], true, 1);
    
    //cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_bdisc_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true,eventW);
    MSPlot[(prefix+"_bdisc_bfBT_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true, Luminosity*(scaleFactor/btagSF)/(EquilumiSF));
    MSPlot[(prefix+"_btagSF_"+decaystring).c_str()]->Fill(btagSF , datasets[d], true, 1);
    
    MSPlot[(prefix+"_nMu_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nMu_bfMuSF_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true, Luminosity*scaleFactor/(EquilumiSF*muonSFtemp));
    MSPlot[(prefix+"_muSF_"+decaystring).c_str()]->Fill(muonSFtemp , datasets[d], true, 1);
    
    MSPlot[(prefix+"_nEl_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nEl_bfElSF_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true, Luminosity*scaleFactor/(EquilumiSF*electronSFtemp));
    MSPlot[(prefix+"_elSF_"+decaystring).c_str()]->Fill(electronSFtemp , datasets[d], true, 1);
    
    if(selectedJets.size() >0){
      // cout << (prefix+"_JetPt_bfJER_"+decaystring).c_str() << endl;
      MSPlot[(prefix+"_JetPt_bfJER_"+decaystring).c_str()]->Fill(jet_Pt_before_JER[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_bfJES_"+decaystring).c_str()]->Fill(jet_Pt_before_JES[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_afJER_"+decaystring).c_str()]->Fill(jet_Pt_after_JER[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_afJES_"+decaystring).c_str()]->Fill(jet_Pt_after_JES[0] , datasets[d], true,eventW);
    }
    MSPlot[(prefix+"_met_bfJES_"+decaystring).c_str()]->Fill(met_before_JES , datasets[d], true,eventW);
    MSPlot[(prefix+"_met_afJES_"+decaystring).c_str()]->Fill(met_after_JES , datasets[d], true,eventW);
    
    
    // vars
    
    MSPlot[(prefix+"_ZbosonMass_"+decaystring).c_str()]->Fill((float) Zboson.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_WbosonMass_"+decaystring).c_str()] ->Fill((float) Wboson.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_mWT_"+decaystring).c_str()]->Fill(mWT, datasets[d], true,eventW);
    MSPlot[(prefix+"_mWT2_"+decaystring).c_str()]->Fill(mWT2, datasets[d], true,eventW);
    MSPlot[(prefix+"_SMTopMass_"+decaystring).c_str()]->Fill((float) SMtop.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_mlb_"+decaystring).c_str()]->Fill((float) (SMbjet+Wlep).M() , datasets[d], true,eventW);
    
    
    if(selectedJets.size()>0) MSPlot[(prefix+"_LeadingJetPt_"+decaystring).c_str()]->Fill((float) selectedJets[0].Pt() , datasets[d], true,eventW);
    MSPlot[(prefix+"_LeadingLepPt_"+decaystring).c_str()]->Fill((float) selectedLeptons[0].Pt() , datasets[d], true,eventW);
    if(selectedJets.size()>1) MSPlot[(prefix+"_2ndLeadingJetPt_"+decaystring).c_str()]->Fill((float) selectedJets[1].Pt() , datasets[d], true,eventW);
    MSPlot[(prefix+"_2ndLeadingLepPt_"+decaystring).c_str()]->Fill((float) selectedLeptons[1].Pt() , datasets[d], true,eventW);
    
    MSPlot[(prefix+"_nJets_"+decaystring).c_str()]->Fill(selectedJets.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVL_"+decaystring).c_str()]->Fill(selectedCSVLJetID.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVM_"+decaystring).c_str()]->Fill(selectedCSVMJetID.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVT_"+decaystring).c_str()]->Fill(selectedCSVTJetID.size() , datasets[d], true,eventW);
  }
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
   // cout << decayChannels[iChan] << endl;
    if(decayChannels[iChan] != -9) continue;
    //cout << "fill all" << endl;
    decaystring = "all";
  
    
    //cout << "nvtx "  << nvtx << endl;
    MSPlot[(prefix+"_NbOfVertices_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true,eventW);
    // MSPlot[(prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true, Luminosity*(scaleFactor/puSF)/(EquilumiSF)); TO FIX
    MSPlot[(prefix+"_puSF_"+decaystring).c_str()]->Fill(puSF , datasets[d], true, 1);
    
    //cout << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_bdisc_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true,eventW);
    MSPlot[(prefix+"_bdisc_bfBT_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true, Luminosity*(scaleFactor/btagSF)/(EquilumiSF));
    MSPlot[(prefix+"_btagSF_"+decaystring).c_str()]->Fill(btagSF , datasets[d], true, 1);
    
    MSPlot[(prefix+"_nMu_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nMu_bfMuSF_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true, Luminosity*scaleFactor/(EquilumiSF*muonSFtemp));
    MSPlot[(prefix+"_muSF_"+decaystring).c_str()]->Fill(muonSFtemp , datasets[d], true, 1);
    
    MSPlot[(prefix+"_nEl_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nEl_bfElSF_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true, Luminosity*scaleFactor/(EquilumiSF*electronSFtemp));
    MSPlot[(prefix+"_elSF_"+decaystring).c_str()]->Fill(electronSFtemp , datasets[d], true, 1);
    
    if(selectedJets.size() >0){
      // cout << (prefix+"_JetPt_bfJER_"+decaystring).c_str() << endl;
      MSPlot[(prefix+"_JetPt_bfJER_"+decaystring).c_str()]->Fill(jet_Pt_before_JER[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_bfJES_"+decaystring).c_str()]->Fill(jet_Pt_before_JES[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_afJER_"+decaystring).c_str()]->Fill(jet_Pt_after_JER[0] , datasets[d], true,eventW);
      MSPlot[(prefix+"_JetPt_afJES_"+decaystring).c_str()]->Fill(jet_Pt_after_JES[0] , datasets[d], true,eventW);
    }
    MSPlot[(prefix+"_met_bfJES_"+decaystring).c_str()]->Fill(met_before_JES , datasets[d], true,eventW);
    MSPlot[(prefix+"_met_afJES_"+decaystring).c_str()]->Fill(met_after_JES , datasets[d], true,eventW);
    
    
    // vars
    
    MSPlot[(prefix+"_ZbosonMass_"+decaystring).c_str()]->Fill((float) Zboson.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_WbosonMass_"+decaystring).c_str()] ->Fill((float) Wboson.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_mWT_"+decaystring).c_str()]->Fill(mWT, datasets[d], true,eventW);
    MSPlot[(prefix+"_mWT2_"+decaystring).c_str()]->Fill(mWT2, datasets[d], true,eventW);
    MSPlot[(prefix+"_SMTopMass_"+decaystring).c_str()]->Fill((float) SMtop.M() , datasets[d], true,eventW);
    MSPlot[(prefix+"_mlb_"+decaystring).c_str()]->Fill((float) (SMbjet+Wlep).M() , datasets[d], true,eventW);
    
    
    if(selectedJets.size()>0) MSPlot[(prefix+"_LeadingJetPt_"+decaystring).c_str()]->Fill((float) selectedJets[0].Pt() , datasets[d], true,eventW);
    MSPlot[(prefix+"_LeadingLepPt_"+decaystring).c_str()]->Fill((float) selectedLeptons[0].Pt() , datasets[d], true,eventW);
    if(selectedJets.size()>1) MSPlot[(prefix+"_2ndLeadingJetPt_"+decaystring).c_str()]->Fill((float) selectedJets[1].Pt() , datasets[d], true,eventW);
    MSPlot[(prefix+"_2ndLeadingLepPt_"+decaystring).c_str()]->Fill((float) selectedLeptons[1].Pt() , datasets[d], true,eventW);
    
    MSPlot[(prefix+"_nJets_"+decaystring).c_str()]->Fill(selectedJets.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVL_"+decaystring).c_str()]->Fill(selectedCSVLJetID.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVM_"+decaystring).c_str()]->Fill(selectedCSVMJetID.size() , datasets[d], true,eventW);
    MSPlot[(prefix+"_nJetsCSVT_"+decaystring).c_str()]->Fill(selectedCSVTJetID.size() , datasets[d], true,eventW);
  
  }
  
    MSPlot[(prefix+"_Decay").c_str()]->Fill(channelInt , datasets[d], true,eventW);
  
  
  //cout << "end plot filling" << endl;
}


