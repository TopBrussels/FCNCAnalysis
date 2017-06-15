
//#define TreeAnalyzer_cxx
//#include "TreeAnalyzer.h"
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
#include "rochester/RoccoR.cc"
#include "TEfficiency.h"
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




struct HighestPt{
  bool operator()(TLorentzVector j1, TLorentzVector j2) const
  {
    return j1.Pt() > j2.Pt();
  }
  
  
};

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH1F*> histo1D_fakevvalidation;
map<string,TH1F*> histo1D_PUSystematics;
map<string,TH1F*> histo1D_ElSystematics;
map<string,TH1F*> histo1D_MuSystematics;
map<string,TH1F*> histo1D_Bcferr1Systematics;
map<string,TH1F*> histo1D_Bcferr2Systematics;
map<string,TH1F*> histo1D_Bhfstats1Systematics;
map<string,TH1F*> histo1D_Bhfstats2Systematics;
map<string,TH1F*> histo1D_Blfstats1Systematics;
map<string,TH1F*> histo1D_Blfstats2Systematics;
map<string,TH1F*> histo1D_BhfSystematics;
map<string,TH1F*> histo1D_BlfSystematics;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCutfl; // TOFIX


map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;


vector < Dataset* > datasets;
vector < Dataset* > datasetsbefore;
std::vector < int>  decayChannels = {0,1,2,3,-9,4,5}; // uuu uue eeu eee all
//std::vector < int>  decayChannels = {-9};
bool firstevent = false;



TFile* triggerEfffile = 0;
string triggerEfffilename = "triggerefficiencies.root";




TFile* charmscalefactorsfile = 0;
string charmscalefactorsfilename = "charmtagefficiencies";

int nbin_Pt_lep0 = 4;
int nbin_Pt_lep1 = 2;
int nbin_Pt_lep2 = 1;
int endbin_Pt_lep0 = 500;
int endbin_Pt_lep1 = 250;
int endbin_Pt_lep2 = 150;
int nbin_charmVSb = 21;
int nbin_charmVSl= 16;


///////////////////////////////////// MVA VARS /////////////////////////////////////////
// Decleration of MVA variables
TFile * tupfile = 0;
TTree* mvatree = 0;
TTree* mvatree_JECup = 0;
TTree* mvatree_JECdown = 0;
TTree* mvatree_JERup = 0;
TTree* mvatree_JERdown = 0;


Int_t MVA_channel = -999.;
Float_t MVA_weight = 1.;
Double_t MVA_weight_nom = 1.;
Double_t MVA_weight_puSF_up = 1.;
Double_t MVA_weight_puSF_down = 1.;
Double_t MVA_weight_electronSF_up =1.;
Double_t MVA_weight_electronSF_down = 1.;
Double_t MVA_weight_muonSF_up = 1.;
Double_t MVA_weight_muonSF_down = 1.;
Double_t MVA_weight_btagSF_cferr1_up = 1.;
Double_t MVA_weight_btagSF_cferr1_down = 1.;
Double_t MVA_weight_btagSF_cferr2_up = 1.;
Double_t MVA_weight_btagSF_cferr2_down = 1.;
Double_t MVA_weight_btagSF_hf_up = 1.;
Double_t MVA_weight_btagSF_hf_down = 1.;
Double_t MVA_weight_btagSF_hfstats1_up = 1.;
Double_t MVA_weight_btagSF_hfstats1_down = 1.;
Double_t MVA_weight_btagSF_hfstats2_up = 1.;
Double_t MVA_weight_btagSF_hfstats2_down = 1.;
Double_t MVA_weight_btagSF_lf_up = 1.;
Double_t MVA_weight_btagSF_lf_down = 1.;
Double_t MVA_weight_btagSF_lfstats1_up =1.;
Double_t MVA_weight_btagSF_lfstats1_down = 1.;
Double_t MVA_weight_btagSF_lfstats2_up = 1.;
Double_t MVA_weight_btagSF_lfstats2_down = 1.;
Double_t MVA_weight_puSF = 1. ;
Double_t MVA_weight_btagSF = 1.;
Double_t MVA_weight_muonSF = 1.;
Double_t MVA_weight_electronSF = 1.;

Int_t         MVA_id1;
Int_t         MVA_id2;
Double_t         MVA_x1;
Double_t         MVA_x2;
Double_t         MVA_q;
Double_t        MVA_weight0, MVA_weight1, MVA_weight2, MVA_weight3, MVA_weight4, MVA_weight5, MVA_weight6, MVA_weight7, MVA_weight8;
Double_t        MVA_hdamp_up;
Double_t        MVA_hdamp_down;

Float_t MVA_region = -999.;
Double_t MVA_EqLumi = -999.;
Double_t MVA_Luminosity = -999.;


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
Int_t MVA_nMuons = -999;
Float_t MVA_NJets_CSVv2T = -999;
Float_t MVA_NJets_CSVv2M = -999;
Float_t MVA_NJets_CSVv2L = -999;
Float_t MVA_nJets = -999;
Float_t MVA_nElectrons = -999;
Float_t MVA_nJets_CharmL = -999;
Float_t MVA_nJets_CharmM = -999;
Float_t MVA_nJets_CharmT = -999;


//SM kinematics
Float_t MVA_mWt = -999.;
Float_t MVA_mWt2 = -999.;
Float_t MVA_SMtop_M = -999.;
Float_t MVA_mlb = -999.;
Float_t MVA_Wboson_M = -999.;

Float_t MVA_dRWlepb = -999.;

Float_t MVA_dPhiWlepb = -999.;

Float_t MVA_Wlep_Charge = -999;
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
Double_t Luminosity = 36000.;
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

Double_t        x1;
Double_t        x2;
Int_t           id1;
Int_t           id2;
Double_t        q;
Double_t        hdamp_up;
Double_t        hdamp_down;

Double_t        btagSFshape;
Double_t        btagSFshape_down_cferr1;
Double_t        btagSFshape_down_cferr2;
Double_t        btagSFshape_down_hf;
Double_t        btagSFshape_down_hfstats1;
Double_t        btagSFshape_down_hfstats2;
Double_t        btagSFshape_down_lf;
Double_t        btagSFshape_down_lfstats1;
Double_t        btagSFshape_down_lfstats2;

Double_t        btagSFshape_up_cferr1;
Double_t        btagSFshape_up_cferr2;
Double_t        btagSFshape_up_hf;
Double_t        btagSFshape_up_hfstats1;
Double_t        btagSFshape_up_hfstats2;
Double_t        btagSFshape_up_lf;
Double_t        btagSFshape_up_lfstats1;
Double_t        btagSFshape_up_lfstats2;

Double_t weight0, weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8;


Int_t           channelInt;
Double_t        nloWeight;
Int_t           run_num;
Long64_t        evt_num;
Int_t           lumi_num;
Int_t           nvtx;
Int_t           npu;
Double_t        puSF;
Double_t        puSF_up;
Double_t        puSF_down;
Int_t           PassedMETFilter;
Int_t           PassedTrigger;
Int_t           PassedTriggerMET;
Int_t           PassedTriggerJET;
Int_t		PassedTriggerNoLogic;
Int_t		PassedTriggerNoLogic2;
Int_t           PassedGoodPV;
Int_t           nbOfLooseElectrons;
Int_t           nElectrons;
Double_t        ElectronSF[3];   //[nElectrons]
Double_t        ElectronSF_up[3];   //[nElectrons]
Double_t        ElectronSF_down[3];   //[nElectrons]
Double_t        pt_electron[3];   //[nElectrons]
Double_t        phi_electron[3];   //[nElectrons]
Double_t        eta_electron[3];   //[nElectrons]
Double_t        eta_superCluster_electron[3];   //[nElectrons]
Double_t        E_electron[3];   //[nElectrons]
Double_t        chargedHadronIso_electron[3];   //[nElectrons]
Double_t        neutralHadronIso_electron[3];   //[nElectrons]
Double_t        photonIso_electron[3];   //[nElectrons]
Double_t        pfIso_electron[3];   //[nElectrons]
Int_t           charge_electron[3];   //[nElectrons]
Double_t        d0_electron[3];   //[nElectrons]
Double_t        d0BeamSpot_electron[3];   //[nElectrons]
Double_t        sigmaIEtaIEta_electron[3];   //[nElectrons]
Double_t        deltaEtaIn_electron[3];   //[nElectrons]
Double_t        deltaPhiIn_electron[3];   //[nElectrons]
Double_t        hadronicOverEm_electron[3];   //[nElectrons]
Int_t           missingHits_electron[3];   //[nElectrons]
Bool_t          passConversion_electron[3];   //[nElectrons]
Bool_t          isId_electron[3];   //[nElectrons]
Bool_t          isIso_electron[3];   //[nElectrons]
Bool_t          isEBEEGap[3];   //[nElectrons]
Int_t           nMuons;
Int_t           nbOfLooseMuons;
Int_t           badmueventmu[3];
Int_t           rejecteventBadPFmuon;
Int_t           badmueventclonemu[3];
Double_t        MuonIDSF[3];   //[nMuons]
Double_t        ptSF_muon[3];
Double_t        MuonTrackSF[3];   //[nMuons]
Double_t        MuonIDSF_up[3];   //[nMuons]
Double_t        MuonIsoSF_up[3];   //[nMuons]
Double_t        MuonIDSF_down[3];   //[nMuons]
Double_t        MuonIsoSF_down[3];   //[nMuons]
Double_t        MuonIsoSF[3];   //[nMuons]
Double_t        MuonTrigSFv2[3];   //[nMuons]
Double_t        MuonTrigSFv3[3];   //[nMuons]
Double_t        pt_muon[3];   //[nMuons]
Double_t        TrackLayers_muon[3];
Double_t        phi_muon[3];   //[nMuons]
Double_t        eta_muon[3];   //[nMuons]
Double_t        E_muon[3];   //[nMuons]
Double_t        chargedHadronIso_muon[3];   //[nMuons]
Double_t        neutralHadronIso_muon[3];   //[nMuons]
Double_t        photonIso_muon[3];   //[nMuons]
Bool_t          isId_muon[3];   //[nMuons]
Bool_t          isIso_muon[3];   //[nMuons]
Double_t        pfIso_muon[3];   //[nMuons]
Int_t           charge_muon[3];   //[nMuons]
Double_t        d0_muon[3];   //[nMuons]
Double_t        d0BeamSpot_muon[3];   //[nMuons]
Int_t           nJets;
Double_t        pt_jet[9];   //[nJets]
Double_t        px_jet[9];   //[nJets]
Double_t        py_jet[9];   //[nJets]
Double_t        pz_jet[9];   //[nJets]
Double_t        phi_jet[9];   //[nJets]
Double_t        eta_jet[9];   //[nJets]
Double_t        E_jet[9];   //[nJets]
Int_t           charge_jet[9];   //[nJets]
Double_t        bdisc_jet[9];   //[nJets]
Double_t        jet_Pt_before_JER[9];   //[nJets]
Double_t        jet_Pt_before_JES[9];   //[nJets]
Double_t        jet_Pt_after_JER[9];   //[nJets]
Double_t        jet_Pt_after_JES[9];   //[nJets]
Double_t        cdiscCvsL_jet[9];   //[nJets]
Double_t        cdiscCvsB_jet[9];   //[nJets]
Int_t           nMCParticles;
Int_t           mc_status[133];   //[nMCParticles]
Int_t           mc_pdgId[133];   //[nMCParticles]
Int_t           mc_mother[133];   //[nMCParticles]
Int_t           mc_granny[133];   //[nMCParticles]
Double_t        mc_pt[133];   //[nMCParticles]
Double_t        mc_phi[133];   //[nMCParticles]
Double_t        mc_eta[133];   //[nMCParticles]
Double_t        mc_E[133];   //[nMCParticles]
Double_t        mc_M[133];   //[nMCParticles]
Bool_t          mc_isLastCopy[133];   //[nMCParticles]
Bool_t          mc_isPromptFinalState[133];   //[nMCParticles]
Bool_t          mc_isHardProcess[133];   //[nMCParticles]
Bool_t          mc_fromHardProcessFinalState[133];   //[nMCParticles]
Double_t        met_Pt;
Double_t        met_Eta;
Double_t        met_Phi;
Double_t        met_Px;
Double_t        met_Py;
Double_t        met_before_JES;
Double_t        met_after_JES;

// List of branches
TBranch         *b_weight0, *b_weight1, *b_weight2, *b_weight3, *b_weight4, *b_weight5, *b_weight6, *b_weight7, *b_weight8; 
TBranch        *b_x1;   //!
TBranch        *b_x2;   //!
TBranch        *b_id1;   //!
TBranch        *b_id2;   //!
TBranch        *b_q;   //!
TBranch        *b_hdamp_up;   //!
TBranch        *b_hdamp_down;   //!
TBranch        *b_btagSFshape;   //!
TBranch        *b_btagSFshape_down_cferr1;   //!
TBranch        *b_btagSFshape_down_cferr2;   //!
TBranch        *b_btagSFshape_down_hf;   //!
TBranch        *b_btagSFshape_down_hfstats1;   //!
TBranch        *b_btagSFshape_down_hfstats2;   //!
TBranch        *b_btagSFshape_down_lf;   //!
TBranch        *b_btagSFshape_down_lfstats1;   //!
TBranch        *b_btagSFshape_down_lfstats2;   //!

TBranch        *b_btagSFshape_up_cferr1;   //!
TBranch        *b_btagSFshape_up_cferr2;   //!
TBranch        *b_btagSFshape_up_hf;   //!
TBranch        *b_btagSFshape_up_hfstats1;   //!
TBranch        *b_btagSFshape_up_hfstats2;   //!
TBranch        *b_btagSFshape_up_lf;   //!
TBranch        *b_btagSFshape_up_lfstats1;   //!
TBranch        *b_btagSFshape_up_lfstats2;   //!

TBranch        *b_channelInt;   //!
TBranch        *b_nloWeight;   //!
TBranch        *b_run_num;   //!
TBranch        *b_evt_num;   //!
TBranch        *b_lumi_num;   //!
TBranch        *b_nvtx;   //!
TBranch        *b_npu;   //!
TBranch        *b_puSF;   //!
TBranch        *b_puSF_down;   //!
TBranch        *b_puSF_up;   //!
TBranch        *b_PassedMETFilter;   //!
TBranch        *b_PassedTrigger;   //!
TBranch        *b_PassedTriggerJET;   //!
TBranch        *b_PassedTriggerMET;
TBranch        *b_PassedTriggerNoLogic;
TBranch		     *b_PassedTriggerNoLogic2;
TBranch        *b_PassedGoodPV;   //!
TBranch         *b_nbOfLooseElectrons;
TBranch        *b_nElectrons;   //!
TBranch        *b_ElectronSF;   //!
TBranch        *b_ElectronSF_up;   //!
TBranch        *b_ElectronSF_down;   //!
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
TBranch        *b_nbOfLooseMuons;
TBranch        *b_nMuons;   //!
TBranch        *b_MuonIDSF;   //!
TBranch        *b_ptSF_muon;
TBranch        *b_MuonTrackSF;
TBranch        *b_MuonIsoSF;   //!
TBranch        *b_MuonIDSF_up;   //!
TBranch        *b_MuonIsoSF_up;   //!
TBranch        *b_MuonIDSF_down;   //!
TBranch        *b_MuonIsoSF_down;   //!
TBranch        *b_MuonTrigSFv2;   //!
TBranch        *b_MuonTrigSFv3;   //!
TBranch        *b_pt_muon;   //!
TBranch        *b_badmueventclonemu;
TBranch        *b_rejecteventBadPFmuon;
TBranch        *b_badmueventmu;
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
TBranch        *b_TrackLayers_muon;
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
TBranch        *b_nMCParticles;   //!
TBranch        *b_mc_status;   //!
TBranch        *b_mc_pdgId;   //!
TBranch        *b_mc_mother;   //!
TBranch        *b_mc_granny;   //!
TBranch        *b_mc_pt;   //!
TBranch        *b_mc_phi;   //!
TBranch        *b_mc_eta;   //!
TBranch        *b_mc_E;   //!
TBranch        *b_mc_M;   //!
TBranch        *b_mc_isLastCopy;   //!
TBranch        *b_mc_isPromptFinalState;   //!
TBranch        *b_mc_isHardProcess;   //!
TBranch        *b_mc_fromHardProcessFinalState;   //!
TBranch        *b_met_Pt;   //!
TBranch        *b_met_Eta;   //!
TBranch        *b_met_Phi;   //!
TBranch        *b_met_Px;   //!
TBranch        *b_met_Py;   //!
TBranch        *b_met_before_JES;   //!
TBranch        *b_met_after_JES;   //!
///////////////////////////////////// HOME MADE FUNCTIONS AND VARS  /////////////////////////////////////////


Double_t met_Pz;
int verbose = 2;


string MakeTimeStamp();
void InitMSPlots(string prefix, vector<int> decayChannels);
void InitMVAMSPlotsSingletop(string prefix,vector<int> decayChannels);
void InitMSPlotsBDT(string prefix, vector<int> decayChannels);
void InitMVAMSPlotsTopPair(string prefix, vector<int> decayChannels);
void InitMVAMSPlotsWZ(string prefix, vector <int> decayChannels);
void InitFakeValidation(string dataSetName, vector <int> decayChannels);
void Init1DPlots(string dataSetName);
void Init2DPlots();
void InitGenInfoPlots(string dataSetName);
void FillGenInfoPlots(string dataSetName);
void InitRecovsGenInfoPlots(string dataSetName);
void FillRecovsGenInfoPlots(string dataSetName, vector<TLorentzVector> selectedElectrons, vector <TLorentzVector> selectedMuons , vector <TLorentzVector> selectedJets);
void FillFakeValidation(string dataSetName, vector <int> decayChannels, bool isData, bool isfakes, bool threelepregion,bool twolepregion);
void Fill1DPlots(string dataSetName, double eventW, bool threelepregion, bool twolepregion);
void InitTree(TTree* tree, bool isData, bool isfakes);
// data from global tree
void ClearMetaData();
void GetMetaData(TTree* tree, bool isData,int Entries, bool isAMC, bool isfakes);
// put everything to default values
void ClearObjects();
void ClearVars();
void ClearMVAVars();
void ClearTLVs();
void ClearMatchingVars();
void ClearMatchingVarsTLV();
void ClearMatchingSampleVars();
void FillGeneralPlots(int d, string prefix, vector<int>decayChannels, bool isData, bool isfakes,bool threelepregion,bool twolepregion);
void FillMVAPlots(int d, string dataSetName, int Region, string prefix, vector<int>decayChannels);
string ConvertIntToString(int nb, bool pad);
void ReconstructObjects(vector<int> selectedJetsID, vector<TLorentzVector> Muons,vector<TLorentzVector> selectedElectrons, vector<TLorentzVector> selectedJets,int Region, bool threelepregion);
void MakeMVAvars(int Region, Double_t scaleFactor);
void createMVAtree(string dataSetName);
void writeMVAtree();
int SMjetCalculator(vector<TLorentzVector> Jets,int verb);
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met, TLorentzVector SMjet);
int FCNCjetCalculator(vector < TLorentzVector>  Jets, TLorentzVector recoZ ,int index, int verb);
int FCNCjetCalculatorCvsBTagger(vector < TLorentzVector>  Jets, int index, int verb);
int FCNCjetCalculatorCvsLTagger(vector < TLorentzVector>  Jets, int index, int verb);
int FCNCjetCalculatorCwp(vector < TLorentzVector>  Jets, std::vector <int> cjetindex, int index, int verb);
std::pair <Double_t,Double_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo);
void LeptonAssigner(vector<TLorentzVector> electrons, vector<TLorentzVector> muons, std::vector<int> electronsCharge,std::vector<int> muonsCharge);
bool MatchingFunction(string dataSetName, vector <TLorentzVector> selectedleptons, vector <TLorentzVector> selectedMuons,vector <TLorentzVector> selectedElectrons, vector <TLorentzVector> selectedJets, bool makePlots, bool debug);
void LeptonMatcher(vector < TLorentzVector> mcParticles, vector <TLorentzVector> Leptons);
void EventSearcher(vector < TLorentzVector> mcParticles, string dataSetName, bool debug);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> > LeptonMatching(vector < TLorentzVector> selectedleptons, vector <TLorentzVector> mcParticles, string dataSetName, bool debug);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> > JetMatching(vector < TLorentzVector> selectedJets, vector <TLorentzVector> mcParticles, string dataSetName, bool debug);
void MatchingEfficiency();
Double_t RochLeptonMatching(TLorentzVector selectedlepton, vector <TLorentzVector> mcParticles, bool isData, double nbtracks, int chargelep, bool isNP, bool debug);
vector<TLorentzVector> selectedMuons;
vector<TLorentzVector> selectedElectrons;
vector<TLorentzVector> selectedLeptons;
vector<TLorentzVector> selectedJets;
vector<TLorentzVector> preselectedJets;
vector<int> selectedElectronsCharge;
vector<int> selectedMuonsCharge;
vector<int> electronID;
vector<int> muonID;
bool Assigned = false;

int WmuIndiceM = -999;
int WelecIndiceM = -999;
int BjetIndiceM = -999;
int CjetIndiceM = -999;
int UjetIndiceM = -999;
int ZelecIndiceM_0 = -999;
int ZelecIndiceM_1 = -999;
int ZmuIndiceM_0 = -999;
int ZmuIndiceM_1 = -999;
Double_t WlepMatched= 0.;
Double_t ZlepMatched = 0.;
Double_t BjetMatched = 0.;
Double_t CjetMatched = 0.;
Double_t CjetMatchedeventCvsLL = 0.;
Double_t CjetMatchedeventCvsLM = 0.;
Double_t CjetMatchedeventCvsLT = 0.;
Double_t CjetMatchedeventCvsBL = 0.;
Double_t CjetMatchedeventCvsBM = 0.;
Double_t CjetMatchedeventCvsBT = 0.;
Double_t CjetMatchedeventL = 0.;
Double_t CjetMatchedeventM = 0.;
Double_t CjetMatchedeventT = 0.;
Double_t CjetMatchedL = 0.;
Double_t CjetMatchedM = 0.;
Double_t CjetMatchedT = 0.;
Double_t CjetMatchedCvsBL = 0.;
Double_t CjetMatchedCvsBM = 0.;
Double_t CjetMatchedCvsBT = 0.;
Double_t CjetMatchedCvsLL = 0.;
Double_t CjetMatchedCvsLM = 0.;
Double_t CjetMatchedCvsLT = 0.;

Double_t UjetMatched = 0.;
Double_t WlepMatchedevent= 0.;
Double_t ZlepMatchedevent = 0.;
Double_t BjetMatchedevent = 0.;
Double_t CjetMatchedevent = 0.;
Double_t UjetMatchedevent = 0.;

TLorentzVector mcpart;
vector <TLorentzVector> mcParticles;
vector <TLorentzVector> mcParticlesroch;
vector <TLorentzVector> partonsrochester;
vector <TLorentzVector> selectedleps;
RoccoR rc("rochester/rcdata.2016.v3");
bool foundTopQ = false;
bool foundAntitopQ = false;
bool foundSMb = false;
bool foundSMmu = false;
bool foundSMel = false;
bool foundSMnuel = false;
bool foundSMnumu = false;
bool foundW = false;
bool foundZmumin = false;
bool foundZmuplus = false;
bool foundZelmin = false;
bool foundZelplus = false;
bool foundZ = false;
bool foundCjet = false;
bool foundUjet = false;
int TopQ_Indice = -999;
int AntiTopQ_Indice = -999;
int SMb_Indice = -999;
int SMmu_Indice = -999;
int SMnuel_Indice = -999;
int SMnumu_Indice = -999;
int SMel_Indice = -999;
int SMW_Indice = -999;
int Zmu_min_Indice = -999;
int Zmu_plus_Indice = -999;
int Zel_min_Indice = -999;
int Zel_plus_Indice = -999;
int Z_Indice = -999;
int Cjet_Indice = -999;
int Ujet_Indice = -999;
bool foundDecay = false;

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
vector <int> selectedJetsID;
Double_t scaleFactor;
Double_t scaleFactor_bfBT;
Double_t scaleFactor_bfELSF;
Double_t scaleFactor_bfMuSF;
Double_t scaleFactor_bfPU;
Double_t scaleFactor_puSF_down;
Double_t scaleFactor_puSF_up;
Double_t scaleFactor_electronSF_down;
Double_t scaleFactor_electronSF_up;
Double_t scaleFactor_muonSF_down;
Double_t scaleFactor_muonSF_up;
Double_t scaleFactor_btagSF_cferr1_down;
Double_t scaleFactor_btagSF_cferr1_up;
Double_t scaleFactor_btagSF_cferr2_down;
Double_t scaleFactor_btagSF_cferr2_up;
Double_t scaleFactor_btagSF_hf_down;
Double_t scaleFactor_btagSF_hf_up;
Double_t scaleFactor_btagSF_hfstats1_down;
Double_t scaleFactor_btagSF_hfstats1_up;
Double_t scaleFactor_btagSF_hfstats2_down;
Double_t scaleFactor_btagSF_hfstats2_up;
Double_t scaleFactor_btagSF_lf_down;
Double_t scaleFactor_btagSF_lf_up;
Double_t scaleFactor_btagSF_lfstats1_down;
Double_t scaleFactor_btagSF_lfstats1_up;
Double_t scaleFactor_btagSF_lfstats2_down;
Double_t scaleFactor_btagSF_lfstats2_up;
Double_t scaleFactor_puSF ;
Double_t scaleFactor_btagSF ;
Double_t scaleFactor_muonSF;
Double_t scaleFactor_electronSF;

Double_t muonSFtemp;
Double_t electronSFtemp;
Double_t EquilumiSF = 1.;
Double_t nloSF = 1.;
Double_t Xsect = 1.;
Double_t mWT;
Double_t mWT2;
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

bool check_matching = false;
bool debugmatching = false;
bool doDilep = false;
bool noTrilep = false;
Double_t tempPx;
Double_t tempPy;
Double_t tempHt;
Double_t tempPx_jet;
Double_t tempPy_jet;
Double_t tempHt_jet;

TLorentzVector tempInvMassObj;
TLorentzVector tempInvMassObj_jet;

string pathOutputdate = "";

///////////////////////////////////// MAIN CODE /////////////////////////////////////////
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
  datasetsbefore.clear();
  
  
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasetsbefore, xmlFile);
  cout << "datasets loaded" <<endl;
  bool makePlots = false;
  bool makeMVAPlots = false;
  bool makeMatchingPlots = false;
  bool makeMVAtree = false;
  bool applyMuonSF = false;
  bool applyMuonSF_down = false;
  bool applyMuonSF_up = false;
  bool applyJEC_up = false;
  bool applyJEC_down = false;
  bool applyJER_up = false;
  bool applyJER_down = false;
  bool applyElectronSF = false;
  bool applyElectronSFup = false;
  bool applyElectronSFdown = false;
  bool applyBTagSF = false; // To fix
  bool applyPUSF = false;
  bool applyPUSF_down = false;
  bool applyPUSF_up = false;
  bool applyNloSF = false;
  bool applyMETfilter = false;
  bool removeBadMu  = true;
  bool doRunGH = false;
  bool doRunBCDEF = false;
  string placeNtup = "singletop/170214";
  int channel = -999;
  datafound = false;
  bool applytrigger = true;
  bool docharmsf = false;
  bool applycharmsf = false;
  bool checktrigger = false;
  bool applytriggerNoLogic = false;
  bool applytriggerNoLogic2 = false;
  bool testing = false;
  int testevent = 100;
  bool dorochester = false;
  bool dofakevalidation = false;
  bool applyTrigSF = false;
  bool doCutflow = false;
  bool MakeSelectionTable = false;
  bool systematicplots = false;
  doDilep = false;
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./Analyzer [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   decay uuu : add decay channel uuu, default value is all decays: uuu / uue / eeu / eee / all = 0 / 1 / 2 / 3 / -9" << endl;
      std::cout << "   debug: set verbose to high value" << endl;
      std::cout << "   applySF: applyPUSF, applyElectronSF, applyMuonsSF, applyMET, applyBtagSF, applyAMC. These can also be set seperatly" << endl;
      std::cout << "   - applyElectronSF UP/DOWN/NOM "<< endl;
      std::cout << "   - applyMuonSF UP/DOWN/NOM "<< endl;
      std::cout << "   - applyPUSF UP/DOWN/NOM "<< endl;
      std::cout << "   - applyJER UP/DOWN" << endl;
      std::cout << "   - applyJEC UP/DOWN" << endl;
      std::cout << "   MakeMVATree: make trees for the MVA" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   MakeMVAPlots: make MVA plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. NtupleMakerOutput/MergedTuples/placeNtup " << endl;
      std::cout << "   Matching: do matching for tZq, FCNC" << endl;
      std::cout << "   debugMatching: do matching for tZq, FCNC" << endl;
      std::cout << "   MakeMatchingPlots: make matching plots" << endl;
      std::cout << "   noTrigger: do not apply trigger " << endl;
      std::cout << "   checkTrigger: check trigger in data " << endl;
      std::cout << "   Test: loop over 1000 events" << endl;
      std::cout << "   withBadMu " << endl;
      std::cout << "   RunGH / RunBF " << endl;
      std::cout << "   doDilep: include dilep plots " << endl;
      std::cout << "   noTrlep: exclude trilep plots " << endl;
      std::cout << "   docharmSF: make charm SF " << endl;
      std::cout << "   applycharmSF: apply charm SF " << endl;
      std::cout << "   rochester: apply rochester" << endl;
      std::cout << "   fakeval: apply fakevalidation" << endl;
      std::cout << "   applyTrigSF" << endl;
      std::cout << "   doCutFlow" << endl;
      std::cout << "   doCutTable" << endl;
       std::cout << "   doSys" << endl;
      return 0;
    }
    if(string(argv[i]).find("doSys")!=std::string::npos) {
      systematicplots= true;
      
    }
    if(string(argv[i]).find("doCutFlow")!=std::string::npos) {
      doCutflow = true;
      
    }
    if(string(argv[i]).find("doCutTable")!=std::string::npos) {
      MakeSelectionTable = true;
      
    }
    if(string(argv[i]).find("applyTrigSF")!=std::string::npos) {
      applyTrigSF = true;
      
    }
    if(string(argv[i]).find("fakeval")!=std::string::npos) {
      dofakevalidation = true;
      
    }
    if(string(argv[i]).find("rochester")!=std::string::npos) {
      dorochester = true;
      
    }
    if(string(argv[i]).find("applyCharmSF")!=std::string::npos) {
      applycharmsf = true;
      docharmsf = false;
    }
    if(string(argv[i]).find("doCharmSF")!=std::string::npos) {
      docharmsf = true;
      applycharmsf = false;
    }
    if(string(argv[i]).find("noTrilep")!=std::string::npos) {
      noTrilep = true;
    }
    if(string(argv[i]).find("doDilep")!=std::string::npos) {
      doDilep = true;
    }
    if(string(argv[i]).find("RunGH")!=std::string::npos) {
      doRunGH = true;
    }
    if(string(argv[i]).find("RunBF")!=std::string::npos) {
      doRunBCDEF = true;
    }
    if(string(argv[i]).find("withBadMu")!=std::string::npos) {
      removeBadMu = false;
    }
    if(string(argv[i]).find("applyJEC")!=std::string::npos) {
      
      i++;
      if(string(argv[i]).find("UP")!=string::npos) applyJEC_up = true;
      else if(string(argv[i]).find("DOWN")!=string::npos) applyJEC_down= true;
      else if(string(argv[i]).find("NOM")!=string::npos) applyJEC_down = applyJEC_up = false;
      else { cout << "argument missing for JEC" << endl; break; }
    }
    if(string(argv[i]).find("applyJER")!=std::string::npos) {
      
      i++;
      if(string(argv[i]).find("UP")!=string::npos) applyJER_up = true;
      else if(string(argv[i]).find("DOWN")!=string::npos) applyJER_down= true;
      else if(string(argv[i]).find("NOM")!=string::npos) applyJER_down = applyJER_up = false;
      else { cout << "argument missing for JER" << endl; break; }
    }
    if(string(argv[i]).find("Testing")!=std::string::npos) {
      testing = true;
      i++;
      testevent = strtol(argv[i], NULL, 10);
    }
    if(string(argv[i]).find("checkTrigger")!=std::string::npos) {
      checktrigger = true;
      applytrigger = false;
    }
    if(string(argv[i]).find("noTrigger")!=std::string::npos) {
      applytrigger = false;
    }
    if(string(argv[i]).find("checkTrigNoLogic")!=std::string::npos) {
      applytriggerNoLogic = true;
    }
    if(string(argv[i]).find("checkNoLogic2")!=std::string::npos) {
      applytriggerNoLogic2 = true;
    }
    if(string(argv[i]).find("MakeMatchingPlots")!=std::string::npos) {
      makeMatchingPlots = true;
    }
    if(string(argv[i]).find("Matching")!=std::string::npos) {
      check_matching = true;
    }
    if(string(argv[i]).find("debugMatching")!=std::string::npos) {
      debugmatching = true;
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
   
    if(string(argv[i]).find("applySF")!=string::npos) {
      applyElectronSF =true;
      applyNloSF = true;
      applyMETfilter = true;
      applyMuonSF = true;
      applyPUSF = true;
      applyBTagSF = true;
      cout << "                applying scalefactors " << endl;
    }
    if(string(argv[i]).find("noElectronSF")!=string::npos) {
      applyElectronSF =false;
         }
    if(string(argv[i]).find("noMuonSF")!=string::npos) {
      applyMuonSF =false;
     
      
    }
    if(string(argv[i]).find("noPUSF")!=std::string::npos) {
      applyPUSF =false;
    }
    if(string(argv[i]).find("noMET")!=string::npos) {
      applyMETfilter =false;
    }
    if(string(argv[i]).find("noBtagSF")!=string::npos) {
      applyBTagSF =false;
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
    if(string(argv[i]).find("MakeMVAPlots")!=string::npos) {
      makeMVAPlots = true;
    }
  }
  
  
  if(!applyNloSF || !applyBTagSF ||! applyElectronSF || !applyMETfilter ||  !applyMuonSF || !dorochester) cout << " WARNING not all booleans set for reweighing" << endl;
  
   if(docharmsf && doDilep) charmscalefactorsfilename = charmscalefactorsfilename + "_" + dateString + "_dilep" ;
  else if(applycharmsf && doDilep) charmscalefactorsfilename = charmscalefactorsfilename + "_170615_1405_dilep" ;
   charmscalefactorsfilename = charmscalefactorsfilename + ".root";
  if(!makePlots) doCutflow = false;
  for (int d = 0; d < datasetsbefore.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasetsbefore[d]->Name();
    
    if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data") !=std::string::npos || dataSetName.find("DATA") !=std::string::npos)
    {
      if(doRunGH && dataSetName.find("RunGH")==std::string::npos){cout << "looking at run GH and this is " << dataSetName << endl;  continue;}
      else if(doRunBCDEF && dataSetName.find("RunBCDEF")==std::string::npos){cout << "looking at run B-F and this is " << dataSetName << endl;  continue;}
      else { Luminosity = datasetsbefore[d]->EquivalentLumi();
        cout << "setting lumi to " << Luminosity << endl;
      }
    }
    if((applyJEC_down || applyJEC_up || applyJER_down || applyJER_up) && dataSetName.find("data")!=std::string::npos){continue;}
    else if((applyJEC_down || applyJEC_up || applyJER_down || applyJER_up) && dataSetName.find("fake")!=std::string::npos){continue;}
    else{
      datasets.push_back(datasetsbefore[d]);
      cout << " looking at " << dataSetName << endl;
    }
    
  }
  if(makePlots){
    firstevent = true;
    //InitMSPlots("control_afterAtLeast1Jet", decayChannels);
    InitMSPlots("control_afterAtLeast1Jet_afterZWindow", decayChannels);
    // InitMSPlots("control_afterAtLeast1Jet_afterZWindow_afterAtLeast1BJet", decayChannels);
    // Init1DPlots();
    MSPlot["cutflow"] = new MultiSamplePlot(datasets, "cutflow", 10, -0.5, 9.5, "Cutflow");
    MSPlot["cutflow_eee"] = new MultiSamplePlot(datasets, "cutflow_eee", 10, -0.5, 9.5, "Cutflow");
    MSPlot["cutflow_eeu"] = new MultiSamplePlot(datasets, "cutflow_eeu", 10, -0.5, 9.5, "Cutflow");
    MSPlot["cutflow_uue"] = new MultiSamplePlot(datasets, "cutflow_uue", 10, -0.5, 9.5, "Cutflow");
    MSPlot["cutflow_uuu"] = new MultiSamplePlot(datasets, "cutflow_uuu", 10, -0.5, 9.5, "Cutflow");
    Init2DPlots();
  }
  vector < string > v_cutflow = {">1l,>0j", "SF pair","lep veto","Z mass",">2l","STSR","TTSR","WZCR"};
  SelectionTable *CutflowTable = new SelectionTable(v_cutflow,datasets);
  SelectionTable *CutflowTable_eee = new SelectionTable(v_cutflow,datasets);
  SelectionTable *CutflowTable_eeu = new SelectionTable(v_cutflow,datasets);
  SelectionTable *CutflowTable_uue = new SelectionTable(v_cutflow,datasets);
  SelectionTable *CutflowTable_uuu = new SelectionTable(v_cutflow,datasets);
  if(MakeSelectionTable){
    CutflowTable->SetLuminosity(Luminosity);
    CutflowTable->SetPrecision(1);
    
    CutflowTable_eee->SetLuminosity(Luminosity);
    CutflowTable_eee->SetPrecision(1);
    
    CutflowTable_eeu->SetLuminosity(Luminosity);
    CutflowTable_eeu->SetPrecision(1);
    
    CutflowTable_uue->SetLuminosity(Luminosity);
    CutflowTable_uue->SetPrecision(1);
    
    CutflowTable_uuu->SetLuminosity(Luminosity);
    CutflowTable_uuu->SetPrecision(1);
  }
  
  if(makeMVAtree && makeMVAPlots) {
   // InitMVAMSPlotsSingletop("singletop", decayChannels);
    InitMVAMSPlotsTopPair("wzcontrol", decayChannels);
    // InitMVAMSPlotsTopPair("ttzcontrol", decayChannels);
    //InitMVAMSPlotsTopPair("toppair", decayChannels);
    
  }
  
  
  /*
   
   for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
   {
   string dataSetName = datasets[d]->Name();
   if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
   {
   Luminosity = datasets[d]->EquivalentLumi();
   datafound = true;
   }
   
   
   
   }
   */
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  bool isData = false;
  bool isfakes = false;
  bool isAMC = false;
  bool isNP = false;
  int nSelectedEntriesST = 0;
  int nSelectedEntriesTT = 0;
  int nSelectedEntriesWZ = 0;
  int nSelectedEntriesTTZ = 0;
  int nSelectedEntriesDilep = 0;
  Double_t nSelectedEntriesDilepweighted = 0;
  Double_t nSelectedEntriesSTweighted = 0.;
  Double_t nSelectedEntriesTTweighted = 0.;
  Double_t nSelectedEntriesWZweighted = 0.;
  Double_t nSelectedEntriesTTZweighted = 0.;
  Long64_t evtID ;
  ofstream myfile;
  ofstream myfiletrigged;
  
  ofstream myfileST;
  ofstream myfileSTtrigged;
  
  ofstream myfileTT;
  ofstream myfileTTtrigged;
  
  ofstream myfileWZ;
  ofstream myfileWZtrigged;
  
  
  TH1::SetDefaultSumw2();
  TH1F* CvsB_Histo_sum = 0;
  TH1F* CvsL_Histo_sum = 0;
  TH1F* CvsB_Histo_data = 0;
  TH1F* CvsL_Histo_data = 0;
  TH1F* SumNormal_cvsb =0 ;
  TH1F* SumNormal_cvsl=0 ;
  TH1F* dataNormal_cvsb = 0;
  TH1F* dataNormal_cvsl =0 ;
  TH1F* charm_SFHisto_cvsb =0 ;
  TH1F* charm_SFHisto_cvsl=0 ;
  TH1F* charm_SFHisto_cvsb_up =0 ;
  TH1F* charm_SFHisto_cvsl_up=0 ;
  TH1F* charm_SFHisto_cvsb_down =0 ;
  TH1F* charm_SFHisto_cvsl_down =0 ;
  
  if(docharmsf){
    CvsB_Histo_sum  = new TH1F("CvsB_Histo_sum", "CvsB_Histo_sum" , nbin_charmVSb,-1, 1);
    CvsL_Histo_sum  = new TH1F("CvsL_Histo_sum", "CvsL_Histo_sum" , nbin_charmVSl,-0.6, 1);
    
    CvsB_Histo_data  = new TH1F("CvsB_Histo_data", "CvsB_Histo_data" , nbin_charmVSb,-1, 1);
    CvsL_Histo_data = new TH1F("CvsL_Histo_data", "CvsL_Histo_data" , nbin_charmVSl,-0.6, 1);
  }
  
  if(applycharmsf){
    charmscalefactorsfile = TFile::Open( charmscalefactorsfilename.c_str(), "READ" );
    charm_SFHisto_cvsb = (TH1F*) (charmscalefactorsfile->Get("charm_SFHisto_cvsb"))->Clone("charm_SFHisto_cvsb");
    charm_SFHisto_cvsl = (TH1F*) (charmscalefactorsfile->Get("charm_SFHisto_cvsl"))->Clone("charm_SFHisto_cvsl");
   
    
  }
  // cout << " nbins " << MuPtSFHisto0->GetNbinsX() << endl;
  int xbinmcharm=-21;


  int endEvent = -5;
  
  
  
  /// Loop over datasets
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    Long64_t  mineventnb = 999999999;
    Long64_t maxeventnb = 0;
    
    firstevent = true;
    ClearMetaData();
    //cout << "meta data cleared" << endl;
    dataSetName = datasets[d]->Name();
    Xsect = datasets[d]->Xsection();
    //if(dataSetName.find("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X")==std::string::npos) continue; // TO FIX
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false;
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data")!= std::string::npos || dataSetName.find("DATA")!= std::string::npos )
    {
      isData = true;
      cout <<" found data " << endl;
    }
    isfakes = false;
    if(dataSetName.find("fake")!=std::string::npos ) {isfakes = true;}
    isAMC= false;
    if ( dataSetName.find("amc")!= std::string::npos || dataSetName.find("AMC") != std::string::npos  )
    {
      isAMC = true;
      //cout << "amc at nlo sample" <<endl;
    }
    isNP = false;
    if (dataSetName.find("FCNC")!=std::string::npos)
    {
      isNP = true;
    }
    if (dataSetName.find("FCNC")==std::string::npos && dataSetName.find("tZq")==std::string::npos)
    {
      check_matching = false;
      
    }
    
    
    if(check_matching){
      ClearMatchingSampleVars();
      if(makeMatchingPlots){
        InitGenInfoPlots(dataSetName);
        InitRecovsGenInfoPlots(dataSetName);
      }
    }
    if(( dataSetName.find("TT")!=std::string::npos || dataSetName.find("WWTo")!=std::string::npos||dataSetName.find("DY")!=std::string::npos || dataSetName.find("Zjets")!=std::string::npos  || dataSetName.find("data")!=std::string::npos || dataSetName.find("fake")!=std::string::npos ) && dofakevalidation  ){
      InitFakeValidation(dataSetName, decayChannels);
    }
    //    if(dataSetName.find("TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut")!=std::string::npos){
    // if(dataSetName.find("WZTo3LNu")!=std::string::npos){
    // if(dataSetName.find("tZq")!=std::string::npos){
    if(dataSetName.find("WZTo3LNu")!=std::string::npos && systematicplots){
      //cout << "init 1D plots" << endl;
      Init1DPlots(dataSetName);
    }
    
    string ntupleFileName = "NtupleMakerOutput/MergedTuples/"+placeNtup+"/"+dataSetName+".root";
    tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "tree";
    string tStatsTreeName = "globaltree";
    string postfix = "";
    if(dataSetName.find("fakeshift")!=std::string::npos) postfix ="_FakeShift";
    if(isData || isfakes) postfix = "";
    if(applyJEC_down) postfix = "_JESdown";
    if(applyJEC_up) postfix = "_JESup";
    if(applyJER_down) postfix = "_JERdown";
    if(applyJER_up) postfix = "_JERup";
    
    tTreeName += postfix;
    tStatsTreeName +=  postfix;
    //cout << "looking at " << tStatsTreeName << " and " << tTreeName << endl;
    /// Get meta data
    tStatsTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tStatsTreeName.c_str());
    globalnEntries = (int) tStatsTree[dataSetName.c_str()]->GetEntries();
    //cout << "getting meta data " << endl;
    GetMetaData(tStatsTree[dataSetName.c_str()], isData, globalnEntries, isAMC, isfakes);
    //cout << "meta data gotten" << endl;
    
    
    /// Get data
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries: " << nEntries << endl;
    
    
    // Set branch addresses and branch pointers
    InitTree(tTree[dataSetName.c_str()], isData, isfakes);
    
    if(makeMVAtree){
      
      TString output_file_name = pathOutputdate+"/MVAtrees/";
      mkdir(output_file_name, 0777);
      output_file_name = pathOutputdate+"/MVAtrees/MVA_tree_" + dataSetName + postfix + ".root";
      cout << "                making MVA tree file with name  " << output_file_name << endl;
      tupfile = new TFile(output_file_name,"RECREATE");
      mvatree = new TTree(("mvatree"+postfix).c_str(), ("mvatree"+postfix).c_str());
      
      createMVAtree(dataSetName);
      
    }
    
    if( (isData || dataSetName.find("WZ")!=std::string::npos) && checktrigger ){
      myfile.open((dataSetName+"eventID.txt").c_str());
      myfiletrigged.open((dataSetName+"eventIDtrigged.txt").c_str()); // sort -u eventID.txt | wc -l  (-u ==> unique lines )
      
      myfileST.open((dataSetName+"eventID_ST.txt").c_str());
      myfileSTtrigged.open((dataSetName+"eventIDtrigged_ST.txt").c_str()); // sort -u eventID.txt | wc -l  (-u ==> unique lines )
      
      myfileTT.open((dataSetName+"eventID_TT.txt").c_str());
      myfileTTtrigged.open((dataSetName+"eventIDtrigged_TT.txt").c_str()); // sort -u eventID.txt | wc -l  (-u ==> unique lines )
      
      myfileWZ.open((dataSetName+"eventID_WZ.txt").c_str());
      myfileWZtrigged.open((dataSetName+"eventIDtrigged_WZ.txt").c_str()); // sort -u eventID.txt | wc -l  (-u ==> unique lines )
      
      
    }
    
    
    endEvent =nEntries;
    if(testing){
      if(endEvent > 1000) endEvent = testevent;
    }
    int istartevt = 0;
    nSelectedEntriesST = 0;
    nSelectedEntriesTT = 0;
    nSelectedEntriesWZ = 0;
    nSelectedEntriesSTweighted = 0.;
    nSelectedEntriesTTweighted = 0.;
    nSelectedEntriesWZweighted = 0.;
    nSelectedEntriesDilep = 0;
    nSelectedEntriesTTZ = 0;
    nSelectedEntriesTTZweighted = 0.;
    nSelectedEntriesDilepweighted = 0.;
    
    // for trig eff
    TH1::SetDefaultSumw2();
    TH1F* histPt_all = 0;
    TH1F* histPt_noTrig_all = 0;
    TH1F* histPt_leadinglep_all = (0);
    TH1F* histPt_2ndleadinglep_all = (0);
    TH1F* histPt_3dleadinglep_all = (0);
    TH1F* histPt_leadinglep_noTrig_all = (0);
    TH1F* histPt_2ndleadinglep_noTrig_all = (0);
    TH1F* histPt_3dleadinglep_noTrig_all = (0);
    
    TH1F* histPt_3mu = 0;
    TH1F* histPt_noTrig_3mu =  0;
    TH1F* histPt_leadinglep_3mu = (0);
    TH1F* histPt_2ndleadinglep_3mu = (0);
    TH1F* histPt_3dleadinglep_3mu = (0);
    TH1F* histPt_leadinglep_noTrig_3mu = (0);
    TH1F* histPt_2ndleadinglep_noTrig_3mu = (0);
    TH1F* histPt_3dleadinglep_noTrig_3mu = (0);
    
    
    TH1F* histPt_3e =  0;
    TH1F* histPt_noTrig_3e = 0;
    TH1F* histPt_leadinglep_3e =  (0);
    TH1F* histPt_2ndleadinglep_3e =  (0);
    TH1F* histPt_3dleadinglep_3e =  (0);
    TH1F* histPt_leadinglep_noTrig_3e =  (0);
    TH1F* histPt_2ndleadinglep_noTrig_3e =  (0);
    TH1F* histPt_3dleadinglep_noTrig_3e  =(0);
    
    TH1F* histPt_2e1mu =0;
    TH1F* histPt_noTrig_2e1mu =0;
    TH1F* histPt_leadinglep_2e1mu =(0);
    TH1F* histPt_2ndleadinglep_2e1mu =(0);
    TH1F* histPt_3dleadinglep_2e1mu =(0);
    TH1F* histPt_leadinglep_noTrig_2e1mu =(0);
    TH1F* histPt_2ndleadinglep_noTrig_2e1mu =(0);
    TH1F* histPt_3dleadinglep_noTrig_2e1mu =(0);
    
    TH1F* histPt_1e2mu = 0;
    TH1F* histPt_noTrig_1e2mu = 0;
    TH1F* histPt_leadinglep_1e2mu = (0);
    TH1F* histPt_2ndleadinglep_1e2mu = (0);
    TH1F* histPt_3dleadinglep_1e2mu = (0);
    TH1F* histPt_leadinglep_noTrig_1e2mu = (0);
    TH1F* histPt_2ndleadinglep_noTrig_1e2mu = (0);
    TH1F* histPt_3dleadinglep_noTrig_1e2mu = (0);
    
    if((isData || dataSetName.find("WZ")!=std::string::npos ) && checktrigger){
      histPt_all = new TH1F(("histPt_all"+ dataSetName).c_str(), ("histPt_all" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_all = new TH1F(("histPt_leadinglep_all"+ dataSetName).c_str(), ("histPt_leadinglep_all" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_all = new TH1F(("histPt_2ndleadinglep_all"+ dataSetName).c_str(), ("histPt_2ndleadinglep_all" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_all = new TH1F(("histPt_3dleadinglep_all"+ dataSetName).c_str(), ("histPt_3dleadinglep_all" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      histPt_noTrig_all = new TH1F(("histPt_noTrig_all"+ dataSetName).c_str(), ("histPt_noTrig_all" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_noTrig_all = new TH1F(("histPt_leadinglep_noTrig_all"+ dataSetName).c_str(), ("histPt_leadinglep_noTrig_all" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_noTrig_all = new TH1F(("histPt_2ndleadinglep_noTrig_all"+ dataSetName).c_str(), ("histPt_2ndleadinglep_noTrig_all" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_noTrig_all = new TH1F(("histPt_3dleadinglep_noTrig_all"+ dataSetName).c_str(), ("histPt_3dleadinglep_noTrig_all" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      
      histPt_3mu = new TH1F(("histPt_3mu"+ dataSetName).c_str(), ("histPt_3mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_3mu = new TH1F(("histPt_leadinglep_3mu"+ dataSetName).c_str(), ("histPt_leadinglep_3mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_3mu = new TH1F(("histPt_2ndleadinglep_3mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_3mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_3mu = new TH1F(("histPt_3dleadinglep_3mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_3mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      histPt_noTrig_3mu = new TH1F(("histPt_noTrig_3mu"+ dataSetName).c_str(), ("histPt_noTrig_3mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_noTrig_3mu = new TH1F(("histPt_leadinglep_noTrig_3mu"+ dataSetName).c_str(), ("histPt_leadinglep_noTrig_3mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_noTrig_3mu = new TH1F(("histPt_2ndleadinglep_noTrig_3mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_noTrig_3mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_noTrig_3mu = new TH1F(("histPt_3dleadinglep_noTrig_3mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_noTrig_3mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      
      
      histPt_3e = new TH1F(("histPt_3e"+ dataSetName).c_str(), ("histPt_3e" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_3e = new TH1F(("histPt_leadinglep_3e"+ dataSetName).c_str(), ("histPt_leadinglep_3e" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_3e = new TH1F(("histPt_2ndleadinglep_3e"+ dataSetName).c_str(), ("histPt_2ndleadinglep_3e" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_3e = new TH1F(("histPt_3dleadinglep_3e"+ dataSetName).c_str(), ("histPt_3dleadinglep_3e" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      histPt_noTrig_3e = new TH1F(("histPt_noTrig_3e"+ dataSetName).c_str(), ("histPt_noTrig_3e" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_noTrig_3e = new TH1F(("histPt_leadinglep_noTrig_3e"+ dataSetName).c_str(), ("histPt_leadinglep_noTrig_3e" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_noTrig_3e = new TH1F(("histPt_2ndleadinglep_noTrig_3e"+ dataSetName).c_str(), ("histPt_2ndleadinglep_noTrig_3e" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_noTrig_3e = new TH1F(("histPt_3dleadinglep_noTrig_3e"+ dataSetName).c_str(), ("histPt_3dleadinglep_noTrig_3e" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      
      
      histPt_2e1mu = new TH1F(("histPt_2e1mu"+ dataSetName).c_str(), ("histPt_2e1mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_2e1mu = new TH1F(("histPt_leadinglep_2e1mu"+ dataSetName).c_str(), ("histPt_leadinglep_2e1mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_2e1mu = new TH1F(("histPt_2ndleadinglep_2e1mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_2e1mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_2e1mu = new TH1F(("histPt_3dleadinglep_2e1mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_2e1mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      histPt_noTrig_2e1mu = new TH1F(("histPt_noTrig_2e1mu"+ dataSetName).c_str(), ("histPt_noTrig_2e1mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_noTrig_2e1mu = new TH1F(("histPt_leadinglep_noTrig_2e1mu"+ dataSetName).c_str(), ("histPt_leadinglep_noTrig_2e1mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_noTrig_2e1mu = new TH1F(("histPt_2ndleadinglep_noTrig_2e1mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_noTrig_2e1mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_noTrig_2e1mu = new TH1F(("histPt_3dleadinglep_noTrig_2e1mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_noTrig_2e1mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      
      
      histPt_1e2mu = new TH1F(("histPt_1e2mu"+ dataSetName).c_str(), ("histPt_1e2mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_1e2mu = new TH1F(("histPt_leadinglep_1e2mu"+ dataSetName).c_str(), ("histPt_leadinglep_1e2mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_1e2mu = new TH1F(("histPt_2ndleadinglep_1e2mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_1e2mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_1e2mu = new TH1F(("histPt_3dleadinglep_1e2mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_1e2mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      histPt_noTrig_1e2mu = new TH1F(("histPt_noTrig_1e2mu"+ dataSetName).c_str(), ("histPt_noTrig_1e2mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_leadinglep_noTrig_1e2mu = new TH1F(("histPt_leadinglep_noTrig_1e2mu"+ dataSetName).c_str(), ("histPt_leadinglep_noTrig_1e2mu" + dataSetName).c_str() , nbin_Pt_lep0, 0., endbin_Pt_lep0);
      histPt_2ndleadinglep_noTrig_1e2mu = new TH1F(("histPt_2ndleadinglep_noTrig_1e2mu"+ dataSetName).c_str(), ("histPt_2ndleadinglep_noTrig_1e2mu" + dataSetName).c_str() , nbin_Pt_lep1, 0., endbin_Pt_lep1);
      histPt_3dleadinglep_noTrig_1e2mu = new TH1F(("histPt_3dleadinglep_noTrig_1e2mu"+ dataSetName).c_str(), ("histPt_3dleadinglep_noTrig_1e2mu" + dataSetName).c_str() , nbin_Pt_lep2, 0., endbin_Pt_lep2);
      
      
    }
    
    
    bool PushBack = true;
    double bdiscrim  = -5.;
    double cbdiscrim = -5.;
    double cldiscrim = -5.;
    
    
    
    for (int ievt = 0; ievt < endEvent; ievt++)
    {
      ClearObjects(); // put everything to default values
      
      if(ievt == 0 ){ firstevent = true;}
      else{ firstevent = false;}
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
      
      // cout << "ievt "  << ievt << endl;
      /// Load event
      tTree[(dataSetName).c_str()]->GetEntry(ievt);
      if(!doDilep && !doCutflow && ((nMuons+nElectrons)!=3)){ continue;  }// making it faster
      if(removeBadMu && rejecteventBadPFmuon){ cout << "removing bad pf muon event " << endl; continue;}
      
      if(isData){
        if(evt_num < mineventnb)  mineventnb = evt_num;
        if(evt_num > maxeventnb)  maxeventnb = evt_num ;
        evtID = evt_num ;
      }
      
      
      if(applyMETfilter && !PassedMETFilter){   continue;}
      if(dataSetName.find("data_MET")!=std::string::npos && !PassedTriggerMET) continue;
      if(dataSetName.find("data_JetHt")!=std::string::npos && !PassedTriggerJET) continue;
      
      if(applytrigger && !PassedTrigger ) continue;
      if(applytriggerNoLogic && !PassedTriggerNoLogic) continue;
      if(applytriggerNoLogic2 && !PassedTriggerNoLogic2) continue;
      
      // fill objects,
      clock_t startFill = clock();
      
      tempPx = 0.;
      tempPy = 0.;
      tempHt = 0.;
      tempInvMassObj.SetPtEtaPhiE(0.,0., 0.,0.);
      
      
      
      muonID.clear();
      for(unsigned int iMu = 0; iMu < nMuons ; iMu++){
        
        if( fabs(eta_muon[iMu]) >= 2.4) {continue; }
        
        muon.Clear();
        
        muon.SetPtEtaPhiE(pt_muon[iMu], eta_muon[iMu], phi_muon[iMu], E_muon[iMu]);
        if(dorochester){
          double ptmu = pt_muon[iMu];
          bool databool = false;
          if( isfakes || isData) databool = true;
          ptmu = RochLeptonMatching(muon, mcParticles, databool, TrackLayers_muon[iMu],charge_muon[iMu], isNP, 0) * pt_muon[iMu] ; //*ptSF_muon[iMu];
          // cout << "rochester " <<  RochLeptonMatching(muon, mcParticles, isData, TrackLayers_muon[iMu],charge_muon[iMu], isNP, 0)  << " pt " << pt_muon[iMu] << endl;
          // muon.SetPtEtaPhiE(ptmu, eta_muon[iMu], phi_muon[iMu], E_muon[iMu]);
          muon.SetPtEtaPhiM(ptmu,eta_muon[iMu], phi_muon[iMu],0.105658);
        }
        if(muon.Pt() < 30. ){ continue; }
        selectedMuons.push_back(muon);
        selectedMuonsCharge.push_back(charge_muon[iMu]);
        selectedLeptons.push_back(muon);
        tempPx = tempPx + muon.Px();
        tempPy = tempPy + muon.Py();
        tempHt = tempHt + muon.Pt();
        tempInvMassObj = tempInvMassObj + muon;
        muonID.push_back(iMu);
      }
      // cout << "nMuons " << nMuons << " selected " << selectedMuons.size() << endl;
      electronID.clear();
      for(unsigned int iEl = 0; iEl < nElectrons ; iEl++){
        if(pt_electron[iEl]<35.){ continue;}
        if(fabs(eta_electron[iEl]) >= 2.1){ continue;}
        electron.Clear();
        electron.SetPtEtaPhiE(pt_electron[iEl], eta_electron[iEl], phi_electron[iEl], E_electron[iEl]);
        selectedElectrons.push_back(electron);
        selectedElectronsCharge.push_back(charge_electron[iEl]);
        selectedLeptons.push_back(electron);
        tempPx = tempPx + electron.Px();
        tempPy = tempPy + electron.Py();
        tempHt = tempHt + electron.Pt();
        tempInvMassObj = tempInvMassObj + electron;
        electronID.push_back(iEl);
      }
      
      
      
      if(!check_matching) sort(selectedLeptons.begin(), selectedLeptons.end(), HighestPt());
      MVA_TotalHt_lep = tempHt,
      MVA_TotalPt_lep = sqrt(tempPx*tempPx + tempPy*tempPx);
      MVA_TotalInvMass_lep = tempInvMassObj.M();
      
      
      tempHt_jet = 0.;
      tempPy_jet = 0.;
      tempPx_jet = 0.;
      PushBack = true;
      tempInvMassObj_jet.Clear();
      bdiscrim  = -5.;
      cbdiscrim = -5.;
      cldiscrim = -5.;
      
      for(unsigned int iJet = 0; iJet < nJets ; iJet++){
        if(pt_jet[iJet] < 30. ) continue;
        if(fabs(eta_jet[iJet]) >= 2.4){ continue;}
        jet.Clear();
        jet.SetPtEtaPhiE(pt_jet[iJet], eta_jet[iJet], phi_jet[iJet], E_jet[iJet]);
        bdiscrim = bdisc_jet[iJet];
        cbdiscrim = cdiscCvsB_jet[iJet];
        cldiscrim = cdiscCvsL_jet[iJet];
        
        if(applycharmsf && !isData){
          int binSF = charm_SFHisto_cvsb->GetXaxis()->FindBin(cbdiscrim);
          double charmSF = charm_SFHisto_cvsb->GetBinContent(binSF);
          cdiscCvsB_jet[iJet] = cdiscCvsB_jet[iJet]*charmSF;
          cbdiscrim = cdiscCvsB_jet[iJet];
          
          binSF = charm_SFHisto_cvsl->GetXaxis()->FindBin(cldiscrim);
          charmSF = charm_SFHisto_cvsl->GetBinContent(binSF);
          cdiscCvsL_jet[iJet] = cdiscCvsL_jet[iJet]*charmSF;
          cldiscrim = cdiscCvsL_jet[iJet];
          
        }
        
        
        
        PushBack = true;
        for(int iM = 0; iM < selectedMuons.size(); iM++){
          if(jet.DeltaR(selectedMuons[iM]) < 0.4) {
            PushBack = false;
            break;
          }
        }
        if(!PushBack) {   continue;}
        for(int iE = 0; iE < selectedElectrons.size(); iE++){
          if( jet.DeltaR(selectedElectrons[iE]) < 0.3) {
            PushBack = false;
            break;
          }
        }
        if(!PushBack){   continue;}
        else{
          //if(isfakes) cout << "pushing back " << iJet << endl;
          selectedJets.push_back(jet);
          selectedJetsID.push_back(iJet);
          tempPx = tempPx + jet.Px();
          tempPy = tempPy + jet.Py();
          tempHt = tempHt + jet.Pt();
          tempInvMassObj = tempInvMassObj + jet;
          
          tempPx_jet = tempPx_jet + jet.Px();
          tempPy_jet = tempPy_jet + jet.Py();
          tempHt_jet = tempHt_jet + jet.Pt();
          tempInvMassObj_jet = tempInvMassObj_jet + jet;
          
          // cout << "iJ " << iJet << " bdisc " << bdiscrim << " WP " << WPb_L << endl;
          if(bdiscrim >  WPb_L){
            selectedCSVLJetID.push_back(iJet);
          }
          else{
            selectednonCSVLJetID.push_back(iJet);
          }
          if(bdiscrim >  WPb_M) selectedCSVMJetID.push_back(iJet);
          else selectednonCSVMJetID.push_back(iJet);
          if(bdiscrim >  WPb_T) selectedCSVTJetID.push_back(iJet);
          else selectednonCSVTJetID.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Loose) selectedCvsBLJetID.push_back(iJet);
          else selectednonCvsBLJetID.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Medium) selectedCvsBMJetID.push_back(iJet);
          else selectednonCvsBMJetID.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Tight) selectedCvsBTJetID.push_back(iJet);
          else selectednonCvsBTJetID.push_back(iJet);
          if(cldiscrim > WPc_CvsL_Loose) selectedCvsLLJetID.push_back(iJet);
          else selectednonCvsLLJetID.push_back(iJet);
          if(cldiscrim > WPc_CvsL_Medium) selectedCvsLMJetID.push_back(iJet);
          else selectednonCvsLMJetID.push_back(iJet);
          if(cldiscrim > WPc_CvsL_Tight) selectedCvsLTJetID.push_back(iJet);
          else selectednonCvsLTJetID.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Loose && cldiscrim > WPc_CvsL_Loose) selectedCharmLJetsindex.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Medium && cldiscrim > WPc_CvsL_Medium) selectedCharmMJetsindex.push_back(iJet);
          if(cbdiscrim > WPc_CvsB_Tight && cldiscrim > WPc_CvsL_Tight) selectedCharmTJetsindex.push_back(iJet);
        }
        
      }
      
      // cout << "before selections " << endl;
      
      // selections
      //if(selectedJetsID.size()>6) continue; // temp fix
      if(selectedJetsID.size() == 0) continue;
      //if(mWT2 > 300.) continue;  // temp fix
      // if(met_Pt < 50. ) continue;
      
      /*
       for( int iJ = 0; iJ < selectedCSVLJetID.size(); iJ++){
       if( (bdisc_jet[selectedCSVLJetID[iJ]] <= WPb_L) && isfakes) cout << "is not loose b !! " << bdiscrim << " <= " << WPb_L << endl;
       }*/
      MVA_TotalHt = tempHt + met_Pt;
      MVA_TotalPt = sqrt(tempPx*tempPx + tempPy*tempPx);
      MVA_TotalInvMass = tempInvMassObj.M();
      
      MVA_TotalHt_jet = tempHt_jet + met_Pt;
      MVA_TotalPt_jet = sqrt(tempPx_jet*tempPx_jet + tempPy_jet*tempPx_jet);
      MVA_TotalInvMass_jet = tempInvMassObj_jet.M();
      
      
      
      // safety
      if(selectedMuons.size() == 3 && selectedElectrons.size() == 0){ channelInt = 0; }
      else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){ channelInt = 1;}
      else if(selectedElectrons.size() == 2 && selectedMuons.size() ==1 ){ channelInt = 2;  }
      else if(selectedElectrons.size() == 3 && selectedMuons.size() == 0){ channelInt = 3;  }
      else if(selectedMuons.size() > 1 && doDilep){ channelInt = 4; }
      else if(selectedElectrons.size() > 1 && doDilep){ channelInt = 5;  }
      else { continue;}
      
      
      
      
      // apply SF
      scaleFactor = 1.;
      scaleFactor_bfBT = 1.;
      scaleFactor_bfELSF = 1.;
      scaleFactor_bfMuSF = 1.;
      scaleFactor_bfPU = 1.;
      scaleFactor_puSF_down=1.;
      scaleFactor_puSF_up=1.;
      scaleFactor_electronSF_down=1.;
      scaleFactor_electronSF_up=1.;
      scaleFactor_muonSF_down=1.;
      scaleFactor_muonSF_up=1.;
      scaleFactor_btagSF_cferr1_down=1.;
      scaleFactor_btagSF_cferr1_up=1.;
      scaleFactor_btagSF_cferr2_down=1.;
      scaleFactor_btagSF_cferr2_up=1.;
      scaleFactor_btagSF_hf_down=1.;
      scaleFactor_btagSF_hf_up=1.;
      scaleFactor_btagSF_hfstats1_down=1.;
      scaleFactor_btagSF_hfstats1_up=1.;
      scaleFactor_btagSF_hfstats2_down=1.;
      scaleFactor_btagSF_hfstats2_up=1.;
      scaleFactor_btagSF_lf_down=1.;
      scaleFactor_btagSF_lf_up=1.;
      scaleFactor_btagSF_lfstats1_down=1.;
      scaleFactor_btagSF_lfstats1_up=1.;
      scaleFactor_btagSF_lfstats2_down=1.;
      scaleFactor_btagSF_lfstats2_up=1.;
      scaleFactor_puSF = 1. ;
      scaleFactor_btagSF = 1.;
      scaleFactor_muonSF = 1.;
      scaleFactor_electronSF = 1.;
      muonSFtemp = 1.;
      electronSFtemp = 1.;
      //puSF = 1.;
      
      
      if (! isData && !isfakes)
      {
        
        if(applyTrigSF){
          if(channelInt == 0) scaleFactor *= 1.;
          if(channelInt == 2) scaleFactor *= 1.0006003602;
          if(channelInt == 1) scaleFactor *= 1.0004001601;
          if(channelInt == 3) scaleFactor *= 0.9541174113;
        }
        if (applyMuonSF) {
          // if(ievt == 2)cout << "                - applying muon factors " << endl;
          muonSFtemp = 1.;
          for(unsigned int iMu = 0;  iMu < nMuons ;  iMu++){
            if(applyMuonSF_up){
              scaleFactor *= MuonIDSF_up[iMu] * MuonIsoSF_up[iMu] ;
              muonSFtemp *= MuonIDSF_up[iMu] * MuonIsoSF_up[iMu] ;
            }
            else if(applyMuonSF_down){
              scaleFactor *= MuonIDSF_down[iMu] * MuonIsoSF_down[iMu] ;
              muonSFtemp *= MuonIDSF_down[iMu] * MuonIsoSF_down[iMu] ;
            }
            else {
              // scaleFactor *= MuonIDSF[iMu] * MuonIsoSF[iMu]*ptSF_muon[iMu] ;
              scaleFactor *= MuonIDSF[iMu] * MuonIsoSF[iMu] * MuonTrackSF[iMu];
              //muonSFtemp *= MuonIDSF[iMu] * MuonIsoSF[iMu]*ptSF_muon[ ;
              muonSFtemp *= MuonIDSF[iMu] * MuonIsoSF[iMu]* MuonTrackSF[iMu];// TO FIX
            }
            scaleFactor_muonSF_up *= MuonIDSF_up[iMu] * MuonIsoSF_up[iMu];// * (MuonTrackSF[iMu]*1.01 );
            scaleFactor_muonSF_down *= MuonIDSF_down[iMu] * MuonIsoSF_down[iMu];// * (MuonTrackSF[iMu]*0.99 ) ;
            scaleFactor_muonSF *= MuonIDSF[iMu] * MuonIsoSF[iMu];
          }
          
        }
        else muonSFtemp = 1.;
        
        if (applyElectronSF) {
          //if(ievt == 2)cout << "                - applying electron factors " << endl;
          electronSFtemp = 1.;
          
          for(unsigned int iEl = 0; iEl < nElectrons ; iEl++){
            if(applyElectronSFup){
              scaleFactor *= ElectronSF_up[iEl] ;
              electronSFtemp *= ElectronSF_up[iEl];
            }
            else if(applyElectronSFdown){
              scaleFactor *= ElectronSF_down[iEl] ;
              electronSFtemp *= ElectronSF_down[iEl];
            }
            else{
              scaleFactor *= ElectronSF[iEl] ;
              electronSFtemp *= ElectronSF[iEl];
            }
            scaleFactor_electronSF_up *= ElectronSF_up[iEl] ;
            scaleFactor_electronSF_down *= ElectronSF_down[iEl] ;
            scaleFactor_electronSF  *= ElectronSF[iEl];
          }
          
        }
        else electronSFtemp = 1.;
        
        if (applyPUSF) {
          scaleFactor *= puSF;
          
          scaleFactor_puSF_up *= puSF_up;
          scaleFactor_puSF_down *= puSF_down;
          scaleFactor_puSF = puSF;
          
          
          
          
          //if(ievt == 2)cout << "                - applying pu factors " << endl;
        }
        else puSF  =1.;
        
        if (applyBTagSF) {
          scaleFactor *= btagSFshape;
          
          scaleFactor_btagSF_cferr1_down*= btagSFshape_down_cferr1  ;
          scaleFactor_btagSF_cferr1_up*= btagSFshape_up_cferr1  ;
          scaleFactor_btagSF_cferr2_down*= btagSFshape_down_cferr2  ;
          scaleFactor_btagSF_cferr2_up*= btagSFshape_up_cferr2  ;
          scaleFactor_btagSF_hf_down*= btagSFshape_down_hf  ;
          scaleFactor_btagSF_hf_up*= btagSFshape_up_hf  ;
          scaleFactor_btagSF_hfstats1_down*= btagSFshape_down_hfstats1  ;
          scaleFactor_btagSF_hfstats1_up*= btagSFshape_up_hfstats1  ;
          scaleFactor_btagSF_hfstats2_down*= btagSFshape_down_hfstats2  ;
          scaleFactor_btagSF_hfstats2_up*= btagSFshape_up_hfstats2  ;
          scaleFactor_btagSF_lf_down*= btagSFshape_down_lf  ;
          scaleFactor_btagSF_lf_up*= btagSFshape_up_lf  ;
          scaleFactor_btagSF_lfstats1_down*= btagSFshape_down_lfstats1  ;
          scaleFactor_btagSF_lfstats1_up*= btagSFshape_up_lfstats1  ;
          scaleFactor_btagSF_lfstats2_down*= btagSFshape_down_lfstats2  ;
          scaleFactor_btagSF_lfstats2_up*= btagSFshape_up_lfstats2  ;
          scaleFactor_btagSF = btagSFshape;
          
        }
        else btagSFshape = 1.;
        
        if (applyNloSF && isAMC) {
          scaleFactor *= nloWeight * nloSF;
        }  // additional SF due to number of events with neg weight!!
        
        /*
        scaleFactor_bfBT = scaleFactor/btagSFshape;
        scaleFactor_bfELSF = scaleFactor / electronSFtemp;
        scaleFactor_bfMuSF = scaleFactor / muonSFtemp;
        scaleFactor_bfPU = scaleFactor / puSF;
        
        //check for nan
        if(scaleFactor_bfBT != scaleFactor_bfBT) scaleFactor_bfBT = 0.;
        if(scaleFactor_bfMuSF != scaleFactor_bfMuSF) scaleFactor_bfMuSF = 0.;
        if(scaleFactor_bfELSF != scaleFactor_bfELSF) scaleFactor_bfELSF= 0.;
        if(scaleFactor_bfPU != scaleFactor_bfPU) scaleFactor_bfPU= 0.;
        
        
        scaleFactor_muonSF_down = ( scaleFactor_muonSF_down * scaleFactor ) / muonSFtemp;
        scaleFactor_muonSF_up = ( scaleFactor_muonSF_up * scaleFactor ) / muonSFtemp;
        scaleFactor_electronSF_down = ( scaleFactor_electronSF_down * scaleFactor) / electronSFtemp;
        scaleFactor_electronSF_up = ( scaleFactor_electronSF_up * scaleFactor) / electronSFtemp;
        scaleFactor_puSF_down = ( scaleFactor_puSF_down * scaleFactor) / puSF;
        scaleFactor_puSF_up = ( scaleFactor_puSF_up * scaleFactor) / puSF;
        scaleFactor_btagSF_cferr1_down= ( btagSFshape_down_cferr1  *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_cferr1_up= ( btagSFshape_up_cferr1 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_cferr2_down= ( btagSFshape_down_cferr2 *scaleFactor)/btagSFshape  ;
        scaleFactor_btagSF_cferr2_up= ( btagSFshape_up_cferr2 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_hf_down= ( btagSFshape_down_hf *scaleFactor)/btagSFshape  ;
        scaleFactor_btagSF_hf_up= ( btagSFshape_up_hf *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_hfstats1_down= ( btagSFshape_down_hfstats1*scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_hfstats1_up= ( btagSFshape_up_hfstats1 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_hfstats2_down= ( btagSFshape_down_hfstats2 *scaleFactor)/btagSFshape  ;
        scaleFactor_btagSF_hfstats2_up= ( btagSFshape_up_hfstats2 *scaleFactor)/btagSFshape  ;
        scaleFactor_btagSF_lf_down= ( btagSFshape_down_lf *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_lf_up= ( btagSFshape_up_lf *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_lfstats1_down= ( btagSFshape_down_lfstats1 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_lfstats1_up= ( btagSFshape_up_lfstats1 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_lfstats2_down= ( btagSFshape_down_lfstats2 *scaleFactor)/btagSFshape ;
        scaleFactor_btagSF_lfstats2_up= ( btagSFshape_up_lfstats2*scaleFactor)/btagSFshape ;*/
        
      }
      else if(isData || isfakes ){
        btagSFshape = 1.;
        scaleFactor = 1.;
        scaleFactor_bfBT = 1.;
        scaleFactor_bfELSF = 1.;
        scaleFactor_bfMuSF = 1.;
        scaleFactor_bfPU = 1.;
        scaleFactor_puSF_down=1.;
        scaleFactor_puSF_up=1.;
        scaleFactor_electronSF_down=1.;
        scaleFactor_electronSF_up=1.;
        scaleFactor_muonSF_down=1.;
        scaleFactor_muonSF_up=1.;
        scaleFactor_btagSF_cferr1_down=1.;
        scaleFactor_btagSF_cferr1_up=1.;
        scaleFactor_btagSF_cferr2_down=1.;
        scaleFactor_btagSF_cferr2_up=1.;
        scaleFactor_btagSF_hf_down=1.;
        scaleFactor_btagSF_hf_up=1.;
        scaleFactor_btagSF_hfstats1_down=1.;
        scaleFactor_btagSF_hfstats1_up=1.;
        scaleFactor_btagSF_hfstats2_down=1.;
        scaleFactor_btagSF_hfstats2_up=1.;
        scaleFactor_btagSF_lf_down=1.;
        scaleFactor_btagSF_lf_up=1.;
        scaleFactor_btagSF_lfstats1_down=1.;
        scaleFactor_btagSF_lfstats1_up=1.;
        scaleFactor_btagSF_lfstats2_down=1.;
        scaleFactor_btagSF_lfstats2_up=1.;
        muonSFtemp = 1.;
        electronSFtemp = 1.;
        puSF = 1.;
        
      }
      
      if(isfakes && ( channelInt == 0 || channelInt == 2)){ scaleFactor *= 0.0237522 * 0.0001 ; }
      if(isfakes && (channelInt== 1 || channelInt == 3)){ scaleFactor *= 0.20771 * 0.0001;  }
      
      
      
      if(doCutflow){
        MSPlot["cutflow"] ->Fill(0. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(0. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(0. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(0. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(0. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
      }
      if(MakeSelectionTable) {
        CutflowTable->Fill(d,0,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) CutflowTable_eee->Fill(d,0,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) CutflowTable_eeu->Fill(d,0,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) CutflowTable_uue->Fill(d,0,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) CutflowTable_uuu->Fill(d,0,scaleFactor*Luminosity/EquilumiSF);
      }
      bool threelepregion = false;
      bool twolepregion = false;
      if(selectedLeptons.size() == 3)  threelepregion = true;
      if(selectedElectrons.size() > 1 || selectedMuons.size() > 1) twolepregion = true;
      if(! threelepregion && ! twolepregion && !docharmsf ) continue;
      if(MakeSelectionTable) {
        CutflowTable->Fill(d,1,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) CutflowTable_eee->Fill(d,1,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) CutflowTable_eeu->Fill(d,1,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) CutflowTable_uue->Fill(d,1,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) CutflowTable_uuu->Fill(d,1,scaleFactor*Luminosity/EquilumiSF);
      }
      if(doCutflow){
        MSPlot["cutflow"] ->Fill(1. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(1. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(1. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(1. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(1. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
      }
      if(selectedElectrons.size() != nbOfLooseElectrons && !docharmsf && threelepregion ) continue;
      if(selectedMuons.size() != nbOfLooseMuons && !docharmsf && threelepregion ) continue;
      if(doCutflow){
        MSPlot["cutflow"] ->Fill(2. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(2. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(2. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(2. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(2. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
      }
      if(MakeSelectionTable) {
        CutflowTable->Fill(d,2,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) CutflowTable_eee->Fill(d,2,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) CutflowTable_eeu->Fill(d,2,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) CutflowTable_uue->Fill(d,2,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) CutflowTable_uuu->Fill(d,2,scaleFactor*Luminosity/EquilumiSF);
      }
      // cout << "in assigner" <<endl;
      
      //cout << "WmuIndiceF " << WmuIndiceF <<" WelecIndiceF "<< WelecIndiceF <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
      
      LeptonAssigner(selectedElectrons, selectedMuons,selectedElectronsCharge ,selectedMuonsCharge);
      //cout << "WmuIndiceF " << WmuIndiceF <<" WelecIndiceF "<< WelecIndiceF <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
      if(!Assigned ) continue;
      //cout << "in reco" << endl;
      
      ReconstructObjects(selectedJetsID, selectedMuons, selectedElectrons, selectedJets, Region, threelepregion);
      
      
      
      // cout << "twolepregion" << " " << twolepregion << " " << "threelepregion" << " " <<  threelepregion << endl;
      if (makePlots)
      {
        //cout << "ievt " << ievt << endl;
        //FillGeneralPlots(d, "control_afterAtLeast1Jet", decayChannels, isData, isfakes, threelepregion, twolepregion);
        //if(dataSetName.find("WZTo3LNu")!=std::string::npos) Fill1DPlots(dataSetName);
       
        
        //if(dataSetName.find("tZq")!=std::string::npos){ Fill1DPlots(dataSetName);}
        
      }
      //cout << "zmass" << endl;
      if(Zboson.M() < 76 || Zboson.M() > 106) continue;
      if(doCutflow){
        MSPlot["cutflow"] ->Fill(3. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(3. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(3. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(3. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(3. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
      }
      if(MakeSelectionTable) {
        CutflowTable->Fill(d,3,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) CutflowTable_eee->Fill(d,3,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) CutflowTable_eeu->Fill(d,3,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) CutflowTable_uue->Fill(d,3,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) CutflowTable_uuu->Fill(d,3,scaleFactor*Luminosity/EquilumiSF);
      }
      if( docharmsf && (selectedElectrons.size()+selectedMuons.size()) > 0){
        if(isData){
          for(int iterJet = 0; iterJet< selectedJetsID.size(); iterJet++){
          CvsB_Histo_data-> Fill(cdiscCvsB_jet[selectedJetsID[iterJet]], scaleFactor*Luminosity/EquilumiSF);
          CvsL_Histo_data-> Fill(cdiscCvsL_jet[selectedJetsID[iterJet]], scaleFactor*Luminosity/EquilumiSF);
          }
                }
        if(!isData && dataSetName.find("NP_overlay")==std::string::npos){
          for(int iterJet = 0; iterJet< selectedJetsID.size(); iterJet++){
            CvsB_Histo_sum-> Fill(cdiscCvsB_jet[selectedJetsID[iterJet]], scaleFactor*Luminosity/EquilumiSF);
            CvsL_Histo_sum-> Fill(cdiscCvsL_jet[selectedJetsID[iterJet]], scaleFactor*Luminosity/EquilumiSF);
          }
        }
      }
      
      
      
      if (makePlots )
      {
        FillGeneralPlots(d, "control_afterAtLeast1Jet_afterZWindow", decayChannels,isData, isfakes, threelepregion, twolepregion);
        
        
      }
      if((dataSetName.find("DY")!=std::string::npos || dataSetName.find("TT")!=std::string::npos || dataSetName.find("WWTo")!=std::string::npos|| dataSetName.find("Zjets")!=std::string::npos  || dataSetName.find("fake")!=std::string::npos || dataSetName.find("data")!=std::string::npos) && dofakevalidation && selectedJetsID.size() > 0 && selectedCSVLJetID.size()>0){
        FillFakeValidation(dataSetName,decayChannels,isData, isfakes, threelepregion, twolepregion);
      }
      if(selectednonCSVLJetID.size()>0 && makePlots ){
        // FillGeneralPlots(d, "control_afterAtLeast1Jet_afterZWindow_afterAtLeast1BJet", decayChannels,isData, isfakes, threelepregion, twolepregion);
        
      }
      
      
      // from here only 3lep analysis !!!!
      if(twolepregion && doDilep){ nSelectedEntriesDilep++; nSelectedEntriesDilepweighted += scaleFactor*Luminosity/EquilumiSF;}
      if((selectedMuons.size()+selectedElectrons.size())!= 3) continue;
      if(doCutflow){
        MSPlot["cutflow"] ->Fill(4. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(4. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(4. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(4. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(4. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
      }
      if(MakeSelectionTable) {
        CutflowTable->Fill(d,4,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 3) CutflowTable_eee->Fill(d,4,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 2) CutflowTable_eeu->Fill(d,4,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 1) CutflowTable_uue->Fill(d,4,scaleFactor*Luminosity/EquilumiSF);
        if(channelInt == 0) CutflowTable_uuu->Fill(d,4,scaleFactor*Luminosity/EquilumiSF);
      }
      if(!threelepregion) cout << "WARNING something went wrong with threelep region" << endl;
      
       if(threelepregion &&dataSetName.find("WZTo3LNu")!=std::string::npos && systematicplots ) Fill1DPlots(dataSetName, Luminosity/EquilumiSF, threelepregion,twolepregion); // FIX EVENTWEIGHT
      
      bool matcher = false;
      if(check_matching) matcher = MatchingFunction(dataSetName, selectedLeptons, selectedMuons, selectedElectrons, selectedJets,makeMatchingPlots, debugmatching);
      if(matcher && debugmatching) cout << " done with matching " << endl;
      
      // Signal regions and background region
      bool selected = false;
      if(selectedJets.size() == 1 && selectedCSVLJetID.size() > 0 && threelepregion){
        Region = 0;
        nSelectedEntriesST++;
        selected = true;
      /*  if(doCutflow){
          MSPlot["cutflow"] ->Fill(5. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(5. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(5. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(5. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(5. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        }*/
        if(MakeSelectionTable) {
          CutflowTable->Fill(d,5,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) CutflowTable_eee->Fill(d,5,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) CutflowTable_eeu->Fill(d,5,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) CutflowTable_uue->Fill(d,5,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) CutflowTable_uuu->Fill(d,5,scaleFactor*Luminosity/EquilumiSF);
        }
      } // ST region
      if(selectedJets.size() > 1 && selectedCSVLJetID.size() > 0 && threelepregion){
        Region = 1;
        nSelectedEntriesTT++;
        selected = true;
       /* if(doCutflow){
          MSPlot["cutflow"] ->Fill(6. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(6. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(6. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(6. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(6. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        }*/
        if(MakeSelectionTable) {
          CutflowTable->Fill(d,6,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) CutflowTable_eee->Fill(d,6,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) CutflowTable_eeu->Fill(d,6,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) CutflowTable_uue->Fill(d,6,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) CutflowTable_uuu->Fill(d,6,scaleFactor*Luminosity/EquilumiSF);
        }
      } // ttbar region
      if(selectedJets.size() >0 && selectedCSVLJetID.size() == 0 && threelepregion){
        Region = 2;
        nSelectedEntriesWZ++;
        selected = true;
        if(doCutflow){
          MSPlot["cutflow"] ->Fill(7. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) MSPlot["cutflow_eee"] ->Fill(7. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) MSPlot["cutflow_eeu"] ->Fill(7. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) MSPlot["cutflow_uue"] ->Fill(7. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) MSPlot["cutflow_uuu"] ->Fill(7. , datasets[d], true,scaleFactor*Luminosity/EquilumiSF);
        }
        if(MakeSelectionTable) {
          CutflowTable->Fill(d,7,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 3) CutflowTable_eee->Fill(d,7,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 2) CutflowTable_eeu->Fill(d,7,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 1) CutflowTable_uue->Fill(d,7,scaleFactor*Luminosity/EquilumiSF);
          if(channelInt == 0) CutflowTable_uuu->Fill(d,7,scaleFactor*Luminosity/EquilumiSF);
        }
      }// WZ control region
      if(selectedJets.size() >1 && selectedCSVLJetID.size() > 1 && threelepregion){
        Region = 3;
        nSelectedEntriesTTZ++;
        selected = true;
      }// ttZ control region
      if(!selected){continue; }
      
      
      
      //if(Region == 2 && (isData || dataSetName.find("WZ")!=std::string::npos)){
      if( (isData || dataSetName.find("WZ")!=std::string::npos) && checktrigger ){
        
        for(int iLep = 0; iLep < selectedLeptons.size() ; iLep++) {
          histPt_noTrig_all->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
          if(PassedTrigger) histPt_all->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
        }
        
        histPt_leadinglep_noTrig_all->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
        histPt_2ndleadinglep_noTrig_all->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
        histPt_3dleadinglep_noTrig_all->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
        
        if(PassedTrigger){
          histPt_leadinglep_all->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_2ndleadinglep_all->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_3dleadinglep_all->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
          
        }
        if(channelInt == 0){
          for(int iLep = 0; iLep < selectedLeptons.size() ; iLep++) {
            histPt_noTrig_3mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
            if(PassedTrigger) histPt_3mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
          }
          
          histPt_leadinglep_noTrig_3mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_2ndleadinglep_noTrig_3mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_3dleadinglep_noTrig_3mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
          
          if(PassedTrigger){
            histPt_leadinglep_3mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_2ndleadinglep_3mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_3dleadinglep_3mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
            
          }
        }
        else if(channelInt == 1){
          for(int iLep = 0; iLep < selectedLeptons.size() ; iLep++) {
            histPt_noTrig_1e2mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
            if(PassedTrigger) histPt_1e2mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
          }
          
          histPt_leadinglep_noTrig_1e2mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_2ndleadinglep_noTrig_1e2mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_3dleadinglep_noTrig_1e2mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
          
          if(PassedTrigger){
            histPt_leadinglep_1e2mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_2ndleadinglep_1e2mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_3dleadinglep_1e2mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
            
          }
        }
        else if(channelInt == 2){
          for(int iLep = 0; iLep < selectedLeptons.size() ; iLep++) {
            histPt_noTrig_2e1mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
            if(PassedTrigger) histPt_2e1mu->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
          }
          
          histPt_leadinglep_noTrig_2e1mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_2ndleadinglep_noTrig_2e1mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_3dleadinglep_noTrig_2e1mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
          
          if(PassedTrigger){
            histPt_leadinglep_2e1mu->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_2ndleadinglep_2e1mu->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_3dleadinglep_2e1mu->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
            
          }
        }
        else if(channelInt == 3){
          for(int iLep = 0; iLep < selectedLeptons.size() ; iLep++) {
            histPt_noTrig_3e->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
            if(PassedTrigger) histPt_3e->Fill(selectedLeptons[iLep].Pt(), scaleFactor*Luminosity/EquilumiSF);
          }
          
          histPt_leadinglep_noTrig_3e->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_2ndleadinglep_noTrig_3e->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
          histPt_3dleadinglep_noTrig_3e->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
          
          if(PassedTrigger){
            histPt_leadinglep_3e->Fill(selectedLeptons[0].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_2ndleadinglep_3e->Fill(selectedLeptons[1].Pt(), scaleFactor*Luminosity/EquilumiSF);
            histPt_3dleadinglep_3e->Fill(selectedLeptons[2].Pt(), scaleFactor*Luminosity/EquilumiSF);
            
          }
        }
      }
      
      
      // if(isfakes || isData) cout << "equilumisf = " << EquilumiSF << " Lumi " << Luminosity << " SF " << scaleFactor << " scaleFactor*Luminosity/EquiLumi " << scaleFactor*Luminosity/datasets[d]->EquivalentLumi() << " equilumi " << datasets[d]->EquivalentLumi() << endl;
      // in MSPlot automatically there is divided by eqlumi, for MC this is 1. but for data this is equal to the lumi
      // for MC: equilumiSF is calculated to fix this factor 1.
      if(!isfakes &&  !isData){
        if(Region == 0 ) nSelectedEntriesSTweighted += scaleFactor*Luminosity/EquilumiSF; //
        if(Region == 1 ) nSelectedEntriesTTweighted += scaleFactor*Luminosity/EquilumiSF;
        if(Region == 2 ) nSelectedEntriesWZweighted += scaleFactor*Luminosity/EquilumiSF;
        if(Region == 3 ) nSelectedEntriesTTZweighted += scaleFactor*Luminosity/EquilumiSF;
      }
      else if(isfakes || isData) {
        if(Region == 0 ) nSelectedEntriesSTweighted += scaleFactor*Luminosity/datasets[d]->EquivalentLumi();
        if(Region == 1 ) nSelectedEntriesTTweighted += scaleFactor*Luminosity/datasets[d]->EquivalentLumi();
        if(Region == 2 ) nSelectedEntriesWZweighted += scaleFactor*Luminosity/datasets[d]->EquivalentLumi();
        if(Region == 3 ) nSelectedEntriesTTZweighted += scaleFactor*Luminosity/datasets[d]->EquivalentLumi();
      }
      if((isData || dataSetName.find("WZ")!=std::string::npos)  && checktrigger ){
        myfile << evt_num << endl;
        if(PassedTrigger)  myfiletrigged << evt_num << endl;
        
        if(Region ==0){
          myfileST << evt_num << endl;
          if(PassedTrigger)  myfileSTtrigged << evt_num << endl;
        }
        else if(Region == 1){
          myfileTT << evt_num << endl;
          if(PassedTrigger)  myfileTTtrigged << evt_num << endl;
        }
        else if(Region == 2){
          myfileWZ << evt_num << endl;
          if(PassedTrigger)  myfileWZtrigged << evt_num << endl;
        }
        
        
      }
      MVA_EqLumi = EquilumiSF;
      MVA_Luminosity = Luminosity;
      if(makeMVAtree ){
        //cout << "ievt " << ievt << endl;
        MakeMVAvars(Region, scaleFactor);
        
      }
      // cout << "region "<< Region << endl;
      /// Make plots
      if (makeMVAPlots )
      {
       // if(Region == 0) FillMVAPlots(d,dataSetName, Region, "singletop", decayChannels);
       // if(Region == 1) FillMVAPlots(d,dataSetName, Region, "toppair", decayChannels);
        if(Region == 2) FillMVAPlots(d,dataSetName, Region, "wzcontrol", decayChannels);
        // if(Region == 3) FillMVAPlots(d,dataSetName, Region, "ttzcontrol", decayChannels);
      }
      
      
      
    } // events
    
    if((isData || dataSetName.find("WZ")!=std::string::npos) && checktrigger){
      if(d == 0)triggerEfffile = TFile::Open( triggerEfffilename.c_str(), "RECREATE" );
      else triggerEfffile = TFile::Open( triggerEfffilename.c_str(), "UPDATE" );
      triggerEfffile->cd();
      histPt_noTrig_all->Write();
      histPt_all->Write();
      histPt_2ndleadinglep_all->Write();
      histPt_3dleadinglep_all->Write();
      histPt_leadinglep_all->Write();
      histPt_2ndleadinglep_noTrig_all->Write();
      histPt_3dleadinglep_noTrig_all->Write();
      histPt_leadinglep_noTrig_all->Write();
      
      histPt_noTrig_3mu->Write();
      histPt_3mu->Write();
      histPt_2ndleadinglep_3mu->Write();
      histPt_3dleadinglep_3mu->Write();
      histPt_leadinglep_3mu->Write();
      histPt_2ndleadinglep_noTrig_3mu->Write();
      histPt_3dleadinglep_noTrig_3mu->Write();
      histPt_leadinglep_noTrig_3mu->Write();
      
      histPt_noTrig_3e->Write();
      histPt_3e->Write();
      histPt_2ndleadinglep_3e->Write();
      histPt_3dleadinglep_3e->Write();
      histPt_leadinglep_3e->Write();
      histPt_2ndleadinglep_noTrig_3e->Write();
      histPt_3dleadinglep_noTrig_3e->Write();
      histPt_leadinglep_noTrig_3e->Write();
      
      histPt_noTrig_2e1mu->Write();
      histPt_2e1mu->Write();
      histPt_2ndleadinglep_2e1mu->Write();
      histPt_3dleadinglep_2e1mu->Write();
      histPt_leadinglep_2e1mu->Write();
      histPt_2ndleadinglep_noTrig_2e1mu->Write();
      histPt_3dleadinglep_noTrig_2e1mu->Write();
      histPt_leadinglep_noTrig_2e1mu->Write();
      
      histPt_noTrig_1e2mu->Write();
      histPt_1e2mu->Write();
      histPt_2ndleadinglep_1e2mu->Write();
      histPt_3dleadinglep_1e2mu->Write();
      histPt_leadinglep_1e2mu->Write();
      histPt_2ndleadinglep_noTrig_1e2mu->Write();
      histPt_3dleadinglep_noTrig_1e2mu->Write();
      histPt_leadinglep_noTrig_1e2mu->Write();
      
      triggerEfffile->Write();
      triggerEfffile->Close();
    }
    
    
    
    if((isData || dataSetName.find("WZ")!=std::string::npos)  && checktrigger){
      myfile.close();
      myfiletrigged.close();
      
      myfileST.close();
      myfileSTtrigged.close();
      
      myfileTT.close();
      myfileTTtrigged.close();
      
      myfileWZ.close();
      myfileWZtrigged.close();
    }
    if(makeMVAtree){
      firstevent = true;
      writeMVAtree();
    }
    /* if(isData || isfakes) {
     nSelectedEntriesSTweighted = nSelectedEntriesST;
     nSelectedEntriesTTweighted = nSelectedEntriesTT;
     nSelectedEntriesWZweighted = nSelectedEntriesWZ;
     }*/
    cout << "                nSelectedEntries ST region: " << nSelectedEntriesST << " weighted " << nSelectedEntriesSTweighted << endl;
    cout << "                nSelectedEntries TT region: " << nSelectedEntriesTT << " weighted " << nSelectedEntriesTTweighted << endl;
    cout << "                nSelectedEntries WZ region: " << nSelectedEntriesWZ  << " weighted " << nSelectedEntriesWZweighted << endl;
    cout << "                nSelectedEntries TTZ region: " << nSelectedEntriesTTZ  << " weighted " << nSelectedEntriesTTZweighted << endl;
    if(doDilep) cout << "                nSelectedEntries dilep region: " << nSelectedEntriesDilep  << " weighted " << nSelectedEntriesDilepweighted << endl;
    cout << endl;
    if(check_matching) MatchingEfficiency();
  } // data
  
  if(applycharmsf){
    delete charm_SFHisto_cvsb;
    delete charm_SFHisto_cvsl;
    delete charm_SFHisto_cvsb_up;
    delete charm_SFHisto_cvsl_up;
    delete charm_SFHisto_cvsb_down;
    delete charm_SFHisto_cvsl_down;
    charmscalefactorsfile->Close();
    delete charmscalefactorsfile;
    
  }
  
  
  
  if(docharmsf){
   
    
    
    charmscalefactorsfile = TFile::Open( charmscalefactorsfilename.c_str(), "RECREATE" );
    
    SumNormal_cvsb = (TH1F*) (CvsB_Histo_sum)->Clone("SumNormal_cvsb");
    SumNormal_cvsl = (TH1F*) (CvsL_Histo_sum)->Clone("SumNormal_cvsl");
    dataNormal_cvsb = (TH1F*) (CvsB_Histo_data)->Clone("dataNormal_cvsb");
    dataNormal_cvsl = (TH1F*) (CvsL_Histo_data)->Clone("dataNormal_cvsl");
    
    SumNormal_cvsb->Scale(1/SumNormal_cvsb->Integral());
    SumNormal_cvsl->Scale(1/SumNormal_cvsl->Integral());
    dataNormal_cvsb->Scale(1/SumNormal_cvsb->Integral());
    dataNormal_cvsl->Scale(1/SumNormal_cvsl->Integral());
    
    charm_SFHisto_cvsb = (TH1F*) dataNormal_cvsb->Clone("charm_SFHisto_cvsb");
    charm_SFHisto_cvsl = (TH1F*) dataNormal_cvsl->Clone("charm_SFHisto_cvsl");
    
    charm_SFHisto_cvsb->Divide((TH1F*)SumNormal_cvsb);
    charm_SFHisto_cvsl->Divide((TH1F*)SumNormal_cvsl);
    
    charmscalefactorsfile->cd();
    
    
    CvsB_Histo_data->Write();
    CvsB_Histo_sum->Write();
    CvsL_Histo_data->Write();
    CvsL_Histo_sum->Write();
    SumNormal_cvsb->Write();
    SumNormal_cvsl->Write();
    dataNormal_cvsb->Write();
    dataNormal_cvsl->Write();
    charm_SFHisto_cvsl->Write();
    charm_SFHisto_cvsb->Write();
    
    
    charm_SFHisto_cvsb_up = (TH1F*) charm_SFHisto_cvsb->Clone("charm_SFHisto_cvsb_up");
    charm_SFHisto_cvsl_up = (TH1F*) charm_SFHisto_cvsl->Clone("charm_SFHisto_cvsl_up");
    charm_SFHisto_cvsb_down = (TH1F*) charm_SFHisto_cvsb->Clone("charm_SFHisto_cvsb_down");
    charm_SFHisto_cvsl_down = (TH1F*) charm_SFHisto_cvsl->Clone("charm_SFHisto_cvsl_down");
    double binerror;
    double bincontent;
    for(int ibin = 1; ibin < charm_SFHisto_cvsb->GetNbinsX(); ibin++){
      binerror = charm_SFHisto_cvsb->GetBinError(ibin);
      bincontent = charm_SFHisto_cvsb->GetBinContent(ibin);
      charm_SFHisto_cvsb_up->SetBinContent(ibin,bincontent + binerror);
      charm_SFHisto_cvsb_down->SetBinContent(ibin,bincontent - binerror);
     
      binerror = charm_SFHisto_cvsl->GetBinError(ibin);
      bincontent = charm_SFHisto_cvsl->GetBinContent(ibin);
      charm_SFHisto_cvsl_up->SetBinContent(ibin,bincontent + binerror);
      charm_SFHisto_cvsl_down->SetBinContent(ibin,bincontent - binerror);
    }
    
    charm_SFHisto_cvsl_up->Write();
    charm_SFHisto_cvsb_up->Write();
    charm_SFHisto_cvsl_down->Write();
    charm_SFHisto_cvsb_down->Write();
    
    
    
    
    charmscalefactorsfile->Write();
    charmscalefactorsfile->Close();
    delete CvsB_Histo_sum;
    delete CvsB_Histo_data;
    delete CvsL_Histo_data;
    delete CvsL_Histo_sum;
    
    delete charmscalefactorsfile;
  }
  

  
  
  if(checktrigger){
    triggerEfffile = TFile::Open( triggerEfffilename.c_str(), "UPDATE" );
    triggerEfffile->cd();
    
    
    gStyle->SetOptStat(0);
    
    TCanvas* tempC = new TCanvas();
    vector <string> plotter = {"histPt_leadinglep", "histPt_2ndleadinglep", "histPt_3dleadinglep", "histPt"};
    vector <string> plotternames = {"leading lepton p_{T} ", "2nd leading lepton p_{T} ", "3d leading lepton p_{T} ", "lepton p_{T} "};
    vector <string> plotterchanels  ={"all", "3mu", "3e", "1e2mu", "2e1mu"};
    vector <string> plotterchan = {"all", "3#mu","3e", "1e2#mu", "2e1#mu"};
    vector <double> xvalues;
    Double_t xl1=0.5, yl1=.3, xl2=xl1+.4, yl2=yl1+.2; //.7
    
    string tempstr = "1e2#mu";
    string tempstri = "1e2mu";
    
    
    for(int ich = 0; ich < plotterchanels.size(); ich++){
      tempstr = plotterchan[ich];
      tempstri = plotterchanels[ich];
      cout << "looking at " << tempstri << " channel" << endl;
      for(int iplot = 0; iplot < plotternames.size(); iplot++){
        TLegend *legtrig = new TLegend(xl1,yl1,xl2,yl2);
        TGraphAsymmErrors* effPt_data_Graph = new TGraphAsymmErrors();
        TGraphAsymmErrors* effPt_MC_Graph = new TGraphAsymmErrors();
        
        effPt_data_Graph->Divide(((TH1F*)triggerEfffile->Get((plotter[iplot]+"_"+tempstri+"data_MET").c_str())),((TH1F*)triggerEfffile->Get((plotter[iplot]+"_noTrig_"+tempstri+"data_MET").c_str())),"cl=0.683 b(1,1) mode ");
        effPt_MC_Graph->Divide(((TH1F*)triggerEfffile->Get((plotter[iplot]+"_"+tempstri+"WZTo3LNu_amc_80X").c_str())),((TH1F*)triggerEfffile->Get((plotter[iplot]+"_noTrig_"+tempstri+"WZTo3LNu_amc_80X").c_str())),"cl=0.683 b(1,1) mode " );
        effPt_data_Graph->SetName("effPt_data_Graph");
        effPt_data_Graph->SetMarkerColor(kBlue);
        effPt_data_Graph->SetLineColor(kBlue);
        effPt_data_Graph->SetLineWidth(2);
        effPt_data_Graph->SetMarkerStyle(2);
        effPt_data_Graph->Write();
        effPt_MC_Graph->SetName("effPt_MC_Graph");
        effPt_MC_Graph->SetMarkerColor(kRed);
        effPt_MC_Graph->SetLineColor(kRed);
        effPt_MC_Graph->SetLineWidth(2);
        effPt_MC_Graph->SetMarkerStyle(2);
        effPt_MC_Graph->Write();
        effPt_data_Graph->SetMaximum(1.5);
        effPt_MC_Graph->SetMaximum(1.5);
        effPt_data_Graph->SetMinimum(0);
        effPt_MC_Graph->SetMinimum(0);
        effPt_MC_Graph->GetXaxis()->SetTitle(plotternames[iplot].c_str());
        effPt_MC_Graph->GetYaxis()->SetTitle("#epsilon");
        effPt_data_Graph->GetXaxis()->SetTitle(plotternames[iplot].c_str());
        effPt_data_Graph->GetYaxis()->SetTitle("#epsilon");
        effPt_MC_Graph->SetTitle(("Trigger efficiency: "+tempstr+" channel").c_str());
        effPt_data_Graph->SetTitle(("Trigger efficiency: "+tempstr+" channel").c_str());
        legtrig->AddEntry(effPt_data_Graph,"Data","AP");
        legtrig->AddEntry(effPt_MC_Graph,"WZ+jets","AP");
        tempC = new TCanvas("triggereff", "triggereff");
        tempC->cd();
        effPt_data_Graph->Draw("AP");
        effPt_MC_Graph->Draw("P,sames");
        legtrig->Draw();
        tempC->Update();
        tempC->SetLogy();
        tempC->Update();
        tempC->Write();
        tempC->SaveAs( ("triggeff_"+tempstri+plotter[iplot]+".png").c_str() );
        
        xvalues.clear();
        for(int xval = 1 ; xval < ((TH1F*)triggerEfffile->Get((plotter[iplot]+"_"+tempstri+"data_MET").c_str()))->GetNbinsX()+1; xval ++){
          xvalues.push_back(((TH1F*)triggerEfffile->Get((plotter[iplot]+"_"+tempstri+"data_MET").c_str()))->GetBinCenter(xval));
          //cout << ((TH1F*)triggerEfffile->Get((plotter[iplot]+"data_MET").c_str()))->GetBinCenter(xval) << endl;
        }
        double x[xvalues.size()], y[xvalues.size()], xl[xvalues.size()], xh[xvalues.size()],yl[xvalues.size()],yh[xvalues.size()];
        for(int ix = 0; ix < xvalues.size(); ix++){
          y[ix] = effPt_data_Graph->Eval(xvalues[ix])/effPt_MC_Graph->Eval(xvalues[ix]);
          x[ix] = xvalues[ix];
          xh[ix] = TMath::Abs(effPt_MC_Graph->GetErrorXhigh(ix));
          xl[ix] = TMath::Abs(effPt_MC_Graph->GetErrorXlow(ix));
          double iyup_MC = TMath::Abs(effPt_MC_Graph->GetErrorYhigh(ix));
          double iydown_MC = TMath::Abs(effPt_MC_Graph->GetErrorYlow(ix));
          double iyup_data = TMath::Abs(effPt_data_Graph->GetErrorYhigh(ix));
          double iydown_data = TMath::Abs(effPt_data_Graph->GetErrorYlow(ix));
          
          double ymc = effPt_MC_Graph->Eval(xvalues[ix]);
          double ydata = effPt_data_Graph->Eval(xvalues[ix]);
          double termoneup = (iyup_data*iyup_data)/(ymc*ymc);
          double termtwoup = (iyup_MC*iyup_MC*ydata*ydata)/(ymc*ymc*ymc*ymc);
          double termonedown = (iydown_data*iydown_data)/(ymc*ymc);
          double termtwodown = (iydown_MC*iydown_MC*ydata*ydata)/(ymc*ymc*ymc*ymc);
          yh[ix] = sqrt(termoneup +   termtwoup);
          yl[ix] = sqrt(termonedown +   termtwodown);
          //cout << "xval: " << x[ix] << " + " << xh[ix] << " - " << xl[ix] << endl;
          //cout << "yval: " << y[ix] << " + " << yh[ix] << " - " << yl[ix] <<  endl;
        }
        
        TGraphAsymmErrors* scalefactors_Graph = new TGraphAsymmErrors(xvalues.size(),x,y,xl,xh,yl,yl);
        tempC = new TCanvas("triggereff", "triggereff");
        tempC->cd();
        scalefactors_Graph->SetTitle(("Trigger ScaleFactors: " + tempstr + " channel").c_str());
        scalefactors_Graph->GetXaxis()->SetTitle(plotternames[iplot].c_str());
        scalefactors_Graph->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
        scalefactors_Graph->SetMarkerColor(kRed);
        scalefactors_Graph->SetLineColor(kRed);
        scalefactors_Graph->SetLineWidth(2);
        scalefactors_Graph->SetMarkerStyle(2);
        scalefactors_Graph->SetMaximum(1.5);
        scalefactors_Graph->SetMinimum(0);
        scalefactors_Graph->Draw("AP");
        tempC->Update();
        tempC->Write();
        tempC->SaveAs( ("SF_trigger_"+tempstri+plotter[iplot]+".png").c_str() );
      }
      // cout << "trigger SF " << histPt_MET_SF->GetBinContent(1) << endl;
    }
    triggerEfffile->Write();
    triggerEfffile->Close();
    // delete triggerEfffile;
    
  }
  
  if(MakeSelectionTable){
    CutflowTable->TableCalculator(false, false, false, false, false, false, false, true); // mergeTT, mergeQCD , mergeW, mergeZ, mergeST, mergeVV, mergeTTV, NPmass -> write NP poverlay entries, nonpromptmc, tth, np zut, np zct
    string CutflowTablestring = "Cutflowtable_table.tex";
    CutflowTable->Write(CutflowTablestring.c_str(), true,true,true,true,true,true,true);  //output, witherr, merged, usebooktabs, addrawnbrs, addeff, add total aff, writelandscape
    
    CutflowTable_eee->TableCalculator(false, false, false, false, false, false, false, true);
    CutflowTablestring = "Cutflowtable_table_eee.tex";
    CutflowTable_eee->Write(CutflowTablestring.c_str(), true,true,true,true,true,true,true);
    
    CutflowTable_eeu->TableCalculator(false, false, false, false, false, false, false, true);
    CutflowTablestring = "Cutflowtable_table_eeu.tex";
    CutflowTable_eeu->Write(CutflowTablestring.c_str(), true,true,true,true,true,true,true);
    
    CutflowTable_uue->TableCalculator(false, false, false, false, false, false, false, true);
    CutflowTablestring = "Cutflowtable_table_uue.tex";
    CutflowTable_uue->Write(CutflowTablestring.c_str(), true,true,true,true,true,true,true);
    
    CutflowTable_uuu->TableCalculator(false, false, false, false, false, false, false, true);
    CutflowTablestring = "Cutflowtable_table_uuu.tex";
    CutflowTable_uuu->Write(CutflowTablestring.c_str(), true,true,true,true,true,true,true);
  }
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  if(makePlots || dofakevalidation || makeMatchingPlots || systematicplots){
    string rootFileName ="NtuplePlots.root";
    string place =pathOutputdate+"/MSPlot/";
    string placeTH1F = pathOutputdate+"/TH1F/";
    string placeTH2F = pathOutputdate+"/TH2F/";
    vector <string> vlabel_chan = {"3#mu", "1e2#mu", "2e1#mu", "3e"};
    if(doDilep) vlabel_chan = {"3#mu", "1e2#mu", "2e1#mu", "3e","2#mu","2e"};
    //vector <string> vlabel_cuts =  {"",">1l,>0j", "SF pair","lep veto","Z mass",">2l","STSR","TTSR","WZCR"};
    mkdir(place.c_str(),0777);
    mkdir(placeTH1F.c_str(),0777);
    mkdir(placeTH2F.c_str(),0777);
    
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutputdate+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    if(makePlots && !makeMatchingPlots){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
        cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        
        //temp->showNumberEntries(showEntriesLegend);
       // temp->setPreliminary(false);
        //temp->setErrorBandFile(errorbandfile, dosystfile);
        //temp->Draw(name,RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal);
        string name = it->first;
        if(!datafound) temp->setDataLumi(Luminosity);
        if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
        if(name.find("eee")!=std::string::npos) temp->setChannel(true, "3e");
        if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "2e1#mu");
        if(name.find("uue")!=std::string::npos) temp->setChannel(true, "1e2#mu");
        if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "3#mu");
        if(name.find("uu")!=std::string::npos && name.find("uuu")==std::string::npos  && name.find("uue")==std::string::npos && doDilep) temp->setChannel(true, "2#mu");
        if(name.find("ee")!=std::string::npos && name.find("eee")==std::string::npos  && name.find("eeu")==std::string::npos && doDilep) temp->setChannel(true, "2e");
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        if(name.find("cutflow")!=std::string::npos) temp->setBins(v_cutflow);
        temp->Draw(name, 1, false, false, false, 10);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        cout << "writing to " << pathOutputdate+"MSPlot" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        temp->Write(fout, name, true, (pathOutputdate+"/MSPlot").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
    }
    if(makeMatchingPlots || makePlots){
      TDirectory* th1dir = fout->mkdir("1D_histograms");
      th1dir->cd();
      gStyle->SetOptStat(1110);
      for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
      {
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
        TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
        tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
      }
    }
    if(dofakevalidation){
      TDirectory* th1dirfakes = fout->mkdir("1D_histograms_fakevalidation");
      th1dirfakes->cd();
      gStyle->SetOptStat(0);
      
      string poststring = "2lep_jets_bjets_zmasswindow";
      string posttitle = " channel - >= 2 lep, >0 j, >0 b, Zmass window";
      vector<string> plotnames = {"ZbosonPt_","ZbosonEta_","ZbosonPhi_","WlepPt_","WlepEta_","WlepPhi_","TrMassW_"};
      vector<string> plottitles = {"p_{T} Z boson (GeV)","#eta Z boson","#phi Z boson","p_{T} l_{W} (GeV)","#eta l_{W}","#phi l_{W}","transv. mass W boson (GeV)"};
      
      /*vector<string> plotnames_u = {"ZbosonPtMu_","ZbosonEtaMu_","ZbosonPhiMu_","WlepPtMu_","WlepEtaMu_","WlepPhiMu_","TrMassWMu_"};
       vector<string> plottitles_u = {"p_{T} Z_{#mu#mu} boson (GeV)","#eta Z_{#mu#mu} boson","#phi Z_{#mu#mu} boson","p_{T} #mu_{W} (GeV)","#eta #mu_{W}","#phi #mu_{W}","transv. mass W_{#mu} boson (GeV)"};
       
       vector<string> plotnames_e = {"ZbosonPtEl_","ZbosonEtaEl_","ZbosonPhiEl_","WlepPtEl_","WlepEtaEl_","WlepPhiEl_","TrMassWEl_"};
       vector<string> plottitles_e = {"p_{T} Z_{ee} boson (GeV)","#eta Z_{ee} boson","#phi Z_{ee} boson","p_{T} e_{W} (GeV)","#eta e_{W}","#phi e_{W}","transv. mass W_{e} boson (GeV)"};
       
       vector<string> channames = {"all","eee","uuu","uue","eeu","ee","uu"};
       vector<string> chantitle = {"all","3e","3#mu","1e2#mu","2e1#mu","2e","2#mu"};
       
       vector<string> regionname = {"2lep","3lep"};
       
       for(int iterv = 0; iterv < 7; iterv++){
       plotnames.push_back(plotnames_u[iterv]);
       plottitles.push_back(plottitles_u[iterv]);
       plotnames.push_back(plotnames_e[iterv]);
       plottitles.push_back(plottitles_e[iterv]);
       }*/
      
      vector<string> channames = {"all","ee","uu"};
      vector<string> chantitle = {"all","2e","2#mu"};
      
      
      string channelst;
      string channeltitle;
      
      vector<string> regionname = {"2lep","3lep"};
      regionname = {"2lep"};
      for(int reg = 0; reg < regionname.size(); reg++){
        
        for(int icha = 0; icha < channames.size(); icha++){
          channelst = regionname[reg]+"_"+channames[icha];
          channeltitle = chantitle[icha];
          for(int iv = 0; iv < plotnames.size() ; iv++){
            
            TH1F* tempDY(0);
            TH1F* tempfake(0);
            TH1F* tempdata(0);
            for (std::map<std::string,TH1F*>::const_iterator it = histo1D_fakevvalidation.begin(); it != histo1D_fakevvalidation.end(); it++)
            {
              TH1F *temp = it->second;
              string name = it->first;
              if(name.find(plotnames[iv].c_str())==std::string::npos) continue;
              if(name.find(channelst.c_str())==std::string::npos) continue;
              
              if(name.find("data")!=std::string::npos){
                if(tempdata == 0) tempdata = (TH1F*) temp;
                else tempdata->Add((TH1F*) temp);
              }
              if(name.find("fake")!=std::string::npos){
                if(tempfake == 0) tempfake = (TH1F*) temp;
                else tempfake->Add((TH1F*) temp);
              }
              if(name.find("DY")!=std::string::npos || name.find("Zjets")!=std::string::npos || dataSetName.find("TT")!=std::string::npos || dataSetName.find("WW")!=std::string::npos ){
                if(tempDY == 0) tempDY = (TH1F*) temp;
                else tempDY->Add((TH1F*) temp);
              }
              
            }
            tempfake->Scale(1./tempfake->Integral());
            tempdata->Scale(1./tempdata->Integral());
            tempDY->Scale(1./tempDY->Integral());
            
            tempfake->SetLineColor(kRed);
            tempDY->SetLineColor(kBlue);
            tempdata->SetLineColor(kBlack);
            tempdata->SetMarkerSize(3);
            TCanvas* Canvasfakes = 0;
            Double_t xl1=0.5, yl1=.7, xl2=xl1+.4, yl2=yl1+.2;
            TLegend *legfakes = new TLegend(xl1,yl1,xl2,yl2);
            
            legfakes->AddEntry(tempdata,"Data","LPF");   // h1 and h2 are histogram pointers
            legfakes->AddEntry(tempfake,"DD non prompt","L");
            legfakes->AddEntry(tempDY,"MC non prompt","L");
            
            double maximum = TMath::Max(TMath::Max(tempfake->GetMaximum(), tempdata->GetMaximum()), tempDY->GetMaximum());
            tempdata->SetMaximum(maximum*1.5);
            
            tempdata->SetTitle((channeltitle+ posttitle).c_str());
            tempdata->GetXaxis()->SetTitle(plottitles[iv].c_str());
            
            Canvasfakes =  TCanvasCreator(tempdata, "tempdata" );//new TCanvas("Canvas_PU","Canvas_PU");
            Canvasfakes->cd();
            
            tempdata->Draw("ep");
            tempfake->Draw("same histo e");
            tempDY->Draw("same histo e");
            legfakes->Draw("");
            Canvasfakes->Update();
            Canvasfakes->SaveAs( ("fakevalidation_"+plotnames[iv]+poststring+channelst+".png").c_str());
            Canvasfakes->SetLogy();
            Canvasfakes->Update();
            Canvasfakes->SaveAs( ("fakevalidation_"+plotnames[iv]+"LogY_"+poststring+channelst+".png").c_str() );
            
            TH1F* ratioFakeVsFake = (TH1F*) tempDY->Clone("ratioFakeVsFake");
            ratioFakeVsFake->Divide((TH1F*) tempfake);
            
            
            double xmaxfakesVsfake = ratioFakeVsFake->GetXaxis()->GetBinLowEdge(ratioFakeVsFake->GetNbinsX()+1);
            double xminfakesVsfake = ratioFakeVsFake->GetXaxis()->GetBinLowEdge(1);
            TLine *line = new TLine(xminfakesVsfake,1, xmaxfakesVsfake,1);
            line->SetLineColor(kGray);
            
            ratioFakeVsFake->SetTitle((channeltitle+posttitle).c_str());
            ratioFakeVsFake->GetYaxis()->SetTitle("MC non prompt/ DY non prompt");
            ratioFakeVsFake->GetXaxis()->SetTitle(plottitles[iv].c_str());
            ratioFakeVsFake->SetMaximum(ratioFakeVsFake->GetMaximum()*1.5);
            Canvasfakes =  TCanvasCreator(ratioFakeVsFake, "ratioFakeVsFake" );//new TCanvas("Canvas_PU","Canvas_PU");
            Canvasfakes->cd();
            ratioFakeVsFake->Draw("histo e");
            line->Draw("same");
            Canvasfakes->Update();
            // Canvasfakes->SaveAs( ("fakevalidation_ratiofakes_"+plotnames[iv]+poststring+channelst+".png").c_str());
            Canvasfakes->SetLogy();
            Canvasfakes->Update();
            // Canvasfakes->SaveAs( ("fakevalidation_ratiofakes_"+plotnames[iv]+"LogY_"+poststring+channelst+".png").c_str());
            
            TH1F* ratioFakeVsData = (TH1F*) tempdata->Clone("ratioFakeVsData");
            ratioFakeVsData->Divide((TH1F*) tempfake);
            ratioFakeVsData->SetMaximum(ratioFakeVsData->GetMaximum()*1.5);
            ratioFakeVsData->GetXaxis()->SetTitle(plottitles[iv].c_str());
            ratioFakeVsData->SetLineColor(kBlue);
            ratioFakeVsData->SetTitle((channeltitle+posttitle).c_str());
            ratioFakeVsData->GetYaxis()->SetTitle("data/DD non prompt");
            Canvasfakes =  TCanvasCreator(ratioFakeVsData, "ratioFakeVsData" );//new TCanvas("Canvas_PU","Canvas_PU");
            Canvasfakes->cd();
            ratioFakeVsData->Draw("histo e");
            line->Draw("same");
            Canvasfakes->Update();
            // Canvasfakes->SaveAs( ("fakevalidation_ratiodata_"+plotnames[iv]+poststring+channelst+".png").c_str());
            Canvasfakes->SetLogy();
            Canvasfakes->Update();
            // Canvasfakes->SaveAs( ("fakevalidation_ratiodata_"+plotnames[iv]+"LogY_"+poststring+channelst+".png").c_str());
            
            legfakes = new TLegend(xl1,yl1,xl2,yl2);
            legfakes->AddEntry(ratioFakeVsData,"DD non prompt/Data","L");   // h1 and h2 are histogram pointers
            legfakes->AddEntry(tempfake,"DD non prompt / MC non prompt","L");
            
            maximum = TMath::Max(ratioFakeVsFake->GetMaximum(),ratioFakeVsData->GetMaximum());
            ratioFakeVsData->SetMaximum(maximum*1.5);
            
            Canvasfakes =  TCanvasCreator(ratioFakeVsData, "ratioFakeVsData" );//new TCanvas("Canvas_PU","Canvas_PU");
            Canvasfakes->cd();
            ratioFakeVsData->SetTitle((channeltitle+posttitle).c_str());
            ratioFakeVsData->GetXaxis()->SetTitle(plottitles[iv].c_str());
            ratioFakeVsData->GetYaxis()->SetTitle("");
            ratioFakeVsData->SetLineColor(kBlue);
            ratioFakeVsFake->SetLineColor(kRed);
            ratioFakeVsData->Draw("histo e");
            ratioFakeVsFake->Draw("same histo e");
            line->Draw("same");
            legfakes->Draw("");
            Canvasfakes->Update();
            Canvasfakes->SaveAs( ("fakevalidation_ratio_"+plotnames[iv]+poststring+channelst+".png").c_str());
            Canvasfakes->SetLogy();
            Canvasfakes->Update();
            Canvasfakes->SaveAs( ("fakevalidation_ratio_"+plotnames[iv]+"LogY_"+poststring+channelst+".png").c_str());
          }
        }
        
      }
    }
    if(systematicplots){
     TDirectory* th1dirsys = fout->mkdir("1D_histograms_sys");
     th1dirsys->cd();
     gStyle->SetOptStat(1110);
     gStyle->SetOptStat(0);
     TH1F *tempnom =0;
     TH1F *tempup = 0;
     TH1F *tempdown = 0;
     int Nnom = 0;
     int Nup = 0;
     int Ndown = 0;
     double max = 0.;
     double max1 = 0.;
     TCanvas* Canvas = 0;
     Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
     TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
     std::map<std::string,TH1F*>::const_iterator nom = histo1D_PUSystematics.begin(), up=histo1D_PUSystematics.begin(), down=histo1D_PUSystematics.begin();
     
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_PUSystematics.begin(); it != histo1D_PUSystematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     leg->AddEntry(tempnom,"Nominal","L");   // h1 and h2 are histogram pointers
     leg->AddEntry(tempup,"Up","L");
     leg->AddEntry(tempdown,"Down","L");
     
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     
     
     tempnom->SetTitle("Nb Of vertices: PU SF");
     Canvas =  TCanvasCreator(tempnom, "Nb Of vertices: PU SF" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("h e");
     tempup->Draw("SAME,h e");
     tempdown->Draw("SAME,h e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"PUSF_nvtx.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"PUSF_nvtx_LogY.png").c_str() );
     
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     // electron SF
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_ElSystematics.begin(); it != histo1D_ElSystematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     
     
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("Pt electrons: El SF");
     Canvas =  TCanvasCreator(tempnom, "Pt leading electron: El SF" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"ELSF_ptelectron.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"ELSF_ptelectron_LogY.png").c_str() );
     
     // muon SF
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_MuSystematics.begin(); it != histo1D_MuSystematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("Pt muons: Mu SF");
     Canvas =  TCanvasCreator(tempnom, "Pt leading muon: Mu SF" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"MUSF_ptmuon.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"MUSF_ptmuon_LogY.png").c_str() );
     
     // btag SF cferr1
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Bcferr1Systematics.begin(); it != histo1D_Bcferr1Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF cferr1");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF cferr1" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdiscferr1.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdiscferr1_LogY.png").c_str() );
     
     // btag SF cferr2
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Bcferr2Systematics.begin(); it != histo1D_Bcferr2Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF cferr2");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF cferr2" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdiscferr2.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdiscferr2_LogY.png").c_str() );
     
     // btag SF hfstats1
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Bhfstats1Systematics.begin(); it != histo1D_Bhfstats1Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF hfstats1");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF hfstats1" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdishfstats1.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdishfstats1_LogY.png").c_str() );
     
     // btag SF hfstats2
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Bhfstats2Systematics.begin(); it != histo1D_Bhfstats2Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF hfstats2");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF hfstats2" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdishfstats2.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdishfstats2_LogY.png").c_str() );
     
     
     // btag SF lfstats1
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Blfstats1Systematics.begin(); it != histo1D_Blfstats1Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF lfstats1");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF lfstats1" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdislfstats1.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdislfstats1_LogY.png").c_str() );
     
     // btag SF lfstats2
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_Blfstats2Systematics.begin(); it != histo1D_Blfstats2Systematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF lfstats2");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF lfstats2" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdislfstats2.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdislfstats2_LogY.png").c_str() );
     
     
     // btag SF hf
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_BhfSystematics.begin(); it != histo1D_BhfSystematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF hf");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF hf" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdishf.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdishf_LogY.png").c_str() );
     
     // btag SF lf
     tempnom =0;
     tempup = 0;
     tempdown = 0;
     
     for (std::map<std::string,TH1F*>::const_iterator it = histo1D_BlfSystematics.begin(); it != histo1D_BlfSystematics.end(); it++)
     {
     string name = it->first;
     if(name.find("nom")!=std::string::npos){
     nom = it;
     if(tempnom == 0) tempnom = nom->second;
     else tempnom->Add(nom->second);
     }
     if(name.find("up")!=std::string::npos){
     up= it;
     if(tempup == 0) tempup = up->second;
     else tempup->Add(up->second);
     }
     if(name.find("down")!=std::string::npos){
     down = it;
     if(tempdown == 0) tempdown = down->second;
     else tempdown->Add(down->second);
     }
     
     }
     Nnom= tempnom->GetNbinsX();
     tempnom->SetBinContent(Nnom,tempnom->GetBinContent(Nnom)+tempnom->GetBinContent(Nnom+1));
     tempnom->SetBinContent(Nnom+1,0);
     tempnom->SetEntries(tempnom->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempnom->Write();
     Nup = tempup->GetNbinsX();
     tempup->SetBinContent(Nup,tempup->GetBinContent(Nup)+tempup->GetBinContent(Nup+1));
     tempup->SetBinContent(Nup+1,0);
     tempup->SetEntries(tempup->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempup->Write();
     Ndown= tempdown->GetNbinsX();
     tempdown->SetBinContent(Ndown,tempdown->GetBinContent(Ndown)+tempdown->GetBinContent(Ndown+1));
     tempdown->SetBinContent(Ndown+1,0);
     tempdown->SetEntries(tempdown->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
     tempdown->Write();
     tempnom->SetLineColor(kRed);
     tempup->SetLineColor(kBlue);
     tempdown->SetLineColor(kGreen+3);
     
     max1 = TMath::Max(tempnom->GetMaximum(), tempup->GetMaximum());
     max = TMath::Max(max1, tempdown->GetMaximum());
     tempnom->SetMaximum(max*1.2);
     tempnom->SetTitle("CSVv2 leading jet: btag SF lf");
     Canvas =  TCanvasCreator(tempnom, "CSVv2: btag SF lf" );//new TCanvas("Canvas_PU","Canvas_PU");
     Canvas->cd();
     tempnom->Draw("histo e");
     tempup->Draw("same histo e");
     tempdown->Draw("same histo e");
     leg->Draw("sames");
     Canvas->SaveAs( (placeTH1F+"BSF_bdislf.png").c_str() );
     Canvas->SetLogy();
     Canvas->Update();
     Canvas->SaveAs( (placeTH1F+"BSF_bdislf_LogY.png").c_str() );
     
     delete tempdown;
     delete tempnom;
     delete tempup;
     delete leg;
    
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
      tempCanvas->SaveAs( (placeTH2F+it->first+".png").c_str() );
    }
    
    fout->Close();
    
    delete fout;
    
  }
  
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

///////////////////////////////////// MVA INPUT /////////////////////////////////////////
void MakeMVAvars(int Region, Double_t scaleFactor){
  clock_t start_sub = clock();
  
  MVA_x1 = x1;
  MVA_x2 = x2;
  MVA_q = q;
  MVA_id1 = id1;
  MVA_id2 = id2;
  MVA_weight0 = weight0;
  MVA_weight1 = weight1;
  MVA_weight2 = weight2;
  MVA_weight3 = weight3;
  MVA_weight4 = weight4;
  MVA_weight5 = weight5;
  MVA_weight6 = weight6;
  MVA_weight7 = weight7;
  MVA_weight8 = weight8;
  MVA_hdamp_down = hdamp_down;
  MVA_hdamp_up = hdamp_up;
  
  
  
  MVA_channel = channelInt;
  MVA_Luminosity = Luminosity;
  MVA_EqLumi = EquilumiSF;
  MVA_weight = static_cast<float>( scaleFactor * Luminosity /EquilumiSF );
  MVA_weight_nom =  scaleFactor ;
  MVA_weight_puSF = scaleFactor_puSF;
  MVA_weight_muonSF = scaleFactor_muonSF;
  MVA_weight_electronSF = scaleFactor_electronSF;
  MVA_weight_btagSF = scaleFactor_btagSF;
  MVA_weight_puSF_up =  scaleFactor_puSF_up ;
  MVA_weight_puSF_down =  scaleFactor_puSF_down ;
  MVA_weight_electronSF_up =  scaleFactor_electronSF_up ;
  MVA_weight_electronSF_down =  scaleFactor_electronSF_down ;
  MVA_weight_muonSF_up =  scaleFactor_muonSF_up ;
  MVA_weight_muonSF_down =  scaleFactor_muonSF_down ;
  MVA_weight_btagSF_cferr1_up =  scaleFactor_btagSF_cferr1_up ;
  MVA_weight_btagSF_cferr1_down =  scaleFactor_btagSF_cferr1_down ;
  MVA_weight_btagSF_cferr2_up =  scaleFactor_btagSF_cferr2_up ;
  MVA_weight_btagSF_cferr2_down =  scaleFactor_btagSF_cferr2_down ;
  MVA_weight_btagSF_hf_up =  scaleFactor_btagSF_hf_up ;
  MVA_weight_btagSF_hf_down =  scaleFactor_btagSF_hf_down ;
  MVA_weight_btagSF_hfstats1_up =  scaleFactor_btagSF_hfstats1_up ;
  MVA_weight_btagSF_hfstats1_down =  scaleFactor_btagSF_hfstats1_down ;
  MVA_weight_btagSF_hfstats2_up =  scaleFactor_btagSF_hfstats2_up ;
  MVA_weight_btagSF_hfstats2_down =  scaleFactor_btagSF_hfstats2_down ;
  MVA_weight_btagSF_lf_up =  scaleFactor_btagSF_lf_up ;
  MVA_weight_btagSF_lf_down =  scaleFactor_btagSF_lf_down ;
  MVA_weight_btagSF_lfstats1_up =  scaleFactor_btagSF_lfstats1_up ;
  MVA_weight_btagSF_lfstats1_down =  scaleFactor_btagSF_lfstats1_down ;
  MVA_weight_btagSF_lfstats2_up =  scaleFactor_btagSF_lfstats2_up ;
  MVA_weight_btagSF_lfstats2_down =  scaleFactor_btagSF_lfstats2_down ;
  
  MVA_region = static_cast<float>( Region);
  
  MVA_lepton0_pt = static_cast<float>( selectedLeptons[0].Pt());
  // cout << "reconstructing " << MVA_lepton0_pt << endl;
  MVA_lepton1_pt= static_cast<float>( selectedLeptons[1].Pt());
  MVA_lepton2_pt= static_cast<float>( selectedLeptons[2].Pt());
  
  MVA_lepton0_eta= static_cast<float>( selectedLeptons[0].Eta());
  MVA_lepton1_eta= static_cast<float>( selectedLeptons[1].Eta());
  MVA_lepton2_eta = static_cast<float>( selectedLeptons[2].Eta());
  MVA_lepton0_phi = static_cast<float>( selectedLeptons[0].Phi());
  MVA_lepton1_phi = static_cast<float>( selectedLeptons[1].Phi());
  MVA_lepton2_phi = static_cast<float>( selectedLeptons[2].Phi());
  
  MVA_nMuons = ( selectedMuons.size());
  MVA_nElectrons = ( selectedElectrons.size());
  // cout << "leptons done" << endl);
  if(selectedJets.size()>0){
    MVA_jet0_pt = static_cast<float>( selectedJets[0].Pt());
    MVA_jet0_eta = static_cast<float>( selectedJets[0].Eta());
    MVA_jet0_phi = static_cast<float>( selectedJets[0].Phi());
  }
  if(selectedJets.size()>1){
    MVA_jet1_pt = static_cast<float>( selectedJets[1].Pt());
    MVA_jet1_eta = static_cast<float>( selectedJets[1].Eta());
    MVA_jet1_phi = static_cast<float>( selectedJets[1].Phi());
  }
  
  MVA_nJets = static_cast<float>( selectedJets.size());
  MVA_NJets_CSVv2L =static_cast<float>( selectedCSVLJetID.size());
  MVA_NJets_CSVv2M= static_cast<float>( selectedCSVMJetID.size());
  MVA_NJets_CSVv2T= static_cast<float>( selectedCSVTJetID.size());
  //cout << "jets done" << endl);
  
  // SM side
  MVA_Wlep_pt = static_cast<float>( Wlep.Pt());
  MVA_Wlep_eta = static_cast<float>( Wlep.Eta());
  MVA_Wlep_phi = static_cast<float>( Wlep.Phi());
  MVA_SMbjet_pt = static_cast<float>( selectedJets[SMjetIndex].Pt());
  MVA_SMbjet_eta = static_cast<float>( selectedJets[SMjetIndex].Eta());
  MVA_SMbjet_phi = static_cast<float>( selectedJets[SMjetIndex].Phi());
  
  MVA_Wboson_pt = static_cast<float>( Wboson.Pt());
  MVA_Wboson_eta = static_cast<float>( Wboson.Eta());
  MVA_Wboson_phi = static_cast<float>( Wboson.Phi());
  MVA_met = static_cast<float>( met_Pt);
  MVA_SMtop_pt = static_cast<float>( SMtop.Pt());
  MVA_SMtop_eta = static_cast<float>( SMtop.Eta());
  MVA_SMtop_phi = static_cast<float>( SMtop.Phi());
  
  //cout << "SM side done" << endl);
  
  // FCNC side
  MVA_Zboson_pt = static_cast<float>( Zboson.Pt());
  MVA_Zboson_eta = static_cast<float>( Zboson.Eta());
  MVA_Zboson_phi = static_cast<float>( Zboson.Phi());
  
  if(selectedJets.size()>1 )
  {
    MVA_LightJet_pt = static_cast<float>( LightJet.Pt());
    MVA_LightJet_eta = static_cast<float>( LightJet.Eta());
    MVA_LightJet_phi = static_cast<float>( LightJet.Phi());
    
    MVA_FCNCtop_pt = static_cast<float>( FCNCtop.Pt());
    MVA_FCNCtop_eta = static_cast<float>( FCNCtop.Eta());
    MVA_FCNCtop_phi =static_cast<float>( FCNCtop.Phi());
    
    
  }
  
  MVA_nJets_CharmL = static_cast<float>( selectedCharmLJetsindex.size());
  MVA_nJets_CharmM = static_cast<float>( selectedCharmMJetsindex.size());
  MVA_nJets_CharmT = static_cast<float>( selectedCharmTJetsindex.size());
  
  // cout << "FCNC side " << endl);
  
  //SM kinematics
  
  MVA_SMtop_M = static_cast<float>(SMtop.M());
  MVA_mlb = static_cast<float>((SMbjet+Wlep).M());
  MVA_Wboson_M = static_cast<float>(Wboson.M());
  
  MVA_dRWlepb = static_cast<float>(Wlep.DeltaR(SMbjet));
  
  MVA_dPhiWlepb = static_cast<float>(Wlep.DeltaPhi(SMbjet));
  
  if(WmuIndiceF != -999) MVA_Wlep_Charge =static_cast<float> ( selectedMuonsCharge[WmuIndiceF]);
  else if(WelecIndiceF != -999) MVA_Wlep_Charge = static_cast<float>( selectedElectronsCharge[WelecIndiceF]);
  MVA_charge_asym = static_cast<float>( MVA_Wlep_Charge*fabs(Wlep.Eta()));
  if(selectedJets.size()>1) MVA_bdiscCSVv2_jet_1 = static_cast<float>( bdisc_jet[selectedJetsID[1]]);
  if(selectedJets.size()>0) MVA_bdiscCSVv2_jet_0 = static_cast<float>(bdisc_jet[selectedJetsID[0]]);
  MVA_CosTheta = static_cast<float>( (CosThetaCalculation(Wlep, metTLV, SMbjet, false)).first);
  MVA_CosTheta_alt = static_cast<float>( (CosThetaCalculation(Wlep, metTLV, SMbjet,false)).second);
  
  //cout << "sm kin" << endl);
  
  // FCNC kinematics
  if(selectedJets.size()>1 ) MVA_FCNCtop_M = static_cast<float>( FCNCtop.M());
  MVA_Zboson_M = static_cast<float>( Zboson.M());
  //cout << "Zboson mass " << Zboson_M << endl);
  
  if(selectedJets.size()>1 ) MVA_dRZc = static_cast<float>( Zboson.DeltaR(LightJet));
  if(selectedJets.size()>1 ) MVA_dPhiZc = static_cast<float>( Zboson.DeltaPhi(LightJet));
  
  if(selectedJets.size()>0) MVA_cdiscCvsB_jet_0 = static_cast<float>( cdiscCvsB_jet[selectedJetsID[0]]);
  if(selectedJets.size()>0) MVA_cdiscCvsL_jet_0 = static_cast<float>( cdiscCvsL_jet[selectedJetsID[0]]);
  if(selectedJets.size()>1) MVA_cdiscCvsB_jet_1 = static_cast<float>( cdiscCvsB_jet[selectedJetsID[1]]);
  if(selectedJets.size()>1)  MVA_cdiscCvsL_jet_1 = static_cast<float>( cdiscCvsL_jet[selectedJetsID[1]]);
  
  //cout << "fcnc kin " << endl);
  
  // interplay
  if(selectedJets.size()>1 ) MVA_dRSMFCNCtop = static_cast<float>( SMtop.DeltaR(FCNCtop));
  MVA_dRZb = static_cast<float>( Zboson.DeltaR(SMbjet));
  if(selectedJets.size()>1 ) MVA_dRWlepc = static_cast<float>( Wlep.DeltaR(LightJet));
  MVA_dRZWlep = static_cast<float>( Zboson.DeltaR(Wlep));
  MVA_dRZSMtop = static_cast<float>( Zboson.DeltaR(SMtop));
  MVA_dPhiZMET = static_cast<float>( Zboson.DeltaPhi(metTLV));
  
  if(selectedJets.size()>1 ) MVA_dPhiSMFCNCtop = static_cast<float>( SMtop.DeltaPhi(FCNCtop));
  MVA_dPhiZb = static_cast<float>( Zboson.DeltaPhi(SMbjet));
  if(selectedJets.size()>1 ) MVA_dPhiWlepc = static_cast<float>( Wlep.DeltaPhi(LightJet));
  MVA_dPhiZWlep = static_cast<float>( Zboson.DeltaPhi(Wlep));
  MVA_dPhiZSMtop = static_cast<float>( Zboson.DeltaPhi(SMtop));
  MVA_m3l = static_cast<float>( (selectedLeptons[0]+selectedLeptons[1]+selectedLeptons[2]).M());
  
  
  
  
  //cout << "interplay done " << endl);
  
  MVA_mWt2 = static_cast<float>(mWT2);
  MVA_mWt = static_cast<float>(mWT);
  
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
void createMVAtree(string dataSetName){
  clock_t start_sub = clock();
  tupfile->cd();
  //cout << "does tree exist? " << !tupfile->GetListOfKeys()->Contains("mvatree") <<endl;
  
  
  // event
  // for pdf unc
  mvatree->Branch("MVA_x1", &MVA_x1, "MVA_x1/D");
  mvatree->Branch("MVA_x2", &MVA_x2, "MVA_x2/D");
  mvatree->Branch("MVA_id1", &MVA_id1, "MVA_id1/I");
  mvatree->Branch("MVA_id2", &MVA_id2, "MVA_id2/I");
  mvatree->Branch("MVA_q", &MVA_q, "MVA_q/D");
  mvatree->Branch("MVA_weight0", &MVA_weight0," MVA_weight0/D");
  mvatree->Branch("MVA_weight1", &MVA_weight1," MVA_weight1/D");
  mvatree->Branch("MVA_weight2", &MVA_weight2," MVA_weight2/D");
  mvatree->Branch("MVA_weight3", &MVA_weight3," MVA_weight3/D");
  mvatree->Branch("MVA_weight4", &MVA_weight4," MVA_weight4/D");
  mvatree->Branch("MVA_weight5", &MVA_weight5," MVA_weight5/D");
  mvatree->Branch("MVA_weight6", &MVA_weight6," MVA_weight6/D");
  mvatree->Branch("MVA_weight7", &MVA_weight7," MVA_weight7/D");
  mvatree->Branch("MVA_weight8", &MVA_weight8," MVA_weight8/D");
  mvatree->Branch("MVA_hdamp_up",&MVA_hdamp_up,"MVA_hdamp_up/D");
  mvatree->Branch("MVA_hdamp_down",&MVA_hdamp_down,"MVA_hdamp_down/D");
  
  mvatree->Branch("MVA_channel", &MVA_channel , "MVA_channel/I");
  mvatree->Branch("MVA_weight", &MVA_weight, "MVA_weight/F");
  mvatree->Branch("MVA_weight_nom", &MVA_weight_nom, "MVA_weight_nom/D");
  mvatree->Branch("MVA_weight_puSF_up", &MVA_weight_puSF_up, "MVA_weight_puSF_up/D");
  mvatree->Branch("MVA_weight_puSF_down", &MVA_weight_puSF_down, "MVA_weight_puSF_down/D");
  mvatree->Branch("MVA_weight_electronSF_up", &MVA_weight_electronSF_up, "MVA_weight_electronSF_up/D");
  mvatree->Branch("MVA_weight_electronSF_down", &MVA_weight_electronSF_down, "MVA_weight_electronSF_down/D");
  mvatree->Branch("MVA_weight_puSF",&MVA_weight_puSF, "MVA_weight_puSF/D");
  mvatree->Branch("MVA_weight_muonSF",&MVA_weight_muonSF, "MVA_weight_muonSF/D");
  mvatree->Branch("MVA_weight_eletcronSF",&MVA_weight_electronSF, "MVA_weight_electronSF/D");
  mvatree->Branch("MVA_weight_btagSF",&MVA_weight_btagSF, "MVA_weight_btagSF/D");
  
  mvatree->Branch("MVA_weight_muonSF_up", &MVA_weight_muonSF_up, "MVA_weight_muonSF_up/D");
  mvatree->Branch("MVA_weight_muonSF_down", &MVA_weight_muonSF_down, "MVA_weight_muonSF_down/D");
  mvatree->Branch("MVA_weight_btagSF_cferr1_up", &MVA_weight_btagSF_cferr1_up, "MVA_weight_btagSF_cferr1_up/D");
  mvatree->Branch("MVA_weight_btagSF_cferr1_down", &MVA_weight_btagSF_cferr1_down, "MVA_weight_btagSF_cferr1_down/D");
  mvatree->Branch("MVA_weight_btagSF_cferr2_up", &MVA_weight_btagSF_cferr2_up, "MVA_weight_btagSF_cferr2_up/D");
  mvatree->Branch("MVA_weight_btagSF_cferr2_down", &MVA_weight_btagSF_cferr2_down, "MVA_weight_btagSF_cferr2_down/D");
  mvatree->Branch("MVA_weight_btagSF_hf_up", &MVA_weight_btagSF_hf_up, "MVA_weight_btagSF_hf_up/D");
  mvatree->Branch("MVA_weight_btagSF_hf_down", &MVA_weight_btagSF_hf_down, "MVA_weight_btagSF_hf_down/D");
  mvatree->Branch("MVA_weight_btagSF_hfstats1_up", &MVA_weight_btagSF_hfstats1_up, "MVA_weight_btagSF_hfstats1_up/D");
  mvatree->Branch("MVA_weight_btagSF_hfstats1_down", &MVA_weight_btagSF_hfstats1_down, "MVA_weight_btagSF_hfstats1_down/D");
  mvatree->Branch("MVA_weight_btagSF_hfstats2_up", &MVA_weight_btagSF_hfstats2_up, "MVA_weight_btagSF_hfstats2_up/D");
  mvatree->Branch("MVA_weight_btagSF_hfstats2_down", &MVA_weight_btagSF_hfstats2_down, "MVA_weight_btagSF_hfstats2_down/D");
  mvatree->Branch("MVA_weight_btagSF_lf_up", &MVA_weight_btagSF_lf_up, "MVA_weight_btagSF_lf_up/D");
  mvatree->Branch("MVA_weight_btagSF_lf_down", &MVA_weight_btagSF_lf_down, "MVA_weight_btagSF_lf_down/D");
  mvatree->Branch("MVA_weight_btagSF_lfstats1_up", &MVA_weight_btagSF_lfstats1_up, "MVA_weight_btagSF_lfstats1_up/D");
  mvatree->Branch("MVA_weight_btagSF_lfstats1_down", &MVA_weight_btagSF_lfstats1_down, "MVA_weight_btagSF_lfstats1_down/D");
  mvatree->Branch("MVA_weight_btagSF_lfstats2_up", &MVA_weight_btagSF_lfstats2_up, "MVA_weight_btagSF_lfstats2_up/D");
  mvatree->Branch("MVA_weight_btagSF_lfstats2_down", &MVA_weight_btagSF_lfstats2_down, "MVA_weight_btagSF_lfstats2_down/D");
  
  
  
  mvatree->Branch("MVA_region", &MVA_region, "MVA_region/F");
  mvatree->Branch("MVA_EqLumi",&MVA_EqLumi,"MVA_EqLumi/D");
  mvatree->Branch("MVA_Luminosity",&MVA_Luminosity,"MVA_Luminosity/D");
  
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
  mvatree->Branch("MVA_nElectrons", &MVA_nElectrons, "MVA_nElectrons/I");
  mvatree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
  mvatree->Branch("MVA_NJets_CSVv2L", &MVA_NJets_CSVv2L, "MVA_NJets_CSVv2L/F");
  mvatree->Branch("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, "MVA_NJets_CSVv2M/F");
  mvatree->Branch("MVA_NJets_CSVv2T", &MVA_NJets_CSVv2T, "MVA_NJets_CSVv2T/F");
  mvatree->Branch("MVA_nMuons", &MVA_nMuons, "MVA_nMuons/I");
  
  
  mvatree->Branch("MVA_met", &MVA_met, "MVA_met/F");
  
  
  //SM kinematics
  mvatree->Branch("MVA_mWt", &MVA_mWt,"MVA_mWt/F");
  mvatree->Branch("MVA_mWt2", &MVA_mWt2,"MVA_mWt2/F");
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
  mvatree->Branch("MVA_bdiscCSVv2_jet_1", &MVA_bdiscCSVv2_jet_1,"MVA_bdiscCSVv2_jet_1/F");
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
void writeMVAtree(){
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

///////////////////////////////////// EventReco /////////////////////////////////////////
void LeptonAssigner(vector<TLorentzVector> electrons, vector<TLorentzVector> muons, std::vector<int> electronsCharge,std::vector<int> muonsCharge){
  clock_t start_sub =  clock();
  
  //  cout << " in assigner " << endl;
  
  Assigned = false;
  if(electronsCharge.size() != electrons.size() || muons.size() != muonsCharge.size()) cout << "ERROR vectors not filled properly" << endl;
  
  
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
  
  Double_t mass01 = 9999.;
  Double_t mass02 = 9999.;
  Double_t mass12 = 9999.;
  
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
    else if(mass12 <= mass01 && mass12 <= mass02 && can12){
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
    else if(mass02 <= mass12 && mass02 <= mass01 && can02){
      Assigned = true;
      ZmuIndiceF_0 = 0; ZmuIndiceF_1=2;WmuIndiceF = 1;
    }
    else if(mass12 <= mass01 && mass12 <= mass02 && can12){
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
void ReconstructObjects(vector<int> selectedJetsID, vector<TLorentzVector> Muons,vector<TLorentzVector> selectedElectrons, vector<TLorentzVector> selectedJets,int Region, bool threelepregion){
  clock_t start_sub =  clock();
  
  
  
  
  //cout << "SMjetIndex "<< SMjetIndex << " selectedJets.size() " << selectedJets.size() << endl;
  SMjetIndex = -999;
  SMjetIndex = SMjetCalculator(selectedJets,verbose);
  if(SMjetIndex > selectedJets.size())  cout << "WARNING SMjetIndex "<< SMjetIndex << endl;
  //SMbjet.SetPtEtaPhiE(selectedJets[SMjetIndex].Pt(),selectedJets[SMjetIndex].Eta(), selectedJets[SMjetIndex].Phi(), selectedJets[SMjetIndex].Energy());
  SMbjet.SetPxPyPzE(selectedJets[SMjetIndex].Px(),selectedJets[SMjetIndex].Py(), selectedJets[SMjetIndex].Pz(), selectedJets[SMjetIndex].Energy());
  // cout << "SMjetIndex new"<< SMjetIndex << endl;
  
  
  bool Wlepfound = false;
  if(WmuIndiceF != -999 && WelecIndiceF != -999) cout << "ERROR: 2 W leptons found" << endl;
  else if(WmuIndiceF!=-999 && threelepregion){
    Wlep.SetPtEtaPhiE(selectedMuons[WmuIndiceF].Pt(),selectedMuons[WmuIndiceF].Eta(), selectedMuons[WmuIndiceF].Phi(), selectedMuons[WmuIndiceF].Energy());
    Wlepfound = true;
  }
  else if(WelecIndiceF!=-999 && threelepregion){
    Wlep.SetPtEtaPhiE(selectedElectrons[WelecIndiceF].Pt(), selectedElectrons[WelecIndiceF].Eta(), selectedElectrons[WelecIndiceF].Phi(), selectedElectrons[WelecIndiceF].Energy());
    Wlepfound = true;
  }
  else if((selectedMuons.size() +selectedElectrons.size()) == 3){ cout << " ERROR: Wlep not found" << endl;}
  else { Wlepfound = false; }
  
  metTLVbf.SetPtEtaPhiE(met_Pt, met_Eta, met_Phi, TMath::Sqrt(met_Px*met_Px+met_Py*met_Py));
  if(Wlepfound){
    metTLV = MetzCalculator(Wlep, metTLVbf, SMbjet);
    //cout << "met" << endl;
    Wboson = metTLV + Wlep;
    SMtop = metTLV + Wlep  + SMbjet;
    
    //cout << "Wboson mass " << Wboson.M() << " SMtop " << SMtop.M() << " mlb " << (Wlep+SMbjet).M() << endl;
    mWT = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_Px)*(Wlep.Px() + met_Px) - (Wlep.Py() + met_Py)* (Wlep.Py() +met_Py));
    Double_t phis = Wlep.Phi() - met_Phi;
    Double_t cosphis = TMath::Cos(phis);
    mWT2 = TMath::Sqrt(2*Wlep.Pt() * met_Pt*(1-cosphis));
  }
  
  if(ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999 && ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999 ) cout << "ERROR: 2 Zbosons found " << endl;
  else if(ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999 ) Zboson = selectedMuons[ZmuIndiceF_0] + selectedMuons[ZmuIndiceF_1];
  else if(ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999 ) Zboson = selectedElectrons[ZelecIndiceF_0] + selectedElectrons[ZelecIndiceF_1];
  else cout << "ERROR: Zboson not found " << endl;
  
  if(selectedJets.size()>1 )
  {
    //cout << "cjet" << endl;
    cjetindex = FCNCjetCalculator(selectedJets,Zboson ,SMjetIndex, 3);
    cjetindex_CvsLtagger = FCNCjetCalculatorCvsLTagger(selectedJets,SMjetIndex, 3);
    cjetindex_CvsBtagger = FCNCjetCalculatorCvsBTagger(selectedJets,SMjetIndex, 3);
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
int SMjetCalculator(vector<TLorentzVector> Jets,int verb){
  clock_t start_sub = clock();
  // cout << "calculating bjet" << endl;
  int index_ = -999 ;
  Double_t tempbdis = -999.;
  
  for(int iJ = 0; iJ < Jets.size() ; iJ++){
    
    if(bdisc_jet[selectedJetsID[iJ]] > tempbdis)
    {
      index_ = iJ;
      tempbdis = bdisc_jet[selectedJetsID[iJ]];
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
void GetMetaData(TTree* tree, bool isData,int Entries, bool isAMC, bool isfakes){
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
   tree->SetBranchAddress("nbEv_Initial_3e =",&nbEv_Initial_2e1mu =,&b_nbEv_Initial_all);
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
  
  if(!isData && !isfakes){
    std::cout.precision (8);
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
/*
 TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met){
 
 clock_t start_sub = clock();
 // cout << "calculating metz" << endl;
 
 Double_t term1 = leptW.Pz() * ( leptW.Px()* v_met.Px() + leptW.Py()*v_met.Py() + pow(80.399, 2)/2.);
 
 Double_t det = pow(leptW.Px() * v_met.Px() + leptW.Py() * v_met.Py() + pow(80.399, 2)/2., 2) - v_met.Pt()*v_met.Pt() * (leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz() );
 
 if(det<0) det=0;
 
 Double_t term2 = leptW.E() * pow(det, 0.5);
 Double_t denom = leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz();
 Double_t sol1 = (term1 - term2) / denom;
 // Double_t sol2 = (term1 + term2) / denom;
 Double_t nu_E = 0;
 
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
 
 
 };*/
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met, TLorentzVector SMjet){
  // Get the z component of the neutrino by constraining W mass
  
  // We assume that the x and y components of the MET are entirely
  // due to the escaping neutrino, and apply WMass constraint in
  // order to extract the z component.
  TLorentzVector neutrino;
  double WMass = 80.4;
  
  // Neutrino momentum components
  double neutrinoPx, neutrinoPy, neutrinoPz, neutrinoE;
  neutrinoPx = v_met.Px();
  neutrinoPy = v_met.Py();
  neutrinoPz = -999;
  neutrinoE  = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(0,2), 0.5);//neglecting neutrino mass
  
  // Lepton momentum components
  double leptonPx, leptonPy, leptonPz, leptonE;
  leptonPx = leptW.Px() ;
  leptonPy =  leptW.Py() ;
  leptonPz =  leptW.Pz() ;
  leptonE  =  leptW.Energy();
  
  // Tranverse momentum of the lepton
  double trmomlep = 4.*( leptonE*leptonE - leptonPz*leptonPz );
  if( trmomlep < 1e-5 ) { neutrinoPz = -1000.; neutrino.SetPxPyPzE( v_met.Px(), v_met.Py(), neutrinoPz, neutrinoE); met_Pz = neutrinoPz; return neutrino;}
  
  //////////////////////////////////////////////////////////////////////////
  // If we assume that the neutrino and the lepton come from the W boson
  // we have:
  //   1) Momentum conservation : p(W) = p(l) + p(nu)
  //   2) Energy conservation : [M(W)]^2 = [p(l) + p(nu)]^2
  //
  // Solving the quadratic equation that yields, we have:
  //
  //                    pz(nu) = [-b +/- sqrt(b*b -4*a*c)]/(2*a)
  //
  // If the solution is complex, we only take the real part:
  //
  //                    pz(nu) = -b/2a
  //
  // Else, we have two solutions: A+B, A-B. Take value giving the top mass
  //////////////////////////////////////////////////////////////////////////
  Double_t TopMass = 172.9;
  double a, b, c;
  
  a = leptonPx*leptonPx + leptonPy*leptonPy;
  
  b = -( WMass*WMass*leptonPz
        + 2*leptonPx*leptonPz*neutrinoPx
        + 2*leptonPy*leptonPz*neutrinoPx );
  
  c = neutrinoPy * ( leptonPx*leptonPx + leptonPz*leptonPz )
  + neutrinoPx * ( leptonPy*leptonPy + leptonPz*leptonPz );
  
  // Square root
  double sqRoot = b*b - 4*a*c;
  
  // Parts of the solution
  double firstPart  = -b/(2*a);
  double secondPart = sqrt(b*b -4*a*c)/(2*a);
  
  // Solutions of the equation
  double neutrinoPz_1, neutrinoPz_2;
  double nu_E1, nu_E2, nu_E;
  TLorentzVector nu1, nu2;
  Double_t InvMass1, InvMass2;
  if( sqRoot<0 ){
    
    neutrinoPz = firstPart;
    
  }
  else{
    
    neutrinoPz_1 = firstPart + secondPart;
    neutrinoPz_2 = firstPart - secondPart;
    
    nu_E1 = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(neutrinoPz_1,2), 0.5);//neglecting neutrino mass
    nu_E2 = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(neutrinoPz_2,2), 0.5);//neglecting neutrino mass
    
    nu1.SetPxPyPzE( v_met.Px(), v_met.Py(), neutrinoPz_1, nu_E1);
    nu1.SetPxPyPzE( v_met.Px(), v_met.Py(), neutrinoPz_2, nu_E2);
    // Take solution closets to top
    InvMass1 = fabs((SMjet+nu1+leptW).M() - TopMass);
    InvMass2 = fabs((SMjet+nu2+leptW).M() - TopMass);
    if( InvMass1 < InvMass2 ) {
      neutrinoPz = neutrinoPz_1;
    }
    else {
      neutrinoPz = neutrinoPz_2;
    }
    
  }
  
  
  //printf("Neutrino pz = %f \n", neutrinoPz);
  
  met_Pz = neutrinoPz;
  nu_E = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(neutrinoPz,2), 0.5);//neglecting neutrino mass
  neutrino.SetPxPyPzE( v_met.Px(), v_met.Py(), neutrinoPz, nu_E);
  return neutrino;
  
}

std::pair <Double_t,Double_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo){
  
  clock_t start_sub = clock();
  // see https://github.com/TopBrussels/TopTreeAnalysis/blob/CMSSW_53X/WHelicities/src/BTagCosThetaCalculation.cc
  // ttbar http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2012/157
  // single top : https://cds.cern.ch/record/1601800
  
  
  Double_t CosTheta = -999.;
  
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
  Double_t ThetaTevatron = ROOT::Math::VectorUtil::Angle( TopWRF, leptWRF );
  CosTheta = -(TMath::Cos(ThetaTevatron));
  //Cos theta is defined as the angle between the lepton and the reversed direction of the top quark, both boosted to the W-boson rest frame.
  //Still reversed direction doesn't need to be calculated since angle between lepton and top and between lepton and reversed top is proportional to theta and Pi-theta.
  //For these the angles the following relation holds: cos(theta) = - cos(Pi-theta)
  // --> Hence the need of the minus sign in the CosTheta definition!!
  
  Double_t ThetaSTAR = ROOT::Math::VectorUtil::Angle( WTRF, leptWRF );
  Double_t CosThetaSTAR= TMath::Cos(ThetaSTAR);
  
  
  
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
  
  
  std::pair <Double_t,Double_t> CosThetaPair(CosTheta, CosThetaSTAR);
  
  
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
  Double_t TempMinMass = 100000.00;
  Double_t TopMass = 172.9;
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
  
  Double_t tempcdis = -9999.;
  int NbInColl = -999;
  
  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      
      if(cdiscCvsB_jet[selectedJetsID[iJ]]>tempcdis){
        NbInColl = iJ;
        tempcdis = cdiscCvsB_jet[selectedJetsID[iJ]];
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
  Double_t tempcdis = -9999.;
  int NbInColl = -999;
  
  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      
      if(cdiscCvsL_jet[selectedJetsID[iJ]]>tempcdis){
        NbInColl = iJ;
        tempcdis = cdiscCvsL_jet[selectedJetsID[iJ]];
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
  Double_t TempMinMass = 100000.00;
  Double_t TopMass = 172.9;
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

///////////////////////////////////// BOOKKEEPING /////////////////////////////////////////
string ConvertIntToString(int Number, int pad){
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

///////////////////////////////////// CLEARING /////////////////////////////////////////
void ClearObjects(){
  //  ClearLeaves(); TO FIX
  ClearTLVs();
  ClearVars();
  
}
void ClearMVAVars(){
  MVA_channel = -999;
  MVA_weight = 1.;
  MVA_weight_nom = 1.;
  MVA_weight_puSF_up = 1.;
  MVA_weight_puSF_down = 1.;
  MVA_weight_electronSF_up =1.;
  MVA_weight_electronSF_down = 1.;
  MVA_weight_muonSF_up = 1.;
  MVA_weight_muonSF_down = 1.;
  MVA_weight_btagSF_cferr1_up = 1.;
  MVA_weight_btagSF_cferr1_down = 1.;
  MVA_weight_btagSF_cferr2_up = 1.;
  MVA_weight_btagSF_cferr2_down = 1.;
  MVA_weight_btagSF_hf_up = 1.;
  MVA_weight_btagSF_hf_down = 1.;
  MVA_weight_btagSF_hfstats1_up = 1.;
  MVA_weight_btagSF_hfstats1_down = 1.;
  MVA_weight_btagSF_hfstats2_up = 1.;
  MVA_weight_btagSF_hfstats2_down = 1.;
  MVA_weight_btagSF_lf_up = 1.;
  MVA_weight_btagSF_lf_down = 1.;
  MVA_weight_btagSF_lfstats1_up =1.;
  MVA_weight_btagSF_lfstats1_down = 1.;
  MVA_weight_btagSF_lfstats2_up = 1.;
  MVA_weight_btagSF_lfstats2_down = 1.;
  
  
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
  MVA_nMuons = -999;
  MVA_NJets_CSVv2T = -999;
  MVA_NJets_CSVv2M = -999;
  MVA_NJets_CSVv2L = -999;
  MVA_nJets = -999;
  MVA_nJets_CharmL = -999;
  MVA_nJets_CharmM = -999;
  MVA_nJets_CharmT = -999;
  MVA_nElectrons = -999;
  
  //SM kinematics
  MVA_mWt = -999.;
  MVA_mWt2 = -999.;
  MVA_SMtop_M = -999.;
  MVA_mlb = -999.;
  MVA_Wboson_M = -999.;
  
  MVA_dRWlepb = -999.;
  
  MVA_dPhiWlepb = -999.;
  
  MVA_Wlep_Charge = -999;
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
  btagSFshape = 1.;
  scaleFactor_bfBT = 1.;
  scaleFactor_bfELSF = 1.;
  scaleFactor_bfMuSF = 1.;
  scaleFactor_bfPU = 1.;
  scaleFactor_puSF_down=1.;
  scaleFactor_puSF_up=1.;
  scaleFactor_electronSF_down=1.;
  scaleFactor_electronSF_up=1.;
  scaleFactor_muonSF_down=1.;
  scaleFactor_muonSF_up=1.;
  scaleFactor_btagSF_cferr1_down=1.;
  scaleFactor_btagSF_cferr1_up=1.;
  scaleFactor_btagSF_cferr2_down=1.;
  scaleFactor_btagSF_cferr2_up=1.;
  scaleFactor_btagSF_hf_down=1.;
  scaleFactor_btagSF_hf_up=1.;
  scaleFactor_btagSF_hfstats1_down=1.;
  scaleFactor_btagSF_hfstats1_up=1.;
  scaleFactor_btagSF_hfstats2_down=1.;
  scaleFactor_btagSF_hfstats2_up=1.;
  scaleFactor_btagSF_lf_down=1.;
  scaleFactor_btagSF_lf_up=1.;
  scaleFactor_btagSF_lfstats1_down=1.;
  scaleFactor_btagSF_lfstats1_up=1.;
  scaleFactor_btagSF_lfstats2_down=1.;
  scaleFactor_btagSF_lfstats2_up=1.;
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
  selectedJetsID.clear();
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
void ClearMatchingVarsTLV(){
  mcpart.Clear();
  mcParticles.clear();
  
}
void ClearMatchingVars(){
  
  foundTopQ = false;
  foundAntitopQ = false;
  foundSMb = false;
  foundSMmu = false;
  foundSMel = false;
  foundSMnumu = false;
  foundSMnuel = false;
  foundW = false;
  foundZmumin = false;
  foundZmuplus = false;
  foundZelmin = false;
  foundZelplus = false;
  foundZ = false;
  foundCjet = false;
  foundUjet = false;
  TopQ_Indice = -999;
  AntiTopQ_Indice = -999;
  SMb_Indice = -999;
  SMmu_Indice = -999;
  SMel_Indice = -999;
  SMW_Indice = -999;
  SMnumu_Indice = -999;
  SMnuel_Indice = -999;
  Zmu_min_Indice = -999;
  Zmu_plus_Indice = -999;
  Zel_min_Indice = -999;
  Zel_plus_Indice = -999;
  Z_Indice = -999;
  Cjet_Indice = -999;
  Ujet_Indice = -999;
  foundDecay = false;
  
  WmuIndiceM = -999;
  WelecIndiceM = -999;
  ZelecIndiceM_0 = -999;
  ZelecIndiceM_1 = -999;
  ZmuIndiceM_0 = -999;
  ZmuIndiceM_1 = -999;
  BjetIndiceM = -999;
  CjetIndiceM = -999;
  UjetIndiceM = -999;
}
void ClearMetaData(){
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
  WPc_CvsL_Tight= 0.;
  
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
void ClearMatchingSampleVars(){
  WlepMatched= 0.;
  ZlepMatched = 0.;
  WlepMatchedevent= 0.;
  ZlepMatchedevent = 0.;
  BjetMatched = 0.;
  BjetMatchedevent = 0.;
  CjetMatched = 0.;
  CjetMatchedevent = 0. ;
  UjetMatched = 0.;
  UjetMatchedevent = 0.;
  
  CjetMatchedeventCvsLL = 0.;
  CjetMatchedeventCvsLM = 0.;
  CjetMatchedeventCvsLT = 0.;
  CjetMatchedeventCvsBL = 0.;
  CjetMatchedeventCvsBM = 0.;
  CjetMatchedeventCvsBT = 0.;
  CjetMatchedeventL = 0.;
  CjetMatchedeventM = 0.;
  CjetMatchedeventT = 0.;
  CjetMatchedL = 0.;
  CjetMatchedM = 0.;
  CjetMatchedT = 0.;
  CjetMatchedCvsBL = 0.;
  CjetMatchedCvsBM = 0.;
  CjetMatchedCvsBT = 0.;
  CjetMatchedCvsLL = 0.;
  CjetMatchedCvsLM = 0.;
  CjetMatchedCvsLT = 0.;
}
///////////////////////////////////// INIT PLOTS /////////////////////////////////////////
void InitFakeValidation(string dataSetName,  vector <int> decayChannels){
  TH1::SetDefaultSumw2();
  string channelstr = "";
  
  vector<string> v_prefixregion = {"2lep", "3lep"};
  if(!doDilep) v_prefixregion = {"3lep"};
  string prefixregion = "";
  for(int iReg = 0; iReg < v_prefixregion.size(); iReg++){
    prefixregion = v_prefixregion[iReg];
    
    for(int iChan =0; iChan < decayChannels.size() ; iChan++){
      
      channelstr = prefixregion + "_";
      if(decayChannels[iChan] == 0) channelstr += "uuu";
      if(decayChannels[iChan] == 1) channelstr += "uue";
      if(decayChannels[iChan] == 2) channelstr += "eeu";
      if(decayChannels[iChan] == 3) channelstr += "eee";
      if(decayChannels[iChan] == -9) channelstr += "all";
      if(decayChannels[iChan] == 4) channelstr +="uu";
      if(decayChannels[iChan] == 5) channelstr +="ee";
      //MSPlot[plotname.c_str()]->setChannel(true, decayChan);
      
      if((decayChannels[iChan] == 4 || decayChannels[iChan] == 5) && prefixregion.find("3lep")!=std::string::npos) continue;
      
      histo1D_fakevvalidation[("ZbosonPt_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPt_"+channelstr+dataSetName).c_str(),"p_{T} Z boson (GeV)",20,0,500);
      histo1D_fakevvalidation[("ZbosonEta_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonEta_"+channelstr+dataSetName).c_str(),"#eta Z boson",6,-3,3);
      histo1D_fakevvalidation[("ZbosonPhi_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPhi_"+channelstr+dataSetName).c_str(),"#phi Z boson",9,-4,4);
      
      histo1D_fakevvalidation[("WlepPt_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPt_"+channelstr+dataSetName).c_str(),"p_{T} l_{W} (GeV)",10,0,450);
      histo1D_fakevvalidation[("WlepEta_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepEta_"+channelstr+dataSetName).c_str(),"#eta l_{W}",6,-3,3);
      histo1D_fakevvalidation[("WlepPhi_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPhi_"+channelstr+dataSetName).c_str(),"#phi l_{W}",8,-4,4);
      histo1D_fakevvalidation[("TrMassW_"+channelstr+dataSetName).c_str()] = new TH1F(("TrMassW_"+channelstr+dataSetName).c_str(),"transv. mass W boson (GeV)",10,15,300);
      
      histo1D_fakevvalidation[("ZbosonPtMu_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPtMu_"+channelstr+dataSetName).c_str(),"p_{T} Z_{#mu #mu} boson (GeV)",10,0,500);
      histo1D_fakevvalidation[("ZbosonEtaMu_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonEtaMu_"+channelstr+dataSetName).c_str(),"#eta Z_{#mu #mu} boson",6,-3,3);
      histo1D_fakevvalidation[("ZbosonPhiMu_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPhiMu_"+channelstr+dataSetName).c_str(),"#phi Z_{#mu #mu} boson",9,-4,4);
      
      histo1D_fakevvalidation[("WlepPtMu_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPtu_"+channelstr+dataSetName).c_str(),"p_{T} #mu_{W} (GeV)",10,0,500);
      histo1D_fakevvalidation[("WlepEtaMu_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepEtau_"+channelstr+dataSetName).c_str(),"#eta #mu_{W}",6,-3,3);
      histo1D_fakevvalidation[("WlepPhiMu_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPhiu_"+channelstr+dataSetName).c_str(),"#phi #mu_{W}",8,-4,4);
      histo1D_fakevvalidation[("TrMassWMu_"+channelstr+dataSetName).c_str()] = new TH1F(("TrMassWu_"+channelstr+dataSetName).c_str(),"transv. mass W_{#mu} boson (GeV)",10,15,300);
      
      histo1D_fakevvalidation[("ZbosonPtEl_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPtEl_"+channelstr+dataSetName).c_str(),"p_{T} Z_{ee} boson (GeV)",10,0,500);
      histo1D_fakevvalidation[("ZbosonEtaEl_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonEtaEl_"+channelstr+dataSetName).c_str(),"#eta Z_{ee} boson",6,-3,3);
      histo1D_fakevvalidation[("ZbosonPhiEl_"+channelstr+dataSetName).c_str()] = new TH1F(("ZbosonPhiEl_"+channelstr+dataSetName).c_str(),"#phi Z_{ee} boson",9,-4,4);
      
      histo1D_fakevvalidation[("WlepPtEl_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPtEl_"+channelstr+dataSetName).c_str(),"p_{T} e_{W} (GeV)",10,0,500);
      histo1D_fakevvalidation[("WlepEtaEl_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepEtaEl_"+channelstr+dataSetName).c_str(),"#eta e_{W}",6,-3,3);
      histo1D_fakevvalidation[("WlepPhiEl_"+channelstr+dataSetName).c_str()] = new TH1F(("WlepPhiEl_"+channelstr+dataSetName).c_str(),"#phi e_{W}",8,-4,4);
      histo1D_fakevvalidation[("TrMassWEl_"+channelstr+dataSetName).c_str()] = new TH1F(("TrMassWEl_"+channelstr+dataSetName).c_str(),"transv. mass W_{e} boson (GeV)",10,15,300);
    }
  }
  
}

void Init1DPlots(string dataSetName){
  TH1::SetDefaultSumw2();
  
  //cout << ("NbVertices_nom_"+dataSetName).c_str() << endl;
  histo1D_PUSystematics[("NbVertices_nom_"+dataSetName).c_str()] = new TH1F(("NbVertices_nom_"+dataSetName).c_str(),"Nb Of Vertices",30,-0.5,59);
  histo1D_ElSystematics[("Pt_electron_nom_"+dataSetName).c_str()] = new TH1F(("Pt_electron_nom_"+dataSetName).c_str(),"Pt leading electron",50,0,600);
  histo1D_MuSystematics[("Pt_muon_nom_"+dataSetName).c_str()] = new TH1F(("Pt_muon_nom_"+dataSetName).c_str(),"Pt leading muon",50,0,600);
  
  histo1D_PUSystematics[("NbVertices_up_"+dataSetName).c_str()] = new TH1F(("NbVertices_up_"+dataSetName).c_str(),"Nb Of Vertices",30,-0.5,59);
  histo1D_ElSystematics[("Pt_electron_up_"+dataSetName).c_str()] = new TH1F(("Pt_electron_up_"+dataSetName).c_str(),"Pt leading electron",50,0,600);
  histo1D_MuSystematics[("Pt_muon_up_"+dataSetName).c_str()] = new TH1F(("Pt_muon_up_"+dataSetName).c_str(),"Pt leading muon",50,0,600);
  
  histo1D_PUSystematics[("NbVertices_down_"+dataSetName).c_str()] = new TH1F(("NbVertices_down_"+dataSetName).c_str(),"Nb Of Vertices",30,-0.5,59);
  histo1D_ElSystematics[("Pt_electron_down_"+dataSetName).c_str()] = new TH1F(("Pt_electron_down_"+dataSetName).c_str(),"Pt leading electron",50,0,600);
  histo1D_MuSystematics[("Pt_muon_down_"+dataSetName).c_str()] = new TH1F(("Pt_muon_down_"+dataSetName).c_str(),"Pt leading muon",50,0,600);
  
  histo1D_Bcferr1Systematics[("Bdis_cferr1_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr1_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bcferr1Systematics[("Bdis_cferr1_up_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr1_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bcferr1Systematics[("Bdis_cferr1_down_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr1_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bcferr2Systematics[("Bdis_cferr2_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr2_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bcferr2Systematics[("Bdis_cferr2_up_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr2_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bcferr2Systematics[("Bdis_cferr2_down_"+dataSetName).c_str()] = new TH1F(("Bdis_cferr2_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  
  histo1D_Bhfstats1Systematics[("Bdis_hfstats1_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats1_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bhfstats1Systematics[("Bdis_hfstats1_up_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats1_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bhfstats1Systematics[("Bdis_hfstats1_down_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats1_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bhfstats2Systematics[("Bdis_hfstats2_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats2_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bhfstats2Systematics[("Bdis_hfstats2_up_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats2_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Bhfstats2Systematics[("Bdis_hfstats2_down_"+dataSetName).c_str()] = new TH1F(("Bdis_hfstats2_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  
  histo1D_Blfstats1Systematics[("Bdis_lfstats1_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats1_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Blfstats1Systematics[("Bdis_lfstats1_up_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats1_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Blfstats1Systematics[("Bdis_lfstats1_down_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats1_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Blfstats2Systematics[("Bdis_lfstats2_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats2_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Blfstats2Systematics[("Bdis_lfstats2_up_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats2_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_Blfstats2Systematics[("Bdis_lfstats2_down_"+dataSetName).c_str()] = new TH1F(("Bdis_lfstats2_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  
  histo1D_BhfSystematics[("Bdis_hf_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_hf_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_BhfSystematics[("Bdis_hf_up_"+dataSetName).c_str()] = new TH1F(("Bdis_hf_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_BhfSystematics[("Bdis_hf_down_"+dataSetName).c_str()] = new TH1F(("Bdis_hf_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_BlfSystematics[("Bdis_lf_nom_"+dataSetName).c_str()] = new TH1F(("Bdis_lf_nom_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_BlfSystematics[("Bdis_lf_up_"+dataSetName).c_str()] = new TH1F(("Bdis_lf_up_"+dataSetName).c_str(),"CSVv2",25,0,1);
  histo1D_BlfSystematics[("Bdis_lf_down_"+dataSetName).c_str()] = new TH1F(("Bdis_lf_down_"+dataSetName).c_str(),"CSVv2",25,0,1);
  
  
  
  
}

void InitMSPlotsBDT(string prefix, vector <int> decayChannels){
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
    
    
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), 60, 0, 60, "BDT");
    
  }
  
  
  
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

void InitMSPlots(string prefix, vector <int> decayChannels){
  clock_t start_sub = clock();
  string decaystring = "";
  // control plots
  vector<string> v_prefixregion = {"2lep", "3lep"};
  if(!doDilep) v_prefixregion = {"3lep"};
  string prefixregion = "";
  bool is3regio = false;
  
  for(int iReg = 0; iReg < v_prefixregion.size(); iReg++){
    prefixregion = v_prefixregion[iReg];
    if(prefixregion.find("3lep")!=std::string::npos){ is3regio = true;}
    for(int iChan =0; iChan < decayChannels.size() ; iChan++){
      
      decaystring = "";
      if(decayChannels[iChan] == 0) decaystring = "uuu";
      if(decayChannels[iChan] == 1) decaystring = "uue";
      if(decayChannels[iChan] == 2) decaystring = "eeu";
      if(decayChannels[iChan] == 3) decaystring = "eee";
      if(decayChannels[iChan] == -9) decaystring = "all";
      if(decayChannels[iChan] == 4) decaystring ="uu";
      if(decayChannels[iChan] == 5) decaystring ="ee";
      //MSPlot[plotname.c_str()]->setChannel(true, decayChan);
      
      if((decayChannels[iChan] == 4 || decayChannels[iChan] == 5) && prefixregion.find("3lep")!=std::string::npos) continue;
      
      MSPlot[(prefixregion+prefix+"_NbOfVertices_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_NbOfVertices_"+decaystring).c_str(), 60, 0, 60, "#  vertices");
      MSPlot[(prefixregion+prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_NbOfVertices_bfPU_"+decaystring).c_str(), 60, 0, 60, "#  vertices before PU reweighing");
      MSPlot[(prefixregion+prefix+"_puSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_puSF_"+decaystring).c_str(), 200, 0, 2, "Pile up SF");
      
      //cout << "init " << (prefixregion+prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
      MSPlot[(prefixregion+prefix+"_bdisc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_bdisc_"+decaystring).c_str(), 21, 0., 1, "CSVv2 discriminant");
      MSPlot[(prefixregion+prefix+"_bdisc_bfBT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_bdisc_bfBT_"+decaystring).c_str(),21, 0., 1, "CSVv2 discriminant before Btag SF");
      MSPlot[(prefixregion+prefix+"_btagSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_btagSF_"+decaystring).c_str(), 80, 0.9, 1.3, "Btag SF");
      
      
      MSPlot[(prefixregion+prefix+"_cvsldisc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_cvsldisc_"+decaystring).c_str(), 16, -0.6., 1, "charm vs loose disc.");
      MSPlot[(prefixregion+prefix+"_cvsbdisc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_cvsbdisc_"+decaystring).c_str(), 21, -1., 1, "charm vs b disc.");
      
      MSPlot[(prefixregion+prefix+"_nMu_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nMu_"+decaystring).c_str(), 10, -0.5, 9.5, "#  Muons");
      MSPlot[(prefixregion+prefix+"_nMu_bfMuSF_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nMu_bfMuSF_"+decaystring).c_str(), 10, -0.5, 9.5, "#  Muons before Muon SF");
      MSPlot[(prefixregion+prefix+"_muSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_muSF_"+decaystring).c_str(), 60, 0.75, 1.05, "Muon SF");
      
      MSPlot[(prefixregion+prefix+"_nEl_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nEl_"+decaystring).c_str(), 10, -0.5, 9.5, "#  Electrons");
      MSPlot[(prefixregion+prefix+"_nEl_bfElSF_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nEl_bfElSF_"+decaystring).c_str(), 10, -0.5, 9.5, "#  Electrons before Electron SF","# e");
      MSPlot[(prefixregion+prefix+"_elSF_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_elSF_"+decaystring).c_str(), 60, 0.8, 1.05, "Electron SF");
      
      MSPlot[(prefixregion+prefix+"_JetPt_bfJER_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_JetPt_bfJER_"+decaystring).c_str(), 20, 0, 500, "Jet p_{T} before JER, after JES ", "GeV");
      MSPlot[(prefixregion+prefix+"_JetPt_bfJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_JetPt_bfJES_"+decaystring).c_str(), 20, 0, 500, "Jet p_{T} before JES, before JER ", "GeV");
      MSPlot[(prefixregion+prefix+"_JetPt_afJER_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_JetPt_afJER_"+decaystring).c_str(), 15, 0, 500, "Jet p_{T}","GeV");
      MSPlot[(prefixregion+prefix+"_JetPt_afJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_JetPt_afJES_"+decaystring).c_str(), 20, 0, 500, "Jet p_{T} after JES, before JER ", "GeV");
      MSPlot[(prefixregion+prefix+"_met_bfJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_bfJES_"+decaystring).c_str(), 40, 0, 400, "E_{T}^{miss} p_{T} before JES ","GeV");
      MSPlot[(prefixregion+prefix+"_met_afJES_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_afJES_"+decaystring).c_str(), 40, 0, 400, "E_{T}^{miss} p_{T}","GeV");
      MSPlot[(prefixregion+prefix+"_met_phi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_phi_"+decaystring).c_str(), 20, -4, 4, "E_{T}^{miss} #phi");
      MSPlot[(prefixregion+prefix+"_met_px_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_px_"+decaystring).c_str(), 25, 0, 50, "E_{T}^{miss} p_{x}");
      MSPlot[(prefixregion+prefix+"_met_py_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_py_"+decaystring).c_str(), 25, 0, 50, "E_{T}^{miss} p_{y}");
      MSPlot[(prefixregion+prefix+"_met_pz_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_pz_"+decaystring).c_str(), 25, 0, 50, "E_{T}^{miss} p_{z}");
      MSPlot[(prefixregion+prefix+"_met_pt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_met_pt_"+decaystring).c_str(), 20, 0, 400, "E_{T}^{miss} p_{T}");
      
      
      // vars
      MSPlot[(prefixregion+prefix+"_ZbosonMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMass_"+decaystring).c_str(), 10, 70, 110, "inv. mass Z boson ","GeV");
      if(prefixregion.find("3lep")!=std::string::npos){
        MSPlot[(prefixregion+prefix+"_WbosonMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WbosonMass_"+decaystring).c_str(), 200, 0, 200, "inv. mass W boson ","GeV");
        MSPlot[(prefixregion+prefix+"_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_mlb_"+decaystring).c_str(), 40, 0, 800, "Inv. Mass l_{W}+b^{SM} ","GeV");
        MSPlot[(prefixregion+prefix+"_SMTopMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_SMTopMass_"+decaystring).c_str(), 100, 0, 500, "inv. mass SM top ","GeV");
        MSPlot[(prefixregion+prefix+"_mWT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_mWT_"+decaystring).c_str(), 25, 0, 300, "Transv. Mass W boson ","GeV");
        MSPlot[(prefixregion+prefix+"_mWT2_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_mWT2_"+decaystring).c_str(), 25, 0, 300, "Transv. Mass W boson ","GeV");
        MSPlot[(prefixregion+prefix+"_3dLeadingLepPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingLepPt_"+decaystring).c_str(), 25, 0, 500, "3d leading lepton p_{T} ","GeV");
        
      }
      MSPlot[(prefixregion+prefix+"_LeadingJetPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingJetPt_"+decaystring).c_str(), 20, 0, 500, "leading jet p_{T} ","GeV");
      MSPlot[(prefixregion+prefix+"_LeadingLepPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingLepPt_"+decaystring).c_str(), 25, 0, 500, "leading lepton p_{T} ","GeV");
      MSPlot[(prefixregion+prefix+"_2ndLeadingJetPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingJetPt_"+decaystring).c_str(), 15, 0, 300, "2nd leading jet p_{T} ","GeV");
      MSPlot[(prefixregion+prefix+"_2ndLeadingLepPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingLepPt_"+decaystring).c_str(), 25, 0, 500, "2nd leading lepton p_{T} ","GeV");
      
      
      
      MSPlot[(prefixregion+prefix+"_nJets_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJets_"+decaystring).c_str(), 10, -0.5, 9.5, "#  Jets");
      MSPlot[(prefixregion+prefix+"_nJetsCSVL_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCSVL_"+decaystring).c_str(), 5, -0.5, 4.5, "#  CSVL");
      MSPlot[(prefixregion+prefix+"_nJetsCSVM_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCSVM_"+decaystring).c_str(), 5, -0.5, 4.5, "#  CSVM");
      MSPlot[(prefixregion+prefix+"_nJetsCSVT_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCSVT_"+decaystring).c_str(), 5, -0.5, 4.5, "#  CSVT");
      MSPlot[(prefixregion+prefix+"_nJetsCharmL_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCharmL_"+decaystring).c_str(), 10, -0.5, 9.5, "#  charm loose jets");
      MSPlot[(prefixregion+prefix+"_nJetsCharmM_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCharmM_"+decaystring).c_str(), 10, -0.5, 9.5, "#  charm medium jets");
      MSPlot[(prefixregion+prefix+"_nJetsCharmT_"+decaystring).c_str()]  = new MultiSamplePlot(datasets, (prefixregion+prefix+"_nJetsCharmT_"+decaystring).c_str(), 10, -0.5, 9.5, "#  charm tight jets");
      // finding Z boson
      
      if(decayChannels[iChan] == 0 || decayChannels[iChan] == -9){
        MSPlot[(prefixregion+prefix+"_3dLeadingMuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingMuPt_"+decaystring).c_str(), 25, 0, 500, "3d leading muon p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_3dLeadingMuIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingMuIso_"+decaystring).c_str(), 10,0,0.5, "3d leading muon rel. iso.");
        MSPlot[(prefixregion+prefix+"_3dLeadingMuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingMuPhi_"+decaystring).c_str(), 20,-4,4, "3d leading muon #phi ");
        MSPlot[(prefixregion+prefix+"_3dLeadingMuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingMuEta_"+decaystring).c_str(), 15,-3,3, "3d leading muon #eta ");
        
        MSPlot[(prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str(), 15,3,3, "#mu_{W} #eta ");
        MSPlot[(prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str(), 20,-4,4, "#mu_{W} #phi ");
        MSPlot[(prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str(), 10,0,250, "#mu_{W}( p_{T} ","GeV");
      }
      if(decayChannels[iChan] == 0 || decayChannels[iChan] == 1 || decayChannels[iChan] == -9 || decayChannels[iChan] == 4 ){
        MSPlot[(prefixregion+prefix+"_2ndLeadingMuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingMuPt_"+decaystring).c_str(), 20, 0, 250, "2nd leading muon p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_2ndLeadingMuIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingMuIso_"+decaystring).c_str(), 10,0,0.5, "2nd leading muon rel. iso. ");
        MSPlot[(prefixregion+prefix+"_2ndLeadingMuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingMuEta_"+decaystring).c_str(), 15,-3,3, "2nd leading muon #eta ");
        
        MSPlot[(prefixregion+prefix+"_ZbosonMuIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMuIso_"+decaystring).c_str(), 10,0,0.5, "#mu_{Z} rel. iso.");
        MSPlot[(prefixregion+prefix+"_ZbosonMudPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMudPhi_"+decaystring).c_str(), 10,-4,4, "#Delta #phi (#mu_{Z},#mu_{Z})");
        MSPlot[(prefixregion+prefix+"_ZbosonMudR_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMudR_"+decaystring).c_str(), 20,0,4, "#Delta R(#mu,#mu) ");
        MSPlot[(prefixregion+prefix+"_ZbosonMassMu_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMassMu_"+decaystring).c_str(), 10, 70, 110, "inv. mass Z_{#mu,#mu} boson ","GeV");
        
        MSPlot[(prefixregion+prefix+"_ZbosonPtMu_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonPtMu_"+decaystring).c_str(), 40,0,500, "Z_{#mu,#mu} boson p_{T} ","GeV");
        
        MSPlot[(prefixregion+prefix+"_2ndLeadingMuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingMuPhi_"+decaystring).c_str(), 20,-4,4, "2nd leading muon #phi");
        
        if(decayChannels[iChan] == 1  ){
          MSPlot[(prefixregion+prefix+"_WlepElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElEta_"+decaystring).c_str(), 30,-6,6, "e_{W} #eta ");
          MSPlot[(prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str(), 20,-4,4, "e_{W} #phi");
          MSPlot[(prefixregion+prefix+"_WlepElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElPt_"+decaystring).c_str(), 25,0,400, "e_{W} p_{T} ","GeV");
        }
      }
      if(decayChannels[iChan] == 0 || decayChannels[iChan] == 1 || decayChannels[iChan] == 2 || decayChannels[iChan] == -9 || decayChannels[iChan] == 4){
        MSPlot[(prefixregion+prefix+"_LeadingMuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingMuPt_"+decaystring).c_str(), 20, 0, 300, "leading muon p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_MuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_MuPt_"+decaystring).c_str(), 20, 0, 300, "muon p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_LeadingMuIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingMuIso_"+decaystring).c_str(), 10,0,0.5, "leading muon rel. iso. ");
        MSPlot[(prefixregion+prefix+"_MuIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_MuIso_"+decaystring).c_str(), 10,0,0.5, "muon rel. iso. ");
        MSPlot[(prefixregion+prefix+"_LeadingMuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingMuPhi_"+decaystring).c_str(), 20,-4,4, "leading muon #phi ");
        MSPlot[(prefixregion+prefix+"_MuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_MuPhi_"+decaystring).c_str(), 20,-4,4, "muon #phi ");
        MSPlot[(prefixregion+prefix+"_LeadingMuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingMuEta_"+decaystring).c_str(), 15,-3,3, "muon #eta ");
        MSPlot[(prefixregion+prefix+"_MuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_MuEta_"+decaystring).c_str(), 15,-3,3, "leading muon #eta ");
      }
      if(decayChannels[iChan] == 3 || decayChannels[iChan] == -9){
        MSPlot[(prefixregion+prefix+"_3dLeadingElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingElPt_"+decaystring).c_str(), 10, 0, 150, "3d leading electron p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_3dLeadingElIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingElIso_"+decaystring).c_str(), 5,0,0.35, "3d leading electron rel. iso.");
        MSPlot[(prefixregion+prefix+"_3dLeadingElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingElPhi_"+decaystring).c_str(), 20,-4,4, "3d leading electron #phi");
        MSPlot[(prefixregion+prefix+"_3dLeadingElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_3dLeadingElEta_"+decaystring).c_str(), 15,-3,3, "3d leading electron #eta ");
        
        
        MSPlot[(prefixregion+prefix+"_WlepElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElEta_"+decaystring).c_str(), 15,-3,3, "e_{W} #eta ");
        MSPlot[(prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str(), 20,-4,4, "e_{W} #phi");
        MSPlot[(prefixregion+prefix+"_WlepElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepElPt_"+decaystring).c_str(), 25,0,500, "e_{W}p_{T} ","GeV");
        
        
        
      }
      if(decayChannels[iChan] == 3 || decayChannels[iChan] == 2 || decayChannels[iChan] == -9 || decayChannels[iChan] == 5){
        MSPlot[(prefixregion+prefix+"_2ndLeadingElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingElPt_"+decaystring).c_str(), 10, 0, 200, "2nd leading electron p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_2ndLeadingElIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingElIso_"+decaystring).c_str(), 5,0,0.35, "2nd leading electron rel. iso.");
        MSPlot[(prefixregion+prefix+"_ZbosonElIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonElIso_"+decaystring).c_str(), 5,0,0.35, "e_{Z}  rel. iso.");
        MSPlot[(prefixregion+prefix+"_2ndLeadingElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingElEta_"+decaystring).c_str(), 15,-3,3, "2nd leading electron #eta ");
        
        MSPlot[(prefixregion+prefix+"_ZbosonEldPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonEldPhi_"+decaystring).c_str(), 20,-4,4, "#Delta #phi(e_{Z},e_{Z})");
        MSPlot[(prefixregion+prefix+"_ZbosonEldR_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonEldR_"+decaystring).c_str(), 20,0,4, "#Delta R(e_{Z},e_{Z})");
        
        MSPlot[(prefixregion+prefix+"_ZbosonPtEl_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonPtEl_"+decaystring).c_str(), 20,0,500, "Z_{ee} boson p_{T}","GeV");
        
        MSPlot[(prefixregion+prefix+"_2ndLeadingElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_2ndLeadingElPhi_"+decaystring).c_str(), 20,-4,4, "2nd leading electron #phi");
        MSPlot[(prefixregion+prefix+"_ZbosonMassEl_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonMassEl_"+decaystring).c_str(), 10, 70, 110, "inv. mass Z_{ee} boson ","GeV");
        
        if(decayChannels[iChan] == 2  ){
          MSPlot[(prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str(), 15,-3,3, "#mu_{W}) #eta ");
          MSPlot[(prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str(), 20,-4,4, "#mu_{W}) #phi");
          MSPlot[(prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str(), 25,0,500, "#mu_{W} p_{T} ","GeV");
        }
        
      }
      if(decayChannels[iChan] == 3 || decayChannels[iChan] == 1 || decayChannels[iChan] == 2|| decayChannels[iChan] == -9 || decayChannels[iChan] == 5){
        MSPlot[(prefixregion+prefix+"_LeadingElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingElPt_"+decaystring).c_str(), 25, 0, 400, "leading electron p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_ElPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ElPt_"+decaystring).c_str(), 25, 0, 500, "electron p_{T} ","GeV");
        MSPlot[(prefixregion+prefix+"_LeadingElIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingElIso_"+decaystring).c_str(), 5,0,0.35, "leading electron rel. iso.");
        MSPlot[(prefixregion+prefix+"_ElIso_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ElIso_"+decaystring).c_str(), 5,0,0.35, "electron rel. iso.");
        MSPlot[(prefixregion+prefix+"_LeadingElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingElPhi_"+decaystring).c_str(), 20,-4,4, "leading electron #phi");
        MSPlot[(prefixregion+prefix+"_ElPhi_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ElPhi_"+decaystring).c_str(), 20,-4,4, "electron #phi");
        MSPlot[(prefixregion+prefix+"_LeadingElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_LeadingElEta_"+decaystring).c_str(), 15,-3,3, "leading electron #eta ");
        MSPlot[(prefixregion+prefix+"_ElEta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ElEta_"+decaystring).c_str(), 15,-3,3, "electron #eta ");
        
        
      }
      
      
      
      MSPlot[(prefixregion+prefix+"_ZbosonPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefixregion+prefix+"_ZbosonPt_"+decaystring).c_str(), 20,0,500, "Z boson p_{T} ","GeV");
      
      
      
      
      
    }
    
    
  }
  MSPlot[(prefix+"_Decay").c_str()]= new MultiSamplePlot(datasets, (prefix+"_Decay").c_str(), 10, -0.5, 9.5, "Decay channel");
  
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
    if(decayChannels[iChan] == 4 || decayChannels[iChan] == 5) continue;
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    MSPlot[(prefix+"_MVA_mWt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mWt_"+decaystring).c_str(),20, 0, 500, "Transv. Mass W boson  ", "GeV");
    MSPlot[(prefix+"_MVA_mWt2_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mWt2_"+decaystring).c_str(),15, 0, 300, "Transv. Mass W boson ","GeV");
  }
  
}
void InitMVAMSPlotsSingletop(string prefix, vector <int> decayChannels){
  clock_t start_sub = clock();
  
  InitMVAMSPlotsWZ(prefix, decayChannels);
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    if(decayChannels[iChan] == 4 || decayChannels[iChan] == 5) continue;
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    MSPlot[ (prefix+"_MVA_channel_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_channel_"+decaystring).c_str(), 5,-0.5, 4.5, "decaymode");
    MSPlot[ (prefix+"_MVA_weight_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_weight_"+decaystring).c_str(), 100,-0.5, 9.5, "eventweight");
    //   cout << "defining " <<  (prefix+"_MVA_lepton0_pt_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_MVA_lepton0_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_pt_"+decaystring).c_str(), 25,0, 500, "leading lepton p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_Zboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_pt_"+decaystring).c_str(), 25,0, 500, "Z boson p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_Zboson_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_eta_"+decaystring).c_str(),30,-6, 6, "Z boson #eta ","GeV");
    MSPlot[(prefix+"_MVA_dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRWlepb_"+decaystring).c_str(),30,0, 6, "#Delta R(l_{W},b)");
    MSPlot[(prefix+"_MVA_dPhiWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiWlepb_"+decaystring).c_str(),20,-4, 4, "#Delta #phi(l_{W},b)");
    MSPlot[(prefix+"_MVA_TotalHt_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_jet_"+decaystring).c_str(),20, 0, 1200, "total jet and E_{T}^{miss} H_{T} ","GeV");
    MSPlot[(prefix+"_MVA_dRZb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZb_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,b)");
    MSPlot[(prefix+"_MVA_dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZWlep_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,l_{W})");
    MSPlot[(prefix+"_MVA_dRZSMtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZSMtop_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,SM top)");
    
    MSPlot[(prefix+"_MVA_dPhiZb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZb_"+decaystring).c_str(),20,-4, 4, "#Delta #phi (Z,b)");
    MSPlot[(prefix+"_MVA_dPhiZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZWlep_"+decaystring).c_str(),20,-4, 4, "#Delta #phi (Z, l_{W})");
    MSPlot[(prefix+"_MVA_dPhiZMET_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZMET_"+decaystring).c_str(),20,-4, 4, "#Delta #phi (Z,E_{T}^{miss})");
    MSPlot[(prefix+"_MVA_dPhiZSMtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZSMtop_"+decaystring).c_str(),20,-4, 4, "#Delta #phi (Z,SM top)");
    MSPlot[(prefix+"_MVA_SMtop_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_eta_"+decaystring).c_str(),60,-6, 6, "SM top #eta");
    
    
    
    MSPlot[(prefix+"_MVA_lepton1_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_pt_"+decaystring).c_str(), 10,0, 250, "2nd leading lepton p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_lepton2_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_pt_"+decaystring).c_str(), 10,0, 250, "3d leading lepton p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_lepton0_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_eta_"+decaystring).c_str(),15,-3, 3, "leading lepton #eta");
    MSPlot[ (prefix+"_MVA_lepton1_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_eta_"+decaystring).c_str(),15,-3, 3, "2nd leading lepton #eta");
    MSPlot[ (prefix+"_MVA_lepton2_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_eta_"+decaystring).c_str(),15,-3,3 , "3d leading lepton #eta");
    MSPlot[ (prefix+"_MVA_lepton0_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton0_phi_"+decaystring).c_str(),20,-4, 4, "leading lepton #phi");
    MSPlot[ (prefix+"_MVA_lepton1_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton1_phi_"+decaystring).c_str(),20,-4, 4, "2nd leading lepton #phi");
    MSPlot[ (prefix+"_MVA_lepton2_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_lepton2_phi_"+decaystring).c_str(),20,-4, 4, "3d leading lepton #phi");
    
    MSPlot[ (prefix+"_MVA_jet0_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_pt_"+decaystring).c_str(), 20,0, 500, "leading jet p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_jet0_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_eta_"+decaystring).c_str(),15,-3, 3, "leading jet #eta");
    MSPlot[ (prefix+"_MVA_jet0_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet0_phi_"+decaystring).c_str(),20,-4, 4, "leading jet #phi");
    
    // SM side
    MSPlot[ (prefix+"_MVA_Wlep_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_pt_"+decaystring).c_str(), 15,0, 300, "l_{W} p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_Wlep_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_eta_"+decaystring).c_str(),15,-3, 3, "l_{W} #eta ");
    MSPlot[ (prefix+"_MVA_Wlep_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_phi_"+decaystring).c_str(),20,-4, 4, "l_{W} #phi ");
    MSPlot[ (prefix+"_MVA_SMbjet_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_pt_"+decaystring).c_str(), 20,0, 400, "b jet p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_SMbjet_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_eta_"+decaystring).c_str(),15,-3, 3, "b jet #eta ");
    MSPlot[ (prefix+"_MVA_SMbjet_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMbjet_phi_"+decaystring).c_str(),20,-4, 4, "b jet #phi ");
    MSPlot[ (prefix+"_MVA_Wboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_pt_"+decaystring).c_str(), 10,0, 400, "W boson p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_Wboson_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_eta_"+decaystring).c_str(),30,-6, 6, "W boson #eta");
    MSPlot[ (prefix+"_MVA_Wboson_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_phi_"+decaystring).c_str(),20,-4, 4, "W boson #phi");
    MSPlot[ (prefix+"_MVA_met_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_met_"+decaystring).c_str(), 20,0, 400, "E_{T}^{miss} p_{T}","GeV");
    MSPlot[ (prefix+"_MVA_SMtop_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_pt_"+decaystring).c_str(), 15,0, 350, "SM top p_{T} ","GeV");
    MSPlot[ (prefix+"_MVA_SMtop_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_phi_"+decaystring).c_str(),20,-4, 4, "SM top #phi");
    
    // FCNC side
    MSPlot[ (prefix+"_MVA_Zboson_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_phi_"+decaystring).c_str(),20,-4, 4, "Z boson #phi");
    
    // nbrs
    MSPlot[ (prefix+"_MVA_nMuons_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nMuons_"+decaystring).c_str(), 10,-0.5, 9.5, "# Muons");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2T_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2T_"+decaystring).c_str(), 10,-0.5, 9.5, "# CSVv2T");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2M_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2M_"+decaystring).c_str(), 10,-0.5, 9.5, "# CSVv2M");
    MSPlot[ (prefix+"_MVA_NJets_CSVv2L_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_NJets_CSVv2L_"+decaystring).c_str(), 10,-0.5, 9.5, "# CSVv2L");
    MSPlot[ (prefix+"_MVA_nJets_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_"+decaystring).c_str(), 10,-0.5, 9.5, "# Jets");
    MSPlot[ (prefix+"_MVA_nElectrons_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nElectrons_"+decaystring).c_str(), 10,-0.5, 9.5, "# electrons");
    
    //SM kinematics
    MSPlot[(prefix+"_MVA_SMtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_SMtop_M_"+decaystring).c_str(),30, 0,300, "inv. mass SM top ","GeV");
    MSPlot[(prefix+"_MVA_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mlb_"+decaystring).c_str(),25, 0, 500, "inv. mass l_{W}b ","GeV");
    MSPlot[(prefix+"_MVA_Wboson_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Wboson_M_"+decaystring).c_str(),100, 0, 100, "inv. mass W boson ","GeV");
    
    
    MSPlot[(prefix+"_MVA_Wlep_Charge_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Wlep_Charge_"+decaystring).c_str(),4, -2, 2, "l_{W} charge");
    MSPlot[(prefix+"_MVA_charge_asym_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_charge_asym_"+decaystring).c_str(),15, -3, 3, "l_{W} charge X |W boson #eta|");
    MSPlot[(prefix+"_MVA_TotalPt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_"+decaystring).c_str(), 50, 0, 2000, "total P_{T} ","GeV");
    MSPlot[(prefix+"_MVA_TotalHt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_"+decaystring).c_str(),20, 0, 1500, "total H_{T} ","GeV");
    MSPlot[(prefix+"_MVA_TotalInvMass_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_"+decaystring).c_str(),25, 0, 2000, "total inv. mass ","GeV");
    MSPlot[(prefix+"_MVA_TotalPt_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_jet_"+decaystring).c_str(), 50, 0, 2000, "total jet and E_{T}^{miss} p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_TotalInvMass_jet_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_jet_"+decaystring).c_str(),50, 0, 2000, "total jet and E_T^{miss} inv. mass ","GeV");
    MSPlot[(prefix+"_MVA_TotalPt_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalPt_lep_"+decaystring).c_str(), 50, 0, 2000, "total lepton p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_TotalHt_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalHt_lep_"+decaystring).c_str(),20, 0, 1000, "total lepton H_{T} ","GeV");
    MSPlot[(prefix+"_MVA_TotalInvMass_lep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_TotalInvMass_lep_"+decaystring).c_str(),20, 0, 1000, "total lepton inv. mass ","GeV");
    
    
    MSPlot[(prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str(),10, 0, 0.6, " leading jet CSVv2 disc.");
    
    // MSPlot[(prefix+"_MVA_CosTheta_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_CosTheta_"+decaystring).c_str(),25, -1, 1, "Cos(#theta *)");
    // MSPlot[(prefix+"_MVA_CosTheta_alt_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_CosTheta_alt_"+decaystring).c_str(),25, -1, 1, "Cos(#theta *)");
    
    // FCNC kinematics
    MSPlot[(prefix+"_MVA_Zboson_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_M_"+decaystring).c_str(),10, 70,110, "inv. mass Z boson ","GeV");
    
    // interplay
    
    MSPlot[(prefix+"_MVA_m3l_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_m3l_"+decaystring).c_str(),20,0, 500, "inv. mass Leptons ","GeV");
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
void InitMVAMSPlotsTopPair(string prefix, vector <int> decayChannels){
  clock_t start_sub = clock();
  
  InitMVAMSPlotsSingletop(prefix, decayChannels);
  
  bool iswzcontrol = false;
  if(prefix.find("wzcontrol")!=std::string::npos ) iswzcontrol = true; 
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    if(decayChannels[iChan] == 4 || decayChannels[iChan] == 5) continue;
    string decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    MSPlot[(prefix+"_MVA_bdiscCSVv2_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_bdiscCSVv2_jet_1_"+decaystring).c_str(),15, 0, 0.6, "2nd leading jet CSVv2");
    MSPlot[(prefix+"_MVA_cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsB_jet_0_"+decaystring).c_str(),10, -1, 1, "leading jet CvsB");
    MSPlot[(prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str(),10, -1, 1, "leading jet CvsL");
    
    
    MSPlot[(prefix+"_MVA_jet1_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_pt_"+decaystring).c_str(), 10,0, 250, "2nd leading jet p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_jet1_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_eta_"+decaystring).c_str(),15,-3, 3, "2nd leading jet #eta ");
    MSPlot[(prefix+"_MVA_jet1_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_jet1_phi_"+decaystring).c_str(),20,-4, 4, "2nd leading jet #phi ");
    
    
    MSPlot[(prefix+"_MVA_LightJet_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_pt_"+decaystring).c_str(), 10,0, 250, "light jet p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_LightJet_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_eta_"+decaystring).c_str(),60,-6, 6, "light jet #eta ");
    MSPlot[(prefix+"_MVA_LightJet_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_LightJet_phi_"+decaystring).c_str(),40,-4, 4, "light jet #phi ");
    
    MSPlot[(prefix+"_MVA_FCNCtop_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_pt_"+decaystring).c_str(), 500,0, 500, "FCNC top p_{T} ","GeV");
    MSPlot[(prefix+"_MVA_FCNCtop_eta_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_eta_"+decaystring).c_str(),60,-6, 6, "FCNC top #eta ");
    MSPlot[(prefix+"_MVA_FCNCtop_phi_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_phi_"+decaystring).c_str(),40,-4, 4, "FCNC top #phi ");
    
    MSPlot[(prefix+"_MVA_nJets_CharmL_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmL_"+decaystring).c_str(), 10,-0.5, 9.5, "# charm loose jets");
    MSPlot[(prefix+"_MVA_nJets_CharmM_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmM_"+decaystring).c_str(), 10,-0.5, 9.5, "# charm medium jets");
    MSPlot[(prefix+"_MVA_nJets_CharmT_"+decaystring).c_str()]=new MultiSamplePlot(datasets, (prefix+"_MVA_nJets_CharmT_"+decaystring).c_str(), 10,-0.5, 9.5, "# charm tight jets");
    
    
    // FCNC kinematics
    MSPlot[(prefix+"_MVA_FCNCtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_FCNCtop_M_"+decaystring).c_str(),300, 0,300, "inv. mass FCNC top ","GeV");
    
    MSPlot[(prefix+"_MVA_dRZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZc_"+decaystring).c_str(),100,-10, 10, "#Delta R(Z,q)");
    MSPlot[(prefix+"_MVA_dPhiZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiZc_"+decaystring).c_str(),40,-4, 4, "#Delta #phi (Z,q)");
    
    MSPlot[(prefix+"_MVA_cdiscCvsB_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsB_jet_1_"+decaystring).c_str(),10, -1, 1, "2nd leading jet CvsB");
    MSPlot[(prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str(),10, -1, 1, "2nd leading jet CvsL");
    
    // interplay
    MSPlot[(prefix+"_MVA_dRSMFCNCtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRSMFCNCtop_"+decaystring).c_str(),100,-10, 10, "#Delta R(SM top,FCNC top)");
    MSPlot[(prefix+"_MVA_dRWlepc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRWlepc_"+decaystring).c_str(),100,-10, 10, "#Delta R(l_{W},q)");
    
    MSPlot[(prefix+"_MVA_dPhiSMFCNCtop_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiSMFCNCtop_"+decaystring).c_str(),40,-4, 4, "#Delta #phi (SM top,FCNC top)");
    MSPlot[(prefix+"_MVA_dPhiWlepc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiWlepc_"+decaystring).c_str(),40,-4, 4, "#Delta #phi (l_{W},q)");
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
void Init2DPlots(){
  TH2::SetDefaultSumw2();
  //histo2D["CosTheta"]= new TH2F("CosTheta", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);
}
void InitGenInfoPlots(string dataSetName){
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  histo1D[("Topmass_Wb_"+dataSetName).c_str()] = new TH1F(("Topmass_Wb_"+dataSetName).c_str(),"inv. mass Wb ",500,0,500);
  histo1D[("pt_Wb_"+dataSetName).c_str()] = new TH1F(("pt_Wb_"+dataSetName).c_str(),"p_{T}(Wb) ",500,0,500);
  histo1D[("eta_Wb_"+dataSetName).c_str()]= new TH1F(("eta_Wb_"+dataSetName).c_str(),"#eta(Wb)",30,-3,3);
  histo1D[("phi_Wb_"+dataSetName).c_str()]= new TH1F(("phi_Wb_"+dataSetName).c_str(),"#phi(Wb)",32,-3.2,3.2);
  histo1D[("dPhi_Wb_"+dataSetName).c_str()]  = new TH1F(("DPhiWb_"+dataSetName).c_str(),"#Delta #Phi (W,b)", 140,-7, 7);
  histo1D[("dR_Wb_"+dataSetName).c_str()] =  new TH1F(("DRWb_"+dataSetName).c_str(),"#Delta R(W,b)", 500,0, 5);
  
  
  histo1D[("Topmass_topquark_"+dataSetName).c_str()] = new TH1F(("Topmass_topquark_"+dataSetName).c_str(),"inv. mass top ",500,0,500);
  histo1D[("pt_topquark_"+dataSetName).c_str()] = new TH1F(("pt_topquark_"+dataSetName).c_str(),"p_{T} top ",500,0,500);
  histo1D[("eta_topquark_"+dataSetName).c_str()]= new TH1F(("eta_topquark_"+dataSetName).c_str(),"#eta top",30,-3,3);
  histo1D[("phi_topquark_"+dataSetName).c_str()]= new TH1F(("phi_topquark_"+dataSetName).c_str(),"#phi top",32,-3.2,3.2);
  
  histo1D[("dR_Wbvstopquark_"+dataSetName).c_str()]  = new TH1F(("DRWbvstopquark_"+dataSetName).c_str(),"#Delta R (Wb,top)", 500,0, 5);
  histo1D[("dPhi_Wbvstopquark_"+dataSetName).c_str()]  = new TH1F(("DPhiWbvstopquark_"+dataSetName).c_str(),"#Delta #phi(Wb,top)", 140,-7, 7);
  
  histo2D[("Topmass_topquarkvsWb_"+dataSetName).c_str()] = new TH2F(("Topmass_topquarkvsWb_"+dataSetName).c_str(),"inv. mass top :Wb:t",500,0,500,500,0,500);
  histo2D[("pt_topquarkvsWb_"+dataSetName).c_str()] = new TH2F(("pt_topquarkvsWb_"+dataSetName).c_str(),"p_{T} :Wb:t",500,0,500,500,0,500);
  histo2D[("eta_topquarkvsWb_"+dataSetName).c_str()]= new TH2F(("eta_topquarkvsWb_"+dataSetName).c_str(),"#eta:Wb:t",30,-3,3,30,-3,3);
  histo2D[("phi_topquarkvsWb_"+dataSetName).c_str()]= new TH2F(("phi_topquarkvsWb_"+dataSetName).c_str(),"#phi:Wb:t",32,-3.2,3.2,32,-3.2,3.2);
  
  histo1D[("dR_Wbtop_"+dataSetName).c_str()]  = new TH1F(("DRWbtop_"+dataSetName).c_str(),"#Delta R (Wb,top)", 500,0, 5);
  histo1D[("dPhi_Wbtop_"+dataSetName).c_str()]  = new TH1F(("DPhiWbtop_"+dataSetName).c_str(),"#Delta #phi(Wb,top)", 140,-7, 7);
  
  histo1D[("mass_ZlepMin_"+dataSetName).c_str()]                                  = new TH1F(("mass_ZlepMin_"+dataSetName).c_str(),"inv. mass l_{Z}^{-} ",250,0,0.5);
  histo1D[("mass_ZlepPlus_"+dataSetName).c_str()]                                  = new TH1F(("mass_ZlepPlus_"+dataSetName).c_str(),"inv. mass l_{Z}^{+} ",250,0,0.5);
  histo1D[("Zmass_Zleptons_"+dataSetName).c_str()]                                  = new TH1F(("Zmass_Zleptons_"+dataSetName).c_str(),"inv. mass (l_{Z}^{-},l_{Z}^{+}) ",200,0,200);
  histo1D[("Zmass_Zboson_"+dataSetName).c_str()]                                  = new TH1F(("Zmass_Zboson_"+dataSetName).c_str(),"inv. mass Z boson ",200,0,200);
  histo2D[("mass_ZlepMinvsZlepPlus_"+dataSetName).c_str()]            = new TH2F(("mass_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"inv. mass leptons ;inv. mass l_{Z}^{-};inv. mass l_{Z}^{+}", 250,0,0.5,250,0,0.5 );
  histo2D[("Zmass_ZbosonvsZleptons_"+dataSetName).c_str()]            = new TH2F(("Zmass_ZbosonbsZleptons_"+dataSetName).c_str(),"inv. mass Z boson ;l_{Z}^{-},l_{Z}^{+};Z boson", 200,0,200,200,0,200 );
  
  histo1D[("pt_ZlepMin_"+dataSetName).c_str()]          = new TH1F(("pt_ZlepMin_"+dataSetName).c_str(),"l_{Z}^{-} p_{T}  ", 200,0,400);
  histo1D[("pt_ZlepPlus_"+dataSetName).c_str()]          = new TH1F(("pt_ZlepPlus_"+dataSetName).c_str(),"l_{Z}^{+} p_{T} ", 200,0,400);
  histo1D[("pt_Zboson_"+dataSetName).c_str()]          = new TH1F(("pt_Zboson_"+dataSetName).c_str(),"Z boson p_{T} ", 200,0,400);
  histo1D[("pt_Zleptons_"+dataSetName).c_str()]          = new TH1F(("pt_Zleptons_"+dataSetName).c_str(),"p_{T} Zleptons", 200,0,400);
  histo2D[("pt_ZlepMinvsZlepPlus_"+dataSetName).c_str()]          = new TH2F(("pt_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"p_{T} Z leptons ;p_{T} lep +;p_{T} lep -", 200,0,400, 200,0,400);
  histo2D[("pt_ZbosonvsZleptons_"+dataSetName).c_str()]          = new TH2F(("pt_ZbosonvsZleptons_"+dataSetName).c_str(),"p_{T} Z boson ;p_{T}(l_{Z}^{-},l_{Z}^{+});p_{T} Z boson", 200,0,400, 200,0,400);
  
  histo1D[("phi_ZlepMin_"+dataSetName).c_str()]          = new TH1F(("phi_ZlepMin_"+dataSetName).c_str(),"#phi lep +",32,-3.2,3.2);
  histo1D[("phi_ZlepPlus_"+dataSetName).c_str()]          = new TH1F(("phi_ZlepPlus_"+dataSetName).c_str(),"#phi lep -",32,-3.2,3.2);
  histo1D[("phi_Zboson_"+dataSetName).c_str()]          = new TH1F(("phi_Zboson_"+dataSetName).c_str(),"#phi Zboson",32,-3.2,3.2);
  histo1D[("phi_Zleptons_"+dataSetName).c_str()]          = new TH1F(("phi_Zleptons_"+dataSetName).c_str(),"#phi Zleptons",32,-3.2,3.2);
  histo2D[("phi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]          = new TH2F(("phi_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"#phi lep;phi lep +;phi lep -",32,-3.2,3.2,32,-3.2,3.2);
  histo2D[("phi_ZbosonvsZleptons_"+dataSetName).c_str()]          = new TH2F(("phiZbosonvsZleptons_"+dataSetName).c_str(),"#phi Z;phi Zleptons;phi Zboson",32,-3.2,3.2,32,-3.2,3.2);
  
  histo1D[("eta_ZlepMin_"+dataSetName).c_str()]          = new TH1F(("eta_ZlepMin_"+dataSetName).c_str(),"l_{Z}^{-} #eta", 30,-3,3);
  histo1D[("eta_ZlepPlus_"+dataSetName).c_str()]          = new TH1F(("eta_ZlepPlus_"+dataSetName).c_str(),"l_{Z}^{+} #eta ", 30,-3,3);
  histo1D[("eta_Zboson_"+dataSetName).c_str()]          = new TH1F(("eta_Zboson_"+dataSetName).c_str(),"#eta Z boson", 30,-3,3);
  histo1D[("eta_Zleptons_"+dataSetName).c_str()]          = new TH1F(("etaZleptons_"+dataSetName).c_str(),"#eta l_{Z} ", 30,-3,3);
  histo2D[("eta_ZlepMinvsZlepPlus_"+dataSetName).c_str()]          = new TH2F(("eta_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"#eta lep;l_{Z}^{+} #eta ;#eta l_{Z}^{-}", 30,-3,3, 30,-3,3);
  histo2D[("eta_ZbosonvsZleptons_"+dataSetName).c_str()]          = new TH2F(("etaZbosonvsZleptons_"+dataSetName).c_str(),"#eta Z;#eta (l_{Z}^{-},l_{Z}^{+});#eta Z boson", 30,-3,3, 30,-3,3);
  
  histo1D[("dR_ZlepMinvsZlepPlus_"+dataSetName).c_str()]          = new TH1F(("dR_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"Distance between the Z boson leptons; #Delta R(l^{+}_{Z}l^{-}_{Z}); # Events", 100,0, 5);
  histo1D[("dPhi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]          = new TH1F(("dPhi_ZlepMinvsZlepPlus_"+dataSetName).c_str(),"#Delta #phi (l_{Z}^{-},l_{Z}^{+})", 140,-7, 7);
  
  
  histo1D[("GenInfo_mWT_"+dataSetName).c_str()] = new TH1F(("GenInfo_mWT_"+dataSetName).c_str(),"GenInfo Transv. Mass W boson ", 100, 0, 400 );
  histo1D[("GenInfo_mWT2_"+dataSetName).c_str()] = new TH1F(("GenInfo_mWT2_"+dataSetName).c_str(),"GenInfo Transv. Mass W boson ", 50, 0, 300 );
  histo2D[("GenInfo_mWTmWT2_"+dataSetName).c_str()] = new TH2F(("GenInfo_mWTmWT2_"+dataSetName).c_str(),"GenInfo Transv. Mass W boson ", 50, 0, 300, 50,0,300);
  
  histo1D[("GenInfo_CosTheta_"+dataSetName).c_str()] = new TH1F(("GenInfo_CosTheta_"+dataSetName).c_str(),"GenInfo Cos(#theta *)",25, -1, 1);
  histo1D[("GenInfo_CosTheta_alt_"+dataSetName).c_str()] = new TH1F(("GenInfo_CosTheta_alt_"+dataSetName).c_str(), " GenInfo Cos(#theta *)",25, -1, 1);
  histo2D[("GenInfo_CosTheta2_"+dataSetName).c_str()] = new TH2F(("GenInfo_CosTheta2_"+dataSetName).c_str(),"GenInfo Cos(#theta *)",25, -1, 1,25,-1,1);
  
  
  
}
void InitRecovsGenInfoPlots(string dataSetName){
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //cout << ("Matched_mlb_"+dataSetName).c_str() << endl;
  histo1D[("Matched_mlb_"+dataSetName).c_str()] = new TH1F(("Matched_mlb_"+dataSetName).c_str(),"inv. mass lb matched ",500,0,500);
  histo1D[("Matched_Zmass_"+dataSetName).c_str()] = new TH1F(("Matched_Zmass_"+dataSetName).c_str(),"inv. mass Z boson matched ",200,0,200);
  
  
  
}

///////////////////////////////////// INIT TREES /////////////////////////////////////////
void InitTree(TTree* tree, bool isData, bool isfakes){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("x1", &x1, &b_x1);
  tree->SetBranchAddress("x2", &x2, &b_x2);
  tree->SetBranchAddress("id1", &id1, &b_id1);
  tree->SetBranchAddress("id2", &id2, &b_id2);
  tree->SetBranchAddress("q", &q, &b_q);
 
    tree->SetBranchAddress("hdamp_up", &hdamp_up, &b_hdamp_up);
    tree->SetBranchAddress("hdamp_down", &hdamp_down, &b_hdamp_down);
  if(!isData && !isfakes){
    tree->SetBranchAddress("weight0", &weight0, &b_weight0);
    tree->SetBranchAddress("weight1", &weight1, &b_weight1);
    tree->SetBranchAddress("weight2", &weight2, &b_weight2);
    tree->SetBranchAddress("weight3", &weight3, &b_weight3);
    tree->SetBranchAddress("weight4", &weight4, &b_weight4);
    tree->SetBranchAddress("weight5", &weight5, &b_weight5);
    tree->SetBranchAddress("weight6", &weight6, &b_weight6);
    tree->SetBranchAddress("weight7", &weight7, &b_weight7);
    tree->SetBranchAddress("weight8", &weight8, &b_weight8);}
  tree->SetBranchAddress("btagSFshape", &btagSFshape, &b_btagSFshape);
  tree->SetBranchAddress("btagSFshape_down_cferr1", &btagSFshape_down_cferr1, &b_btagSFshape_down_cferr1);
  tree->SetBranchAddress("btagSFshape_down_cferr2", &btagSFshape_down_cferr2, &b_btagSFshape_down_cferr2);
  tree->SetBranchAddress("btagSFshape_down_hf", &btagSFshape_down_hf, &b_btagSFshape_down_hf);
  tree->SetBranchAddress("btagSFshape_down_hfstats1", &btagSFshape_down_hfstats1, &b_btagSFshape_down_hfstats1);
  tree->SetBranchAddress("btagSFshape_down_hfstats2", &btagSFshape_down_hfstats2, &b_btagSFshape_down_hfstats2);
  tree->SetBranchAddress("btagSFshape_down_lf", &btagSFshape_down_lf, &b_btagSFshape_down_lf);
  tree->SetBranchAddress("btagSFshape_down_lfstats1", &btagSFshape_down_lfstats1, &b_btagSFshape_down_lfstats1);
  tree->SetBranchAddress("btagSFshape_down_lfstats2", &btagSFshape_down_lfstats2, &b_btagSFshape_down_lfstats2);
  tree->SetBranchAddress("btagSFshape_up_cferr1", &btagSFshape_up_cferr1, &b_btagSFshape_up_cferr1);
  tree->SetBranchAddress("btagSFshape_up_cferr2", &btagSFshape_up_cferr2, &b_btagSFshape_up_cferr2);
  tree->SetBranchAddress("btagSFshape_up_hf", &btagSFshape_up_hf, &b_btagSFshape_up_hf);
  tree->SetBranchAddress("btagSFshape_up_hfstats1", &btagSFshape_up_hfstats1, &b_btagSFshape_up_hfstats1);
  tree->SetBranchAddress("btagSFshape_up_hfstats2", &btagSFshape_up_hfstats2, &b_btagSFshape_up_hfstats2);
  tree->SetBranchAddress("btagSFshape_up_lf", &btagSFshape_up_lf, &b_btagSFshape_up_lf);
  tree->SetBranchAddress("btagSFshape_up_lfstats1", &btagSFshape_up_lfstats1, &b_btagSFshape_up_lfstats1);
  tree->SetBranchAddress("btagSFshape_up_lfstats2", &btagSFshape_up_lfstats2, &b_btagSFshape_up_lfstats2);
  tree->SetBranchAddress("channelInt", &channelInt, &b_channelInt);
  tree->SetBranchAddress("nloWeight", &nloWeight, &b_nloWeight);
  tree->SetBranchAddress("run_num", &run_num, &b_run_num);
  tree->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
  tree->SetBranchAddress("lumi_num", &lumi_num, &b_lumi_num);
  tree->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
  tree->SetBranchAddress("npu", &npu, &b_npu);
  tree->SetBranchAddress("puSF", &puSF, &b_puSF);
  tree->SetBranchAddress("puSF_up", &puSF_up, &b_puSF_up);
  tree->SetBranchAddress("puSF_down", &puSF_down, &b_puSF_down);
  tree->SetBranchAddress("PassedMETFilter", &PassedMETFilter, &b_PassedMETFilter);
  tree->SetBranchAddress("PassedTriggerMET", &PassedTriggerMET, &b_PassedTriggerMET);
  tree->SetBranchAddress("PassedTrigger", &PassedTrigger, &b_PassedTrigger);
  tree->SetBranchAddress("PassedTriggerJET", &PassedTriggerJET, &b_PassedTriggerJET);
  tree->SetBranchAddress("PassedTriggerNoLogic", &PassedTriggerNoLogic, &b_PassedTriggerNoLogic);
  tree->SetBranchAddress("PassedTriggerNoLogic2", &PassedTriggerNoLogic2, &b_PassedTriggerNoLogic2);
  tree->SetBranchAddress("PassedGoodPV", &PassedGoodPV, &b_PassedGoodPV);
  tree->SetBranchAddress("nbOfLooseElectrons", &nbOfLooseElectrons, &b_nbOfLooseElectrons);
  tree->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
  tree->SetBranchAddress("ElectronSF", ElectronSF, &b_ElectronSF);
  tree->SetBranchAddress("ElectronSF_up", ElectronSF_up, &b_ElectronSF_up);
  tree->SetBranchAddress("ElectronSF_down", ElectronSF_down, &b_ElectronSF_down);
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
  tree->SetBranchAddress("nbOfLooseMuons", &nbOfLooseMuons, &b_nbOfLooseMuons);
  tree->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
  tree->SetBranchAddress("MuonIDSF", MuonIDSF, &b_MuonIDSF);
  tree->SetBranchAddress("ptSF_muon",ptSF_muon,&b_ptSF_muon);
  tree->SetBranchAddress("MuonTrackSF", MuonTrackSF, &b_MuonTrackSF);
  tree->SetBranchAddress("MuonIsoSF", MuonIsoSF, &b_MuonIsoSF);
  tree->SetBranchAddress("MuonIDSF_up", MuonIDSF_up, &b_MuonIDSF_up);
  tree->SetBranchAddress("MuonIsoSF_up", MuonIsoSF_up, &b_MuonIsoSF_up);
  tree->SetBranchAddress("MuonIDSF_down", MuonIDSF_down, &b_MuonIDSF_down);
  tree->SetBranchAddress("MuonIsoSF_down", MuonIsoSF_down, &b_MuonIsoSF_down);
  tree->SetBranchAddress("MuonTrigSFv2", MuonTrigSFv2, &b_MuonTrigSFv2);
  
  //tree->SetBranchAddress("rejecteventBadPFmuon",rejecteventBadPFmuon, &b_rejecteventBadPFmuon);
  tree->SetBranchAddress("MuonTrigSFv3", MuonTrigSFv3, &b_MuonTrigSFv3);
  tree->SetBranchAddress("badmueventmu", badmueventmu, &b_badmueventmu);
  tree->SetBranchAddress("badmueventclonemu", badmueventclonemu, &b_badmueventclonemu);
  tree->SetBranchAddress("pt_muon", pt_muon, &b_pt_muon);
  tree->SetBranchAddress("TrackLayers_muon",TrackLayers_muon, &b_TrackLayers_muon);
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
  if(!isData && !isfakes){
    tree->SetBranchAddress("nMCParticles", &nMCParticles, &b_nMCParticles);
    tree->SetBranchAddress("mc_status", mc_status, &b_mc_status);
    tree->SetBranchAddress("mc_pdgId", mc_pdgId, &b_mc_pdgId);
    tree->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
    tree->SetBranchAddress("mc_granny", mc_granny, &b_mc_granny);
    tree->SetBranchAddress("mc_pt", mc_pt, &b_mc_pt);
    tree->SetBranchAddress("mc_phi", mc_phi, &b_mc_phi);
    tree->SetBranchAddress("mc_eta", mc_eta, &b_mc_eta);
    tree->SetBranchAddress("mc_E", mc_E, &b_mc_E);
    tree->SetBranchAddress("mc_M", mc_M, &b_mc_M);
    tree->SetBranchAddress("mc_isLastCopy", mc_isLastCopy, &b_mc_isLastCopy);
    tree->SetBranchAddress("mc_isPromptFinalState", mc_isPromptFinalState, &b_mc_isPromptFinalState);
    tree->SetBranchAddress("mc_isHardProcess", mc_isHardProcess, &b_mc_isHardProcess);
    tree->SetBranchAddress("mc_fromHardProcessFinalState", mc_fromHardProcessFinalState, &b_mc_fromHardProcessFinalState);
  }
  tree->SetBranchAddress("met_Pt", &met_Pt, &b_met_Pt);
  tree->SetBranchAddress("met_Eta", &met_Eta, &b_met_Eta);
  tree->SetBranchAddress("met_Phi", &met_Phi, &b_met_Phi);
  tree->SetBranchAddress("met_Px", &met_Px, &b_met_Px);
  tree->SetBranchAddress("met_Py", &met_Py, &b_met_Py);
  tree->SetBranchAddress("met_before_JES", &met_before_JES, &b_met_before_JES);
  tree->SetBranchAddress("met_after_JES", &met_after_JES, &b_met_after_JES)  ;
  
}




void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData, bool isfakes, bool threelepregion,bool twolepregion){
  // cout << "fill plots" << endl;
  string decaystring = "";
  Double_t eventW = 1.;
  // if(isData  ) scaleFactor = 1.;
  eventW = Luminosity/EquilumiSF;
  //if(isfakes) eventW *= 0.0001;
  
  vector<string> v_prefixregion = {"2lep", "3lep"};
  if(!doDilep) v_prefixregion = {"3lep"};
  string prefixregion = "";
  for(int iReg = 0; iReg < v_prefixregion.size(); iReg++){
    prefixregion = v_prefixregion[iReg];
    if(prefixregion.find("3lep")!=std::string::npos && !threelepregion) continue;
    if(prefixregion.find("2lep")!=std::string::npos && !twolepregion) continue;
    for(int iChan =0; iChan < decayChannels.size() ; iChan++){
      decaystring = "";
      //  cout << decayChannels[iChan] << " " << channelInt << " " <<  threelepregion << " " << twolepregion << endl;
      // if(decayChannels[iChan] == -9) continue;;
      //cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
      if(decayChannels[iChan] != channelInt && decayChannels[iChan] != -9) continue;
      //if(decayChannels[iChan] == -9)cout << decayChannels[iChan] << " " << channelInt << " passed " << endl;
      if(decayChannels[iChan] == 0) decaystring = "uuu";
      if(decayChannels[iChan] == 1) decaystring = "uue";
      if(decayChannels[iChan] == 2) decaystring = "eeu";
      if(decayChannels[iChan] == 3) decaystring = "eee";
      if(decayChannels[iChan] == 4) decaystring = "uu";
      if(decayChannels[iChan] == 5) decaystring = "ee";
      if(decayChannels[iChan] == -9) decaystring = "all";
      
      if((decayChannels[iChan] == 4 || decayChannels[iChan] == 5) && prefixregion.find("3lep")!=std::string::npos) continue;
      //cout << "nvtx "  << nvtx << endl;
      MSPlot[(prefixregion+prefix+"_NbOfVertices_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true, eventW*scaleFactor_bfPU);
      MSPlot[(prefixregion+prefix+"_puSF_"+decaystring).c_str()]->Fill(puSF , datasets[d], true, 1);
      
      MSPlot[(prefixregion+prefix+"_btagSF_"+decaystring).c_str()]->Fill(btagSFshape , datasets[d], true, 1);
      
      //cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
      for(int iterJ = 0 ; iterJ < selectedJets.size(); iterJ++){
        MSPlot[(prefixregion+prefix+"_bdisc_"+decaystring).c_str()]->Fill(bdisc_jet[iterJ] , datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_bdisc_bfBT_"+decaystring).c_str()]->Fill(bdisc_jet[iterJ] , datasets[d], true,eventW*scaleFactor_bfBT);
        
        MSPlot[(prefixregion+prefix+"_cvsldisc_"+decaystring).c_str()]->Fill(cdiscCvsL_jet[iterJ], datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_cvsbdisc_"+decaystring).c_str()]->Fill(cdiscCvsB_jet[iterJ], datasets[d], true,eventW*scaleFactor);
      }
      MSPlot[(prefixregion+prefix+"_nMu_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nMu_bfMuSF_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true, eventW*scaleFactor_bfMuSF);
      MSPlot[(prefixregion+prefix+"_muSF_"+decaystring).c_str()]->Fill(muonSFtemp , datasets[d], true, 1);
      
      MSPlot[(prefixregion+prefix+"_nEl_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nEl_bfElSF_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true, eventW*scaleFactor_bfELSF);
      MSPlot[(prefixregion+prefix+"_elSF_"+decaystring).c_str()]->Fill(electronSFtemp , datasets[d], true, 1);
      
      if(selectedJets.size() >0){
        // cout << (prefix+"_JetPt_bfJER_"+decaystring).c_str() << endl;
        MSPlot[(prefixregion+prefix+"_JetPt_bfJER_"+decaystring).c_str()]->Fill(jet_Pt_before_JER[0] , datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_JetPt_bfJES_"+decaystring).c_str()]->Fill(jet_Pt_before_JES[0] , datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_JetPt_afJER_"+decaystring).c_str()]->Fill(jet_Pt_after_JER[0] , datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_JetPt_afJES_"+decaystring).c_str()]->Fill(jet_Pt_after_JES[0] , datasets[d], true,eventW*scaleFactor);
      }
      MSPlot[(prefixregion+prefix+"_met_bfJES_"+decaystring).c_str()]->Fill(met_before_JES , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_afJES_"+decaystring).c_str()]->Fill(met_after_JES , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_px_"+decaystring).c_str()]->Fill(met_Px , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_py_"+decaystring).c_str()]->Fill(met_Py , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_phi_"+decaystring).c_str()]->Fill(met_Phi , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_pz_"+decaystring).c_str()]->Fill(met_Pz , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_met_pt_"+decaystring).c_str()]->Fill(met_Pt , datasets[d], true,eventW*scaleFactor);
      
      // vars
      
      MSPlot[(prefixregion+prefix+"_ZbosonMass_"+decaystring).c_str()]->Fill(Zboson.M() , datasets[d], true,eventW*scaleFactor);
      if(prefixregion.find("3lep")!=std::string::npos){
        MSPlot[(prefixregion+prefix+"_WbosonMass_"+decaystring).c_str()] ->Fill(Wboson.M() , datasets[d], true,eventW*scaleFactor);
        
        MSPlot[(prefixregion+prefix+"_mWT_"+decaystring).c_str()]->Fill(mWT, datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_mWT2_"+decaystring).c_str()]->Fill(mWT2, datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_SMTopMass_"+decaystring).c_str()]->Fill( SMtop.M() , datasets[d], true,eventW*scaleFactor);
        MSPlot[(prefixregion+prefix+"_mlb_"+decaystring).c_str()]->Fill((SMbjet+Wlep).M() , datasets[d], true,eventW*scaleFactor);
      }
      
      if(selectedJets.size()>0) MSPlot[(prefixregion+prefix+"_LeadingJetPt_"+decaystring).c_str()]->Fill( selectedJets[0].Pt() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_LeadingLepPt_"+decaystring).c_str()]->Fill(selectedLeptons[0].Pt() , datasets[d], true,eventW*scaleFactor);
      if(selectedJets.size()>1) MSPlot[(prefixregion+prefix+"_2ndLeadingJetPt_"+decaystring).c_str()]->Fill(selectedJets[1].Pt() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_2ndLeadingLepPt_"+decaystring).c_str()]->Fill(selectedLeptons[1].Pt() , datasets[d], true,eventW*scaleFactor);
      
      MSPlot[(prefixregion+prefix+"_nJets_"+decaystring).c_str()]->Fill(selectedJets.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCSVL_"+decaystring).c_str()]->Fill(selectedCSVLJetID.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCSVM_"+decaystring).c_str()]->Fill(selectedCSVMJetID.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCSVT_"+decaystring).c_str()]->Fill(selectedCSVTJetID.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCharmL_"+decaystring).c_str()]->Fill(selectedCharmLJetsindex.size() , datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCharmM_"+decaystring).c_str()]->Fill(selectedCharmMJetsindex.size(), datasets[d], true,eventW*scaleFactor);
      MSPlot[(prefixregion+prefix+"_nJetsCharmT_"+decaystring).c_str()]->Fill(selectedCharmTJetsindex.size(), datasets[d], true,eventW*scaleFactor);
      
      
      
      if(decayChannels[iChan] != 0 && decayChannels[iChan]!=4 ){ // not uuu or uu
        if(selectedElectrons.size()>2 && (decayChannels[iChan]==3 || decayChannels[iChan] ==-9)  ){
          MSPlot[(prefixregion+prefix+"_3dLeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[2].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_3dLeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[2].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_3dLeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[2]] , datasets[d], true,eventW*scaleFactor);
          //cout << "_3dLeadingElEta_" << eta_electron[electronID[2]] << endl;
          MSPlot[(prefixregion+prefix+"_3dLeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[2]] , datasets[d], true,eventW*scaleFactor);
          
          if(WelecIndiceF != -999){
            MSPlot[(prefixregion+prefix+"_WlepElEta_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Eta(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Phi(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepElPt_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Pt(), datasets[d], true,eventW*scaleFactor);
          }
          
        }
        if(selectedElectrons.size()>1 && (decayChannels[iChan]==3 ||decayChannels[iChan]==2 || decayChannels[iChan]==5 || decayChannels[iChan] ==-9) ){
          //cout << "filling " << (prefixregion+prefix+"_2ndLeadingElEta_"+decaystring).c_str() << endl;
          MSPlot[(prefixregion+prefix+"_2ndLeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[1].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[1].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[1]], datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[1]] , datasets[d], true,eventW*scaleFactor);
          if( ZelecIndiceF_1!= -999 && ZelecIndiceF_0 != -999  ){
            MSPlot[(prefixregion+prefix+"_ZbosonElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[ZelecIndiceF_0]] + pfIso_electron[electronID[ZelecIndiceF_1]]  , datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_ZbosonEldPhi_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(selectedElectrons[ZelecIndiceF_0],selectedElectrons[ZelecIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonEldR_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(selectedElectrons[ZelecIndiceF_0],selectedElectrons[ZelecIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonPtEl_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonMassEl_"+decaystring).c_str()]->Fill(Zboson.M(),datasets[d], true,eventW*scaleFactor) ;
            
          }
          
          if(selectedMuons.size()>0 && (decayChannels[iChan]==2 || decayChannels[iChan] ==-9)){
            if(WmuIndiceF != -999){
              MSPlot[(prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Eta(), datasets[d], true,eventW*scaleFactor);
              MSPlot[(prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Phi(), datasets[d], true,eventW*scaleFactor);
              MSPlot[(prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Pt(), datasets[d], true,eventW*scaleFactor);
            }
          }
        }
        if(selectedElectrons.size()>0 && (decayChannels[iChan]==3 ||decayChannels[iChan]==2 ||decayChannels[iChan]==1 || decayChannels[iChan]==5 || decayChannels[iChan] ==-9)){
          MSPlot[(prefixregion+prefix+"_LeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[0].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[0].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[0]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[0]] , datasets[d], true,eventW*scaleFactor);
          
        }
        for(int iter = 0; iter < selectedElectrons.size(); iter++){
          MSPlot[(prefixregion+prefix+"_ElPt_"+decaystring).c_str()]->Fill(selectedElectrons[iter].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_ElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[iter]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_ElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[iter].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_ElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[iter]] , datasets[d], true,eventW*scaleFactor);
        }
        
        
      }
      if(decayChannels[iChan] != 3 && decayChannels[iChan]!= 5){ // not eee or  ee
        if(selectedMuons.size()>2 && (decayChannels[iChan] == 0 || decayChannels[iChan] ==-9)){
          MSPlot[(prefixregion+prefix+"_3dLeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[2].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_3dLeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[2]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_3dLeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[2].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_3dLeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[2].Eta() , datasets[d], true,eventW*scaleFactor);
          
          
          if(WmuIndiceF != -999){
            MSPlot[(prefixregion+prefix+"_WlepMuEta_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Eta(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Phi(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepMuPt_"+decaystring).c_str()]->Fill(selectedMuons[WmuIndiceF].Pt(), datasets[d], true,eventW*scaleFactor);
          }
          
        }
        if(selectedMuons.size()>1 && (decayChannels[iChan] == 0 || decayChannels[iChan] == 1 || decayChannels[iChan]==4 || decayChannels[iChan] ==-9)){
          MSPlot[(prefixregion+prefix+"_2ndLeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[1].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[1]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[1].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_2ndLeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[1].Eta() , datasets[d], true,eventW*scaleFactor);
          
          if( ZmuIndiceF_1!= -999 && ZmuIndiceF_0 != -999  ){
            MSPlot[(prefixregion+prefix+"_ZbosonMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[ZmuIndiceF_0]]+pfIso_muon[muonID[ZmuIndiceF_1]] , datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_ZbosonMudPhi_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(selectedMuons[ZmuIndiceF_0],selectedMuons[ZmuIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonMudR_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(selectedMuons[ZmuIndiceF_0],selectedMuons[ZmuIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonPtMu_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
            MSPlot[(prefixregion+prefix+"_ZbosonMassMu_"+decaystring).c_str()]->Fill(Zboson.M(),datasets[d], true,eventW*scaleFactor) ;
          }
          if(WelecIndiceF != -999 && selectedElectrons.size() > 0 && (decayChannels[iChan]==1 || decayChannels[iChan] ==-9)){
            MSPlot[(prefixregion+prefix+"_WlepElEta_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Eta(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Phi(), datasets[d], true,eventW*scaleFactor);
            MSPlot[(prefixregion+prefix+"_WlepElPt_"+decaystring).c_str()]->Fill(selectedElectrons[WelecIndiceF].Pt(), datasets[d], true,eventW*scaleFactor);
          }
        }
        if(selectedMuons.size()>0 && ( decayChannels[iChan] == 0 || decayChannels[iChan] == 1 || decayChannels[iChan] == 2 || decayChannels[iChan]==4 || decayChannels[iChan] ==-9)){
          MSPlot[(prefixregion+prefix+"_LeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[0].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[0]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[0].Pt() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_LeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[0].Eta() , datasets[d], true,eventW*scaleFactor);
          
        }
        for(int iter = 0; iter < selectedMuons.size(); iter++){
          MSPlot[(prefixregion+prefix+"_MuPhi_"+decaystring).c_str()]->Fill(selectedMuons[iter].Phi() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_MuEta_"+decaystring).c_str()]->Fill(selectedMuons[iter].Eta() , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_MuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[iter]] , datasets[d], true,eventW*scaleFactor);
          MSPlot[(prefixregion+prefix+"_MuPt_"+decaystring).c_str()]->Fill(selectedMuons[iter].Pt() , datasets[d], true,eventW*scaleFactor);
        }
        
        
      }
      
      MSPlot[(prefixregion+prefix+"_ZbosonPt_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
      
      
    }
    /*
     for(int iChan =0; iChan < decayChannels.size() ; iChan++){
     // cout << decayChannels[iChan] << endl;
     if(decayChannels[iChan] != -9) continue;
     //cout << "fill all" << endl;
     decaystring = "all";
     
     
     //cout << "nvtx "  << nvtx << endl;
     MSPlot[(prefix+"_NbOfVertices_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_NbOfVertices_bfPU_"+decaystring).c_str()]->Fill(nvtx , datasets[d], true, eventW*scaleFactor_bfPU);
     MSPlot[(prefix+"_puSF_"+decaystring).c_str()]->Fill(puSF , datasets[d], true, 1);
     
     //cout << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
     MSPlot[(prefix+"_bdisc_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_bdisc_bfBT_"+decaystring).c_str()]->Fill(bdisc_jet[0] , datasets[d], true, eventW*scaleFactor_bfBT);
     MSPlot[(prefix+"_btagSF_"+decaystring).c_str()]->Fill(btagSFshape , datasets[d], true, 1);
     
     MSPlot[(prefix+"_nMu_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_nMu_bfMuSF_"+decaystring).c_str()]->Fill(selectedMuons.size() , datasets[d], true,eventW*scaleFactor_bfMuSF);
     MSPlot[(prefix+"_muSF_"+decaystring).c_str()]->Fill(muonSFtemp , datasets[d], true, 1);
     
     MSPlot[(prefix+"_nEl_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_nEl_bfElSF_"+decaystring).c_str()]->Fill(selectedElectrons.size() , datasets[d], true, eventW*scaleFactor_bfELSF);
     MSPlot[(prefix+"_elSF_"+decaystring).c_str()]->Fill(electronSFtemp , datasets[d], true, 1);
     
     if(selectedJets.size() >0){
     // cout << (prefix+"_JetPt_bfJER_"+decaystring).c_str() << endl;
     MSPlot[(prefix+"_JetPt_bfJER_"+decaystring).c_str()]->Fill(jet_Pt_before_JER[0] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_JetPt_bfJES_"+decaystring).c_str()]->Fill(jet_Pt_before_JES[0] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_JetPt_afJER_"+decaystring).c_str()]->Fill(jet_Pt_after_JER[0] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_JetPt_afJES_"+decaystring).c_str()]->Fill(jet_Pt_after_JES[0] , datasets[d], true,eventW*scaleFactor);
     }
     MSPlot[(prefix+"_met_bfJES_"+decaystring).c_str()]->Fill(met_before_JES , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_met_afJES_"+decaystring).c_str()]->Fill(met_after_JES , datasets[d], true,eventW*scaleFactor);
     
     
     // vars
     
     MSPlot[(prefix+"_ZbosonMass_"+decaystring).c_str()]->Fill( Zboson.M() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_WbosonMass_"+decaystring).c_str()] ->Fill(Wboson.M() , datasets[d], true,eventW*scaleFactor);
     if((selectedElectrons.size() +selectedMuons.size()) > 2 ){
     MSPlot[(prefix+"_mWT_"+decaystring).c_str()]->Fill(mWT, datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_mWT2_"+decaystring).c_str()]->Fill(mWT2, datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_SMTopMass_"+decaystring).c_str()]->Fill( SMtop.M() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_mlb_"+decaystring).c_str()]->Fill((SMbjet+Wlep).M() , datasets[d], true,eventW*scaleFactor);
     }
     
     if(selectedJets.size()>0) MSPlot[(prefix+"_LeadingJetPt_"+decaystring).c_str()]->Fill( selectedJets[0].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingLepPt_"+decaystring).c_str()]->Fill( selectedLeptons[0].Pt() , datasets[d], true,eventW*scaleFactor);
     if(selectedJets.size()>1) MSPlot[(prefix+"_2ndLeadingJetPt_"+decaystring).c_str()]->Fill(selectedJets[1].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingLepPt_"+decaystring).c_str()]->Fill(selectedLeptons[1].Pt() , datasets[d], true,eventW*scaleFactor);
     
     MSPlot[(prefix+"_nJets_"+decaystring).c_str()]->Fill(selectedJets.size() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_nJetsCSVL_"+decaystring).c_str()]->Fill(selectedCSVLJetID.size() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_nJetsCSVM_"+decaystring).c_str()]->Fill(selectedCSVMJetID.size() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_nJetsCSVT_"+decaystring).c_str()]->Fill(selectedCSVTJetID.size() , datasets[d], true,eventW*scaleFactor);
     
     
     
     if(selectedElectrons.size()>2  ){
     MSPlot[(prefix+"_3dLeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[2].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[2].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[2]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[2]] , datasets[d], true,eventW*scaleFactor);
     
     }
     if(selectedElectrons.size()>1){
     MSPlot[(prefix+"_2ndLeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[1].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[1].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[1]], datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[1]] , datasets[d], true,eventW*scaleFactor);
     
     if( ZelecIndiceF_1!= -999 && ZelecIndiceF_0 != -999  ){
     MSPlot[(prefix+"_ZbosonElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[ZelecIndiceF_0]] + pfIso_electron[electronID[ZelecIndiceF_1]]  , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_ZbosonEldPhi_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(selectedElectrons[ZelecIndiceF_0],selectedElectrons[ZelecIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
     MSPlot[(prefix+"_ZbosonEldR_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(selectedElectrons[ZelecIndiceF_0],selectedElectrons[ZelecIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
     MSPlot[(prefix+"_ZbosonPtEl_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
     
     }
     
     
     }
     if(selectedElectrons.size()>0 ){
     MSPlot[(prefix+"_LeadingElPt_"+decaystring).c_str()]->Fill(selectedElectrons[0].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[0].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[0]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[0]] , datasets[d], true,eventW*scaleFactor);
     }
     
     for(int iter = 0; iter < selectedElectrons.size(); iter++){
     MSPlot[(prefix+"_ElPt_"+decaystring).c_str()]->Fill(selectedElectrons[iter].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_ElEta_"+decaystring).c_str()]->Fill(eta_electron[electronID[iter]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_ElPhi_"+decaystring).c_str()]->Fill(selectedElectrons[iter].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_ElIso_"+decaystring).c_str()]->Fill(pfIso_electron[electronID[iter]] , datasets[d], true,eventW*scaleFactor);
     }
     
     
     
     if(selectedMuons.size()>2 ){
     MSPlot[(prefix+"_3dLeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[2].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[2]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[2].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_3dLeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[2].Eta() , datasets[d], true,eventW*scaleFactor);
     }
     if(selectedMuons.size()>1 ){
     MSPlot[(prefix+"_2ndLeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[1].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[1]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[1].Pt() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_2ndLeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[1].Eta() , datasets[d], true,eventW*scaleFactor);
     
     if( ZmuIndiceF_1!= -999 && ZmuIndiceF_0 != -999  ){
     MSPlot[(prefix+"_ZbosonMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[ZmuIndiceF_0]]+pfIso_muon[muonID[ZmuIndiceF_1]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_ZbosonMudPhi_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(selectedMuons[ZmuIndiceF_0],selectedMuons[ZmuIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
     MSPlot[(prefix+"_ZbosonMudR_"+decaystring).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(selectedMuons[ZmuIndiceF_0],selectedMuons[ZmuIndiceF_1]),datasets[d], true,eventW*scaleFactor) ;
     MSPlot[(prefix+"_ZbosonPtMu_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
     }
     }
     if(selectedMuons.size()>0 ){
     MSPlot[(prefix+"_LeadingMuEta_"+decaystring).c_str()]->Fill(selectedMuons[0].Eta() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingMuPhi_"+decaystring).c_str()]->Fill(selectedMuons[0].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingMuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[0]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_LeadingMuPt_"+decaystring).c_str()]->Fill(selectedMuons[0].Pt() , datasets[d], true,eventW*scaleFactor);
     
     }
     for(int iter = 0; iter < selectedMuons.size(); iter++){
     MSPlot[(prefix+"_MuPhi_"+decaystring).c_str()]->Fill(selectedMuons[iter].Phi() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_MuEta_"+decaystring).c_str()]->Fill(selectedMuons[iter].Eta() , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_MuIso_"+decaystring).c_str()]->Fill(pfIso_muon[muonID[iter]] , datasets[d], true,eventW*scaleFactor);
     MSPlot[(prefix+"_MuPt_"+decaystring).c_str()]->Fill(selectedMuons[iter].Pt() , datasets[d], true,eventW*scaleFactor);
     }
     
     
     MSPlot[(prefix+"_ZbosonPt_"+decaystring).c_str()]->Fill(Zboson.Pt(),datasets[d], true,eventW*scaleFactor) ;
     
     }*/
  }
  MSPlot[(prefix+"_Decay").c_str()]->Fill(channelInt , datasets[d], true,eventW*scaleFactor);
  
  
  
  
  
  //cout << "end plot filling" << endl;
}
void FillFakeValidation(string dataSetName, vector <int> decayChannels, bool isData, bool isfakes, bool threelepregion,bool twolepregion){
  // cout << "fill plots" << endl;
  string decaystr= "";
  Double_t eventW = 1.;
  // if(isData  ) scaleFactor = 1.;
  eventW = Luminosity/EquilumiSF;
  
  vector<string> v_prefixregion = {"2lep", "3lep"};
  if(!doDilep) v_prefixregion = {"3lep"};
  string prefixregion = "";
  for(int iReg = 0; iReg < v_prefixregion.size(); iReg++){
    prefixregion = v_prefixregion[iReg];
    if(prefixregion.find("3lep")!=std::string::npos && !threelepregion) continue;
    if(prefixregion.find("2lep")!=std::string::npos && !twolepregion) continue;
    for(int iChan =0; iChan < decayChannels.size() ; iChan++){
      decaystr = prefixregion + "_";
      //  cout << decayChannels[iChan] << " " << channelInt << " " <<  threelepregion << " " << twolepregion << endl;
      // if(decayChannels[iChan] == -9) continue;;
      //cout << decayChannels[iChan] << " " << channelInt << " " << (prefix+"_bdisc_bfBT_"+decaystring).c_str() << endl;
      if(decayChannels[iChan] != channelInt && decayChannels[iChan] != -9) continue;
      //if(decayChannels[iChan] == -9)cout << decayChannels[iChan] << " " << channelInt << " passed " << endl;
      if(decayChannels[iChan] == 0) decaystr += "uuu";
      if(decayChannels[iChan] == 1) decaystr += "uue";
      if(decayChannels[iChan] == 2) decaystr += "eeu";
      if(decayChannels[iChan] == 3) decaystr += "eee";
      if(decayChannels[iChan] == 4) decaystr += "uu";
      if(decayChannels[iChan] == 5) decaystr += "ee";
      if(decayChannels[iChan] == -9) decaystr += "all";
      
      if((decayChannels[iChan] == 4 || decayChannels[iChan] == 5) && prefixregion.find("3lep")!=std::string::npos) continue;
      
      histo1D_fakevvalidation[("ZbosonPt_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Pt(), eventW*scaleFactor);
      histo1D_fakevvalidation[("ZbosonEta_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Eta(), eventW*scaleFactor);
      histo1D_fakevvalidation[("ZbosonPhi_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Phi(), eventW*scaleFactor);
      if(WelecIndiceF != -999 && selectedElectrons.size() > 0 ){
        histo1D_fakevvalidation[("WlepPt_"+decaystr+dataSetName).c_str()] ->Fill(selectedElectrons[WelecIndiceF].Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepEta_"+decaystr+dataSetName).c_str()]->Fill(selectedElectrons[WelecIndiceF].Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepPhi_"+decaystr+dataSetName).c_str()]->Fill(selectedElectrons[WelecIndiceF].Phi(), eventW*scaleFactor);
      }
      if(WmuIndiceF != -999 && selectedMuons.size() > 0 ){
        histo1D_fakevvalidation[("WlepPt_"+decaystr+dataSetName).c_str()] ->Fill(selectedMuons[WmuIndiceF].Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepEta_"+decaystr+dataSetName).c_str()]->Fill(selectedMuons[WmuIndiceF].Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepPhi_"+decaystr+dataSetName).c_str()]->Fill(selectedMuons[WmuIndiceF].Phi(), eventW*scaleFactor);
      }
      histo1D_fakevvalidation[("TrMassW_"+decaystr+dataSetName).c_str()]->Fill(mWT, eventW*scaleFactor);
      
      if(ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999 &&selectedMuons.size() > 1 ){
        histo1D_fakevvalidation[("ZbosonPtMu_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("ZbosonEtaMu_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("ZbosonPhiMu_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Phi(), eventW*scaleFactor);
      }
      if(WmuIndiceF != -999 && selectedMuons.size() > 0 ){
        histo1D_fakevvalidation[("WlepPtMu_"+decaystr+dataSetName).c_str()] ->Fill(selectedMuons[WmuIndiceF].Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepEtaMu_"+decaystr+dataSetName).c_str()]->Fill(selectedMuons[WmuIndiceF].Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepPhiMu_"+decaystr+dataSetName).c_str()]->Fill(selectedMuons[WmuIndiceF].Phi(), eventW*scaleFactor);
        histo1D_fakevvalidation[("TrMassWMu_"+decaystr+dataSetName).c_str()]->Fill(mWT, eventW*scaleFactor);
      }
      if(ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999 && selectedElectrons.size() > 1 ){
        histo1D_fakevvalidation[("ZbosonPtEl_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("ZbosonEtaEl_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("ZbosonPhiEl_"+decaystr+dataSetName).c_str()] ->Fill(Zboson.Phi(), eventW*scaleFactor);
      }
      if(WelecIndiceF != -999 && selectedElectrons.size() > 0 ){
        histo1D_fakevvalidation[("WlepPtEl_"+decaystr+dataSetName).c_str()] ->Fill(selectedElectrons[WelecIndiceF].Pt(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepEtaEl_"+decaystr+dataSetName).c_str()]->Fill(selectedElectrons[WelecIndiceF].Eta(), eventW*scaleFactor);
        histo1D_fakevvalidation[("WlepPhiEl_"+decaystr+dataSetName).c_str()]->Fill(selectedElectrons[WelecIndiceF].Phi(), eventW*scaleFactor);
        histo1D_fakevvalidation[("TrMassWEl_"+decaystr+dataSetName).c_str()]->Fill(mWT, eventW*scaleFactor);
      }
      
      
    }
  }
}
void Fill1DPlots(string dataSetName, double eventW, bool twolepregion, bool threelepregion){
  //Double_t eventW = 1.;
  //eventW = Luminosity/EquilumiSF;
  string prefix = "";
  //if(twolepregion) prefix = "2lep_";
  //if(threelepregion) prefix = "3lep_";
  
  // cout << "muonSF nom : " << eventW*scaleFactor << " up " << eventW*scaleFactor_muonSF_up << " down " << eventW*scaleFactor_puSF_down << endl;
  
  histo1D_PUSystematics[(prefix+"NbVertices_nom_"+dataSetName).c_str()]->Fill(nvtx, eventW*scaleFactor);
  if(selectedElectrons.size()>0) histo1D_ElSystematics[(prefix+"Pt_electron_nom_"+dataSetName).c_str()]->Fill(selectedElectrons[0].Pt(), eventW*scaleFactor);
  if(selectedMuons.size()>0) histo1D_MuSystematics[(prefix+"Pt_muon_nom_"+dataSetName).c_str()] ->Fill(selectedMuons[0].Pt(), eventW*scaleFactor);
  histo1D_Bcferr1Systematics[(prefix+"Bdis_cferr1_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_Bcferr2Systematics[(prefix+"Bdis_cferr2_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_Bhfstats1Systematics[(prefix+"Bdis_hfstats1_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_Bhfstats2Systematics[(prefix+"Bdis_hfstats2_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_Blfstats1Systematics[(prefix+"Bdis_lfstats1_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_Blfstats2Systematics[(prefix+"Bdis_lfstats2_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_BhfSystematics[(prefix+"Bdis_hf_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  histo1D_BlfSystematics[(prefix+"Bdis_lf_nom_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor);
  
  
  
  
  histo1D_PUSystematics[(prefix+"NbVertices_up_"+dataSetName).c_str()]->Fill(nvtx, eventW*scaleFactor*scaleFactor_puSF_up/scaleFactor_puSF);
  if(selectedElectrons.size()>0) histo1D_ElSystematics[(prefix+"Pt_electron_up_"+dataSetName).c_str()]->Fill(selectedElectrons[0].Pt(), eventW*scaleFactor_electronSF_up*scaleFactor/scaleFactor_electronSF);
  if(selectedMuons.size()>0) histo1D_MuSystematics[(prefix+"Pt_muon_up_"+dataSetName).c_str()]->Fill(selectedMuons[0].Pt(), eventW*scaleFactor_muonSF_up*scaleFactor/scaleFactor_muonSF);
  histo1D_Bcferr1Systematics[(prefix+"Bdis_cferr1_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_cferr1_up*scaleFactor/scaleFactor_btagSF);
  histo1D_Bcferr2Systematics[(prefix+"Bdis_cferr2_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_cferr2_up*scaleFactor/scaleFactor_btagSF);
  histo1D_Bhfstats1Systematics[(prefix+"Bdis_hfstats1_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hfstats1_up*scaleFactor/scaleFactor_btagSF);
  histo1D_Bhfstats2Systematics[(prefix+"Bdis_hfstats2_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hfstats2_up*scaleFactor/scaleFactor_btagSF);
  histo1D_Blfstats1Systematics[(prefix+"Bdis_lfstats1_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lfstats1_up*scaleFactor/scaleFactor_btagSF);
  histo1D_Blfstats2Systematics[(prefix+"Bdis_lfstats2_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lfstats2_up*scaleFactor/scaleFactor_btagSF);
  histo1D_BhfSystematics[(prefix+"Bdis_hf_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hf_up*scaleFactor/scaleFactor_btagSF);
  histo1D_BlfSystematics[(prefix+"Bdis_lf_up_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lf_up*scaleFactor/scaleFactor_btagSF);
  
  
  histo1D_PUSystematics[(prefix+"NbVertices_down_"+dataSetName).c_str()]->Fill(nvtx, eventW*scaleFactor*scaleFactor_puSF_down/scaleFactor_puSF);
  if(selectedElectrons.size()>0) histo1D_ElSystematics[(prefix+"Pt_electron_down_"+dataSetName).c_str()]->Fill(selectedElectrons[0].Pt(), eventW*scaleFactor_electronSF_down*scaleFactor/scaleFactor_electronSF);
  if(selectedMuons.size()>0) histo1D_MuSystematics[(prefix+"Pt_muon_down_"+dataSetName).c_str()]->Fill(selectedMuons[0].Pt(), eventW*scaleFactor_muonSF_down*scaleFactor/scaleFactor_muonSF);
  histo1D_Bcferr1Systematics[(prefix+"Bdis_cferr1_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_cferr1_down*scaleFactor/scaleFactor_btagSF);
  histo1D_Bcferr2Systematics[(prefix+"Bdis_cferr2_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_cferr2_down*scaleFactor/scaleFactor_btagSF);
  histo1D_Bhfstats1Systematics[(prefix+"Bdis_hfstats1_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hfstats1_down*scaleFactor/scaleFactor_btagSF);
  histo1D_Bhfstats2Systematics[(prefix+"Bdis_hfstats2_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hfstats2_down*scaleFactor/scaleFactor_btagSF);
  histo1D_Blfstats1Systematics[(prefix+"Bdis_lfstats1_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lfstats1_down*scaleFactor/scaleFactor_btagSF);
  histo1D_Blfstats2Systematics[(prefix+"Bdis_lfstats2_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lfstats2_down*scaleFactor/scaleFactor_btagSF);
  histo1D_BhfSystematics[(prefix+"Bdis_hf_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_hf_down*scaleFactor/scaleFactor_btagSF);
  histo1D_BlfSystematics[(prefix+"Bdis_lf_down_"+dataSetName).c_str()]->Fill(bdisc_jet[0],eventW*scaleFactor_btagSF_lf_down*scaleFactor/scaleFactor_btagSF);
  
}


void FillGenInfoPlots(string dataSetName){
  if( foundSMb && foundW){
    histo1D[("Topmass_Wb_"+dataSetName).c_str()]->Fill((mcParticles[SMb_Indice] + mcParticles[SMW_Indice]).M());
    histo1D[("pt_Wb_"+dataSetName).c_str()]->Fill((mcParticles[SMb_Indice] + mcParticles[SMW_Indice]).Pt());
    histo1D[("eta_Wb_"+dataSetName).c_str()]->Fill((mcParticles[SMb_Indice] + mcParticles[SMW_Indice]).Eta());
    histo1D[("phi_Wb_"+dataSetName).c_str()]->Fill((mcParticles[SMb_Indice] + mcParticles[SMW_Indice]).Phi());
    histo1D[("dPhi_Wb_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticles[SMb_Indice],mcParticles[SMW_Indice]));
    histo1D[("dR_Wb_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticles[SMb_Indice],mcParticles[SMW_Indice]));
    
  }
  
  if( foundSMb && foundW && foundAntitopQ){
    //cout << "in tbar" << endl;
    histo1D[("Topmass_topquark_"+dataSetName).c_str()]->Fill(mcParticles[AntiTopQ_Indice].M());
    histo1D[("pt_topquark_"+dataSetName).c_str()]->Fill(mcParticles[AntiTopQ_Indice].Pt());
    histo1D[("eta_topquark_"+dataSetName).c_str()]->Fill(mcParticles[AntiTopQ_Indice].Eta());
    histo1D[("phi_topquark_"+dataSetName).c_str()]->Fill(mcParticles[AntiTopQ_Indice].Phi());
    
    histo2D[("Topmass_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[AntiTopQ_Indice].M() );
    histo2D[("pt_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[AntiTopQ_Indice].Pt() );
    histo2D[("phi_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[AntiTopQ_Indice].Phi() );
    histo2D[("eta_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[AntiTopQ_Indice].Eta() );
    
    histo1D[("dPhi_Wbvstopquark_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticles[SMb_Indice]+mcParticles[SMW_Indice],mcParticles[AntiTopQ_Indice]));
    histo1D[("dR_Wbvstopquark_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticles[SMb_Indice]+mcParticles[SMW_Indice],mcParticles[AntiTopQ_Indice]));
    
  }
  if( foundSMb && foundW && foundTopQ){
    //cout << "in top" << endl;
    histo1D[("Topmass_topquark_"+dataSetName).c_str()]->Fill(mcParticles[TopQ_Indice].M());
    histo1D[("pt_topquark_"+dataSetName).c_str()]->Fill(mcParticles[TopQ_Indice].Pt());
    histo1D[("eta_topquark_"+dataSetName).c_str()]->Fill(mcParticles[TopQ_Indice].Eta());
    histo1D[("phi_topquark_"+dataSetName).c_str()]->Fill(mcParticles[TopQ_Indice].Phi());
    
    histo1D[("dPhi_Wbvstopquark_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticles[SMb_Indice]+mcParticles[SMW_Indice],mcParticles[TopQ_Indice]));
    histo1D[("dR_Wbvstopquark_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticles[SMb_Indice]+mcParticles[SMW_Indice],mcParticles[TopQ_Indice]));
    
    
    histo2D[("Topmass_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[TopQ_Indice].M() );
    histo2D[("pt_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[TopQ_Indice].Pt() );
    histo2D[("phi_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[TopQ_Indice].Phi() );
    histo2D[("eta_topquarkvsWb_"+dataSetName).c_str()]->Fill((mcParticles[SMW_Indice] + mcParticles[SMb_Indice]).M(),mcParticles[TopQ_Indice].Eta() );
    
  }
  
  
  
  //FCNC TOP
  if( (foundZmumin && foundZmuplus)){
    histo1D[("Zmass_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).M());
    histo1D[("pt_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Pt());
    histo1D[("eta_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Eta());
    histo1D[("phi_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Phi());
    histo2D[("pt_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Pt(), mcParticles[Zmu_plus_Indice].Pt());
    histo2D[("phi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Phi(), mcParticles[Zmu_plus_Indice].Phi());
    histo2D[("eta_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Eta(), mcParticles[Zmu_plus_Indice].Eta());
    histo1D[("dPhi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticles[Zmu_min_Indice],mcParticles[Zmu_plus_Indice]));
    histo1D[("dR_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticles[Zmu_min_Indice],mcParticles[Zmu_plus_Indice]));
    histo1D[("mass_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].M());
    histo1D[("mass_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_plus_Indice].M());
    histo1D[("pt_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Pt());
    histo1D[("pt_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_plus_Indice].Pt());
    histo1D[("eta_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Eta());
    histo1D[("eta_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_plus_Indice].Eta());
    histo1D[("phi_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].Phi());
    histo1D[("phi_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_plus_Indice].Phi());
    histo2D[("mass_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zmu_min_Indice].M(), mcParticles[Zmu_plus_Indice].M());
    
  }
  if( (foundZelmin && foundZelplus)){
    histo1D[("Zmass_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).M());
    histo1D[("pt_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Pt());
    histo1D[("eta_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Eta());
    histo1D[("phi_Zleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Phi());
    histo2D[("pt_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Pt(), mcParticles[Zel_plus_Indice].Pt());
    histo2D[("phi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Phi(), mcParticles[Zel_plus_Indice].Phi());
    histo2D[("eta_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Eta(), mcParticles[Zel_plus_Indice].Eta());
    histo1D[("dPhi_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticles[Zel_min_Indice],mcParticles[Zel_plus_Indice]));
    histo1D[("dR_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticles[Zel_min_Indice],mcParticles[Zel_plus_Indice]));
    
    histo1D[("mass_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].M());
    histo1D[("mass_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_plus_Indice].M());
    histo1D[("pt_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Pt());
    histo1D[("pt_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_plus_Indice].Pt());
    histo1D[("eta_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Eta());
    histo1D[("eta_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_plus_Indice].Eta());
    histo1D[("phi_ZlepMin_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].Phi());
    histo1D[("phi_ZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_plus_Indice].Phi());
    histo2D[("mass_ZlepMinvsZlepPlus_"+dataSetName).c_str()]->Fill(mcParticles[Zel_min_Indice].M(), mcParticles[Zel_plus_Indice].M());
    
  }
  
  if( foundZ){
    histo1D[("Zmass_Zboson_"+dataSetName).c_str()]->Fill(mcParticles[Z_Indice].M());
    histo1D[("pt_Zboson_"+dataSetName).c_str()]->Fill(mcParticles[Z_Indice].Pt());
    histo1D[("phi_Zboson_"+dataSetName).c_str()]->Fill(mcParticles[Z_Indice].Phi());
    histo1D[("eta_Zboson_"+dataSetName).c_str()]->Fill(mcParticles[Z_Indice].Eta());
    
  }
  if( foundZ && (foundZmumin && foundZmuplus)){
    histo2D[("Zmass_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).M(),mcParticles[Z_Indice].M() );
    histo2D[("pt_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Pt(),mcParticles[Z_Indice].Pt() );
    histo2D[("phi_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Pt(),mcParticles[Z_Indice].Phi() );
    histo2D[("eta_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zmu_min_Indice] + mcParticles[Zmu_plus_Indice]).Pt(),mcParticles[Z_Indice].Eta() );
  }
  
  if( (foundZelmin && foundZelplus) && foundZ ){
    histo2D[("Zmass_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).M(),mcParticles[Z_Indice].M() );
    histo2D[("pt_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Pt(),mcParticles[Z_Indice].Pt() );
    histo2D[("phi_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Pt(),mcParticles[Z_Indice].Phi() );
    histo2D[("eta_ZbosonvsZleptons_"+dataSetName).c_str()]->Fill((mcParticles[Zel_min_Indice] + mcParticles[Zel_plus_Indice]).Pt(),mcParticles[Z_Indice].Eta() );
    
  }
  
  Double_t genmWT;
  Double_t genphis ;
  Double_t gencosphis;
  Double_t genmWT2 ;
  Double_t genCosTheta ;
  Double_t genCosTheta_alt ;
  
  
  if(SMmu_Indice != -999 && SMnumu_Indice != -999) {
    genmWT= TMath::Sqrt((mcParticles[SMmu_Indice].Pt() + mcParticles[SMnumu_Indice].Pt())*(mcParticles[SMmu_Indice].Pt() +mcParticles[SMnumu_Indice].Pt())-(mcParticles[SMmu_Indice].Px() + mcParticles[SMnumu_Indice].Px())*(mcParticles[SMmu_Indice].Px() + mcParticles[SMnumu_Indice].Px()) - (mcParticles[SMmu_Indice].Py() + mcParticles[SMnumu_Indice].Py())* (mcParticles[SMmu_Indice].Py() +mcParticles[SMnumu_Indice].Py()));
    genphis = mcParticles[SMmu_Indice].Phi() - met_Phi;
    gencosphis = TMath::Cos(genphis);
    genmWT2 = TMath::Sqrt(2*mcParticles[SMmu_Indice].Pt() * mcParticles[SMnumu_Indice].Pt()*(1-gencosphis));
    
    if(SMb_Indice!=-999) genCosTheta =  CosThetaCalculation(mcParticles[SMmu_Indice], mcParticles[SMnumu_Indice], mcParticles[SMb_Indice], false).first;
    if(SMb_Indice!=-999) genCosTheta_alt =  CosThetaCalculation(mcParticles[SMmu_Indice], mcParticles[SMnumu_Indice], mcParticles[SMb_Indice], false).second;
    
  }
  if(SMel_Indice != -999 && SMnuel_Indice != -999) {
    genmWT= TMath::Sqrt((mcParticles[SMel_Indice].Pt() + mcParticles[SMnuel_Indice].Pt())*(mcParticles[SMel_Indice].Pt() +mcParticles[SMnuel_Indice].Pt())-(mcParticles[SMel_Indice].Px() + mcParticles[SMnuel_Indice].Px())*(mcParticles[SMel_Indice].Px() + mcParticles[SMnuel_Indice].Px()) - (mcParticles[SMel_Indice].Py() + mcParticles[SMnuel_Indice].Py())* (mcParticles[SMel_Indice].Py() +mcParticles[SMnuel_Indice].Py()));
    genphis = mcParticles[SMel_Indice].Phi() - met_Phi;
    gencosphis = TMath::Cos(genphis);
    genmWT2 = TMath::Sqrt(2*mcParticles[SMel_Indice].Pt() * mcParticles[SMnuel_Indice].Pt()*(1-gencosphis));
    
    if(SMb_Indice!=-999) genCosTheta =  CosThetaCalculation(mcParticles[SMel_Indice], mcParticles[SMnuel_Indice], mcParticles[SMb_Indice], false).first;
    if(SMb_Indice!=-999) genCosTheta_alt =  CosThetaCalculation(mcParticles[SMel_Indice], mcParticles[SMnuel_Indice], mcParticles[SMb_Indice], false).second;
  }
  
  
  histo1D[("GenInfo_mWT_"+dataSetName).c_str()]->Fill(genmWT);
  histo1D[("GenInfo_mWT2_"+dataSetName).c_str()] ->Fill(genmWT2);
  histo2D[("GenInfo_mWTmWT2_"+dataSetName).c_str()]->Fill(genmWT,genmWT2);
  
  histo1D[("GenInfo_CosTheta_"+dataSetName).c_str()] ->Fill(genCosTheta);
  histo1D[("GenInfo_CosTheta_alt_"+dataSetName).c_str()] ->Fill(genCosTheta_alt);
  histo2D[("GenInfo_CosTheta2_"+dataSetName).c_str()] ->Fill(genCosTheta, genCosTheta_alt);
  
  
}
void FillRecovsGenInfoPlots(string dataSetName, vector<TLorentzVector> selectedElectrons, vector <TLorentzVector> selectedMuons , vector <TLorentzVector> selectedJets){
  
  //cout << ("Matched_mlb_"+dataSetName).c_str() << endl;
  cout << "WmuIndiceM "<< WmuIndiceM <<" BjetIndiceM "<< BjetIndiceM <<" selectedJets "<< selectedJets.size() <<" selectedMuons "<< selectedMuons.size() << endl;
  if(WmuIndiceM != -999 && BjetIndiceM != -999)  histo1D[("Matched_mlb_"+dataSetName).c_str()]->Fill((selectedJets[BjetIndiceM] + selectedMuons[WmuIndiceM]).M());
  if(WelecIndiceM != -999 && BjetIndiceM != -999)  histo1D[("Matched_mlb_"+dataSetName).c_str()]->Fill((selectedJets[BjetIndiceM] + selectedElectrons[WelecIndiceM]).M());
  
  if(ZmuIndiceM_1 != -999 && ZmuIndiceM_0 != -999) histo1D[("Matched_Zmass_"+dataSetName).c_str()]->Fill((selectedMuons[ZmuIndiceM_0] + selectedMuons[ZmuIndiceF_1]).M());
  if(ZelecIndiceM_1 != -999 && ZelecIndiceM_0 != -999) histo1D[("Matched_Zmass_"+dataSetName).c_str()]->Fill((selectedElectrons[ZelecIndiceM_0] + selectedElectrons[ZelecIndiceF_1]).M());
  
  
}
void FillMVAPlots(int d, string dataSetName, int Region, string prefix, vector<int> decayChannels ){
  clock_t start_sub = clock();
  //cout << "in MVA plots" << endl;
  string sregion;
  sregion = prefix;
  string decaystring = "";
  
  //if(dataSetName.find("fake")!=std::string::npos) scaleFactor *= 0.0001;
  
  //cout <<  "region " << sregion << endl;
  for(int iChan = 0;iChan < decayChannels.size() ;iChan++){
    // cout << "chan " << decayChannels[iChan] << " chan in evt " << MVA_channel << endl;
    if(decayChannels[iChan]!= -9 && decayChannels[iChan] != MVA_channel){
      
      // cout << "continuing " <<decayChannels[iChan] << " is not " << MVA_channel << endl;
      continue;
    }
    
    decaystring= "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    else if(decayChannels[iChan] == 1) decaystring = "uue";
    else if(decayChannels[iChan] == 2) decaystring = "eeu";
    else if(decayChannels[iChan] == 3) decaystring = "eee";
    else  if(decayChannels[iChan] == 4 || decayChannels[iChan] == 5 ) continue;
    
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << decaystring << endl;
    MSPlot[(sregion+"_MVA_Zboson_eta_"+decaystring).c_str()]->Fill(MVA_Zboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_channel_"+decaystring).c_str()]->Fill(MVA_channel, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_weight_"+decaystring).c_str()]->Fill(MVA_weight, datasets[d], true, 1);
    MSPlot[(sregion+"_MVA_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    //  cout << "filling " << (sregion+"_MVA_lepton0_pt_"+decaystring).c_str() << endl;
    MSPlot[(sregion+"_MVA_lepton0_pt_"+decaystring).c_str()]->Fill(MVA_lepton0_pt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton1_pt_"+decaystring).c_str()]->Fill(MVA_lepton1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_SMtop_eta_"+decaystring).c_str()]->Fill(MVA_SMtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton2_pt_"+decaystring).c_str()]->Fill(MVA_lepton2_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton0_eta_"+decaystring).c_str()]->Fill(MVA_lepton0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton1_eta_"+decaystring).c_str()]->Fill(MVA_lepton1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton2_eta_"+decaystring).c_str()]->Fill(MVA_lepton2_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton0_phi_"+decaystring).c_str()]->Fill(MVA_lepton0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton1_phi_"+decaystring).c_str()]->Fill(MVA_lepton1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_lepton2_phi_"+decaystring).c_str()]->Fill(MVA_lepton2_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_jet0_pt_"+decaystring).c_str()]->Fill(MVA_jet0_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_jet0_eta_"+decaystring).c_str()]->Fill(MVA_jet0_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_jet0_phi_"+decaystring).c_str()]->Fill(MVA_jet0_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJets.size()>1) MSPlot[(sregion+"_MVA_jet1_pt_"+decaystring).c_str()]->Fill(MVA_jet1_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJets.size()>1) MSPlot[(sregion+"_MVA_jet1_eta_"+decaystring).c_str()]->Fill(MVA_jet1_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJets.size()>1) MSPlot[(sregion+"_MVA_jet1_phi_"+decaystring).c_str()]->Fill(MVA_jet1_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    // SM side
    MSPlot[(sregion+"_MVA_Wlep_pt_"+decaystring).c_str()]->Fill(MVA_Wlep_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wlep_eta_"+decaystring).c_str()]->Fill(MVA_Wlep_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wlep_phi_"+decaystring).c_str()]->Fill(MVA_Wlep_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_SMbjet_pt_"+decaystring).c_str()]->Fill(MVA_SMbjet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_SMbjet_eta_"+decaystring).c_str()]->Fill(MVA_SMbjet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_SMbjet_phi_"+decaystring).c_str()]->Fill(MVA_SMbjet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wboson_pt_"+decaystring).c_str()]->Fill(MVA_Wboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wboson_eta_"+decaystring).c_str()]->Fill(MVA_Wboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wboson_phi_"+decaystring).c_str()]->Fill(MVA_Wboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_met_"+decaystring).c_str()]->Fill(MVA_met, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_SMtop_pt_"+decaystring).c_str()]->Fill(MVA_SMtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_SMtop_phi_"+decaystring).c_str()]->Fill(MVA_SMtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    // FCNC side
    // MSPlot[(sregion+"_MVA_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Zboson_eta_"+decaystring).c_str()]->Fill(MVA_Zboson_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Zboson_phi_"+decaystring).c_str()]->Fill(MVA_Zboson_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    // nbrs
    MSPlot[(sregion+"_MVA_nMuons_"+decaystring).c_str()]->Fill(MVA_nMuons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_NJets_CSVv2T_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2T, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_NJets_CSVv2M_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_NJets_CSVv2L_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2L, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_nJets_"+decaystring).c_str()]->Fill(MVA_nJets, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_nElectrons_"+decaystring).c_str()]->Fill(MVA_nElectrons, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    
    //SM kinematics
    MSPlot[(sregion+"_MVA_SMtop_M_"+decaystring).c_str()]->Fill(MVA_SMtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_mlb_"+decaystring).c_str()]->Fill(MVA_mlb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Wboson_M_"+decaystring).c_str()]->Fill(MVA_Wboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_dRWlepb_"+decaystring).c_str()]->Fill(MVA_dRWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_dPhiWlepb_"+decaystring).c_str()]->Fill(MVA_dPhiWlepb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_Wlep_Charge_"+decaystring).c_str()]->Fill(MVA_Wlep_Charge , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalPt_"+decaystring).c_str()]->Fill(MVA_TotalPt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalHt_"+decaystring).c_str()]->Fill(MVA_TotalHt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalInvMass_"+decaystring).c_str()]->Fill( MVA_TotalInvMass, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_TotalPt_jet_"+decaystring).c_str()]->Fill(MVA_TotalPt_jet , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalHt_jet_"+decaystring).c_str()]->Fill(MVA_TotalHt_jet , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalInvMass_jet_"+decaystring).c_str()]->Fill( MVA_TotalInvMass_jet, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_TotalPt_lep_"+decaystring).c_str()]->Fill(MVA_TotalPt_lep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalHt_lep_"+decaystring).c_str()]->Fill(MVA_TotalHt_lep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_TotalInvMass_lep_"+decaystring).c_str()]->Fill( MVA_TotalInvMass_lep, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
   if(selectedJetsID.size()>0) MSPlot[(sregion+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
   if(selectedJetsID.size()>1)  MSPlot[(sregion+"_MVA_bdiscCSVv2_jet_1_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    // MSPlot[(sregion+"_MVA_CosTheta_"+decaystring).c_str()]->Fill(MVA_CosTheta , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    //  MSPlot[(sregion+"_MVA_CosTheta_alt_"+decaystring).c_str()]->Fill(MVA_CosTheta_alt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
   if(selectedJetsID.size()>1) MSPlot[(sregion+"_MVA_cdiscCvsB_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJetsID.size()>1) MSPlot[(sregion+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_1 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJetsID.size()>0) MSPlot[(sregion+"_MVA_cdiscCvsB_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_0 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    if(selectedJetsID.size()>0) MSPlot[(sregion+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_0 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    
    
    // FCNC kinematics
    
    MSPlot[(sregion+"_MVA_dRZb_"+decaystring).c_str()]->Fill(MVA_dRZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_dRZWlep_"+decaystring).c_str()]->Fill(MVA_dRZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_dRZSMtop_"+decaystring).c_str()]->Fill(MVA_dRZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    MSPlot[(sregion+"_MVA_dPhiZb_"+decaystring).c_str()]->Fill(MVA_dPhiZb , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_dPhiZWlep_"+decaystring).c_str()]->Fill(MVA_dPhiZWlep , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_dPhiZMET_"+decaystring).c_str()]->Fill(MVA_dPhiZMET , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_dPhiZSMtop_"+decaystring).c_str()]->Fill(MVA_dPhiZSMtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_m3l_"+decaystring).c_str()]->Fill(MVA_m3l , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    MSPlot[(sregion+"_MVA_mWt_"+decaystring).c_str()]->Fill(MVA_mWt , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_mWt2_"+decaystring).c_str()]->Fill(MVA_mWt2 , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    MSPlot[(sregion+"_MVA_nJets_CharmL_"+decaystring).c_str()]->Fill(MVA_nJets_CharmL, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_nJets_CharmM_"+decaystring).c_str()]->Fill(MVA_nJets_CharmM, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_nJets_CharmT_"+decaystring).c_str()]->Fill(MVA_nJets_CharmT, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    MSPlot[(sregion+"_MVA_Zboson_M_"+decaystring).c_str()]->Fill(MVA_Zboson_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
    
    
    if(Region == 1){
      MSPlot[(sregion+"_MVA_LightJet_pt_"+decaystring).c_str()]->Fill(MVA_LightJet_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_LightJet_eta_"+decaystring).c_str()]->Fill(MVA_LightJet_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_LightJet_phi_"+decaystring).c_str()]->Fill(MVA_LightJet_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"_MVA_FCNCtop_pt_"+decaystring).c_str()]->Fill(MVA_FCNCtop_pt, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_FCNCtop_eta_"+decaystring).c_str()]->Fill(MVA_FCNCtop_eta, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_FCNCtop_phi_"+decaystring).c_str()]->Fill(MVA_FCNCtop_phi, datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      
      MSPlot[(sregion+"_MVA_FCNCtop_M_"+decaystring).c_str()]->Fill(MVA_FCNCtop_M , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
      MSPlot[(sregion+"_MVA_dRZc_"+decaystring).c_str()]->Fill(MVA_dRZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_dPhiZc_"+decaystring).c_str()]->Fill(MVA_dPhiZc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
     
      // interplay
      MSPlot[(sregion+"_MVA_dRSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dRSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_dRWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_dPhiSMFCNCtop_"+decaystring).c_str()]->Fill(MVA_dPhiSMFCNCtop , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      MSPlot[(sregion+"_MVA_dPhiWlepc_"+decaystring).c_str()]->Fill(MVA_dRWlepc , datasets[d], true, Luminosity*scaleFactor/EquilumiSF);
      
     
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

///////////////////////////////////// MATCHING /////////////////////////////////////////
bool MatchingFunction(string dataSetName, vector <TLorentzVector> Leptons, vector <TLorentzVector> selectedMuons, vector<TLorentzVector> selectedElectrons , vector <TLorentzVector> selectedJets, bool makePlots, bool debug){
  ClearMatchingVars(); // to do with each new dataset
  ClearMatchingVarsTLV();
  int MuonSize = selectedMuons.size();
  for (int iMC = 0; iMC < nMCParticles; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC],mc_phi[iMC],mc_E[iMC]);
    mcParticles.push_back(mcpart);
  }
  if(mcParticles.size() != nMCParticles){cout << "ERROR mcP not filled correctly" << endl;  }
  
  if(debug) cout << "looking at " << mcParticles.size() << " mcParticles" << endl;
  
  EventSearcher(mcParticles, dataSetName, debug);
  
  if(debug) cout << "found decay " << foundDecay << endl;
  
  if((foundSMmu || foundSMel) &&
     ((foundZelmin && foundZelplus) || (foundZmuplus && foundZmumin)) &&
     foundSMb ){
    foundDecay = true;
  }
  if(debug) cout << "found decay " << foundDecay << endl;
  
  if(!foundDecay){ return false;  }
  
  if(makePlots){
    // InitGenInfoPlots(dataSetName);
    FillGenInfoPlots(dataSetName);
  }
  
  bool matchedEvent = true;
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> > OutputLeptonMatching;
  OutputLeptonMatching =  LeptonMatching(Leptons, mcParticles, dataSetName, debug);
  //cout << "WmuIndiceM "<< WmuIndiceM << endl;
  for(unsigned int iPart = 0 ; iPart < (OutputLeptonMatching.second).size(); iPart++){
    //cout << "(OutputLeptonMatching.second)[iPart] " << (OutputLeptonMatching.second)[iPart] << endl;
    if((OutputLeptonMatching.second)[iPart].find("SMmu")!=string::npos){ WmuIndiceM = (OutputLeptonMatching.first)[iPart].first ;
      //cout << "OutputLeptonMatching.first)[iPart].first " << (OutputLeptonMatching.first)[iPart].first << endl;
    }
    if((OutputLeptonMatching.second)[iPart].find("SMel")!=string::npos){ WelecIndiceM = (OutputLeptonMatching.first)[iPart].first - MuonSize;  if(debug){ cout << "found SM el" << endl; }}
    if((OutputLeptonMatching.second)[iPart].find("Zmumin")!=string::npos){ ZmuIndiceM_0 = ((OutputLeptonMatching.first)[iPart].first);  if(debug){ cout << "found FCNC mu" << endl; } }
    if((OutputLeptonMatching.second)[iPart].find("Zelmin")!=string::npos){ ZelecIndiceM_0 = (OutputLeptonMatching.first)[iPart].first - MuonSize; if(debug){ cout << "found FCNC el" << endl; }}
    if((OutputLeptonMatching.second)[iPart].find("Zmuplus")!=string::npos){ ZmuIndiceM_1 = ((OutputLeptonMatching.first)[iPart].first  ); if(debug){ cout << "found FCNC mu" << endl; }}
    if((OutputLeptonMatching.second)[iPart].find("Zelplus")!=string::npos){ZelecIndiceM_1 = (OutputLeptonMatching.first)[iPart].first - MuonSize; if(debug){ cout << "found FCNC el" << endl; }}
  }
  
  
  if(WmuIndiceF == WmuIndiceM && WmuIndiceM != -999 ) WlepMatched++;
  else if(WelecIndiceF == WelecIndiceM && WelecIndiceM != -999 ) WlepMatched++;
  else matchedEvent = false;
  if(ZmuIndiceF_0 == ZmuIndiceM_0 && ZmuIndiceF_1 == ZmuIndiceM_1 && ZmuIndiceM_1!= -999 && ZmuIndiceM_0 != -999  ) ZlepMatched++;
  else if(ZmuIndiceF_1== ZmuIndiceM_0 && ZmuIndiceF_0 == ZmuIndiceM_1 && ZmuIndiceM_1!= -999 && ZmuIndiceM_0 != -999 ) ZlepMatched++;
  else matchedEvent = false;
  if(ZelecIndiceM_0== ZelecIndiceF_0 && ZelecIndiceF_1 == ZelecIndiceM_1 && ZelecIndiceM_1!= -999 && ZelecIndiceM_0 != -999) ZlepMatched++;
  else if(ZelecIndiceM_1== ZelecIndiceF_0 && ZelecIndiceF_0 == ZelecIndiceM_1 && ZelecIndiceM_1!= -999 && ZelecIndiceM_0 != -999 ) ZlepMatched++;
  else matchedEvent = false;
  
  if(WmuIndiceM != -999 || WelecIndiceM != -999 ) WlepMatchedevent++;;
  if( (ZmuIndiceM_0 != -999 && ZmuIndiceM_1 != -999) || (ZelecIndiceM_0 != -999 && ZelecIndiceM_1 != -999)) ZlepMatchedevent++;
  
  if(debug){
    cout << "******* numbers for matching *********" << endl;
    cout << "WmuIndiceM "<< WmuIndiceM <<" WmuIndiceF "<< WmuIndiceF << endl;
    cout << "WelecIndiceM "<< WelecIndiceM <<" WelecIndiceF "<< WelecIndiceF << endl;
    cout << "ZelecIndiceM_1 "<< ZelecIndiceM_1 <<" ZelecIndiceM_0 "<< ZelecIndiceM_0 << endl;
    cout << "ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
    cout << "ZmuIndiceM_0 "<< ZmuIndiceM_0 <<" ZmuIndiceM_1 "<< ZmuIndiceM_1 << endl;
    cout << "ZmuIndiceF_0 "<< ZmuIndiceF_0 <<" ZmuIndiceF_1 "<< ZmuIndiceF_1 << endl;
    cout << "WlepMatched "<< WlepMatched <<" WlepMatchedevent "<< WlepMatchedevent << endl;
    cout << "ZlepMatched "<< ZlepMatched <<" ZlepMatchedevent "<< ZlepMatchedevent << endl;
    
  }
  
  if(debug)cout << "**** OutputJetMatching begin ****" << endl;
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> > OutputJetMatching;
  OutputJetMatching =  JetMatching(selectedJets, mcParticles, dataSetName, true);
  if(debug) cout << "**** OutputJetMatching end ****" << endl;
  
  for(unsigned int iPart = 0 ; iPart < (OutputJetMatching.second).size(); iPart++){
    if((OutputJetMatching.second)[iPart].find("SMb")!=string::npos){ BjetIndiceM = (OutputJetMatching.first)[iPart].first  ; }
    if((OutputJetMatching.second)[iPart].find("FCNCc")!=string::npos){ CjetIndiceM = (OutputJetMatching.first)[iPart].first  ;}
    if((OutputJetMatching.second)[iPart].find("FCNCu")!=string::npos){ UjetIndiceM = (OutputJetMatching.first)[iPart].first  ;}
  }
  
  if(cjetindex == CjetIndiceM && CjetIndiceM != -999) CjetMatched++;
  if(cjetindex_Cloose == CjetIndiceM && CjetIndiceM != -999) CjetMatchedL++;
  if(cjetindex_Cmedium == CjetIndiceM && CjetIndiceM != -999) CjetMatchedM++;
  if(cjetindex_Ctight == CjetIndiceM && CjetIndiceM != -999) CjetMatchedT++;
  if(cjetindex_CvsBtagger == CjetIndiceM && CjetIndiceM != -999) CjetMatchedCvsBL++;
  if(cjetindex_CvsLtagger == CjetIndiceM && CjetIndiceM != -999) CjetMatchedCvsLL++;
  
  if(SMjetIndex == BjetIndiceM && BjetIndiceM != -999) BjetMatched++;
  else matchedEvent = false;
  if(cjetindex == UjetIndiceM && UjetIndiceM != -999) UjetMatched++;
  
  if(CjetIndiceM != -999 && cjetindex != -999) CjetMatchedevent++;
  if(CjetIndiceM != -999 && cjetindex_Cloose != -999) CjetMatchedeventL++;
  if(CjetIndiceM != -999 && cjetindex_Cmedium != -999) CjetMatchedeventM++;
  if(CjetIndiceM != -999 && cjetindex_Ctight != -999) CjetMatchedeventT++;
  if(CjetIndiceM != -999 && cjetindex_CvsBtagger != -999) CjetMatchedeventCvsBL++;
  if(CjetIndiceM != -999 && cjetindex_CvsLtagger != -999) CjetMatchedeventCvsLL++;
  
  
  if(BjetIndiceM != -999 && SMjetIndex != -999) BjetMatchedevent++;
  if(UjetIndiceM != -999 && cjetindex != -999) UjetMatchedevent++;
  
  //if(makePlots) FillRecovsGenInfoPlots(dataSetName, selectedElectrons, selectedMuons , selectedJets); // TO FIX
  
  return true;
}
void MatchingEfficiency(){
  
  cout << "                LEPTON MATCHING" << endl;
  cout << "                  - W lepton matching: " << WlepMatched << " out of " << WlepMatchedevent << " or " << (WlepMatched/WlepMatchedevent)*100. << " % " << endl;
  cout << "                  - Z lepton matching: " << ZlepMatched << " out of " << ZlepMatchedevent << " or " << (ZlepMatched/ZlepMatchedevent)*100. << " % " << endl;
  cout << "                  - Lepton matched : " << ZlepMatched + WlepMatched << " out of " << WlepMatchedevent + ZlepMatchedevent << " or " << ((WlepMatched+ZlepMatched)/(WlepMatchedevent+ZlepMatchedevent))*100. << " % " << endl;
  
  cout << "                JET MATCHTING" << endl;
  cout << "                  - B jet matching: " << BjetMatched << " out of " << BjetMatchedevent << " or " << (BjetMatched/BjetMatchedevent)*100. << " % " << endl;
  cout << "                  - C jet matching based on mass req: " << CjetMatched << " out of " << CjetMatchedevent << " or " << (CjetMatched/CjetMatchedevent)*100. << " % " << endl;
  cout << "                  - C jet matching based on CvsL tagger req: " << CjetMatchedCvsLL << " out of " << CjetMatchedeventCvsLL << " or " << (CjetMatchedCvsLL/CjetMatchedeventCvsLL)*100. << " % " << endl;
  cout << "                  - C jet matching based on CvsB tagger req: " << CjetMatchedCvsBL << " out of " << CjetMatchedeventCvsBL << " or " << (CjetMatchedCvsBL/CjetMatchedeventCvsBL)*100. << " % " << endl;
  cout << "                  - C jet matching based on WP L, Pt req: " << CjetMatchedL << " out of " << CjetMatchedeventL << " or " << (CjetMatchedL/CjetMatchedeventL)*100. << " % " << endl;
  cout << "                  - C jet matching based on WP M, Pt req: " << CjetMatchedM << " out of " << CjetMatchedeventM << " or " << (CjetMatchedM/CjetMatchedeventM)*100. << " % " << endl;
  cout << "                  - C jet matching based on WP T, Pt req: " << CjetMatchedT << " out of " << CjetMatchedeventT << " or " << (CjetMatchedT/CjetMatchedeventT)*100. << " % " << endl;
  cout << "                  - U jet matching based in mass req: " << UjetMatched << " out of " << UjetMatchedevent << " or " << (UjetMatched/UjetMatchedevent)*100. << " % " << endl;
};





pair< vector< pair<unsigned int, unsigned int>>, vector <string> > LeptonMatching(vector < TLorentzVector> selectedleptons, vector <TLorentzVector> mcParticles, string dataSetName, bool debug){
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  for (unsigned int iMC = 0; iMC < mcParticles.size(); iMC++)
  {
    if ( (mc_status[iMC] > 1 && mc_status[iMC] <= 20) || mc_status[iMC] >= 30 ) continue;  /// Final state particle or particle from hardest process
    if( abs(mc_pdgId[iMC]) ==  13 || abs(mc_pdgId[iMC]) ==  11 ){
      partons.push_back(mcParticles[iMC]);
      partonID.push_back(iMC);
    } // leptons
  }
  
  if(debug) cout << "nb of partons to be matched " << partons.size() << endl;
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; //
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedleptons,2,true,true,0.1 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedleptons collection
    // i is nb in partons
    // partonID contains place in mcParticles vector
  }
  
  
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    // cout << "mc_pdgId[partonID[partonIDnb]] " << mc_pdgId[partonID[partonIDnb]] << " mc_mother[partonID[partonIDnb]] " << mc_mother[partonID[partonIDnb]] << " mc_granny[partonID[partonIDnb]] " << endl;
    
    //SM
    if( abs(mc_pdgId[partonID[partonIDnb]]) ==  13 && abs(mc_mother[partonID[partonIDnb]])  == 24 && abs(mc_granny[partonID[partonIDnb]])  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
      NPair.push_back("SMmu");
      //cout << "i " << i << " SM mu = JetPartonPair[i].first " << JetPartonPair[i].first << endl;
      
    } // mu from W from t
    if( abs(mc_pdgId[partonID[partonIDnb]]) ==  11 && abs(mc_mother[partonID[partonIDnb]])  == 24 && abs(mc_granny[partonID[partonIDnb]])  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
      NPair.push_back("SMel");
    } // el from W from t
    
    if(dataSetName.find("TT_FCNC") ==std::string::npos ){
      //FCNC
      if( mc_pdgId[partonID[partonIDnb]] ==  13 && abs(mc_mother[partonID[partonIDnb]])  == 23 ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zmumin");
      } // mu- from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  11 && abs(mc_mother[partonID[partonIDnb]])  == 23  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zelmin");
      } // el- from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  -13 && abs(mc_mother[partonID[partonIDnb]])  == 23  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zmuplus");
      } // mu+ from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  -11 && abs(mc_mother[partonID[partonIDnb]])  == 23  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zelplus");
      } // el+ from Z from t
    }
    if(dataSetName.find("TT_FCNC")!=std::string::npos ){
      if( mc_pdgId[partonID[partonIDnb]] ==  13 && abs(mc_mother[partonID[partonIDnb]])  == 23&& abs(mc_granny[partonID[partonIDnb]])  == 6  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zmumin");
      } // mu- from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  11 && abs(mc_mother[partonID[partonIDnb]])  == 23 && abs(mc_granny[partonID[partonIDnb]])  == 6  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zelmin");
      } // el- from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  -13 && abs(mc_mother[partonID[partonIDnb]])  == 23 && abs(mc_granny[partonID[partonIDnb]])  == 6  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zmuplus");
      } // mu+ from Z from t
      if( mc_pdgId[partonID[partonIDnb]] ==  -11 && abs(mc_mother[partonID[partonIDnb]])  == 23  && abs(mc_granny[partonID[partonIDnb]])  == 6 ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("Zelplus");
      } // el+ from Z from t
      
    }
    
  }
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
  
}
pair< vector< pair<unsigned int, unsigned int>>, vector <string> > JetMatching(vector < TLorentzVector> selectedJets, vector <TLorentzVector> mcParticles, string dataSetName, bool debug){
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  debug = false;
  if(debug) cout << "begin jet matching" << endl;
  for (unsigned int iMC = 0; iMC < mcParticles.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mc_status[iMC] << "  pdgId: " << setw(3) << mc_pdgId[iMC] << "  Mother: " << setw(4) << mc_mother[iMC] << "  Granny: " << setw(4) << mc_granny[iMC] << endl;
    
    //cout << "iMC " << iMC << endl;
    if ( (mc_status[iMC] > 1 && mc_status[iMC] <= 20) || mc_status[iMC] >= 30 ) continue;  /// Final state particle or particle from hardest process
    //cout << "iMC through " << iMC << endl;
    
    if( abs(mc_pdgId[iMC]) <=  5 || abs(mc_pdgId[iMC]) ==  21 ){
      //cout << "pushing back" << endl;
      partons.push_back(mcParticles[iMC]);
      partonID.push_back(iMC);
    } //jets and gluons
    
    // cout << "iMC done " << iMC << endl;
    
  }
  
  if(debug) cout << "partons to be matched: " << partons.size() <<endl;
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedJets,2,true,true,0.3 );
  
  if(debug)  cout << "jet parton matching done " << endl;
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedJets collection
    // i is nb in partons
    // partonID contains place in mcParticles vector
  }
  if(debug) cout << "pushing back" << endl;
  
  /*
   if(JetPartonPair.size() < 1 && dataSetName.find("TT")==std::string::npos){
   pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
   return returnVectorEr;
   }
   else if(JetPartonPair.size() < 2 && dataSetName.find("TT")!=std::string::npos){
   pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
   return returnVectorEr;
   }
   */
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; //
  NPair.clear();
  
  if(debug) cout << "get pairs" << endl;
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    if( abs(mc_pdgId[partonID[partonIDnb]]) ==  5 && abs(mc_mother[partonID[partonIDnb]])  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
      NPair.push_back("SMb");
    } // b from t
    
    //FCNC
    if(dataSetName.find("TT")!=std::string::npos){
      if( abs(mc_pdgId[partonID[partonIDnb]]) ==  2  && abs(mc_mother[partonID[partonIDnb]])  == 6 ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("FCNCu");
      } // q from t
      if( abs(mc_pdgId[partonID[partonIDnb]]) ==  4 && abs(mc_mother[partonID[partonIDnb]])  == 6 ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("FCNCc");
      } // q from t
    }
    else {
      
      if( abs(mc_pdgId[partonID[partonIDnb]]) ==  2  ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("FCNCu");
      } // q from t
      if( abs(mc_pdgId[partonID[partonIDnb]]) ==  4 ){
        PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,partonID[partonIDnb]));
        NPair.push_back("FCNCc");
      } // q from t
    }
  }
  
  if(debug) cout << "return pairs" << endl;
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
  
  
}
void EventSearcher(vector < TLorentzVector> mcParticles, string dataSetName, bool debug){
  
  for (unsigned int iMC = 0; iMC < mcParticles.size(); iMC++)
  {
    if (debug)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mc_status[iMC] << "  pdgId: " << setw(3) << mc_pdgId[iMC]<< "  Mother: " << setw(4) << mc_mother[iMC] << "  Granny: " << setw(4) << mc_granny[iMC] << "  Pt: " << setw(7) << left << mc_pt[iMC] << "  Eta: " << mc_eta[iMC] << endl;
    
    
    if ( (mc_status[iMC] > 1 && mc_status[iMC] <= 20) || mc_status[iMC] >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mc_pdgId[iMC]== 6 ){
      TopQ_Indice = iMC; foundTopQ = true;
    }
    else if( mc_pdgId[iMC]== -6 ){ AntiTopQ_Indice = iMC; foundAntitopQ = true; }
    
    //SM
    else if( abs( mc_pdgId[iMC]) ==  5 && abs(mc_mother[iMC])  == 6 ){
      if(mc_status[iMC] == 23) {SMb_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMb_Indice == -999) {SMb_Indice = iMC;}
      foundSMb = true;
      
    } // b from t
    else if( mc_pdgId[iMC]==  13 && mc_mother[iMC]  == -24 && mc_granny[iMC]  == -6 ){
      if(mc_status[iMC] == 23) {SMmu_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMmu_Indice == -999) {SMmu_Indice = iMC;}
      foundSMmu = true;
      //cout << "found mu from tbar" << endl;
      
    }
    else if( abs(mc_pdgId[iMC])==  14 && mc_mother[iMC]  == -24 && mc_granny[iMC]  == -6 ){
      if(mc_status[iMC] == 23) {SMnumu_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMnumu_Indice == -999) {SMnumu_Indice = iMC;}
      foundSMnumu = true;
      
    }
    else if( abs(mc_pdgId[iMC])==  12 && mc_mother[iMC]  == -24 && mc_granny[iMC]  == -6 ){
      if(mc_status[iMC] == 23) {SMnuel_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMnuel_Indice== -999) {SMnuel_Indice = iMC;}
      foundSMnuel = true;
    }
    else if( mc_pdgId[iMC]==  -13 && mc_mother[iMC]  == 24 && mc_granny[iMC]  == 6 ){
      if(mc_status[iMC] == 23) {SMmu_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMmu_Indice == -999) {SMmu_Indice = iMC;}
      foundSMmu = true;
      
    } // mu+ from W+ from top
    else if( mc_pdgId[iMC]==  11 && mc_mother[iMC]  == -24 && mc_granny[iMC]  == -6 ){
      if(mc_status[iMC] == 23) {SMel_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMel_Indice== -999) {SMel_Indice = iMC;}
      foundSMel = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mc_pdgId[iMC]==  -11 && mc_mother[iMC]  == 24 && mc_granny[iMC]  == 6 ){
      if(mc_status[iMC] == 23) {SMel_Indice = iMC;  }
      else if( mc_status[iMC] != 23 && SMel_Indice== -999) {SMel_Indice = iMC;}
      foundSMel = true;
      
    } // el+ from W+ from top
    else if( mc_pdgId[iMC]==  24 && mc_mother[iMC]  == 6   ){
      if(mc_status[iMC] == 23) {SMW_Indice= iMC; }
      else if( mc_status[iMC] != 23 && SMW_Indice== -999) {SMW_Indice = iMC; }
      foundW = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mc_pdgId[iMC]==  -24 && mc_mother[iMC]  == -6  ){
      if(mc_status[iMC] == 23) {SMW_Indice= iMC; }
      else if( mc_status[iMC] != 23 && SMW_Indice== -999) {SMW_Indice = iMC;}
      foundW = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    // u jet
    if(dataSetName.find("TT")==std::string::npos){
      //FCNC
      if( mc_pdgId[iMC]==  13 && mc_mother[iMC]  == 23  ){
        if(mc_status[iMC] == 23) {Zmu_min_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zmu_min_Indice== -999) {Zmu_min_Indice = iMC;}
        foundZmumin = true;
        
        // cout << "mu from Z from tbar found  " << endl;
      } // mu - from Z  from tbar
      else if( mc_pdgId[iMC]==  -13 && mc_mother[iMC]  == 23  ){
        if(mc_status[iMC] == 23) {Zmu_plus_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zmu_plus_Indice== -999) {Zmu_plus_Indice= iMC;}
        foundZmuplus = true;
        
        //cout << "mu from Z from tbar found  " << endl;
      } // mu + from Z  from tbar
      
      else if( abs( mc_pdgId[iMC]) ==  4  ){
        if(mc_status[iMC] == 23) {Cjet_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Cjet_Indice == -999) {Cjet_Indice = iMC;}
        foundCjet = true;
        
      } // c jet
      else if( abs( mc_pdgId[iMC]) ==  2  ){
        if(mc_status[iMC] == 23) {Ujet_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Ujet_Indice == -999) {Ujet_Indice = iMC;}
        foundUjet = true;
        
      }
      else if( mc_pdgId[iMC]==  11 && mc_mother[iMC]  == 23   ){
        if(mc_status[iMC] == 23) {Zel_min_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zel_min_Indice== -999) {Zel_min_Indice = iMC;}
        foundZelmin = true;
        //cout << "el from Z from tbar found  " << endl;
      } // el - from Z  from tbar
      else if( mc_pdgId[iMC]==  -11 && mc_mother[iMC]  == 23  ){
        if(mc_status[iMC] == 23) {Zel_plus_Indice= iMC;  }
        else if( mc_status[iMC] != 23 && Zel_plus_Indice== -999) {Zel_plus_Indice = iMC;}
        foundZelplus = true;
        
        // cout << "el from Z from tbar found  " << endl;
      } // el + from Z  from tbar
      else if( mc_pdgId[iMC]==  23 ){
        if(mc_status[iMC] == 23) {Z_Indice= iMC;  }
        else if( mc_status[iMC] != 23 && Z_Indice== -999) {Z_Indice = iMC; }
        foundZ = true;
      } //Z from top with mu decay
      
    }
    
    else if(dataSetName.find("TT")!=std::string::npos){
      
      //FCNC
      if( mc_pdgId[iMC]==  13 && mc_mother[iMC]  == 23 && abs(mc_granny[iMC])  == 6 ){
        if(mc_status[iMC] == 23) {Zmu_min_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zmu_min_Indice== -999) {Zmu_min_Indice = iMC;}
        foundZmumin = true;
        
        // cout << "mu from Z from tbar found  " << endl;
      } // mu - from Z  from tbar
      else if( mc_pdgId[iMC]==  -13 && mc_mother[iMC]  == 23  && abs(mc_granny[iMC])  == 6){
        if(mc_status[iMC] == 23) {Zmu_plus_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zmu_plus_Indice== -999) {Zmu_plus_Indice= iMC;}
        foundZmuplus = true;
        
        //cout << "mu from Z from tbar found  " << endl;
      } // mu + from Z  from tbar
      else if( abs( mc_pdgId[iMC]) ==  4 && abs(mc_mother[iMC])  == 6 ){
        if(mc_status[iMC] == 23) {Cjet_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Cjet_Indice == -999) {Cjet_Indice = iMC;}
        foundCjet = true;
        
      } // c jet
      else if( abs( mc_pdgId[iMC]) ==  2  && abs(mc_mother[iMC])  == 6){
        if(mc_status[iMC] == 23) {Ujet_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Ujet_Indice == -999) {Ujet_Indice = iMC;}
        foundUjet = true;
        
      }
      
      else if( mc_pdgId[iMC]==  11 && mc_mother[iMC]  == 23  && abs(mc_granny[iMC])  == 6 ){
        if(mc_status[iMC] == 23) {Zel_min_Indice = iMC;  }
        else if( mc_status[iMC] != 23 && Zel_min_Indice== -999) {Zel_min_Indice = iMC;}
        foundZelmin = true;
        //cout << "el from Z from tbar found  " << endl;
      } // el - from Z  from tbar
      else if( mc_pdgId[iMC]==  -11 && mc_mother[iMC]  == 23 && abs(mc_granny[iMC])  == 6 ){
        if(mc_status[iMC] == 23) {Zel_plus_Indice= iMC;  }
        else if( mc_status[iMC] != 23 && Zel_plus_Indice== -999) {Zel_plus_Indice = iMC;}
        foundZelplus = true;
        
        // cout << "el from Z from tbar found  " << endl;
      } // el + from Z  from tbar
      else if( mc_pdgId[iMC]==  23&& abs(mc_mother[iMC])  == 6 ){
        if(mc_status[iMC] == 23) {Z_Indice= iMC;  }
        else if( mc_status[iMC] != 23 && Z_Indice== -999) {Z_Indice = iMC; }
        foundZ = true;
      } //Z from top
      
      
    }
    
  }
  if(debug){
    cout << " ****** event searcher ******" << endl;
    cout << " foundZelmin " << foundZelmin <<" foundZelplus "<< foundZelplus <<" foundZmumin "<< foundZmumin <<" foundZmuplus "<< foundZmuplus <<endl;
    cout << "foundSMb "<< foundSMb <<" foundSMel "<< foundSMel <<" foundSMmu "<< foundSMmu << endl;
    
  }
  
}
Double_t RochLeptonMatching(TLorentzVector selectedlepton, vector <TLorentzVector> mcParticles, bool isData, double nbtracks, int chargelep, bool isNP, bool debug){
  if(debug) cout << "in rochester " << endl;
  
  mcParticlesroch.clear();
  if(!isData){
    for (int iMC = 0; iMC < nMCParticles; iMC++)
    {
      mcpart.Clear();
      mcpart.SetPtEtaPhiE(mc_pt[iMC], mc_eta[iMC],mc_phi[iMC],mc_E[iMC]);
      mcParticlesroch.push_back(mcpart);
    }
    if(mcParticlesroch.size() != nMCParticles){cout << "ERROR mcP roch not filled correctly" << endl;  }
    
  }
  
  if(debug) cout << "nb of mc particles " << mcParticlesroch.size() << endl;
  
  
  //vector <TLorentzVector> partonsrochester;
  partonsrochester.clear();
  double rocsf = 1.;
  
  
  if(isData){
    rocsf = rc.kScaleDT(chargelep, selectedlepton.Pt(), selectedlepton.Eta(), selectedlepton.Phi(), 0, 0);
    return rocsf;
  }
  else if(isNP){
    rocsf = rc.kScaleAndSmearMC(chargelep,selectedlepton.Pt(),selectedlepton.Eta(), selectedlepton.Phi(), nbtracks, gRandom->Rndm(), gRandom->Rndm(), 0, 0);
    return rocsf;
  }
  else{
    if(debug) cout << "loop over mc particles " << endl;
    for (unsigned int iMC = 0; iMC < mcParticlesroch.size(); iMC++)
    {
      if ( (mc_status[iMC] > 1 && mc_status[iMC] <= 20) || mc_status[iMC] >= 30 ) continue;  /// Final state particle or particle from hardest process
      if( abs(mc_pdgId[iMC]) ==  13  ){
        partonsrochester.push_back(mcParticlesroch[iMC]);
      } // muons
    }
    //vector <TLorentzVector> selectedleps;
    selectedleps.clear();
    selectedleps.push_back(selectedlepton);
    
    JetPartonMatching matchingToollep = JetPartonMatching(partonsrochester, selectedleps,2,true,true,0.1 );
    
    if(debug)  cout << "mu parton matching done " << endl;
    if(debug)  cout << matchingToollep.getNumberOfAvailableCombinations() << endl;
    
    if (matchingToollep.getNumberOfAvailableCombinations() != 1)
      cerr << "matching.getNumberOfAvailableCombinations() = " << matchingToollep.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
    
    
    /// Fill match in JetPartonPair;
    vector< pair<unsigned int, unsigned int> > lepPartonPair; // First one is jet number, second one is mcParticle number
    lepPartonPair.clear();
    
    
    for (unsigned int i = 0; i < partonsrochester.size(); i++)
    {
      int matchedlepNumber = matchingToollep.getMatchForParton(i, 0);
      //cout << "matchedlepnb " <<matchedlepNumber << endl;
      if (matchedlepNumber > -1)
        lepPartonPair.push_back( pair<unsigned int, unsigned int> (matchedlepNumber, i) );
      // matched lep number is nb in selectedlep collection
      // i is nb in partons
      // partonID contains place in mcParticles vector
    }
    if(debug) cout << "pushing back partons" << endl;
    if(debug) cout << "nb of leptonparton pairs " << lepPartonPair.size() << endl;
    if(lepPartonPair.size()>0){
      for (unsigned int i = 0; i < lepPartonPair.size(); i++)
      {
        unsigned int ipart = lepPartonPair[i].second; // place in mcParticles vector
        // unsigned int particlenb = lepPartonPair[i].first;  // place in selectedLeptons vector
        
        double deltaR = sqrt(pow(selectedlepton.Eta()-partonsrochester[ipart].Eta(),2) + pow(selectedlepton.Phi()-partonsrochester[ipart].Phi(),2));
        if(deltaR<0.1){
          rocsf = rc.kScaleFromGenMC(chargelep,selectedlepton.Pt(),selectedlepton.Eta(), selectedlepton.Phi(), nbtracks, partonsrochester[ipart].Pt(), gRandom->Rndm(), 0, 0);
        }
        else {
          rocsf = rc.kScaleAndSmearMC(chargelep,selectedlepton.Pt(),selectedlepton.Eta(), selectedlepton.Phi(), nbtracks, gRandom->Rndm(), gRandom->Rndm(), 0, 0);
        }
      }
    }
    else{
      if(debug) cout << "didn't found match" << endl;
      rocsf = rc.kScaleAndSmearMC(chargelep,selectedlepton.Pt(),selectedlepton.Eta(), selectedlepton.Phi(), nbtracks, gRandom->Rndm(), gRandom->Rndm(), 0, 0);
    }
  }
  
  
  
  if(debug) cout << "in rochester sf is: "  << rocsf << endl;
  return rocsf;
  
}




