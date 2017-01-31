#include "include/Jet.h"
#include "include/Electron.h"
#include "include/Muon.h"
#include "include/Truth.h"
#include "include/Event.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

#define NONE 10E+8

float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
float GetDeltaPhi(float phi1,float phi2);
float GetDeltaEta(float eta1,float eta2);

int main(int argc, char *argv[])
{
   std::vector<Jet> *jets = new std::vector<Jet>();
   std::vector<Electron> *electrons = new std::vector<Electron>();
   std::vector<Muon> *muons = new std::vector<Muon>();
   std::vector<Truth> *truth = new std::vector<Truth>();
   std::vector<Event> *event = new std::vector<Event>();
   
   std::ifstream infile(argv[1]);
   std::string fname;
   if( infile.good() )
     {
	getline(infile,fname);
     }   
   infile.close();
   
   std::string foutName = std::string(argv[2])+".root";

   TFile *f = TFile::Open(fname.c_str());
   
   TTree *tr = (TTree*)f->Get("Nt");
   
   tr->SetBranchAddress("Truth",&truth);
   tr->SetBranchAddress("Jet",&jets);
   tr->SetBranchAddress("Event",&event);
   tr->SetBranchAddress("ElectronTight",&electrons);
   tr->SetBranchAddress("MuonTight",&muons);

   TFile *fout = new TFile(foutName.c_str(),"RECREATE");

   bool Pass;
   // Generator
   float TopLepWLepPt, TopLepWLepEta, TopLepWLepPhi, TopLepWLepE;
   int TopLepWLepId;
   float TopLepWNuPt, TopLepWNuEta, TopLepWNuPhi, TopLepWNuE;
   int TopLepWNuId;
   float TopLepPt, TopLepEta, TopLepPhi, TopLepE, TopLepM;
   float TopLepBJetPt, TopLepBJetEta, TopLepBJetPhi, TopLepBJetE;
   float TopLepWPt, TopLepWEta, TopLepWPhi, TopLepWE, TopLepWM;
   float HiggsPt, HiggsEta, HiggsPhi, HiggsE, HiggsM;
   float HiggsBJet1Pt, HiggsBJet1Eta, HiggsBJet1Phi, HiggsBJet1E;
   float HiggsBJet2Pt, HiggsBJet2Eta, HiggsBJet2Phi, HiggsBJet2E;
   // Reconstruction
   float MetRecPx, MetRecPy;
   float TopLepWLepRecPt, TopLepWLepRecEta, TopLepWLepRecPhi, TopLepWLepRecE;
   float TopLepRecPt, TopLepRecEta, TopLepRecPhi, TopLepRecE, TopLepRecM;
   float TopLepBJetRecPt, TopLepBJetRecEta, TopLepBJetRecPhi, TopLepBJetRecE, TopLepBJetRecCSVv2;
   float HiggsRecPt, HiggsRecEta, HiggsRecPhi, HiggsRecE, HiggsRecM;
   float HiggsBJet1RecPt, HiggsBJet1RecEta, HiggsBJet1RecPhi, HiggsBJet1RecE, HiggsBJet1RecCSVv2;
   float HiggsBJet2RecPt, HiggsBJet2RecEta, HiggsBJet2RecPhi, HiggsBJet2RecE, HiggsBJet2RecCSVv2;
   std::vector<float> OtherJetRecPt, OtherJetRecEta, OtherJetRecPhi, OtherJetRecE, OtherJetRecCSVv2;
   // Transfer functions
   float dMetPx, dMetPy;
   float dTopLepBJetPx, dTopLepBJetPy, dTopLepBJetPz, dTopLepBJetE;
   float dHiggsBJet1Px, dHiggsBJet1Py, dHiggsBJet1Pz, dHiggsBJet1E;
   float dHiggsBJet2Px, dHiggsBJet2Py, dHiggsBJet2Pz, dHiggsBJet2E;
   float dElecPx, dElecPy, dElecPz, dElecE;
   float dMuonPx, dMuonPy, dMuonPz, dMuonE;
   // Additional kinematic variables
   float HiggsBJet1HiggsBJet2Dr, HiggsTopLepWNuDr, HiggsTopLepWLepDr;
   float TopLepHiggsDr, TopLepWTopLepBJetDr;
   int TopLepWLepCharge;
   
   TTree *trGEN = new TTree("trGEN","trGEN");

   trGEN->Branch("Pass",&Pass,"Pass/O");
   
   // Generator
   trGEN->Branch("TopLepWLepPt",&TopLepWLepPt,"TopLepWLepPt/F");
   trGEN->Branch("TopLepWLepEta",&TopLepWLepEta,"TopLepWLepEta/F");
   trGEN->Branch("TopLepWLepPhi",&TopLepWLepPhi,"TopLepWLepPhi/F");
   trGEN->Branch("TopLepWLepE",&TopLepWLepE,"TopLepWLepE/F");
   trGEN->Branch("TopLepWLepId",&TopLepWLepId,"TopLepWLepId/I");
   trGEN->Branch("TopLepWNuPt",&TopLepWNuPt,"TopLepWNuPt/F");
   trGEN->Branch("TopLepWNuEta",&TopLepWNuEta,"TopLepWNuEta/F");
   trGEN->Branch("TopLepWNuPhi",&TopLepWNuPhi,"TopLepWNuPhi/F");
   trGEN->Branch("TopLepWNuE",&TopLepWNuE,"TopLepWNuE/F");
   trGEN->Branch("TopLepWNuId",&TopLepWNuId,"TopLepWNuId/I");
   trGEN->Branch("TopLepPt",&TopLepPt,"TopLepPt/F");
   trGEN->Branch("TopLepEta",&TopLepEta,"TopLepEta/F");
   trGEN->Branch("TopLepPhi",&TopLepPhi,"TopLepPhi/F");
   trGEN->Branch("TopLepE",&TopLepE,"TopLepE/F");
   trGEN->Branch("TopLepM",&TopLepM,"TopLepM/F");
   trGEN->Branch("TopLepBJetPt",&TopLepBJetPt,"TopLepBJetPt/F");
   trGEN->Branch("TopLepBJetEta",&TopLepBJetEta,"TopLepBJetEta/F");
   trGEN->Branch("TopLepBJetPhi",&TopLepBJetPhi,"TopLepBJetPhi/F");
   trGEN->Branch("TopLepBJetE",&TopLepBJetE,"TopLepBJetE/F");
   trGEN->Branch("TopLepWPt",&TopLepWPt,"TopLepWPt/F");
   trGEN->Branch("TopLepWEta",&TopLepWEta,"TopLepWEta/F");
   trGEN->Branch("TopLepWPhi",&TopLepWPhi,"TopLepWPhi/F");
   trGEN->Branch("TopLepWE",&TopLepWE,"TopLepWE/F");
   trGEN->Branch("TopLepWM",&TopLepWM,"TopLepWM/F");
   trGEN->Branch("HiggsPt",&HiggsPt,"HiggsPt/F");
   trGEN->Branch("HiggsEta",&HiggsEta,"HiggsEta/F");
   trGEN->Branch("HiggsPhi",&HiggsPhi,"HiggsPhi/F");
   trGEN->Branch("HiggsE",&HiggsE,"HiggsE/F");
   trGEN->Branch("HiggsM",&HiggsM,"HiggsM/F");
   trGEN->Branch("HiggsBJet1Pt",&HiggsBJet1Pt,"HiggsBJet1Pt/F");
   trGEN->Branch("HiggsBJet1Eta",&HiggsBJet1Eta,"HiggsBJet1Eta/F");
   trGEN->Branch("HiggsBJet1Phi",&HiggsBJet1Phi,"HiggsBJet1Phi/F");
   trGEN->Branch("HiggsBJet1E",&HiggsBJet1E,"HiggsBJet1E/F");
   trGEN->Branch("HiggsBJet2Pt",&HiggsBJet2Pt,"HiggsBJet2Pt/F");
   trGEN->Branch("HiggsBJet2Eta",&HiggsBJet2Eta,"HiggsBJet2Eta/F");
   trGEN->Branch("HiggsBJet2Phi",&HiggsBJet2Phi,"HiggsBJet2Phi/F");
   trGEN->Branch("HiggsBJet2E",&HiggsBJet2E,"HiggsBJet2E/F");

   // Reconstruction
   trGEN->Branch("MetRecPx",&MetRecPx,"MetRecPx/F");
   trGEN->Branch("MetRecPy",&MetRecPy,"MetRecPy/F");
   trGEN->Branch("TopLepWLepRecPt",&TopLepWLepRecPt,"TopLepWLepRecPt/F");
   trGEN->Branch("TopLepWLepRecEta",&TopLepWLepRecEta,"TopLepWLepRecEta/F");
   trGEN->Branch("TopLepWLepRecPhi",&TopLepWLepRecPhi,"TopLepWLepRecPhi/F");
   trGEN->Branch("TopLepWLepRecE",&TopLepWLepRecE,"TopLepWLepRecE/F");
   trGEN->Branch("TopLepRecPt",&TopLepRecPt,"TopLepRecPt/F");
   trGEN->Branch("TopLepRecEta",&TopLepRecEta,"TopLepRecEta/F");
   trGEN->Branch("TopLepRecPhi",&TopLepRecPhi,"TopLepRecPhi/F");
   trGEN->Branch("TopLepRecE",&TopLepRecE,"TopLepRecE/F");
   trGEN->Branch("TopLepRecM",&TopLepRecM,"TopLepRecM/F");
   trGEN->Branch("TopLepBJetRecPt",&TopLepBJetRecPt,"TopLepBJetRecPt/F");
   trGEN->Branch("TopLepBJetRecEta",&TopLepBJetRecEta,"TopLepBJetRecEta/F");
   trGEN->Branch("TopLepBJetRecPhi",&TopLepBJetRecPhi,"TopLepBJetRecPhi/F");
   trGEN->Branch("TopLepBJetRecE",&TopLepBJetRecE,"TopLepBJetRecE/F");
   trGEN->Branch("TopLepBJetRecCSVv2",&TopLepBJetRecCSVv2,"TopLepBJetRecCSVv2/F");
   trGEN->Branch("HiggsRecPt",&HiggsRecPt,"HiggsRecPt/F");
   trGEN->Branch("HiggsRecEta",&HiggsRecEta,"HiggsRecEta/F");
   trGEN->Branch("HiggsRecPhi",&HiggsRecPhi,"HiggsRecPhi/F");
   trGEN->Branch("HiggsRecE",&HiggsRecE,"HiggsRecE/F");
   trGEN->Branch("HiggsRecM",&HiggsRecM,"HiggsRecM/F");
   trGEN->Branch("HiggsBJet1RecPt",&HiggsBJet1RecPt,"HiggsBJet1RecPt/F");
   trGEN->Branch("HiggsBJet1RecEta",&HiggsBJet1RecEta,"HiggsBJet1RecEta/F");
   trGEN->Branch("HiggsBJet1RecPhi",&HiggsBJet1RecPhi,"HiggsBJet1RecPhi/F");
   trGEN->Branch("HiggsBJet1RecE",&HiggsBJet1RecE,"HiggsBJet1RecE/F");
   trGEN->Branch("HiggsBJet1RecCSVv2",&HiggsBJet1RecCSVv2,"HiggsBJet1RecCSVv2/F");
   trGEN->Branch("HiggsBJet2RecPt",&HiggsBJet2RecPt,"HiggsBJet2RecPt/F");
   trGEN->Branch("HiggsBJet2RecEta",&HiggsBJet2RecEta,"HiggsBJet2RecEta/F");
   trGEN->Branch("HiggsBJet2RecPhi",&HiggsBJet2RecPhi,"HiggsBJet2RecPhi/F");
   trGEN->Branch("HiggsBJet2RecE",&HiggsBJet2RecE,"HiggsBJet2RecE/F");
   trGEN->Branch("HiggsBJet2RecCSVv2",&HiggsBJet2RecCSVv2,"HiggsBJet2RecCSVv2/F");
   trGEN->Branch("OtherJetRecPt","std::vector<float>",&OtherJetRecPt);
   trGEN->Branch("OtherJetRecEta","std::vector<float>",&OtherJetRecEta);
   trGEN->Branch("OtherJetRecPhi","std::vector<float>",&OtherJetRecPhi);
   trGEN->Branch("OtherJetRecE","std::vector<float>",&OtherJetRecE);
   trGEN->Branch("OtherJetRecCSVv2","std::vector<float>",&OtherJetRecCSVv2);
   
   // Transfer functions
   trGEN->Branch("dMetPx",&dMetPx,"dMetPx/F");
   trGEN->Branch("dMetPy",&dMetPy,"dMetPy/F");
   trGEN->Branch("dTopLepBJetPx",&dTopLepBJetPx,"dTopLepBJetPx/F");
   trGEN->Branch("dTopLepBJetPy",&dTopLepBJetPy,"dTopLepBJetPy/F");
   trGEN->Branch("dTopLepBJetPz",&dTopLepBJetPz,"dTopLepBJetPz/F");
   trGEN->Branch("dTopLepBJetE",&dTopLepBJetE,"dTopLepBJetE/F");
   trGEN->Branch("dHiggsBJet1Px",&dHiggsBJet1Px,"dHiggsBJet1Px/F");
   trGEN->Branch("dHiggsBJet1Py",&dHiggsBJet1Py,"dHiggsBJet1Py/F");
   trGEN->Branch("dHiggsBJet1Pz",&dHiggsBJet1Pz,"dHiggsBJet1Pz/F");
   trGEN->Branch("dHiggsBJet1E",&dHiggsBJet1E,"dHiggsBJet1E/F"); 
   trGEN->Branch("dHiggsBJet2Px",&dHiggsBJet2Px,"dHiggsBJet2Px/F");
   trGEN->Branch("dHiggsBJet2Py",&dHiggsBJet2Py,"dHiggsBJet2Py/F");
   trGEN->Branch("dHiggsBJet2Pz",&dHiggsBJet2Pz,"dHiggsBJet2Pz/F");
   trGEN->Branch("dHiggsBJet2E",&dHiggsBJet2E,"dHiggsBJet2E/F"); 
   trGEN->Branch("dElecPx",&dElecPx,"dElecPx/F");
   trGEN->Branch("dElecPy",&dElecPy,"dElecPy/F");
   trGEN->Branch("dElecPz",&dElecPz,"dElecPz/F");
   trGEN->Branch("dElecE",&dElecE,"dElecE/F");
   trGEN->Branch("dMuonPx",&dMuonPx,"dMuonPx/F");
   trGEN->Branch("dMuonPy",&dMuonPy,"dMuonPy/F");
   trGEN->Branch("dMuonPz",&dMuonPz,"dMuonPz/F");
   trGEN->Branch("dMuonE",&dMuonE,"dMuonE/F");

   // Additional kinematic variables
   
   trGEN->Branch("HiggsBJet1HiggsBJet2Dr",&HiggsBJet1HiggsBJet2Dr,"HiggsBJet1HiggsBJet2Dr/F");
   trGEN->Branch("HiggsTopLepWNuDr",&HiggsTopLepWNuDr,"HiggsTopLepWNuDr/F");
   trGEN->Branch("HiggsTopLepWLepDr",&HiggsTopLepWLepDr,"HiggsTopLepWLepDr/F");
   trGEN->Branch("TopLepHiggsDr",&TopLepHiggsDr,"TopLepHiggsDr/F");
   trGEN->Branch("TopLepWTopLepBJetDr",&TopLepWTopLepBJetDr,"TopLepWTopLepBJetDr/F");
   trGEN->Branch("TopLepWLepCharge",&TopLepWLepCharge,"TopLepWLepCharge/I");
   
   int nev = tr->GetEntries();
   std::cout << "Total number of events = " << nev << std::endl;
   
   for(int i=0;i<nev;i++)
     {
	Pass = false;
	TopLepWLepPt = NONE; TopLepWLepEta = NONE; TopLepWLepPhi = NONE; TopLepWLepE = NONE;
	TopLepWLepId = NONE;
	TopLepWNuPt = NONE; TopLepWNuEta = NONE; TopLepWNuPhi = NONE; TopLepWNuE = NONE;
	TopLepWNuId = NONE;
	TopLepPt = NONE; TopLepEta = NONE; TopLepPhi = NONE; TopLepE = NONE; TopLepM = NONE;
	TopLepBJetPt = NONE; TopLepBJetEta = NONE; TopLepBJetPhi = NONE; TopLepBJetE = NONE;
	TopLepWPt = NONE; TopLepWEta = NONE; TopLepWPhi = NONE; TopLepWE = NONE; TopLepWM = NONE;
	HiggsPt = NONE; HiggsEta = NONE; HiggsPhi = NONE; HiggsE = NONE; HiggsM = NONE;
	HiggsBJet1Pt = NONE; HiggsBJet1Eta = NONE; HiggsBJet1Phi = NONE; HiggsBJet1E = NONE;
	HiggsBJet2Pt = NONE; HiggsBJet2Eta = NONE; HiggsBJet2Phi = NONE; HiggsBJet2E = NONE;   

	MetRecPx = NONE; MetRecPy = NONE;
	TopLepWLepRecPt = NONE; TopLepWLepRecEta = NONE; TopLepWLepRecPhi = NONE; TopLepWLepRecE = NONE;
	TopLepRecPt = NONE; TopLepRecEta = NONE; TopLepRecPhi = NONE; TopLepRecE = NONE; TopLepRecM = NONE;
	TopLepBJetRecPt = NONE; TopLepBJetRecEta = NONE; TopLepBJetRecPhi = NONE; TopLepBJetRecE = NONE; TopLepBJetRecCSVv2 = NONE;
	HiggsRecPt = NONE; HiggsRecEta = NONE; HiggsRecPhi = NONE; HiggsRecE = NONE; HiggsRecM = NONE;
	HiggsBJet1RecPt = NONE; HiggsBJet1RecEta = NONE; HiggsBJet1RecPhi = NONE; HiggsBJet1RecE = NONE; HiggsBJet1RecCSVv2 = NONE;
	HiggsBJet2RecPt = NONE; HiggsBJet2RecEta = NONE; HiggsBJet2RecPhi = NONE; HiggsBJet2RecE = NONE; HiggsBJet2RecCSVv2 = NONE;
	OtherJetRecPt.clear(); OtherJetRecEta.clear(); OtherJetRecPhi.clear(); OtherJetRecE.clear(); OtherJetRecCSVv2.clear();
   
	dMetPx = NONE; dMetPy = NONE;
	dTopLepBJetPx = NONE; dTopLepBJetPy = NONE; dTopLepBJetPz = NONE; dTopLepBJetE = NONE;
	dHiggsBJet1Px = NONE; dHiggsBJet1Py = NONE; dHiggsBJet1Pz = NONE; dHiggsBJet1E = NONE;
	dHiggsBJet2Px = NONE; dHiggsBJet2Py = NONE; dHiggsBJet2Pz = NONE; dHiggsBJet2E = NONE;
	dElecPx; dElecPy = NONE; dElecPz = NONE; dElecE = NONE;
	dMuonPx; dMuonPy = NONE; dMuonPz = NONE; dMuonE = NONE;

	HiggsBJet1HiggsBJet2Dr = NONE; HiggsTopLepWNuDr = NONE; HiggsTopLepWLepDr = NONE;
	TopLepHiggsDr = NONE; TopLepWTopLepBJetDr = NONE;
	TopLepWLepCharge = NONE;
	
	int Truth_idxTopLep = -1;
	int Truth_idxTopLepBJet = -1;
	int Truth_idxTopLepW = -1;
	int Truth_idxTopLepWLep = -1;
	int Truth_idxTopLepWNu = -1;	
	int Truth_idxHiggs = -1;
	int Truth_idxHiggsBJet1 = -1;
	int Truth_idxHiggsBJet2 = -1;
	
	tr->GetEntry(i);

	int NElec = electrons->size();
	int NMuon = muons->size();
	int NJet = jets->size();
	
	std::vector<int> plabel = truth->at(0).mc_truth_label();
	
	int nplabel = plabel.size();

	for(int j=0;j<nplabel;j++)
	  {
	     int plab = plabel[j];
	     
	     if( plab == 2 || plab == 3 ) Truth_idxTopLep = j;
	     
	     if( plab == 21 || plab == 31 ) Truth_idxTopLepBJet = j;
	     if( plab == 20 || plab == 30 ) Truth_idxTopLepW = j;
	     if( plab == 220 || plab == 330 ) Truth_idxTopLepWLep = j;
	     if( plab == 221 || plab == 331 ) Truth_idxTopLepWNu = j;
	     
	     if( plab == 1 ) Truth_idxHiggs = j;
	     if( plab == 10 ) Truth_idxHiggsBJet1 = j;
	     if( plab == 11 ) Truth_idxHiggsBJet2 = j;
	  }
	
	bool pass = (Truth_idxHiggsBJet1 >= 0 && Truth_idxTopLepBJet >= 0 &&
		     Truth_idxTopLepWLep >= 0);

	if( !pass )
	  {
	     std::cout << "These are not TopHLepbb events" << std::endl;
	     
	     std::cout << "TopLep=" << Truth_idxTopLep << std::endl;
	     std::cout << "TopLepBJet=" << Truth_idxTopLepBJet << std::endl;
	     std::cout << "TopLepW=" << Truth_idxTopLepW << std::endl;
	     std::cout << "TopLepWLep=" << Truth_idxTopLepWLep << std::endl;
	     std::cout << "TopLepWNu=" << Truth_idxTopLepWNu << std::endl;
	     std::cout << "Higgs=" << Truth_idxHiggs << std::endl;
	     std::cout << "HiggsBJet1=" << Truth_idxHiggsBJet1 << std::endl;
	     std::cout << "HiggsBJet2=" << Truth_idxHiggsBJet2 << std::endl;
	     
	     exit(1);
	  }	

	TopLepPt = truth->at(0).mc_truth_pt()[Truth_idxTopLep];
	TopLepEta = truth->at(0).mc_truth_eta()[Truth_idxTopLep];
	TopLepPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLep];
	TopLepE = truth->at(0).mc_truth_E()[Truth_idxTopLep];
	     
	TopLepWPt = truth->at(0).mc_truth_pt()[Truth_idxTopLepW];
	TopLepWEta = truth->at(0).mc_truth_eta()[Truth_idxTopLepW];
	TopLepWPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLepW];
	TopLepWE = truth->at(0).mc_truth_E()[Truth_idxTopLepW];

	TopLepWLepPt = truth->at(0).mc_truth_pt()[Truth_idxTopLepWLep];
	TopLepWLepEta = truth->at(0).mc_truth_eta()[Truth_idxTopLepWLep];
	TopLepWLepPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLepWLep];
	TopLepWLepE = truth->at(0).mc_truth_E()[Truth_idxTopLepWLep];
	TopLepWLepId = truth->at(0).mc_truth_id()[Truth_idxTopLepWLep];
	if( TopLepWLepId < 0 ) TopLepWLepCharge = 1;
	else TopLepWLepCharge = -1;

	TopLepWNuPt = truth->at(0).mc_truth_pt()[Truth_idxTopLepWNu];
	TopLepWNuEta = truth->at(0).mc_truth_eta()[Truth_idxTopLepWNu];
	TopLepWNuPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLepWNu];
	TopLepWNuE = truth->at(0).mc_truth_E()[Truth_idxTopLepWNu];
	TopLepWNuId = truth->at(0).mc_truth_id()[Truth_idxTopLepWNu];
	
	TopLepBJetPt = truth->at(0).mc_truth_pt()[Truth_idxTopLepBJet];
	TopLepBJetEta = truth->at(0).mc_truth_eta()[Truth_idxTopLepBJet];
	TopLepBJetPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLepBJet];
	TopLepBJetE = truth->at(0).mc_truth_E()[Truth_idxTopLepBJet];

	HiggsPt = truth->at(0).mc_truth_pt()[Truth_idxHiggs];
	HiggsEta = truth->at(0).mc_truth_eta()[Truth_idxHiggs];
	HiggsPhi = truth->at(0).mc_truth_phi()[Truth_idxHiggs];
	HiggsE = truth->at(0).mc_truth_E()[Truth_idxHiggs];	
	
	HiggsBJet1Pt = truth->at(0).mc_truth_pt()[Truth_idxHiggsBJet1];
	HiggsBJet1Eta = truth->at(0).mc_truth_eta()[Truth_idxHiggsBJet1];
	HiggsBJet1Phi = truth->at(0).mc_truth_phi()[Truth_idxHiggsBJet1];
	HiggsBJet1E = truth->at(0).mc_truth_E()[Truth_idxHiggsBJet1];
	
	HiggsBJet2Pt = truth->at(0).mc_truth_pt()[Truth_idxHiggsBJet2];
	HiggsBJet2Eta = truth->at(0).mc_truth_eta()[Truth_idxHiggsBJet2];
	HiggsBJet2Phi = truth->at(0).mc_truth_phi()[Truth_idxHiggsBJet2];
	HiggsBJet2E = truth->at(0).mc_truth_E()[Truth_idxHiggsBJet2];
	     
	TLorentzVector *Truth_TopLep = new TLorentzVector();
	Truth_TopLep->SetPtEtaPhiE(TopLepPt,TopLepEta,TopLepPhi,TopLepE);
	TopLepM = Truth_TopLep->M();

	TLorentzVector *Truth_TopLepW = new TLorentzVector();
	Truth_TopLepW->SetPtEtaPhiE(TopLepWPt,TopLepWEta,TopLepWPhi,TopLepWE);
	TopLepWM = Truth_TopLepW->M();
	
	TLorentzVector *Truth_TopLepWNu = new TLorentzVector();
	Truth_TopLepWNu->SetPtEtaPhiE(TopLepWNuPt,TopLepWNuEta,TopLepWNuPhi,TopLepWNuE);

	TLorentzVector *Truth_TopLepWLep = new TLorentzVector();
	Truth_TopLepWLep->SetPtEtaPhiE(TopLepWLepPt,TopLepWLepEta,TopLepWLepPhi,TopLepWLepE);

	TLorentzVector *Truth_TopLepBJet = new TLorentzVector();
	Truth_TopLepBJet->SetPtEtaPhiE(TopLepBJetPt,TopLepBJetEta,TopLepBJetPhi,TopLepBJetE);

	TLorentzVector *Truth_Higgs = new TLorentzVector();
	Truth_Higgs->SetPtEtaPhiE(HiggsPt,HiggsEta,HiggsPhi,HiggsE);
	HiggsM = Truth_Higgs->M();

	TLorentzVector *Truth_HiggsBJet1 = new TLorentzVector();
	Truth_HiggsBJet1->SetPtEtaPhiE(HiggsBJet1Pt,HiggsBJet1Eta,HiggsBJet1Phi,HiggsBJet1E);
	
	TLorentzVector *Truth_HiggsBJet2 = new TLorentzVector();
	Truth_HiggsBJet2->SetPtEtaPhiE(HiggsBJet2Pt,HiggsBJet2Eta,HiggsBJet2Phi,HiggsBJet2E);
	      	
	TLorentzVector *Rec_TopLep = new TLorentzVector();
	TLorentzVector *Rec_TopLepWLep = new TLorentzVector();
	TLorentzVector *Rec_TopLepBJet = new TLorentzVector();
	TLorentzVector *Rec_Higgs = new TLorentzVector();
	TLorentzVector *Rec_HiggsBJet1 = new TLorentzVector();
	TLorentzVector *Rec_HiggsBJet2 = new TLorentzVector();
	
	// MET
	
	TLorentzVector MET = *Truth_TopLepWNu;
	float TruthMETPt = MET.Pt();
	float TruthMETPhi = MET.Phi();
	float TruthMETPx = TruthMETPt*cos(TruthMETPhi);
	float TruthMETPy = TruthMETPt*sin(TruthMETPhi);
	     
	float RecMETPt = event->at(0).metpt();
	float RecMETPhi = event->at(0).metphi();
	MetRecPx = RecMETPt*cos(RecMETPhi);
        MetRecPy = RecMETPt*sin(RecMETPhi);
	     
	dMetPx = (TruthMETPx-MetRecPx);
	dMetPy = (TruthMETPy-MetRecPy);

	// TopLepBJet
	     
	float drMinTopLepBJet = NONE;
	int idxMinTopLepBJet = -1;
	     
	for(int ij=0;ij<NJet;ij++)
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();
	     
	     float drTopLepBJet = GetDeltaR(eta,phi,TopLepBJetEta,TopLepBJetPhi);
	     if( drTopLepBJet < drMinTopLepBJet && drTopLepBJet < 0.4 && 
		 fabs(pt-TopLepBJetPt)/std::max(pt,TopLepBJetPt) < 0.8 )
	       {
		  drMinTopLepBJet = drTopLepBJet;
		  idxMinTopLepBJet = ij;
	       }		  
	  }	     
	     
	TopLepBJetRecPt = (idxMinTopLepBJet >= 0) ? jets->at(idxMinTopLepBJet).pt() : NONE;
	TopLepBJetRecEta = (idxMinTopLepBJet >= 0) ? jets->at(idxMinTopLepBJet).eta() : NONE;
	TopLepBJetRecPhi = (idxMinTopLepBJet >= 0) ? jets->at(idxMinTopLepBJet).phi() : NONE;
	TopLepBJetRecE = (idxMinTopLepBJet >= 0) ? jets->at(idxMinTopLepBJet).E() : NONE;
	TopLepBJetRecCSVv2 = (idxMinTopLepBJet >= 0) ? jets->at(idxMinTopLepBJet).CSVv2() : NONE;
	
	if( idxMinTopLepBJet >= 0 )
	  Rec_TopLepBJet->SetPtEtaPhiE(TopLepBJetRecPt,
				       TopLepBJetRecEta,
				       TopLepBJetRecPhi,
				       TopLepBJetRecE);
	
	dTopLepBJetPx = (idxMinTopLepBJet >= 0) ? Truth_TopLepBJet->Px()-Rec_TopLepBJet->Px() : NONE;
	dTopLepBJetPy = (idxMinTopLepBJet >= 0) ? Truth_TopLepBJet->Py()-Rec_TopLepBJet->Py() : NONE;
	dTopLepBJetPz = (idxMinTopLepBJet >= 0) ? Truth_TopLepBJet->Pz()-Rec_TopLepBJet->Pz() : NONE;
	dTopLepBJetE = (idxMinTopLepBJet >= 0) ? Truth_TopLepBJet->E()-Rec_TopLepBJet->E() : NONE;

	// HiggsBJet1, HiggsBJet2
	
	float drMinHiggsBJet1 = NONE;
	int idxMinHiggsBJet1 = -1;
	float drMinHiggsBJet2 = NONE;
	int idxMinHiggsBJet2 = -1;

	for(int ij=0;ij<NJet;ij++)
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();
	     
	     float drHiggsBJet1 = GetDeltaR(eta,phi,HiggsBJet1Eta,HiggsBJet1Phi);
	     if( drHiggsBJet1 < drMinHiggsBJet1 && drHiggsBJet1 < 0.4 && 
		 fabs(pt-HiggsBJet1Pt)/std::max(pt,HiggsBJet1Pt) < 0.8 )
	       {
		  drMinHiggsBJet1 = drHiggsBJet1;
		  idxMinHiggsBJet1 = ij;
	       }		  
	     float drHiggsBJet2 = GetDeltaR(eta,phi,HiggsBJet2Eta,HiggsBJet2Phi);
	     if( drHiggsBJet2 < drMinHiggsBJet2 && drHiggsBJet2 < 0.4 && 
		 fabs(pt-HiggsBJet2Pt)/std::max(pt,HiggsBJet2Pt) < 0.8 )
	       {
		  drMinHiggsBJet2 = drHiggsBJet2;
		  idxMinHiggsBJet2 = ij;
	       }		  
	  }	     

	HiggsBJet1RecPt = (idxMinHiggsBJet1 >= 0) ? jets->at(idxMinHiggsBJet1).pt() : NONE;
	HiggsBJet1RecEta = (idxMinHiggsBJet1 >= 0) ? jets->at(idxMinHiggsBJet1).eta() : NONE;
	HiggsBJet1RecPhi = (idxMinHiggsBJet1 >= 0) ? jets->at(idxMinHiggsBJet1).phi() : NONE;
	HiggsBJet1RecE = (idxMinHiggsBJet1 >= 0) ? jets->at(idxMinHiggsBJet1).E() : NONE;
	HiggsBJet1RecCSVv2 = (idxMinHiggsBJet1 >= 0) ? jets->at(idxMinHiggsBJet1).CSVv2() : NONE;

	HiggsBJet2RecPt = (idxMinHiggsBJet2 >= 0) ? jets->at(idxMinHiggsBJet2).pt() : NONE;
	HiggsBJet2RecEta = (idxMinHiggsBJet2 >= 0) ? jets->at(idxMinHiggsBJet2).eta() : NONE;
	HiggsBJet2RecPhi = (idxMinHiggsBJet2 >= 0) ? jets->at(idxMinHiggsBJet2).phi() : NONE;
	HiggsBJet2RecE = (idxMinHiggsBJet2 >= 0) ? jets->at(idxMinHiggsBJet2).E() : NONE;
	HiggsBJet2RecCSVv2 = (idxMinHiggsBJet2 >= 0) ? jets->at(idxMinHiggsBJet2).CSVv2() : NONE;
	
	if( idxMinHiggsBJet1 >= 0 )
	  Rec_HiggsBJet1->SetPtEtaPhiE(HiggsBJet1RecPt,
				       HiggsBJet1RecEta,
				       HiggsBJet1RecPhi,
				       HiggsBJet1RecE);

	if( idxMinHiggsBJet2 >= 0 )
	  Rec_HiggsBJet2->SetPtEtaPhiE(HiggsBJet2RecPt,
				       HiggsBJet2RecEta,
				       HiggsBJet2RecPhi,
				       HiggsBJet2RecE);
	
	dHiggsBJet1Px = (idxMinHiggsBJet1 >= 0) ? Truth_HiggsBJet1->Px()-Rec_HiggsBJet1->Px() : NONE;
	dHiggsBJet1Py = (idxMinHiggsBJet1 >= 0) ? Truth_HiggsBJet1->Py()-Rec_HiggsBJet1->Py() : NONE;
	dHiggsBJet1Pz = (idxMinHiggsBJet1 >= 0) ? Truth_HiggsBJet1->Pz()-Rec_HiggsBJet1->Pz() : NONE;
	dHiggsBJet1E = (idxMinHiggsBJet1 >= 0) ? Truth_HiggsBJet1->E()-Rec_HiggsBJet1->E() : NONE;

	dHiggsBJet2Px = (idxMinHiggsBJet2 >= 0) ? Truth_HiggsBJet2->Px()-Rec_HiggsBJet2->Px() : NONE;
	dHiggsBJet2Py = (idxMinHiggsBJet2 >= 0) ? Truth_HiggsBJet2->Py()-Rec_HiggsBJet2->Py() : NONE;
	dHiggsBJet2Pz = (idxMinHiggsBJet2 >= 0) ? Truth_HiggsBJet2->Pz()-Rec_HiggsBJet2->Pz() : NONE;
	dHiggsBJet2E = (idxMinHiggsBJet2 >= 0) ? Truth_HiggsBJet2->E()-Rec_HiggsBJet2->E() : NONE;

	*Rec_Higgs = *Rec_HiggsBJet1+*Rec_HiggsBJet2;

	if( idxMinHiggsBJet1 >= 0 && idxMinHiggsBJet2 >= 0 )
	  {		  
	     HiggsRecPt = Rec_Higgs->Pt();
	     HiggsRecEta = Rec_Higgs->PseudoRapidity();
	     HiggsRecPhi = Rec_Higgs->Phi();
	     HiggsRecE = Rec_Higgs->E();
	     HiggsRecM = Rec_Higgs->M();
	  }

	// Leptons

	bool LisElecTruth = (abs(TopLepWLepId) == 11);
	
	float drMinTopLepWLep_Elec = NONE;
	int idxMinTopLepWLep_Elec = -1;
	float drMinTopLepWLep_Muon = NONE;
	int idxMinTopLepWLep_Muon = -1;
	     
	for(int ie=0;ie<NElec;ie++)
	  {
	     float pt = electrons->at(ie).pt();
	     float eta = electrons->at(ie).eta();
	     float phi = electrons->at(ie).phi();
		  
	     float drTopLepWLep = GetDeltaR(eta,phi,TopLepWLepEta,TopLepWLepPhi);
	     if( drTopLepWLep < drMinTopLepWLep_Elec && drTopLepWLep < 0.1 && fabs(pt/TopLepWLepPt-1.) < 0.5 && LisElecTruth )
	       {
		  drMinTopLepWLep_Elec = drTopLepWLep;
		  idxMinTopLepWLep_Elec = ie;
	       }		  
	  }	     
	     
	for(int im=0;im<NMuon;im++)
	  {
	     float pt = muons->at(im).pt();
	     float eta = muons->at(im).eta();
	     float phi = muons->at(im).phi();
	     
	     float drTopLepWLep = GetDeltaR(eta,phi,TopLepWLepEta,TopLepWLepPhi);
	     if( drTopLepWLep < drMinTopLepWLep_Muon && drTopLepWLep < 0.1 && fabs(pt/TopLepWLepPt-1.) < 0.5 && !LisElecTruth )
	       {
		  drMinTopLepWLep_Muon = drTopLepWLep;
		  idxMinTopLepWLep_Muon = im;
	       }		  
	  }	     
	
	bool LisElec = (drMinTopLepWLep_Elec < drMinTopLepWLep_Muon) ? 1 : 0;
	     
	if( LisElec )
	  {	     
	     TopLepWLepRecPt = (idxMinTopLepWLep_Elec >= 0) ? electrons->at(idxMinTopLepWLep_Elec).pt() : NONE;
	     TopLepWLepRecEta = (idxMinTopLepWLep_Elec >= 0) ? electrons->at(idxMinTopLepWLep_Elec).eta() : NONE;
	     TopLepWLepRecPhi = (idxMinTopLepWLep_Elec >= 0) ? electrons->at(idxMinTopLepWLep_Elec).phi() : NONE;
	     TopLepWLepRecE = (idxMinTopLepWLep_Elec >= 0) ? electrons->at(idxMinTopLepWLep_Elec).E() : NONE;
	  }	
	else
	  {	     
	     TopLepWLepRecPt = (idxMinTopLepWLep_Muon >= 0) ? muons->at(idxMinTopLepWLep_Muon).pt() : NONE;
	     TopLepWLepRecEta = (idxMinTopLepWLep_Muon >= 0) ? muons->at(idxMinTopLepWLep_Muon).eta() : NONE;
	     TopLepWLepRecPhi = (idxMinTopLepWLep_Muon >= 0) ? muons->at(idxMinTopLepWLep_Muon).phi() : NONE;
	     TopLepWLepRecE = (idxMinTopLepWLep_Muon >= 0) ? muons->at(idxMinTopLepWLep_Muon).E() : NONE;
	  }	
	
	if( (idxMinTopLepWLep_Elec >= 0 && LisElec) || (idxMinTopLepWLep_Muon >= 0 && !LisElec) )
	  Rec_TopLepWLep->SetPtEtaPhiE(TopLepWLepRecPt,
				       TopLepWLepRecEta,
				       TopLepWLepRecPhi,
				       TopLepWLepRecE);
	
	if( LisElec && LisElecTruth )
	  {
	     dElecPx = (idxMinTopLepWLep_Elec >= 0) ? Truth_TopLepWLep->Px()-Rec_TopLepWLep->Px() : NONE;
	     dElecPy = (idxMinTopLepWLep_Elec >= 0) ? Truth_TopLepWLep->Py()-Rec_TopLepWLep->Py() : NONE;
	     dElecPz = (idxMinTopLepWLep_Elec >= 0) ? Truth_TopLepWLep->Pz()-Rec_TopLepWLep->Pz() : NONE;
	     dElecE = (idxMinTopLepWLep_Elec >= 0) ? Truth_TopLepWLep->E()-Rec_TopLepWLep->E() : NONE;
	  }
	
	if( !LisElec && !LisElecTruth )
	  {		  
	     dMuonPx = (idxMinTopLepWLep_Muon >= 0) ? Truth_TopLepWLep->Px()-Rec_TopLepWLep->Px() : NONE;
	     dMuonPy = (idxMinTopLepWLep_Muon >= 0) ? Truth_TopLepWLep->Py()-Rec_TopLepWLep->Py() : NONE;
	     dMuonPz = (idxMinTopLepWLep_Muon >= 0) ? Truth_TopLepWLep->Pz()-Rec_TopLepWLep->Pz() : NONE;
	     dMuonE = (idxMinTopLepWLep_Muon >= 0) ? Truth_TopLepWLep->E()-Rec_TopLepWLep->E() : NONE;
	  }
	
	*Rec_TopLep = *Rec_TopLepWLep+*Truth_TopLepWNu+*Rec_TopLepBJet;
	if( idxMinTopLepBJet >= 0 && 
	    ((idxMinTopLepWLep_Elec >= 0 && LisElec) ||
		(idxMinTopLepWLep_Muon >= 0 && !LisElec))
	  )
	  {
	     TopLepRecPt = Rec_TopLep->Pt();
	     TopLepRecEta = Rec_TopLep->PseudoRapidity();
	     TopLepRecPhi = Rec_TopLep->Phi();
	     TopLepRecE = Rec_TopLep->E();
	     TopLepRecM = Rec_TopLep->M();
	  }	

	// Additional kinematic variables
	HiggsBJet1HiggsBJet2Dr = GetDeltaR(HiggsBJet1Eta,HiggsBJet1Phi,HiggsBJet2Eta,HiggsBJet2Phi);
	HiggsTopLepWNuDr = GetDeltaR(HiggsEta,HiggsPhi,TopLepWNuEta,TopLepWNuPhi);
	HiggsTopLepWLepDr = GetDeltaR(HiggsEta,HiggsPhi,TopLepWLepEta,TopLepWLepPhi);
	TopLepHiggsDr = GetDeltaR(TopLepEta,TopLepPhi,HiggsEta,HiggsPhi);
        TopLepWTopLepBJetDr = GetDeltaR(TopLepWEta,TopLepWPhi,TopLepBJetEta,TopLepBJetPhi);
	
	for(int ij=0;ij<NJet;ij++) 
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();
	     float E = jets->at(ij).E();
	     float CSVv2 = jets->at(ij).CSVv2();

	     if( idxMinTopLepBJet != ij &&
		 idxMinHiggsBJet1 != ij && idxMinHiggsBJet2 != ij )
	       {
		  OtherJetRecPt.push_back(pt);
		  OtherJetRecEta.push_back(eta);
		  OtherJetRecPhi.push_back(phi);
		  OtherJetRecE.push_back(E);
		  OtherJetRecCSVv2.push_back(CSVv2);
	       }	     
	  }	

	delete Truth_TopLep;
	delete Truth_TopLepW;
	delete Truth_TopLepWNu;
	delete Truth_TopLepWLep;
	delete Truth_TopLepBJet;
	delete Truth_Higgs;
	delete Truth_HiggsBJet1;
	delete Truth_HiggsBJet2;

	delete Rec_TopLep;
	delete Rec_TopLepWLep;
	delete Rec_TopLepBJet;
	delete Rec_Higgs;
	delete Rec_HiggsBJet1;
	delete Rec_HiggsBJet2;

	Pass = 	
	  ( idxMinHiggsBJet1 >= 0 && idxMinHiggsBJet2 >= 0 && 
	    idxMinTopLepBJet >= 0 &&
	    ((idxMinTopLepWLep_Elec >= 0 && LisElec) || (idxMinTopLepWLep_Muon >= 0 && !LisElec)) &&
	    idxMinHiggsBJet1 != idxMinHiggsBJet2 &&
	    idxMinHiggsBJet1 != idxMinTopLepBJet &&
	    idxMinTopLepBJet != idxMinHiggsBJet2 );
	
	trGEN->Fill();
     }   

   int nevGEN = trGEN->GetEntries();
   std::cout << "Saved events = " << nevGEN << std::endl;
   
   fout->Write();
   fout->Close();
   
   f->Close();
}

float GetDeltaPhi(float phi1,float phi2)
{
   float DeltaPhi = phi2 - phi1;
//   float DeltaPhi = TMath::Abs(phi2 - phi1);
//   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return DeltaPhi;
}

float GetDeltaEta(float eta1,float eta2)
{   
   return (eta2-eta1);
}

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{   
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
