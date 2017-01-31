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
   float TopHadPt, TopHadEta, TopHadPhi, TopHadE, TopHadM;
   float TopHadBJetPt, TopHadBJetEta, TopHadBJetPhi, TopHadBJetE;
   float TopLepWPt, TopLepWEta, TopLepWPhi, TopLepWE, TopLepWM;
   float TopHadWPt, TopHadWEta, TopHadWPhi, TopHadWE, TopHadWM;
   float TopHadWNonBJet1Pt, TopHadWNonBJet1Eta, TopHadWNonBJet1Phi, TopHadWNonBJet1E;
   float TopHadWNonBJet2Pt, TopHadWNonBJet2Eta, TopHadWNonBJet2Phi, TopHadWNonBJet2E;
   // Reconstruction
   float MetRecPx, MetRecPy;
   float TopLepWLepRecPt, TopLepWLepRecEta, TopLepWLepRecPhi, TopLepWLepRecE;
   float TopLepRecPt, TopLepRecEta, TopLepRecPhi, TopLepRecE, TopLepRecM;
   float TopLepBJetRecPt, TopLepBJetRecEta, TopLepBJetRecPhi, TopLepBJetRecE, TopLepBJetRecCSVv2;
   float TopHadRecPt, TopHadRecEta, TopHadRecPhi, TopHadRecE, TopHadRecM;
   float TopHadBJetRecPt, TopHadBJetRecEta, TopHadBJetRecPhi, TopHadBJetRecE, TopHadBJetRecCSVv2;
   float TopHadWRecPt, TopHadWRecEta, TopHadWRecPhi, TopHadWRecE, TopHadWRecM;
   float TopHadWNonBJet1RecPt, TopHadWNonBJet1RecEta, TopHadWNonBJet1RecPhi, TopHadWNonBJet1RecE, TopHadWNonBJet1RecCSVv2;
   float TopHadWNonBJet2RecPt, TopHadWNonBJet2RecEta, TopHadWNonBJet2RecPhi, TopHadWNonBJet2RecE, TopHadWNonBJet2RecCSVv2;
   std::vector<float> OtherJetRecPt, OtherJetRecEta, OtherJetRecPhi, OtherJetRecE, OtherJetRecCSVv2;
   // Transfer functions
   float dMetPx, dMetPy;
   float dTopLepBJetPx, dTopLepBJetPy, dTopLepBJetPz, dTopLepBJetE;
   float dTopHadWNonBJet1Px, dTopHadWNonBJet1Py, dTopHadWNonBJet1Pz, dTopHadWNonBJet1E;
   float dTopHadWNonBJet2Px, dTopHadWNonBJet2Py, dTopHadWNonBJet2Pz, dTopHadWNonBJet2E;
   float dTopHadBJetPx, dTopHadBJetPy, dTopHadBJetPz, dTopHadBJetE;
   float dElecPx, dElecPy, dElecPz, dElecE;
   float dMuonPx, dMuonPy, dMuonPz, dMuonE;
   // Additional kinematic variables
   float TopHadWNonBJet1TopHadWNonBJet2Dr, TopHadWTopLepWNuDr, TopHadWTopLepWLepDr, TopHadWTopHadBJetDr, TopLepWLepTopHadBJetDr;
   float TopLepTopHadDr, TopLepTopHadWDr, TopHadTopHadWDr, TopLepWTopLepBJetDr;
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
   trGEN->Branch("TopHadPt",&TopHadPt,"TopHadPt/F");
   trGEN->Branch("TopHadEta",&TopHadEta,"TopHadEta/F");
   trGEN->Branch("TopHadPhi",&TopHadPhi,"TopHadPhi/F");
   trGEN->Branch("TopHadE",&TopHadE,"TopHadE/F");
   trGEN->Branch("TopHadM",&TopHadM,"TopHadM/F");
   trGEN->Branch("TopHadBJetPt",&TopHadBJetPt,"TopHadBJetPt/F");
   trGEN->Branch("TopHadBJetEta",&TopHadBJetEta,"TopHadBJetEta/F");
   trGEN->Branch("TopHadBJetPhi",&TopHadBJetPhi,"TopHadBJetPhi/F");
   trGEN->Branch("TopHadBJetE",&TopHadBJetE,"TopHadBJetE/F");
   trGEN->Branch("TopLepWPt",&TopLepWPt,"TopLepWPt/F");
   trGEN->Branch("TopLepWEta",&TopLepWEta,"TopLepWEta/F");
   trGEN->Branch("TopLepWPhi",&TopLepWPhi,"TopLepWPhi/F");
   trGEN->Branch("TopLepWE",&TopLepWE,"TopLepWE/F");
   trGEN->Branch("TopLepWM",&TopLepWM,"TopLepWM/F");
   trGEN->Branch("TopHadWPt",&TopHadWPt,"TopHadWPt/F");
   trGEN->Branch("TopHadWEta",&TopHadWEta,"TopHadWEta/F");
   trGEN->Branch("TopHadWPhi",&TopHadWPhi,"TopHadWPhi/F");
   trGEN->Branch("TopHadWE",&TopHadWE,"TopHadWE/F");
   trGEN->Branch("TopHadWM",&TopHadWM,"TopHadWM/F");
   trGEN->Branch("TopHadWNonBJet1Pt",&TopHadWNonBJet1Pt,"TopHadWNonBJet1Pt/F");
   trGEN->Branch("TopHadWNonBJet1Eta",&TopHadWNonBJet1Eta,"TopHadWNonBJet1Eta/F");
   trGEN->Branch("TopHadWNonBJet1Phi",&TopHadWNonBJet1Phi,"TopHadWNonBJet1Phi/F");
   trGEN->Branch("TopHadWNonBJet1E",&TopHadWNonBJet1E,"TopHadWNonBJet1E/F");
   trGEN->Branch("TopHadWNonBJet2Pt",&TopHadWNonBJet2Pt,"TopHadWNonBJet2Pt/F");
   trGEN->Branch("TopHadWNonBJet2Eta",&TopHadWNonBJet2Eta,"TopHadWNonBJet2Eta/F");
   trGEN->Branch("TopHadWNonBJet2Phi",&TopHadWNonBJet2Phi,"TopHadWNonBJet2Phi/F");
   trGEN->Branch("TopHadWNonBJet2E",&TopHadWNonBJet2E,"TopHadWNonBJet2E/F");

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
   trGEN->Branch("TopHadRecPt",&TopHadRecPt,"TopHadRecPt/F");
   trGEN->Branch("TopHadRecEta",&TopHadRecEta,"TopHadRecEta/F");
   trGEN->Branch("TopHadRecPhi",&TopHadRecPhi,"TopHadRecPhi/F");
   trGEN->Branch("TopHadRecE",&TopHadRecE,"TopHadRecE/F");
   trGEN->Branch("TopHadRecM",&TopHadRecM,"TopHadRecM/F");
   trGEN->Branch("TopHadBJetRecPt",&TopHadBJetRecPt,"TopHadBJetRecPt/F");
   trGEN->Branch("TopHadBJetRecEta",&TopHadBJetRecEta,"TopHadBJetRecEta/F");
   trGEN->Branch("TopHadBJetRecPhi",&TopHadBJetRecPhi,"TopHadBJetRecPhi/F");
   trGEN->Branch("TopHadBJetRecE",&TopHadBJetRecE,"TopHadBJetRecE/F");
   trGEN->Branch("TopHadBJetRecCSVv2",&TopHadBJetRecCSVv2,"TopHadBJetRecCSVv2/F");
   trGEN->Branch("TopHadWRecPt",&TopHadWRecPt,"TopHadWRecPt/F");
   trGEN->Branch("TopHadWRecEta",&TopHadWRecEta,"TopHadWRecEta/F");
   trGEN->Branch("TopHadWRecPhi",&TopHadWRecPhi,"TopHadWRecPhi/F");
   trGEN->Branch("TopHadWRecE",&TopHadWRecE,"TopHadWRecE/F");
   trGEN->Branch("TopHadWRecM",&TopHadWRecM,"TopHadWRecM/F");
   trGEN->Branch("TopHadWNonBJet1RecPt",&TopHadWNonBJet1RecPt,"TopHadWNonBJet1RecPt/F");
   trGEN->Branch("TopHadWNonBJet1RecEta",&TopHadWNonBJet1RecEta,"TopHadWNonBJet1RecEta/F");
   trGEN->Branch("TopHadWNonBJet1RecPhi",&TopHadWNonBJet1RecPhi,"TopHadWNonBJet1RecPhi/F");
   trGEN->Branch("TopHadWNonBJet1RecE",&TopHadWNonBJet1RecE,"TopHadWNonBJet1RecE/F");
   trGEN->Branch("TopHadWNonBJet1RecCSVv2",&TopHadWNonBJet1RecCSVv2,"TopHadWNonBJet1RecCSVv2/F");
   trGEN->Branch("TopHadWNonBJet2RecPt",&TopHadWNonBJet2RecPt,"TopHadWNonBJet2RecPt/F");
   trGEN->Branch("TopHadWNonBJet2RecEta",&TopHadWNonBJet2RecEta,"TopHadWNonBJet2RecEta/F");
   trGEN->Branch("TopHadWNonBJet2RecPhi",&TopHadWNonBJet2RecPhi,"TopHadWNonBJet2RecPhi/F");
   trGEN->Branch("TopHadWNonBJet2RecE",&TopHadWNonBJet2RecE,"TopHadWNonBJet2RecE/F");
   trGEN->Branch("TopHadWNonBJet2RecCSVv2",&TopHadWNonBJet2RecCSVv2,"TopHadWNonBJet2RecCSVv2/F");
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
   trGEN->Branch("dTopHadWNonBJet1Px",&dTopHadWNonBJet1Px,"dTopHadWNonBJet1Px/F");
   trGEN->Branch("dTopHadWNonBJet1Py",&dTopHadWNonBJet1Py,"dTopHadWNonBJet1Py/F");
   trGEN->Branch("dTopHadWNonBJet1Pz",&dTopHadWNonBJet1Pz,"dTopHadWNonBJet1Pz/F");
   trGEN->Branch("dTopHadWNonBJet1E",&dTopHadWNonBJet1E,"dTopHadWNonBJet1E/F"); 
   trGEN->Branch("dTopHadWNonBJet2Px",&dTopHadWNonBJet2Px,"dTopHadWNonBJet2Px/F");
   trGEN->Branch("dTopHadWNonBJet2Py",&dTopHadWNonBJet2Py,"dTopHadWNonBJet2Py/F");
   trGEN->Branch("dTopHadWNonBJet2Pz",&dTopHadWNonBJet2Pz,"dTopHadWNonBJet2Pz/F");
   trGEN->Branch("dTopHadWNonBJet2E",&dTopHadWNonBJet2E,"dTopHadWNonBJet2E/F"); 
   trGEN->Branch("dTopHadBJetPx",&dTopHadBJetPx,"dTopHadBJetPx/F");
   trGEN->Branch("dTopHadBJetPy",&dTopHadBJetPy,"dTopHadBJetPy/F");
   trGEN->Branch("dTopHadBJetPz",&dTopHadBJetPz,"dTopHadBJetPz/F");
   trGEN->Branch("dTopHadBJetE",&dTopHadBJetE,"dTopHadBJetE/F");
   trGEN->Branch("dElecPx",&dElecPx,"dElecPx/F");
   trGEN->Branch("dElecPy",&dElecPy,"dElecPy/F");
   trGEN->Branch("dElecPz",&dElecPz,"dElecPz/F");
   trGEN->Branch("dElecE",&dElecE,"dElecE/F");
   trGEN->Branch("dMuonPx",&dMuonPx,"dMuonPx/F");
   trGEN->Branch("dMuonPy",&dMuonPy,"dMuonPy/F");
   trGEN->Branch("dMuonPz",&dMuonPz,"dMuonPz/F");
   trGEN->Branch("dMuonE",&dMuonE,"dMuonE/F");

   // Additional kinematic variables
   
   trGEN->Branch("TopHadWNonBJet1TopHadWNonBJet2Dr",&TopHadWNonBJet1TopHadWNonBJet2Dr,"TopHadWNonBJet1TopHadWNonBJet2Dr/F");
   trGEN->Branch("TopHadWTopLepWNuDr",&TopHadWTopLepWNuDr,"TopHadWTopLepWNuDr/F");
   trGEN->Branch("TopHadWTopLepWLepDr",&TopHadWTopLepWLepDr,"TopHadWTopLepWLepDr/F");
   trGEN->Branch("TopHadWTopHadBJetDr",&TopHadWTopHadBJetDr,"TopHadWTopHadBJetDr/F");
   trGEN->Branch("TopLepWLepTopHadBJetDr",&TopLepWLepTopHadBJetDr,"TopLepWLepTopHadBJetDr/F");
   trGEN->Branch("TopLepTopHadDr",&TopLepTopHadDr,"TopLepTopHadDr/F");
   trGEN->Branch("TopLepTopHadWDr",&TopLepTopHadWDr,"TopLepTopHadWDr/F");
   trGEN->Branch("TopHadTopHadWDr",&TopHadTopHadWDr,"TopHadTopHadWDr/F");
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
	TopHadPt = NONE; TopHadEta = NONE; TopHadPhi = NONE; TopHadE = NONE; TopHadM = NONE;
	TopHadBJetPt = NONE; TopHadBJetEta = NONE; TopHadBJetPhi = NONE; TopHadBJetE = NONE;
	TopLepWPt = NONE; TopLepWEta = NONE; TopLepWPhi = NONE; TopLepWE = NONE; TopLepWM = NONE;
	TopHadWPt = NONE; TopHadWEta = NONE; TopHadWPhi = NONE; TopHadWE = NONE; TopHadWM = NONE;
	TopHadWNonBJet1Pt = NONE; TopHadWNonBJet1Eta = NONE; TopHadWNonBJet1Phi = NONE; TopHadWNonBJet1E = NONE;
	TopHadWNonBJet2Pt = NONE; TopHadWNonBJet2Eta = NONE; TopHadWNonBJet2Phi = NONE; TopHadWNonBJet2E = NONE;

	MetRecPx = NONE; MetRecPy = NONE;
	TopLepWLepRecPt = NONE; TopLepWLepRecEta = NONE; TopLepWLepRecPhi = NONE; TopLepWLepRecE = NONE;
	TopLepRecPt = NONE; TopLepRecEta = NONE; TopLepRecPhi = NONE; TopLepRecE = NONE; TopLepRecM = NONE;
	TopLepBJetRecPt = NONE; TopLepBJetRecEta = NONE; TopLepBJetRecPhi = NONE; TopLepBJetRecE = NONE; TopLepBJetRecCSVv2 = NONE;
	TopHadRecPt = NONE; TopHadRecEta = NONE; TopHadRecPhi = NONE; TopHadRecE = NONE; TopHadRecM = NONE;
	TopHadBJetRecPt = NONE; TopHadBJetRecEta = NONE; TopHadBJetRecPhi = NONE; TopHadBJetRecE = NONE; TopHadBJetRecCSVv2 = NONE;
	TopHadWRecPt = NONE; TopHadWRecEta = NONE; TopHadWRecPhi = NONE; TopHadWRecE = NONE; TopHadWRecM = NONE;
	TopHadWNonBJet1RecPt = NONE; TopHadWNonBJet1RecEta = NONE; TopHadWNonBJet1RecPhi = NONE; TopHadWNonBJet1RecE = NONE; TopHadWNonBJet1RecCSVv2 = NONE;
	TopHadWNonBJet2RecPt = NONE; TopHadWNonBJet2RecEta = NONE; TopHadWNonBJet2RecPhi = NONE; TopHadWNonBJet2RecE = NONE; TopHadWNonBJet2RecCSVv2 = NONE;
	OtherJetRecPt.clear(); OtherJetRecEta.clear(); OtherJetRecPhi.clear(); OtherJetRecE.clear(); OtherJetRecCSVv2.clear();
   
	dMetPx = NONE; dMetPy = NONE;
	dTopLepBJetPx = NONE; dTopLepBJetPy = NONE; dTopLepBJetPz = NONE; dTopLepBJetE = NONE;
	dTopHadWNonBJet1Px = NONE; dTopHadWNonBJet1Py = NONE; dTopHadWNonBJet1Pz = NONE; dTopHadWNonBJet1E = NONE;
	dTopHadWNonBJet2Px = NONE; dTopHadWNonBJet2Py = NONE; dTopHadWNonBJet2Pz = NONE; dTopHadWNonBJet2E = NONE;
	dTopHadBJetPx = NONE; dTopHadBJetPy = NONE; dTopHadBJetPz = NONE; dTopHadBJetE = NONE;
	dElecPx; dElecPy = NONE; dElecPz = NONE; dElecE = NONE;
	dMuonPx; dMuonPy = NONE; dMuonPz = NONE; dMuonE = NONE;

	TopHadWNonBJet1TopHadWNonBJet2Dr = NONE; TopHadWTopLepWNuDr = NONE; TopHadWTopLepWLepDr = NONE; TopHadWTopHadBJetDr = NONE; TopLepWLepTopHadBJetDr = NONE;
	TopLepTopHadDr = NONE; TopLepTopHadWDr = NONE; TopHadTopHadWDr = NONE; TopLepWTopLepBJetDr = NONE;
	TopLepWLepCharge = NONE;
	
	int Truth_idxTopLep = -1;
	int Truth_idxTopHad = -1;
	int Truth_idxTopLepBJet = -1;
	int Truth_idxTopHadBJet = -1;
	int Truth_idxTopLepW = -1;
	int Truth_idxTopLepWLep = -1;
	int Truth_idxTopLepWNu = -1;	
	int Truth_idxTopHadW = -1;
	int Truth_idxTopHadWNonBJet1 = -1;
	int Truth_idxTopHadWNonBJet2 = -1;
	
	tr->GetEntry(i);

	int NElec = electrons->size();
	int NMuon = muons->size();
	int NJet = jets->size();
	
	std::vector<int> plabel = truth->at(0).mc_truth_label();
	
	int nplabel = plabel.size();
	bool isTopLep2 = 0;
	for(int j=0;j<nplabel;j++)
	  {
	     int plab = plabel[j];
	     if( plab == 220 ) isTopLep2 = 1;
	  }	

	for(int j=0;j<nplabel;j++)
	  {
	     int plab = plabel[j];
	     
	     if( plab == 2 && isTopLep2 ) Truth_idxTopLep = j;
	     if( plab == 2 && !isTopLep2 ) Truth_idxTopHad = j;

	     if( plab == 3 && isTopLep2 ) Truth_idxTopHad = j;
	     if( plab == 3 && !isTopLep2 ) Truth_idxTopLep = j;
	     
	     if( plab == 21 && isTopLep2 ) Truth_idxTopLepBJet = j;
	     if( plab == 21 && !isTopLep2 ) Truth_idxTopHadBJet = j;

	     if( plab == 31 && isTopLep2 ) Truth_idxTopHadBJet = j;
	     if( plab == 31 && !isTopLep2 ) Truth_idxTopLepBJet = j;
	     
	     if( plab == 20 && isTopLep2 ) Truth_idxTopLepW = j;
	     if( plab == 20 && !isTopLep2 ) Truth_idxTopHadW = j;

	     if( plab == 30 && isTopLep2 ) Truth_idxTopHadW = j;
	     if( plab == 30 && !isTopLep2 ) Truth_idxTopLepW = j;
	     
	     if( plab == 220 || plab == 330 ) Truth_idxTopLepWLep = j;
	     if( plab == 221 || plab == 331 ) Truth_idxTopLepWNu = j;
	     
	     if( plab == 222 || plab == 332 ) Truth_idxTopHadWNonBJet1 = j;
	     if( plab == 223 || plab == 333 ) Truth_idxTopHadWNonBJet2 = j;
	  }

	if( Truth_idxTopLepWLep < 0 || Truth_idxTopHadWNonBJet1 < 0 ) continue;
	
	bool pass = (Truth_idxTopHadWNonBJet1 >= 0 && Truth_idxTopLepBJet >= 0 &&
		     Truth_idxTopHadBJet >= 0 && Truth_idxTopLepWLep >= 0);

	if( !pass )
	  {
	     std::cout << "These are not TopTopLepHad events" << std::endl;
	     
	     std::cout << "TopLep=" << Truth_idxTopLep << std::endl;
	     std::cout << "TopHad=" << Truth_idxTopHad << std::endl;
	     std::cout << "TopLepBJet=" << Truth_idxTopLepBJet << std::endl;
	     std::cout << "TopHadBJet=" << Truth_idxTopHadBJet << std::endl;
	     std::cout << "TopLepW=" << Truth_idxTopLepW << std::endl;
	     std::cout << "TopLepWLep=" << Truth_idxTopLepWLep << std::endl;
	     std::cout << "TopLepWNu=" << Truth_idxTopLepWNu << std::endl;
	     std::cout << "TopHadW=" << Truth_idxTopHadW << std::endl;
	     std::cout << "TopHadWNonBJet1=" << Truth_idxTopHadWNonBJet1 << std::endl;
	     std::cout << "TopHadWNonBJet2=" << Truth_idxTopHadWNonBJet2 << std::endl;
	     
	     exit(1);
	  }

	TopLepPt = truth->at(0).mc_truth_pt()[Truth_idxTopLep];
	TopLepEta = truth->at(0).mc_truth_eta()[Truth_idxTopLep];
	TopLepPhi = truth->at(0).mc_truth_phi()[Truth_idxTopLep];
	TopLepE = truth->at(0).mc_truth_E()[Truth_idxTopLep];

	TopHadPt = truth->at(0).mc_truth_pt()[Truth_idxTopHad];
	TopHadEta = truth->at(0).mc_truth_eta()[Truth_idxTopHad];
	TopHadPhi = truth->at(0).mc_truth_phi()[Truth_idxTopHad];
	TopHadE = truth->at(0).mc_truth_E()[Truth_idxTopHad];
	     
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
	
	TopHadBJetPt = truth->at(0).mc_truth_pt()[Truth_idxTopHadBJet];
	TopHadBJetEta = truth->at(0).mc_truth_eta()[Truth_idxTopHadBJet];
	TopHadBJetPhi = truth->at(0).mc_truth_phi()[Truth_idxTopHadBJet];
	TopHadBJetE = truth->at(0).mc_truth_E()[Truth_idxTopHadBJet];

	TopHadWPt = truth->at(0).mc_truth_pt()[Truth_idxTopHadW];
	TopHadWEta = truth->at(0).mc_truth_eta()[Truth_idxTopHadW];
	TopHadWPhi = truth->at(0).mc_truth_phi()[Truth_idxTopHadW];
	TopHadWE = truth->at(0).mc_truth_E()[Truth_idxTopHadW];	
	
	TopHadWNonBJet1Pt = truth->at(0).mc_truth_pt()[Truth_idxTopHadWNonBJet1];
	TopHadWNonBJet1Eta = truth->at(0).mc_truth_eta()[Truth_idxTopHadWNonBJet1];
	TopHadWNonBJet1Phi = truth->at(0).mc_truth_phi()[Truth_idxTopHadWNonBJet1];
	TopHadWNonBJet1E = truth->at(0).mc_truth_E()[Truth_idxTopHadWNonBJet1];
	
	TopHadWNonBJet2Pt = truth->at(0).mc_truth_pt()[Truth_idxTopHadWNonBJet2];
	TopHadWNonBJet2Eta = truth->at(0).mc_truth_eta()[Truth_idxTopHadWNonBJet2];
	TopHadWNonBJet2Phi = truth->at(0).mc_truth_phi()[Truth_idxTopHadWNonBJet2];
	TopHadWNonBJet2E = truth->at(0).mc_truth_E()[Truth_idxTopHadWNonBJet2];
	     
	TLorentzVector *Truth_TopLep = new TLorentzVector();
	Truth_TopLep->SetPtEtaPhiE(TopLepPt,TopLepEta,TopLepPhi,TopLepE);
	TopLepM = Truth_TopLep->M();
	
	TLorentzVector *Truth_TopHad = new TLorentzVector();
	Truth_TopHad->SetPtEtaPhiE(TopHadPt,TopHadEta,TopHadPhi,TopHadE);
	TopHadM = Truth_TopHad->M();
	
	TLorentzVector *Truth_TopLepW = new TLorentzVector();
	Truth_TopLepW->SetPtEtaPhiE(TopLepWPt,TopLepWEta,TopLepWPhi,TopLepWE);
	TopLepWM = Truth_TopLepW->M();
	
	TLorentzVector *Truth_TopLepWNu = new TLorentzVector();
	Truth_TopLepWNu->SetPtEtaPhiE(TopLepWNuPt,TopLepWNuEta,TopLepWNuPhi,TopLepWNuE);

	TLorentzVector *Truth_TopLepWLep = new TLorentzVector();
	Truth_TopLepWLep->SetPtEtaPhiE(TopLepWLepPt,TopLepWLepEta,TopLepWLepPhi,TopLepWLepE);

	TLorentzVector *Truth_TopLepBJet = new TLorentzVector();
	Truth_TopLepBJet->SetPtEtaPhiE(TopLepBJetPt,TopLepBJetEta,TopLepBJetPhi,TopLepBJetE);
	
	TLorentzVector *Truth_TopHadBJet = new TLorentzVector();
	Truth_TopHadBJet->SetPtEtaPhiE(TopHadBJetPt,TopHadBJetEta,TopHadBJetPhi,TopHadBJetE);

	TLorentzVector *Truth_TopHadW = new TLorentzVector();
	Truth_TopHadW->SetPtEtaPhiE(TopHadWPt,TopHadWEta,TopHadWPhi,TopHadWE);
	TopHadWM = Truth_TopHadW->M();

	TLorentzVector *Truth_TopHadWNonBJet1 = new TLorentzVector();
	Truth_TopHadWNonBJet1->SetPtEtaPhiE(TopHadWNonBJet1Pt,TopHadWNonBJet1Eta,TopHadWNonBJet1Phi,TopHadWNonBJet1E);
	
	TLorentzVector *Truth_TopHadWNonBJet2 = new TLorentzVector();
	Truth_TopHadWNonBJet2->SetPtEtaPhiE(TopHadWNonBJet2Pt,TopHadWNonBJet2Eta,TopHadWNonBJet2Phi,TopHadWNonBJet2E);
	      	
	TLorentzVector *Rec_TopLep = new TLorentzVector();
	TLorentzVector *Rec_TopHad = new TLorentzVector();
	TLorentzVector *Rec_TopLepWLep = new TLorentzVector();
	TLorentzVector *Rec_TopLepBJet = new TLorentzVector();
	TLorentzVector *Rec_TopHadBJet = new TLorentzVector();
	TLorentzVector *Rec_TopHadW = new TLorentzVector();
	TLorentzVector *Rec_TopHadWNonBJet1 = new TLorentzVector();
	TLorentzVector *Rec_TopHadWNonBJet2 = new TLorentzVector();
	
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

	// TopHadBJet
	     
	float drMinTopHadBJet = NONE;
	int idxMinTopHadBJet = -1;
	     
	for(int ij=0;ij<NJet;ij++)
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();
	     
	     float drTopHadBJet = GetDeltaR(eta,phi,TopHadBJetEta,TopHadBJetPhi);
	     if( drTopHadBJet < drMinTopHadBJet && drTopHadBJet < 0.4 && 
		 fabs(pt-TopHadBJetPt)/std::max(pt,TopHadBJetPt) < 0.8 )
	       {
		  drMinTopHadBJet = drTopHadBJet;
		  idxMinTopHadBJet = ij;
	       }		  
	  }	     
	     
	TopHadBJetRecPt = (idxMinTopHadBJet >= 0) ? jets->at(idxMinTopHadBJet).pt() : NONE;
	TopHadBJetRecEta = (idxMinTopHadBJet >= 0) ? jets->at(idxMinTopHadBJet).eta() : NONE;
	TopHadBJetRecPhi = (idxMinTopHadBJet >= 0) ? jets->at(idxMinTopHadBJet).phi() : NONE;
	TopHadBJetRecE = (idxMinTopHadBJet >= 0) ? jets->at(idxMinTopHadBJet).E() : NONE;
	TopHadBJetRecCSVv2 = (idxMinTopHadBJet >= 0) ? jets->at(idxMinTopHadBJet).CSVv2() : NONE;

	if( idxMinTopHadBJet >= 0 )
	  Rec_TopHadBJet->SetPtEtaPhiE(TopHadBJetRecPt,
				       TopHadBJetRecEta,
				       TopHadBJetRecPhi,
				       TopHadBJetRecE);
	
	dTopHadBJetPx = (idxMinTopHadBJet >= 0) ? Truth_TopHadBJet->Px()-Rec_TopHadBJet->Px() : NONE;
	dTopHadBJetPy = (idxMinTopHadBJet >= 0) ? Truth_TopHadBJet->Py()-Rec_TopHadBJet->Py() : NONE;
	dTopHadBJetPz = (idxMinTopHadBJet >= 0) ? Truth_TopHadBJet->Pz()-Rec_TopHadBJet->Pz() : NONE;
	dTopHadBJetE = (idxMinTopHadBJet >= 0) ? Truth_TopHadBJet->E()-Rec_TopHadBJet->E() : NONE;
	
	// TopHadWNonBJet1, TopHadWNonBJet2
	
	float drMinTopHadWNonBJet1 = NONE;
	int idxMinTopHadWNonBJet1 = -1;
	float drMinTopHadWNonBJet2 = NONE;
	int idxMinTopHadWNonBJet2 = -1;

	for(int ij=0;ij<NJet;ij++)
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();

	     float drTopHadWNonBJet1 = GetDeltaR(eta,phi,TopHadWNonBJet1Eta,TopHadWNonBJet1Phi);
	     if( drTopHadWNonBJet1 < drMinTopHadWNonBJet1 && drTopHadWNonBJet1 < 0.4 && 
		 fabs(pt-TopHadWNonBJet1Pt)/std::max(pt,TopHadWNonBJet1Pt) < 0.8 )
	       {
		  drMinTopHadWNonBJet1 = drTopHadWNonBJet1;
		  idxMinTopHadWNonBJet1 = ij;
	       }		  
	     float drTopHadWNonBJet2 = GetDeltaR(eta,phi,TopHadWNonBJet2Eta,TopHadWNonBJet2Phi);
	     if( drTopHadWNonBJet2 < drMinTopHadWNonBJet2 && drTopHadWNonBJet2 < 0.4 && 
		 fabs(pt-TopHadWNonBJet2Pt)/std::max(pt,TopHadWNonBJet2Pt) < 0.8 )
	       {
		  drMinTopHadWNonBJet2 = drTopHadWNonBJet2;
		  idxMinTopHadWNonBJet2 = ij;
	       }
	  }

	TopHadWNonBJet1RecPt = (idxMinTopHadWNonBJet1 >= 0) ? jets->at(idxMinTopHadWNonBJet1).pt() : NONE;
	TopHadWNonBJet1RecEta = (idxMinTopHadWNonBJet1 >= 0) ? jets->at(idxMinTopHadWNonBJet1).eta() : NONE;
	TopHadWNonBJet1RecPhi = (idxMinTopHadWNonBJet1 >= 0) ? jets->at(idxMinTopHadWNonBJet1).phi() : NONE;
	TopHadWNonBJet1RecE = (idxMinTopHadWNonBJet1 >= 0) ? jets->at(idxMinTopHadWNonBJet1).E() : NONE;
	TopHadWNonBJet1RecCSVv2 = (idxMinTopHadWNonBJet1 >= 0) ? jets->at(idxMinTopHadWNonBJet1).CSVv2() : NONE;

	TopHadWNonBJet2RecPt = (idxMinTopHadWNonBJet2 >= 0) ? jets->at(idxMinTopHadWNonBJet2).pt() : NONE;
	TopHadWNonBJet2RecEta = (idxMinTopHadWNonBJet2 >= 0) ? jets->at(idxMinTopHadWNonBJet2).eta() : NONE;
	TopHadWNonBJet2RecPhi = (idxMinTopHadWNonBJet2 >= 0) ? jets->at(idxMinTopHadWNonBJet2).phi() : NONE;
	TopHadWNonBJet2RecE = (idxMinTopHadWNonBJet2 >= 0) ? jets->at(idxMinTopHadWNonBJet2).E() : NONE;
	TopHadWNonBJet2RecCSVv2 = (idxMinTopHadWNonBJet2 >= 0) ? jets->at(idxMinTopHadWNonBJet2).CSVv2() : NONE;
	
	if( idxMinTopHadWNonBJet1 >= 0 )
	  Rec_TopHadWNonBJet1->SetPtEtaPhiE(TopHadWNonBJet1RecPt,
					    TopHadWNonBJet1RecEta,
					    TopHadWNonBJet1RecPhi,
					    TopHadWNonBJet1RecE);
	
	if( idxMinTopHadWNonBJet2 >= 0 )
	  Rec_TopHadWNonBJet2->SetPtEtaPhiE(TopHadWNonBJet2RecPt,
					    TopHadWNonBJet2RecEta,
					    TopHadWNonBJet2RecPhi,
					    TopHadWNonBJet2RecE);
	
	dTopHadWNonBJet1Px = (idxMinTopHadWNonBJet1 >= 0) ? Truth_TopHadWNonBJet1->Px()-Rec_TopHadWNonBJet1->Px() : NONE;
	dTopHadWNonBJet1Py = (idxMinTopHadWNonBJet1 >= 0) ? Truth_TopHadWNonBJet1->Py()-Rec_TopHadWNonBJet1->Py() : NONE;
	dTopHadWNonBJet1Pz = (idxMinTopHadWNonBJet1 >= 0) ? Truth_TopHadWNonBJet1->Pz()-Rec_TopHadWNonBJet1->Pz() : NONE;
	dTopHadWNonBJet1E = (idxMinTopHadWNonBJet1 >= 0) ? Truth_TopHadWNonBJet1->E()-Rec_TopHadWNonBJet1->E() : NONE;

	dTopHadWNonBJet2Px = (idxMinTopHadWNonBJet2 >= 0) ? Truth_TopHadWNonBJet2->Px()-Rec_TopHadWNonBJet2->Px() : NONE;
	dTopHadWNonBJet2Py = (idxMinTopHadWNonBJet2 >= 0) ? Truth_TopHadWNonBJet2->Py()-Rec_TopHadWNonBJet2->Py() : NONE;
	dTopHadWNonBJet2Pz = (idxMinTopHadWNonBJet2 >= 0) ? Truth_TopHadWNonBJet2->Pz()-Rec_TopHadWNonBJet2->Pz() : NONE;
	dTopHadWNonBJet2E = (idxMinTopHadWNonBJet2 >= 0) ? Truth_TopHadWNonBJet2->E()-Rec_TopHadWNonBJet2->E() : NONE;

	*Rec_TopHadW = *Rec_TopHadWNonBJet1+*Rec_TopHadWNonBJet2;

	if( idxMinTopHadWNonBJet1 >= 0 && idxMinTopHadWNonBJet2 >= 0 )
	  {		  
	     TopHadWRecPt = Rec_TopHadW->Pt();
	     TopHadWRecEta = Rec_TopHadW->PseudoRapidity();
	     TopHadWRecPhi = Rec_TopHadW->Phi();
	     TopHadWRecE = Rec_TopHadW->E();
	     TopHadWRecM = Rec_TopHadW->M();
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
	     if( drTopLepWLep < drMinTopLepWLep_Elec && drTopLepWLep < 0.1 && LisElecTruth &&
		 fabs(pt-TopLepWLepPt)/std::max(pt,TopLepWLepPt) < 0.33 )
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
	     if( drTopLepWLep < drMinTopLepWLep_Muon && drTopLepWLep < 0.1 && !LisElecTruth &&
		 fabs(pt-TopLepWLepPt)/std::max(pt,TopLepWLepPt) < 0.33 )
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
	  
	*Rec_TopHad = *Rec_TopHadW+*Rec_TopHadBJet;
	if( idxMinTopHadWNonBJet1 >= 0 && idxMinTopHadWNonBJet2 >= 0 && idxMinTopHadBJet >= 0 )
	  {	     
	     TopHadRecPt = Rec_TopHad->Pt();
	     TopHadRecEta = Rec_TopHad->PseudoRapidity();
	     TopHadRecPhi = Rec_TopHad->Phi();
	     TopHadRecE = Rec_TopHad->E();
	     TopHadRecM = Rec_TopHad->M();
	  }	

	// Additional kinematic variables
	TopHadWNonBJet1TopHadWNonBJet2Dr = GetDeltaR(TopHadWNonBJet1Eta,TopHadWNonBJet1Phi,TopHadWNonBJet2Eta,TopHadWNonBJet2Phi);
	TopHadWTopLepWNuDr = GetDeltaR(TopHadWEta,TopHadWPhi,TopLepWNuEta,TopLepWNuPhi);
	TopHadWTopLepWLepDr = GetDeltaR(TopHadWEta,TopHadWPhi,TopLepWLepEta,TopLepWLepPhi);
	TopHadWTopHadBJetDr = GetDeltaR(TopHadWEta,TopHadWPhi,TopHadBJetEta,TopHadBJetPhi);
	TopLepWLepTopHadBJetDr = GetDeltaR(TopLepWLepEta,TopLepWLepPhi,TopHadBJetEta,TopHadBJetPhi);
	TopLepTopHadDr = GetDeltaR(TopLepEta,TopLepPhi,TopHadEta,TopHadPhi);
	TopLepTopHadWDr = GetDeltaR(TopLepEta,TopLepPhi,TopHadWEta,TopHadWPhi);
	TopHadTopHadWDr = GetDeltaR(TopHadEta,TopHadPhi,TopHadWEta,TopHadWPhi);
        TopLepWTopLepBJetDr = GetDeltaR(TopLepWEta,TopLepWPhi,TopLepBJetEta,TopLepBJetPhi);

	for(int ij=0;ij<NJet;ij++)
	  {
	     if( !jets->at(ij).isTight() ) continue;
	     
	     float pt = jets->at(ij).pt();
	     float eta = jets->at(ij).eta();
	     float phi = jets->at(ij).phi();
	     float E = jets->at(ij).E();
	     float CSVv2 = jets->at(ij).CSVv2();

	     if( idxMinTopHadWNonBJet1 != ij && idxMinTopHadWNonBJet2 != ij &&
		 idxMinTopHadBJet != ij && idxMinTopLepBJet != ij )
	       {
		  OtherJetRecPt.push_back(pt);
		  OtherJetRecEta.push_back(eta);
		  OtherJetRecPhi.push_back(phi);
		  OtherJetRecE.push_back(E);
		  OtherJetRecCSVv2.push_back(CSVv2);
	       }
	  }	
	
	delete Truth_TopLep;
	delete Truth_TopHad;
	delete Truth_TopLepW;
	delete Truth_TopLepWNu;
	delete Truth_TopLepWLep;
	delete Truth_TopLepBJet;
	delete Truth_TopHadBJet;
	delete Truth_TopHadW;
	delete Truth_TopHadWNonBJet1;
	delete Truth_TopHadWNonBJet2;

	delete Rec_TopLep;
	delete Rec_TopHad;
	delete Rec_TopLepWLep;
	delete Rec_TopLepBJet;
	delete Rec_TopHadBJet;
	delete Rec_TopHadW;
	delete Rec_TopHadWNonBJet1;
	delete Rec_TopHadWNonBJet2;
	
	Pass = 
	  ( idxMinTopHadWNonBJet1 >= 0 && idxMinTopHadWNonBJet2 >= 0 &&
	    idxMinTopLepBJet >= 0 && idxMinTopHadBJet >= 0 &&
	    ((idxMinTopLepWLep_Elec >= 0 && LisElec) || (idxMinTopLepWLep_Muon >= 0 && !LisElec)) &&
	    idxMinTopHadWNonBJet1 != idxMinTopHadWNonBJet2 &&
	    idxMinTopHadWNonBJet1 != idxMinTopHadBJet &&
	    idxMinTopHadWNonBJet1 != idxMinTopLepBJet &&
	    idxMinTopLepBJet != idxMinTopHadBJet &&
	    idxMinTopLepBJet != idxMinTopHadWNonBJet2 &&
	    idxMinTopHadBJet != idxMinTopHadWNonBJet2 );
	  
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
