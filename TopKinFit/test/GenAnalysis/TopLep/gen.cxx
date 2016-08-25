#include "include/Jet.h"
#include "include/Electron.h"
#include "include/Muon.h"
#include "include/Truth.h"
#include "include/Event.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

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
   tr->SetBranchAddress("Electron",&electrons);
   tr->SetBranchAddress("Muon",&muons);

   TFile *fout = new TFile(foutName.c_str(),"RECREATE");

   float TopWLPt;
   float TopWLEta;
   float Top1M, Top2M;
   float TopW1M, TopW2M;
   float dMetPx, dMetPy;
   float dBJet1Px, dBJet1Py, dBJet1Pz, dBJet1E;
   float dBJet2Px, dBJet2Py, dBJet2Pz, dBJet2E;
   float dNonBJet1Px, dNonBJet1Py, dNonBJet1Pz, dNonBJet1E;
   float dNonBJet2Px, dNonBJet2Py, dNonBJet2Pz, dNonBJet2E;
   float dElecPx, dElecPy, dElecPz, dElecE;
   float dMuonPx, dMuonPy, dMuonPz, dMuonE;
   
   float TopWHadRecM;
   float TopLepRecM, TopHadRecM;
   
   TTree *trGEN = new TTree("trGEN","trGEN");
   
   trGEN->Branch("TopWLPt",&TopWLPt,"TopWLPt/F");
   trGEN->Branch("TopWLEta",&TopWLEta,"TopWLEta/F");
   trGEN->Branch("Top1M",&Top1M,"Top1M/F");
   trGEN->Branch("Top2M",&Top2M,"Top2M/F");
   trGEN->Branch("TopLepRecM",&TopLepRecM,"TopLepRecM/F");
   trGEN->Branch("TopHadRecM",&TopHadRecM,"TopHadRecM/F");
   trGEN->Branch("TopWHadRecM",&TopWHadRecM,"TopWHadRecM/F");
   trGEN->Branch("TopW1M",&TopW1M,"TopW1M/F");
   trGEN->Branch("TopW2M",&TopW2M,"TopW2M/F");
   trGEN->Branch("dMetPx",&dMetPx,"dMetPx/F");
   trGEN->Branch("dMetPy",&dMetPy,"dMetPy/F");
   trGEN->Branch("dBJet1Px",&dBJet1Px,"dBJet1Px/F");
   trGEN->Branch("dBJet1Py",&dBJet1Py,"dBJet1Py/F");
   trGEN->Branch("dBJet1Pz",&dBJet1Pz,"dBJet1Pz/F");
   trGEN->Branch("dBJet1E",&dBJet1E,"dBJet1E/F");
   trGEN->Branch("dBJet2Px",&dBJet2Px,"dBJet2Px/F");
   trGEN->Branch("dBJet2Py",&dBJet2Py,"dBJet2Py/F");
   trGEN->Branch("dBJet2Pz",&dBJet2Pz,"dBJet2Pz/F");
   trGEN->Branch("dBJet2E",&dBJet2E,"dBJet2E/F"); 
   trGEN->Branch("dNonBJet1Px",&dNonBJet1Px,"dNonBJet1Px/F");
   trGEN->Branch("dNonBJet1Py",&dNonBJet1Py,"dNonBJet1Py/F");
   trGEN->Branch("dNonBJet1Pz",&dNonBJet1Pz,"dNonBJet1Pz/F");
   trGEN->Branch("dNonBJet1E",&dNonBJet1E,"dNonBJet1E/F");
   trGEN->Branch("dNonBJet2Px",&dNonBJet2Px,"dNonBJet2Px/F");
   trGEN->Branch("dNonBJet2Py",&dNonBJet2Py,"dNonBJet2Py/F");
   trGEN->Branch("dNonBJet2Pz",&dNonBJet2Pz,"dNonBJet2Pz/F");
   trGEN->Branch("dNonBJet2E",&dNonBJet2E,"dNonBJet2E/F");
   trGEN->Branch("dElecPx",&dElecPx,"dElecPx/F");
   trGEN->Branch("dElecPy",&dElecPy,"dElecPy/F");
   trGEN->Branch("dElecPz",&dElecPz,"dElecPz/F");
   trGEN->Branch("dElecE",&dElecE,"dElecE/F");
   trGEN->Branch("dMuonPx",&dMuonPx,"dMuonPx/F");
   trGEN->Branch("dMuonPy",&dMuonPy,"dMuonPy/F");
   trGEN->Branch("dMuonPz",&dMuonPz,"dMuonPz/F");
   trGEN->Branch("dMuonE",&dMuonE,"dMuonE/F");
   
   int nev = tr->GetEntries();
   std::cout << "Total number of events = " << nev << std::endl;
   
   for(int i=0;i<nev;i++)
     {
	int Truth_idxTop1 = -1;
	int Truth_idxTop2 = -1;
	int Truth_idxTopB1 = -1;
	int Truth_idxTopB2 = -1;
	int Truth_idxTopW1 = -1;
	int Truth_idxTopW2 = -1;
	int Truth_idxTopWL1 = -1;
	int Truth_idxTopWL2 = -1;
	int Truth_idxTopWNu1 = -1;
	int Truth_idxTopWNu2 = -1;
	int Truth_idxTopWJ11 = -1;
	int Truth_idxTopWJ21 = -1;
	int Truth_idxTopWJ12 = -1;
	int Truth_idxTopWJ22 = -1;
	
	tr->GetEntry(i);

	std::vector<int> plabel = truth->at(0).mc_truth_label();
	
	int nplabel = plabel.size();
	for(int j=0;j<nplabel;j++)
	  {
	     int plab = plabel[j];
	     
	     if( plab == 2 ) Truth_idxTop1 = j;
	     if( plab == 3 ) Truth_idxTop2 = j;
	     if( plab == 20 ) Truth_idxTopB1 = j;
	     if( plab == 30 ) Truth_idxTopB2 = j;
	     if( plab == 21 ) Truth_idxTopW1 = j;
	     if( plab == 31 ) Truth_idxTopW2 = j;
	     if( plab == 210 ) Truth_idxTopWL1 = j;
	     if( plab == 310 ) Truth_idxTopWL2 = j;
	     if( plab == 211 ) Truth_idxTopWNu1 = j;
	     if( plab == 311 ) Truth_idxTopWNu2 = j;
	     if( plab == 223 ) Truth_idxTopWJ11 = j;
	     if( plab == 224 ) Truth_idxTopWJ21 = j;
	     if( plab == 323 ) Truth_idxTopWJ12 = j;
	     if( plab == 324 ) Truth_idxTopWJ22 = j;
	  }
	
	bool pass = (Truth_idxTopB1 >= 0 && Truth_idxTopB2 >= 0 && 
		     ((Truth_idxTopWL1 >= 0 && !(Truth_idxTopWL2 >= 0)) || (!(Truth_idxTopWL1 >= 0) && Truth_idxTopWL2 >= 0)) &&
		     ((Truth_idxTopWJ11 >= 0 && Truth_idxTopWJ21 >= 0) ||
			 (Truth_idxTopWJ12 >= 0 && Truth_idxTopWJ22 >= 0)));	     
	
	if( pass )
	  {
	     float Truth_Top1_Pt = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTop1] : -666.;
	     float Truth_Top1_Eta = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTop1] : -666.;
	     float Truth_Top1_Phi = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTop1] : -666.;
	     float Truth_Top1_E = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTop1] : -666.;
	     int Truth_Top1_Id = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTop1] : -666.;

	     float Truth_Top2_Pt = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTop2] : -666.;
	     float Truth_Top2_Eta = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTop2] : -666.;
	     float Truth_Top2_Phi = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTop2] : -666.;
	     float Truth_Top2_E = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTop2] : -666.;
	     int Truth_Top2_Id = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTop2] : -666.;

	     float Truth_TopWL_Pt = -666.;
	     float Truth_TopWL_Eta = -666.;
	     float Truth_TopWL_Phi = -666.;
	     float Truth_TopWL_E = -666.;
	     int Truth_TopWL_Id = -666.;

	     float Truth_TopWNu_Pt = -666.;
	     float Truth_TopWNu_Eta = -666.;
	     float Truth_TopWNu_Phi = -666.;
	     float Truth_TopWNu_E = -666.;
	     int Truth_TopWNu_Id = -666.;

	     float Truth_TopWJ1_Pt = -666.;
	     float Truth_TopWJ1_Eta = -666.;
	     float Truth_TopWJ1_Phi = -666.;
	     float Truth_TopWJ1_E = -666.;
	     int Truth_TopWJ1_Id = -666.;

	     float Truth_TopWJ2_Pt = -666.;
	     float Truth_TopWJ2_Eta = -666.;
	     float Truth_TopWJ2_Phi = -666.;
	     float Truth_TopWJ2_E = -666.;
	     int Truth_TopWJ2_Id = -666.;
	     
	     if( Truth_idxTopWL1 >= 0 )
	       {	
		  Truth_TopWL_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWL1];
		  Truth_TopWL_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWL1];
		  Truth_TopWL_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWL1];
		  Truth_TopWL_E = truth->at(0).mc_truth_E()[Truth_idxTopWL1];
		  Truth_TopWL_Id = truth->at(0).mc_truth_id()[Truth_idxTopWL1];
	       }
	     else if( Truth_idxTopWL2 >= 0 )
	       {
		  Truth_TopWL_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWL2];
		  Truth_TopWL_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWL2];
		  Truth_TopWL_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWL2];
		  Truth_TopWL_E = truth->at(0).mc_truth_E()[Truth_idxTopWL2];
		  Truth_TopWL_Id = truth->at(0).mc_truth_id()[Truth_idxTopWL2];
	       }

	     if( Truth_idxTopWNu1 >= 0 )
	       {		  
		  Truth_TopWNu_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWNu1];
		  Truth_TopWNu_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWNu1];
		  Truth_TopWNu_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWNu1];
		  Truth_TopWNu_E = truth->at(0).mc_truth_E()[Truth_idxTopWNu1];
		  Truth_TopWNu_Id = truth->at(0).mc_truth_id()[Truth_idxTopWNu1];
	       }
	     else if( Truth_idxTopWNu2 >= 0 )
	       {
		  Truth_TopWNu_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWNu2];
		  Truth_TopWNu_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWNu2];
		  Truth_TopWNu_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWNu2];
		  Truth_TopWNu_E = truth->at(0).mc_truth_E()[Truth_idxTopWNu2];
		  Truth_TopWNu_Id = truth->at(0).mc_truth_id()[Truth_idxTopWNu2];
	       }	     
	     
	     float Truth_TopW1_Pt = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_Eta = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_Phi = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_E = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopW1] : -666.;
	     int Truth_TopW1_Id = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopW1] : -666.;

	     float Truth_TopW2_Pt = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_Eta = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_Phi = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_E = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopW2] : -666.;
	     int Truth_TopW2_Id = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopW2] : -666.;
	     
	     float Truth_TopB1_Pt = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_Eta = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_Phi = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_E = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopB1] : -666.;
		  
	     float Truth_TopB2_Pt = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_Eta = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_Phi = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_E = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopB2] : -666.;

	     if( Truth_idxTopWJ11 >= 0 )
	       {		  
		  Truth_TopWJ1_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWJ11];
		  Truth_TopWJ1_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWJ11];
		  Truth_TopWJ1_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWJ11];
		  Truth_TopWJ1_E = truth->at(0).mc_truth_E()[Truth_idxTopWJ11];
		  Truth_TopWJ1_Id = truth->at(0).mc_truth_id()[Truth_idxTopWJ11];
	       }
	     else if( Truth_idxTopWJ12 >= 0 )
	       {
		  Truth_TopWJ1_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWJ12];
		  Truth_TopWJ1_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWJ12];
		  Truth_TopWJ1_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWJ12];
		  Truth_TopWJ1_E = truth->at(0).mc_truth_E()[Truth_idxTopWJ12];
		  Truth_TopWJ1_Id = truth->at(0).mc_truth_id()[Truth_idxTopWJ12];
	       }	     	     	     

	     if( Truth_idxTopWJ21 >= 0 )
	       {		  
		  Truth_TopWJ2_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWJ21];
		  Truth_TopWJ2_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWJ21];
		  Truth_TopWJ2_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWJ21];
		  Truth_TopWJ2_E = truth->at(0).mc_truth_E()[Truth_idxTopWJ21];
		  Truth_TopWJ2_Id = truth->at(0).mc_truth_id()[Truth_idxTopWJ21];
	       }
	     else if( Truth_idxTopWJ22 >= 0 )
	       {
		  Truth_TopWJ2_Pt = truth->at(0).mc_truth_pt()[Truth_idxTopWJ22];
		  Truth_TopWJ2_Eta = truth->at(0).mc_truth_eta()[Truth_idxTopWJ22];
		  Truth_TopWJ2_Phi = truth->at(0).mc_truth_phi()[Truth_idxTopWJ22];
		  Truth_TopWJ2_E = truth->at(0).mc_truth_E()[Truth_idxTopWJ22];
		  Truth_TopWJ2_Id = truth->at(0).mc_truth_id()[Truth_idxTopWJ22];
	       }	     	     	     
	     
	     TLorentzVector *Top1 = new TLorentzVector();
	     Top1->SetPtEtaPhiE(Truth_Top1_Pt,Truth_Top1_Eta,Truth_Top1_Phi,Truth_Top1_E);

	     TLorentzVector *Top2 = new TLorentzVector();
	     Top2->SetPtEtaPhiE(Truth_Top2_Pt,Truth_Top2_Eta,Truth_Top2_Phi,Truth_Top2_E);

	     TLorentzVector *TopW1 = new TLorentzVector();
	     TopW1->SetPtEtaPhiE(Truth_TopW1_Pt,Truth_TopW1_Eta,Truth_TopW1_Phi,Truth_TopW1_E);

	     TLorentzVector *TopW2 = new TLorentzVector();
	     TopW2->SetPtEtaPhiE(Truth_TopW2_Pt,Truth_TopW2_Eta,Truth_TopW2_Phi,Truth_TopW2_E);

	     TLorentzVector *TopWJ1 = new TLorentzVector();
	     TopWJ1->SetPtEtaPhiE(Truth_TopWJ1_Pt,Truth_TopWJ1_Eta,Truth_TopWJ1_Phi,Truth_TopWJ1_E);

	     TLorentzVector *TopWJ2 = new TLorentzVector();
	     TopWJ2->SetPtEtaPhiE(Truth_TopWJ2_Pt,Truth_TopWJ2_Eta,Truth_TopWJ2_Phi,Truth_TopWJ2_E);
	     
	     TLorentzVector *TopWNu = new TLorentzVector();
	     TopWNu->SetPtEtaPhiE(Truth_TopWNu_Pt,Truth_TopWNu_Eta,Truth_TopWNu_Phi,Truth_TopWNu_E);

	     TLorentzVector *TopWL = new TLorentzVector();
	     TopWL->SetPtEtaPhiE(Truth_TopWL_Pt,Truth_TopWL_Eta,Truth_TopWL_Phi,Truth_TopWL_E);

	     TLorentzVector *TopB1 = new TLorentzVector();
	     TopB1->SetPtEtaPhiE(Truth_TopB1_Pt,Truth_TopB1_Eta,Truth_TopB1_Phi,Truth_TopB1_E);

	     TLorentzVector *TopB2 = new TLorentzVector();
	     TopB2->SetPtEtaPhiE(Truth_TopB2_Pt,Truth_TopB2_Eta,Truth_TopB2_Phi,Truth_TopB2_E);
	     
	     TLorentzVector MET = *TopWNu;
	     float METPt = MET.Pt();
	     float METPhi = MET.Phi();
	     float METPx = METPt*cos(METPhi);
	     float METPy = METPt*sin(METPhi);
	     
	     float metPt = event->at(0).metpt();
	     float metPhi = event->at(0).metphi();
	     float metPx = metPt*cos(metPhi);
	     float metPy = metPt*sin(metPhi);
	     
	     dMetPx = (METPx-metPx);
	     dMetPy = (METPy-metPy);

	     // match b jets
	     
	     float drMinTopB1 = 10E+10;
	     int idxMinTopB1 = -1;
	     float drMinTopB2 = 10E+10;
	     int idxMinTopB2 = -1;
	     
	     int NJets = jets->size();
	     for(int ij=0;ij<NJets;ij++)
	       {
		  float pt = jets->at(ij).pt();
		  float eta = jets->at(ij).eta();
		  float phi = jets->at(ij).phi();

		  float drTopB1 = GetDeltaR(eta,phi,Truth_TopB1_Eta,Truth_TopB1_Phi);
		  if( drTopB1 < drMinTopB1 && drTopB1 < 0.4 && fabs(pt/Truth_TopB1_Pt-1.) < 0.5 )
		    {
		       drMinTopB1 = drTopB1;
		       idxMinTopB1 = ij;
		    }		  
		  float drTopB2 = GetDeltaR(eta,phi,Truth_TopB2_Eta,Truth_TopB2_Phi);
		  if( drTopB2 < drMinTopB2 && drTopB2 < 0.4 && fabs(pt/Truth_TopB2_Pt-1.) < 0.5 ) 
		    {
		       drMinTopB2 = drTopB2;
		       idxMinTopB2 = ij;
		    }		  
	       }	     
	     
	     TLorentzVector *TopB1Rec = new TLorentzVector();
	     TLorentzVector *TopB2Rec = new TLorentzVector();
	     
	     if( idxMinTopB1 >= 0 )
	       TopB1Rec->SetPtEtaPhiE(jets->at(idxMinTopB1).pt(),
				      jets->at(idxMinTopB1).eta(),
				      jets->at(idxMinTopB1).phi(),
				      jets->at(idxMinTopB1).E());

	     if( idxMinTopB2 >= 0 )
	       TopB2Rec->SetPtEtaPhiE(jets->at(idxMinTopB2).pt(),
				      jets->at(idxMinTopB2).eta(),
				      jets->at(idxMinTopB2).phi(),
				      jets->at(idxMinTopB2).E());
	     
	     dBJet1Px = (idxMinTopB1 >= 0) ? TopB1->Px()-TopB1Rec->Px() : 10E+10;
	     dBJet1Py = (idxMinTopB1 >= 0) ? TopB1->Py()-TopB1Rec->Py() : 10E+10;
	     dBJet1Pz = (idxMinTopB1 >= 0) ? TopB1->Pz()-TopB1Rec->Pz() : 10E+10;
	     dBJet1E = (idxMinTopB1 >= 0) ? TopB1->E()-TopB1Rec->E() : 10E+10;
	     
	     dBJet2Px = (idxMinTopB2 >= 0) ? TopB2->Px()-TopB2Rec->Px() : 10E+10;
	     dBJet2Py = (idxMinTopB2 >= 0) ? TopB2->Py()-TopB2Rec->Py() : 10E+10;
	     dBJet2Pz = (idxMinTopB2 >= 0) ? TopB2->Pz()-TopB2Rec->Pz() : 10E+10;
	     dBJet2E = (idxMinTopB2 >= 0) ? TopB2->E()-TopB2Rec->E() : 10E+10;

	     // match W jets
	     
	     float drMinTopWJ1 = 10E+10;
	     int idxMinTopWJ1 = -1;
	     float drMinTopWJ2 = 10E+10;
	     int idxMinTopWJ2 = -1;
	     
	     for(int ij=0;ij<NJets;ij++)
	       {
		  float pt = jets->at(ij).pt();
		  float eta = jets->at(ij).eta();
		  float phi = jets->at(ij).phi();

		  float drTopWJ1 = GetDeltaR(eta,phi,Truth_TopWJ1_Eta,Truth_TopWJ1_Phi);
		  if( drTopWJ1 < drMinTopWJ1 && drTopWJ1 < 0.4 && fabs(pt/Truth_TopWJ1_Pt-1.) < 0.5 )
		    {
		       drMinTopWJ1 = drTopWJ1;
		       idxMinTopWJ1 = ij;
		    }		  
		  float drTopWJ2 = GetDeltaR(eta,phi,Truth_TopWJ2_Eta,Truth_TopWJ2_Phi);
		  if( drTopWJ2 < drMinTopWJ2 && drTopWJ2 < 0.4 && fabs(pt/Truth_TopWJ2_Pt-1.) < 0.5 ) 
		    {
		       drMinTopWJ2 = drTopWJ2;
		       idxMinTopWJ2 = ij;
		    }		  
	       }	     
	     
	     TLorentzVector *TopWJ1Rec = new TLorentzVector();
	     TLorentzVector *TopWJ2Rec = new TLorentzVector();
	     
	     if( idxMinTopWJ1 >= 0 )
	       TopWJ1Rec->SetPtEtaPhiE(jets->at(idxMinTopWJ1).pt(),
				       jets->at(idxMinTopWJ1).eta(),
				       jets->at(idxMinTopWJ1).phi(),
				       jets->at(idxMinTopWJ1).E());

	     if( idxMinTopWJ2 >= 0 )
	       TopWJ2Rec->SetPtEtaPhiE(jets->at(idxMinTopWJ2).pt(),
				       jets->at(idxMinTopWJ2).eta(),
				       jets->at(idxMinTopWJ2).phi(),
				       jets->at(idxMinTopWJ2).E());
	     
	     dNonBJet1Px = (idxMinTopWJ1 >= 0) ? TopWJ1->Px()-TopWJ1Rec->Px() : 10E+10;
	     dNonBJet1Py = (idxMinTopWJ1 >= 0) ? TopWJ1->Py()-TopWJ1Rec->Py() : 10E+10;
	     dNonBJet1Pz = (idxMinTopWJ1 >= 0) ? TopWJ1->Pz()-TopWJ1Rec->Pz() : 10E+10;
	     dNonBJet1E = (idxMinTopWJ1 >= 0) ? TopWJ1->E()-TopWJ1Rec->E() : 10E+10;
	     
	     dNonBJet2Px = (idxMinTopWJ2 >= 0) ? TopWJ2->Px()-TopWJ2Rec->Px() : 10E+10;
	     dNonBJet2Py = (idxMinTopWJ2 >= 0) ? TopWJ2->Py()-TopWJ2Rec->Py() : 10E+10;
	     dNonBJet2Pz = (idxMinTopWJ2 >= 0) ? TopWJ2->Pz()-TopWJ2Rec->Pz() : 10E+10;
	     dNonBJet2E = (idxMinTopWJ2 >= 0) ? TopWJ2->E()-TopWJ2Rec->E() : 10E+10;

	     TopWHadRecM = 10E+10;

	     TLorentzVector TopWHadRec = *TopWJ1Rec+*TopWJ2Rec;
	     TLorentzVector TopHadRec;
	     
	     if( idxMinTopWJ1 >= 0 && idxMinTopWJ2 >= 0 )
	       {		  
		  TopWHadRecM = TopWHadRec.M();
	       }	     

	     delete TopWJ1Rec;
	     delete TopWJ2Rec;
	     
	     // match leptons

	     bool LisElecTruth = (abs(Truth_TopWL_Id) == 11);
	     
	     float drMinTopWL_Elec = 10E+10;
	     int idxMinTopWL_Elec = -1;
	     float drMinTopWL_Muon = 10E+10;
	     int idxMinTopWL_Muon = -1;
	     
	     int NElec = electrons->size();
	     int NMuon = muons->size();
	     
	     for(int ie=0;ie<NElec;ie++)
	       {
		  float pt = electrons->at(ie).pt();
		  float eta = electrons->at(ie).eta();
		  float phi = electrons->at(ie).phi();
		  
		  float drTopWL = GetDeltaR(eta,phi,Truth_TopWL_Eta,Truth_TopWL_Phi);
		  if( drTopWL < drMinTopWL_Elec && drTopWL < 0.1 && fabs(pt/Truth_TopWL_Pt-1.) < 0.5 && LisElecTruth )
		    {
		       drMinTopWL_Elec = drTopWL;
		       idxMinTopWL_Elec = ie;
		    }		  
	       }	     
	     
	     for(int im=0;im<NMuon;im++)
	       {
		  float pt = muons->at(im).pt();
		  float eta = muons->at(im).eta();
		  float phi = muons->at(im).phi();

		  float drTopWL = GetDeltaR(eta,phi,Truth_TopWL_Eta,Truth_TopWL_Phi);
		  if( drTopWL < drMinTopWL_Muon && drTopWL < 0.1 && fabs(pt/Truth_TopWL_Pt-1.) < 0.5 && !LisElecTruth )
		    {
		       drMinTopWL_Muon = drTopWL;
		       idxMinTopWL_Muon = im;
		    }		  
	       }	     

	     bool LisElec = (drMinTopWL_Elec < drMinTopWL_Muon) ? 1 : 0;
	     
	     TLorentzVector *TopWLRec = new TLorentzVector();
	     
	     if( idxMinTopWL_Elec >= 0 && LisElec )
	       TopWLRec->SetPtEtaPhiE(electrons->at(idxMinTopWL_Elec).pt(),
				      electrons->at(idxMinTopWL_Elec).eta(),
				      electrons->at(idxMinTopWL_Elec).phi(),
				      electrons->at(idxMinTopWL_Elec).E());
	     
	     if( idxMinTopWL_Muon >= 0 && !LisElec )
	       TopWLRec->SetPtEtaPhiE(muons->at(idxMinTopWL_Muon).pt(),
				      muons->at(idxMinTopWL_Muon).eta(),
				      muons->at(idxMinTopWL_Muon).phi(),
				      muons->at(idxMinTopWL_Muon).E());

	     dElecPx = 10E+10;
	     dElecPy = 10E+10;
	     dElecPz = 10E+10;
	     dElecE = 10E+10;

	     dMuonPx = 10E+10;
	     dMuonPy = 10E+10;
	     dMuonPz = 10E+10;
	     dMuonE = 10E+10;
	     
	     if( LisElec && LisElecTruth )
	       {
		  dElecPx = (idxMinTopWL_Elec >= 0) ? TopWL->Px()-TopWLRec->Px() : 10E+10;
		  dElecPy = (idxMinTopWL_Elec >= 0) ? TopWL->Py()-TopWLRec->Py() : 10E+10;
		  dElecPz = (idxMinTopWL_Elec >= 0) ? TopWL->Pz()-TopWLRec->Pz() : 10E+10;
		  dElecE = (idxMinTopWL_Elec >= 0) ? TopWL->E()-TopWLRec->E() : 10E+10;
	       }

	     if( !LisElec && !LisElecTruth )
	       {		  
		  dMuonPx = (idxMinTopWL_Muon >= 0) ? TopWL->Px()-TopWLRec->Px() : 10E+10;
		  dMuonPy = (idxMinTopWL_Muon >= 0) ? TopWL->Py()-TopWLRec->Py() : 10E+10;
		  dMuonPz = (idxMinTopWL_Muon >= 0) ? TopWL->Pz()-TopWLRec->Pz() : 10E+10;
		  dMuonE = (idxMinTopWL_Muon >= 0) ? TopWL->E()-TopWLRec->E() : 10E+10;
	       }
	     
	     delete TopWLRec;
	     
	     // additional kinematic variables
	     
	     Top1M = Top1->M();
	     Top2M = Top2->M();
	     
	     TopW1M = TopW1->M();
	     TopW2M = TopW2->M();
	     
	     TopWLPt = Truth_TopWL_Pt;
	     TopWLEta = Truth_TopWL_Eta;
	     
	     TopLepRecM = 10E+10;
	     TopHadRecM = 10E+10;
	     
	     if( idxMinTopB1 >= 0 && idxMinTopB2 >= 0 )
	       {		  
		  if( Truth_idxTopWL1 >= 0 )
		    {
		       TopLepRecM = (*TopW1+*TopB1Rec).M();
		       if( idxMinTopWJ1 >= 0 && idxMinTopWJ2 >= 0 ) TopHadRecM = (TopWHadRec+*TopB2Rec).M();
		    }
		  else
		    {
		       TopLepRecM = (*TopW2+*TopB2Rec).M();
		       if( idxMinTopWJ1 >= 0 && idxMinTopWJ2 >= 0 ) TopHadRecM = (TopWHadRec+*TopB1Rec).M();
		    }
	       }	     
		  
	     delete Top1;
	     delete Top2;

	     delete TopW1;
	     delete TopW2;
	     
	     delete TopB1;
	     delete TopB2;

	     delete TopWNu;

	     delete TopB1Rec;
	     delete TopB2Rec;
	     
	     trGEN->Fill();
	  }	     
     }   
   
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
