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

   float TopWL1Pt, TopWL2Pt;
   float TopWL1Eta, TopWL2Eta;   
   float TopWL1TopWL2Dr;
   float Top1Top2Dr, Top1Top2DPhi, Top1Top2DEta;
   float Top1M, Top2M;
   float TopW1M, TopW2M;
   float TopW1TopB1Dr, TopW2TopB2Dr;
   float dMetPx, dMetPy;
   float dBJet1Px, dBJet1Py, dBJet1Pz, dBJet1E;
   float dBJet2Px, dBJet2Py, dBJet2Pz, dBJet2E;
   float dNonBJet1Px, dNonBJet1Py, dNonBJet1Pz, dNonBJet1E;
   float dNonBJet2Px, dNonBJet2Py, dNonBJet2Pz, dNonBJet2E;
   float dElec1Px, dElec1Py, dElec1Pz, dElec1E;
   float dElec2Px, dElec2Py, dElec2Pz, dElec2E;
   float dMuon1Px, dMuon1Py, dMuon1Pz, dMuon1E;
   float dMuon2Px, dMuon2Py, dMuon2Pz, dMuon2E;
   
   TTree *trGEN = new TTree("trGEN","trGEN");
   
   trGEN->Branch("TopWL1Pt",&TopWL1Pt,"TopWL1Pt/F");
   trGEN->Branch("TopWL1Eta",&TopWL1Eta,"TopWL1Eta/F");
   trGEN->Branch("TopWL2Pt",&TopWL2Pt,"TopWL2Pt/F");
   trGEN->Branch("TopWL2Eta",&TopWL2Eta,"TopWL2Eta/F");
   trGEN->Branch("TopWL1TopWL2Dr",&TopWL1TopWL2Dr,"TopWL1TopWL2Dr/F");
   trGEN->Branch("Top1Top2Dr",&Top1Top2Dr,"Top1Top2Dr/F");
   trGEN->Branch("Top1Top2DEta",&Top1Top2DEta,"Top1Top2DEta/F");
   trGEN->Branch("Top1Top2DPhi",&Top1Top2DPhi,"Top1Top2DPhi/F");
   trGEN->Branch("Top1M",&Top1M,"Top1M/F");
   trGEN->Branch("Top2M",&Top2M,"Top2M/F");
   trGEN->Branch("TopW1M",&TopW1M,"TopW1M/F");
   trGEN->Branch("TopW2M",&TopW2M,"TopW2M/F");
   trGEN->Branch("TopW1TopB1Dr",&TopW1TopB1Dr,"TopW1TopB1Dr/F");
   trGEN->Branch("TopW2TopB2Dr",&TopW2TopB2Dr,"TopW2TopB2Dr/F");
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
   trGEN->Branch("dElec1Px",&dElec1Px,"dElec1Px/F");
   trGEN->Branch("dElec1Py",&dElec1Py,"dElec1Py/F");
   trGEN->Branch("dElec1Pz",&dElec1Pz,"dElec1Pz/F");
   trGEN->Branch("dElec1E",&dElec1E,"dElec1E/F");
   trGEN->Branch("dElec2Px",&dElec2Px,"dElec2Px/F");
   trGEN->Branch("dElec2Py",&dElec2Py,"dElec2Py/F");
   trGEN->Branch("dElec2Pz",&dElec2Pz,"dElec2Pz/F");
   trGEN->Branch("dElec2E",&dElec2E,"dElec2E/F");
   trGEN->Branch("dMuon1Px",&dMuon1Px,"dMuon1Px/F");
   trGEN->Branch("dMuon1Py",&dMuon1Py,"dMuon1Py/F");
   trGEN->Branch("dMuon1Pz",&dMuon1Pz,"dMuon1Pz/F");
   trGEN->Branch("dMuon1E",&dMuon1E,"dMuon1E/F");
   trGEN->Branch("dMuon2Px",&dMuon2Px,"dMuon2Px/F");
   trGEN->Branch("dMuon2Py",&dMuon2Py,"dMuon2Py/F");
   trGEN->Branch("dMuon2Pz",&dMuon2Pz,"dMuon2Pz/F");
   trGEN->Branch("dMuon2E",&dMuon2E,"dMuon2E/F");
   
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
		     Truth_idxTopWL1 >= 0 && Truth_idxTopWL2 >= 0 &&
		     Truth_idxTopWNu1 >= 0 && Truth_idxTopWNu2 >= 0);
	
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

//	     std::cout << Truth_Top1_Id << " " << Truth_Top2_Id << std::endl;
	     
	     float Truth_TopWL1_Pt = (Truth_idxTopWL1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopWL1] : -666.;
	     float Truth_TopWL1_Eta = (Truth_idxTopWL1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopWL1] : -666.;
	     float Truth_TopWL1_Phi = (Truth_idxTopWL1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopWL1] : -666.;
	     float Truth_TopWL1_E = (Truth_idxTopWL1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopWL1] : -666.;
	     int Truth_TopWL1_Id = (Truth_idxTopWL1 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopWL1] : -666.;

	     float Truth_TopWL2_Pt = (Truth_idxTopWL2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopWL2] : -666.;
	     float Truth_TopWL2_Eta = (Truth_idxTopWL2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopWL2] : -666.;
	     float Truth_TopWL2_Phi = (Truth_idxTopWL2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopWL2] : -666.;
	     float Truth_TopWL2_E = (Truth_idxTopWL2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopWL2] : -666.;
	     int Truth_TopWL2_Id = (Truth_idxTopWL2 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopWL2] : -666.;

	     float Truth_TopWNu1_Pt = (Truth_idxTopWNu1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopWNu1] : -666.;
	     float Truth_TopWNu1_Eta = (Truth_idxTopWNu1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopWNu1] : -666.;
	     float Truth_TopWNu1_Phi = (Truth_idxTopWNu1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopWNu1] : -666.;
	     float Truth_TopWNu1_E = (Truth_idxTopWNu1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopWNu1] : -666.;
	     int Truth_TopWNu1_Id = (Truth_idxTopWNu1 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopWNu1] : -666.;

	     float Truth_TopWNu2_Pt = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_Eta = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_Phi = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_E = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopWNu2] : -666.;
	     int Truth_TopWNu2_Id = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_id()[Truth_idxTopWNu2] : -666.;
	     
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
	     
	     TLorentzVector *Top1 = new TLorentzVector();
	     Top1->SetPtEtaPhiE(Truth_Top1_Pt,Truth_Top1_Eta,Truth_Top1_Phi,Truth_Top1_E);

	     TLorentzVector *Top2 = new TLorentzVector();
	     Top2->SetPtEtaPhiE(Truth_Top2_Pt,Truth_Top2_Eta,Truth_Top2_Phi,Truth_Top2_E);

	     TLorentzVector *TopW1 = new TLorentzVector();
	     TopW1->SetPtEtaPhiE(Truth_TopW1_Pt,Truth_TopW1_Eta,Truth_TopW1_Phi,Truth_TopW1_E);

	     TLorentzVector *TopW2 = new TLorentzVector();
	     TopW2->SetPtEtaPhiE(Truth_TopW2_Pt,Truth_TopW2_Eta,Truth_TopW2_Phi,Truth_TopW2_E);

	     TLorentzVector *TopWNu1 = new TLorentzVector();
	     TopWNu1->SetPtEtaPhiE(Truth_TopWNu1_Pt,Truth_TopWNu1_Eta,Truth_TopWNu1_Phi,Truth_TopWNu1_E);

	     TLorentzVector *TopWNu2 = new TLorentzVector();
	     TopWNu2->SetPtEtaPhiE(Truth_TopWNu2_Pt,Truth_TopWNu2_Eta,Truth_TopWNu2_Phi,Truth_TopWNu2_E);

	     TLorentzVector *TopWL1 = new TLorentzVector();
	     TopWL1->SetPtEtaPhiE(Truth_TopWL1_Pt,Truth_TopWL1_Eta,Truth_TopWL1_Phi,Truth_TopWL1_E);

	     TLorentzVector *TopWL2 = new TLorentzVector();
	     TopWL2->SetPtEtaPhiE(Truth_TopWL2_Pt,Truth_TopWL2_Eta,Truth_TopWL2_Phi,Truth_TopWL2_E);

	     TLorentzVector *TopB1 = new TLorentzVector();
	     TopB1->SetPtEtaPhiE(Truth_TopB1_Pt,Truth_TopB1_Eta,Truth_TopB1_Phi,Truth_TopB1_E);

	     TLorentzVector *TopB2 = new TLorentzVector();
	     TopB2->SetPtEtaPhiE(Truth_TopB2_Pt,Truth_TopB2_Eta,Truth_TopB2_Phi,Truth_TopB2_E);
	     
	     TLorentzVector MET = *TopWNu1+*TopWNu2;
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

	     delete TopB1Rec;
	     delete TopB2Rec;

	     // match leptons

	     bool L1isElecTruth = (abs(Truth_TopWL1_Id) == 11);
	     bool L2isElecTruth = (abs(Truth_TopWL2_Id) == 11);
	     
	     float drMinTopWL1_Elec = 10E+10;
	     int idxMinTopWL1_Elec = -1;
	     float drMinTopWL1_Muon = 10E+10;
	     int idxMinTopWL1_Muon = -1;
	     float drMinTopWL2_Elec = 10E+10;
	     int idxMinTopWL2_Elec = -1;
	     float drMinTopWL2_Muon = 10E+10;
	     int idxMinTopWL2_Muon = -1;
	     
	     int NElec = electrons->size();
	     int NMuon = muons->size();
	     
	     for(int ie=0;ie<NElec;ie++)
	       {
		  float pt = electrons->at(ie).pt();
		  float eta = electrons->at(ie).eta();
		  float phi = electrons->at(ie).phi();
		  
		  float drTopWL1 = GetDeltaR(eta,phi,Truth_TopWL1_Eta,Truth_TopWL1_Phi);
		  if( drTopWL1 < drMinTopWL1_Elec && drTopWL1 < 0.1 && fabs(pt/Truth_TopWL1_Pt-1.) < 0.5 && L1isElecTruth )
		    {
		       drMinTopWL1_Elec = drTopWL1;
		       idxMinTopWL1_Elec = ie;
		    }		  
		  float drTopWL2 = GetDeltaR(eta,phi,Truth_TopWL2_Eta,Truth_TopWL2_Phi);
		  if( drTopWL2 < drMinTopWL2_Elec && drTopWL2 < 0.1 && fabs(pt/Truth_TopWL2_Pt-1.) < 0.5 && L2isElecTruth ) 
		    {
		       drMinTopWL2_Elec = drTopWL2;
		       idxMinTopWL2_Elec = ie;
		    }		  
	       }	     
	     
	     for(int im=0;im<NMuon;im++)
	       {
		  float pt = muons->at(im).pt();
		  float eta = muons->at(im).eta();
		  float phi = muons->at(im).phi();

		  float drTopWL1 = GetDeltaR(eta,phi,Truth_TopWL1_Eta,Truth_TopWL1_Phi);
		  if( drTopWL1 < drMinTopWL1_Muon && drTopWL1 < 0.1 && fabs(pt/Truth_TopWL1_Pt-1.) < 0.5 && !L1isElecTruth )
		    {
		       drMinTopWL1_Muon = drTopWL1;
		       idxMinTopWL1_Muon = im;
		    }		  
		  float drTopWL2 = GetDeltaR(eta,phi,Truth_TopWL2_Eta,Truth_TopWL2_Phi);
		  if( drTopWL2 < drMinTopWL2_Muon && drTopWL2 < 0.1 && fabs(pt/Truth_TopWL2_Pt-1.) < 0.5 && !L2isElecTruth )
		    {
		       drMinTopWL2_Muon = drTopWL2;
		       idxMinTopWL2_Muon = im;
		    }		  
	       }	     

	     bool L1isElec = (drMinTopWL1_Elec < drMinTopWL1_Muon) ? 1 : 0;
	     bool L2isElec = (drMinTopWL2_Elec < drMinTopWL2_Muon) ? 1 : 0;
	     
	     TLorentzVector *TopWL1Rec = new TLorentzVector();
	     TLorentzVector *TopWL2Rec = new TLorentzVector();
	     
	     if( idxMinTopWL1_Elec >= 0 && L1isElec )
	       TopWL1Rec->SetPtEtaPhiE(electrons->at(idxMinTopWL1_Elec).pt(),
				       electrons->at(idxMinTopWL1_Elec).eta(),
				       electrons->at(idxMinTopWL1_Elec).phi(),
				       electrons->at(idxMinTopWL1_Elec).E());

	     if( idxMinTopWL1_Muon >= 0 && !L1isElec )
	       TopWL1Rec->SetPtEtaPhiE(muons->at(idxMinTopWL1_Muon).pt(),
				       muons->at(idxMinTopWL1_Muon).eta(),
				       muons->at(idxMinTopWL1_Muon).phi(),
				       muons->at(idxMinTopWL1_Muon).E());
	     
	     if( idxMinTopWL2_Elec >= 0 && L2isElec )
	       TopWL2Rec->SetPtEtaPhiE(electrons->at(idxMinTopWL2_Elec).pt(),
				       electrons->at(idxMinTopWL2_Elec).eta(),
				       electrons->at(idxMinTopWL2_Elec).phi(),
				       electrons->at(idxMinTopWL2_Elec).E());

	     if( idxMinTopWL2_Muon >= 0 && !L2isElec )
	       TopWL2Rec->SetPtEtaPhiE(muons->at(idxMinTopWL2_Muon).pt(),
				       muons->at(idxMinTopWL2_Muon).eta(),
				       muons->at(idxMinTopWL2_Muon).phi(),
				       muons->at(idxMinTopWL2_Muon).E());

	     dElec1Px = 10E+10;
	     dElec1Py = 10E+10;
	     dElec1Pz = 10E+10;
	     dElec1E = 10E+10;

	     dElec2Px = 10E+10;
	     dElec2Py = 10E+10;
	     dElec2Pz = 10E+10;
	     dElec2E = 10E+10;

	     dMuon1Px = 10E+10;
	     dMuon1Py = 10E+10;
	     dMuon1Pz = 10E+10;
	     dMuon1E = 10E+10;

	     dMuon2Px = 10E+10;
	     dMuon2Py = 10E+10;
	     dMuon2Pz = 10E+10;
	     dMuon2E = 10E+10;
	     
	     if( L1isElec && L1isElecTruth )
	       {
		  dElec1Px = (idxMinTopWL1_Elec >= 0) ? TopWL1->Px()-TopWL1Rec->Px() : 10E+10;
		  dElec1Py = (idxMinTopWL1_Elec >= 0) ? TopWL1->Py()-TopWL1Rec->Py() : 10E+10;
		  dElec1Pz = (idxMinTopWL1_Elec >= 0) ? TopWL1->Pz()-TopWL1Rec->Pz() : 10E+10;
		  dElec1E = (idxMinTopWL1_Elec >= 0) ? TopWL1->E()-TopWL1Rec->E() : 10E+10;
	       }
	     if( L2isElec && L2isElecTruth )
	       {		  
		  dElec2Px = (idxMinTopWL2_Elec >= 0) ? TopWL2->Px()-TopWL2Rec->Px() : 10E+10;
		  dElec2Py = (idxMinTopWL2_Elec >= 0) ? TopWL2->Py()-TopWL2Rec->Py() : 10E+10;
		  dElec2Pz = (idxMinTopWL2_Elec >= 0) ? TopWL2->Pz()-TopWL2Rec->Pz() : 10E+10;
		  dElec2E = (idxMinTopWL2_Elec >= 0) ? TopWL2->E()-TopWL2Rec->E() : 10E+10;
	       }	     

	     if( !L1isElec && !L1isElecTruth )
	       {		  
		  dMuon1Px = (idxMinTopWL1_Muon >= 0) ? TopWL1->Px()-TopWL1Rec->Px() : 10E+10;
		  dMuon1Py = (idxMinTopWL1_Muon >= 0) ? TopWL1->Py()-TopWL1Rec->Py() : 10E+10;
		  dMuon1Pz = (idxMinTopWL1_Muon >= 0) ? TopWL1->Pz()-TopWL1Rec->Pz() : 10E+10;
		  dMuon1E = (idxMinTopWL1_Muon >= 0) ? TopWL1->E()-TopWL1Rec->E() : 10E+10;
	       }
	     if( !L2isElec && !L2isElecTruth )
	       {		  
		  dMuon2Px = (idxMinTopWL2_Muon >= 0) ? TopWL2->Px()-TopWL2Rec->Px() : 10E+10;
		  dMuon2Py = (idxMinTopWL2_Muon >= 0) ? TopWL2->Py()-TopWL2Rec->Py() : 10E+10;
		  dMuon2Pz = (idxMinTopWL2_Muon >= 0) ? TopWL2->Pz()-TopWL2Rec->Pz() : 10E+10;
		  dMuon2E = (idxMinTopWL2_Muon >= 0) ? TopWL2->E()-TopWL2Rec->E() : 10E+10;
	       }	     
	     
	     delete TopWL1Rec;
	     delete TopWL2Rec;
	     
	     // additional kinematic variables
	     
	     Top1Top2Dr = GetDeltaR(Truth_Top1_Eta,Truth_Top1_Phi,
				    Truth_Top2_Eta,Truth_Top2_Phi);
	     
	     Top1M = Top1->M();
	     Top2M = Top2->M();

	     TopW1M = TopW1->M();
	     TopW2M = TopW2->M();
	     
	     Top1Top2DEta = GetDeltaEta(Truth_Top1_Eta,Truth_Top2_Eta);
	     
	     Top1Top2DPhi = GetDeltaPhi(Truth_Top1_Phi,Truth_Top2_Phi);
	     
	     TopWL1TopWL2Dr = GetDeltaR(Truth_TopWL1_Eta,Truth_TopWL1_Phi,
					Truth_TopWL2_Eta,Truth_TopWL2_Phi);
	     
	     TopWL1Pt = Truth_TopWL1_Pt;
	     TopWL1Eta = Truth_TopWL1_Eta;	     
	     
	     TopWL2Pt = Truth_TopWL2_Pt;
	     TopWL2Eta = Truth_TopWL2_Eta;
	     
	     TopW1TopB1Dr = GetDeltaR(Truth_TopW1_Eta,Truth_TopW1_Phi,
				      Truth_TopB1_Eta,Truth_TopB1_Phi);

	     TopW2TopB2Dr = GetDeltaR(Truth_TopW2_Eta,Truth_TopW2_Phi,
				      Truth_TopB2_Eta,Truth_TopB2_Phi);
	     
	     delete Top1;
	     delete Top2;

	     delete TopW1;
	     delete TopW2;
	     
	     delete TopB1;
	     delete TopB2;

	     delete TopWNu1;
	     delete TopWNu2;
	     
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
