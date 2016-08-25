#include "kinfit.h"

#include "include/Jet.h"
#include "include/Electron.h"
#include "include/Muon.h"
#include "include/Truth.h"
#include "include/Event.h"

#include "TH2D.h"

float GetDeltaR(float eta1,float phi1,float eta2,float phi2);

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

   int nToys = atoi(argv[3]);
   
   std::cout << "Run " << nToys << " toys" << std::endl;
   
   KINFIT::kfit *kf = new KINFIT::kfit();

   kf->Init(TOPLEP);
   
   std::string pdfFileName = "/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/KinFit/test/GenAnalysis/TopLep/pdf.root";
   kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopWM_Fit");
   kf->SetPDF("TopMass",pdfFileName.c_str(),"TopM_Fit");
   kf->SetPDF("MetPx",pdfFileName.c_str(),"dMetPx_Gaus");
   kf->SetPDF("MetPy",pdfFileName.c_str(),"dMetPy_Gaus");
   kf->SetPDF("BJetPx",pdfFileName.c_str(),"dBJetPx_Fit");
   kf->SetPDF("BJetPy",pdfFileName.c_str(),"dBJetPy_Fit");
   kf->SetPDF("BJetPz",pdfFileName.c_str(),"dBJetPz_Fit");
   kf->SetPDF("BJetE",pdfFileName.c_str(),"dBJetE_Fit");
   kf->SetPDF("ElecPx",pdfFileName.c_str(),"dElecPx_Fit");
   kf->SetPDF("ElecPy",pdfFileName.c_str(),"dElecPy_Fit");
   kf->SetPDF("ElecPz",pdfFileName.c_str(),"dElecPz_Fit");
   kf->SetPDF("ElecE",pdfFileName.c_str(),"dElecE_Fit");
   kf->SetPDF("MuonPx",pdfFileName.c_str(),"dMuonPx_Fit");
   kf->SetPDF("MuonPy",pdfFileName.c_str(),"dMuonPy_Fit");
   kf->SetPDF("MuonPz",pdfFileName.c_str(),"dMuonPz_Fit");
   kf->SetPDF("MuonE",pdfFileName.c_str(),"dMuonE_Fit");

   kf->SetWMassBW(80.4,2.1);
   kf->SetTopMass(172.5);
   kf->SetNToy(nToys);
   kf->SetNWRMS(4);
   kf->SetNMetRMS(4);

   kf->SetNBJetPxRMS(3);
   kf->SetNBJetPyRMS(3);
   kf->SetNBJetPzRMS(3);
   kf->SetNBJetERMS(3);

   kf->SetNElecPxRMS(3);
   kf->SetNElecPyRMS(3);
   kf->SetNElecPzRMS(3);
   kf->SetNElecERMS(3);

   kf->SetNMuonPxRMS(3);
   kf->SetNMuonPyRMS(3);
   kf->SetNMuonPzRMS(3);
   kf->SetNMuonERMS(3);

   TFile *f = TFile::Open(fname.c_str());
   
   TTree *tr = (TTree*)f->Get("Nt");
   
   tr->SetBranchAddress("Jet",&jets);
   tr->SetBranchAddress("Electron",&electrons);
   tr->SetBranchAddress("Muon",&muons);
   tr->SetBranchAddress("Truth",&truth);
   tr->SetBranchAddress("Event",&event);

   TFile *fout = new TFile(foutName.c_str(),"RECREATE");

   float P0Disc;
   float P0Prob;
   
   bool P0MatchBJet;
   bool P0RecoBJet;
   bool P0MatchTopWL;
   bool P0RecoTopWL;
   
   float P0Term0, P0Term1;

   float NuPx;
   float NuPy;
   float NuPz;
   float WMassGen;
   float WMass;
   float MetXMeas, MetYMeas;
   float MetXFit, MetYFit;

   float NuGenPx;
   float NuGenPy;
   float NuGenPz;
   
   float TopP, TopPt, TopEta, TopM, TopRap;
   float WP, WPt, WEta, WPhi, WRap;

   float TopGenP, TopGenPt, TopGenEta, TopGenM, TopGenRap;
   float WGenP, WGenPt, WGenEta, WGenPhi, WGenRap;
   
   double tim;
   
   TTree *trKFIT = new TTree("trKFIT","trKFIT");
   
   trKFIT->Branch("P0Disc",&P0Disc,"P0Disc/F");
   trKFIT->Branch("P0Prob",&P0Prob,"P0Prob/F");
   
   trKFIT->Branch("P0Term0",&P0Term0,"P0Term0/F");
   trKFIT->Branch("P0Term1",&P0Term1,"P0Term1/F");
   
   trKFIT->Branch("P0MatchBJet",&P0MatchBJet,"P0MatchBJet/O");
   trKFIT->Branch("P0RecoBJet",&P0RecoBJet,"P0RecoBJet/O");
   
   trKFIT->Branch("P0MatchTopWL",&P0MatchTopWL,"P0MatchTopWL/O");
   trKFIT->Branch("P0RecoTopWL",&P0RecoTopWL,"P0RecoTopWL/O");
   
   trKFIT->Branch("NuPx",&NuPx,"NuPx/F");
   trKFIT->Branch("NuPy",&NuPy,"NuPy/F");
   trKFIT->Branch("NuPz",&NuPz,"NuPz/F");
   trKFIT->Branch("WMassGen",&WMassGen,"WMassGen/F");
   trKFIT->Branch("WMass",&WMass,"WMass/F");
   trKFIT->Branch("MetXMeas",&MetXMeas,"MetXMeas/F");
   trKFIT->Branch("MetYMeas",&MetYMeas,"MetYMeas/F");
   trKFIT->Branch("MetXFit",&MetXFit,"MetXFit/F");
   trKFIT->Branch("MetYFit",&MetYFit,"MetYFit/F");

   trKFIT->Branch("NuGenPx",&NuGenPx,"NuGenPx/F");
   trKFIT->Branch("NuGenPy",&NuGenPy,"NuGenPy/F");
   trKFIT->Branch("NuGenPz",&NuGenPz,"NuGenPz/F");
   
   trKFIT->Branch("TopPt",&TopPt,"TopPt/F");
   trKFIT->Branch("TopP",&TopP,"TopP/F");
   trKFIT->Branch("TopEta",&TopEta,"TopEta/F");
   trKFIT->Branch("TopRap",&TopRap,"TopRap/F");
   trKFIT->Branch("TopM",&TopM,"TopM/F");

   trKFIT->Branch("TopGenPt",&TopGenPt,"TopGenPt/F");
   trKFIT->Branch("TopGenP",&TopGenP,"TopGenP/F");
   trKFIT->Branch("TopGenEta",&TopGenEta,"TopGenEta/F");
   trKFIT->Branch("TopGenRap",&TopGenRap,"TopGenRap/F");
   trKFIT->Branch("TopGenM",&TopGenM,"TopGenM/F");

   trKFIT->Branch("WPt",&WPt,"WPt/F");
   trKFIT->Branch("WP",&WP,"WP/F");
   trKFIT->Branch("WEta",&WEta,"WEta/F");
   trKFIT->Branch("WPhi",&WPhi,"WPhi/F");
   trKFIT->Branch("WRap",&WRap,"WRap/F");

   trKFIT->Branch("WGenPt",&WGenPt,"WGenPt/F");
   trKFIT->Branch("WGenP",&WGenP,"WGenP/F");
   trKFIT->Branch("WGenEta",&WGenEta,"WGenEta/F");
   trKFIT->Branch("WGenPhi",&WGenPhi,"WGenPhi/F");
   trKFIT->Branch("WGenRap",&WGenRap,"WGenRap/F");
   
   trKFIT->Branch("tim",&tim,"tim/D");
   
   int nev = tr->GetEntries();
   std::cout << "Total number of events = " << nev << std::endl;

   int nMax = -1;
//   int nMax = 100;   
   
   for(int i=0;i<nev;i++)
     {
	if( nMax >=0 && i > nMax ) break;

	std::vector<float> BJetPt;
	std::vector<float> BJetEta;
	std::vector<float> BJetPhi;
	std::vector<float> BJetE;

	std::vector<float> NonBJetPt;
	std::vector<float> NonBJetEta;
	std::vector<float> NonBJetPhi;
	std::vector<float> NonBJetE;

	std::vector<float> ElectronPt;
	std::vector<float> ElectronEta;
	std::vector<float> ElectronPhi;
	std::vector<float> ElectronE;

	std::vector<float> MuonPt;
	std::vector<float> MuonEta;
	std::vector<float> MuonPhi;
	std::vector<float> MuonE;
	
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
	
	int nElectrons = electrons->size();
	int nMuons = muons->size();
	int nJets = jets->size();

	bool isL1 = (Truth_idxTopWL1 >= 0);
	
	bool pass = (Truth_idxTopB1 >= 0 && Truth_idxTopB2 >= 0 &&
		     ((Truth_idxTopWL1 >= 0 && !(Truth_idxTopWL2 >= 0)) || (!(Truth_idxTopWL1 >= 0) && Truth_idxTopWL2 >= 0)) &&
		     ((Truth_idxTopWJ11 >= 0 && Truth_idxTopWJ21 >= 0) ||
			 (Truth_idxTopWJ12 >= 0 && Truth_idxTopWJ22 >= 0)));
	
	bool passReco = 0;
	
	if( pass )
	  {
	     tim = 10E+10;
	     clock_t tStart = clock();
	     
	     // reco-level selection

	     // BJets et NonBJets
	     for(int ij=0;ij<nJets;ij++)
	       {		  	     
		  float jpt = jets->at(ij).pt();
		  float jeta = jets->at(ij).eta();
		  float jphi = jets->at(ij).phi();
		  float jE = jets->at(ij).E();

		  if( jpt < 20. ) continue;
		  
		  float btag = jets->at(ij).CSVv2();
		  bool isBTag = 0;
		  if( btag > 0.800 ) isBTag = 1;
//		  if( btag > 0.460 ) isBTag = 1;

		  if( isBTag )
		    {		       
		       BJetPt.push_back(jpt);
		       BJetEta.push_back(jeta);
		       BJetPhi.push_back(jphi);
		       BJetE.push_back(jE);
		    }	
		  else
		    {
		       NonBJetPt.push_back(jpt);
		       NonBJetEta.push_back(jeta);
		       NonBJetPhi.push_back(jphi);
		       NonBJetE.push_back(jE);
		    }		  
	       }	     

	     // Electrons
	     for(int ie=0;ie<nElectrons;ie++)
	       {		  	     
		  float ept = electrons->at(ie).pt();
		  float eeta = electrons->at(ie).eta();
		  float ephi = electrons->at(ie).phi();
		  float eE = electrons->at(ie).E();

		  if( ept < 20. ) continue;
		  
		  ElectronPt.push_back(ept);
		  ElectronEta.push_back(eeta);
		  ElectronPhi.push_back(ephi);
		  ElectronE.push_back(eE);
	       }

	     // Muons
	     for(int im=0;im<nMuons;im++)
	       {		  	     
		  float mpt = muons->at(im).pt();
		  float meta = muons->at(im).eta();
		  float mphi = muons->at(im).phi();
		  float mE = muons->at(im).E();

		  if( mpt < 20. ) continue;
		  
		  MuonPt.push_back(mpt);
		  MuonEta.push_back(meta);
		  MuonPhi.push_back(mphi);
		  MuonE.push_back(mE);
	       }	     

	     // MET
	     float metpt = event->at(0).metpt();
	     float metphi = event->at(0).metphi();
	     float metpx = metpt*cos(metphi);
	     float metpy = metpt*sin(metphi);
	     float metE = sqrt(metpx*metpx+metpy*metpy);
	     
	     passReco = (BJetPt.size() >= 2 &&
			 (ElectronPt.size()+MuonPt.size()) == 1);

	     if( !passReco ) continue;
	     
	     std::cout << "Event #" << i << std::endl;
	     
	     kf->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
	     kf->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
	     
	     kf->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
	     kf->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
	     
	     kf->SetMet(metpx,metpy);
	     
	     kf->Run();
	     
	     P0Disc = kf->GetDisc();
	     int nPerm = kf->GetNPerm();
	     std::cout << "Permutations=" << nPerm << std::endl;
	     std::cout << "Nbjet=" << BJetPt.size() << " Njet=" << NonBJetPt.size() <<
	       " Nlep=" << ElectronPt.size()+MuonPt.size() << std::endl;
	     std::cout << "Chi2=" << P0Disc << std::endl;
	     
	     P0Prob = 10E+10;

	     P0Term0 = 10E+10;
	     P0Term1 = 10E+10;
	     
	     P0MatchBJet = 0;
	     P0RecoBJet = 0;

	     P0MatchTopWL = 0;
	     P0RecoTopWL = 0;
	     
	     int nDF = 0;
	     for(int ip=0;ip<nPerm;ip++)
	       {
		  if(ip == 0)
		    {
		       P0Term0 = kf->GetDiscTerm(ip,0);
		       P0Term1 = kf->GetDiscTerm(ip,1);
		       P0Prob = TMath::Prob(P0Disc,nDF);
		    }		  
	       }		  

	     NuPx = kf->GetNuPx(0,0);
	     NuPy = kf->GetNuPy(0,0);
	     NuPz = kf->GetNuPz(0,0);
	     MetXFit = kf->GetMetX(0);
	     MetYFit = kf->GetMetY(0);
	     MetXMeas = metpx;
	     MetYMeas = metpy;
	     WMass = kf->GetWMass(0,0);
	     
	     TopPt = kf->GetTopPt(0,0);
	     TopP = kf->GetTopP(0,0);
	     TopEta = kf->GetTopEta(0,0);
	     TopM = kf->GetTopMass(0,0);
	     TopRap = kf->GetTopRap(0,0);
	     
	     WPt = kf->GetWPt(0,0);
	     WP = kf->GetWP(0,0);
	     WEta = kf->GetWEta(0,0);
	     WPhi = kf->GetWPhi(0,0);
	     WRap = kf->GetWRap(0,0);	     
	     
	     float Truth_Top_Pt = (isL1) ? truth->at(0).mc_truth_pt()[Truth_idxTop1] : truth->at(0).mc_truth_pt()[Truth_idxTop2];
	     float Truth_Top_Eta = (isL1) ? truth->at(0).mc_truth_eta()[Truth_idxTop1] : truth->at(0).mc_truth_eta()[Truth_idxTop2];
	     float Truth_Top_Phi = (isL1) ? truth->at(0).mc_truth_phi()[Truth_idxTop1] : truth->at(0).mc_truth_phi()[Truth_idxTop2];
	     float Truth_Top_E = (isL1) ? truth->at(0).mc_truth_E()[Truth_idxTop1] : truth->at(0).mc_truth_E()[Truth_idxTop2];	     
	     
	     float Truth_TopB_Pt = (isL1) ? truth->at(0).mc_truth_pt()[Truth_idxTopB1] : truth->at(0).mc_truth_pt()[Truth_idxTopB2];
	     float Truth_TopB_Eta = (isL1) ? truth->at(0).mc_truth_eta()[Truth_idxTopB1] : truth->at(0).mc_truth_eta()[Truth_idxTopB2];
	     float Truth_TopB_Phi = (isL1) ? truth->at(0).mc_truth_phi()[Truth_idxTopB1] : truth->at(0).mc_truth_phi()[Truth_idxTopB2];
	     float Truth_TopB_E = (isL1) ? truth->at(0).mc_truth_E()[Truth_idxTopB1] : truth->at(0).mc_truth_E()[Truth_idxTopB2];	     
	     
	     float Truth_TopW_Pt = (isL1) ? truth->at(0).mc_truth_pt()[Truth_idxTopW1] : truth->at(0).mc_truth_pt()[Truth_idxTopW2];
	     float Truth_TopW_Eta = (isL1) ? truth->at(0).mc_truth_eta()[Truth_idxTopW1] : truth->at(0).mc_truth_eta()[Truth_idxTopW2];
	     float Truth_TopW_Phi = (isL1) ? truth->at(0).mc_truth_phi()[Truth_idxTopW1] : truth->at(0).mc_truth_phi()[Truth_idxTopW2];
	     float Truth_TopW_E = (isL1) ? truth->at(0).mc_truth_E()[Truth_idxTopW1] : truth->at(0).mc_truth_E()[Truth_idxTopW2];	     
	     
	     float Truth_TopWL_Pt = (isL1) ? truth->at(0).mc_truth_pt()[Truth_idxTopWL1] : truth->at(0).mc_truth_pt()[Truth_idxTopWL2];
	     float Truth_TopWL_Eta = (isL1) ? truth->at(0).mc_truth_eta()[Truth_idxTopWL1] : truth->at(0).mc_truth_eta()[Truth_idxTopWL2];
	     float Truth_TopWL_Phi = (isL1) ? truth->at(0).mc_truth_phi()[Truth_idxTopWL1] : truth->at(0).mc_truth_phi()[Truth_idxTopWL2];
	     float Truth_TopWL_E = (isL1) ? truth->at(0).mc_truth_E()[Truth_idxTopWL1] : truth->at(0).mc_truth_E()[Truth_idxTopWL2];
	     int Truth_TopWL_Id = (isL1) ? truth->at(0).mc_truth_id()[Truth_idxTopWL1] : truth->at(0).mc_truth_id()[Truth_idxTopWL2];

	     float Truth_TopWNu_Pt = (isL1) ? truth->at(0).mc_truth_pt()[Truth_idxTopWNu1] : truth->at(0).mc_truth_pt()[Truth_idxTopWNu2];
	     float Truth_TopWNu_Eta = (isL1) ? truth->at(0).mc_truth_eta()[Truth_idxTopWNu1] : truth->at(0).mc_truth_eta()[Truth_idxTopWNu2];
	     float Truth_TopWNu_Phi = (isL1) ? truth->at(0).mc_truth_phi()[Truth_idxTopWNu1] : truth->at(0).mc_truth_phi()[Truth_idxTopWNu2];
	     float Truth_TopWNu_E = (isL1) ? truth->at(0).mc_truth_E()[Truth_idxTopWNu1] : truth->at(0).mc_truth_E()[Truth_idxTopWNu2];
	     int Truth_TopWNu_Id = (isL1) ? truth->at(0).mc_truth_id()[Truth_idxTopWNu1] : truth->at(0).mc_truth_id()[Truth_idxTopWNu2];
	     
	     TLorentzVector *TopWNu = new TLorentzVector();

	     TopWNu->SetPtEtaPhiE(Truth_TopWNu_Pt,Truth_TopWNu_Eta,Truth_TopWNu_Phi,Truth_TopWNu_E);
	     
	     NuGenPx = TopWNu->Px();
	     NuGenPy = TopWNu->Py();
	     NuGenPz = TopWNu->Pz();

	     delete TopWNu;

	     TLorentzVector *Top = new TLorentzVector();

	     Top->SetPtEtaPhiE(Truth_Top_Pt,Truth_Top_Eta,Truth_Top_Phi,Truth_Top_E);
	     
	     TopGenPt = Top->Pt();
	     TopGenP = Top->P();
	     TopGenEta = Top->PseudoRapidity();
	     TopGenRap = Top->Rapidity();
	     TopGenM = Top->M();

	     delete Top;

	     TLorentzVector *TopW = new TLorentzVector();

	     TopW->SetPtEtaPhiE(Truth_TopW_Pt,Truth_TopW_Eta,Truth_TopW_Phi,Truth_TopW_E);
	     
	     WMassGen = TopW->M();
	     
	     WGenP = TopW->P();
	     WGenPt = TopW->Pt();
	     WGenEta = TopW->PseudoRapidity();
	     WGenRap = TopW->Rapidity();
	     WGenPhi = TopW->Phi();
	     
	     delete TopW;
	     
	     for(int ip=0;ip<nPerm;ip++)
	       {
		  if( kf->GetDisc(ip) > 10E+9 ) continue;
		  
		  // b jets
		       
		  bool matchBJet = 0;
		  
		  int bjetIdx = kf->GetIndex(BJET_TOPLEP,ip);
		  
		  float bjetEta = BJetEta[bjetIdx];
		  float bjetPhi = BJetPhi[bjetIdx];
		  
		  float dr_TruthTopB_BJet = GetDeltaR(Truth_TopB_Eta,Truth_TopB_Phi,bjetEta,bjetPhi);
		  
		  if( dr_TruthTopB_BJet < 0.4 ) matchBJet = 1;
		  
		  bool recoBJet = 0;
		  
		  float dr_TruthTopB_min = 666.;
		  
		  for(int ij=0;ij<BJetPt.size();ij++)
		    {
		       float jeta = BJetEta[ij];
		       float jphi = BJetPhi[ij];
		       
		       float dr_TruthTopB = GetDeltaR(Truth_TopB_Eta,Truth_TopB_Phi,jeta,jphi);
		       if( dr_TruthTopB < dr_TruthTopB_min ) dr_TruthTopB_min = dr_TruthTopB;
		    }			   
		  
		  if( dr_TruthTopB_min < 0.4 ) recoBJet = 1;
		  
		  // Top lepton
		       
		  bool matchTopWL = 0;
		  
		  bool topWLisElectron = 1;
		  int topWLIdx = kf->GetIndex(ELECTRON_TOPLEP,ip);
		  if( topWLIdx < 0 )
		    {
		       topWLIdx = kf->GetIndex(MUON_TOPLEP,ip);
		       topWLisElectron = 0;
		    }		       
		  
		  float topWLEta = (topWLisElectron) ? ElectronEta[topWLIdx] : MuonEta[topWLIdx];
		  float topWLPhi = (topWLisElectron) ? ElectronPhi[topWLIdx] : MuonPhi[topWLIdx];
		  
		  float dr_TruthTopWL_TopWL = GetDeltaR(Truth_TopWL_Eta,Truth_TopWL_Phi,topWLEta,topWLPhi);
		  
		  if( dr_TruthTopWL_TopWL < 0.1 && dr_TruthTopB_BJet < 0.4 ) matchTopWL = 1;
		  
		  bool recoTopWL = 0;
		  
		  float dr_TruthTopWL_min = 666.;
		  
		  if( topWLisElectron )
		    {			    
		       for(int ie=0;ie<ElectronPt.size();ie++)
			 {
			    float eeta = ElectronEta[ie];
			    float ephi = ElectronPhi[ie];
			    
			    float dr_TruthTopWL = GetDeltaR(Truth_TopWL_Eta,Truth_TopWL_Phi,eeta,ephi);
			    if( dr_TruthTopWL < dr_TruthTopWL_min ) dr_TruthTopWL_min = dr_TruthTopWL;
			 }
		    }		       
		  else
		    {
		       for(int im=0;im<MuonPt.size();im++)
			 {
			    float meta = MuonEta[im];
			    float mphi = MuonPhi[im];
			    
			    float dr_TruthTopWL = GetDeltaR(Truth_TopWL_Eta,Truth_TopWL_Phi,meta,mphi);
			    if( dr_TruthTopWL < dr_TruthTopWL_min ) dr_TruthTopWL_min = dr_TruthTopWL;
			 }
		    }		       
		  
		  if( dr_TruthTopWL_min < 0.1 ) recoTopWL = 1;
		  
		  if( ip == 0 )
		    {
		       P0MatchBJet = matchBJet;		       
		       P0RecoBJet = recoBJet;
		       
		       P0MatchTopWL = matchTopWL;		       
		       P0RecoTopWL = recoTopWL;
		    }		
	       }	     
		  
	     tim = (double)(clock() - tStart)/CLOCKS_PER_SEC;
	     
	     trKFIT->Fill();
	  }	
     }   
   
   fout->Write();
   fout->Close();
   
   f->Close();
}

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{   
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
