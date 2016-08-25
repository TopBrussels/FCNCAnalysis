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

   kf->Init(TOPTOPLEPLEP);
   
   std::string pdfFileName = "/home-pbs/kskovpen/ttH/CMSSW_7_6_3_patch1/src/KinFit/test/GenAnalysis/TopTopLepLep/pdf.root";
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
   
   bool P0MatchBJet1, P0MatchBJet2;   
   bool P0RecoBJet1, P0RecoBJet2;
   bool P0MatchTopWL1, P0MatchTopWL2;
   bool P0RecoTopWL1, P0RecoTopWL2;
   
   float P0Term0, P0Term1;

   float Nu1Px, Nu2Px;
   float Nu1Py, Nu2Py;
   float Nu1Pz, Nu2Pz;
   float WMassGen1, WMassGen2;
   float WMass1, WMass2;
   float MetXMeas, MetYMeas;
   float MetXFit, MetYFit;

   float NuGen1Px, NuGen2Px;
   float NuGen1Py, NuGen2Py;
   float NuGen1Pz, NuGen2Pz;
   
   float DrTopTop, MTopTop, PTopTop, PtTopTop, EtaTopTop, PhiTopTop, RapTopTop;
   float Top1P, Top1Pt, Top1Eta, Top1M, Top1Rap;
   float Top2P, Top2Pt, Top2Eta, Top2M, Top2Rap;
   float W1P, W1Pt, W1Eta, W1Phi, W1Rap;
   float W2P, W2Pt, W2Eta, W2Phi, W2Rap;

   float DrTopTopGen, MTopTopGen, PTopTopGen, PtTopTopGen, EtaTopTopGen, PhiTopTopGen, RapTopTopGen;
   float TopGen1P, TopGen1Pt, TopGen1Eta, TopGen1M, TopGen1Rap;
   float TopGen2P, TopGen2Pt, TopGen2Eta, TopGen2M, TopGen2Rap;
   float WGen1P, WGen1Pt, WGen1Eta, WGen1Phi, WGen1Rap;
   float WGen2P, WGen2Pt, WGen2Eta, WGen2Phi, WGen2Rap;
   
   double tim;
   
   TTree *trKFIT = new TTree("trKFIT","trKFIT");
   
   trKFIT->Branch("P0Disc",&P0Disc,"P0Disc/F");
   trKFIT->Branch("P0Prob",&P0Prob,"P0Prob/F");
   
   trKFIT->Branch("P0Term0",&P0Term0,"P0Term0/F");
   trKFIT->Branch("P0Term1",&P0Term1,"P0Term1/F");
   
   trKFIT->Branch("P0MatchBJet1",&P0MatchBJet1,"P0MatchBJet1/O");
   trKFIT->Branch("P0MatchBJet2",&P0MatchBJet2,"P0MatchBJet2/O");
   trKFIT->Branch("P0RecoBJet1",&P0RecoBJet1,"P0RecoBJet1/O");
   trKFIT->Branch("P0RecoBJet2",&P0RecoBJet2,"P0RecoBJet2/O");

   trKFIT->Branch("P0MatchTopWL1",&P0MatchTopWL1,"P0MatchTopWL1/O");
   trKFIT->Branch("P0MatchTopWL2",&P0MatchTopWL2,"P0MatchTopWL2/O");
   trKFIT->Branch("P0RecoTopWL1",&P0RecoTopWL1,"P0RecoTopWL1/O");
   trKFIT->Branch("P0RecoTopWL2",&P0RecoTopWL2,"P0RecoTopWL2/O");
   
   trKFIT->Branch("Nu1Px",&Nu1Px,"Nu1Px/F");
   trKFIT->Branch("Nu2Px",&Nu2Px,"Nu2Px/F");
   trKFIT->Branch("Nu1Py",&Nu1Py,"Nu1Py/F");
   trKFIT->Branch("Nu2Py",&Nu2Py,"Nu2Py/F");
   trKFIT->Branch("Nu1Pz",&Nu1Pz,"Nu1Pz/F");
   trKFIT->Branch("Nu2Pz",&Nu2Pz,"Nu2Pz/F");
   trKFIT->Branch("WMassGen1",&WMassGen1,"WMassGen1/F");
   trKFIT->Branch("WMassGen2",&WMassGen2,"WMassGen2/F");
   trKFIT->Branch("WMass1",&WMass1,"WMass1/F");
   trKFIT->Branch("WMass2",&WMass2,"WMass2/F");
   trKFIT->Branch("MetXMeas",&MetXMeas,"MetXMeas/F");
   trKFIT->Branch("MetYMeas",&MetYMeas,"MetYMeas/F");
   trKFIT->Branch("MetXFit",&MetXFit,"MetXFit/F");
   trKFIT->Branch("MetYFit",&MetYFit,"MetYFit/F");

   trKFIT->Branch("NuGen1Px",&NuGen1Px,"NuGen1Px/F");
   trKFIT->Branch("NuGen2Px",&NuGen2Px,"NuGen2Px/F");
   trKFIT->Branch("NuGen1Py",&NuGen1Py,"NuGen1Py/F");
   trKFIT->Branch("NuGen2Py",&NuGen2Py,"NuGen2Py/F");
   trKFIT->Branch("NuGen1Pz",&NuGen1Pz,"NuGen1Pz/F");
   trKFIT->Branch("NuGen2Pz",&NuGen2Pz,"NuGen2Pz/F");
   
   trKFIT->Branch("Top1Pt",&Top1Pt,"Top1Pt/F");
   trKFIT->Branch("Top1P",&Top1P,"Top1P/F");
   trKFIT->Branch("Top1Eta",&Top1Eta,"Top1Eta/F");
   trKFIT->Branch("Top1Rap",&Top1Rap,"Top1Rap/F");
   trKFIT->Branch("Top1M",&Top1M,"Top1M/F");
   trKFIT->Branch("Top2Pt",&Top2Pt,"Top2Pt/F");
   trKFIT->Branch("Top2P",&Top2P,"Top2P/F");
   trKFIT->Branch("Top2Eta",&Top2Eta,"Top2Eta/F");
   trKFIT->Branch("Top2Rap",&Top2Rap,"Top2Rap/F");
   trKFIT->Branch("Top2M",&Top2M,"Top2M/F");

   trKFIT->Branch("TopGen1Pt",&TopGen1Pt,"TopGen1Pt/F");
   trKFIT->Branch("TopGen1P",&TopGen1P,"TopGen1P/F");
   trKFIT->Branch("TopGen1Eta",&TopGen1Eta,"TopGen1Eta/F");
   trKFIT->Branch("TopGen1Rap",&TopGen1Rap,"TopGen1Rap/F");
   trKFIT->Branch("TopGen1M",&TopGen1M,"TopGen1M/F");
   trKFIT->Branch("TopGen2Pt",&TopGen2Pt,"TopGen2Pt/F");
   trKFIT->Branch("TopGen2P",&TopGen2P,"TopGen2P/F");
   trKFIT->Branch("TopGen2Eta",&TopGen2Eta,"TopGen2Eta/F");
   trKFIT->Branch("TopGen2Rap",&TopGen2Rap,"TopGen2Rap/F");
   trKFIT->Branch("TopGen2M",&TopGen2M,"TopGen2M/F");

   trKFIT->Branch("W1Pt",&W1Pt,"W1Pt/F");
   trKFIT->Branch("W1P",&W1P,"W1P/F");
   trKFIT->Branch("W1Eta",&W1Eta,"W1Eta/F");
   trKFIT->Branch("W1Phi",&W1Phi,"W1Phi/F");
   trKFIT->Branch("W1Rap",&W1Rap,"W1Rap/F");

   trKFIT->Branch("W2Pt",&W2Pt,"W2Pt/F");
   trKFIT->Branch("W2P",&W2P,"W2P/F");
   trKFIT->Branch("W2Eta",&W2Eta,"W2Eta/F");
   trKFIT->Branch("W2Phi",&W2Phi,"W2Phi/F");
   trKFIT->Branch("W2Rap",&W2Rap,"W2Rap/F");

   trKFIT->Branch("WGen1Pt",&WGen1Pt,"WGen1Pt/F");
   trKFIT->Branch("WGen1P",&WGen1P,"WGen1P/F");
   trKFIT->Branch("WGen1Eta",&WGen1Eta,"WGen1Eta/F");
   trKFIT->Branch("WGen1Phi",&WGen1Phi,"WGen1Phi/F");
   trKFIT->Branch("WGen1Rap",&WGen1Rap,"WGen1Rap/F");

   trKFIT->Branch("WGen2Pt",&WGen2Pt,"WGen2Pt/F");
   trKFIT->Branch("WGen2P",&WGen2P,"WGen2P/F");
   trKFIT->Branch("WGen2Eta",&WGen2Eta,"WGen2Eta/F");
   trKFIT->Branch("WGen2Phi",&WGen2Phi,"WGen2Phi/F");
   trKFIT->Branch("WGen2Rap",&WGen2Rap,"WGen2Rap/F");
   
   trKFIT->Branch("PtTopTop",&PtTopTop,"PtTopTop/F");
   trKFIT->Branch("PTopTop",&PTopTop,"PTopTop/F");
   trKFIT->Branch("EtaTopTop",&EtaTopTop,"EtaTopTop/F");
   trKFIT->Branch("PhiTopTop",&PhiTopTop,"PhiTopTop/F");
   trKFIT->Branch("DrTopTop",&DrTopTop,"DrTopTop/F");
   trKFIT->Branch("MTopTop",&MTopTop,"MTopTop/F");
   trKFIT->Branch("RapTopTop",&RapTopTop,"RapTopTop/F");

   trKFIT->Branch("PtTopTopGen",&PtTopTopGen,"PtTopTopGen/F");
   trKFIT->Branch("PTopTopGen",&PTopTopGen,"PTopTopGen/F");
   trKFIT->Branch("EtaTopTopGen",&EtaTopTopGen,"EtaTopTopGen/F");
   trKFIT->Branch("PhiTopTopGen",&PhiTopTopGen,"PhiTopTopGen/F");
   trKFIT->Branch("DrTopTopGen",&DrTopTopGen,"DrTopTopGen/F");
   trKFIT->Branch("MTopTopGen",&MTopTopGen,"MTopTopGen/F");
   trKFIT->Branch("RapTopTopGen",&RapTopTopGen,"RapTopTopGen/F");
   
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
	  }
	
	int nElectrons = electrons->size();
	int nMuons = muons->size();
	int nJets = jets->size();

	bool pass = (Truth_idxTopB1 >= 0 && Truth_idxTopB2 >= 0 && Truth_idxTopWL1 >= 0 && Truth_idxTopWL2 >= 0 &&
		     Truth_idxTopWNu1 >= 0 && Truth_idxTopWNu2 >= 0);
	
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
//		  if( btag > 0.800 ) isBTag = 1;
		  if( btag > 0.460 ) isBTag = 1;

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
			 (ElectronPt.size()+MuonPt.size()) == 2);

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
	     
	     P0MatchBJet1 = 0;
	     P0MatchBJet2 = 0;
	     P0RecoBJet1 = 0;
	     P0RecoBJet2 = 0;

	     P0MatchTopWL1 = 0;
	     P0MatchTopWL2 = 0;
	     P0RecoTopWL1 = 0;
	     P0RecoTopWL2 = 0;
	     
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
	     
	     Nu1Px = kf->GetNuPx(0,0);
	     Nu2Px = kf->GetNuPx(0,1);
	     Nu1Py = kf->GetNuPy(0,0);
	     Nu2Py = kf->GetNuPy(0,1);
	     Nu1Pz = kf->GetNuPz(0,0);
	     Nu2Pz = kf->GetNuPz(0,1);
	     MetXFit = kf->GetMetX(0);
	     MetYFit = kf->GetMetY(0);
	     MetXMeas = metpx;
	     MetYMeas = metpy;
	     WMass1 = kf->GetWMass(0,0);
	     WMass2 = kf->GetWMass(0,1);
	     
	     DrTopTop = kf->GetDrTopTop(0);
	     PtTopTop = kf->GetPtTopTop(0);
	     PTopTop = kf->GetPTopTop(0);
	     EtaTopTop = kf->GetEtaTopTop(0);
	     PhiTopTop = kf->GetPhiTopTop(0);
	     MTopTop = kf->GetMTopTop(0);
	     RapTopTop = kf->GetRapTopTop(0);
	     
	     Top1Pt = kf->GetTopPt(0,0);
	     Top1P = kf->GetTopP(0,0);
	     Top1Eta = kf->GetTopEta(0,0);
	     Top1M = kf->GetTopMass(0,0);
	     Top1Rap = kf->GetTopRap(0,0);

	     Top2Pt = kf->GetTopPt(0,1);
	     Top2P = kf->GetTopP(0,1);
	     Top2Eta = kf->GetTopEta(0,1);
	     Top2M = kf->GetTopMass(0,1);
	     Top2Rap = kf->GetTopRap(0,1);

	     W1Pt = kf->GetWPt(0,0);
	     W1P = kf->GetWP(0,0);
	     W1Eta = kf->GetWEta(0,0);
	     W1Phi = kf->GetWPhi(0,0);
	     W1Rap = kf->GetWRap(0,0);

	     W2Pt = kf->GetWPt(0,1);
	     W2P = kf->GetWP(0,1);
	     W2Eta = kf->GetWEta(0,1);
	     W2Phi = kf->GetWPhi(0,1);
	     W2Rap = kf->GetWRap(0,1);
	     
	     float Truth_Top1_Pt = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTop1] : -666.;
	     float Truth_Top1_Eta = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTop1] : -666.;
	     float Truth_Top1_Phi = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTop1] : -666.;
	     float Truth_Top1_E = (Truth_idxTop1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTop1] : -666.;
	     
	     float Truth_Top2_Pt = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTop2] : -666.;
	     float Truth_Top2_Eta = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTop2] : -666.;
	     float Truth_Top2_Phi = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTop2] : -666.;
	     float Truth_Top2_E = (Truth_idxTop2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTop2] : -666.;
	     
	     float Truth_TopB1_Pt = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_Eta = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_Phi = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopB1] : -666.;
	     float Truth_TopB1_E = (Truth_idxTopB1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopB1] : -666.;
	     
	     float Truth_TopB2_Pt = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_Eta = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_Phi = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopB2] : -666.;
	     float Truth_TopB2_E = (Truth_idxTopB2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopB2] : -666.;

	     float Truth_TopW1_Pt = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_Eta = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_Phi = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopW1] : -666.;
	     float Truth_TopW1_E = (Truth_idxTopW1 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopW1] : -666.;
		  
	     float Truth_TopW2_Pt = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_Eta = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_Phi = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopW2] : -666.;
	     float Truth_TopW2_E = (Truth_idxTopW2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopW2] : -666.;
	     
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
		  
	     float Truth_TopWNu2_Pt = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_pt()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_Eta = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_eta()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_Phi = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_phi()[Truth_idxTopWNu2] : -666.;
	     float Truth_TopWNu2_E = (Truth_idxTopWNu2 >= 0) ? truth->at(0).mc_truth_E()[Truth_idxTopWNu2] : -666.;

	     TLorentzVector *TopWNu1 = new TLorentzVector();
	     TLorentzVector *TopWNu2 = new TLorentzVector();

	     TopWNu1->SetPtEtaPhiE(Truth_TopWNu1_Pt,Truth_TopWNu1_Eta,Truth_TopWNu1_Phi,Truth_TopWNu1_E);
	     TopWNu2->SetPtEtaPhiE(Truth_TopWNu2_Pt,Truth_TopWNu2_Eta,Truth_TopWNu2_Phi,Truth_TopWNu2_E);
	     
	     NuGen1Px = TopWNu1->Px();
	     NuGen1Py = TopWNu1->Py();
	     NuGen1Pz = TopWNu1->Pz();

	     NuGen2Px = TopWNu2->Px();
	     NuGen2Py = TopWNu2->Py();
	     NuGen2Pz = TopWNu2->Pz();
	     
	     delete TopWNu1;
	     delete TopWNu2;

	     TLorentzVector *Top1 = new TLorentzVector();
	     TLorentzVector *Top2 = new TLorentzVector();

	     Top1->SetPtEtaPhiE(Truth_Top1_Pt,Truth_Top1_Eta,Truth_Top1_Phi,Truth_Top1_E);
	     Top2->SetPtEtaPhiE(Truth_Top2_Pt,Truth_Top2_Eta,Truth_Top2_Phi,Truth_Top2_E);
	     
	     TopGen1Pt = Top1->Pt();
	     TopGen1P = Top1->P();
	     TopGen1Eta = Top1->PseudoRapidity();
	     TopGen1Rap = Top1->Rapidity();
	     TopGen1M = Top1->M();

	     TopGen2Pt = Top2->Pt();
	     TopGen2P = Top2->P();
	     TopGen2Eta = Top2->PseudoRapidity();
	     TopGen2Rap = Top2->Rapidity();
	     TopGen2M = Top2->M();
	     
	     PtTopTopGen = (*Top1+*Top2).Pt();
	     PTopTopGen = (*Top1+*Top2).P();
	     EtaTopTopGen = (*Top1+*Top2).PseudoRapidity();
	     RapTopTopGen = (*Top1+*Top2).Rapidity();
	     PhiTopTopGen = (*Top1+*Top2).Phi();
	     MTopTopGen = (*Top1+*Top2).M();
	     DrTopTopGen = Top1->DeltaR(*Top2);
	     
	     delete Top1;
	     delete Top2;

	     TLorentzVector *TopW1 = new TLorentzVector();
	     TLorentzVector *TopW2 = new TLorentzVector();

	     TopW1->SetPtEtaPhiE(Truth_TopW1_Pt,Truth_TopW1_Eta,Truth_TopW1_Phi,Truth_TopW1_E);
	     TopW2->SetPtEtaPhiE(Truth_TopW2_Pt,Truth_TopW2_Eta,Truth_TopW2_Phi,Truth_TopW2_E);
	     
	     WMassGen1 = TopW1->M();
	     WMassGen2 = TopW2->M();
	     
	     WGen1P = TopW1->P();
	     WGen1Pt = TopW1->Pt();
	     WGen1Eta = TopW1->PseudoRapidity();
	     WGen1Rap = TopW1->Rapidity();
	     WGen1Phi = TopW1->Phi();

	     WGen2P = TopW2->P();
	     WGen2Pt = TopW2->Pt();
	     WGen2Eta = TopW2->PseudoRapidity();
	     WGen2Rap = TopW2->Rapidity();
	     WGen2Phi = TopW2->Phi();
	     
	     delete TopW1;
	     delete TopW2;
	     
	     for(int ip=0;ip<nPerm;ip++)
	       {
		  if( kf->GetDisc(ip) > 10E+9 ) continue;
		  
		  // b jets
		       
		  bool matchBJet1 = 0;
		  bool matchBJet2 = 0;
		  
		  int bjet1Idx = kf->GetIndex(BJET1_TOPTOPLEPLEP,ip);
		  int bjet2Idx = kf->GetIndex(BJET2_TOPTOPLEPLEP,ip);
		  
		  float bjet1Eta = BJetEta[bjet1Idx];
		  float bjet1Phi = BJetPhi[bjet1Idx];
		  
		  float bjet2Eta = BJetEta[bjet2Idx];
		  float bjet2Phi = BJetPhi[bjet2Idx];
		  
		  float dr_TruthTopB1_BJet1 = GetDeltaR(Truth_TopB1_Eta,Truth_TopB1_Phi,bjet1Eta,bjet1Phi);
		  float dr_TruthTopB2_BJet1 = GetDeltaR(Truth_TopB2_Eta,Truth_TopB2_Phi,bjet1Eta,bjet1Phi);
		  
		  float dr_TruthTopB1_BJet2 = GetDeltaR(Truth_TopB1_Eta,Truth_TopB1_Phi,bjet2Eta,bjet2Phi);
		  float dr_TruthTopB2_BJet2 = GetDeltaR(Truth_TopB2_Eta,Truth_TopB2_Phi,bjet2Eta,bjet2Phi);
		  
		  if( dr_TruthTopB1_BJet1 < 0.4 || dr_TruthTopB2_BJet1 < 0.4 ) matchBJet1 = 1;
		  if( dr_TruthTopB1_BJet2 < 0.4 || dr_TruthTopB2_BJet2 < 0.4 ) matchBJet2 = 1;
		  
		  bool recoBJet1 = 0;
		  bool recoBJet2 = 0;
		  
		  float dr_TruthTopB1_min = 666.;
		  float dr_TruthTopB2_min = 666.;
		  
		  for(int ij=0;ij<BJetPt.size();ij++)
		    {
		       float jeta = BJetEta[ij];
		       float jphi = BJetPhi[ij];
		       
		       float dr_TruthTopB1 = GetDeltaR(Truth_TopB1_Eta,Truth_TopB1_Phi,jeta,jphi);
		       if( dr_TruthTopB1 < dr_TruthTopB1_min ) dr_TruthTopB1_min = dr_TruthTopB1;
		       
		       float dr_TruthTopB2 = GetDeltaR(Truth_TopB2_Eta,Truth_TopB2_Phi,jeta,jphi);
		       if( dr_TruthTopB2 < dr_TruthTopB2_min ) dr_TruthTopB2_min = dr_TruthTopB2;
		    }			   
		  
		  if( dr_TruthTopB1_min < 0.4 ) recoBJet1 = 1;
		  if( dr_TruthTopB2_min < 0.4 ) recoBJet2 = 1;
		  
		  // Top leptons
		       
		  bool matchTopWL1 = 0;
		  bool matchTopWL2 = 0;
		  
		  bool topWL1isElectron = 1;
		  int topWL1Idx = kf->GetIndex(ELECTRON1_TOPTOPLEPLEP,ip);
		  if( topWL1Idx < 0 )
		    {
		       topWL1Idx = kf->GetIndex(MUON1_TOPTOPLEPLEP,ip);
		       topWL1isElectron = 0;
		    }		       
		  
		  bool topWL2isElectron = 1;
		  int topWL2Idx = kf->GetIndex(ELECTRON2_TOPTOPLEPLEP,ip);
		  if( topWL2Idx < 0 ) 
		    {
		       topWL2Idx = kf->GetIndex(MUON2_TOPTOPLEPLEP,ip);
		       topWL2isElectron = 0;
		    }
		  
		  float topWL1Eta = (topWL1isElectron) ? ElectronEta[topWL1Idx] : MuonEta[topWL1Idx];
		  float topWL1Phi = (topWL1isElectron) ? ElectronPhi[topWL1Idx] : MuonPhi[topWL1Idx];
		  
		  float topWL2Eta = (topWL2isElectron) ? ElectronEta[topWL2Idx] : MuonEta[topWL2Idx];
		  float topWL2Phi = (topWL2isElectron) ? ElectronPhi[topWL2Idx] : MuonPhi[topWL2Idx];
		  
		  float dr_TruthTopWL1_TopWL1 = GetDeltaR(Truth_TopWL1_Eta,Truth_TopWL1_Phi,topWL1Eta,topWL1Phi);
		  float dr_TruthTopWL2_TopWL1 = GetDeltaR(Truth_TopWL2_Eta,Truth_TopWL2_Phi,topWL1Eta,topWL1Phi);
		  
		  float dr_TruthTopWL1_TopWL2 = GetDeltaR(Truth_TopWL1_Eta,Truth_TopWL1_Phi,topWL2Eta,topWL2Phi);
		  float dr_TruthTopWL2_TopWL2 = GetDeltaR(Truth_TopWL2_Eta,Truth_TopWL2_Phi,topWL2Eta,topWL2Phi);
		  
		  if( (dr_TruthTopWL1_TopWL1 < 0.1 || dr_TruthTopWL1_TopWL2 < 0.1) &&
		      (dr_TruthTopB1_BJet1 < 0.4 || dr_TruthTopB1_BJet2 < 0.4) ) matchTopWL1 = 1;
		  if( (dr_TruthTopWL2_TopWL1 < 0.1 || dr_TruthTopWL2_TopWL2 < 0.1) &&
		      (dr_TruthTopB2_BJet1 < 0.4 || dr_TruthTopB2_BJet2 < 0.4) ) matchTopWL2 = 1;
		  
		  bool recoTopWL1 = 0;
		  bool recoTopWL2 = 0;
		  
		  float dr_TruthTopWL1_min = 666.;
		  float dr_TruthTopWL2_min = 666.;
		  
		  if( topWL1isElectron )
		    {			    
		       for(int ie=0;ie<ElectronPt.size();ie++)
			 {
			    float eeta = ElectronEta[ie];
			    float ephi = ElectronPhi[ie];
			    
			    float dr_TruthTopWL1 = GetDeltaR(Truth_TopWL1_Eta,Truth_TopWL1_Phi,eeta,ephi);
			    if( dr_TruthTopWL1 < dr_TruthTopWL1_min ) dr_TruthTopWL1_min = dr_TruthTopWL1;
			 }
		    }		       
		  else
		    {
		       for(int im=0;im<MuonPt.size();im++)
			 {
			    float meta = MuonEta[im];
			    float mphi = MuonPhi[im];
			    
			    float dr_TruthTopWL1 = GetDeltaR(Truth_TopWL1_Eta,Truth_TopWL1_Phi,meta,mphi);
			    if( dr_TruthTopWL1 < dr_TruthTopWL1_min ) dr_TruthTopWL1_min = dr_TruthTopWL1;
			 }
		    }		       
		  
		  if( topWL2isElectron )
		    {			    
		       for(int ie=0;ie<ElectronPt.size();ie++)
			 {
			    float eeta = ElectronEta[ie];
			    float ephi = ElectronPhi[ie];
			    
			    float dr_TruthTopWL2 = GetDeltaR(Truth_TopWL2_Eta,Truth_TopWL2_Phi,eeta,ephi);
			    if( dr_TruthTopWL2 < dr_TruthTopWL2_min ) dr_TruthTopWL2_min = dr_TruthTopWL2;
			 }
		    }		       
		  else
		    {
		       for(int im=0;im<MuonPt.size();im++)
			 {
			    float meta = MuonEta[im];
			    float mphi = MuonPhi[im];
			    
			    float dr_TruthTopWL2 = GetDeltaR(Truth_TopWL2_Eta,Truth_TopWL2_Phi,meta,mphi);
			    if( dr_TruthTopWL2 < dr_TruthTopWL2_min ) dr_TruthTopWL2_min = dr_TruthTopWL2;
			 }
		    }		       
		  
		  if( dr_TruthTopWL1_min < 0.1 ) recoTopWL1 = 1;
		  if( dr_TruthTopWL2_min < 0.1 ) recoTopWL2 = 1;
		  
		  if( ip == 0 )
		    {
		       P0MatchBJet1 = matchBJet1;
		       P0MatchBJet2 = matchBJet2;
		       
		       P0RecoBJet1 = recoBJet1;
		       P0RecoBJet2 = recoBJet2;
		       
		       P0MatchTopWL1 = matchTopWL1;
		       P0MatchTopWL2 = matchTopWL2;
		       
		       P0RecoTopWL1 = recoTopWL1;
		       P0RecoTopWL2 = recoTopWL2;
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
