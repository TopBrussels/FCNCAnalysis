#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);

void setTDRStyle();
double BW(double* x,double* par);

// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookHowToFit

void draw_gen()
{
   gROOT->SetBatch();

   setTDRStyle();
   
   TFile *f = TFile::Open("run_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8/data.root");

   float Top1M, Top2M;
   float TopW1M, TopW2M;
   float TopWHadRecM, TopHadRecM, TopLepRecM;
   float dMetPx, dMetPy;
   float dBJet1Px, dBJet1Py, dBJet1Pz, dBJet1E;
   float dBJet2Px, dBJet2Py, dBJet2Pz, dBJet2E;
   float dElecPx, dElecPy, dElecPz, dElecE;
   float dMuonPx, dMuonPy, dMuonPz, dMuonE;
   float dNonBJet1Px, dNonBJet1Py, dNonBJet1Pz, dNonBJet1E;
   float dNonBJet2Px, dNonBJet2Py, dNonBJet2Pz, dNonBJet2E;

   TTree *tr = (TTree*)f->Get("trGEN");
   tr->SetBranchAddress("TopWHadRecM",&TopWHadRecM);
   tr->SetBranchAddress("TopHadRecM",&TopHadRecM);
   tr->SetBranchAddress("TopLepRecM",&TopLepRecM);
   tr->SetBranchAddress("Top1M",&Top1M);
   tr->SetBranchAddress("Top2M",&Top2M);
   tr->SetBranchAddress("TopW1M",&TopW1M);
   tr->SetBranchAddress("TopW2M",&TopW2M);
   tr->SetBranchAddress("dMetPx",&dMetPx);
   tr->SetBranchAddress("dMetPy",&dMetPy);
   tr->SetBranchAddress("dBJet1Px",&dBJet1Px);
   tr->SetBranchAddress("dBJet1Py",&dBJet1Py);
   tr->SetBranchAddress("dBJet1Pz",&dBJet1Pz);
   tr->SetBranchAddress("dBJet1E",&dBJet1E);
   tr->SetBranchAddress("dBJet2Px",&dBJet2Px);
   tr->SetBranchAddress("dBJet2Py",&dBJet2Py);
   tr->SetBranchAddress("dBJet2Pz",&dBJet2Pz);
   tr->SetBranchAddress("dBJet2E",&dBJet2E);
   tr->SetBranchAddress("dElecPx",&dElecPx);
   tr->SetBranchAddress("dElecPy",&dElecPy);
   tr->SetBranchAddress("dElecPz",&dElecPz);
   tr->SetBranchAddress("dElecE",&dElecE);
   tr->SetBranchAddress("dMuonPx",&dMuonPx);
   tr->SetBranchAddress("dMuonPy",&dMuonPy);
   tr->SetBranchAddress("dMuonPz",&dMuonPz);
   tr->SetBranchAddress("dMuonE",&dMuonE);
   tr->SetBranchAddress("dNonBJet1Px",&dNonBJet1Px);
   tr->SetBranchAddress("dNonBJet1Py",&dNonBJet1Py);
   tr->SetBranchAddress("dNonBJet1Pz",&dNonBJet1Pz);
   tr->SetBranchAddress("dNonBJet1E",&dNonBJet1E);
   tr->SetBranchAddress("dNonBJet2Px",&dNonBJet2Px);
   tr->SetBranchAddress("dNonBJet2Py",&dNonBJet2Py);
   tr->SetBranchAddress("dNonBJet2Pz",&dNonBJet2Pz);
   tr->SetBranchAddress("dNonBJet2E",&dNonBJet2E);

   int nHist = 0;
   
   TH1D *h[1000];
   std::string hName[1000];
   std::string hLab[1000];

   TFile *fOut = new TFile("pdf.root","RECREATE");
   
   h[nHist] = new TH1D("h_TopM","h_TopM",100,150.,190.);
   h[nHist]->Sumw2();
   hName[nHist] = "TopM";
   hLab[nHist] = "m(t) [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_TopWM","h_TopWM",100,60.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "TopWM";
   hLab[nHist] = "m(W) [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dMetPx","h_dMetPx",100,-100.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMetPx";
   hLab[nHist] = "MetPx_{gen} - MetPx_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dMetPy","h_dMetPy",100,-100.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMetPy";
   hLab[nHist] = "MetPy_{gen} - MetPy_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dBJetPx","h_dBJetPx",100,-100.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "dBJetPx";
   hLab[nHist] = "BJetPx_{gen} - BJetPx_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dBJetPy","h_dBJetPy",100,-100.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "dBJetPy";
   hLab[nHist] = "BJetPy_{gen} - BJetPy_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dBJetPz","h_dBJetPz",100,-100.,100.);
   h[nHist]->Sumw2();
   hName[nHist] = "dBJetPz";
   hLab[nHist] = "BJetPz_{gen} - BJetPz_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dBJetE","h_dBJetE",100,-200.,200.);
   h[nHist]->Sumw2();
   hName[nHist] = "dBJetE";
   hLab[nHist] = "BJetE_{gen} - BJetE_{reco} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_dElecPx","h_dElecPx",100,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "dElecPx";
   hLab[nHist] = "ElecPx_{gen} - ElecPx_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dElecPy","h_dElecPy",100,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "dElecPy";
   hLab[nHist] = "ElecPy_{gen} - ElecPy_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dElecPz","h_dElecPz",100,-6.,6.);
   h[nHist]->Sumw2();
   hName[nHist] = "dElecPz";
   hLab[nHist] = "ElecPz_{gen} - ElecPz_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dElecE","h_dElecE",100,-6.,6.);
   h[nHist]->Sumw2();
   hName[nHist] = "dElecE";
   hLab[nHist] = "ElecE_{gen} - ElecE_{reco} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_dMuonPx","h_dMuonPx",100,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMuonPx";
   hLab[nHist] = "MuonPx_{gen} - MuonPx_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dMuonPy","h_dMuonPy",100,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMuonPy";
   hLab[nHist] = "MuonPy_{gen} - MuonPy_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dMuonPz","h_dMuonPz",100,-6.,6.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMuonPz";
   hLab[nHist] = "MuonPz_{gen} - MuonPz_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dMuonE","h_dMuonE",100,-6.,6.);
   h[nHist]->Sumw2();
   hName[nHist] = "dMuonE";
   hLab[nHist] = "MuonE_{gen} - MuonE_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dNonBJetPx","h_dNonBJetPx",100,-70.,70.);
   h[nHist]->Sumw2();
   hName[nHist] = "dNonBJetPx";
   hLab[nHist] = "NonBJetPx_{gen} - NonBJetPx_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dNonBJetPy","h_dNonBJetPy",100,-70.,70.);
   h[nHist]->Sumw2();
   hName[nHist] = "dNonBJetPy";
   hLab[nHist] = "NonBJetPy_{gen} - NonBJetPy_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dNonBJetPz","h_dNonBJetPz",100,-70.,70.);
   h[nHist]->Sumw2();
   hName[nHist] = "dNonBJetPz";
   hLab[nHist] = "NonBJetPz_{gen} - NonBJetPz_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_dNonBJetE","h_dNonBJetE",100,-70.,70.);
   h[nHist]->Sumw2();
   hName[nHist] = "dNonBJetE";
   hLab[nHist] = "NonBJetE_{gen} - NonBJetE_{reco} [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_TopWHadRecM","h_TopWHadRecM",100,20.,140.);
   h[nHist]->Sumw2();
   hName[nHist] = "TopWHadRecM";
   hLab[nHist] = "m(W) [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_TopHadRecM","h_TopHadRecM",100,0.,300.);
   h[nHist]->Sumw2();
   hName[nHist] = "TopHadRecM";
   hLab[nHist] = "m(t) [GeV]";
   nHist++;

   h[nHist] = new TH1D("h_TopLepRecM","h_TopLepRecM",100,80.,260.);
   h[nHist]->Sumw2();
   hName[nHist] = "TopLepRecM";
   hLab[nHist] = "m(t) [GeV]";
   nHist++;
   
   int nev = tr->GetEntries();   

   TRandom3 *rnd = new TRandom3();
   
   for(int i=0;i<nev;i++)
     {
	tr->GetEntry(i);
	
	for(int ih=0;ih<nHist;ih++)
	  {
	     if( hName[ih] == "TopM" ) {h[ih]->Fill(Top1M);h[ih]->Fill(Top2M);}
	     else if( hName[ih] == "TopWM" ) {h[ih]->Fill(TopW1M);h[ih]->Fill(TopW2M);}
	     else if( hName[ih] == "dMetPx" ) h[ih]->Fill(dMetPx);
	     else if( hName[ih] == "dMetPy" ) h[ih]->Fill(dMetPy);
	     else if( hName[ih] == "dBJetPx" ) {h[ih]->Fill(dBJet1Px);h[ih]->Fill(dBJet2Px);}
	     else if( hName[ih] == "dBJetPy" ) {h[ih]->Fill(dBJet1Py);h[ih]->Fill(dBJet2Py);}
	     else if( hName[ih] == "dBJetPz" ) {h[ih]->Fill(dBJet1Pz);h[ih]->Fill(dBJet2Pz);}
	     else if( hName[ih] == "dBJetE" ) {h[ih]->Fill(dBJet1E);h[ih]->Fill(dBJet2E);}
	     else if( hName[ih] == "dElecPx" ) {h[ih]->Fill(dElecPx);}
	     else if( hName[ih] == "dElecPy" ) {h[ih]->Fill(dElecPy);}
	     else if( hName[ih] == "dElecPz" ) {h[ih]->Fill(dElecPz);}
	     else if( hName[ih] == "dElecE" ) {h[ih]->Fill(dElecE);}
	     else if( hName[ih] == "dMuonPx" ) {h[ih]->Fill(dMuonPx);}
	     else if( hName[ih] == "dMuonPy" ) {h[ih]->Fill(dMuonPy);}
	     else if( hName[ih] == "dMuonPz" ) {h[ih]->Fill(dMuonPz);}
	     else if( hName[ih] == "dMuonE" ) {h[ih]->Fill(dMuonE);}
	     else if( hName[ih] == "dNonBJetPx" ) {h[ih]->Fill(dNonBJet1Px); h[ih]->Fill(dNonBJet2Px);}
	     else if( hName[ih] == "dNonBJetPy" ) {h[ih]->Fill(dNonBJet1Py); h[ih]->Fill(dNonBJet2Py);}
	     else if( hName[ih] == "dNonBJetPz" ) {h[ih]->Fill(dNonBJet1Pz); h[ih]->Fill(dNonBJet2Pz);}
	     else if( hName[ih] == "dNonBJetE" ) {h[ih]->Fill(dNonBJet1E); h[ih]->Fill(dNonBJet2E);}
	     else if( hName[ih] == "TopWHadRecM" ) {h[ih]->Fill(TopWHadRecM);}
	     else if( hName[ih] == "TopHadRecM" ) {h[ih]->Fill(TopHadRecM);}
	     else if( hName[ih] == "TopLepRecM" ) {h[ih]->Fill(TopLepRecM);}
	  }	
     }
   
   delete rnd;

   // Plots
   
   TCanvas *c1 = new TCanvas();
   c1->Draw();
   c1->cd();

   TPad *c1_1;
   
   gStyle->SetHistTopMargin(0);

   for(int i=0;i<nHist;i++)
     {	
//	addbin(h[i]);

	h[i]->Scale(1./h[i]->Integral());
	h[i]->Scale(1./h[i]->GetMaximum());
   
	h[i]->SetLineWidth(2);	
	h[i]->SetLineColor(kRed);
	h[i]->SetMarkerColor(kRed);
	h[i]->SetMarkerStyle(20);
	h[i]->Draw("hist e1");
	h[i]->GetXaxis()->SetTitle(hLab[i].c_str());
//	h[i]->GetYaxis()->SetTitle("Normalized to unity");
	
	float max = h[i]->GetMaximum();
	
	h[i]->SetMaximum(1.2*max);
	
	if( hName[i] == "TopWM" || hName[i] == "TopM" )
	  {
	     std::string funcName = hName[i]+"_Fit";
	     TF1 *func = new TF1(funcName.c_str(),BW,h[i]->GetXaxis()->GetBinLowEdge(1),
				 h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()),3);
	     
	     double mean = 80.4;
	     if( hName[i] == "TopM" ) mean = 173.0;
	     double sigma = 2.0;
	     if( hName[i] == "TopM" ) sigma = 2.0;
	     
	     func->SetParameter(0,mean);
	     func->SetParName(0,"mean");	     
	     func->SetParameter(1,sigma);
	     func->SetParName(1,"sigma");
	     func->SetParameter(2,1.0);
	     func->SetParName(2,"constant");
	     
	     func->FixParameter(2,1.);
	     
	     h[i]->Fit(funcName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();

	     std::string funcGausName = hName[i]+"_Gaus";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     funcGaus->Write();
	  }	

	if( hName[i] == "dMetPx" || hName[i] == "dMetPy" )
	  {
	     double mean = 0.;
	     double sigma = 20.;

	     std::string funcGausName = hName[i]+"_Gaus";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     funcGaus->FixParameter(0,1);
	     h[i]->Fit(funcGausName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcGausName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();
	  }	

	if( hName[i] == "dBJetPx" || hName[i] == "dBJetPy" || hName[i] == "dBJetPz" || hName[i] == "dBJetE" ||
	    hName[i] == "dNonBJetPx" || hName[i] == "dNonBJetPy" || hName[i] == "dNonBJetPz" || hName[i] == "dNonBJetE" )
	  {
	     double mean = 0.;
	     double sigma = 20.;

	     std::string funcGausName = hName[i]+"_Fit";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)+gaus(3)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     if( hName[i] == "dBJetPz" ) funcGaus->FixParameter(0,0.7);
	     if( hName[i] == "dNonBJetPx" || hName[i] == "dNonBJetPy" ) funcGaus->FixParameter(0,0.5);
	     if( hName[i] == "dNonBJetPz" ) funcGaus->FixParameter(0,0.2);

	     funcGaus->SetParameter(3,0.3);
	     funcGaus->SetParameter(4,mean);
	     funcGaus->SetParameter(5,sigma*2.);
	     funcGaus->FixParameter(3,0.3);	     
	     if( hName[i] == "dBJetPz" ) funcGaus->FixParameter(3,0.3);
	     if( hName[i] == "dNonBJetPx" || hName[i] == "dNonBJetPy" ) funcGaus->FixParameter(3,0.5);
	     if( hName[i] == "dNonBJetPz" ) funcGaus->FixParameter(3,0.8);
	     h[i]->Fit(funcGausName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcGausName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();
	  }	

	if( hName[i] == "dElecPx" || hName[i] == "dElecPy" || hName[i] == "dElecPz" || hName[i] == "dElecE" )
	  {
	     double mean = 0.;
	     double sigma = 1.;

	     std::string funcGausName = hName[i]+"_Fit";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)+gaus(3)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     funcGaus->FixParameter(0,0.7);
	     if( hName[i] == "dElecPz" ) funcGaus->FixParameter(0,0.8);

	     funcGaus->SetParameter(3,0.3);
	     funcGaus->SetParameter(4,mean);
	     funcGaus->SetParameter(5,sigma*2.);
	     funcGaus->FixParameter(3,0.3);
	     if( hName[i] == "dElecPz" ) funcGaus->FixParameter(3,0.2);
	     h[i]->Fit(funcGausName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcGausName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();
	  }	

	if( hName[i] == "dMuonPx" || hName[i] == "dMuonPy" || hName[i] == "dMuonPz" || hName[i] == "dMuonE" )
	  {
	     double mean = 0.;
	     double sigma = 1.;

	     std::string funcGausName = hName[i]+"_Fit";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)+gaus(3)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     funcGaus->FixParameter(0,0.7);
	     if( hName[i] == "dMuonPz" ) funcGaus->FixParameter(0,0.8);

	     funcGaus->SetParameter(3,0.3);
	     funcGaus->SetParameter(4,mean);
	     funcGaus->SetParameter(5,sigma*2.);
	     funcGaus->FixParameter(3,0.3);
	     if( hName[i] == "dMuonPz" ) funcGaus->FixParameter(3,0.2);
	     h[i]->Fit(funcGausName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcGausName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();
	  }	

	if( hName[i] == "TopWHadRecM" || hName[i] == "TopHadRecM" || hName[i] == "TopLepRecM" )
	  {
	     double mean = 80.;
	     double sigma = 10.;
	     
	     if( hName[i] == "TopHadRecM" || hName[i] == "TopLepRecM" )
	       {
		  mean = 170;
		  sigma = 20;
	       }	     

	     std::string funcGausName = hName[i]+"_Fit";
	     TF1 *funcGaus = new TF1(funcGausName.c_str(),"gaus(0)+gaus(3)",h[i]->GetXaxis()->GetBinLowEdge(1),
				     h[i]->GetXaxis()->GetBinUpEdge(h[i]->GetXaxis()->GetNbins()));
	     funcGaus->SetParameter(0,1);
	     funcGaus->SetParameter(1,mean);
	     funcGaus->SetParameter(2,sigma);
	     if( hName[i] != "TopLepRecM" ) funcGaus->FixParameter(0,0.7);

	     funcGaus->SetParameter(3,0.3);
	     if( hName[i] == "TopLepRecM" ) funcGaus->SetParameter(3,0.7);
	     if( hName[i] == "TopLepRecM" ) funcGaus->SetParameter(4,150);
	     funcGaus->SetParameter(4,mean);
	     funcGaus->SetParameter(5,sigma*2.);

	     h[i]->Fit(funcGausName.c_str(),"QR");
	     TF1 *fit = h[i]->GetFunction(funcGausName.c_str());
	     fit->SetLineColor(1);
	     fit->Draw("same");
	     fit->Write();
	  }	
	
	std::string figName = "pics/"+hName[i]+".eps";
	c1->Print(figName.c_str());
	c1->Clear();
     }   

   fOut->Write();
   fOut->Close();
   
   gApplication->Terminate();
}

void addbin(TH1D *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   // Set overflow and underflow bins to 0
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

void tdrGrid(bool gridOn) {
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(450); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
   tdrStyle->SetStatX(0.9);
   tdrStyle->SetStatY(0.9);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.17);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.06);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.00);
  tdrStyle->SetTitleYOffset(1.00);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

double BW(double* x,double* par)
{
   // mWmean = par[0]
   // GammaW = par[1]
   
   float f = par[2]*par[0]*par[0]*par[1]*par[1]/(pow(x[0]*x[0]-par[0]*par[0],2)+par[0]*par[0]*par[1]*par[1]);
   
   return f;
}
