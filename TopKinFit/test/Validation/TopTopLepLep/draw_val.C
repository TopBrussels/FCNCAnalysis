#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);

void setTDRStyle();

void draw_val()
{
   gROOT->SetBatch();
   
   setTDRStyle();
   
   TFile *fSig = TFile::Open("runTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8/data.root");

   float P0Disc;
   
   float Nu1Px, Nu1Py, Nu1Pz;
   float Nu2Px, Nu2Py, Nu2Pz;

   float NuGen1Px, NuGen1Py, NuGen1Pz;
   float NuGen2Px, NuGen2Py, NuGen2Pz;
   
   float MTopTop, PTopTop, PtTopTop, EtaTopTop, PhiTopTop, DrTopTop;
   float MTopTopGen, PTopTopGen, PtTopTopGen, EtaTopTopGen, PhiTopTopGen, DrTopTopGen;

   float Top1P, Top1Pt, Top1Eta, Top1M;
   float Top2P, Top2Pt, Top2Eta, Top2M;
   float W1P, W1Pt, W1Eta, W1Phi;
   float W2P, W2Pt, W2Eta, W2Phi;

   float TopGen1P, TopGen1Pt, TopGen1Eta, TopGen1M;
   float TopGen2P, TopGen2Pt, TopGen2Eta, TopGen2M;
   float WGen1P, WGen1Pt, WGen1Eta, WGen1Phi;
   float WGen2P, WGen2Pt, WGen2Eta, WGen2Phi;
   
   float WMass1, WMass2;
   float WMassGen1, WMassGen2;
   
   TTree *trSig = (TTree*)fSig->Get("trKFIT");

   trSig->SetBranchAddress("P0Disc",&P0Disc);
   
   trSig->SetBranchAddress("Nu1Px",&Nu1Px);
   trSig->SetBranchAddress("Nu1Py",&Nu1Py);
   trSig->SetBranchAddress("Nu1Pz",&Nu1Pz);
   trSig->SetBranchAddress("Nu2Px",&Nu2Px);
   trSig->SetBranchAddress("Nu2Py",&Nu2Py);
   trSig->SetBranchAddress("Nu2Pz",&Nu2Pz);

   trSig->SetBranchAddress("NuGen1Px",&NuGen1Px);
   trSig->SetBranchAddress("NuGen1Py",&NuGen1Py);
   trSig->SetBranchAddress("NuGen1Pz",&NuGen1Pz);
   trSig->SetBranchAddress("NuGen2Px",&NuGen2Px);
   trSig->SetBranchAddress("NuGen2Py",&NuGen2Py);
   trSig->SetBranchAddress("NuGen2Pz",&NuGen2Pz);
   
   trSig->SetBranchAddress("MTopTop",&MTopTop);
   trSig->SetBranchAddress("PtTopTop",&PtTopTop);
   trSig->SetBranchAddress("PTopTop",&PTopTop);
   trSig->SetBranchAddress("EtaTopTop",&EtaTopTop);
   trSig->SetBranchAddress("PhiTopTop",&PhiTopTop);
   trSig->SetBranchAddress("DrTopTop",&DrTopTop);
   
   trSig->SetBranchAddress("MTopTopGen",&MTopTopGen);
   trSig->SetBranchAddress("PtTopTopGen",&PtTopTopGen);
   trSig->SetBranchAddress("PTopTopGen",&PTopTopGen);
   trSig->SetBranchAddress("EtaTopTopGen",&EtaTopTopGen);
   trSig->SetBranchAddress("PhiTopTopGen",&PhiTopTopGen);
   trSig->SetBranchAddress("DrTopTopGen",&DrTopTopGen);

   trSig->SetBranchAddress("Top1P",&Top1P);
   trSig->SetBranchAddress("Top1Pt",&Top1Pt);
   trSig->SetBranchAddress("Top1Eta",&Top1Eta);
   trSig->SetBranchAddress("Top1M",&Top1M);

   trSig->SetBranchAddress("Top2P",&Top2P);
   trSig->SetBranchAddress("Top2Pt",&Top2Pt);
   trSig->SetBranchAddress("Top2Eta",&Top2Eta);
   trSig->SetBranchAddress("Top2M",&Top2M);

   trSig->SetBranchAddress("TopGen1P",&TopGen1P);
   trSig->SetBranchAddress("TopGen1Pt",&TopGen1Pt);
   trSig->SetBranchAddress("TopGen1Eta",&TopGen1Eta);
   trSig->SetBranchAddress("TopGen1M",&TopGen1M);

   trSig->SetBranchAddress("TopGen2P",&TopGen2P);
   trSig->SetBranchAddress("TopGen2Pt",&TopGen2Pt);
   trSig->SetBranchAddress("TopGen2Eta",&TopGen2Eta);
   trSig->SetBranchAddress("TopGen2M",&TopGen2M);

   trSig->SetBranchAddress("W1P",&W1P);
   trSig->SetBranchAddress("W1Pt",&W1Pt);
   trSig->SetBranchAddress("W1Eta",&W1Eta);
   trSig->SetBranchAddress("W1Phi",&W1Phi);

   trSig->SetBranchAddress("W2P",&W2P);
   trSig->SetBranchAddress("W2Pt",&W2Pt);
   trSig->SetBranchAddress("W2Eta",&W2Eta);
   trSig->SetBranchAddress("W2Phi",&W2Phi);

   trSig->SetBranchAddress("WGen1P",&WGen1P);
   trSig->SetBranchAddress("WGen1Pt",&WGen1Pt);
   trSig->SetBranchAddress("WGen1Eta",&WGen1Eta);
   trSig->SetBranchAddress("WGen1Phi",&WGen1Phi);

   trSig->SetBranchAddress("WGen2P",&WGen2P);
   trSig->SetBranchAddress("WGen2Pt",&WGen2Pt);
   trSig->SetBranchAddress("WGen2Eta",&WGen2Eta);
   trSig->SetBranchAddress("WGen2Phi",&WGen2Phi);

   trSig->SetBranchAddress("WMass1",&WMass1);
   trSig->SetBranchAddress("WMass2",&WMass2);
   
   trSig->SetBranchAddress("WMassGen1",&WMassGen1);
   trSig->SetBranchAddress("WMassGen2",&WMassGen2);
   
   int nHist = 0;
   
   TH1D *hSig[1000];
   TH1D *hSigGen[1000];
   std::string hName[1000];
   std::string hLab[1000];
   
   // Nu1
   hSig[nHist] = new TH1D("hSig_PxNu1","hSig_PxNu1",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PxNu1","hSigGen_PxNu1",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu1Px";
   hLab[nHist] = "p_{x} (#nu_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PyNu1","hSig_PyNu1",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PyNu1","hSigGen_PyNu1",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu1Py";
   hLab[nHist] = "p_{y} (#nu_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PzNu1","hSig_PzNu1",30,-500.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PzNu1","hSigGen_PzNu1",30,-500.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu1Pz";
   hLab[nHist] = "p_{z} (#nu_{1}) [GeV]";
   nHist++;

   // Nu2
   hSig[nHist] = new TH1D("hSig_PxNu2","hSig_PxNu2",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PxNu2","hSigGen_PxNu2",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu2Px";
   hLab[nHist] = "p_{x} (#nu_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PyNu2","hSig_PyNu2",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PyNu2","hSigGen_PyNu2",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu2Py";
   hLab[nHist] = "p_{y} (#nu_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PzNu2","hSig_PzNu2",30,-500.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PzNu2","hSigGen_PzNu2",30,-500.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Nu2Pz";
   hLab[nHist] = "p_{z} (#nu_{2}) [GeV]";
   nHist++;

   // TopTop
   
   hSig[nHist] = new TH1D("hSig_MTopTop","hSig_MTopTop",30,0.,2000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_MTopTop","hSigGen_MTopTop",30,0.,2000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "MTopTop";
   hLab[nHist] = "m(t#bar{t}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PtTopTop","hSig_PtTopTop",30,0.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PtTopTop","hSigGen_PtTopTop",30,0.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "PtTopTop";
   hLab[nHist] = "p_{T}(t#bar{t}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PTopTop","hSig_PTopTop",30,0.,2000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PTopTop","hSigGen_PTopTop",30,0.,2000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "PTopTop";
   hLab[nHist] = "p(t#bar{t}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_EtaTopTop","hSig_EtaTopTop",30,-6.,6.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_EtaTopTop","hSigGen_EtaTopTop",30,-6.,6.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "EtaTopTop";
   hLab[nHist] = "#eta(t#bar{t})";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PhiTopTop","hSig_PhiTopTop",30,-3.5,3.5);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PhiTopTop","hSigGen_PhiTopTop",30,-3.5,3.5);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "PhiTopTop";
   hLab[nHist] = "#phi(t#bar{t})";
   nHist++;

   hSig[nHist] = new TH1D("hSig_DrTopTop","hSig_DrTopTop",30,0.,6.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_DrTopTop","hSigGen_DrTopTop",30,0.,6.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "DrTopTop";
   hLab[nHist] = "#Delta R (t,#bar{t})";
   nHist++;

   // W
   
   hSig[nHist] = new TH1D("hSig_WMass1","hSig_WMass1",30,60.,110.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_WMass1","hSigGen_WMass1",30,60.,110.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "WMass1";
   hLab[nHist] = "m(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WMass2","hSig_WMass2",30,60.,110.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_WMass2","hSigGen_WMass2",30,60.,110.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "WMass2";
   hLab[nHist] = "m(W_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W1P","hSig_W1P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_W1P","hSigGen_W1P",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "W1P";
   hLab[nHist] = "p(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W2P","hSig_W2P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_W2P","hSigGen_W2P",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "W2P";
   hLab[nHist] = "p(W_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W1Pt","hSig_W1Pt",30,0.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_W1Pt","hSigGen_W1Pt",30,0.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "W1Pt";
   hLab[nHist] = "p_{T}(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W2Pt","hSig_W2Pt",30,0.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_W2Pt","hSigGen_W2Pt",30,0.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "W2Pt";
   hLab[nHist] = "p_{T}(W_{2}) [GeV]";
   nHist++;

   // Top
   
   hSig[nHist] = new TH1D("hSig_Top1M","hSig_Top1M",30,150.,200.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top1M","hSigGen_Top1M",30,150.,200.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top1M";
   hLab[nHist] = "m(t_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top2M","hSig_Top2M",30,150.,200.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top2M","hSigGen_Top2M",30,150.,200.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top2M";
   hLab[nHist] = "m(t_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top1P","hSig_Top1P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top1P","hSigGen_Top1P",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top1P";
   hLab[nHist] = "p(t_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top2P","hSig_Top2P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top2P","hSigGen_Top2P",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top2P";
   hLab[nHist] = "p(t_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top1Eta","hSig_Top1Eta",30,-5.,5.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top1Eta","hSigGen_Top1Eta",30,-5.,5.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top1Eta";
   hLab[nHist] = "#eta(t_{1})";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top2Eta","hSig_Top2Eta",30,-5.,5.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_Top2Eta","hSigGen_Top2Eta",30,-5.,5.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "Top2Eta";
   hLab[nHist] = "#eta(t_{2})";
   nHist++;
   
   int nevSig = trSig->GetEntries();   
   
   for(int i=0;i<nevSig;i++)
     {
	trSig->GetEntry(i);

	if( P0Disc > 10E+9 ) continue;
//	if( P0Disc > 1 ) continue;
	
	for(int ih=0;ih<nHist;ih++)
	  {	     
	     if( hName[ih] == "Nu1Px" ) 
	       {
		  hSig[ih]->Fill(Nu1Px);
		  hSigGen[ih]->Fill(NuGen1Px);
	       }	     
	     else if( hName[ih] == "Nu1Py" ) 
	       {
		  hSig[ih]->Fill(Nu1Py);
		  hSigGen[ih]->Fill(NuGen1Py);
	       }	     
	     else if( hName[ih] == "Nu1Pz" ) 
	       {
		  hSig[ih]->Fill(Nu1Pz);
		  hSigGen[ih]->Fill(NuGen1Pz);
	       }	     
	     else if( hName[ih] == "Nu2Px" ) 
	       {
		  hSig[ih]->Fill(Nu2Px);
		  hSigGen[ih]->Fill(NuGen2Px);
	       }	     
	     else if( hName[ih] == "Nu2Py" ) 
	       {
		  hSig[ih]->Fill(Nu2Py);
		  hSigGen[ih]->Fill(NuGen2Py);
	       }	     
	     else if( hName[ih] == "Nu2Pz" ) 
	       {
		  hSig[ih]->Fill(Nu2Pz);
		  hSigGen[ih]->Fill(NuGen2Pz);
	       }	     

	     else if( hName[ih] == "MTopTop" )
	       {
		  hSig[ih]->Fill(MTopTop);
		  hSigGen[ih]->Fill(MTopTopGen);
	       }	     
	     else if( hName[ih] == "PtTopTop" )
	       {
		  hSig[ih]->Fill(PtTopTop);
		  hSigGen[ih]->Fill(PtTopTopGen);
	       }	     
	     else if( hName[ih] == "PTopTop" )
	       {
		  hSig[ih]->Fill(PTopTop);
		  hSigGen[ih]->Fill(PTopTopGen);
	       }	     
	     else if( hName[ih] == "EtaTopTop" )
	       {
		  hSig[ih]->Fill(EtaTopTop);
		  hSigGen[ih]->Fill(EtaTopTopGen);
	       }	     
	     else if( hName[ih] == "PhiTopTop" )
	       {
		  hSig[ih]->Fill(PhiTopTop);
		  hSigGen[ih]->Fill(PhiTopTopGen);
	       }	     
	     else if( hName[ih] == "DrTopTop" )
	       {
		  hSig[ih]->Fill(DrTopTop);
		  hSigGen[ih]->Fill(DrTopTopGen);
	       }	     

	     else if( hName[ih] == "WMass1" )
	       {
		  hSig[ih]->Fill(WMass1);
		  hSigGen[ih]->Fill(WMassGen1);
	       }	     
	     else if( hName[ih] == "WMass2" )
	       {
		  hSig[ih]->Fill(WMass2);
		  hSigGen[ih]->Fill(WMassGen2);
	       }	     
	     else if( hName[ih] == "W1P" )
	       {
		  hSig[ih]->Fill(W1P);
		  hSigGen[ih]->Fill(WGen1P);
	       }	     
	     else if( hName[ih] == "W2P" )
	       {
		  hSig[ih]->Fill(W2P);
		  hSigGen[ih]->Fill(WGen2P);
	       }	     
	     else if( hName[ih] == "W1Pt" )
	       {
		  hSig[ih]->Fill(W1Pt);
		  hSigGen[ih]->Fill(WGen1Pt);
	       }	     
	     else if( hName[ih] == "W2Pt" )
	       {
		  hSig[ih]->Fill(W2Pt);
		  hSigGen[ih]->Fill(WGen2Pt);
	       }	     

	     else if( hName[ih] == "Top1M" )
	       {
		  hSig[ih]->Fill(Top1M);
		  hSigGen[ih]->Fill(TopGen1M);
	       }	     
	     else if( hName[ih] == "Top2M" )
	       {
		  hSig[ih]->Fill(Top2M);
		  hSigGen[ih]->Fill(TopGen2M);
	       }	     
	     else if( hName[ih] == "Top1P" )
	       {
		  hSig[ih]->Fill(Top1P);
		  hSigGen[ih]->Fill(TopGen1P);
	       }	     
	     else if( hName[ih] == "Top2P" )
	       {
		  hSig[ih]->Fill(Top2P);
		  hSigGen[ih]->Fill(TopGen2P);
	       }	     
	     else if( hName[ih] == "Top1Eta" )
	       {
		  hSig[ih]->Fill(Top1Eta);
		  hSigGen[ih]->Fill(TopGen1Eta);
	       }	     
	     else if( hName[ih] == "Top2Eta" )
	       {
		  hSig[ih]->Fill(Top2Eta);
		  hSigGen[ih]->Fill(TopGen2Eta);
	       }	     
	  }	
     }

   // Plots
   
   TCanvas *c1 = new TCanvas();
   c1->Draw();
   c1->cd();

   TPad *c1_1;
   
   gStyle->SetHistTopMargin(0);

   for(int i=0;i<nHist;i++)
     {		
	addbin(hSig[i]);
	addbin(hSigGen[i]);

	hSig[i]->Scale(1./hSig[i]->Integral());
	hSigGen[i]->Scale(1./hSigGen[i]->Integral());
   
	hSig[i]->SetLineWidth(2);	
	hSig[i]->SetLineColor(kRed);
	hSig[i]->SetMarkerColor(kRed);
	hSig[i]->SetMarkerStyle(20);

	hSigGen[i]->SetLineWidth(2);	
	hSigGen[i]->SetLineColor(kBlue);
	hSigGen[i]->SetMarkerColor(kBlue);
	hSigGen[i]->SetMarkerStyle(22);
	
	hSig[i]->Draw("hist e1");
	hSigGen[i]->Draw("hist e1 same");
	
	hSig[i]->GetXaxis()->SetTitle(hLab[i].c_str());
	hSig[i]->GetYaxis()->SetTitle("Normalized to unity");
   
	float maxSig = hSig[i]->GetMaximum();
	float maxSigGen = hSigGen[i]->GetMaximum();
	float max = (maxSig > maxSigGen) ? maxSig : maxSigGen;
	
	hSig[i]->SetMaximum(1.2*max);      

        TLegend *leg = new TLegend(0.70,0.90,0.90,0.70);
	leg->SetFillColor(253);
	leg->SetBorderSize(0);
	
	leg->AddEntry(hSig[i],"KinFit","lp");
	leg->AddEntry(hSigGen[i],"MC Truth","lp");
	
	leg->Draw();
	
	std::string figName = "pics/"+hName[i]+"Val.eps";
	c1->Print(figName.c_str());
	c1->Clear();
     }   
   
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
