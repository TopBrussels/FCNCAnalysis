#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);
void addbin2D(TH2D *h);

void setTDRStyle();

void draw_valDiff()
{
   gROOT->SetBatch();

   setTDRStyle();
   
   TFile *fSig = TFile::Open("runTEST_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8/data.root");

   double tim;
   
   float P0Disc;
   
   float NuPx, NuPy, NuPz;

   float NuGenPx, NuGenPy, NuGenPz;

   float TopP, TopPt, TopEta, TopM, TopRap;
   float WP, WPt, WEta, WPhi, WRap;

   float TopGenP, TopGenPt, TopGenEta, TopGenM, TopGenRap;
   float WGenP, WGenPt, WGenEta, WGenPhi, WGenRap;
   
   float WMass;
   float WMassGen;
   
   TTree *trSig = (TTree*)fSig->Get("trKFIT");

   trSig->SetBranchAddress("tim",&tim);
   
   trSig->SetBranchAddress("P0Disc",&P0Disc);
   
   trSig->SetBranchAddress("NuPx",&NuPx);
   trSig->SetBranchAddress("NuPy",&NuPy);
   trSig->SetBranchAddress("NuPz",&NuPz);

   trSig->SetBranchAddress("NuGenPx",&NuGenPx);
   trSig->SetBranchAddress("NuGenPy",&NuGenPy);
   trSig->SetBranchAddress("NuGenPz",&NuGenPz);

   trSig->SetBranchAddress("TopP",&TopP);
   trSig->SetBranchAddress("TopPt",&TopPt);
   trSig->SetBranchAddress("TopEta",&TopEta);
   trSig->SetBranchAddress("TopM",&TopM);
   trSig->SetBranchAddress("TopRap",&TopRap);

   trSig->SetBranchAddress("TopGenP",&TopGenP);
   trSig->SetBranchAddress("TopGenPt",&TopGenPt);
   trSig->SetBranchAddress("TopGenEta",&TopGenEta);
   trSig->SetBranchAddress("TopGenM",&TopGenM);
   trSig->SetBranchAddress("TopGenRap",&TopGenRap);

   trSig->SetBranchAddress("WP",&WP);
   trSig->SetBranchAddress("WPt",&WPt);
   trSig->SetBranchAddress("WEta",&WEta);
   trSig->SetBranchAddress("WPhi",&WPhi);
   trSig->SetBranchAddress("WRap",&WRap);

   trSig->SetBranchAddress("WGenP",&WGenP);
   trSig->SetBranchAddress("WGenPt",&WGenPt);
   trSig->SetBranchAddress("WGenEta",&WGenEta);
   trSig->SetBranchAddress("WGenPhi",&WGenPhi);
   trSig->SetBranchAddress("WGenRap",&WGenRap);

   trSig->SetBranchAddress("WMass",&WMass);
   
   trSig->SetBranchAddress("WMassGen",&WMassGen);
   
   int nHist = 0;
   int nHistArr = 0;
   int nHist2D = 0;
   
   TH1D *hSig[1000];
   std::string hName[1000];
   std::string hLab[1000];

   TH2D *hSig2D[1000];
   std::string hName2D[1000];
   std::string hLab2DX[1000];
   std::string hLab2DY[1000];
   
   TH1D *hSigArr[1000][100];
   std::string hNameArr[1000];
   std::string hLabX[1000];
   std::string hLabY[1000];

   hSig[nHist] = new TH1D("hSig_tim","hSig_tim",60,0.,5.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "tim";
   hLab[nHist] = "Computation time [s]";
   nHist++;
   
   // Nu1
/*   hSig[nHist] = new TH1D("hSig_PxNu1","hSig_PxNu1",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu1Px";
   hLab[nHist] = "p_{x} (#nu_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PyNu1","hSig_PyNu1",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu1Py";
   hLab[nHist] = "p_{y} (#nu_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PzNu1","hSig_PzNu1",30,-500.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu1Pz";
   hLab[nHist] = "p_{z} (#nu_{1}) [GeV]";
   nHist++;

   // Nu2
   hSig[nHist] = new TH1D("hSig_PxNu2","hSig_PxNu2",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu2Px";
   hLab[nHist] = "p_{x} (#nu_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PyNu2","hSig_PyNu2",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu2Py";
   hLab[nHist] = "p_{y} (#nu_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PzNu2","hSig_PzNu2",30,-500.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Nu2Pz";
   hLab[nHist] = "p_{z} (#nu_{2}) [GeV]";
   nHist++;*/

   // W
   
/*   hSig[nHist] = new TH1D("hSig_WMass1","hSig_WMass1",30,60.,110.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMass1";
   hLab[nHist] = "m(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WMass2","hSig_WMass2",30,60.,110.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMass2";
   hLab[nHist] = "m(W_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W1P","hSig_W1P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "W1P";
   hLab[nHist] = "p(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W2P","hSig_W2P",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "W2P";
   hLab[nHist] = "p(W_{2}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W1Pt","hSig_W1Pt",30,0.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "W1Pt";
   hLab[nHist] = "p_{T}(W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_W2Pt","hSig_W2Pt",30,0.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "W2Pt";
   hLab[nHist] = "p_{T}(W_{2}) [GeV]";
   nHist++;*/

   // Top
   
   hSig[nHist] = new TH1D("hSig_TopP","hSig_TopP",60,-500.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "TopP";
   hLab[nHist] = "p(t)^{gen}-p(t)^{rec} [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_TopPt","hSig_TopPt",100,-500.,500.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "TopPt";
   hLab[nHist] = "p_{T}(t)^{gen}-p_{T}(t)^{rec} [GeV]";
   nHist++;

   for(int ib=0;ib<10;ib++)
     {	
	std::string hname = "hSigArr_TopPt_"+std::string(Form("%d",ib));
	hSigArr[nHistArr][ib] = new TH1D(hname.c_str(),hname.c_str(),100,-500.,500.);
	hSigArr[nHistArr][ib]->Sumw2();
     }   
   hNameArr[nHistArr] = "TopPtArr";
   hLabX[nHistArr] = "p_{T}(t)^{gen} [GeV]";
   hLabY[nHistArr] = "p_{T}(t)^{gen}-p_{T}(t)^{rec} [GeV]";
   nHistArr++;
   
   hSig[nHist] = new TH1D("hSig_TopEta","hSig_TopEta",60,-3.,3.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "TopEta";
   hLab[nHist] = "#eta(t)^{gen}-#eta(t)^{rec}";
   nHist++;

   hSig[nHist] = new TH1D("hSig_TopRap","hSig_TopRap",100,-3.,3.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "TopRap";
   hLab[nHist] = "y(t)^{gen}-y(t)^{rec}";
   nHist++;

   hSig[nHist] = new TH1D("hSig_TopM","hSig_TopM",100,-200.,200.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "TopM";
   hLab[nHist] = "m(t)^{gen}-m(t)^{rec}";
   nHist++;
   
   for(int ib=0;ib<10;ib++)
     {	
	std::string hname = "hSigArr_TopRap_"+std::string(Form("%d",ib));
	hSigArr[nHistArr][ib] = new TH1D(hname.c_str(),hname.c_str(),100,-3.,3.);
	hSigArr[nHistArr][ib]->Sumw2();
     }   
   hNameArr[nHistArr] = "TopRapArr";
   hLabX[nHistArr] = "y(t)^{gen}";
   hLabY[nHistArr] = "y(t)^{gen}-y(t)^{rec}";
   nHistArr++;

   hSig2D[nHist2D] = new TH2D("hSig2D_TopPt","hSig2D_TopPt",20,0.,500.,20,-500.,500.);
   hSig2D[nHist2D]->Sumw2();
   hName2D[nHist2D] = "TopPt2D";
   hLab2DX[nHist2D] = "p_{T}(t)^{gen} [GeV]";
   hLab2DY[nHist2D] = "p_{T}(t)^{gen}-p_{T}(t)^{rec} [GeV]";
   nHist2D++;

   hSig2D[nHist2D] = new TH2D("hSig2D_TopRap","hSig2D_TopRap",20,-3.,3.,20,-3.,3.);
   hSig2D[nHist2D]->Sumw2();
   hName2D[nHist2D] = "TopRap2D";
   hLab2DX[nHist2D] = "y(t)^{gen}";
   hLab2DY[nHist2D] = "y(t)^{gen}-y(t)^{rec}";
   nHist2D++;
   
   int nevSig = trSig->GetEntries();   
   
   float nSol = 0;
   
   for(int i=0;i<nevSig;i++)
     {
	trSig->GetEntry(i);

	if( P0Disc > 10E+9 ) continue;
//	if( P0Disc > 0.5 ) continue;
	
	nSol++;
	
	for(int ih=0;ih<nHist;ih++)
	  {	     
	     if( hName[ih] == "tim" )
	       {
		  hSig[ih]->Fill(tim);
	       }	     

	     else if( hName[ih] == "TopP" )
	       {
		  hSig[ih]->Fill(TopGenP-TopP);
	       }	     
	     else if( hName[ih] == "TopPt" )
	       {
		  hSig[ih]->Fill(TopGenPt-TopPt);
	       }	     
	     else if( hName[ih] == "TopEta" )
	       {
		  hSig[ih]->Fill(TopGenEta-TopEta);
	       }	     
	     else if( hName[ih] == "TopRap" )
	       {
		  hSig[ih]->Fill(TopGenRap-TopRap);
	       }	     
	     else if( hName[ih] == "TopM" )
	       {
		  hSig[ih]->Fill(TopGenM-TopM);
	       }	     
	  }
	
	for(int ih=0;ih<nHistArr;ih++)
	  {	     
	     if( hNameArr[ih] == "TopPtArr" )
	       {
		  int iBin = -1;
		  if( TopGenPt < 30. ) continue;
		  else if( TopGenPt > 30. && TopGenPt < 120. ) iBin = 0;
		  else if( TopGenPt > 120. && TopGenPt < 250. ) iBin = 1;
		  else if( TopGenPt > 250. && TopGenPt < 350. ) iBin = 2;
		  else if( TopGenPt > 350. && TopGenPt < 500. ) iBin = 3;
		  else continue;
		  
		  hSigArr[ih][iBin]->Fill(TopGenPt-TopPt);
	       }	     
	     else if( hNameArr[ih] == "TopRapArr" )
	       {
		  int iBin = -1;
		  if( fabs(TopGenRap) > 2.5 ) continue;
		  else if( TopGenRap > -2.5 && TopGenRap < -1.5 ) iBin = 0;
		  else if( TopGenRap > -1.5 && TopGenRap < -1.0 ) iBin = 1;
		  else if( TopGenRap > -1.0 && TopGenRap < -0.5 ) iBin = 2;
		  else if( TopGenRap > -0.5 && TopGenRap < 0 ) iBin = 3;
		  else if( TopGenRap > 0 && TopGenRap < 0.5 ) iBin = 4;
		  else if( TopGenRap > 0.5 && TopGenRap < 1.0 ) iBin = 5;
		  else if( TopGenRap > 1.0 && TopGenRap < 1.5 ) iBin = 6;
		  else if( TopGenRap > 1.5 && TopGenRap < 2.5 ) iBin = 7;
		  else continue;
		  
		  hSigArr[ih][iBin]->Fill(TopGenRap-TopRap);
	       }	     
	  }
	
	for(int ih=0;ih<nHist2D;ih++)
	  {	     
	     if( hName2D[ih] == "TopPt2D" )
	       {
		  hSig2D[ih]->Fill(TopPt,TopGenPt-TopPt);
	       }	     
	     else if( hName2D[ih] == "TopRap2D" )
	       {
		  hSig2D[ih]->Fill(TopRap,TopGenRap-TopRap);
	       }	     
	  }
	
     }

   std::cout << nSol/float(nevSig)*100. << "% of events were resolved" << std::endl;
   
   // Plots
   
   TCanvas *c1 = new TCanvas();
   c1->Draw();
   c1->cd();

   TPad *c1_1;
   
   gStyle->SetHistTopMargin(0);

   for(int i=0;i<nHist;i++)
     {		
	addbin(hSig[i]);

//	hSig[i]->Scale(1./hSig[i]->Integral());

	if( hName[i] == "tim" ) std::cout << "Time/event = " << hSig[i]->GetMean() << " +- " << hSig[i]->GetRMS() << " s " << std::endl;
	
	hSig[i]->SetLineWidth(2);	
	hSig[i]->SetLineColor(kRed);
//	hSig[i]->SetMarkerColor(kRed);
//	hSig[i]->SetMarkerStyle(20);
	hSig[i]->SetMarkerSize(0.);

	hSig[i]->Draw("hist e1");
	
	hSig[i]->GetXaxis()->SetTitle(hLab[i].c_str());
	hSig[i]->GetYaxis()->SetTitle("Number of events");
   
	float maxSig = hSig[i]->GetMaximum();
	
	hSig[i]->SetMaximum(1.2*maxSig);

        TLegend *leg = new TLegend(0.70,0.90,0.90,0.70);
	leg->SetFillColor(253);
	leg->SetBorderSize(0);
	
	leg->AddEntry(hSig[i],"KinFit","lp");
	
//	leg->Draw();
	
	std::string figName = "pics/"+hName[i]+"Val.eps";
	c1->Print(figName.c_str());
	c1->Clear();
     }   

   for(int i=0;i<nHistArr;i++)
     {
	float x[100];
	float y[100];
	float ex[100];
	float ey[100];
	TGraphErrors *hMean;
	TGraphErrors *hRMS;
	int nBins = -1;
	
	// Mean
	
	if( hNameArr[i] == "TopPtArr" )
	  {	
	     nBins = 4;

	     x[0] = 75;
	     x[1] = 185;
	     x[2] = 300;
	     x[3] = 425;

	     ex[0] = 45;
	     ex[1] = 65;
	     ex[2] = 50;
	     ex[3] = 75;
	  }
	else if( hNameArr[i] == "TopRapArr" )
	  {
	     nBins = 8;
	     
	     x[0] = -2.0;
	     x[1] = -1.25;
	     x[2] = -0.75;
	     x[3] = -0.25;
	     x[4] = 0.25;
	     x[5] = 0.75;
	     x[6] = 1.25;
	     x[7] = 2.0;

	     ex[0] = 0.5;
	     ex[1] = 0.25;
	     ex[2] = 0.25;
	     ex[3] = 0.25;
	     ex[4] = 0.25;
	     ex[5] = 0.25;
	     ex[6] = 0.25;
	     ex[7] = 0.5;
	  }
	
	for(int ib=0;ib<nBins;ib++)
	  {	     
	     addbin(hSigArr[i][ib]);
	     y[ib] = hSigArr[i][ib]->GetMean();
	     ey[ib] = hSigArr[i][ib]->GetMeanError();
	  }	
	
	hMean = new TGraphErrors(nBins,x,y,ex,ey);
	
	hMean->SetLineWidth(2);	
	hMean->SetLineColor(kRed);
	hMean->SetMarkerColor(kRed);
	hMean->SetMarkerStyle(20);

	hMean->Draw("AP");
	
	std::string hLabYMod = hLabY[i] + " (Mean)";
	
	hMean->GetXaxis()->SetTitle(hLabX[i].c_str());
	hMean->GetYaxis()->SetTitle(hLabYMod.c_str());

	TLine *line = new TLine(hMean->GetXaxis()->GetBinLowEdge(1),0.,
				hMean->GetXaxis()->GetBinUpEdge(hMean->GetXaxis()->GetNbins()),0.);
	
	line->SetLineColor(kBlack);
	line->SetLineStyle(2);
	line->SetLineWidth(2);
	line->Draw();
	
	std::string figName = "pics/"+hNameArr[i]+"MeanVal.eps";
	c1->Print(figName.c_str());
	c1->Clear();
	delete hMean;
	
	// RMS
	
	for(int ib=0;ib<nBins;ib++)
	  {	     
	     y[ib] = hSigArr[i][ib]->GetRMS();
	     ey[ib] = hSigArr[i][ib]->GetRMSError();
	  }	
	
	hRMS = new TGraphErrors(nBins,x,y,ex,ey);
	
	hRMS->SetLineWidth(2);	
	hRMS->SetLineColor(kRed);
	hRMS->SetMarkerColor(kRed);
	hRMS->SetMarkerStyle(20);

	hRMS->Draw("AP");
	
	hLabYMod = hLabY[i] + " (RMS)";
	
	hRMS->GetXaxis()->SetTitle(hLabX[i].c_str());
	hRMS->GetYaxis()->SetTitle(hLabYMod.c_str());
	
	figName = "pics/"+hNameArr[i]+"RMSVal.eps";
	c1->Print(figName.c_str());
	c1->Clear();
	delete hRMS;
     }   

   for(int i=0;i<nHist2D;i++)
     {		
	addbin2D(hSig2D[i]);
   
	hSig2D[i]->SetLineWidth(2);	
	hSig2D[i]->SetLineColor(kRed);

	hSig2D[i]->Draw("BOX");
	
	hSig2D[i]->GetXaxis()->SetTitle(hLab2DX[i].c_str());
	hSig2D[i]->GetYaxis()->SetTitle(hLab2DY[i].c_str());
	
	std::string figName = "pics/"+hName2D[i]+"Val.eps";
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

void addbin2D(TH2D *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   Int_t y_nbins = h->GetYaxis()->GetNbins();
   
   for(int iy=1;iy<y_nbins+1;iy++)
     {	
	h->SetBinContent(1,iy,h->GetBinContent(0,iy)+h->GetBinContent(1,iy));
	h->SetBinError(1,iy,TMath::Sqrt(pow(h->GetBinError(0,iy),2)+pow(h->GetBinError(1,iy),2)));
	h->SetBinContent(x_nbins,iy,h->GetBinContent(x_nbins,iy)+h->GetBinContent(x_nbins+1,iy));
	h->SetBinError(x_nbins,iy,TMath::Sqrt(pow(h->GetBinError(x_nbins,iy),2)+
					      pow(h->GetBinError(x_nbins+1,iy),2)));
	// Set overflow and underflow bins to 0
	h->SetBinContent(0,iy,0.);
	h->SetBinError(0,iy,0.);
	h->SetBinContent(x_nbins+1,iy,0.);
	h->SetBinError(x_nbins+1,iy,0.);
     }   

   for(int ix=1;ix<x_nbins+1;ix++)
     {	
	h->SetBinContent(ix,1,h->GetBinContent(ix,0)+h->GetBinContent(ix,1));
	h->SetBinError(ix,1,TMath::Sqrt(pow(h->GetBinError(ix,0),2)+pow(h->GetBinError(ix,1),2)));
	h->SetBinContent(ix,y_nbins,h->GetBinContent(ix,y_nbins)+h->GetBinContent(ix,y_nbins+1));
	h->SetBinError(ix,y_nbins,TMath::Sqrt(pow(h->GetBinError(ix,y_nbins),2)+
					      pow(h->GetBinError(ix,y_nbins+1),2)));
	// Set overflow and underflow bins to 0
	h->SetBinContent(ix,0,0.);
	h->SetBinError(ix,0,0.);
	h->SetBinContent(ix,y_nbins+1,0.);
	h->SetBinError(ix,y_nbins+1,0.);
     }   
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
