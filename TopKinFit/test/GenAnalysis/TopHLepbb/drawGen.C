#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);

void setTDRStyle();

void drawGen()
{
   gROOT->SetBatch();

   setTDRStyle();
   
   TChain ch("trGEN");
   ch.Add("runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_hut-MadGraph5-pythia8/data.root");
   ch.Add("runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_hct-MadGraph5-pythia8/data.root");

   float TopLepPt, TopLepEta, TopLepPhi, TopLepE;
   float HiggsPt, HiggsEta, HiggsPhi, HiggsE;
   float TopLepWLepPt, TopLepWLepEta, TopLepWLepPhi, TopLepWLepE;
   float TopLepBJetPt, TopLepBJetEta, TopLepBJetPhi, TopLepBJetE;
   float HiggsBJet1Pt, HiggsBJet1Eta, HiggsBJet1Phi, HiggsBJet1E;
   float HiggsBJet2Pt, HiggsBJet2Eta, HiggsBJet2Phi, HiggsBJet2E;

   ch.SetBranchAddress("TopLepPt",&TopLepPt);
   ch.SetBranchAddress("TopLepEta",&TopLepEta);
   ch.SetBranchAddress("TopLepPhi",&TopLepPhi);
   ch.SetBranchAddress("TopLepE",&TopLepE);
   ch.SetBranchAddress("HiggsPt",&HiggsPt);
   ch.SetBranchAddress("HiggsEta",&HiggsEta);
   ch.SetBranchAddress("HiggsPhi",&HiggsPhi);
   ch.SetBranchAddress("HiggsE",&HiggsE);
   ch.SetBranchAddress("TopLepWLepPt",&TopLepWLepPt);
   ch.SetBranchAddress("TopLepWLepEta",&TopLepWLepEta);
   ch.SetBranchAddress("TopLepWLepPhi",&TopLepWLepPhi);
   ch.SetBranchAddress("TopLepWLepE",&TopLepWLepE);
   ch.SetBranchAddress("TopLepBJetPt",&TopLepBJetPt);
   ch.SetBranchAddress("TopLepBJetEta",&TopLepBJetEta);
   ch.SetBranchAddress("TopLepBJetPhi",&TopLepBJetPhi);
   ch.SetBranchAddress("TopLepBJetE",&TopLepBJetE);
   ch.SetBranchAddress("HiggsBJet1Pt",&HiggsBJet1Pt);
   ch.SetBranchAddress("HiggsBJet1Eta",&HiggsBJet1Eta);
   ch.SetBranchAddress("HiggsBJet1Phi",&HiggsBJet1Phi);
   ch.SetBranchAddress("HiggsBJet1E",&HiggsBJet1E);
   ch.SetBranchAddress("HiggsBJet2Pt",&HiggsBJet2Pt);
   ch.SetBranchAddress("HiggsBJet2Eta",&HiggsBJet2Eta);
   ch.SetBranchAddress("HiggsBJet2Phi",&HiggsBJet2Phi);
   ch.SetBranchAddress("HiggsBJet2E",&HiggsBJet2E);

   int nHist = 0;
   
   TH1D *h[1000];
   std::string hName[1000];
   std::string hLab[1000];

   h[nHist] = new TH1D("h_GenTopLepPt","h_GenTopLepPt",50,0.,500.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepPt";
   hLab[nHist] = "p_{T} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_GenTopLepEta","h_GenTopLepEta",50,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepEta";
   hLab[nHist] = "#eta";
   nHist++;

   h[nHist] = new TH1D("h_GenHiggsPt","h_GenHiggsPt",50,0.,500.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenHiggsPt";
   hLab[nHist] = "p_{T} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_GenHiggsEta","h_GenHiggsEta",50,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenHiggsEta";
   hLab[nHist] = "#eta";
   nHist++;

   h[nHist] = new TH1D("h_GenTopLepWLepPt","h_GenTopLepWLepPt",50,0.,300.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepWLepPt";
   hLab[nHist] = "p_{T} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_GenTopLepWLepEta","h_GenTopLepWLepEta",50,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepWLepEta";
   hLab[nHist] = "#eta";
   nHist++;

   h[nHist] = new TH1D("h_GenTopLepBJetPt","h_GenTopLepBJetPt",50,0.,300.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepBJetPt";
   hLab[nHist] = "p_{T} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_GenTopLepBJetEta","h_GenTopLepBJetEta",50,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenTopLepBJetEta";
   hLab[nHist] = "#eta";
   nHist++;

   h[nHist] = new TH1D("h_GenHiggsBJetPt","h_GenHiggsBJetPt",50,0.,300.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenHiggsBJetPt";
   hLab[nHist] = "p_{T} [GeV]";
   nHist++;
   
   h[nHist] = new TH1D("h_GenHiggsBJetEta","h_GenHiggsBJetEta",50,-5.,5.);
   h[nHist]->Sumw2();
   hName[nHist] = "GenHiggsBJetEta";
   hLab[nHist] = "#eta";
   nHist++;
   
   int nev = ch.GetEntries();   
   
   for(int i=0;i<nev;i++)
     {
	ch.GetEntry(i);
	
	TLorentzVector TopLep; TopLep.SetPtEtaPhiE(TopLepPt,TopLepEta,TopLepPhi,TopLepE);
	TLorentzVector Higgs; Higgs.SetPtEtaPhiE(HiggsPt,HiggsEta,HiggsPhi,HiggsE);
	TLorentzVector TopLepWLep; TopLepWLep.SetPtEtaPhiE(TopLepWLepPt,TopLepWLepEta,TopLepWLepPhi,TopLepWLepE);
	TLorentzVector TopLepBJet; TopLepBJet.SetPtEtaPhiE(TopLepBJetPt,TopLepBJetEta,TopLepBJetPhi,TopLepBJetE);
	TLorentzVector HiggsBJet1; HiggsBJet1.SetPtEtaPhiE(HiggsBJet1Pt,HiggsBJet1Eta,HiggsBJet1Phi,HiggsBJet1E);
	TLorentzVector HiggsBJet2; HiggsBJet2.SetPtEtaPhiE(HiggsBJet2Pt,HiggsBJet2Eta,HiggsBJet2Phi,HiggsBJet2E);
	
	for(int ih=0;ih<nHist;ih++)
	  {
	     if( hName[ih] == "GenTopLepPt" ) h[ih]->Fill(TopLepPt);
	     else if( hName[ih] == "GenTopLepEta" ) h[ih]->Fill(TopLepEta);
	     else if( hName[ih] == "GenHiggsPt" ) h[ih]->Fill(HiggsPt);
	     else if( hName[ih] == "GenHiggsEta" ) h[ih]->Fill(HiggsEta);
	     else if( hName[ih] == "GenTopLepWLepPt" ) h[ih]->Fill(TopLepWLepPt);
	     else if( hName[ih] == "GenTopLepWLepEta" ) h[ih]->Fill(TopLepWLepEta);
	     else if( hName[ih] == "GenTopLepBJetPt" ) h[ih]->Fill(TopLepBJetPt);
	     else if( hName[ih] == "GenTopLepBJetEta" ) h[ih]->Fill(TopLepBJetEta);
	     else if( hName[ih] == "GenHiggsBJetPt" ) {h[ih]->Fill(HiggsBJet1Pt);h[ih]->Fill(HiggsBJet2Pt);}
	     else if( hName[ih] == "GenHiggsBJetEta" ) {h[ih]->Fill(HiggsBJet1Eta);h[ih]->Fill(HiggsBJet2Eta);}
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
	addbin(h[i]);

	h[i]->Scale(1./h[i]->Integral());
	h[i]->Scale(1./h[i]->GetMaximum());
   
	h[i]->SetLineWidth(2);	
	h[i]->SetLineColor(kRed);
	h[i]->SetMarkerColor(kRed);
	h[i]->SetMarkerStyle(20);
	h[i]->Draw("hist e1");
	h[i]->GetXaxis()->SetTitle(hLab[i].c_str());
	h[i]->GetYaxis()->SetTitle("Normalized to unity");
	
	float max = h[i]->GetMaximum();
	
	h[i]->SetMaximum(1.2*max);
	
	h[i]->Draw("hist e1");
	
	std::string figName = "pics/"+hName[i]+".eps";	
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
