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
   
   float NuPx, NuPy, NuPz;

   float NuGenPx, NuGenPy, NuGenPz;

   float TopP, TopPt, TopEta, TopM;
   float WP, WPt, WEta, WPhi;

   float TopGenP, TopGenPt, TopGenEta, TopGenM;
   float WGenP, WGenPt, WGenEta, WGenPhi;
   
   float WMass;
   float WMassGen;
   
   TTree *trSig = (TTree*)fSig->Get("trKFIT");

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

   trSig->SetBranchAddress("TopGenP",&TopGenP);
   trSig->SetBranchAddress("TopGenPt",&TopGenPt);
   trSig->SetBranchAddress("TopGenEta",&TopGenEta);
   trSig->SetBranchAddress("TopGenM",&TopGenM);

   trSig->SetBranchAddress("WP",&WP);
   trSig->SetBranchAddress("WPt",&WPt);
   trSig->SetBranchAddress("WEta",&WEta);
   trSig->SetBranchAddress("WPhi",&WPhi);

   trSig->SetBranchAddress("WGenP",&WGenP);
   trSig->SetBranchAddress("WGenPt",&WGenPt);
   trSig->SetBranchAddress("WGenEta",&WGenEta);
   trSig->SetBranchAddress("WGenPhi",&WGenPhi);

   trSig->SetBranchAddress("WMass",&WMass);
   
   trSig->SetBranchAddress("WMassGen",&WMassGen);
   
   int nHist = 0;
   
   TH1D *hSig[1000];
   TH1D *hSigGen[1000];
   std::string hName[1000];
   std::string hLab[1000];
   
   // Nu1
   hSig[nHist] = new TH1D("hSig_PxNu","hSig_PxNu",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PxNu","hSigGen_PxNu",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "NuPx";
   hLab[nHist] = "p_{x} (#nu) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PyNu","hSig_PyNu",30,-300.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PyNu","hSigGen_PyNu",30,-300.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "NuPy";
   hLab[nHist] = "p_{y} (#nu) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_PzNu","hSig_PzNu",30,-500.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_PzNu","hSigGen_PzNu",30,-500.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "NuPz";
   hLab[nHist] = "p_{z} (#nu) [GeV]";
   nHist++;

   // W
   
   hSig[nHist] = new TH1D("hSig_WMass","hSig_WMass",30,60.,110.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_WMass","hSigGen_WMass",30,60.,110.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "WMass";
   hLab[nHist] = "m(W) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WP","hSig_WP",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_WP","hSigGen_WP",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "WP";
   hLab[nHist] = "p(W) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WPt","hSig_WPt",30,0.,500.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_WPt","hSigGen_WPt",30,0.,500.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "WPt";
   hLab[nHist] = "p_{T}(W) [GeV]";
   nHist++;

   // Top
   
   hSig[nHist] = new TH1D("hSig_TopM","hSig_TopM",30,50.,300.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_TopM","hSigGen_TopM",30,50.,300.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "TopM";
   hLab[nHist] = "m(t) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_TopP","hSig_TopP",30,0.,1000.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_TopP","hSigGen_TopP",30,0.,1000.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "TopP";
   hLab[nHist] = "p(t) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_TopEta","hSig_TopEta",30,-5.,5.);
   hSig[nHist]->Sumw2();
   hSigGen[nHist] = new TH1D("hSigGen_TopEta","hSigGen_TopEta",30,-5.,5.);
   hSigGen[nHist]->Sumw2();
   hName[nHist] = "TopEta";
   hLab[nHist] = "#eta(t)";
   nHist++;
   
   int nevSig = trSig->GetEntries();   
   
   for(int i=0;i<nevSig;i++)
     {
	trSig->GetEntry(i);

	if( P0Disc > 10E+9 ) continue;
//	if( P0Disc > 1 ) continue;
	
	for(int ih=0;ih<nHist;ih++)
	  {	     
	     if( hName[ih] == "NuPx" )
	       {
		  hSig[ih]->Fill(NuPx);
		  hSigGen[ih]->Fill(NuGenPx);
	       }	     
	     else if( hName[ih] == "NuPy" ) 
	       {
		  hSig[ih]->Fill(NuPy);
		  hSigGen[ih]->Fill(NuGenPy);
	       }	     
	     else if( hName[ih] == "NuPz" ) 
	       {
		  hSig[ih]->Fill(NuPz);
		  hSigGen[ih]->Fill(NuGenPz);
	       }	     

	     else if( hName[ih] == "WMass" )
	       {
		  hSig[ih]->Fill(WMass);
		  hSigGen[ih]->Fill(WMassGen);
	       }	     
	     else if( hName[ih] == "WP" )
	       {
		  hSig[ih]->Fill(WP);
		  hSigGen[ih]->Fill(WGenP);
	       }	     
	     else if( hName[ih] == "WPt" )
	       {
		  hSig[ih]->Fill(WPt);
		  hSigGen[ih]->Fill(WGenPt);
	       }	     

	     else if( hName[ih] == "TopM" )
	       {
		  hSig[ih]->Fill(TopM);
		  hSigGen[ih]->Fill(TopGenM);
	       }	     
	     else if( hName[ih] == "TopP" )
	       {
		  hSig[ih]->Fill(TopP);
		  hSigGen[ih]->Fill(TopGenP);
	       }	     
	     else if( hName[ih] == "TopEta" )
	       {
		  hSig[ih]->Fill(TopEta);
		  hSigGen[ih]->Fill(TopGenEta);
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
