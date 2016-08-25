#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);

void setTDRStyle();

void draw_test()
{
   gROOT->SetBatch();
   
   setTDRStyle();
   
   TFile *fSig = TFile::Open("runTEST5K_MERGED/TT_TuneCUETP8M1_13TeV-powheg-pythia8/data.root");

   float P0Disc, P1Disc, P2Disc, P3Disc, P4Disc;
   float P0Term0, P0Term1;
   
   bool P0MatchBJet1, P0MatchBJet2;
   bool P0RecoBJet1, P0RecoBJet2;

   bool P0MatchTopWL1, P0MatchTopWL2;
   bool P0RecoTopWL1, P0RecoTopWL2;
   
   float Nu1Px, Nu1Py, Nu1Pz;
   float Nu2Px, Nu2Py, Nu2Pz;
   float WMass1, WMass2;
   float WMassGen1, WMassGen2;
   float Top1M, Top2M;
   float MetXMeas, MetYMeas;
   float MetXFit, MetYFit;

   TTree *trSig = (TTree*)fSig->Get("trKFIT");
   trSig->SetBranchAddress("P0Disc",&P0Disc);
   trSig->SetBranchAddress("P0Term0",&P0Term0);
   trSig->SetBranchAddress("P0Term1",&P0Term1);
   trSig->SetBranchAddress("P0MatchBJet1",&P0MatchBJet1);
   trSig->SetBranchAddress("P0MatchBJet2",&P0MatchBJet2);
   trSig->SetBranchAddress("P0RecoBJet1",&P0RecoBJet1);
   trSig->SetBranchAddress("P0RecoBJet2",&P0RecoBJet2);
   trSig->SetBranchAddress("P0MatchTopWL1",&P0MatchTopWL1);
   trSig->SetBranchAddress("P0MatchTopWL2",&P0MatchTopWL2);
   trSig->SetBranchAddress("P0RecoTopWL1",&P0RecoTopWL1);
   trSig->SetBranchAddress("P0RecoTopWL2",&P0RecoTopWL2);
   trSig->SetBranchAddress("Nu1Px",&Nu1Px);
   trSig->SetBranchAddress("Nu1Py",&Nu1Py);
   trSig->SetBranchAddress("Nu1Pz",&Nu1Pz);
   trSig->SetBranchAddress("Nu2Px",&Nu2Px);
   trSig->SetBranchAddress("Nu2Py",&Nu2Py);
   trSig->SetBranchAddress("Nu2Pz",&Nu2Pz);
   trSig->SetBranchAddress("WMass1",&WMass1);
   trSig->SetBranchAddress("WMass2",&WMass2);
   trSig->SetBranchAddress("WMassGen1",&WMassGen1);
   trSig->SetBranchAddress("WMassGen2",&WMassGen2);
   trSig->SetBranchAddress("Top1M",&Top1M);
   trSig->SetBranchAddress("Top2M",&Top2M);
   trSig->SetBranchAddress("MetXMeas",&MetXMeas);
   trSig->SetBranchAddress("MetYMeas",&MetYMeas);
   trSig->SetBranchAddress("MetXFit",&MetXFit);
   trSig->SetBranchAddress("MetYFit",&MetYFit);

   int nHist = 0;
   
   TH1D *hSig[1000];
   std::string hName[1000];
   std::string hLab[1000];
   
   hSig[nHist] = new TH1D("hSig_P0Disc","hSig_P0Disc",100,0.,30.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "P0Disc";
   hLab[nHist] = "D";
   nHist++;

   // Term
   hSig[nHist] = new TH1D("hSig_P0Term0","hSig_P0Term0",25,0.,1.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "P0Term0";
   hLab[nHist] = "#chi^{2} (Wmass1)";
   nHist++;

   hSig[nHist] = new TH1D("hSig_P0Term1","hSig_P0Term1",25,0.,1.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "P0Term1";
   hLab[nHist] = "#chi^{2} (Wmass1)";
   nHist++;
   
   // Nu1
   hSig[nHist] = new TH1D("hSig_PxNu1","hSig_PxNu1",30,-300.,300.);
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
   nHist++;

   // WMass
   hSig[nHist] = new TH1D("hSig_WMass1","hSig_WMass1",30,75.,90.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMass1";
   hLab[nHist] = "mass (W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WMass2","hSig_WMass2",30,75.,90.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMass2";
   hLab[nHist] = "mass (W_{2}) [GeV]";
   nHist++;

   // WMassGen
   hSig[nHist] = new TH1D("hSig_WMassGen1","hSig_WMassGen1",30,75.,90.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMassGen1";
   hLab[nHist] = "generated mass (W_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_WMassGen2","hSig_WMassGen2",30,75.,90.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "WMassGen2";
   hLab[nHist] = "generated mass (W_{2}) [GeV]";
   nHist++;

   // TopMass
   hSig[nHist] = new TH1D("hSig_Top1M","hSig_Top1M",30,120.,240.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Top1M";
   hLab[nHist] = "mass (top_{1}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_Top2M","hSig_Top2M",30,120.,240.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "Top2M";
   hLab[nHist] = "mass (top_{2}) [GeV]";
   nHist++;

   // MET
   hSig[nHist] = new TH1D("hSig_dMETx","hSig_dMETx",30,-5.,5.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "dMETx";
   hLab[nHist] = "#Delta (EMiss^{measured}_{x}-EMiss^{fit}_{x}) [GeV]";
   nHist++;

   hSig[nHist] = new TH1D("hSig_dMETy","hSig_dMETy",30,-5.,5.);
   hSig[nHist]->Sumw2();
   hName[nHist] = "dMETy";
   hLab[nHist] = "#Delta (EMiss^{measured}_{y}-EMiss^{fit}_{y}) [GeV]";
   nHist++;
   
   int nevSig = trSig->GetEntries();   
   
   for(int i=0;i<nevSig;i++)
     {
	trSig->GetEntry(i);

	for(int ih=0;ih<nHist;ih++)
	  {	     
	     if( hName[ih] == "P0Disc" ) hSig[ih]->Fill(P0Disc);
	     else if( hName[ih] == "P0Term0" ) hSig[ih]->Fill(P0Term0);
	     else if( hName[ih] == "P0Term1" ) hSig[ih]->Fill(P0Term1);
	     else if( hName[ih] == "Nu1Px" ) hSig[ih]->Fill(Nu1Px);
	     else if( hName[ih] == "Nu1Py" ) hSig[ih]->Fill(Nu1Py);
	     else if( hName[ih] == "Nu1Pz" ) hSig[ih]->Fill(Nu1Pz);
	     else if( hName[ih] == "Nu2Px" ) hSig[ih]->Fill(Nu2Px);
	     else if( hName[ih] == "Nu2Py" ) hSig[ih]->Fill(Nu2Py);
	     else if( hName[ih] == "Nu2Pz" ) hSig[ih]->Fill(Nu2Pz);
	     else if( hName[ih] == "WMass1" ) hSig[ih]->Fill(WMass1);
	     else if( hName[ih] == "WMass2" ) hSig[ih]->Fill(WMass2);
	     else if( hName[ih] == "WMassGen1" ) hSig[ih]->Fill(WMassGen1);
	     else if( hName[ih] == "WMassGen2" ) hSig[ih]->Fill(WMassGen2);
	     else if( hName[ih] == "Top1M" ) hSig[ih]->Fill(Top1M);
	     else if( hName[ih] == "Top2M" ) hSig[ih]->Fill(Top2M);
	     else if( hName[ih] == "dMETx" ) hSig[ih]->Fill(MetXMeas-MetXFit);
	     else if( hName[ih] == "dMETy" ) hSig[ih]->Fill(MetYMeas-MetYFit);
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

	hSig[i]->Scale(1./hSig[i]->Integral());   
   
	hSig[i]->SetLineWidth(2);	
	hSig[i]->SetLineColor(kRed);
	hSig[i]->SetMarkerColor(kRed);
	hSig[i]->SetMarkerStyle(20);
	
	hSig[i]->Draw("hist e1");
	
	hSig[i]->GetXaxis()->SetTitle(hLab[i].c_str());
	hSig[i]->GetYaxis()->SetTitle("Normalized to unity");
   
//	if( hName[i] == "P0Disc" ) c1->SetLogx(1);
	
	float max = hSig[i]->GetMaximum();
	
	hSig[i]->SetMaximum(1.2*max);      
	
	std::string figName = "pics/"+hName[i]+".eps";
	c1->Print(figName.c_str());
	c1->Clear();
	
//	if( hName[i] == "P0Disc" ) c1->SetLogx(0);
     }   

   // Tables

   float nP0MatchBJets = 0.;
   float nP0RecoBJets = 0.;

   float nP0MatchTopWLs = 0.;
   float nP0RecoTopWLs = 0.;
   
   for(int i=0;i<nevSig;i++)
     {
	trSig->GetEntry(i);
	
	if( P0RecoBJet1 && P0RecoBJet2 && P0MatchBJet1 && P0MatchBJet2 ) nP0MatchBJets++;
	if( P0RecoBJet1 && P0RecoBJet2 ) nP0RecoBJets++;

	if( P0RecoTopWL1 && P0RecoTopWL2 && P0MatchTopWL1 && P0MatchTopWL2 ) nP0MatchTopWLs++;
	if( P0RecoTopWL1 && P0RecoTopWL2 ) nP0RecoTopWLs++;
     }   
   
   float effP0MatchBJets = nP0MatchBJets/nP0RecoBJets;

   float effP0MatchTopWLs = nP0MatchTopWLs/nP0RecoTopWLs;

   std::string headerDoc = "\\documentclass[a4paper]{article} \n \
\\usepackage[english]{babel} \n \
\\usepackage{graphicx} \n \
\\usepackage{color} \n \
\\begin{document}";
   
   std::string headerTab = "\
\\begin{table}[hbt!] \n \
\\small \n \
\\begin{center} \n \
\\resizebox{16cm}{!}{ \n \
\\begin{tabular}{|c|c|c|c|c|} \n \
\\hline";
   
   std::string row1 = "Permutation id & 1 \\\\";
   std::string row2 = "b jets \\epsilon, \\% & ";
   std::string row3 = "b jets + top leptons \\epsilon, \\% & ";
   
   stringstream sss("");
   sss.precision(2);
   
   sss.str("");
   sss << effP0MatchBJets;
   std::string effP0MatchBJetsStr = sss.str();

   sss.str("");
   sss << effP0MatchTopWLs;
   std::string effP0MatchTopWLsStr = sss.str();
   
   row2 += effP0MatchBJetsStr + " \\\\ ";

   row3 += effP0MatchTopWLsStr + " \\\\ ";
   
   std::string footer = "\n \
\\hline \
\\end{tabular} \n \
} \n \
\\end{center} \n \
\\end{table} \n \
\\end{document}";
   
   std::ofstream ftable("table/matching.tex");
   ftable << headerDoc << std::endl;
   ftable << headerTab << std::endl;
   ftable << row1 << std::endl;
   ftable << row2 << std::endl;
   ftable << row3 << std::endl;
   ftable << footer << std::endl;
   ftable.close();  
   
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
