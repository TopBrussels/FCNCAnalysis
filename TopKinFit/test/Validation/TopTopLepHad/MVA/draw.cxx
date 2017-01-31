#include "PlotStyle.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"

#include <vector>
#include <iostream>
#include <algorithm>

int main()
{
   SetPlotStyle();

   std::string v[6] = {"Truth","HighestCSVv2","CSVv2L","CSVv2M","CSVv2T","All"};

   for(int iv=0;iv<6;iv++)
     {	
	std::string fname = "TMVAPartReco"+v[iv]+".root";
	
	TFile f(fname.c_str());

	TH1D *h_sig_train = (TH1D*)f.Get("Method_BDT/BDT/MVA_BDT_Train_S");
	TH1D *h_bkg_train = (TH1D*)f.Get("Method_BDT/BDT/MVA_BDT_Train_B");
	
	TH1D *h_sig_test = (TH1D*)f.Get("Method_BDT/BDT/MVA_BDT_S");
	TH1D *h_bkg_test = (TH1D*)f.Get("Method_BDT/BDT/MVA_BDT_B");
	
	TCanvas *c1 = new TCanvas();
   
	h_sig_train->SetLineColor(kRed-4);
	h_sig_train->SetMarkerSize(0.7);
	h_sig_train->SetMarkerStyle(20);
	h_sig_train->SetMarkerColor(kRed-4);
	h_sig_test->SetLineColor(kRed-4);
	h_sig_test->SetMarkerSize(0.0);
	//   h_sig_test->SetFillColor(kRed-4);
	//   h_sig_test->SetFillStyle(3354);
	
	h_bkg_train->SetLineColor(kBlue-6);
	h_bkg_train->SetMarkerSize(0.7);
	h_bkg_train->SetMarkerStyle(22);
	h_bkg_train->SetMarkerColor(kBlue-6);
	h_bkg_test->SetLineColor(kBlue-6);
	h_bkg_test->SetMarkerSize(0.0);
	h_bkg_test->SetFillColor(kBlue-10);
	
	gStyle->SetHatchesSpacing(2.0);
	
	h_bkg_test->Draw("hist e1");
	h_sig_test->Draw("hist e1 same");
	
	h_sig_train->Draw("e1 same");
	h_bkg_train->Draw("e1 same");
	
	double max_sig_test = h_sig_test->GetMaximum();
	double max_bkg_test = h_bkg_test->GetMaximum();
	double max_sig_train = h_sig_train->GetMaximum();
	double max_bkg_train = h_bkg_train->GetMaximum();
	
	double max_test = std::max(max_sig_test,max_bkg_test);
	double max_train = std::max(max_sig_train,max_bkg_train);
	double max = std::max(max_test,max_train);
	
	h_bkg_test->SetMaximum(1.2*max);
	
	h_bkg_test->GetXaxis()->SetTitle("MVA output");
	h_bkg_test->GetYaxis()->SetTitle("(1/N) dN/dx");

	TLegend *leg0 = new TLegend(0.20,0.92,0.42,0.70);
	leg0->SetFillColor(253);
	leg0->SetBorderSize(0);
	leg0->AddEntry(h_sig_train,"Signal (train)","l");
	leg0->AddEntry(h_bkg_train,"Background (train)","l");
	leg0->AddEntry(h_sig_test,"Signal (test)","f");
	leg0->AddEntry(h_bkg_test,"Background (test)","f");
	leg0->Draw();
   
	std::string disc = "pics/disc"+v[iv]+".eps";
	c1->Print(disc.c_str());
	c1->Clear();   
   
	TH1D *h_S = (TH1D*)h_bkg_test->Clone("h_S");
	h_S->Clear();
	float maxSign = -1.;
	int nBins = h_S->GetXaxis()->GetNbins();
	for(int i=1;i<=nBins;i++)
	  {
	     double errSig;
	     double errBkg;
	     float intSig = h_sig_test->IntegralAndError(i,nBins,errSig);
	     float intBkg = h_bkg_test->IntegralAndError(i,nBins,errBkg);
	     float sign = intSig/sqrt(intBkg+intSig);
	     h_S->SetBinContent(i,sign);
	     h_S->SetBinError(i,0.);
	     if( sign > maxSign ) maxSign = sign;
	  }
	
	h_S->GetYaxis()->SetTitle("S/#sqrt{S+B}");
	h_S->Draw("hist e1");
	h_S->SetMaximum(maxSign*1.2);
	h_S->SetMinimum(0.0);
	//   h_S->SetMaximum(3.9);
	//   h_S->SetMinimum(3.8);
	//   h_S->SetMaximum(3.7);
	//   h_S->SetMinimum(3.6);
	std::string sign = "pics/sign"+v[iv]+".eps";
	c1->Print(sign.c_str());
	c1->Clear();

	TH1D *h_CorrelationMatrixS = (TH1D*)f.Get("CorrelationMatrixS");
	h_CorrelationMatrixS->Draw("COLZ");
	std::string corS = "pics/corS"+v[iv]+".eps";
	c1->Print(corS.c_str());
	c1->Clear();
   
	TH1D *h_CorrelationMatrixB = (TH1D*)f.Get("CorrelationMatrixB");
	h_CorrelationMatrixB->Draw("COLZ");
	std::string corB = "pics/corB"+v[iv]+".eps";
	c1->Print(corB.c_str());
	c1->Clear();
   
	f.Close();
     }   
}
