#include "PlotStyle.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <sstream>
#include <ctime>

#include <vector>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>

std::string intToStr (int number);
std::string floatToStr (float number);

using namespace std;
int main(int argc, char *argv[])
{

    if(argc < 4)
    {
        std::cout << "INVALID number of arguments. The necessary arguments are: " << std::endl;
        std::cout << "    int baseline_bjets             = strtol(argv[1], NULL,10);" << std::endl;
        std::cout << "    int baseline_jets                 = strtol(argv[2], NULL,10);" << std::endl;
        std::cout << "    std::string SignalSample            = argv[3];" << std::endl;
        std::cout << "    std::string channel            = argv[3];" << std::endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    std::string SignalSample  = argv[3];//Valid arguments are: SThut, SThct, TThct, TThut
    std::string channel            = argv[4];

    SetPlotStyle();
   
    std::string category = "b"+intToStr(baseline_bjets)+"j"+intToStr(baseline_jets);
    std::string TrainingName = "CombTraining_" + SignalSample + channel + "_" +  category;//Example: Training_SThut_El_b3j3
    std::string outputpics = "pics/"+TrainingName+"/";
    mkdir("pics",0777);
    mkdir(outputpics.c_str(),0777);

   TFile f(("../weights/"+TrainingName+".root").c_str());

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

   TLegend *leg0 = new TLegend(0.82,0.92,0.995,0.40);
   leg0->SetFillColor(253);
   leg0->SetBorderSize(0);
   leg0->AddEntry(h_sig_train,"Signal (train)","l");
   leg0->AddEntry(h_bkg_train,"Background (train)","l");
   leg0->AddEntry(h_sig_test,"Signal (test)","f");
   leg0->AddEntry(h_bkg_test,"Background (test)","f");
   leg0->Draw();
   
   c1->Print((outputpics+"disc.eps").c_str());
   c1->Clear();   
   
   TH1D *h_S = (TH1D*)h_bkg_test->Clone("h_S");
   h_S->Clear();
   float maxSign = -1.;
   float xmaxSign = -1.;
   int nBins = h_S->GetXaxis()->GetNbins();
   for(int i=1;i<=nBins;i++)
     {
	double errSig;
	double errBkg;
	float intSig = h_sig_test->IntegralAndError(i,nBins,errSig);
	float intBkg = h_bkg_test->IntegralAndError(i,nBins,errBkg);
//	float sign = intSig/sqrt(intBkg+intSig);
	float sign = intSig/intBkg;
	if(sign <= 0 || sign != sign) sign = 0;
	if(intBkg == 0 && intSig != 0) sign = 50;
//	cout << "sign S/B = " << sign << endl;
	h_S->SetBinContent(i,sign);
	h_S->SetBinError(i,0.);
	if( sign > maxSign )
	{
	    float binWidth = (h_S->GetXaxis()->GetXmax()-h_S->GetXaxis()->GetXmin())/nBins;
	    maxSign = sign;
	    xmaxSign = h_S->GetXaxis()->GetXmin()+i*binWidth;
	}
     }
   
//   h_S->GetYaxis()->SetTitle("S/#sqrt{S+B}");
   h_S->GetYaxis()->SetTitle("S/B");
   h_S->GetXaxis()->SetTitle(("MVA disc. (max sign:"+floatToStr(xmaxSign)+")").c_str());
   h_S->Draw("hist e1");
   h_S->SetMaximum(maxSign*1.2);
   h_S->SetMinimum(0.0);
//   h_S->SetMaximum(3.9);
//   h_S->SetMinimum(3.8);
//   h_S->SetMaximum(3.7);
//   h_S->SetMinimum(3.6);
//   c1->Print((outputpics+"sign.eps").c_str());
   c1->Print((outputpics+"sign_SoverB.eps").c_str());
   c1->Clear();




   TH1D *h_Roc = (TH1D*)h_bkg_test->Clone("h_Roc");
   h_Roc->Clear();
   float minDist =999.;
   float xminDist = -1.;
//   int nBins = h_Roc->GetXaxis()->GetNbins();
   
   float TotalBckg = h_bkg_test->Integral();
   float TotalSig = h_sig_test->Integral();
   
   for(int i=1;i<=nBins;i++)
  {
	      double errSig;
	      double errBkg;
	      float intSig = h_sig_test->Integral(i,nBins);
	      float intBkg = h_bkg_test->Integral(i,nBins);
	      
	      float SigEff = intSig/TotalSig;
	      float BkgRej = 1-(intBkg/TotalBckg);
	      
	      float dist = std::sqrt((1-SigEff)*(1-SigEff) + (1-BkgRej)*(1-BkgRej));
	      h_Roc->SetBinContent(i,dist);
	      h_Roc->SetBinError(i,0.);
	      if( dist < minDist )
	      {
	          float binWidth = (h_Roc->GetXaxis()->GetXmax()-h_Roc->GetXaxis()->GetXmin())/nBins;
	          minDist = dist;
	          xminDist = h_Roc->GetXaxis()->GetXmin()+i*binWidth;
	      }
     }
   
   h_Roc->GetYaxis()->SetTitle("#sqrt(2) - #sqrt(SigEff^2 + BkgRej^2)");
   h_Roc->GetXaxis()->SetTitle(("MVA disc. (optimal ROC:"+floatToStr(xminDist)+")").c_str());
cout << SignalSample << " " << channel << ": ROC best cut is " << xminDist << endl;
   h_Roc->Draw("hist e1");
   h_Roc->SetMaximum(1.2);
//   h_Roc->SetMinimum(0.);
//   h_S->SetMaximum(3.9);
//   h_S->SetMinimum(3.8);
//   h_S->SetMaximum(3.7);
//   h_S->SetMinimum(3.6);
//   c1->Print((outputpics+"sign.eps").c_str());
   c1->Print((outputpics+"ROC_dist.eps").c_str());
   c1->Clear();







   TH1D *h_CorrelationMatrixS = (TH1D*)f.Get("CorrelationMatrixS");
   h_CorrelationMatrixS->Draw("COLZ");
   c1->Print((outputpics+"corS.eps").c_str());
   c1->Clear();
   
   TH1D *h_CorrelationMatrixB = (TH1D*)f.Get("CorrelationMatrixB");
   h_CorrelationMatrixB->Draw("COLZ");
   c1->Print((outputpics+"corB.eps").c_str());
   c1->Clear();

  std::vector<std::string> MVAvars;


    MVAvars.push_back("MVA_TOPTOPLEPHAD");
    MVAvars.push_back("LepCharge");
    MVAvars.push_back("MVA_TOPTOPLEPHBB");
    MVAvars.push_back("MVA_TOPHLEPBB");
    MVAvars.push_back("MVA_TOPHLEPBB");
	  MVAvars.push_back("HiggsMass_TOPHLEPBB");
	  MVAvars.push_back("HiggsMass_TOPHLEPBB");
	  MVAvars.push_back("HiggsEta_TOPHLEPBB");
	  MVAvars.push_back("HiggsEta_TOPHLEPBB");
	  MVAvars.push_back("TopLepMass_TOPHLEPBB");
	  MVAvars.push_back("TopLepMass_TOPHLEPBB");
    MVAvars.push_back("TopLepPt_TOPHLEPBB");
    MVAvars.push_back("TopLepPt_TOPHLEPBB");
    MVAvars.push_back("TopLepEta_TOPHLEPBB");
    MVAvars.push_back("TopLepEta_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB");
    MVAvars.push_back("TopLepHiggsDr_TOPHLEPBB");
    MVAvars.push_back("TopLepHiggsDr_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPHLEPBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPHLEPBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPHLEPBB");
    MVAvars.push_back("TopHadMass_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepMass_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepTopHadDr_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepBJetCSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadBJetCSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("HiggsMass_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepMass_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepHiggsDr_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetCSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD");


    for(int i_vars = 0; i_vars < MVAvars.size(); i_vars++)
    {   
        TH1F *HistoNorm_S = 0;
        TH1F *HistoNorm_B = 0;
	      TH1F *histo_S( (TH1F*) f.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Signal").c_str()));
	      TH1F *histo_B( (TH1F*) f.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Background").c_str()));

        if(histo_S)
        {
          HistoNorm_S = new TH1F(*histo_S);
        }
        else
        {
            std::cout << "Input histo doesn't exist" << std::endl;
            continue;
        }
        if(histo_B)
        {
          HistoNorm_B = new TH1F(*histo_B);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;

        Double_t norm_S = 1;
        Double_t scale_S = norm_S/(HistoNorm_S->Integral());
        HistoNorm_S->Scale(scale_S);
        Double_t norm_B = 1;
        Double_t scale_B = norm_B/(HistoNorm_B->Integral());
        HistoNorm_B->Scale(scale_B);
        
        HistoNorm_S->GetXaxis()->SetTitle(MVAvars[i_vars].c_str());
        HistoNorm_S->SetLineColor(4);
        HistoNorm_B->SetLineColor(2);
        HistoNorm_S->SetFillColor(4);
        HistoNorm_B->SetFillColor(2);
        HistoNorm_S->SetFillStyle(3004);
        HistoNorm_B->SetFillStyle(3005);

        float maxhist= HistoNorm_S->GetMaximum();
        if (HistoNorm_B->GetMaximum() > maxhist) maxhist = HistoNorm_B->GetMaximum();
        
        maxhist = maxhist*1.4;
       
        HistoNorm_S->SetMaximum(maxhist);
        HistoNorm_B->SetMaximum(maxhist);

        TCanvas *c1 = new TCanvas();
        c1->cd();
        HistoNorm_S->Draw("hist");
        HistoNorm_B->Draw("hist same");

        TLegend *leg = new TLegend(0.7,0.75,0.98,0.95);
         leg->AddEntry(HistoNorm_S,"Signal","f");
         leg->AddEntry(HistoNorm_B,"Background","f");
         leg->Draw();



        c1->SaveAs((outputpics+MVAvars[i_vars]+"_Norm.png").c_str());

      delete c1;
    }

   
   f.Close();
}


std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}

// function that converts an int into a string
string floatToStr (float number)
{
  	ostringstream buff;
  	buff<<number;
  	return buff.str();
}

