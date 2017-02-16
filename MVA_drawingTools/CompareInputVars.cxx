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
    std::string TrainingName_kevin = "Training_" + SignalSample + channel + "_" +  category;//Example: Training_SThut_El_b3j3
    std::string TrainingName_kirill = "TMVA_";
    if(SignalSample.find("SThut")!=string::npos) TrainingName_kirill += "HutST_all_" + category;
    else if(SignalSample.find("TThut")!=string::npos) TrainingName_kirill += "HutTT_all_" + category;
    else if(SignalSample.find("SThct")!=string::npos) TrainingName_kirill += "HctST_all_" + category;
    else if(SignalSample.find("TThct")!=string::npos) TrainingName_kirill += "HctTT_all_" + category;
    std::string outputpics = "pics_KirillKevin/"+TrainingName_kevin+"/";
    mkdir("pics_KirillKevin",0777);
    mkdir(outputpics.c_str(),0777);

   TFile f_kevin(("../weights_CSVv2M/"+TrainingName_kevin+".root").c_str());
   TFile f_kirill(("/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/tHFCNC/NtupleAnalyzer/test/MVA/"+TrainingName_kirill+".root").c_str());

   TH1D *h_sig_train_kevin = (TH1D*)f_kevin.Get("Method_BDT/BDT/MVA_BDT_Train_S");
   TH1D *h_sig_train_kirill = (TH1D*)f_kirill.Get("Method_BDT/BDT/MVA_BDT_Train_S");
   TH1D *h_bkg_train_kevin = (TH1D*)f_kevin.Get("Method_BDT/BDT/MVA_BDT_Train_B");
   TH1D *h_bkg_train_kirill = (TH1D*)f_kirill.Get("Method_BDT/BDT/MVA_BDT_Train_B");

   TCanvas *c1 = new TCanvas();
   h_sig_train_kevin->SetLineColor(kRed-4);
   h_sig_train_kevin->SetMarkerSize(0.7);
   h_sig_train_kevin->SetMarkerStyle(20);
   h_sig_train_kevin->SetMarkerColor(kRed-4);
   h_sig_train_kirill->SetLineColor(kRed-4);
   h_sig_train_kirill->SetMarkerSize(0.0);
//   h_sig_train_kirill->SetFillColor(kRed-4);
//   h_sig_train_kirill->SetFillStyle(3354);

   h_bkg_train_kevin->SetLineColor(kBlue-6);
   h_bkg_train_kevin->SetMarkerSize(0.7);
   h_bkg_train_kevin->SetMarkerStyle(22);
   h_bkg_train_kevin->SetMarkerColor(kBlue-6);
   h_bkg_train_kirill->SetLineColor(kBlue-6);
   h_bkg_train_kirill->SetMarkerSize(0.0);
   h_bkg_train_kirill->SetFillColor(kBlue-10);

   gStyle->SetHatchesSpacing(2.0);
   
   h_bkg_train_kirill->Draw("hist e1");
   h_sig_train_kirill->Draw("hist e1 same");
   
   h_sig_train_kevin->Draw("e1 same");
   h_bkg_train_kevin->Draw("e1 same");

   double max_sig_test = h_sig_train_kirill->GetMaximum();
   double max_bkg_test = h_bkg_train_kirill->GetMaximum();
   double max_sig_train = h_sig_train_kevin->GetMaximum();
   double max_bkg_train = h_bkg_train_kevin->GetMaximum();
   
   double max_test = std::max(max_sig_test,max_bkg_test);
   double max_train = std::max(max_sig_train,max_bkg_train);
   double max = std::max(max_test,max_train);
   
   h_bkg_train_kirill->SetMaximum(1.2*max);
   
   h_bkg_train_kirill->GetXaxis()->SetTitle("MVA output");
   h_bkg_train_kirill->GetYaxis()->SetTitle("(1/N) dN/dx");

   TLegend *leg0 = new TLegend(0.82,0.92,0.995,0.40);
   leg0->SetFillColor(253);
   leg0->SetBorderSize(0);
   leg0->AddEntry(h_sig_train_kevin,"Signal (kevin)","l");
   leg0->AddEntry(h_bkg_train_kevin,"Background (kevin)","l");
   leg0->AddEntry(h_sig_train_kirill,"Signal (kirill)","f");
   leg0->AddEntry(h_bkg_train_kirill,"Background (kirill)","f");
   leg0->Draw();
   
   c1->Print((outputpics+"disc.png").c_str());
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
        TH1F *h_sig_train__kevin = 0;
        TH1F *h_sig_train__kirill = 0;
        TH1F *h_bkg_train__kevin = (TH1F*) f_kevin.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Background").c_str());
        TH1F *h_bkg_train__kirill = (TH1F*) f_kirill.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Background").c_str());
	      TH1F *histo_S( (TH1F*) f_kevin.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Signal").c_str()));
	      TH1F *histo_B( (TH1F*) f_kirill.Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Signal").c_str()));

        if(histo_S)
        {
          h_sig_train__kevin = new TH1F(*histo_S);
        }
        else
        {
            std::cout << "Input histo doesn't exist" << std::endl;
            continue;
        }
        if(histo_B)
        {
          h_sig_train__kirill = new TH1F(*histo_B);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;

         TCanvas *c0 = new TCanvas();
				 c0->cd();
         h_sig_train__kevin->SetLineColor(kRed-4);
         h_sig_train__kevin->SetMarkerSize(0.7);
         h_sig_train__kevin->SetMarkerStyle(20);
         h_sig_train__kevin->SetMarkerColor(kRed-4);
         h_sig_train__kirill->SetLineColor(kRed-4);
         h_sig_train__kirill->SetMarkerSize(0.0);
      //   h_sig_train__kirill->SetFillColor(kRed-4);
      //   h_sig_train__kirill->SetFillStyle(3354);

         h_bkg_train__kevin->SetLineColor(kBlue-6);
         h_bkg_train__kevin->SetMarkerSize(0.7);
         h_bkg_train__kevin->SetMarkerStyle(22);
         h_bkg_train__kevin->SetMarkerColor(kBlue-6);
         h_bkg_train__kirill->SetLineColor(kBlue-6);
         h_bkg_train__kirill->SetMarkerSize(0.0);
         h_bkg_train__kirill->SetFillColor(kBlue-10);
				 
				 h_bkg_train__kirill->Scale(1./h_bkg_train__kirill->Integral());
				 h_bkg_train__kevin->Scale(1./h_bkg_train__kevin->Integral());
				 h_sig_train__kirill->Scale(1./h_sig_train__kirill->Integral());
				 h_sig_train__kevin->Scale(1./h_sig_train__kevin->Integral());

         gStyle->SetHatchesSpacing(2.0);
         
         h_bkg_train__kirill->Draw("hist e1");
         h_sig_train__kirill->Draw("hist e1 same");
         
         h_sig_train__kevin->Draw("e1 same");
         h_bkg_train__kevin->Draw("e1 same");

         double max_sig_test = h_sig_train__kirill->GetMaximum();
         double max_bkg_test = h_bkg_train__kirill->GetMaximum();
         double max_sig_train_ = h_sig_train__kevin->GetMaximum();
         double max_bkg_train_ = h_bkg_train__kevin->GetMaximum();
         
         double max_test = std::max(max_sig_test,max_bkg_test);
         double max_train_ = std::max(max_sig_train_,max_bkg_train_);
         double max = std::max(max_test,max_train_);
         
         h_bkg_train__kirill->SetMaximum(1.2*max);
         
         h_bkg_train__kirill->GetXaxis()->SetTitle(MVAvars[i_vars].c_str());
         h_bkg_train__kirill->GetYaxis()->SetTitle("norm. events");

         TLegend *leg = new TLegend(0.82,0.92,0.995,0.40);
         leg->SetFillColor(253);
         leg->SetBorderSize(0);
         leg->AddEntry(h_sig_train__kevin,"Signal (kevin)","l");
         leg->AddEntry(h_bkg_train__kevin,"Background (kevin)","l");
         leg->AddEntry(h_sig_train__kirill,"Signal (kirill)","f");
         leg->AddEntry(h_bkg_train__kirill,"Background (kirill)","f");
         leg->Draw();

        c0->SaveAs((outputpics+MVAvars[i_vars]+"_Norm.png").c_str());

      delete c0;
    }
   
   f_kevin.Close();
   f_kirill.Close();
}


std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}

