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

    SetPlotStyle();
   
    std::string outputpics = "pics_TopKinFit/";
    mkdir("pics_TopKinFit",0777);
    mkdir(outputpics.c_str(),0777);
    
    string location = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_25/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/";

   TFile f1((location+"TMVAFullRecoTruth.root").c_str());
   std::string label1 = "Truth";
   TFile f2((location+"TMVAFullRecoHighestCSVv2.root").c_str());
   std::string label2 = "HighestCSVv2";
   TFile f3((location+"TMVAFullRecoCSVv2L.root").c_str());
   std::string label3 = "CSVv2L";
   TFile f4((location+"TMVAFullRecoCSVv2M.root").c_str());
   std::string label4 = "CSVv2M";
   TFile f5((location+"TMVAFullRecoCSVv2T.root").c_str());
   std::string label5 = "CSVv2T";
   TFile f6((location+"TMVAFullRecoAll.root").c_str());
   std::string label6 = "All";

        TH1F *HistoNorm_1 = 0;
        TH1F *HistoNorm_2 = 0;
        TH1F *HistoNorm_3 = 0;
        TH1F *HistoNorm_4 = 0;
        TH1F *HistoNorm_5 = 0;
        TH1F *HistoNorm_6 = 0;
	      TH1F *histo_1( (TH1F*) f1.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_2( (TH1F*) f2.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_3( (TH1F*) f3.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_4( (TH1F*) f4.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_5( (TH1F*) f5.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_6( (TH1F*) f6.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));

        if(histo_1)
        {
          HistoNorm_1 = new TH1F(*histo_1);
        }
        else
        {
            std::cout << "Input histo doesn't exist" << std::endl;
        }
        if(histo_2)
        {
          HistoNorm_2 = new TH1F(*histo_2);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;
        if(histo_3)
        {
          HistoNorm_3 = new TH1F(*histo_3);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;
        if(histo_4)
        {
          HistoNorm_4 = new TH1F(*histo_4);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;
        if(histo_5)
        {
          HistoNorm_5 = new TH1F(*histo_5);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;
        if(histo_6)
        {
          HistoNorm_6 = new TH1F(*histo_6);
        }
        else std::cout << "Input histo doesn't exist" << std::endl;
/*
cout <<  label1 << " " << HistoNorm_1->Integral() << endl;
cout <<  label2 << " " << HistoNorm_2->Integral() << endl;
cout <<  label3 << " " << HistoNorm_3->Integral() << endl;
*/
/*
        Double_t norm_S = 1;
        Double_t scale_S = norm_S/(HistoNorm_1->Integral());
        HistoNorm_1->Scale(1/HistoNorm_1->Integral());
        Double_t norm_B = 1;
        Double_t scale_B = norm_B/(HistoNorm_2->Integral());
        HistoNorm_2->Scale(1/HistoNorm_2->Integral());
        Double_t norm_C = 1;
        Double_t scale_C = norm_C/(HistoNorm_3->Integral());
        HistoNorm_3->Scale(scale_C);
*/

        
        HistoNorm_1->SetTitle("Full reconstruction");
        HistoNorm_1->GetXaxis()->SetRangeUser(0.1,0.95);
        HistoNorm_1->GetYaxis()->SetRangeUser(0.3,1.2);
        HistoNorm_1->SetTitle("");
        HistoNorm_1->SetLineColor(1);
        HistoNorm_2->SetLineColor(2);
        HistoNorm_3->SetLineColor(3);
        HistoNorm_4->SetLineColor(4);
        HistoNorm_5->SetLineColor(6);
        HistoNorm_6->SetLineColor(7);
        HistoNorm_1->SetLineWidth(2);
        HistoNorm_2->SetLineWidth(2);
        HistoNorm_3->SetLineWidth(2);
        HistoNorm_4->SetLineWidth(2);
        HistoNorm_5->SetLineWidth(2);
        HistoNorm_6->SetLineWidth(2);
//        HistoNorm_1->SetFillColor(4);
//        HistoNorm_2->SetFillColor(2);
//        HistoNorm_1->SetFillStyle(3004);
//        HistoNorm_2->SetFillStyle(3005);
/*
        float maxhist= HistoNorm_1->GetMaximum();
        if (HistoNorm_2->GetMaximum() > maxhist) maxhist = HistoNorm_2->GetMaximum();
        if (HistoNorm_3->GetMaximum() > maxhist) maxhist = HistoNorm_3->GetMaximum();
        
        maxhist = maxhist*1.4;
       
        HistoNorm_1->SetMaximum(maxhist);
        HistoNorm_2->SetMaximum(maxhist);
        HistoNorm_3->SetMaximum(maxhist);
*/
        TCanvas *c1 = new TCanvas();
        c1->cd();
        HistoNorm_1->Draw("hist");
        HistoNorm_2->Draw("hist same");
        HistoNorm_3->Draw("hist same");
        HistoNorm_4->Draw("hist same");
        HistoNorm_5->Draw("hist same");
        HistoNorm_6->Draw("hist same");

        TLegend *leg = new TLegend(0.7,0.65,0.98,0.95);
         leg->AddEntry(HistoNorm_1,label1.c_str(),"f");
         leg->AddEntry(HistoNorm_2,label2.c_str(),"f");
         leg->AddEntry(HistoNorm_3,label3.c_str(),"f");
         leg->AddEntry(HistoNorm_4,label4.c_str(),"f");
         leg->AddEntry(HistoNorm_5,label5.c_str(),"f");
         leg->AddEntry(HistoNorm_6,label6.c_str(),"f");
         leg->Draw();

//        c1->SetLogy();
//        c1->SetLogx();

        c1->SaveAs((outputpics+"ROCCurves_TopKinFit_FullReco.eps").c_str());

      delete c1;

   
   f1.Close();
   f2.Close();
}


std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}

