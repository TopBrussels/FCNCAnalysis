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
    std::string TrainingName = "Training_" + SignalSample + channel + "_" +  category;//Example: Training_SThut_El_b3j3
    std::string outputpics = "pics/"+TrainingName+"/";
    mkdir("pics",0777);
    mkdir(outputpics.c_str(),0777);

   TFile f1(("../weights_CSVv2M/"+TrainingName+".root").c_str());
   std::string label1 = "Normal";
   TFile f2(("../weights___reducedvariables/"+TrainingName+".root").c_str());
   std::string label2 = "Reduced vars";

        TH1F *HistoNorm_S = 0;
        TH1F *HistoNorm_B = 0;
	      TH1F *histo_S( (TH1F*) f1.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));
	      TH1F *histo_B( (TH1F*) f2.Get("Method_BDT/BDT/MVA_BDT_rejBvsS"));

        if(histo_S)
        {
          HistoNorm_S = new TH1F(*histo_S);
        }
        else
        {
            std::cout << "Input histo doesn't exist" << std::endl;
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
        
//        HistoNorm_S->GetXaxis()->SetTitle(MVAvars[i_vars].c_str());
        HistoNorm_S->SetLineColor(4);
        HistoNorm_B->SetLineColor(2);
        HistoNorm_S->SetLineWidth(2);
        HistoNorm_B->SetLineWidth(2);
//        HistoNorm_S->SetFillColor(4);
//        HistoNorm_B->SetFillColor(2);
//        HistoNorm_S->SetFillStyle(3004);
//        HistoNorm_B->SetFillStyle(3005);

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
         leg->AddEntry(HistoNorm_S,label1.c_str(),"f");
         leg->AddEntry(HistoNorm_B,label2.c_str(),"f");
         leg->Draw();



        c1->SaveAs((outputpics+TrainingName+"ROCCurves_.png").c_str());

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

