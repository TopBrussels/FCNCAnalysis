//
//  Combine2DHistos.cpp
//  
//
//  Created by Shimaa on 06/09/16.
//
//

#include <stdio.h>

#include <stdio.h>
#include <vector>
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"

#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include "TStyle.h"
#include "TH2.h"
#include "TKey.h"
#include "setTDRStyle.C"
#include <cmath>
#include "TPaveText.h"
#include "TROOT.h"
#include "TRint.h"
#include "TSystemDirectory.h"


void Combine2DHistos(string channel = "Elec_Elec")
//void Draw2DPlots(string channel = "Elec_Elec")
{
    string runDate = "Test_Adding2DHistos_23Aug";
    string chan = "";
    chan += channel;
    const char *dirname;
    if(chan.find("Elec_Elec") == 0) dirname = "../Merged_Histos/ElecElec/";
    TString pathTohistoDir= dirname+runDate+"/";
    cout << " Hello 1 " << endl;
    //    TString Varhistos[] = {"2L_Nb_jets_vs_CSVLbjets",
    //        "2L_Nb_jets_vs_CSVMbjets",
    //        "2L_Nb_jets_vs_CSVTbjets"
    //    };
    
    TString histosName = "2L_Nb_jets_vs_CSVLbjets_ElEl_";
    
    cout << " Hello 2 " << endl;
    
    TString Samples_rootFiles[]= {
        pathTohistoDir+"merged_DataRun2015D.root",
        pathTohistoDir+"merged_DYJets10-50_amc.root",
        pathTohistoDir+"merged_DYJetsToLL_M-50toInf_Madgraph.root",
        pathTohistoDir+"merged_TTJets_madgraphMLM.root"//,
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        //        pathTohistoDir+"",
        
        
        
    };
    cout << " Hello 3 " << endl;
    
    TString LegendNames[] = {
        "DataRun2015D",
        "DY10-50_amc",
        "DYJetsM-50_MG",
        "TTJets_MG-MLM"//,
        
    };
    
    size_t Histos_Nb = 3;
    size_t SamplesNb = 4;
    
    
    vector<TFile*> Root_Files(0);
    for (size_t isample = 0; isample<SamplesNb; isample++)
    {
        cout << Samples_rootFiles[isample] << endl;
        Root_Files.push_back((TFile*)TFile::Open(Samples_rootFiles[isample]));
    }
    
 cout << " Hello 4 " << endl;
    vector<TH2F*> Histos(0);
    for (unsigned int ihisto = 0 ; ihisto < Root_Files.size(); ihisto++)
    {
        auto histptr = (TH2F*)Root_Files[ihisto]->Get(histosName);
        cout<< histosName.Data()<<endl;
        cout << histptr << endl;
        Histos.push_back(histptr);
        
    }
    cout <<"number of histos  = " << Histos.size() <<endl;
   
   
    TCanvas* c = new TCanvas("c","Overlay Canvas",800,800);
    TLegend* l = new TLegend(0.35,0.15,0.90,0.35);
    
    TH2F* CombHisto = new TH2F();
    TH2F* Temphisto = new TH2F();
    c->cd();
    
    int colors[] = {kRed,kBlue,kGreen,kMagenta};
    
    for (unsigned int i = 0 ; i<Histos.size(); i++)
    {
        cout << " Hello 4 " << endl;
        //Histos[i]->SetMarkerColor(colors[i]);
        Histos[i]->SetFillColor(i);
        cout << " Hello 5 " << endl;
        Histos[i]->SetMarkerStyle(1);
        cout << " Hello 5 " << endl;
        //CombHisto->Add(Histos[i]);
        if (i==0) {
            Histos[i]->Draw();
        }else{
            Histos[i]->Draw("same");
        }
        cout << " Hello 6 " << endl;
        l->AddEntry(Histos[i],LegendNames[i],"l");
    }
    
    //CombHisto->SetTitle(";#Jets; #CSVLbJets");
    //CombHisto->Draw();
    l->Draw("same");
    c->Update();
    c->SaveAs(histosName + "_2DCombinHisto.png");
                         
    
    
    
}

