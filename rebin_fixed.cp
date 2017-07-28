//
//  rebin.cpp
//  
//
//  Created by Isis Van Parijs on 16/07/2017.
//
//

#include <stdio.h>
#include "TFile.h"
#include "TString.h"



void rebin(){
  
  TString inputfilename = "MVAoutput/outputtemplates/Zct_singletop/Reader_Zct_singletop";
  
  TFile* inputfile = new TFile((inputfilename+".root").Data(),"READ");
  TFile* outputfile = new TFile((inputfilename+"_out.root").Data(),"RECREATE");
  if (!inputfile->IsOpen()) {
    cout << "can 't open " << inputfilename << endl;
    exit(1) ;
  }
  
  TList* listofHistograms = inputfile->GetListOfKeys() ;
  if (!listofHistograms) { cout << " No keys found in file " << inputfilename << endl;  ; exit(1) ; }
  TIter next(listofHistograms) ;
  TKey* key ;
  TObject* obj ;
  TH1F* tempHisto;
  double entries_per_bins = 0;
  int nbins = -5;
  
  tempHisto = (TH1F*)inputfile->Get("Zct_BDT_singletop_uuu_WZTo3LNu_amc_80X");
  //inputfile->GetObject(obj->GetTitle(),tempHisto);
  nbins = tempHisto->GetNbinsX();
  cout << "nbins before " << tempHisto->GetNbinsX() << endl;
  entries_per_bins = tempHisto->GetEntries() /tempHisto->GetNbinsX();

  entries_per_bins = 20;
  cout << "Asking for " << entries_per_bins << " entries per bin" << endl;
  
  while ( (key = (TKey*)next() )) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
        && (!obj->InheritsFrom("TH2"))
        && (!obj->InheritsFrom("TH1"))
        ) {
      cout<< "Object "<< obj->GetName() << " is not 1D or 2D histogram :  will not be converted"<< endl;
    }
    inputfile->GetObject(obj->GetTitle(),tempHisto);
    cout << "Histo name: "<< obj->GetName() <<" title: " << obj->GetTitle() << endl;
    
    Double_t *xq = new Double_t[nbins+1];
    Double_t *content = new Double_t[nbins+1];
    xq[0] = tempHisto->GetXaxis()->GetXmin();
    content[0] = tempHisto->GetBinContent(0); // underflow
    xq[nbins] = tempHisto->GetXaxis()->GetXmax();
    content[nbins] = tempHisto->GetBinContent(nbins); // overflow
    
    double entries_so_far = 0;
    int new_i = 1;
    for(int i = 1; i < nbins; i++) {
      entries_so_far += tempHisto->GetBinContent(i);
      cout << entries_so_far << " / " << entries_per_bins << endl;
      if (entries_so_far >= entries_per_bins || (i==(nbins-1))){
        xq[ new_i ] = tempHisto->GetXaxis()->GetBinUpEdge( i );
        content[new_i] = entries_so_far;
        new_i++;
        entries_so_far = 0;
      }
    }
    
    TH1F* historebinned = (TH1F*) tempHisto->Rebin(new_i-1,"historebinned",xq);
    historebinned->SetName(tempHisto->GetName());
    historebinned->SetTitle(tempHisto->GetTitle());
    
    cout << "nbins after " << tempHisto->GetNbinsX() << endl;
    outputfile->cd();
    historebinned->Write(tempHisto->GetTitle());
    
    
  }
  
}