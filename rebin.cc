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
    cout << "nbins before " << tempHisto->GetNbinsX() << endl;
    
    float newNbOfBins = 5;
    tempHisto->Rebin(newNbOfBins);
    cout << "nbins after " << tempHisto->GetNbinsX() << endl;
    outputfile->cd();
    tempHisto->Write();
    
   ;
  }
  
}