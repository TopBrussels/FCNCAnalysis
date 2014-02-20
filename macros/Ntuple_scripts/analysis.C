#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Riostream.h"
#include <vector>
#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"


void analysis::Loop(){
  // running default loop:
  myAnalysis(0,45,0);
}

void analysis::myAnalysis(int nsel, int mode, bool silent)
{



  char plotName[300];

  
  
  if (nsel == 0)                	{sprintf(plotName,"TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root");  }
 
  
  
  
  
  char newRootFile[300];
  sprintf(newRootFile,"../results/an_%d.root", mode);
  
 
 
 
  TFile f_var(newRootFile, "UPDATE");
  
  
  if(!silent){
    std::cout << "[Info:] results root file " << newRootFile << std::endl;
  }
  
  
  //////////
  char title[300];
  
 
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {  // start eventloop
    

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
   /* 
    if (lumi != lum && nsel != 666 && mode !=3){
      if (jentry == 0)std::cout << "[Warning:] This tree was made with a different luminosity (" << lum << ") than " << lumi << std::endl;
      //xlWeight*=(lumi/lum);
    }
    
  */
  }// event loop.
  
  
  /*
  
  //Get the results
  
  if (!silent){ 
    cout << "------------------------------------------" << endl;
    cout << "[Results:] " << plotName <<  endl;
    cout << "------------------------------------------" << endl;  
    for (int i = 2; i < 9; i++){
      if (i == 2) cout << " leptons: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 3) cout << " inv. mass: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 4) cout << " met: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 5) cout << " jet: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 6) cout << " jet_bt: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
      if (i == 7) cout << " ht: " <<  histo->GetBinContent(i) << " +/-  " <<  histo->GetBinError(i)  << endl;
    }
    
   
    
    
  }
  */
  f_var.Write(newRootFile,TFile::kOverwrite);
  f_var.Close();
}
