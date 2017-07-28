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



void rebin_fixed(string coupling= "Zut", string region = "singletop", string temp = "BDT"){
    TH1::SetDefaultSumw2();
  
 string inputfilename = "MVAoutput/outputtemplates/" + coupling +"_"+ region +"/Reader_"+coupling+"_"+region;
  
  if(temp.find("MTW")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/MTWtemplate/Reader_Zut_MTW";
  TFile* inputfile = new TFile((inputfilename+".root").c_str(),"READ");
  TFile* outputfile = new TFile((inputfilename+"_out.root").c_str(),"RECREATE");
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
  
  vector <string> histograms;
  vector <string> channels;
  channels.push_back("uuu");
  channels.push_back("uue");
  channels.push_back("eeu");
  channels.push_back("eee");
  
  vector <string> thesystlist;
  thesystlist.push_back(""); // nominal
  /*thesystlist.push_back("puSFDown");
  thesystlist.push_back("electronSFDown");
  thesystlist.push_back("muonSFDown");
  thesystlist.push_back("btagSFDown");
  thesystlist.push_back("btagSF_cferr1Down");
  thesystlist.push_back("btagSF_cferr2Down");
  thesystlist.push_back("btagSF_hfDown");
  thesystlist.push_back("btagSF_hfstats1Down");
  thesystlist.push_back("btagSF_hfstats2Down");
  thesystlist.push_back("btagSF_lfDown");
  thesystlist.push_back("btagSF_lfstats1Down");
  thesystlist.push_back("btagSF_lfstats2Down");
  thesystlist.push_back("JERUp");
  thesystlist.push_back("JESUp");
  thesystlist.push_back("puSFUp");
  thesystlist.push_back("electronSFUp");
  thesystlist.push_back("muonSFUp");
  thesystlist.push_back("btagSF_cferr1Up");
  thesystlist.push_back("btagSF_cferr2Up");
  thesystlist.push_back("btagSF_hfUp");
  thesystlist.push_back("btagSF_hfstats1Up");
  thesystlist.push_back("btagSF_hfstats2Up");
  thesystlist.push_back("btagSF_lfUp");
  thesystlist.push_back("btagSF_lfstats1Up");
  thesystlist.push_back("btagSF_lfstats2Up");
  thesystlist.push_back("JERDown");
  thesystlist.push_back("JESDown");*/


  const int nProcess =19;
  string processName[nProcess] = {"FakeMu","FakeEl","NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct","NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct","NP_overlay_ST_FCNC_zct","WZTo3LNu_amc","tZq_amc","tHq","TTWJetsToLNu_amc","TTZToLLNuNu_amc","ttHToNonbb","ttHTobb","ZZTo4L","WZZ_amc","ZZZ_amc","tWll","STtW_top","STtW_atop"};
  /*
  //empHisto = (TH1F*)inputfile->Get("Zct_BDT_singletop_uuu_WZTo3LNu_amc_80X");
  //inputfile->GetObject(obj->GetTitle(),tempHisto);
  nbins = tempHisto->GetNbinsX();
  cout << "nbins before " << tempHisto->GetNbinsX() << endl;
  entries_per_bins = tempHisto->GetEntries() /tempHisto->GetNbinsX();

  entries_per_bins = 20;
  cout << "Asking for " << entries_per_bins << " entries per bin" << endl;
  */
  bool foundfirstfakes = false;
  TH1F* tempHistoFake = 0;
  
  while ( (key = (TKey*)next() )) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
        && (!obj->InheritsFrom("TH2"))
        && (!obj->InheritsFrom("TH1"))
        ) {
      cout<< "Object "<< obj->GetName() << " is not 1D or 2D histogram :  will not be converted"<< endl;
    }
    inputfile->GetObject(obj->GetTitle(),tempHisto);
    //cout << "Histo name: "<< obj->GetName() <<" title: " << obj->GetTitle() << endl;
    histograms.push_back(obj->GetName());
    
    tempHisto = (TH1F*)inputfile->Get(obj->GetTitle());
    outputfile->cd();
    tempHisto->Write(obj->GetTitle());
    
    string namehist = (obj->GetTitle());
    if(namehist.find("Fake")!=std::string::npos){
      if(!foundfirstfakes){
        tempHistoFake = (TH1F*)inputfile->Get(obj->GetTitle());
      }
      else{
        
          tempHistoFake->Add((TH1F*)inputfile->Get(obj->GetTitle()));
        
      }
    }
    
    /*
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
    */
    
  }
  
  //TH1F* histoST = new TH1F("histoST","histoST",tempHisto->GetNbinsX(), tempHisto->GetXaxis()->GetXmin(), tempHisto->GetXaxis()->GetXmax());
  TH1F* STtW_top;
  TH1F* triboson;
  TH1F* ttH;
  TH1F* WZ;
  TH1F* FakeEl_nom;
  TH1F* FakeMu_nom;
  TH1F* FakeEl_check;
  TH1F* FakeMu_check;
  
  
 
  
  if(temp.find("MTW")==std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        cout << "sys " << thesystlist[iSys] << endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        
        if(channel.find("eee")!=std::string::npos || channel.find("uue")!=std::string::npos){
        FakeEl_nom = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeEl_80X"+systematic).c_str())->Clone("FakeEl_nom");
        FakeEl_check = (TH1F*)inputfile->Get(("check"+coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeEl_80X"+systematic).c_str())->Clone("FakeEl_check");
        FakeEl_nom->Scale(FakeEl_check->Integral()/FakeEl_nom->Integral());
        FakeEl_nom->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeEl_80X"+systematic).c_str());
        FakeEl_nom->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeEl_80X"+systematic).c_str());
        }
        else{
        FakeMu_nom = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeMu_80X"+systematic).c_str())->Clone("FakeMu_nom");
        FakeMu_check = (TH1F*)inputfile->Get(("check"+coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeMu_80X"+systematic).c_str())->Clone("FakeMu_check");
        FakeMu_nom->Scale(FakeMu_check->Integral()/FakeMu_nom->Integral());
        FakeMu_nom->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeMu_80X"+systematic).c_str());
        FakeMu_nom->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_FakeMu_80X"+systematic).c_str());
        }
        
        
        WZ= (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str())->Clone("WZ");
        WZ->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str());
        WZ->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZTo3LNu_80X"+systematic).c_str());
        
        
        for(int iBin = 0; iBin < WZ->GetNbinsX(); iBin++){
          if(WZ->GetBinContent(iBin) < 1) WZ ->SetBinContent(iBin,1.);
        }
        
        
        
        STtW_top= (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_STtW_top_80X"+systematic).c_str())->Clone("STtW_top");
        STtW_top->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_STtW_atop_80X"+systematic).c_str()));
        STtW_top->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_80X"+systematic).c_str());
        STtW_top->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_80X"+systematic).c_str());
        
        triboson = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZZZ_amc_80X"+systematic).c_str())->Clone("triboson");
        triboson->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZZ_amc_80X"+systematic).c_str()));
        triboson->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_triboson_80X"+systematic).c_str());
        triboson->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_triboson_80X"+systematic).c_str());
        
        
        ttH = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ttHToNonbb_80X"+systematic).c_str())->Clone("triboson");
        if(temp.find("MTW")==std::string::npos) ttH->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ttHTobb_80X"+systematic).c_str()));
        ttH->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ttH_80X"+systematic).c_str());
        ttH->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ttH_80X"+systematic).c_str());
        
        
        outputfile->cd();
        
        if(channel.find("eee")!=std::string::npos || channel.find("uue")!=std::string::npos){FakeEl_nom->Write();}
        else FakeMu_nom->Write();
        STtW_top->Write();
        WZ->Write();
        triboson->Write();
        ttH->Write();
      }
    }
  }
    else if(temp.find("MTW")!=std::string::npos){
      for(int iChan = 0; iChan < channels.size(); iChan++){
        string channel = channels[iChan];
        
        
        for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
          cout << "sys " << thesystlist[iSys] << endl;
          string systematic = "_" + thesystlist[iSys];
          if(iSys == 0) systematic = "";
          
          WZ= (TH1F*)inputfile->Get((temp+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str())->Clone("WZ");
          WZ->SetTitle((temp+"_"+ channel +"_WZ_80X"+systematic).c_str());
          WZ->SetName((temp+"_"+ channel +"_Z_80X"+systematic).c_str());
          
          
          STtW_top= (TH1F*)inputfile->Get((temp+"_"+ channel +"_STtW_top_80X"+systematic).c_str())->Clone("STtW_top");
          STtW_top->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_STtW_atop_80X"+systematic).c_str()));
          STtW_top->SetTitle((temp+"_"+ channel +"_ST_80X"+systematic).c_str());
          STtW_top->SetName((temp+"_"+ channel +"_ST_80X"+systematic).c_str());
          
          triboson = (TH1F*)inputfile->Get((temp+"_"+ channel +"_ZZZ_amc_80X"+systematic).c_str())->Clone("triboson");
          triboson->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_WZZ_amc_80X"+systematic).c_str()));
          triboson->SetTitle((temp+"_"+ channel +"_triboson_80X"+systematic).c_str());
          triboson->SetName((temp+"_"+ channel +"_triboson_80X"+systematic).c_str());
          
          
          ttH = (TH1F*)inputfile->Get((temp+"_"+ channel +"_ttHToNonbb_80X"+systematic).c_str())->Clone("triboson");
          if(temp.find("MTW")==std::string::npos) ttH->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_ttHTobb_80X"+systematic).c_str()));
          ttH->SetTitle((temp+"_"+ channel +"_ttH_80X"+systematic).c_str());
          ttH->SetName((temp+"_"+ channel +"_ttH_80X"+systematic).c_str());
          
          
          outputfile->cd();
          STtW_top->Write();
          WZ->Write();
          triboson->Write();
          ttH->Write();
        }
      }
    }

  
}