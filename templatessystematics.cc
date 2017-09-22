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



void templates(string coupling= "Zut", string region = "toppair", string temp = "TTTTZ"){
  TH1::SetDefaultSumw2();
  
  string inputfilename = "MVAoutput/outputtemplates/" + coupling +"_"+ region +"/Reader_"+coupling+"_"+region;
 
  if(temp.find("MTW")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/MTWtemplate/Reader_Zut_MTW";
  if(temp.find("STTTZ")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/STTTZtemplate/Reader_Zut_TTZ";
  if(temp.find("TTTTZ")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/TTTTZtemplate/Reader_Zut_TTZ";
  TFile* inputfile = new TFile((inputfilename+".root").c_str(),"UPDATE");
  if (!inputfile->IsOpen()) {
    cout << "can 't open " << inputfilename << endl;
    exit(1) ;
  }
  if(temp.find("STTTZ")!=std::string::npos  || temp.find("TTTTZ")!=std::string::npos) temp = "Zmass";
  vector <string> channels;
  channels.push_back("uuu");
  channels.push_back("uue");
  channels.push_back("eeu");
  channels.push_back("eee");
  
  vector <string> thesystlist;
  thesystlist.push_back(""); // nominal
  thesystlist.push_back("puSFDown");
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
   thesystlist.push_back("JESDown");
  
  
  const int nProcess =10;
  string processName[nProcess] = {"tHq","TTWJetsToLNu_amc","ttHToNonbb","ttHTobb","WZZ_amc","ZZZ_amc","tWll","STtW_top","STtW_atop"};
  
  
  TH1F* other;
  TH1F* ST_Zct;
  TH1F* ST_Zut;
  TH1F* TT_Zct;
  TH1F* TT_Zut;
  if(temp.find("BDT")!=std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        cout << "sys " << thesystlist[iSys] << endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        
        
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          cout << (coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str() << endl;
          if(iProcess == 0){ other = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str())->Clone("other");}
          else{
            other->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }
        other->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_other_80X"+systematic).c_str());
        other->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_other_80X"+systematic).c_str());
        
        if(coupling.find("Zct")!=std::string::npos){
          ST_Zct = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_ST_FCNC_zct_80X"+systematic).c_str())->Clone("ST_Zct");
          ST_Zct->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zct_80X"+systematic).c_str());
          ST_Zct->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zct_80X"+systematic).c_str());
          
          TT_Zct = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zct_80X"+systematic).c_str());
          TT_Zct->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zct_80X"+systematic).c_str());
        }
        
        if(coupling.find("Zut")!=std::string::npos){
          ST_Zut = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_ST_FCNC_zut_80X"+systematic).c_str())->Clone("ST_Zut");
          ST_Zut->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zut_80X"+systematic).c_str());
          ST_Zut->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zut_80X"+systematic).c_str());
          
          
          TT_Zut = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zut_80X"+systematic).c_str());
          TT_Zut->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zut_80X"+systematic).c_str());
        }
        
        inputfile->cd();
        other->Write();
        if(coupling.find("Zct")!=std::string::npos)ST_Zct->Write();
        if(coupling.find("Zut")!=std::string::npos)ST_Zut->Write();
        if(coupling.find("Zut")!=std::string::npos)TT_Zut->Write();
       if(coupling.find("Zct")!=std::string::npos) TT_Zct->Write();
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
        
        
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          cout << (temp+"_"+ channel +"_"+processName[iProcess]+ "_80X"+systematic).c_str() << endl;
          if(processName[iProcess].find("ttHTobb")!=std::string::npos ) continue;
          if(iProcess == 0){ other = (TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+ "_80X"+systematic).c_str())->Clone("other");}
          else{
            other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }
        //cout <<  "out " <<endl;
        // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str()));
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_tZq_amc_80X"+systematic).c_str()));
       // cout <<  "out " <<endl;
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str()));
       // cout <<  "out " <<endl;
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str()));
      //  cout <<  "out " <<endl;
        other->SetTitle((temp+"_"+ channel +"_other_80X"+systematic).c_str());
        other->SetName((temp+"_"+ channel +"_other_80X"+systematic).c_str());
        
        
          ST_Zct = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_ST_FCNC_zct_80X"+systematic).c_str())->Clone("ST_Zct");
          ST_Zct->SetTitle((temp+"_"+  channel +"_ST_Zct_80X"+systematic).c_str());
          ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X"+systematic).c_str());
          
          TT_Zct = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X"+systematic).c_str());
          TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X"+systematic).c_str());
        
        
      
          ST_Zut = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_ST_FCNC_zut_80X"+systematic).c_str())->Clone("ST_Zut");
          ST_Zut->SetTitle((temp+"_"+  channel +"_ST_Zut_80X"+systematic).c_str());
          ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X"+systematic).c_str());
          
          
          TT_Zut = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X"+systematic).c_str());
          TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X"+systematic).c_str());
        
        
        inputfile->cd();
        other->Write();
        ST_Zct->Write();
        ST_Zut->Write();
        TT_Zut->Write();
         TT_Zct->Write();

        
        inputfile->cd();
        other->Write();
      }
    }
  }
  else if(temp.find("Zmass")!=std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        cout << "sys " << thesystlist[iSys] << endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        
        
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          if(iProcess == 0){ other = (TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+ "_80X"+systematic).c_str())->Clone("other");}
          else{
            other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }
       /* other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_tZq_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str()));*/
        
        
        other->SetTitle((temp+"_"+ channel +"_other_80X"+systematic).c_str());
        other->SetName((temp+"_"+ channel +"_other_80X"+systematic).c_str());
        
       
          ST_Zct = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_ST_FCNC_zct_80X"+systematic).c_str())->Clone("ST_Zct");
          ST_Zct->SetTitle((temp+"_"+  channel +"_ST_Zct_80X"+systematic).c_str());
          ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X"+systematic).c_str());
          
          TT_Zct = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X"+systematic).c_str());
          TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X"+systematic).c_str());
            ST_Zut = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_ST_FCNC_zut_80X"+systematic).c_str())->Clone("ST_Zut");
          ST_Zut->SetTitle((temp+"_"+  channel +"_ST_Zut_80X"+systematic).c_str());
          ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X"+systematic).c_str());
          
          
          TT_Zut = (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfile->Get((temp+"_"+  channel +"_NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X"+systematic).c_str());
          TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X"+systematic).c_str());
        
        inputfile->cd();
        other->Write();
        ST_Zct->Write();
        ST_Zut->Write();
        TT_Zut->Write();
         TT_Zct->Write();
        
        inputfile->cd();
        other->Write();
      }
    }
  }
  
    
    
  }
