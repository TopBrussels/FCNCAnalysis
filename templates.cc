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
#include "TH1.h"
#include <iostream>

using namespace std;

void templates(string coupling= "Zut", string region = "singletop", string temp = "BDT"){
  TH1::SetDefaultSumw2();
  
  string inputfilename = "MVAoutput/outputtemplates/" + coupling +"_"+ region +"/Reader_"+coupling+"_"+region;
 
  bool isSingle = false;
  if(temp.find("MTW")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/MTWtemplate/Reader_Zut_MTW";
  if(temp.find("STTTZ")!=std::string::npos){ inputfilename ="MVAoutput/outputtemplates/STTTZtemplate/Reader_Zut_TTZ"; isSingle = true;}
  if(temp.find("TTTTZ")!=std::string::npos) inputfilename ="MVAoutput/outputtemplates/TTTTZtemplate/Reader_Zut_TTZ";
  TFile* inputfile = new TFile((inputfilename+".root").c_str(),"UPDATE");
  if (!inputfile->IsOpen()) {
    std::cout << "can 't open " << inputfilename << std::endl;
    exit(1) ;
  }
  cout << "input filename " <<(inputfilename+".root").c_str()<< endl;
  string outputfilename = "MVAoutput/outputtemplates/" + coupling +"_"+ region +"/Reader_"+coupling+"_"+region;
  
  if(temp.find("MTW")!=std::string::npos) outputfilename ="MVAoutput/outputtemplates/MTWtemplate/Reader_Zut_MTW";
  if(temp.find("STTTZ")!=std::string::npos) outputfilename ="MVAoutput/outputtemplates/STTTZtemplate/Reader_Zut_TTZ";
  if(temp.find("TTTTZ")!=std::string::npos) outputfilename ="MVAoutput/outputtemplates/TTTTZtemplate/Reader_Zut_TTZ";

 /* TFile* outputfile = new TFile((outputfilename+"_new.root").c_str(),"RECREATE");
  if (!outputfile->IsOpen()) {
    std::cout << "can 't open " << outputfilename<< std::endl;
    exit(1) ;
  }*/
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
 // const int nbin = 8;
 // double xbins_BDT_singletop[nbin] =  {-1.,-0.2,0.,0.4,0.6,0.7,0.8,1.};
  
  
  
  
  TH1F* other;
  TH1F* ST_Zct;
  TH1F* ST_Zut;
  TH1F* TT_Zct;
  TH1F* TT_Zut;
  TH1F* signal;
  TH1F* bkg;
  TH1F* total;
  TH1F* wfake;
  TH1F* data;
  TH1F* zfake;
  
  TH1F* wz_new;
  TH1F* zz_new;
  TH1F* tzq_new;
  TH1F* ttz_new;
  TH1F* data_new;
  TH1F* other_new;
  TH1F* ST_Zct_new;
  TH1F* ST_Zut_new;
  TH1F* TT_Zct_new;
  TH1F* TT_Zut_new;
  TH1F* wfake_new;
  TH1F* bkg_new;
  TH1F* zfake_new;
  if(temp.find("BDT")!=std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      // determine binning
    /*  double targetbincontent = 5.;
      std::vector<double> vbinx;
      std::vector<double> vbinc;
      

       for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
        std::cout << (coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str() << std::endl;
        if(iProcess == 0){ total = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str())->Clone("other");}
        else{
          total->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
        }
         
      }
      vbinx.push_back(total->GetBinLowEdge(0));
      
      
      
      vbinx.push_back(total->GetBinUpperEdge(total->GetNbinsX()));
      */
      
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        std::cout << "sys " << thesystlist[iSys] << std::endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        
        
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          std::cout << (coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str() << std::endl;
          if(iProcess == 0){ other = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str())->Clone("other");}
          else{
            other->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }
        bkg = (TH1F*) other->Clone("TotalBkg");
        bkg->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"tZq_amc"+"_80X"+systematic).c_str())->Clone("tZq"));
        bkg->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"WZTo3LNu_amc"+"_80X"+systematic).c_str())->Clone("WZ"));
        bkg->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"ZZTo4L"+"_80X"+systematic).c_str())->Clone("ZZ"));
        bkg->Add((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"TTZToLLNuNu_amc"+"_80X"+systematic).c_str())->Clone("ttZ"));
        std::cout << "got all bkgs" <<std::endl;
        
        other_new = (TH1F*) other->Clone(""); //Rebin(nbin-1,"other",xbins_BDT_singletop);
        other_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_other_80X"+systematic).c_str());
        other_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_other_80X"+systematic).c_str());
        bkg_new = (TH1F*) bkg->Clone(""); //Rebin(nbin-1,"bkg",xbins_BDT_singletop);
        bkg_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalBackground_80X"+systematic).c_str());
        bkg_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalBackground_80X"+systematic).c_str());
        std::cout << "got all bkgs" <<std::endl;
        if(coupling.find("Zct")!=std::string::npos){
          ST_Zct = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_ST_FCNC_zct_80X"+systematic).c_str())->Clone("ST_Zct");
          ST_Zct_new = (TH1F*) ST_Zct->Clone(""); //Rebin(nbin-1,"ST_Zct",xbins_BDT_singletop);
          ST_Zct_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zct_80X"+systematic).c_str());
          ST_Zct_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zct_80X"+systematic).c_str());
          
          TT_Zct = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X"+systematic).c_str())->Clone("TT_Zct"));
          
          TT_Zct_new = (TH1F*) TT_Zct->Clone("");//Rebin(nbin-1,"TT_Zct",xbins_BDT_singletop);
          TT_Zct_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zct_80X"+systematic).c_str());
          TT_Zct_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zct_80X"+systematic).c_str());
        }
        std::cout << "got all sig" <<std::endl;
        if(coupling.find("Zut")!=std::string::npos){
          ST_Zut = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_ST_FCNC_zut_80X"+systematic).c_str())->Clone("ST_Zut");
          ST_Zut_new = (TH1F*) ST_Zut->Clone("");//Rebin(nbin-1,"ST_Zut",xbins_BDT_singletop);
          ST_Zut_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zut_80X"+systematic).c_str());
          ST_Zut_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ST_Zut_80X"+systematic).c_str());
          
          
          TT_Zut = (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X"+systematic).c_str())->Clone("TT_Zut"));
          TT_Zut_new = (TH1F*) TT_Zut->Clone(""); //Rebin(nbin-1,"TT_Zut",xbins_BDT_singletop);
          TT_Zut_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zut_80X"+systematic).c_str());
          TT_Zut_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TT_Zut_80X"+systematic).c_str());
        }
        std::cout << "got all sig" <<std::endl;
        if(iSys == 0){
        if( channel.find("uuu")!=std::string::npos || channel.find("eeu")!=std::string::npos) wfake = ((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"WFakeMu"+"_80X"+systematic).c_str())->Clone("wfake"));
        if( channel.find("uue")!=std::string::npos || channel.find("eee")!=std::string::npos) wfake = ((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"WFakeEl"+"_80X"+systematic).c_str())->Clone("wfake"));
        
          wfake_new = (TH1F*) wfake->Clone(""); //Rebin(nbin-1,"wfake",xbins_BDT_singletop);
        if( channel.find("uuu")!=std::string::npos || channel.find("eeu")!=std::string::npos){
        wfake_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WFakeMu_80X"+systematic).c_str());
        wfake_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WFakeMu_80X"+systematic).c_str());
        }
        if( channel.find("uue")!=std::string::npos || channel.find("eee")!=std::string::npos){
          wfake_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WFakeEl_80X"+systematic).c_str());
          wfake_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WFakeEl_80X"+systematic).c_str());
        }
        
        if( channel.find("uuu")!=std::string::npos || channel.find("uue")!=std::string::npos) zfake = ((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"ZFakeMu"+"_80X"+systematic).c_str())->Clone("zfake"));
        if( channel.find("eeu")!=std::string::npos || channel.find("eee")!=std::string::npos) zfake = ((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"ZFakeEl"+"_80X"+systematic).c_str())->Clone("zfake"));
        
          zfake_new = (TH1F*) zfake->Clone(""); //Rebin(nbin-1,"zake",xbins_BDT_singletop);
        if( channel.find("uuu")!=std::string::npos || channel.find("uue")!=std::string::npos){
          zfake_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZFakeMu_80X"+systematic).c_str());
          zfake_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZFakeMu_80X"+systematic).c_str());
        }
        if( channel.find("eeu")!=std::string::npos || channel.find("eee")!=std::string::npos){
          zfake_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZFakeEl_80X"+systematic).c_str());
          zfake_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZFakeEl_80X"+systematic).c_str());
        }
        
        
       
          data = ((TH1F*)inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"data").c_str())->Clone("data"));
          data_new = (TH1F*) data->Clone(""); //Rebin(nbin-1,"fake",xbins_BDT_singletop);
          data_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_data").c_str());
          data_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_data").c_str());
        }
        
        if(coupling.find("Zct")!=std::string::npos){ signal = (TH1F*) ST_Zct_new->Clone("TotalSig"); signal->Add((TH1F*) TT_Zct_new->Clone(""));
         
        }
        else if(coupling.find("Zut")!=std::string::npos){ signal = (TH1F*) ST_Zut_new->Clone("TotalSig"); signal->Add((TH1F*)TT_Zut_new->Clone(""));}
        signal->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalSignal_80X"+systematic).c_str());
        signal->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalSignal_80X"+systematic).c_str());
        
        total = (TH1F*) bkg_new->Clone("Total");
        total ->Add((TH1F*)signal->Clone(""));
        total->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalProcesses_80X"+systematic).c_str());
        total->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TotalProcesses_80X"+systematic).c_str());
        
        ttz_new = (TH1F*)((TH1F*) inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"TTZToLLNuNu_amc"+"_80X"+systematic).c_str()))->Clone(""); //Rebin(nbin-1,"ttz",xbins_BDT_singletop);
        ttz_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str());
        ttz_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str());
        
        wz_new = (TH1F*)((TH1F*) inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"WZTo3LNu_amc"+"_80X"+systematic).c_str()))->Clone(""); //Rebin(nbin-1,"wz",xbins_BDT_singletop);
        wz_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str());
        wz_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str());
        
        tzq_new =(TH1F*) ((TH1F*) inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"tZq_amc"+"_80X"+systematic).c_str()))->Clone(""); //Rebin(nbin-1,"tzq",xbins_BDT_singletop);
        tzq_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_tZq_amc_80X"+systematic).c_str());
        tzq_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_tZq_amc_80X"+systematic).c_str());
        
        zz_new =(TH1F*) ((TH1F*) inputfile->Get((coupling+"_"+temp+"_"+ region+"_"+ channel +"_"+"ZZTo4L"+"_80X"+systematic).c_str()))->Clone(""); //Rebin(nbin-1,"zz",xbins_BDT_singletop);
        zz_new->SetTitle((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str());
        zz_new->SetName((coupling+"_"+temp+"_"+ region+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str());
        
        
       // outputfile->cd();
        inputfile->cd();
        other_new->Write();
        if(coupling.find("Zct")!=std::string::npos)ST_Zct_new->Write();
        if(coupling.find("Zut")!=std::string::npos)ST_Zut_new->Write();
        if(coupling.find("Zut")!=std::string::npos)TT_Zut_new->Write();
       if(coupling.find("Zct")!=std::string::npos) TT_Zct_new->Write();
        bkg_new->Write();
        signal->Write();
        total->Write();
        if(iSys == 0 )wfake_new->Write();
        if(iSys == 0 )zfake_new->Write();
        tzq_new->Write();
        zz_new->Write();
        wz_new->Write();
        ttz_new->Write();
        if(iSys == 0 ) data_new->Write();
      }
    }
  }
  else if(temp.find("MTW")!=std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        std::cout << "sys " << thesystlist[iSys] << std::endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        if(iSys == 0){
          
          data_new = ((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+"data").c_str())->Clone("data"));
          data_new->SetTitle((temp+"_"+ channel +"_data").c_str());
          data_new->SetName((temp+"_"+ channel +"_data").c_str());
          
        }
        
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          std::cout << (temp+"_"+ channel +"_"+processName[iProcess]+ "_80X"+systematic).c_str() << std::endl;
          if(processName[iProcess].find("ttHTobb")!=std::string::npos ) continue;
          if(iProcess == 0){ other = (TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+ "_80X"+systematic).c_str())->Clone("other");}
          else{
            other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }
        //std::cout <<  "out " <<std::endl;
        // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str()));
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_tZq_amc_80X"+systematic).c_str()));
       // std::cout <<  "out " <<std::endl;
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str()));
       // std::cout <<  "out " <<std::endl;
       // other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str()));
      //  std::cout <<  "out " <<std::endl;
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
        
        
       // outputfile->cd();
         inputfile->cd();
        other->Write();
        ST_Zct->Write();
        ST_Zut->Write();
        TT_Zut->Write();
        TT_Zct->Write();
        if(iSys == 0) data_new->Write();

      }
    }
  }
  else if(temp.find("Zmass")!=std::string::npos){
    for(int iChan = 0; iChan < channels.size(); iChan++){
      string channel = channels[iChan];
      
      
      for(int iSys = 0; iSys < thesystlist.size() ; iSys++){
        std::cout << "sys " << thesystlist[iSys] << std::endl;
        string systematic = "_" + thesystlist[iSys];
        if(iSys == 0) systematic = "";
        if(iSys == 0){
          
          data_new = ((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+"data").c_str())->Clone("data"));
          data_new->SetTitle((temp+"_"+ channel +"_data").c_str());
          data_new->SetName((temp+"_"+ channel +"_data").c_str());
          
        }
        for(int iProcess = 0 ; iProcess < nProcess-1 ; iProcess++){
          if(iProcess == 0 ){
             cout << (temp+"_"+channel +"_"+processName[iProcess]+"_80X"+systematic).c_str() << endl;
            other = (TH1F*)inputfile->Get((temp+"_"+channel+"_"+processName[iProcess]+"_80X"+systematic).c_str())->Clone("other");}
          else{
            cout << (temp+"_"+channel+"_"+processName[iProcess]+"_80X"+systematic).c_str() << endl;
            other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_"+processName[iProcess]+"_80X"+systematic).c_str()));
          }
        }

       /* other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_WZTo3LNu_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_tZq_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_TTZToLLNuNu_amc_80X"+systematic).c_str()));
        other->Add((TH1F*)inputfile->Get((temp+"_"+ channel +"_ZZTo4L_80X"+systematic).c_str()));*/
        cout << "checking" << endl;
        
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
        
       // outputfile->cd();
         inputfile->cd();
        other->Write();
        ST_Zct->Write();
        ST_Zut->Write();
        TT_Zut->Write();
         TT_Zct->Write();
        if(iSys == 0) data_new->Write();
      }
    }
  }
  
  string inputfilenamePDF = "PDFs/root/";
  cout << "PDFs" << endl;
  vector <string> pdfsys = {"PDFEnvelopeDown", "PDFEnvelopeUp"};
  for(int isys = 0; isys < 2; isys++){
    string syst = pdfsys[isys];
    
    for(int iCHan = 0; iCHan < channels.size(); iCHan++){
       string channel = channels[iCHan];
      if(temp.find("BDT")!=std::string::npos && region.find("singletop")!=std::string::npos){
        if(coupling.find("Zct")!=std::string::npos){
          /*
         
          //cout << "Get PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_BDT_singletop_Zctttz.root" << endl;
          TFile* inputfilePDF_STNP = new TFile("PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_BDT_singletop_Zctttz.root","READ");
          //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
          ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
          ST_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
          
          
          //cout << "Get PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_singletop_Zctttz.root" << endl;
          TFile* inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_singletop_Zctttz.root" ,"READ");
          //cout << "Get PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_singletop_Zctttz.root" << endl;
          TFile* inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_singletop_Zctttz.root" ,"READ");
          
          TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          TT_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zct->Clone("TotalSig"); signal->Add((TH1F*) TT_Zct->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          
          TFile* inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_BDT_singletop_Zctttz.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          /*
          TFile* inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_BDT_singletop_Zctttz.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ttZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_zz_new = new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_BDT_singletop_Zctttz.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ZZ");
         zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_BDT_singletop_Zctttz.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("tZq");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          
         // outputfile->cd();
           inputfile->cd();
         // ST_Zct->Write();
         // signal->Write();
        //  TT_Zct->Write();
          wz_new->Write();
          zz_new->Write();
        //  ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
        //  inputfilePDF_ttz_new->Close();
        //  inputfilePDF_TT->Close();
        //  inputfilePDF_STNP->Close();
         // inputfilePDF_aTT->Close();
        }
        else if(coupling.find("Zut")!=std::string::npos){
          
          
          /*  TFile* inputfilePDF_STNP = new TFile("PDFs/root/NP_overlay_ST_FCNC_zut_80XPDFhistograms_BDT_singletop_Zutttz.root","READ");
           ST_Zut = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
          ST_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XPDFhistograms_BDT_singletop_Zutttz.root" ,"READ");
          TFile* inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XPDFhistograms_BDT_singletop_Zutttz.root " ,"READ");
          
          TT_Zut = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          TT_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zut->Clone("TotalSig"); signal->Add((TH1F*) TT_Zut->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_BDT_singletop_Zutttz.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
         /*
          TFile* inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_BDT_singletop_Zutttz.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_zz_new = new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_BDT_singletop_Zutttz.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_BDT_singletop_Zutttz.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          //outputfile->cd();
           inputfile->cd();
         /* ST_Zut->Write();
          signal->Write();
          TT_Zut->Write();*/
          wz_new->Write();
          zz_new->Write();
        //  ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
         // inputfilePDF_ttz_new->Close();
         /* inputfilePDF_TT->Close();
          inputfilePDF_STNP->Close();
          inputfilePDF_aTT->Close();*/
        }
        

        
      }
      else if(temp.find("BDT")!=std::string::npos && region.find("toppair")!=std::string::npos){
        if(coupling.find("Zct")!=std::string::npos){
          
          /*
          //cout << "Get PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_BDT_toppair_Zctttz.root" << endl;
          TFile* inputfilePDF_STNP = new TFile("PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_BDT_toppair_Zctttz.root","READ");
          //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
          ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
          ST_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
          
          
          //cout << "Get PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_toppair_Zctttz.root" << endl;
          TFile* inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_toppair_Zctttz.root" ,"READ");
          //cout << "Get PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_toppair_Zctttz.root" << endl;
          TFile* inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_BDT_toppair_Zctttz.root" ,"READ");
          
          TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          TT_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zct->Clone("TotalSig"); signal->Add((TH1F*) TT_Zct->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
           TFile* inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_BDT_toppair_Zctttz.root","READ");
           wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
         /*
          TFile* inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_BDT_toppair_Zctttz.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_zz_new = new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_BDT_toppair_Zctttz.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_BDT_toppair_Zctttz.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          //outputfile->cd();
           inputfile->cd();
         /* ST_Zct->Write();
          signal->Write();*/
         // TT_Zct->Write();
           wz_new ->Write();
          zz_new->Write();
          //ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
          //inputfilePDF_ttz_new->Close();
          /*inputfilePDF_TT->Close();
          inputfilePDF_STNP->Close();
          inputfilePDF_aTT->Close();*/
        }
        else if(coupling.find("Zut")!=std::string::npos){
          
          /*
          TFile* inputfilePDF_STNP = new TFile("PDFs/root/NP_overlay_ST_FCNC_zut_80XPDFhistograms_BDT_toppair_Zutttz.root","READ");
          ST_Zut = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
          ST_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XPDFhistograms_BDT_toppair_Zutttz.root" ,"READ");
          TFile* inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XPDFhistograms_BDT_toppair_Zutttz.root " ,"READ");
          
          TT_Zut = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          TT_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zut->Clone("TotalSig"); signal->Add((TH1F*) TT_Zut->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_BDT_toppair_Zutttz.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
          /*
          TFile* inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_BDT_toppair_Zutttz.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_zz_new = new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_BDT_toppair_Zutttz.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_BDT_toppair_Zutttz.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
         // outputfile->cd();
           inputfile->cd();
        /*  ST_Zut->Write();
          signal->Write();
          TT_Zut->Write();*/
           wz_new ->Write();
          zz_new->Write();
          //ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
         // inputfilePDF_ttz_new->Close();
        /*  inputfilePDF_TT->Close();
          inputfilePDF_STNP->Close();
          inputfilePDF_aTT->Close();*/
        }
        

        
      }
      else if(temp.find("MTW")!=std::string::npos){
        //cout << "Get PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_MTW_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_STNP = new TFile("PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
        ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zct->SetTitle((temp+"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
        ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
        
        
        //cout << "Get PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_MTW_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_MTW_toppair_Zutttz.root" ,"READ");
        //cout << "Get PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_MTW_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_MTW_toppair_Zutttz.root" ,"READ");
        
        TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
        TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
        TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        
        TFile* inputfilePDF_STNPu = new TFile("PDFs/root/NP_overlay_ST_FCNC_zut_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        ST_Zut = (TH1F*)(inputfilePDF_STNPu->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zut->SetTitle((temp+"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
        ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_TTu = new TFile("PDFs/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XPDFhistograms_MTW_toppair_Zutttz.root" ,"READ");
        TFile* inputfilePDF_aTTu= new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XPDFhistograms_MTW_toppair_Zutttz.root " ,"READ");
        
        TT_Zut = (TH1F*) inputfilePDF_TTu->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
        TT_Zut->Add( (TH1F*)inputfilePDF_aTTu->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
        TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        wz_new->SetTitle((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        wz_new->SetName((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        cout << (temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str() << endl;
        
        TFile* inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        ttz_new->SetTitle((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        ttz_new->SetName((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        
        TFile* inputfilePDF_zz_new = new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        zz_new->SetTitle((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        zz_new->SetName((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        
        TFile* inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_MTW_toppair_Zutttz.root","READ");
        tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        tzq_new->SetTitle((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        tzq_new->SetName((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        
        
        
        //outputfile->cd();
         inputfile->cd();
        ST_Zct->Write();
         TT_Zct->Write();
        wz_new ->Write();
        zz_new->Write();
        ttz_new->Write();
        tzq_new->Write();
        ST_Zut->Write();
        TT_Zut->Write();
        inputfilePDF_tzq_new ->Close();
        inputfilePDF_WZ->Close();
        inputfilePDF_zz_new->Close();
        inputfilePDF_ttz_new->Close();
        inputfilePDF_TTu->Close();
        inputfilePDF_TT->Close();
        inputfilePDF_STNPu->Close();
        inputfilePDF_STNP->Close();
        inputfilePDF_aTTu->Close();
        inputfilePDF_aTT->Close();

      }
      else if(temp.find("Zmass")!=std::string::npos){
        //cout << "Get PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_STNP ;
        if(!isSingle) inputfilePDF_STNP= new TFile("PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else if(isSingle) inputfilePDF_STNP= new TFile("PDFs/root/NP_overlay_ST_FCNC_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
        ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zct->SetTitle((temp+"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
        ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
        
        
        //cout << "Get PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_TT;
        if(!isSingle) inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root" ,"READ");
        else inputfilePDF_TT = new TFile("PDFs/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_singletop_Zutttz.root" ,"READ");
        //cout << "Get PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root" << endl;
        TFile* inputfilePDF_aTT;
        if(!isSingle) inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_toppair_Zutttz.root" ,"READ");
        else inputfilePDF_aTT = new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XPDFhistograms_Zmass_singletop_Zutttz.root" ,"READ");
        
        TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
        TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
        TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        
        TFile* inputfilePDF_STNPu;
        if(!isSingle) inputfilePDF_STNPu = new TFile("PDFs/root/NP_overlay_ST_FCNC_zut_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else inputfilePDF_STNPu = new TFile("PDFs/root/NP_overlay_ST_FCNC_zut_80XPDFhistograms_Zmass_singletop_Zutttz.root","READ");
        ST_Zut = (TH1F*)(inputfilePDF_STNPu->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zut->SetTitle((temp+"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
        ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_TTu;
        if(!isSingle) inputfilePDF_TTu = new TFile("PDFs/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XPDFhistograms_Zmass_toppair_Zutttz.root" ,"READ");
        else  inputfilePDF_TTu = new TFile("PDFs/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XPDFhistograms_Zmass_singletop_Zutttz.root" ,"READ");
        TFile* inputfilePDF_aTTu;
        if(!isSingle) inputfilePDF_aTTu= new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XPDFhistograms_Zmass_toppair_Zutttz.root " ,"READ");
        else inputfilePDF_aTTu= new TFile("PDFs/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XPDFhistograms_Zmass_singletop_Zutttz.root " ,"READ");
        
        TT_Zut = (TH1F*) inputfilePDF_TTu->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
        TT_Zut->Add( (TH1F*)inputfilePDF_aTTu->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
        TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_WZ;
        if(!isSingle) inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else inputfilePDF_WZ = new TFile("PDFs/root/WZTo3LNu_amc_80XPDFhistograms_Zmass_singletop_Zutttz.root","READ");
        wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        wz_new->SetTitle((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        wz_new->SetName((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_ttz_new ;
        if(!isSingle) inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else inputfilePDF_ttz_new = new TFile("PDFs/root/TTZToLLNuNu_amc_80XPDFhistograms_Zmass_singletop_Zutttz.root","READ");
        ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        ttz_new->SetTitle((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        ttz_new->SetName((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        
        TFile* inputfilePDF_zz_new;
        if(!isSingle) inputfilePDF_zz_new= new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else inputfilePDF_zz_new= new TFile("PDFs/root/ZZTo4L_80XPDFhistograms_Zmass_singletop_Zutttz.root","READ");
        zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        zz_new->SetTitle((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        zz_new->SetName((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        
        TFile* inputfilePDF_tzq_new;
        if(!isSingle) inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_Zmass_toppair_Zutttz.root","READ");
        else  inputfilePDF_tzq_new = new TFile("PDFs/root/tZq_amc_80XPDFhistograms_Zmass_singletop_Zutttz.root","READ");
        tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        tzq_new->SetTitle((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        tzq_new->SetName((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        
        
        
        //outputfile->cd();
        inputfile->cd();
        ST_Zct->Write();
        TT_Zct->Write();
        wz_new ->Write();
        zz_new->Write();
        ttz_new->Write();
        tzq_new->Write();
        ST_Zut->Write();
        TT_Zut->Write();
        
        inputfilePDF_tzq_new ->Close();
        inputfilePDF_WZ->Close();
        inputfilePDF_zz_new->Close();
        inputfilePDF_ttz_new->Close();
        inputfilePDF_TTu->Close();
        inputfilePDF_TT->Close();
        inputfilePDF_STNPu->Close();
        inputfilePDF_STNP->Close();
        inputfilePDF_aTTu->Close();
        inputfilePDF_aTT->Close();
      }
    }
  }



   cout << "Reno/Facts" << endl;
  vector <string> renfactsys = {"RenFactEnvelopeDown", "RenFactEnvelopeUp"};
  for(int isys = 0; isys < 2; isys++){
    string syst = renfactsys[isys];
    
    for(int iCHan = 0; iCHan < channels.size(); iCHan++){
      string channel = channels[iCHan];
      if(temp.find("BDT")!=std::string::npos && region.find("singletop")!=std::string::npos){
        if(coupling.find("Zct")!=std::string::npos){
          
          /*
          //cout << "Get RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_BDT_singletop_Zct.root" << endl;
          TFile* inputfilePDF_STNP = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_BDT_singletop_Zct.root","READ");
          //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
          ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
          ST_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
          
          
          //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_singletop_Zct.root" << endl;
          TFile* inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_singletop_Zct.root" ,"READ");
          //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_singletop_Zct.root" << endl;
          TFile* inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_singletop_Zct.root" ,"READ");
          
          TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          TT_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zct->Clone("TotalSig"); signal->Add((TH1F*) TT_Zct->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          
          TFile* inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_BDT_singletop_Zct.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
          TFile* inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_BDT_singletop_Zct.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ttZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_"+syst).c_str());
          
          TFile* inputfilePDF_zz_new = new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_BDT_singletop_Zct.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ZZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_BDT_singletop_Zct.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("tZq");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          
          //outputfile->cd();
           inputfile->cd();
         // ST_Zct->Write();
         // signal->Write();
         // TT_Zct->Write();
          wz_new->Write();
          zz_new->Write();
          ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
          inputfilePDF_ttz_new->Close();
         // inputfilePDF_TT->Close();
         // inputfilePDF_STNP->Close();
         // inputfilePDF_aTT->Close();
        }
        else if(coupling.find("Zut")!=std::string::npos){
          /*
          
          TFile* inputfilePDF_STNP = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zut_80XRenoFactohistograms_BDT_singletop_Zut.root","READ");
          ST_Zut = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
          ST_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XRenoFactohistograms_BDT_singletop_Zut.root" ,"READ");
          TFile* inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XRenoFactohistograms_BDT_singletop_Zut.root " ,"READ");
          
          TT_Zut = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          TT_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zut->Clone("TotalSig"); signal->Add((TH1F*) TT_Zut->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_BDT_singletop_Zut.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
          TFile* inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_BDT_singletop_Zut.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          
          TFile* inputfilePDF_zz_new = new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_BDT_singletop_Zut.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_BDT_singletop_Zut.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          //outputfile->cd();
           inputfile->cd();
        /*  ST_Zut->Write();
          signal->Write();
          TT_Zut->Write();*/
          wz_new->Write();
          zz_new->Write();
          ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
          inputfilePDF_ttz_new->Close();
        /*  inputfilePDF_TT->Close();
          inputfilePDF_STNP->Close();
          inputfilePDF_aTT->Close();*/
        }
        
        
        
      }
      else if(temp.find("BDT")!=std::string::npos && region.find("toppair")!=std::string::npos){
        if(coupling.find("Zct")!=std::string::npos){
          /*
          
          //cout << "Get RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_BDT_toppair_Zct.root" << endl;
          TFile* inputfilePDF_STNP = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_BDT_toppair_Zct.root","READ");
          //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
          ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
          ST_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
          
          
          //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_toppair_Zct.root" << endl;
          TFile* inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_toppair_Zct.root" ,"READ");
          //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_toppair_Zct.root" << endl;
          TFile* inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_BDT_toppair_Zct.root" ,"READ");
          
          TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
          TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
          TT_Zct->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          TT_Zct->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zct->Clone("TotalSig"); signal->Add((TH1F*) TT_Zct->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          */
          TFile* inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_BDT_toppair_Zct.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_BDT_toppair_Zct.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          
          TFile* inputfilePDF_zz_new = new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_BDT_toppair_Zct.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          TFile* inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_BDT_toppair_Zct.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
          
          //outputfile->cd();
           inputfile->cd();
        //  ST_Zct->Write();
        //signal->Write();
         // TT_Zct->Write();
          wz_new ->Write();
          zz_new->Write();
          ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
          inputfilePDF_ttz_new->Close();
         // inputfilePDF_TT->Close();
         // inputfilePDF_STNP->Close();
         // inputfilePDF_aTT->Close();
        }
        else if(coupling.find("Zut")!=std::string::npos){
          
          
       /*   TFile* inputfilePDF_STNP = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zut_80XRenoFactohistograms_BDT_toppair_Zut.root","READ");
          ST_Zut = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
          //  //cout << "cloned" << endl;
          ST_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
          ST_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XRenoFactohistograms_BDT_toppair_Zut.root" ,"READ");
          TFile* inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XRenoFactohistograms_BDT_toppair_Zut.root " ,"READ");
          
          TT_Zut = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
          TT_Zut->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
          TT_Zut->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          TT_Zut->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
          
          signal = (TH1F*) ST_Zut->Clone("TotalSig"); signal->Add((TH1F*) TT_Zut->Clone(""));
          signal->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
          signal->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TotalSignal_80X_"+syst).c_str());
        */
          TFile* inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_BDT_toppair_Zut.root","READ");
          wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          wz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          wz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_BDT_toppair_Zut.root","READ");
          ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          ttz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          ttz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
          
          TFile* inputfilePDF_zz_new = new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_BDT_toppair_Zut.root","READ");
          zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          zz_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          zz_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
          
          
          TFile* inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_BDT_toppair_Zut.root","READ");
          tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
          tzq_new->SetTitle((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          tzq_new->SetName((coupling+"_"+temp+"_"+ region +"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
          
          
         // outputfile->cd();
           inputfile->cd();
         // ST_Zut->Write();
         // signal->Write();
         // TT_Zut->Write();
          wz_new ->Write();
          zz_new->Write();
          ttz_new->Write();
          tzq_new->Write();
          inputfilePDF_tzq_new ->Close();
          inputfilePDF_WZ->Close();
          inputfilePDF_zz_new->Close();
          inputfilePDF_ttz_new->Close();
         // inputfilePDF_TT->Close();
        //  inputfilePDF_STNP->Close();
       //   inputfilePDF_aTT->Close();
        }
        
       
        
      }
      else if(temp.find("MTW")!=std::string::npos){
        //cout << "Get RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_MTW_toppair_Zut.root" << endl;
        TFile* inputfilePDF_STNP = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
        ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zct->SetTitle((temp+"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
        ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
        
        
        //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_MTW_toppair_Zut.root" << endl;
        TFile* inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_MTW_toppair_Zut.root" ,"READ");
        //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_MTW_toppair_Zut.root" << endl;
        TFile* inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_MTW_toppair_Zut.root" ,"READ");
        
        TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
        TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
        TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        
        TFile* inputfilePDF_STNPu = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zut_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        ST_Zut = (TH1F*)(inputfilePDF_STNPu->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zut->SetTitle((temp+"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
        ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_TTu = new TFile("RenoFact/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XRenoFactohistograms_MTW_toppair_Zut.root" ,"READ");
        TFile* inputfilePDF_aTTu= new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XRenoFactohistograms_MTW_toppair_Zut.root " ,"READ");
        
        TT_Zut = (TH1F*) inputfilePDF_TTu->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
        TT_Zut->Add( (TH1F*)inputfilePDF_aTTu->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
        TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        wz_new->SetTitle((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        wz_new->SetName((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        ttz_new->SetTitle((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        ttz_new->SetName((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        
        TFile* inputfilePDF_zz_new = new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        zz_new->SetTitle((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        zz_new->SetName((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        
        TFile* inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_MTW_toppair_Zut.root","READ");
        tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        tzq_new->SetTitle((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        tzq_new->SetName((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        
        
        
       // outputfile->cd();
         inputfile->cd();
        ST_Zct->Write();
        TT_Zct->Write();
        wz_new ->Write();
        zz_new->Write();
        ttz_new->Write();
        tzq_new->Write();
        ST_Zut->Write();
        TT_Zut->Write();
        
        inputfilePDF_tzq_new ->Close();
        inputfilePDF_WZ->Close();
        inputfilePDF_zz_new->Close();
        inputfilePDF_ttz_new->Close();
        inputfilePDF_TTu->Close();
        inputfilePDF_TT->Close();
        inputfilePDF_STNPu->Close();
        inputfilePDF_STNP->Close();
        inputfilePDF_aTTu->Close();
        inputfilePDF_aTT->Close();
        
      }
      else if(temp.find("Zmass")!=std::string::npos){
        //cout << "Get RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root" << endl;
        TFile* inputfilePDF_STNP ;
        if(!isSingle) inputfilePDF_STNP= new TFile("RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else if(isSingle) inputfilePDF_STNP= new TFile("RenoFact/root/NP_overlay_ST_FCNC_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        //cout << "histo " << ("NP_overlay_ST_FCNC_zct_80X_"+temp +"_" +channel+"_"+syst).c_str() << endl;
        ST_Zct = (TH1F*)(inputfilePDF_STNP->Get(("NP_overlay_ST_FCNC_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zct->SetTitle((temp+"_"+ channel +"_ST_Zct_80X_"+syst).c_str());
        ST_Zct->SetName((temp+"_"+  channel +"_ST_Zct_80X_"+syst).c_str());
        
        
        //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root" << endl;
        TFile* inputfilePDF_TT;
        if(!isSingle) inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root" ,"READ");
        else inputfilePDF_TT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_singletop_Zut.root" ,"READ");
        //cout << "Get RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root" << endl;
        TFile* inputfilePDF_aTT;
        if(!isSingle) inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_toppair_Zut.root" ,"READ");
        else inputfilePDF_aTT = new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80XRenoFactohistograms_Zmass_singletop_Zut.root" ,"READ");
        
        TT_Zct = (TH1F*) inputfilePDF_TT->Get(("NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct");
        TT_Zct->Add( (TH1F*)inputfilePDF_aTT->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zct"));
        TT_Zct->SetTitle((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        TT_Zct->SetName((temp+"_"+  channel +"_TT_Zct_80X_"+syst).c_str());
        
        TFile* inputfilePDF_STNPu;
        if(!isSingle) inputfilePDF_STNPu = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zut_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else inputfilePDF_STNPu = new TFile("RenoFact/root/NP_overlay_ST_FCNC_zut_80XRenoFactohistograms_Zmass_singletop_Zut.root","READ");
        ST_Zut = (TH1F*)(inputfilePDF_STNPu->Get(("NP_overlay_ST_FCNC_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("ST_Zct"));
        //  //cout << "cloned" << endl;
        ST_Zut->SetTitle((temp+"_"+ channel +"_ST_Zut_80X_"+syst).c_str());
        ST_Zut->SetName((temp+"_"+  channel +"_ST_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_TTu;
         if(!isSingle) inputfilePDF_TTu = new TFile("RenoFact/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XRenoFactohistograms_Zmass_toppair_Zut.root" ,"READ");
        else  inputfilePDF_TTu = new TFile("RenoFact/root/NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80XRenoFactohistograms_Zmass_singletop_Zut.root" ,"READ");
        TFile* inputfilePDF_aTTu;
        if(!isSingle) inputfilePDF_aTTu= new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XRenoFactohistograms_Zmass_toppair_Zut.root " ,"READ");
        else inputfilePDF_aTTu= new TFile("RenoFact/root/NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80XRenoFactohistograms_Zmass_singletop_Zut.root " ,"READ");
        
        TT_Zut = (TH1F*) inputfilePDF_TTu->Get(("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut");
        TT_Zut->Add( (TH1F*)inputfilePDF_aTTu->Get(("NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("TT_Zut"));
        TT_Zut->SetTitle((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        TT_Zut->SetName((temp+"_"+  channel +"_TT_Zut_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_WZ;
        if(!isSingle) inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else inputfilePDF_WZ = new TFile("RenoFact/root/WZTo3LNu_amc_80XRenoFactohistograms_Zmass_singletop_Zut.root","READ");
        wz_new = (TH1F*) inputfilePDF_WZ->Get(("WZTo3LNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        wz_new->SetTitle((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        wz_new->SetName((temp+"_"+  channel +"_WZTo3LNu_amc_80X_"+syst).c_str());
        
        
        TFile* inputfilePDF_ttz_new ;
        if(!isSingle) inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else inputfilePDF_ttz_new = new TFile("RenoFact/root/TTZToLLNuNu_amc_80XRenoFactohistograms_Zmass_singletop_Zut.root","READ");
        ttz_new = (TH1F*) inputfilePDF_ttz_new->Get(("TTZToLLNuNu_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        ttz_new->SetTitle((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        ttz_new->SetName((temp+"_"+  channel +"_TTZToLLNuNu_amc_80X_"+syst).c_str());
        
        TFile* inputfilePDF_zz_new;
        if(!isSingle) inputfilePDF_zz_new= new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else inputfilePDF_zz_new= new TFile("RenoFact/root/ZZTo4L_80XRenoFactohistograms_Zmass_singletop_Zut.root","READ");
        zz_new = (TH1F*) inputfilePDF_zz_new->Get(("ZZTo4L_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        zz_new->SetTitle((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        zz_new->SetName((temp+"_"+  channel +"_ZZTo4L_80X_"+syst).c_str());
        
        TFile* inputfilePDF_tzq_new;
        if(!isSingle) inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_Zmass_toppair_Zut.root","READ");
        else  inputfilePDF_tzq_new = new TFile("RenoFact/root/tZq_amc_80XRenoFactohistograms_Zmass_singletop_Zut.root","READ");
        tzq_new = (TH1F*) inputfilePDF_tzq_new->Get(("tZq_amc_80X_"+temp +"_"+channel+"_"+syst).c_str())->Clone("WZ");
        tzq_new->SetTitle((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        tzq_new->SetName((temp+"_"+  channel +"_tZq_amc_80X_"+syst).c_str());
        
        
        
        //outputfile->cd();
         inputfile->cd();
         ST_Zct->Write();
         TT_Zct->Write();
        wz_new ->Write();
        zz_new->Write();
        ttz_new->Write();
        tzq_new->Write();
        ST_Zut->Write();
        TT_Zut->Write();
        
        inputfilePDF_tzq_new ->Close();
        inputfilePDF_WZ->Close();
        inputfilePDF_zz_new->Close();
        inputfilePDF_ttz_new->Close();
        inputfilePDF_TTu->Close();
        inputfilePDF_TT->Close();
        inputfilePDF_STNPu->Close();
        inputfilePDF_STNP->Close();
        inputfilePDF_aTTu->Close();
        inputfilePDF_aTT->Close();
        
      }
    }
  }
  inputfile->Write();

  inputfile->Close();

  
  


}
