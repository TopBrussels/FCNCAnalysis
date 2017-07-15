//Isis Van Parijs

#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
//#include "inputs.h"


using namespace std;



void datacards(TString  coupling = "Zct", TString channel = "3mu", TString region = "ST" ){
  TString regions[3] = {"ST","TT","WZ"};
  TString couplings[2] = {"Zct","Zut"};
  TString channels[4] = {"3mu","3e","2e1mu","1e2mu"};
  
  for(int iRegion = 0 ; iRegion < 3; iRegion++)
  {
    
    region = regions[iRegion];
    for(int iCoupling = 0; iCoupling < 2; iCoupling++){
      
      coupling = couplings[iCoupling];
      if(coupling.Contains("Zct") && region.Contains("WZ")) continue;
      for( int iChannel =0 ; iChannel < 4 ; iChannel++){
        channel = channels[iChannel];
        
        string inputfilename = "./datacards/rootfiles/";
        TString outputfilename = "./datacardcreator/datacard_"+channel+"_"+region+"_" +coupling+".txt";
        cout << "Writing to " << outputfilename << endl;
        int modes = 0;
        TString mode = region + coupling + "_" + channel;
        if(region.Contains("WZ")) mode = region + "_" +channel;
        
        if(mode.Contains("TTZct")){ inputfilename += "Zct_toppair/Reader_Zct_toppair.root";}
        else if(mode.Contains("TTZut")){ inputfilename +=  "Zut_toppair/Reader_Zut_toppair.root";}
        else if(mode.Contains("STZct")){ inputfilename +=  "Zct_singletop/Reader_Zct_singletop.root";}
        else if(mode.Contains("STZut")){ inputfilename += "Zut_singletop/Reader_Zut_singletop.root";}
        else if(mode.Contains("WZ")){     inputfilename += "MTWtemplate/Reader_Zut_MTW.root";   }
        
        ofstream salida(outputfilename.Data());
        
        
        TFile *_file0 = TFile::Open(inputfilename.c_str());
        
        const int nProcess = 17;
        TString processName[nProcess] = {"FakeEl","NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct","NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct","NP_overlay_ST_FCNC_zct","WZTo3LNu_amc","tZq_amc","tHq","TTWJetsToLNu_amc","TTZToLLNuNu_amc","ttHToNonbb","ttHTobb","ZZTo4L","WZZ_amc","ZZZ_amc","tWll","STtW_top","STtW_atop"};
        int processNb[nProcess] = {1,0,-1,-2,2,3,4,5,6,7,8,9,10,10,11,11,11};
        if(channel.Contains("3e") || channel.Contains("1e2mu")) processName[0] = "FakeEl";
        else if (channel.Contains("3mu") || channel.Contains("2e1mu")) processName[0] = "FakeMu";
        
        if(coupling.Contains("Zct")){
          processName[1] = "NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct";
          processName[2] = "NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct";
          processName[3] = "NP_overlay_ST_FCNC_zct";
        }
        else if(coupling.Contains("Zut")){
          processName[1] = "NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut";
          processName[2] = "NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut";
          processName[3] = "NP_overlay_ST_FCNC_zut";
        }
        
        const int nShapeSystematics = 13;
        TString ShapeSystematics[nShapeSystematics] =  {"puSF","electronSF","muonSF","btagSF_cferr1","btagSF_cferr2","btagSF_hf","btagSF_hfstats1","btagSF_hfstats2","btagSF_lf","btagSF_lfstats1","btagSF_lfstats2","JES","JER"};
        const int nRateSystematics = 9;
        TString RateSystematics[nRateSystematics] ={"Lumi","TrigRate","WZrate","ttZrate","tZqrate", "STrate", "ttHrate","ZZrate","FakeRate"};
        double RateSystematicsNb[nRateSystematics] ={1.025,1.01,1.3,1.3,1.3,1.3,1.3,1.3,1.5};
        if(mode.Contains("3mu") || mode.Contains("2e1mu")) RateSystematicsNb[1] = 1.01;
        else if(mode.Contains("3e") || mode.Contains("1e2mu")) RateSystematicsNb[1] = 1.05;
        
        salida << "imax  * number of bins" << endl;
        salida << "jmax  * number of processes minus one" << endl;
        salida << "kmax  * number of nuisance parameters" << endl;
        salida << "--------------------------------------------------------------------------------" << endl;
        salida << "shapes *         " << mode << "  " << "../" << inputfilename ;
        
        if(mode.Contains("WZ_3mu")) salida << "       MTW_uuu_$PROCESS_80X  MTW_uuu_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("WZ_3e")) salida << "       MTW_eee_$PROCESS_80X  MTW_eee_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("WZ_1e2mu")) salida << "       MTW_uue_$PROCESS_80X  MTW_uue_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("WZ_2e1mu")) salida << "       MTW_eeu_$PROCESS_80X  MTW_eeu_$PROCESS_80X_$SYSTEMATIC" << endl;
        
        else if(mode.Contains("STZct_3mu")) salida << "       Zct_BDT_singletop_uuu_$PROCESS_80X  Zct_BDT_singletop_uuu_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZct_3e")) salida << "       Zct_BDT_singletop_eee_$PROCESS_80X  Zct_BDT_singletop_eee_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZct_1e2mu")) salida << "       Zct_BDT_singletop_uue_$PROCESS_80X  Zct_BDT_singletop_uue_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZct_2e1mu")) salida << "       Zct_BDT_singletop_eeu_$PROCESS_80X  Zct_BDT_singletop_eeu_$PROCESS_80X_$SYSTEMATIC" << endl;
        
        else if(mode.Contains("STZut_3mu")) salida << "       Zut_BDT_singletop_uuu_$PROCESS_80X  Zut_BDT_singletop_uuu_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZut_3e")) salida << "       Zut_BDT_singletop_eee_$PROCESS_80X  Zut_BDT_singletop_eee_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZut_1e2mu")) salida << "       Zut_BDT_singletop_uue_$PROCESS_80X  Zut_BDT_singletop_uue_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("STZut_2e1mu")) salida << "       Zut_BDT_singletop_eeu_$PROCESS_80X  Zut_BDT_singletop_eeu_$PROCESS_80X_$SYSTEMATIC" << endl;
        
        else if(mode.Contains("TTZct_3mu")) salida << "       Zct_BDT_toppair_uuu_$PROCESS_80X  Zct_BDT_toppair_uuu_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZct_3e")) salida << "       Zct_BDT_toppair_eee_$PROCESS_80X  Zct_BDT_toppair_eee_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZct_1e2mu")) salida << "       Zct_BDT_toppair_uue_$PROCESS_80X  Zct_BDT_toppair_uue_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZct_2e1mu")) salida << "       Zct_BDT_toppair_eeu_$PROCESS_80X  Zct_BDT_toppair_eeu_$PROCESS_80X_$SYSTEMATIC" << endl;
        
        else if(mode.Contains("TTZut_3mu")) salida << "       Zut_BDT_toppair_uuu_$PROCESS_80X  Zut_BDT_toppair_uuu_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZut_3e")) salida << "       Zut_BDT_toppair_eee_$PROCESS_80X  Zut_BDT_toppair_eee_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZut_1e2mu")) salida << "       Zut_BDT_toppair_uue_$PROCESS_80X  Zut_BDT_toppair_uue_$PROCESS_80X_$SYSTEMATIC" << endl;
        else if(mode.Contains("TTZut_2e1mu")) salida << "       Zut_BDT_toppair_eeu_$PROCESS_80X  Zut_BDT_toppair_eeu_$PROCESS_80X_$SYSTEMATIC" << endl;
        
        salida << "shapes data_obs  " << mode << "  " << inputfilename ;
        if(mode.Contains("WZ_3mu")) salida << "       MTW_uuu_data" << endl;
        else if(mode.Contains("WZ_3e")) salida << "       MTW_eee_data" << endl;
        else if(mode.Contains("WZ_1e2mu")) salida << "       MTW_uue_data" << endl;
        else if(mode.Contains("WZ_2e1mu")) salida << "       MTW_eeu_data" << endl;
        
        else if(mode.Contains("STZct_3mu")) salida << "       Zct_BDT_singletop_uuu_data" << endl;
        else if(mode.Contains("STZct_3e")) salida << "       Zct_BDT_singletop_eee_data" << endl;
        else if(mode.Contains("STZct_1e2mu")) salida << "       Zct_BDT_singletop_uue_data" << endl;
        else if(mode.Contains("STZct_2e1mu")) salida << "       Zct_BDT_singletop_eeu_data" << endl;
        
        else if(mode.Contains("STZut_3mu")) salida << "       Zut_BDT_singletop_uuu_data" << endl;
        else if(mode.Contains("STZut_3e")) salida << "       Zut_BDT_singletop_eee_data" << endl;
        else if(mode.Contains("STZut_1e2mu")) salida << "       Zut_BDT_singletop_uue_data" << endl;
        else if(mode.Contains("STZut_2e1mu")) salida << "       Zut_BDT_singletop_eeu_data" << endl;
        
        else if(mode.Contains("TTZct_3mu")) salida << "       Zct_BDT_toppair_uuu_data" << endl;
        else if(mode.Contains("TTZct_3e")) salida << "       Zct_BDT_toppair_eee_data" << endl;
        else if(mode.Contains("TTZct_1e2mu")) salida << "       Zct_BDT_toppair_uue_data" << endl;
        else if(mode.Contains("TTZct_2e1mu")) salida << "       Zct_BDT_toppair_eeu_data" << endl;
        
        else if(mode.Contains("TTZut_3mu")) salida << "       Zut_BDT_toppair_uuu_data" << endl;
        else if(mode.Contains("TTZut_3e")) salida << "       Zut_BDT_toppair_eee_data" << endl;
        else if(mode.Contains("TTZut_1e2mu")) salida << "       Zut_BDT_toppair_uue_data" << endl;
        else if(mode.Contains("TTZut_2e1mu")) salida << "       Zut_BDT_toppair_eeu_data" << endl;
        salida << "--------------------------------------------------------------------------------" << endl;
        
        
        vector <int> ProcessToSkip;
        int NbOfProcessToFill = 0;
        cout << "Checking of empty histograms " << endl ;
        TH1F* tempHisto;
        for(int iProcess = 0 ; iProcess < nProcess; iProcess++){
          
          if(mode.Contains("STZct_3mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_singletop_uuu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("STZut_3mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_singletop_uuu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZct_3mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_toppair_uuu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZut_3mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_toppair_uuu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("WZ_3mu")) tempHisto = (TH1F*) _file0->Get("MTW_uuu_"+processName[iProcess]+"_80X");
          
          if(mode.Contains("STZct_3e")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_singletop_eee_"+processName[iProcess]+"_80X");
          else if(mode.Contains("STZut_3e")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_singletop_eee_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZct_3e")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_toppair_eee_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZut_3e")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_toppair_eee_"+processName[iProcess]+"_80X");
          else if(mode.Contains("WZ_3e")) tempHisto = (TH1F*) _file0->Get("MTW_eee_"+processName[iProcess]+"_80X");
          
          
          if(mode.Contains("STZct_2e1mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_singletop_eeu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("STZut_2e1mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_singletop_eeu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZct_2e1mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_toppair_eeu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZut_2e1mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_toppair_eeu_"+processName[iProcess]+"_80X");
          else if(mode.Contains("WZ_2e1mu")) tempHisto = (TH1F*) _file0->Get("MTW_eeu_"+processName[iProcess]+"_80X");
          
          
          if(mode.Contains("STZct_1e2mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_singletop_uue_"+processName[iProcess]+"_80X");
          else if(mode.Contains("STZut_1e2mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_singletop_uue_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZct_1e2mu")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_toppair_uue_"+processName[iProcess]+"_80X");
          else if(mode.Contains("TTZut_1e2mu")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_toppair_uue_"+processName[iProcess]+"_80X");
          else if(mode.Contains("WZ_1e2mu")) tempHisto = (TH1F*) _file0->Get("MTW_uue_"+processName[iProcess]+"_80X");
          
          
          
          if(tempHisto){          cout << "* looking at " << tempHisto->GetTitle() << endl;
          if(tempHisto->GetEntries() ==0){
            ProcessToSkip.push_back(iProcess);
            cout <<"  empty " << endl;
          }
          else NbOfProcessToFill++;
          }
          else {
            cout << "* looking at " << processName[iProcess]<< " (does not exist)" << endl;
            ProcessToSkip.push_back(iProcess);
            cout <<"  empty " << endl;
          }
        }
        
        cout << "  --> Found " << ProcessToSkip.size() << " empty histogram(s)" << endl;
        
        
        
        salida << "bin          " << mode << endl ;
        salida << "observation  " << "-1" << endl;
        salida << "--------------------------------------------------------------------------------" << endl;
        salida << "bin          " ;
        
        bool SkipProcess = false;
        for(int iProcess =0; iProcess < nProcess ; iProcess++)
        {
          SkipProcess = false;
          for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
            if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
          }
          if(!SkipProcess){
            salida << mode << "           " ;
          }
        }
        salida << endl;
        salida << "process      " ;
        
        for(int iProcess =0; iProcess < nProcess ; iProcess++)
        {
          SkipProcess = false;
          for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
            if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
          }
          if(!SkipProcess){
            salida << processName[iProcess] << "           " ;
          }
        }
        salida << endl;
        salida << "process      " ;
        for(int iProcess =0; iProcess < nProcess ; iProcess++)
        {
          SkipProcess = false;
          for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
            if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
          }
          if(!SkipProcess){
            salida << processNb[iProcess] << "           " ;
          }
        }
        salida << endl;
        salida << "rate         " ;
        for(int iProcess =0; iProcess < nProcess ; iProcess++)
        {
          SkipProcess = false;
          for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
            if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
          }
          if(!SkipProcess){
            salida <<  "-1           " ;
          }
        }
        salida << endl;
        salida << endl;
        
        
        for(int iShapeSys =0; iShapeSys < nShapeSystematics ; iShapeSys++){
          salida << ShapeSystematics[iShapeSys]  <<   "       shape          " ;
          for(int iProcess =0; iProcess < nProcess ; iProcess++)
          {
            SkipProcess = false;
            for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
              if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
            }
            if(!SkipProcess && processName[iProcess].Contains("Fake")){
              salida << "-" << "           " ;
            }
            else if(!SkipProcess && processName[iProcess].Contains("STtW_atop") && ShapeSystematics[iShapeSys].Contains("JES") && region.Contains("ST")  && channel.Contains("1e2mu")){
              salida << "-" << "           " ;
            }
            else if(!SkipProcess){
              salida << 1 << "           " ;
            }
          }
          salida << endl;
        }
        for(int iRateSys =0; iRateSys < nRateSystematics ; iRateSys++){
          salida << RateSystematics[iRateSys]  <<   "       lnN          " ;
          for(int iProcess =0; iProcess < nProcess ; iProcess++)
          {
            SkipProcess = false;
            for(int iSkip = 0; iSkip < ProcessToSkip.size(); iSkip++){
              if(iProcess == ProcessToSkip[iSkip]){ SkipProcess = true; break;}
            }
            if(!SkipProcess && processName[iProcess].Contains("Fake")){
              salida << "-" << "           " ;
            }
            else if(!SkipProcess && iRateSys < 2){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("WZ") && processName[iProcess].Contains("WZTo3") && mode.Contains("WZ")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("Fake") && processName[iProcess].Contains("Fake") && mode.Contains("WZ")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("ttZ") && processName[iProcess].Contains("TTZ")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("ZZ") && processName[iProcess].Contains("ZZT")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("ttH") && processName[iProcess].Contains("ttH")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("ST") && (processName[iProcess].Contains("STtW") || processName[iProcess].Contains("tWll") )){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(RateSystematics[iRateSys].Contains("tZq") && processName[iProcess].Contains("tZq")){
              salida << RateSystematicsNb[iRateSys] << "           " ;
            }
            else if(!SkipProcess){
              salida << "-" << "           " ;
            }
          }
          salida << endl;
        }
        salida << endl;
        salida << "--------------------------------------------------------------------------------" << endl;
        if(mode.Contains("3mu") || mode.Contains("2e1mu") ) salida << "RateParamMu   rateParam    " << mode << "    " << processName[0] <<  "     1    " <<  endl;
        if(mode.Contains("3e") || mode.Contains("1e2mu") ) salida << "RateParamEl   rateParam    " << mode << "    " << processName[0] <<  "     1    " << endl;
        salida << "RateParamWZ   rateParam    " << mode << "    " << processName[4] <<  "     1    " << endl;
        
      }
    }
  }
}


