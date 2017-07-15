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
double normalization(double nevents, double xsec, double lumi){
  
  if (nevents !=0) return lumi*xsec/nevents;
  else return 1;
  
}

double precision(double error){
  
  int precisionValue;
  double factErr = 0;
  int iN = 0;
  if (error == 0 || error >= 1) precisionValue = 0;
  else if (error < 1) {
    iN = 0;
    factErr = 0;
    /*while (factErr < 1){
      //factErr = error*(10**iN);
      iN++;
    }*/
    precisionValue = iN-1;
    precisionValue = 0;
  }
  
  if (factErr > 9.5) precisionValue-=1;
  
  return precisionValue;
  
}

double errorEfi (double efi, double errTotal, double totalEvents, double error, double number){
  
  if (totalEvents !=0 && number != 0) return efi*((errTotal/totalEvents)+(error/number));
  else return 0;
  
}

double efficiency (double finalevents, double totalevents){
  
  if (totalevents !=0) return finalevents*100/totalevents;
  else return 0;
  
}





void datacards(TString  coupling = "Zct", TString channel = "3mu", TString region = "ST" ){
  
  double lumi = 35822.43741;
  bool dolandscape = false;
  string inputfilename = "";
  TString outputfilename = "datacard_"+channel+"_"+region+"_" +coupling+".txt";
  cout << "Writing to " << outputfilename << endl;
  int modes = 0;
  TString mode = region + coupling + "_" + channel;
  if(region.Contains("WZ")) mode = region + "_" +channel;
  if(mode.Contains("TTZct")){ inputfilename += "Reader_Zct_toppair.root";}
  else if(mode.Contains("TTZut")){ inputfilename +=  "Reader_Zut_toppair.root";}
  else if(mode.Contains("STZct")){ inputfilename +=  "Reader_Zct_singletop.root";}
  else if(mode.Contains("STZut")){ inputfilename += "Reader_Zut_singletop.root";}
  else if(mode.Contains("WZ")){     inputfilename += "Reader_Zut_MTW.root";   }

  ofstream salida(outputfilename.Data());

  
  TFile *_file0 = TFile::Open(inputfilename.c_str());
  
  const int nProcess = 17;
  TString processName[nProcess] = {"FakeEl","NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct","NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct","NP_overlay_ST_FCNC_zct","WZTo3LNu_amc","tZq_amc","tHq","TTWJetsToLNu_amc","TTZToLLNuNu_amc","ttHToNonbb","ttHTobb","ZZTo4L","WZZ_amc","ZZZ_amc","tWll","STtW_top","STtW_atop"};
  int processNb[nProcess] = {1,0,-1,-2,2,3,4,5,6,7,8,9,10,11,12,13,14};
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
  const int nRateSystematics = 8;
  TString RateSystematics[nRateSystematics] ={"Lumi","TrigRate","WZrate","ttZrate","tZqrate", "STrate", "ttHrate","ZZrate"};
  double RateSystematicsNb[nRateSystematics] ={1.025,1.01,1.3,1.3,1.3,1.3,1.3,1.3};
  if(mode.Contains("3mu") || mode.Contains("2e1mu")) RateSystematicsNb[1] = 1.01;
  else if(mode.Contains("3e") || mode.Contains("1e2mu")) RateSystematicsNb[1] = 1.05;
  
  salida << "imax  * number of bins" << endl;
  salida << "jmax  * number of processes minus one" << endl;
  salida << "kmax  * number of nuisance parameters" << endl;
  salida << "--------------------------------------------------------------------------------" << endl;
  salida << "shapes *         " << mode << "  " << inputfilename ;
  
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
    else if(mode.Contains("STZut")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_singletop_uuu_"+processName[iProcess]+"_80X");
    else if(mode.Contains("TTZct")) tempHisto = (TH1F*) _file0->Get("Zct_BDT_toppair_uuu_"+processName[iProcess]+"_80X");
    else if(mode.Contains("TTZut")) tempHisto = (TH1F*) _file0->Get("Zut_BDT_toppair_uuu_"+processName[iProcess]+"_80X");
    
    cout << "* looking at " << tempHisto->GetTitle() << endl;
    if(tempHisto->GetEntries() ==0){
      ProcessToSkip.push_back(iProcess);
      cout <<"  empty " << endl;
    }
    else NbOfProcessToFill++;
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
      if( processName[iProcess].Contains("Fake")){
        salida << "-" << "           " ;
      }
      if(!SkipProcess){
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
  if(mode.Contains("3mu") || mode.Contains("2e1mu") ) salida << "RateParamMu   rateParam    " << mode << "    " << processName[0] << endl;
  if(mode.Contains("3e") || mode.Contains("1e2mu") ) salida << "RateParamEl   rateParam    " << mode << "    " << processName[0] << endl;
  salida << "RateParamWZ   rateParam    " << mode << "    " << processName[4] << endl;
  
  /*
  TH1F*  h[nCuts][3][2];
  for(int i=1; i<nCuts+1; i++){
    cout << "adding " << cutLabel[i] << endl;
    for(int j = 0 ; j < 3 ; j++){
     h[i][j][0] = (TH1F*) _file0->Get("hist_BDT_"+cutLabel[i] + "_" + values[j]+"_sig");
     h[i][j][1] = (TH1F*) _file0->Get("hist_BDT_"+cutLabel[i] + "_" + values[j]+"_bkg");
    }
  }
  
  cout << "filling vector " << endl;
  double vectorValue[nCuts][3][2][3]; // sys nom process entries
  
  for (int i = 1; i < nCuts+1; i++){
    for(int k =0; k < 3; k++){
      for (int j = 0; j < nProcess; j++){
        //cout << "i " << i << " j " << j << " k " << k << endl;
        //cout << h[i][k][j]->Integral() << endl;
       // cout << "looking at cut " << cutLabel[i] << " for process " << processName[j] << endl;
        vectorValue[i][k][j][0] = h[i][k][j]->Integral();
        vectorValue[i][k][j][1] = 1; //precision(h[j]->GetBinError(i));
        vectorValue[i][k][j][2] = sqrt(h[i][k][j]->Integral());
        //cout << vectorValue[i][k][j][0] << " pm " <<vectorValue[i][k][j][2]<< endl;
      }
    }
  }
  
  
  
  cout << "write table" << endl;
  
  salida << "\\documentclass[a4paper,12pt]{article}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  if(dolandscape)  salida << "  \\begin{landscape}" << endl;
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "   \\caption{The effect of the systematic uncertainties on the event yield in the BDT for the " ;
  if(mode.Contains("TTZct")) salida << " \\TTSR involving the tZc vertex" ;
  else if(mode.Contains("TTZut")) salida << " \\TTSR involving the tZu vertex" ;
  else  if(mode.Contains("STZct")) salida << " \\STSR  involving the tZc vertex" ;
  else if(mode.Contains("STZut")) salida << " \\STSR involving the tZu vertex" ;
  salida << ". The nominal event yield (nom), and the plus one sigma (up) and minus one sigma (down) yields are given for each systematic. Non prompt lepton backgrounds are not included in the yields.}" << endl;
  salida << "  \\begin{tabular} {l";
  for(int i = 0; i < 2 ; i++)
      {
  salida << "|c" ;
  }
  salida << "|c}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i <3; i++){
    cout << "values name " << values[i] << endl;
    salida << " & ";
    salida << values[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int i=1; i  < nCuts+1; i++){
    salida << "\\multicolumn{" << nval+1 << "}{c}{" << cutLabelName[i] << "} \\\\ \\hline " << endl;
    for(int k=0; k < nProcess; k++){
    salida   << processLabel[k] ;
     for (int j = 0; j < 3; j++){ // values  vectorValue[nCuts][3][2][3]; // sys nom process entries
      if (vectorValue[i][j][k][0] == 0 ){
        salida << " & $-$ " ;
      } else {
        salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[i][j][k][1]) << vectorValue[i][j][k][0] ;
        salida << " $\\pm $"  << setprecision(vectorValue[i][j][k][1])<< vectorValue[i][j][k][2];
      }
    }
    salida <<  " \\\\  " << endl;
    }
    if(i!= nCuts) salida << "\\hline" << endl;
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  if(dolandscape)  salida << "  \\end{landscape}" << endl;
  salida << endl;
  salida << endl;
  

  salida << "\\end{document}" << endl;
  */
}


