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





void tableseff_jesjer_sys(TString mode = "STZct" ){
  
  double lumi = 35822.43741;
  bool dolandscape = false;
  string sysratiosfilename = "sysratio_";
  int modes = 0;
  if(mode.Contains("TTZct")){ modes = 0; sysratiosfilename += "TTZct";}
  else if(mode.Contains("TTZut")){ modes = 1; sysratiosfilename += "TTZut";}
  else if(mode.Contains("STZct")){ modes = 2; sysratiosfilename += "STZct";}
  else if(mode.Contains("STZut")){ modes = 3; sysratiosfilename += "STZut";}
  sysratiosfilename += "_jerjes.root";
  
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/systableeffjesjer_%d_%fpb.tex", modes, lumi);
  ofstream salida(myTexFile);
  
  
  
  TFile *_file0 = TFile::Open(sysratiosfilename.c_str());
  
  const int nProcess = 2;
  TString processName[nProcess] = { "sig", "bkg"};
  TString processLabel[nProcess] = { "Signal","Background"};
  const int nCuts = 2;
  
  const int nval = 3;
  TString cutLabel[3] =  {"blank","JES","JER"};
  TString cutLabelName[3] =  {"blank","Jet Energy Scale","Jet Energy Resolution"};
  
  TString values[3] = {"nom","up","down"};
 // TString valuesname[3] = {"Nominal","sig","sig"}; //"Nominal","$+1\\sigma$","$-1\\sigma$"};
  TString valuesname[nval] = { "\\textbf{Signal}","\\textbf{Background}","\\textbf{Background}"};
  TH1F*  h[nCuts][3][2];
  for(int i=1; i<nCuts+1; i++){
    cout << "adding " << cutLabel[i] << endl;
    for(int j = 0 ; j < 3 ; j++){
     h[i][j][0] = (TH1F*) _file0->Get("hist_BDT_"+cutLabel[i] + "_" + values[j]+"_sig");
     h[i][j][1] = (TH1F*) _file0->Get("hist_BDT_"+cutLabel[i] + "_" + values[j]+"_bkg");
    }
  }
  
  cout << "filling vector " << endl;
  double pvectorValue[nCuts][3][2][3]; // sys nom process entries
  
  for (int i = 1; i < nCuts+1; i++){
    for(int k =0; k < 3; k++){
      for (int j = 0; j < nProcess; j++){
        //cout << "i " << i << " j " << j << " k " << k << endl;
        //cout << h[i][k][j]->Integral() << endl;
       // cout << "looking at cut " << cutLabel[i] << " for process " << processName[j] << endl;
        pvectorValue[i][k][j][0] = h[i][k][j]->Integral();
        pvectorValue[i][k][j][1] = 1; //precision(h[j]->GetBinError(i));
        pvectorValue[i][k][j][2] = sqrt(h[i][k][j]->Integral());
        //cout << vectorValue[i][k][j][0] << " pm " <<vectorValue[i][k][j][2]<< endl;
      }
    }
  }
  
  double vectorValue[nCuts][3][2][3]; // sys nom process entries
  
  for (int i = 1; i < nCuts+1; i++){
    for(int k =0; k < 3; k++){
      for (int j = 0; j < nProcess; j++){
        //cout << "i " << i << " j " << j << " k " << k << endl;
        //cout << h[i][k][j]->Integral() << endl;
        // cout << "looking at cut " << cutLabel[i] << " for process " << processName[j] << endl;
        vectorValue[i][k][j][0] = ((pvectorValue[i][k][j][0] - pvectorValue[i][0][j][0]) /(pvectorValue[i][0][j][0]))*100;
        vectorValue[i][k][j][1] = 1; //precision(h[j]->GetBinError(i));
       // vectorValue[i][k][j][2] = sqrt(h[i][k][j]->Integral());
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
  salida << "   \\caption{The relative effect of the systematic uncertainties on the event yield in the BDT for the " ;
  if(mode.Contains("TTZct")) salida << " \\TTSR involving the tZc vertex" ;
  else if(mode.Contains("TTZut")) salida << " \\TTSR involving the tZu vertex" ;
  else  if(mode.Contains("STZct")) salida << " \\STSR  involving the tZc vertex" ;
  else if(mode.Contains("STZut")) salida << " \\STSR involving the tZu vertex" ;
  salida << ". The relative yields for plus one sigma (up) and minus one sigma (down) are given for each systematic. Non prompt lepton backgrounds are not included in the yields. The CSVv2 yields are obtained by quadratically adding the different contributions.}" << endl;
  salida << "  \\begin{tabular} {l";
  for(int i = 1; i < 2 ; i++)
      {
  salida << "|c" ;
  }
  salida << "|c}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 1; i <3; i++){
    cout << "values name " << values[i] << endl;
    salida << " & ";
    salida << values[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int i=1; i  < nCuts+1; i++){ // no btag
    salida << "\\multicolumn{" << nval << "}{c}{" << cutLabelName[i] << "} \\\\ \\hline " << endl;
    for(int k=0; k < nProcess; k++){
    salida   << processLabel[k] ;
     for (int j = 1; j < 3; j++){ // values  vectorValue[nCuts][3][2][3]; // sys nom process entries
      if (vectorValue[i][j][k][0] == 0 ){
        salida << " & $-$ " ;
      } else {
        salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[i][j][k][1]) << vectorValue[i][j][k][0] << "\\%";
       // salida << " $\\pm $"  << setprecision(vectorValue[i][j][k][1])<< vectorValue[i][j][k][2];
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
  
}


