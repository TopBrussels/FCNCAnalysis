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





void tables(int mode = 0){
  
  double lumi = 35822.43741;
  bool dolandscape = true;
  
  
  char myTexFile[300];
  sprintf(myTexFile,"tables/table_%d_%fpb.tex", mode, lumi);
  ofstream salida(myTexFile);
  
  char myRootFile[300];
  sprintf(myRootFile,"cutflowtables.root");
  TFile *_file0 = TFile::Open(myRootFile);
  
  const int nProcess = 9;
  TString processName[nProcess] = { "FCNCZut", "FCNCZct","fake","ttZ", "WZ", "tZq", "other", "data", "mc"};
  TString processLabel[nProcess] = { "\\textbf{FCNC tZu}","\\textbf{non prompt lepton}", "\\textbf{FCNC tZc}", "\\textbf{\\ttZ}", "\\textbf{\\WZ+jets}", "\\textbf{\\tZq}","\\textbf{$Other$}", "\\textbf{data}", "\\textbf{$MC$ (no fakes)}"};
  const int nCuts = 9;
  
  TString cutLabel[10] =  {"blank","$>1l,>0j, <6j, m_{T}^{W} < 300$", "SF pair","lep veto","Z mass","$>2l$","$\\Delta R (l_{W},b) \\leq 2.5$","\\STSR","\\TTSR","\\WZCR"};
  
  TH1F*  h[nProcess];
  for(int i=0; i<nProcess; i++){
    cout << "adding " << processName[i] << endl;
    if (i == nProcess-1) {
     h[i] =  (TH1F*)h[4]->Clone();
    // h[i]->Add(h[1]);
     //h[i]->Add(h[2]);
     //h[i]->Add(h[3]);
     //h[i]->Add(h[3]);
     //h[i]->Add(h[5]);
     //h[i]->Add(h[6]);
     //h[i]->Add(h[7]);
     }
    else{
     
    if(mode == 0)  h[i] = (TH1F*) _file0->Get("CutflowTableHisto_uuu__"+processName[i]);
    else if(mode == 1)  h[i] = (TH1F*) _file0->Get("CutflowTableHisto_uue__"+processName[i]);
    else if(mode == 2)  h[i] = (TH1F*) _file0->Get("CutflowTableHisto_eeu__"+processName[i]);
    else if(mode == 3)  h[i] = (TH1F*) _file0->Get("CutflowTableHisto_eee__"+processName[i]);
    else if(mode == 9)  h[i] = (TH1F*) _file0->Get("CutflowTableHisto__"+processName[i]);
    }
    // Lepton ID and HLT SF
    // if (mode == 0 && i != 8) h[i]->Scale(0.97713);
    // if (mode == 1 && i != 8) h[i]->Scale(0.910067);
    // if (mode == 2 && i != 8) h[i]->Scale(0.945736);
  }
  
  cout << "filling vector " << endl;
  double vectorValue[nProcess][nCuts][3];
  
  
  for (int i = 1; i < nCuts+1; i++){
    for (int j = 0; j < nProcess; j++){
      cout << h[j]->GetBinContent(i) << endl;
      cout << "looking at cut " << cutLabel[i] << " for process " << processName[j] << endl;
      vectorValue[j][i][0] = h[j]->GetBinContent(i);
      vectorValue[j][i][1] = 1; //precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i);
      cout << vectorValue[j][i][0] << " pm " << vectorValue[j][i][2] << endl;
    }
  }
  cout << "write table" << endl;
  
  salida << "\\documentclass[a4paper,12pt]{article}" << endl;
  salida << "\\begin{document}" << endl;
  salida << endl;
  salida << endl;
  
  if(dolandscape)  salida << "  \\begin{landscapce}" << endl;
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < nProcess; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int i=1; i  < nCuts+1; i++){
    salida <<  cutLabel[i] ;
    for (int j = 0; j < nProcess; j++){
      /*if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
        salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
      } else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){*/
      if (i != 0 && vectorValue[j][i][0] == 0 ){
        salida << " & $-$ " ;
      } else {
        salida << " & " << std::setiosflags(std::ios::fixed) << setprecision(vectorValue[j][i][1]) << vectorValue[j][i][0] ;
        salida << " $\\pm $"  << setprecision(vectorValue[j][i][1])<< vectorValue[j][i][2];
      }
    }
    salida <<  " \\\\  " << endl;
  }
  salida << "   \\hline " << endl;
  salida << "  \\end{tabular}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  if(dolandscape)  salida << "  \\end{landscapce}" << endl;
  salida << endl;
  salida << endl;
  

  salida << "\\end{document}" << endl;
  
}


