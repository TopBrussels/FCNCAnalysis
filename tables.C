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





void tables(int mode = 9){
  
  double lumi = 35822.43741;
  bool dolandscape = true;
  
   TH1::SetDefaultSumw2();
  char myTexFile[300];
  sprintf(myTexFile,"tables/table_%d_%fpb.tex", mode, lumi);
  ofstream salida(myTexFile);
  
  char myRootFile[300];
  sprintf(myRootFile,"cutflowtables.root");
  TFile *_file0 = TFile::Open(myRootFile);
  
  TString channelname = "";
  if(mode == 0)  channelname = "uuu_";
  else if(mode == 1)  channelname = "uue_";
  else if(mode == 2)  channelname = "eeu_";
  else if(mode == 3)  channelname = "eee_";
  else if(mode == 9)  channelname = "";
  const int nProcess = 23;
  TString processName[24] {"NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut_80X","NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zut_80X","NP_overlay_ST_FCNC_zut_80X", "NP_overlay_TT_FCNC-T2ZJ_aTleptonic_ZToll_kappa_zct_80X","NP_overlay_TT_FCNC-aT2ZJ_Tleptonic_ZToll_kappa_zct_80X","NP_overlay_ST_FCNC_zct_80X","DYJets50_80X","TTJets_powheg_80X","fake_80X","WZTo3LNu_amc_80X","tZq_amc_80X","TTZToLLNuNu_amc_80X","ZZTo4L_80X",
      "STtW_atop_80X","STtW_top_80X","tWll_80X","ZZZ_amc_80X","WZZ_amc_80X", "TTWJetsToLNu_amc_80X", "tHq_80X", "ttHToNonbb_80X", "ttHTobb_80X","data"};
 // TString processName[nProcess] = { "FCNCZut", "FCNCZct","fake","ttZ", "WZ", "tZq", "other", "data", "mc"};
  TString processLabel[14] = { "\\textbf{FCNC tZu}", "\\textbf{FCNC tZc}","\\textbf{DY amc}" ,"\\textbf{ttbar}","\\textbf{DD fake}", "\\textbf{WZ+jets}", "\\textbf{tZq}","\\textbf{ttZ}", "\\textbf{ZZ}", "\\textbf{other}","\\textbf{data}",  "\\textbf{MC (no fakes)}", "\\textbf{MC}", "\\textbf{MC +DDf}"};
  
   //TString processLabel[12] = { "\\textbf{FCNC tZu tt}", "\\textbf{FCNC tZu}", "\\textbf{FCNC tZc}","\\textbf{DY}" ,"\\textbf{ttbar}","\\textbf{non prompt lepton}", "\\textbf{WZ+jets}", "\\textbf{tZq}","\\textbf{ttZ}", "\\textbf{ZZ}", "\\textbf{other}","\\textbf{data}",  "\\textbf{MC (no fakes)}"};
  const int nCuts = 11;
  
  TString cutLabel[11] =  {"Z mass","$>2l$","STSR","TTSR","WZCR","TTCR","STCR","$Tr_{WZCR \\rightarrow STSR}$","$Tr_{WZCR \\rightarrow TTSR}$" ,"$Tr_{TTCR \\rightarrow TTSR}$", "$Tr_{STCR \\rightarrow STSR}$"};
  
  TH1F*  h[14];
  
  
  for(int i=0; i< 23; i++){
    cout << "adding " <<i << " " <<  processName[i] << endl;
    
    if(i == 0){ h[0] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); cout << "new " << endl;} // Zut
    else if(i  == 1 || i == 2){ h[0]->Add(((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone(""))); cout << "adding" << endl;}
    else if(i == 3){ h[1] = ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); cout << "new " << endl;} // Zct
    else if(i  == 4 || i == 5) h[1]->Add(((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")));
    else if(i == 6) { h[2] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); cout << "new " << endl;} // DY
    else if(i == 7){  h[3] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); }//TT
    else if(i == 8)  h[4] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); // fake
    else if(i == 9)  h[5] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); //WZ
    else if(i == 10)  h[6] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); // tZq
    else if(i == 11)  h[7] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); // ttZ
    else if(i == 12)  h[8] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); //ZZ
    else if(i == 13) h[9] =   ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); //ST,
    else if(i < 22 && i > 13) h[9]->Add( ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone(""))); // ST, twll; ZZZ; WZZ; ttW, tHq , ttH
    else if( i == 22){ h[10] =  ((TH1F*) _file0->Get("CutflowTableHisto_"+channelname+processName[i])->Clone("")); cout << "data" << endl; } // data
    
    
    
    // Lepton ID and HLT SF
    // if (mode == 0 && i != 8) h[i]->Scale(0.97713);
    // if (mode == 1 && i != 8) h[i]->Scale(0.910067);
    // if (mode == 2 && i != 8) h[i]->Scale(0.945736);
  }

  h[11] = (TH1F*) h[0]->Clone("");
  
   h[11]->Add(h[5]);
   h[11]->Add(h[6]);
   h[11]->Add(h[7]);
   h[11]->Add(h[8]);
   h[11]->Add(h[9]);
  
  h[13] = (TH1F*) h[11]->Clone("");
  h[12] = (TH1F*) h[11]->Clone("");
  h[12]->Add(h[2]);
  h[12]->Add(h[3]);
  
  h[13]->Add(h[4]);

  
  cout << "filling vector " << endl;
  double vectorValue[14][nCuts][3];
  
  
  for (int i = 0; i < nCuts-4; i++){
    //cout << "looking at cut " << cutLabel[i] << endl;
    for (int j = 0; j < 14; j++){
     // cout << "j " << j << endl;
     // cout << h[j]->GetBinContent(i+1) << endl;
      cout << "looking at cut " << cutLabel[i] << " for process " << processLabel[j] << endl;
      vectorValue[j][i][0] = h[j]->GetBinContent(i+1);
      vectorValue[j][i][1] = 1; //precision(h[j]->GetBinError(i));
      vectorValue[j][i][2] = h[j]->GetBinError(i+1);
      cout << vectorValue[j][i][0] << " pm " << vectorValue[j][i][2] << endl;
    }
  }
  
  for (int j = 0; j < 14; j++){
    // cout << "j " << j << endl;
    // cout << h[j]->GetBinContent(i+1) << endl;
    cout << "looking at cut " << cutLabel[7] << " for process " << processLabel[j] << endl;
    
    cout << h[j]->GetBinContent(7+1) <<" over " <<h[j]->GetBinContent(8+1)<<endl;
    vectorValue[j][7][0] = h[j]->GetBinContent(7+1)/h[j]->GetBinContent(8+1)*h[j]->GetBinContent(9+1);
    vectorValue[j][7][1] = 2; //precision(h[j]->GetBinError(i));
    double a = h[j]->GetBinContent(7+1);
    double sa = h[j]->GetBinError(7+1);
    double b= h[j]->GetBinContent(8+1);
    double sb = h[j]->GetBinError(8+1);
    double c = h[j]->GetBinContent(9+1);
    double sc = h[j]->GetBinError(9+1);
    vectorValue[j][7][2] = TMath::Sqrt(pow(sa,2.)*pow(c/b,2)+pow(a*c/pow(b,2.),2.)*pow(sb,2.)+pow(a/b,2.)*pow(sc,2.));
    cout << vectorValue[j][7][0] << " pm " << vectorValue[j][7][2] << endl;
    
    vectorValue[j][8][0] = h[j]->GetBinContent(7+1)/h[j]->GetBinContent(8+1) ;
    vectorValue[j][8][1] = 2; //precision(h[j]->GetBinError(i));
     a = h[j]->GetBinContent(7+1);
     sa = h[j]->GetBinError(7+1);
     b= h[j]->GetBinContent(8+1);
     sb = h[j]->GetBinError(8+1);
    vectorValue[j][8][2] = TMath::Sqrt(pow(sa,2.)*pow(1./b,2)+pow(a/pow(b,2.),2.)*pow(sb,2.));
    cout << vectorValue[j][8][0] << " pm " << vectorValue[j][8][2] << endl;
   
    if(j == 3){
    vectorValue[j][9][0] = 15./(60.-15.);
    vectorValue[j][9][1] = 2; //precision(h[j]->GetBinError(i));
    vectorValue[j][9][2] = 0.;
    //cout << vectorValue[j][8][0] << " pm " << vectorValue[j][8][2] << endl;

    vectorValue[j][10][0] = 15./(60.-15.);
    vectorValue[j][10][1] = 2; //precision(h[j]->GetBinError(i));
    vectorValue[j][10][2] = 0;
   // cout << vectorValue[j][9][0] << " pm " << vectorValue[j][9][2] << endl;
   // cout << h[j]->GetBinContent(6+1) << endl;
    }
    else{
      vectorValue[j][9][0] = 0.;
      vectorValue[j][9][1] = 1.;
      vectorValue[j][9][2] = 0.;

      vectorValue[j][10][0] = 0.;
      vectorValue[j][10][1] = 1.;
      vectorValue[j][10][2] = 0.;
    }
      
  }
  
  cout << "write table" << endl;
  
  salida << "\\documentclass[a4paper,12pt]{article}" << endl;
  salida << "\\usepackage[a4paper,bindingoffset=0.2in,left=1in,right=1in,top=1in,bottom=0.5in,[footskip=.25in]{geometry}" << endl;
  salida << "\\usepackage{lscape}" <<endl;
  salida << "\\begin{document}" << endl;
  salida << "\\thispagestyle{empty}" << endl;
  salida << endl;
  salida << endl;
  /*
  if(dolandscape)  salida << "  \\begin{landscape}" << endl;
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 0; i < 12; i++){
    salida << " & " << processLabel[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int i=2; i  < nCuts; i++){
    salida <<  cutLabel[i] ;
    for (int j = 0; j < 12; j++){
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
  if(dolandscape)  salida << "  \\end{landscape}" << endl;
  salida << endl;
  salida << endl;
  */
  if(dolandscape)  salida << "  \\begin{landscape}" << endl;
   salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l| ";
  for(int i = 2; i < nCuts; i++){
  salida << "c|";
  }
    salida << "}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 2; i < nCuts; i++){
    salida << " & " << cutLabel[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int j=0; j  < 14; j++){
    salida <<  processLabel[j] ;
    for (int i = 2; i < nCuts; i++){
      /*if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
       salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
       } else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){*/
      if ((i != 0 && vectorValue[j][i][0] == 0 )|| vectorValue[j][i][0]!= vectorValue[j][i][0] ){
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
  salida << " \\caption{Event yields: " ;
    
  if(mode == 9) salida << " All channel.";
  else if(mode == 0)  salida << " 3mu channel.";
  else if(mode == 1) salida << " 1e2mu channel.";
  else if(mode == 2) salida << " 2e1mu channel.";
  else if(mode == 3) salida << " 3e channel.";
  salida << "The last three columns represent the transfer factors that have to be applied.}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;

   salida << endl;
  salida << endl;

   if(dolandscape)  salida << "  \\end{landscape}" << endl;
  
  
  
  
  ///////////////////
  cutLabel[7] =  "STSR (WZCR)";
  cutLabel[8] =  "TTSR (WZCR)";
  cutLabel[9] ="TTSR (TTCR)";
  cutLabel[10] = "STSR (STCR)";
 
  
  for (int j = 0; j < 14; j++){
    // cout << "j " << j << endl;
    // cout << h[j]->GetBinContent(i+1) << endl;
    cout << "looking at cut " << cutLabel[7] << " for process " << processLabel[j] << endl;
    
     vectorValue[j][7][0] = (h[j]->GetBinContent(7+1)/h[j]->GetBinContent(8+1))*h[j]->GetBinContent(9+1); // 1jet
    vectorValue[j][7][1] = 2; //precision(h[j]->GetBinError(i));
    double a = h[j]->GetBinContent(7+1);
    double sa = h[j]->GetBinError(7+1);
    double b= h[j]->GetBinContent(8+1);
    double sb = h[j]->GetBinError(8+1);
    double c = h[j]->GetBinContent(9+1);
    double sc = h[j]->GetBinError(9+1);
    vectorValue[j][7][2] = TMath::Sqrt(pow(sa,2.)*pow(c/b,2)+pow(a*c/pow(b,2.),2.)*pow(sb,2.) + pow(sc,2.)*pow(a/b,2.));
    cout << vectorValue[j][7][0] << " pm " << vectorValue[j][7][2] << endl;
    
    vectorValue[j][8][0] = (h[j]->GetBinContent(7+1)/h[j]->GetBinContent(8+1))*h[j]->GetBinContent(4+1); // >0jet
    vectorValue[j][8][1] = 2; //precision(h[j]->GetBinError(i));
     a = h[j]->GetBinContent(7+1);
     sa = h[j]->GetBinError(7+1);
     b= h[j]->GetBinContent(8+1);
     sb = h[j]->GetBinError(8+1);
     c = h[j]->GetBinContent(4+1);
     sc = h[j]->GetBinError(4+1);
    vectorValue[j][8][2] = TMath::Sqrt(pow(sa,2.)*pow(c/b,2)+pow(a*c/pow(b,2.),2.)*pow(sb,2.) + pow(sc,2.)*pow(a/b,2.));
    cout << vectorValue[j][8][0] << " pm " << vectorValue[j][8][2] << endl;
    
    if(j == 3){
      vectorValue[j][9][0] = h[j]->GetBinContent(5+1)*15./(60.-15.);
      vectorValue[j][9][1] = 2; //precision(h[j]->GetBinError(i));
      vectorValue[j][9][2] = h[j]->GetBinError(5+1)*15./(60.-15.);
      cout << vectorValue[j][8][0] << " pm " << vectorValue[j][8][2] << endl;
      
      vectorValue[j][10][0] = h[j]->GetBinContent(6+1)*15./(60.-15.);
      vectorValue[j][10][1] = 2; //precision(h[j]->GetBinError(i));
      vectorValue[j][10][2] = h[j]->GetBinError(6+1)*15./(60.-15.);
      cout << vectorValue[j][9][0] << " pm " << vectorValue[j][9][2] << endl;
      // cout << h[j]->GetBinContent(6+1) << endl;
    }
    else{
      vectorValue[j][9][0] = 0.;
      vectorValue[j][9][1] = 1.;
      vectorValue[j][9][2] = 0.;
      
      vectorValue[j][10][0] = 0.;
      vectorValue[j][10][1] = 1.;
      vectorValue[j][10][2] = 0.;
    }
    
  }
  
  cout << "write table" << endl;
  
   salida << "\\thispagestyle{empty}" << endl;
  salida << endl;
  salida << endl;

  if(dolandscape)  salida << "  \\begin{landscape}" << endl;
  salida << "  \\begin{table}" << endl;
  salida << "  \\begin{center}" << endl;
  salida << "  \\begin{tabular} {|l| ";
  for(int i = 2; i < nCuts; i++){
    salida << "c|";
  }
  salida << "}" << endl;
  salida << "  \\hline " << endl;
  for (int i = 2; i < nCuts; i++){
    salida << " & " << cutLabel[i] ;
  }
  salida << "  \\\\ " << endl;
  salida << "  \\hline " << endl;
  
  for (int j=0; j  < 14; j++){
    salida <<  processLabel[j] ;
    for (int i = 2; i < nCuts; i++){
      /*if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] != 0){
       salida << " & $\\leq$ " << setprecision(vectorValue[j][i-1][1]) << 2*vectorValue[j][i-1][2];
       } else if (i != 0 && vectorValue[j][i][0] == 0 && vectorValue[j][i-1][0] == 0){*/
      if ((i != 0 && vectorValue[j][i][0] == 0 )|| vectorValue[j][i][0]!= vectorValue[j][i][0] ){
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
  salida << " \\caption{Event yields for " ;
  
  if(mode == 9) salida << " All channel.";
  else if(mode == 0)  salida << " 3mu channel.";
  else if(mode == 1) salida << " 1e2mu channel.";
  else if(mode == 2) salida << " 2e1mu channel.";
  else if(mode == 3) salida << " 3e channel.";
  salida << "The last tree columns represent the predicted event yield after applying the transfer factors. The regions between brakets correspond to the region from which the prediction is made.}" << endl;
  salida << "  \\end{center}" << endl;
  salida << "  \\end{table}" << endl;
  
  salida << endl;
  salida << endl;
  
  if(dolandscape)  salida << "  \\end{landscape}" << endl;

  salida << "\\end{document}" << endl;
  
}


