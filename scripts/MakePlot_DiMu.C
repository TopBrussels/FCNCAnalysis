#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include "TStyle.h"
#include "TH2.h"
#include "TKey.h"
#include "setTDRStyle.C"
#include <cmath>
#include "THStack.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TRint.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"



void MakePlot_DiMu(string channel = "Mu_Mu", bool plotData = true) {


 
 string chan = ""; 
 chan += channel; 
 chan += " channel";
 chan += ", 2.16 fb^{-1}, 2015D";  
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
 gStyle->SetErrorX(0);
 setTDRStyle();
 gROOT->SetBatch(1);

 TPaveText* labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
 labelcms2->SetTextAlign(12);
 labelcms2->SetTextSize(0.04);
 labelcms2->SetFillColor(kWhite);
 labelcms2->AddText(chan.c_str());
 labelcms2->SetBorderSize(0);

 TPaveText* labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
 labelcms->SetTextAlign(12);
 labelcms->SetTextSize(0.04);
 labelcms->SetFillColor(kWhite);
 labelcms->AddText("CMS Preliminary, #sqrt{s} = 13 TeV");
 labelcms->SetBorderSize(0);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(600);
  gStyle->SetLabelFont(18,"");
  
  gStyle->SetTitleXOffset(1.2);//1.5
  gStyle->SetTitleYOffset(1.2);//1.7


  

 cout << " **** GETTING FILES *** " << endl;
 cout << " - Rootfiles: " << endl; 
  // Get all rootfiles
  vector<TString> listrootfiles; 
  listrootfiles.clear(); 
  const char *dirname; 
  if(chan.find("Elec_Elec") == 0) dirname = "../OutPutHistos/merged_ElEl_AllSamples/";
  else if(chan.find("Mu_Mu") == 0) dirname = "../OutPutHistos/merged_MuMu_AllSamples/";
  else if(chan.find("emu") == 0) dirname = "../_ElMu/";
  const char *ext=".root";
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
     while ((file=(TSystemFile*)next())) {
       //cout << "in while " << endl; 
       fname = file->GetName();
       // cout << fname << endl; 
       //if(fname.EndsWith(ext)) cout << "is a root file" << endl; 
       //if (!file->IsDirectory()) cout << "not a dir" << endl; 
       if (!file->IsDirectory() && fname.EndsWith(ext)) {
           cout << fname << endl;
           listrootfiles.push_back(dirname+fname);
        }
     }
  }
  else cout << " no files found " << endl; 


 /// Get ProcessNames
 vector<TString> Vmyprocess; 
 Vmyprocess.clear(); 
 cout << " - samples " << endl; 
 for(unsigned int k =0; k<listrootfiles.size(); k++)
 { 
  TString Proc = listrootfiles[k]; 
 // cout << "1" << endl; 
  TObjArray *oProc = Proc.Tokenize("_");
 // cout << "2" << endl;
 // for(unsigned int o = 0; o < oProc->GetEntries(); o++)
 // {
 //    cout << ((TObjString *)(oProc->At(3)))->String() << endl; 
 // } 
//  TString ProcSample =  ((TObjString *)(oProc->At(4)))->String(); 
  TString ProcSample =  ((TObjString *)(oProc->At(3)))->String();
//  cout << "3" << endl;
  TObjArray *oProcSample = ProcSample.Tokenize("."); 
//  cout << "4" << endl;
  Vmyprocess.push_back(((TObjString *)(oProcSample->At(0)))->String());  
  cout << ((TObjString *)(oProcSample->At(0)))->String() << endl; 
 }
 // Set the colors 
 vector<Color_t> color; 
 color.clear(); 
 for(unsigned int i = 0; i<Vmyprocess.size(); i++) {
    if(Vmyprocess[i].Contains("DYJets50")){ color.push_back(kMagenta);}
    //if(Vmyprocess[i].Contains("DYJets10-50")){ color.push_back(kBlue);}
    if(Vmyprocess[i].Contains("Data")){ color.push_back(kBlack);}
    if(Vmyprocess[i].Contains("TTJets")){ color.push_back(kCyan);}
    //if(Vmyprocess[i].Contains("TTJets")){ color.push_back(kGreen);}
   // if(Vmyprocess[i].Contains("DYJets50")){ color.push_back(kBlue-3);}
   // if(Vmyprocess[i].Contains("ST")){ color.push_back(kViolet);}
   // if(Vmyprocess[i].Contains("WJets")){ color.push_back(kCyan+3);} 
    if(Vmyprocess[i].Contains("WZ")){ color.push_back(kGreen);}
    //if(Vmyprocess[i].Contains("ZZ")){ color.push_back(kMagenta+2);}
    //if(Vmyprocess[i].Contains("tZq")){ color.push_back(kBlue); }
   // if(Vmyprocess[i].Contains("WW")){ color.push_back(kCyan+2);}
   /* if(Vmyprocess[i].Contains("ZZ")){ color.push_back(kGreen-2);}    
    if(Vmyprocess[i].Contains("ttbar")){ color.push_back(kGreen);}
    if(Vmyprocess[i].Contains("Zjets1050")){ color.push_back(kBlue-2);}
    if(Vmyprocess[i].Contains("Zjets50")){ color.push_back(kBlue-3);} 
    if(Vmyprocess[i].Contains("ttZ")){ color.push_back(kMagenta-3);} 
    if(Vmyprocess[i].Contains("Synch")){ color.push_back(kMagenta-3);}
    if(Vmyprocess[i].Contains("ttW")){ color.push_back(kCyan);}
    if(Vmyprocess[i].Contains("STs")){ color.push_back(kCyan+3);} 
    if(Vmyprocess[i].Contains("STt") && !Vmyprocess[i].Contains("STtW")){ color.push_back(kCyan+4);}
    if(Vmyprocess[i].Contains("STtW")){ color.push_back(kCyan+5);}
*/ } 
  // list all histograms in the rootfiles  
  TFile *f1 = new TFile(listrootfiles[0]);

  vector<string> listHisto;
  listHisto.clear(); 
  vector<string> listTitle;
  listTitle.clear();

  if (!f1->IsOpen()) {
    printf("<E> Cannot open input file %s\n",".") ;
    exit(1) ;
  }
  cout << " - Histograms: " << endl; 
  
  TList* list = f1->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
      
  while (( key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
//         && (!obj->InheritsFrom("TH2"))
	 && (!obj->InheritsFrom("TH1")) 
       ) {
      printf("<W> Object %s is not 1D or 2D histogram : "
             "will not be converted\n",obj->GetName()) ;
    }
    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());
   listHisto.push_back(obj->GetName());
   listTitle.push_back(obj->GetTitle()); 
  }
/*  cout << " *** CHECKS *** " << endl;  
  for(unsigned int i = 0 ; i < Vmyprocess.size(); i++){
    cout << listrootfiles[i] << " " << Vmyprocess[i] << " " << color[i] << endl; 

  }
*/
  cout << " **** DONE GETTING FILES, FILLING THStack *** " << endl; 
  const int sizeRF= listrootfiles.size(); 
  const int sizeH = listHisto.size(); 
  TH1F* h[sizeH][sizeRF]; 
  TH1F* histo[sizeH]; 
  THStack* hStack[sizeH]; 
  TH1F* addH;
  TH1F* hdata;
//  vector<string> Vmyprocess = {"WZ","tZq","data"};
  
  for(unsigned int iVar = 0; iVar < listHisto.size(); iVar++)
  {
    TLegend* leg = new TLegend(0.7,0.7,0.96,0.96);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    hStack[iVar] = new THStack(listHisto[iVar].c_str(), listTitle[iVar].c_str());
    int datanb = -1;
    unsigned int counter = 1; 
    for(unsigned int iProcess = 0; iProcess < listrootfiles.size(); iProcess++)
    {
         
     TString myprocess = Vmyprocess[iProcess];
     TString myrootfile = listrootfiles[iProcess];
     TFile *_file0 = TFile::Open(myrootfile);
     h[iVar][iProcess] = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
   //  if(! h[iVar][iProcess]) continue; 
     if(myprocess.CompareTo("Data")){
       leg->AddEntry(h[iVar][iProcess], Vmyprocess[iProcess],"f");
       h[iVar][iProcess]->SetFillColor(color[iProcess]);
       h[iVar][iProcess]->SetLineColor(kBlack);
       h[iVar][iProcess]->SetLineWidth(1);
       hStack[iVar]->Add(h[iVar][iProcess]);
       if(iProcess < counter){
           histo[iVar] = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
           histo[iVar]->Sumw2();
       }
       else{
           addH = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
           addH->Sumw2(); 
           histo[iVar]->Add(addH,1);
      }
     }
     else{
         datanb = iProcess;
         h[iVar][datanb]->SetMarkerStyle(20);
         h[iVar][datanb]->SetMarkerSize(0.5);
         h[iVar][datanb]->SetLineWidth(1);
         h[iVar][datanb]->SetMarkerColor(kBlack);
         h[iVar][datanb]->SetLineColor(kBlack);
         hdata = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
         cout << " data? " << Vmyprocess[iProcess] << endl; 
     }
    }
    
   
   double max = 0; 
   TCanvas *c1 = new TCanvas();

   if(plotData)
   {
    cout << " PLOTTING DATA" <<endl; 
    max = TMath::Max(hStack[iVar]->GetMaximum(), h[iVar][datanb]->GetMaximum());
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    //pad1->SetLogy(1);
    pad1->Draw();
    pad1->cd();   
    hStack[iVar]->Draw("histo");
    hStack[iVar]->SetMaximum(max * 1.6);
    //hStack[iVar]->SetMinimum(1);
    hStack[iVar]->GetXaxis()->SetTitle(listTitle[iVar].c_str());
    hStack[iVar]->GetYaxis()->SetTitle("events / lumi fb^{-1}");    
    hStack[iVar]->GetYaxis()->CenterTitle(); 
    
    h[iVar][datanb]->Draw("e,same");
    leg->AddEntry(h[iVar][datanb], Vmyprocess[datanb],"p");   
    
    labelcms->Draw();
    labelcms2->Draw(); 
    leg->Draw(); 
    c1->cd(); 
 
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();
    hdata->Sumw2();
    histo[iVar]->Divide(hdata);
    histo[iVar]->SetMarkerStyle(21);
    histo[iVar]->SetTitle("");
    histo[iVar]->GetXaxis()->SetTitle(listTitle[iVar].c_str());
    histo[iVar]->GetXaxis()->SetTitleSize(0.08);
    histo[iVar]->GetYaxis()->SetTitleSize(0.04);
    histo[iVar]->GetYaxis()->SetTitle("MC/Data");
    histo[iVar]->GetYaxis()->CenterTitle(); 
    histo[iVar]->SetMarkerStyle(2);
    histo[iVar]->Draw("ep");
    histo[iVar]->Draw(); 
    c1->cd(); 
 
    string pngname = "../OutPutHistos/1DPlot/DiMu/"; 
    pngname += listHisto[iVar];  
    pngname += "Ratio.png"; 
    string pdfname = "../OutPutHistos/1DPlot/DiMu/";
    pdfname += listHisto[iVar];
    pdfname += "Ratio.pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str());

    pad1->SetLogy(1); 
    string pngnamelogy = "../OutPutHistos/1DPlot/DiMu/";
    pngnamelogy += listHisto[iVar];
    pngnamelogy += "Ratio_logy.png";
    string pdfnamelogy = "../OutPutHistos/1DPlot/DiMu/";
    pdfnamelogy += listHisto[iVar];
    pdfnamelogy += "Ratio_logy.pdf";
    c1->SaveAs(pngnamelogy.c_str());
    c1->SaveAs(pdfnamelogy.c_str());   
  }
  else{
   cout << " NOT PLOTTING DATA " << endl; 
    max  = hStack[iVar]->GetMaximum();
    hStack[iVar]->Draw("histo");
    hStack[iVar]->SetMaximum(max * 1.6);
    //hStack[iVar]->SetMinimum(1);
    hStack[iVar]->GetXaxis()->SetTitle(listTitle[iVar].c_str());
    hStack[iVar]->GetYaxis()->SetTitle("events / lumi fb^{-1}");
    hStack[iVar]->GetYaxis()->CenterTitle();

    labelcms->Draw();
    labelcms2->Draw();
    leg->Draw();

    string pngname = "../OutPutHistos/1DPlot/DiMu/";
    pngname += listHisto[iVar];
    pngname += ".png";
    string pdfname = "../OutPutHistos/1DPlot/DiMu/";
    pdfname += listHisto[iVar];
    pdfname += ".pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str()); 

    c1->SetLogy(); 
    string pngnamelogy = "../OutPutHistos/1DPlot/DiMu/";
    pngnamelogy += listHisto[iVar];
    pngnamelogy += "logy.png";
    string pdfnamelogy = "../OutPutHistos/1DPlot/DiMu/";
    pdfnamelogy += listHisto[iVar];
    pdfnamelogy += "logy.pdf";
    c1->SaveAs(pngnamelogy.c_str());
    c1->SaveAs(pdfnamelogy.c_str()); 



   }
  }
  cout << " **** DONE with THStack *** " << endl;
  

}
