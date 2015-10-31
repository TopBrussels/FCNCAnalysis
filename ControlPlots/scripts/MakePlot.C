#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
//#include "TStyle.h"
#include "TH2.h"
#include "TKey.h"
#include "setTDRStyle.C"
#include <cmath>
#include "THStack.h"

void MakePlot(string channel = "ee", bool plotData = false) {


 
 string chan = ""; 
 chan += channel; 
 chan += " channel";  
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
 gStyle->SetErrorX(0);
 setTDRStyle();
 gROOT->SetBatch(1);

 labelcms2  = new TPaveText(0.12,0.85,0.5,0.88,"NDCBR");
 labelcms2->SetTextAlign(12);
 labelcms2->SetTextSize(0.04);
 labelcms2->SetFillColor(kWhite);
 labelcms2->AddText(chan.c_str());
 labelcms2->SetBorderSize(0);

 labelcms  = new TPaveText(0.12,0.88,0.5,0.94,"NDCBR");
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
  if(chan.find("ee") == 0) dirname = "../_ElEl_allSamples/";
  else if(chan.find("mumu") == 0) dirname = "../_MuMu_allSamples/";
  else if(chan.find("emu") == 0) dirname = "../_ElMu_allSamples/";
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
 for(int k =0; k<listrootfiles.size(); k++)
 { 
  TString Proc = listrootfiles[k]; 
  TObjArray *oProc = Proc.Tokenize("_"); 
  TString ProcSample =  ((TObjString *)(oProc->At(4)))->String(); 
  TObjArray *oProcSample = ProcSample.Tokenize("."); 
  Vmyprocess.push_back(((TObjString *)(oProcSample->At(0)))->String());  
  cout << ((TObjString *)(oProcSample->At(0)))->String() << endl; 
 }
 // Set the colors 
 vector<Color_t> color; 
 color.clear(); 
 for(int i = 0; i<Vmyprocess.size(); i++) {
    if(Vmyprocess[i].CompareTo("WZ")) color.push_back(kMagenta); 
    if(Vmyprocess[i].CompareTo("tZq")) color.push_back(kBlue); 
    if(Vmyprocess[i].CompareTo("data")) color.push_back(kBlack);
 } 
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


  cout << " **** DONE GETTING FILES, FILLING THStack *** " << endl; 
  const int sizeRF= listrootfiles.size(); 
  const int sizeH = listHisto.size(); 
//  Color_t color[3] = {kMagenta, kBlue, kBlack};
  TH1F* h[sizeH][sizeRF]; 
  TH1F* histo[sizeH]; 
  THStack* hStack[sizeH]; 
  TH1F* addH;
  TH1F* hdata;
//  vector<string> Vmyprocess = {"WZ","tZq","data"};
  
  for(int iVar = 0; iVar < listHisto.size(); iVar++)
  {
    leg = new TLegend(0.7,0.7,0.96,0.96);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);

    hStack[iVar] = new THStack(listHisto[iVar].c_str(), listTitle[iVar].c_str());
    int datanb = -1;
    int counter = 1; 
    for(int iProcess = 0; iProcess < listrootfiles.size(); iProcess++)
    {
         
     TString myprocess = Vmyprocess[iProcess];
     TString myrootfile = listrootfiles[iProcess];
     TFile *_file0 = TFile::Open(myrootfile);
     h[iVar][iProcess] = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
     if(myprocess.CompareTo("data")!=0){
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
         h[iVar][iProcess]->SetMarkerStyle(20);
         h[iVar][iProcess]->SetMarkerSize(0.5);
         h[iVar][iProcess]->SetLineWidth(1);
         h[iVar][iProcess]->SetMarkerColor(kBlack);
         h[iVar][iProcess]->SetLineColor(kBlack);
         hdata = (TH1F*) _file0->Get((listHisto[iVar]).c_str());
     }
    }
    
   
   double max = 0; 
   TCanvas *c1 = new TCanvas();

   if(plotData)
   {
    max = TMath::Max(hStack[iVar]->GetMaximum(), h[iVar][datanb]->GetMaximum());
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    //pad1->SetLogy(1);
    pad1->Draw();
    pad1->cd();     
    hStack[iVar]->Draw("histo");
    hStack[iVar]->SetMaximum(max * 1.2);
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
    histo[iVar]->GetYaxis()->SetTitle("MC/Data");
    histo[iVar]->GetYaxis()->CenterTitle(); 
    histo[iVar]->SetMarkerStyle(2);
    histo[iVar]->Draw("ep");
    histo[iVar]->Draw(); 
    c1->cd(); 
 
    string pngname = "../1DPlot/"; 
    pngname += listHisto[iVar];  
    pngname += "Ratio.png"; 
    string pdfname = "../1DPlot/";
    pdfname += listHisto[iVar];
    pdfname += "Ratio.pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str());
  }
  else{
    max  = hStack[iVar]->GetMaximum();
    hStack[iVar]->Draw("histo");
    hStack[iVar]->SetMaximum(max * 1.2);
    //hStack[iVar]->SetMinimum(1);
    hStack[iVar]->GetXaxis()->SetTitle(listTitle[iVar].c_str());
    hStack[iVar]->GetYaxis()->SetTitle("events / lumi fb^{-1}");
    hStack[iVar]->GetYaxis()->CenterTitle();

    labelcms->Draw();
    labelcms2->Draw();
    leg->Draw();

    string pngname = "../1DPlot/";
    pngname += listHisto[iVar];
    pngname += ".png";
    string pdfname = "../1DPlot/";
    pdfname += listHisto[iVar];
    pdfname += ".pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str()); 

    c1->SetLogy(); 
    string pngnamelogy = "../1DPlot/";
    pngnamelogy += listHisto[iVar];
    pngnamelogy += "logy.png";
    string pdfnamelogy = "../1DPlot/";
    pdfnamelogy += listHisto[iVar];
    pdfnamelogy += "logy.pdf";
    c1->SaveAs(pngnamelogy.c_str());
    c1->SaveAs(pdfnamelogy.c_str()); 



   }
  }
  cout << " **** DONE with THStack *** " << endl;
  

}
