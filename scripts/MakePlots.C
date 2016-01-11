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
#include <sys/stat.h>



void MakePlots(string channel = "Mu", string date = "11_1_2016", bool plotData = false) {

 string dirname_ = "../Merged/ControlPlots_"+channel+"/ControlPlots_"+date+"/";
 const char * dirname = dirname_.c_str();
 cout << "Directory name:  " << dirname << endl;
 
 string chan = ""; 
 chan += channel; 
 chan += " channel";  
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
  const char *ext=".root";
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
     while ((file=(TSystemFile*)next())) {
       fname = file->GetName();
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
  TObjArray *oProc = Proc.Tokenize("_");
  TString ProcSample =  ((TObjString *)(oProc->At(9)))->String();
  TObjArray *oProcSample = ProcSample.Tokenize("."); 
  Vmyprocess.push_back(((TObjString *)(oProcSample->At(0)))->String());
  cout << ((TObjString *)(oProcSample->At(0)))->String() << endl; 
 }
 // Set the colors 
 vector<Color_t> color; 
 color.clear(); 
 for(unsigned int i = 0; i<Vmyprocess.size(); i++) {
cout << "DEBUG: " << Vmyprocess[i] << endl;
    if(Vmyprocess[i].Contains("diboson")){ color.push_back(kMagenta);}
    if(Vmyprocess[i].Contains("tZq")){ color.push_back(kBlue); }
    if(Vmyprocess[i].Contains("data")){ color.push_back(kBlack);}
    if(Vmyprocess[i].Contains("ttV")){ color.push_back(kCyan);}
    if(Vmyprocess[i].Contains("TTJets")){ color.push_back(kGreen); cout << "color filled" << endl;}
    if(Vmyprocess[i].Contains("Zjets")){ color.push_back(kBlue-3);}
    if(Vmyprocess[i].Contains("ST")){ color.push_back(kCyan+5);}
/*    if(Vmyprocess[i].Contains("WZ")){ color.push_back(kMagenta);} 
    if(Vmyprocess[i].Contains("WW")){ color.push_back(kMagenta+1);}
    if(Vmyprocess[i].Contains("ZZ")){ color.push_back(kMagenta+2);}
    if(Vmyprocess[i].Contains("tZq")){ color.push_back(kBlue); }
    if(Vmyprocess[i].Contains("data")){ color.push_back(kBlack);}
    if(Vmyprocess[i].Contains("ZZ")){ color.push_back(kGreen-2);}    
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
	 && (!obj->InheritsFrom("TH1")) 
       ) {
      printf("<W> Object %s is not 1D or 2D histogram : "
             "will not be converted\n",obj->GetName()) ;
    }
    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());
   listHisto.push_back(obj->GetName());
   listTitle.push_back(obj->GetTitle()); 
  }
  cout << " *** CHECKS *** " << endl;  
  for(unsigned int i = 0 ; i < Vmyprocess.size(); i++){
    cout << listrootfiles[i] << " " << Vmyprocess[i] << " " << color[i] << endl; 
  }

  cout << " **** DONE GETTING FILES, FILLING THStack *** " << endl; 
  const int sizeRF= listrootfiles.size(); 
  const int sizeH = listHisto.size(); 
  TH1F* h[sizeH][sizeRF]; 
  TH1F* histo[sizeH]; 
  THStack* hStack[sizeH]; 
  TH1F* addH;
  TH1F* hdata;
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
     if(myprocess.CompareTo("data")){
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

   string PlotsDir = "../Plots";
   string PlotsDir_ = PlotsDir+"/Plots_"+channel;
   string PlotsDir__ = PlotsDir_+"/Plots_"+date+"/";
   int mkdirstatus = mkdir(PlotsDir.c_str(),0777);
   mkdirstatus = mkdir(PlotsDir_.c_str(),0777);
   mkdirstatus = mkdir(PlotsDir__.c_str(),0777);


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
 
    string pngname = PlotsDir__; 
    pngname += listHisto[iVar];  
    pngname += "Ratio.png"; 
    string pdfname = PlotsDir__;
    pdfname += listHisto[iVar];
    pdfname += "Ratio.pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str());

    pad1->SetLogy(1); 
    string pngnamelogy = PlotsDir__;
    pngnamelogy += listHisto[iVar];
    pngnamelogy += "Ratio_logy.png";
    string pdfnamelogy = PlotsDir__;
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

    string pngname = PlotsDir__;
    pngname += listHisto[iVar];
    pngname += ".png";
    string pdfname = PlotsDir__;
    pdfname += listHisto[iVar];
    pdfname += ".pdf";
    c1->SaveAs(pngname.c_str());
    c1->SaveAs(pdfname.c_str()); 

    c1->SetLogy(); 
    string pngnamelogy = PlotsDir__;
    pngnamelogy += listHisto[iVar];
    pngnamelogy += "logy.png";
    string pdfnamelogy = PlotsDir__;
    pdfnamelogy += listHisto[iVar];
    pdfnamelogy += "logy.pdf";
    c1->SaveAs(pngnamelogy.c_str());
    c1->SaveAs(pdfnamelogy.c_str()); 



   }
  }
  cout << " **** DONE with THStack *** " << endl;
  

}
