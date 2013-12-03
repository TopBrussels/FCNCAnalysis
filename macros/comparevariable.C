#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "tdrstyle.C"
#include "TString.h"

void comparevariable(const TString & variable = "jet_Pt", const TString & title = "Leading Jet p_{T} (GeV)" ){

  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();

  TFile * f = new TFile("../data/FCNC_selection_3L.root");

  const int ndataset = 17;
  TH1F * h[ndataset];
  h[0] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToWW_WToLNuL_HctL",variable.Data()));
  h[1] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToWW_WToLNuL_HctR",variable.Data()));
  h[2] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL",variable.Data()));
  h[3] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR",variable.Data()));
  h[4] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL",variable.Data()));
  h[5] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR",variable.Data()));
  h[6] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToLL_HctL",variable.Data()));
  h[7] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToLL_HctR",variable.Data()));
  h[8] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL",variable.Data()));
  h[9] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR",variable.Data()));
  h[10] = (TH1F *) f->Get(Form("Histos1D/%s_TTJetsTocZbW",variable.Data()));
  h[11] = (TH1F *) f->Get(Form("Histos1D/%s_ttbar",variable.Data()));
  h[12] = (TH1F *) f->Get(Form("Histos1D/%s_ttz",variable.Data()));
  h[13] = (TH1F *) f->Get(Form("Histos1D/%s_ttw",variable.Data()));
  h[14] = (TH1F *) f->Get(Form("Histos1D/%s_Zjets",variable.Data()));
  h[15] = (TH1F *) f->Get(Form("Histos1D/%s_zz",variable.Data()));
  h[16] = (TH1F *) f->Get(Form("Histos1D/%s_wz",variable.Data()));
  
  TCanvas * c = new TCanvas("c","c",800,600);
  h[0]->Rebin(10);
  h[0]->Scale(1.0/h[0]->Integral());
  h[0]->Draw();
  h[0]->SetMaximum(0.3);
  h[0]->SetTitle("");
  h[0]->GetYaxis()->SetTitle("Normalized Entries");
  h[0]->GetXaxis()->SetTitle(Form("%s",title.Data()));
  h[0]->SetLineWidth(2);
  h[0]->SetStats(0);
  for(int i = 1 ; i < ndataset; i++){
    h[i]->Rebin(10);
    h[i]->Scale(1.0/h[i]->Integral());
    h[i]->SetLineWidth(2);
    h[i]->Draw("same");
    if( i < 9) {
      h[i]->SetLineColor(i+1);
    }else if( i < 11){
      h[i]->SetLineColor(i+2);
    }else{
      h[i]->SetLineStyle(2);
      h[i]->SetLineColor(i-10);
    }
  } 

  TLegend *l = new TLegend(0.42,0.50,0.90,0.90);
  for(int i = 0 ; i < ndataset; i++){
    l->AddEntry(h[i],h[i]->GetName(),"L");
  }
  l->SetTextSize(0.02);
  l->SetFillColor(0);
  l->SetLineColor(0);
  l->Draw();


  c->Print(Form("%s.png", variable.Data()));

}



