void DrawNorm()
{


//  string VariableName = "MVA_BDT";
//  string VariableName = "InclCharges_Hjets";
//  string VariableName = "InclCharges_SMb_Lep";
//  string VariableName = "Charges_Hjets";
//  string VariableName = "Charges_SMb_Lep";
//  string VariableName = "CvsL_Hjet1";
//  string VariableName = "CvsL_Hjet2";
//  string VariableName = "CvsL_SMTopJet";
//  string VariableName = "CvsB_Hjet1";
//  string VariableName = "CvsB_Hjet2";
//  string VariableName = "CvsB_SMTopJet";
//  string VariableName = "bDisc_Hjet1";
//  string VariableName = "bDisc_Hjet2";
//  string VariableName = "bDisc_SMTopJet";
//  string VariableName = "InclCharges_jetTop_jetAntiTop";
//  string VariableName = "inclCharges_JetFCNH_Lep";
  string VariableName = "Charges_jetTop_jetAntiTop";
//  string VariableName = "Charges_JetFCNH_Lep";
//  string VariableName = "CvsL_FCNHJet";
//  string VariableName = "CvsB_FCNHJet";
//  string VariableName = "bDisc_FCNHJet";
  
//  gSystem->Load("$ROOTSYS/test/libEvent");

  TFile *infile= new TFile("../MVA/TrainFilesTrainedJetCombMVA_El_NP_overlay_TTtoTHToBB-1L-Kappa-hct.root","READ");
  TFile *fout = new TFile("NormHistos.root","RECREATE");

 
  TH1F *HistoNorm_S = 0;
  TH1F *HistoNorm_B = 0;
	TH1F *histo_S( (TH1F*) infile->Get(("Method_BDT/BDT/"+VariableName+"__Signal").c_str()));
	TH1F *histo_B( (TH1F*) infile->Get(("Method_BDT/BDT/"+VariableName+"__Background").c_str()));

  if(histo_S)
  {
    HistoNorm_S = new TH1F(*histo_S);
  }
  else cout << "Input histo doesn't exist" << endl;
  if(histo_B)
  {
    HistoNorm_B = new TH1F(*histo_B);
  }
  else cout << "Input histo doesn't exist" << endl;

  Double_t norm_S = 1;
  Double_t scale_S = norm_S/(HistoNorm_S->Integral());
  HistoNorm_S->Scale(scale_S);
  Double_t norm_B = 1;
  Double_t scale_B = norm_B/(HistoNorm_B->Integral());
  HistoNorm_B->Scale(scale_B);
  
  HistoNorm_S->GetXaxis()->SetTitle(VariableName.c_str());
  HistoNorm_S->SetLineColor(4);
  HistoNorm_B->SetLineColor(2);
  HistoNorm_S->SetFillColor(4);
  HistoNorm_B->SetFillColor(2);
  HistoNorm_S->SetFillStyle(3004);
  HistoNorm_B->SetFillStyle(3005);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  HistoNorm_S->Draw("hist");
  HistoNorm_B->Draw("hist same");

   leg = new TLegend(0.7,0.75,0.98,0.95);
   leg->AddEntry(HistoNorm_S,"Signal","f");
   leg->AddEntry(HistoNorm_B,"Background","f");
   leg->Draw();



  c1->SaveAs((VariableName+"_Norm.png").c_str());
  
  fout->Write();
  

}
