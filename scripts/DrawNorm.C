void DrawNorm()
{


//  string VariableName = "MVA_BDT";
  vector<std::string> MVAvars;


      MVAvars.push_back("SumCharge_Hjets");
      MVAvars.push_back("SumCharge_TopJets");
      MVAvars.push_back("SumCharge_FCNHJetLep");
      MVAvars.push_back("CvsL_Hjet1");
      MVAvars.push_back("CvsL_Hjet2");
      MVAvars.push_back("CvsL_SMb");
      MVAvars.push_back("CvsL_FCNHjet");
      MVAvars.push_back("CvsB_Hjet1");
      MVAvars.push_back("CvsB_Hjet2");
      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("CvsB_FCNHjet");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("TransvLepTopmass");
      MVAvars.push_back("HadTopmass");
      MVAvars.push_back("DR_H_HadTop");
      MVAvars.push_back("DPhi_H_LepTop");
      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HadTopPt");
      MVAvars.push_back("JetCombBDT");

//  gSystem->Load("$ROOTSYS/test/libEvent");

  TFile *infile= new TFile("../MVA/TrainFiles/EventMVA_SThypo_4Vars_3B1J_El_NP_overlay_ST_tHToBB_1L_Kappa_hct_19vars.root","READ");
  TFile *fout = new TFile("NormHistos.root","RECREATE");

  for(int i_vars = 0; i_vars < MVAvars.size(); i_vars++)
  {   
      TH1F *HistoNorm_S = 0;
      TH1F *HistoNorm_B = 0;
	    TH1F *histo_S( (TH1F*) infile->Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Signal").c_str()));
	    TH1F *histo_B( (TH1F*) infile->Get(("Method_BDT/BDT/"+MVAvars[i_vars]+"__Background").c_str()));

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
      
      HistoNorm_S->GetXaxis()->SetTitle(MVAvars[i_vars].c_str());
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

      TLegend *leg = new TLegend(0.7,0.75,0.98,0.95);
       leg->AddEntry(HistoNorm_S,"Signal","f");
       leg->AddEntry(HistoNorm_B,"Background","f");
       leg->Draw();



      c1->SaveAs((MVAvars[i_vars]+"_Norm.png").c_str());

    delete c1;
  }

      TH1F *HistoNorm_S = 0;
      TH1F *HistoNorm_B = 0;
	    TH1F *histo_S( (TH1F*) infile->Get("Method_BDT/BDT/MVA_BDT_S"));
	    TH1F *histo_B( (TH1F*) infile->Get("Method_BDT/BDT/MVA_BDT_B"));

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
      
      HistoNorm_S->GetXaxis()->SetTitle("BDT");
      HistoNorm_S->SetLineColor(4);
      HistoNorm_B->SetLineColor(2);
      HistoNorm_S->SetFillColor(4);
      HistoNorm_B->SetFillColor(2);
      HistoNorm_S->SetFillStyle(3004);
      HistoNorm_B->SetFillStyle(3005);

      TCanvas *c1 = new TCanvas();
      c1->cd();
      HistoNorm_B->Draw("hist");
      HistoNorm_S->Draw("hist same");

      TLegend *leg = new TLegend(0.7,0.75,0.98,0.95);
       leg->AddEntry(HistoNorm_S,"Signal","f");
       leg->AddEntry(HistoNorm_B,"Background","f");
       leg->Draw();



      c1->SaveAs("BDT_Norm.png");

    delete c1;

  fout->Write();
  

}
