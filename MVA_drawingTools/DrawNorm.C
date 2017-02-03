void DrawNorm()
{


//  string VariableName = "MVA_BDT";
  vector<std::string> MVAvars;


    MVAvars.push_back("MVA_TOPTOPLEPHAD");
    MVAvars.push_back("MVA_TOPTOPLEPHBB");
    MVAvars.push_back("MVA_TOPHLEPBB_hut");
    MVAvars.push_back("MVA_TOPHLEPBB_hct");
        //Variables for signal/background training
	  MVAvars.push_back("HiggsMass_TOPHLEPBB");
	  MVAvars.push_back("HiggsMass_TOPHLEPBB");
	  MVAvars.push_back("HiggsEta_TOPHLEPBB");
	  MVAvars.push_back("HiggsEta_TOPHLEPBB");
	  MVAvars.push_back("TopLepMass_TOPHLEPBB");
	  MVAvars.push_back("TopLepMass_TOPHLEPBB");
    MVAvars.push_back("TopLepPt_TOPHLEPBB");
    MVAvars.push_back("TopLepPt_TOPHLEPBB");
    MVAvars.push_back("TopLepEta_TOPHLEPBB");
    MVAvars.push_back("TopLepEta_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB");
    MVAvars.push_back("TopLepHiggsDr_TOPHLEPBB");
    MVAvars.push_back("TopLepHiggsDr_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPHLEPBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPHLEPBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPHLEPBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPHLEPBB");
    MVAvars.push_back("TopHadMass_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepMass_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepTopHadDr_TOPTOPLEPHAD");
    MVAvars.push_back("TopLepBJetCSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadBJetCSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD");
    MVAvars.push_back("HiggsMass_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepMass_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepHiggsDr_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet1CSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("HiggsBJet2CSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopLepBJetCSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetCSVv2_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB");
    MVAvars.push_back("TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD");
    MVAvars.push_back("TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD");

//  gSystem->Load("$ROOTSYS/test/libEvent");

  TFile *infile= new TFile("../weights/Training_TThut_All_b2j4.root","READ");
//  TFile *infile= new TFile("../weights/Training_TThut_All_b2j3.root","READ");
//  TFile *infile= new TFile("../weights/Training_SThut_All_b2j3.root","READ");
//  TFile *infile= new TFile("../weights/Training_SThct_All_b2j4.root","READ");
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
      else
      {
          cout << "Input histo doesn't exist" << endl;
          continue;
      }
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

      float maxhist= HistoNorm_S->GetMaximum();
      if (HistoNorm_B->GetMaximum() > maxhist) maxhist = HistoNorm_B->GetMaximum();
      
      maxhist = maxhist*1.5;
     
      HistoNorm_S->SetMaximum(maxhist);
      HistoNorm_B->SetMaximum(maxhist);

      TCanvas *c1 = new TCanvas();
      c1->cd();
      HistoNorm_S->Draw("hist");
      HistoNorm_B->Draw("hist same");

      TLegend *leg = new TLegend(0.7,0.75,0.98,0.95);
       leg->AddEntry(HistoNorm_S,"Signal","f");
       leg->AddEntry(HistoNorm_B,"Background","f");
       leg->Draw();



      c1->SaveAs(("pics/"+MVAvars[i_vars]+"_Norm.png").c_str());

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

      float maxhist= HistoNorm_S->GetMaximum();
      if (HistoNorm_B->GetMaximum() > maxhist) maxhist = HistoNorm_B->GetMaximum();
      
      maxhist = maxhist*1.5;
     
      HistoNorm_S->SetMaximum(maxhist);
      HistoNorm_B->SetMaximum(maxhist);


      TCanvas *c1 = new TCanvas();
      c1->cd();
      HistoNorm_B->Draw("hist");
      HistoNorm_S->Draw("hist same");

      TLegend *leg = new TLegend(0.7,0.75,0.98,0.95);
       leg->AddEntry(HistoNorm_S,"Signal","f");
       leg->AddEntry(HistoNorm_B,"Background","f");
       leg->Draw();



      c1->SaveAs("pics/BDT_Norm.png");

    delete c1;

  fout->Write();
  

}
