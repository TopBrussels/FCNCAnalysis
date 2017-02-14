map<string,TH1F*> histo1D;


void DrawNorm()
{
  TFile *sample1= new TFile("../weights/Training_SThut_All_b2j3.root","READ");
  std::string label_sample1 = "ST-hct FullSim";
  TFile *sample2= new TFile("/user/kskovpen/analysis/tHFCNC_OLD/CMSSW_8_0_12/src/tHFCNC/NtupleAnalyzer/test/MVA/TMVA_HutST_all_b2j3.root","READ");
  std::string label_sample2 = "ST-hct FastSim";

  string treename = "InputVariables_Id";
  TDirectory *tree_sample1 = (TDirectory*)sample1->Get(treename.c_str());
cout << tree_sample1->GetName() << endl;
  TDirectory *tree_sample2 = (TDirectory*)sample2->Get(treename.c_str());
  
  string condition = "";

  vector<std::string> vars;
  vector<int> nbins;
  vector<double> xmin;
  vector<double> xmax;



  //Defining the variables we want to plot with the nbins, xmin and xmax
  vars.push_back("HiggsBJet1CSVv2_TOPHLEPBB__Signal_Id");
  nbins.push_back(50);
  xmin.push_back(-1);
  xmax.push_back(1);

  TFile *fout = new TFile("NormHistos_Comp2Samples_FullSimFastSim.root","RECREATE");

  for(int i_vars = 0; i_vars < vars.size(); i_vars++)
  {   
cout << "HERE 1" << endl;
      tree_sample1->cd();
cout << "HERE 2: " << vars[i_vars] << endl;
      histo1D[(vars[i_vars]+"_sample1").c_str()] = (TH1F*) tree_sample1->Get(vars[i_vars].c_str())->Clone();
cout << "HERE 3" << endl;
      tree_sample2->cd();
      histo1D[(vars[i_vars]+"_sample2").c_str()] = (TH1F*) tree_sample2->Get(vars[i_vars].c_str())->Clone();
      gStyle->SetOptStat(kFALSE);
//      tree_sample1->Draw((vars[i_vars]+">>h_"+vars[i_vars]+"_sample1").c_str(),condition.c_str());
//      tree_sample2->Draw((vars[i_vars]+">>h_"+vars[i_vars]+"_sample2").c_str(),condition.c_str());


      if(!histo1D[(vars[i_vars]+"_sample1").c_str()] || !histo1D[(vars[i_vars]+"_sample2").c_str()])
      {
          cout << "Input histo doesn't exist" << endl;
          continue;
      }

//      histo1D[(vars[i_vars]+"_sample1").c_str()] = tmp_hist_sample1->Clone();
//      histo1D[(vars[i_vars]+"_sample2").c_str()] = tmp_hist_sample2->Clone();

      Double_t norm_sample1 = 1;
      Double_t scale_sample1 = norm_sample1/(histo1D[(vars[i_vars]+"_sample1").c_str()]->Integral());
      histo1D[(vars[i_vars]+"_sample1").c_str()]->Scale(scale_sample1);
      Double_t norm_sample2 = 1;
      Double_t scale_sample2 = norm_sample2/(histo1D[(vars[i_vars]+"_sample2")]->Integral());
      histo1D[(vars[i_vars]+"_sample2")]->Scale(scale_sample2);

      float maxhist= histo1D[(vars[i_vars]+"_sample1").c_str()]->GetMaximum();
      if (histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum() > maxhist) maxhist = histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum();
      
      maxhist = maxhist*1.3;
     
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetMaximum(maxhist);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetMaximum(maxhist);

      
      histo1D[(vars[i_vars]+"_sample1").c_str()]->GetXaxis()->SetTitle(vars[i_vars].c_str());
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetLineColor(4);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetLineColor(2);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetFillColor(4);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetFillColor(2);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetFillStyle(3004);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetFillStyle(3005);


      TH1F* histo_ratio;
      histo_ratio = (TH1F*) histo1D[(vars[i_vars]+"_sample1")]->Clone();
//      histo_ratio->Sumw2();
      histo_ratio->SetName("histo_ratio");
      histo_ratio->SetTitle("");
    
      histo_ratio->Divide(histo1D[(vars[i_vars]+"_sample2")]);


      TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,552);
      c1->SetFillColor(10);
      c1->  cd();   
    

      TPad* canvas_1 = new TPad("canvas_1", "canvas_1",0,0.25,1.0,0.98);
      canvas_1 ->Draw();
      canvas_1 ->cd();

      histo1D[(vars[i_vars]+"_sample1").c_str()]->Draw("hist");
      histo1D[(vars[i_vars]+"_sample2").c_str()]->Draw("hist same");

      TLegend* leg =  new TLegend(0.65,0.60,0.85,0.85,NULL,"brNDC");
      leg->AddEntry(histo1D[(vars[i_vars]+"_sample1").c_str()],label_sample1.c_str(),"f");
      leg->AddEntry(histo1D[(vars[i_vars]+"_sample2").c_str()],label_sample2.c_str(),"f");
      leg->SetFillColor(0);
      leg->Draw();




      double titleoffsety=0.2;
      double titlesizex=0.17;
      double titlesizey=0.1;
      double labelsizex=0.14;
    
      c1->cd();  


      TPad* canvas_2 = new TPad("canvas_2", "canvas_2",0,0.,1.0,0.32);
      canvas_2->Draw();
      canvas_2->cd();
      gPad->SetBottomMargin(0.375);
      gPad->SetGridy();
    
      histo_ratio->SetMarkerStyle(20);
      histo_ratio->SetMarkerSize(0.75);
      histo_ratio->SetLineWidth(2);
   
      histo_ratio->GetYaxis()->SetTitle("Blue/Red");
      histo_ratio->GetXaxis()->SetTitle((vars[i_vars]).c_str());
      histo_ratio->GetYaxis()->SetNdivisions( 505 );

      histo_ratio->SetMarkerStyle(20);
      histo_ratio->GetXaxis()->SetLabelSize( labelsizex);
      histo_ratio->GetXaxis()->SetTitleSize( titlesizex );
      histo_ratio->GetYaxis()->SetTitleSize(titlesizey);
      histo_ratio->GetYaxis()->SetTitleOffset(titleoffsety);
      histo_ratio->SetMarkerSize(1.);
      histo_ratio->GetYaxis()->SetNdivisions(5);
      histo_ratio->SetTitle("");

      histo_ratio->SetMinimum(0.7);
      histo_ratio->SetMaximum(1.3);
      histo_ratio->Draw("P");

      if(vars[i_vars] == "I_nJets")
      {
          for(int iBin = 0; iBin < histo_ratio->GetNbinsX()+1; iBin++)
          {
              cout << "else if (n_jets == " << iBin-1 << ") weight = " << histo_ratio->GetBinContent(iBin) << ";" << endl;
          }
      }
      c1->cd();  

      c1->SaveAs((vars[i_vars]+"_Norm.png").c_str());

    delete c1;
  }

  fout->Write();
  

}
