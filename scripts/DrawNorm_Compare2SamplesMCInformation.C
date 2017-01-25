map<string,TH1F*> histo1D;


void DrawNorm_Compare2SamplesMCInformation()
{
  TFile *sample1= new TFile("../MC_ComparisonOutput/MC_Comparison_NP_overlay_ST_tHToBB_1L_Kappa_hct-Official.root","READ");
  std::string label_sample1 = "ST-hct Official";
  TFile *sample2= new TFile("../MC_ComparisonOutput/MC_Comparison_NP_overlay_ST_tHToBB_1L_Kappa_hct-Private.root","READ");
  std::string label_sample2 = "ST-hct Private";

/*  string treename = "eventTree/MCParticles";
  TTree *tree_sample1 = (TTree*)sample1->Get(treename.c_str());
  TTree *tree_sample2 = (TTree*)sample2->Get(treename.c_str());
*/
  string conidition = "";

  vector<std::string> vars;
  int nbins = 100;
  vector<double> xmin;
  vector<double> xmax;



        vars.push_back("eta_Higgs");
        xmin.push_back(-3);
        xmax.push_back(3);
        
        vars.push_back("eta_Top");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("eta_AntiTop");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("phi_Higgs");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("phi_Top");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("phi_AntiTop");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("pt_Higgs");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("pt_Top");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("pt_AntiTop");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("DeltaR_Higgs_Top");
        xmin.push_back(2);
        xmax.push_back(5.5);

        vars.push_back("DeltaR_Higgs_AntiTop");
        xmin.push_back(2);
        xmax.push_back(5.5);
/*
        vars.push_back("eta_u_TopMother");
        xmin.push_back(-3);
        xmax.push_back(3);
*/
        vars.push_back("eta_b_TopMother");
        xmin.push_back(-3);
        xmax.push_back(3);
/*
        vars.push_back("eta_c_TopMother");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("phi_u_TopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);
*/
        vars.push_back("phi_b_TopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);
/*
        vars.push_back("phi_c_TopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("pt_u_TopMother");
        xmin.push_back(0);
        xmax.push_back(500);
*/
        vars.push_back("pt_b_TopMother");
        xmin.push_back(0);
        xmax.push_back(500);
/*
        vars.push_back("pt_c_TopMother");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("eta_u_AntiTopMother");
        xmin.push_back(-3);
        xmax.push_back(3);
*/
        vars.push_back("eta_b_AntiTopMother");
        xmin.push_back(-3);
        xmax.push_back(3);
/*
        vars.push_back("eta_c_AntiTopMother");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("phi_u_AntiTopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);
*/
        vars.push_back("phi_b_AntiTopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);
/*
        vars.push_back("phi_c_AntiTopMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("pt_u_AntiTopMother");
        xmin.push_back(0);
        xmax.push_back(500);
*/
        vars.push_back("pt_b_AntiTopMother");
        xmin.push_back(0);
        xmax.push_back(500);
/*
        vars.push_back("pt_c_AntiTopMother");
        xmin.push_back(0);
        xmax.push_back(500);
*/
        vars.push_back("pt_b_HiggsMother");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("eta_b_HiggsMother");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("phi_b_HiggsMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("pt_antib_HiggsMother");
        xmin.push_back(0);
        xmax.push_back(500);

        vars.push_back("eta_antib_HiggsMother");
        xmin.push_back(-3);
        xmax.push_back(3);

        vars.push_back("phi_antib_HiggsMother");
        xmin.push_back(-3.2);
        xmax.push_back(3.2);

        vars.push_back("DeltaR_bFromHiggs_Higgs");
        xmin.push_back(0.);
        xmax.push_back(5.5);

        vars.push_back("DeltaR_antibFromHiggs_Higgs");
        xmin.push_back(0.);
        xmax.push_back(5.5);

        vars.push_back("DeltaR_b_antib_FromHiggs");
        xmin.push_back(0.);
        xmax.push_back(5.5);

  TFile *fout = new TFile("NormHistos_Comp2Samples_MCinfo.root","RECREATE");


  for(int i_vars = 0; i_vars < vars.size(); i_vars++)
  {   
      histo1D[(vars[i_vars]+"_sample1").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample1").c_str(),vars[i_vars].c_str(),nbins,xmin[i_vars],xmax[i_vars]);
      histo1D[(vars[i_vars]+"_sample2").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample2").c_str(),vars[i_vars].c_str(),nbins,xmin[i_vars],xmax[i_vars]);
      gStyle->SetOptStat(kFALSE);

      if(!histo1D[(vars[i_vars]+"_sample1").c_str()] || !histo1D[(vars[i_vars]+"_sample2").c_str()])
      {
          cout << "Input histo doesn't exist" << endl;
          continue;
      }

      histo1D[(vars[i_vars]+"_sample1").c_str()] = (TH1F*) sample1->Get(vars[i_vars].c_str());
      histo1D[(vars[i_vars]+"_sample2").c_str()] = (TH1F*) sample2->Get(vars[i_vars].c_str());

      Double_t norm_sample1 = 1;
      Double_t scale_sample1 = norm_sample1/(histo1D[(vars[i_vars]+"_sample1").c_str()]->Integral());
      histo1D[(vars[i_vars]+"_sample1").c_str()]->Scale(scale_sample1);
      Double_t norm_sample2 = 1;
      Double_t scale_sample2 = norm_sample2/(histo1D[(vars[i_vars]+"_sample2")]->Integral());
      histo1D[(vars[i_vars]+"_sample2")]->Scale(scale_sample2);

      float maxhist= histo1D[(vars[i_vars]+"_sample1").c_str()]->GetMaximum();
      if (histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum() > maxhist) maxhist = histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum();
      
      maxhist = maxhist*1.5;
     
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

      TLegend* leg =  new TLegend(0.75,0.65,0.9,0.9,NULL,"brNDC");
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

      c1->cd();  

      c1->SaveAs(("plots_MCInfoComp/"+vars[i_vars]+"_Norm.png").c_str());
      c1->SaveAs(("plots_MCInfoComp/"+vars[i_vars]+"_Norm.eps").c_str());
      c1->SaveAs(("plots_MCInfoComp/"+vars[i_vars]+"_Norm.pdf").c_str());

    delete c1;
  }

  fout->Write();
  

}
