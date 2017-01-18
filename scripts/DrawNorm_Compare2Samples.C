map<string,TH1F*> histo1D;


void DrawNorm_Compare2Samples()
{
  TFile *sample1= new TFile("../Merged/Ntuples_El/Ntuples_12_1_2017/FCNC_1L3B__Run2_TopTree_Study_NP_overlay_ST_tHToBB_1L_Kappa_hct.root","READ");
  std::string label_sample1 = "ST-hct FullSim";
  TFile *sample2= new TFile("/user/kderoove/FCNC/TopTreeFramework_Run2/CMSSW_7_6_5/src/TopBrussels/FCNCAnalysis/Merged/Ntuples_Mu/Ntuples_20_7_2016/FCNC_1L3B__Run2_TopTree_Study_NP_overlay_ST_tHToBB_1L_Kappa_hct.root","READ");
  std::string label_sample2 = "ST-hct FastSim";

  string treename = "ObjectVarsTree";
  TTree *tree_sample1 = (TTree*)sample1->Get(treename.c_str());
  TTree *tree_sample2 = (TTree*)sample2->Get(treename.c_str());

  vector<std::string> vars;
  vector<int> nbins;
  vector<double> xmin;
  vector<double> xmax;



  //Defining the variables we want to plot with the nbins, xmin and xmax
  vars.push_back("I_nvtx");
  nbins.push_back(51);
  xmin.push_back(-0.5);
  xmax.push_back(50.5);
  
  vars.push_back("I_npu");
  nbins.push_back(51);
  xmin.push_back(-0.5);
  xmax.push_back(50.5);
  
  vars.push_back("I_nJets_CSVL");
  nbins.push_back(11);
  xmin.push_back(-0.5);
  xmax.push_back(10.5);
  
  vars.push_back("I_nJets_CSVM");
  nbins.push_back(11);
  xmin.push_back(-0.5);
  xmax.push_back(10.5);
  
  vars.push_back("I_nJets_CSVT");
  nbins.push_back(11);
  xmin.push_back(-0.5);
  xmax.push_back(10.5);
  
  vars.push_back("I_nJets");
  nbins.push_back(11);
  xmin.push_back(-0.5);
  xmax.push_back(10.5); 
  
  vars.push_back("pt_lepton");
  nbins.push_back(50);
  xmin.push_back(20.);
  xmax.push_back(300.); 
  
  vars.push_back("eta_lepton");
  nbins.push_back(50);
  xmin.push_back(-2.5);
  xmax.push_back(2.5); 
  
  vars.push_back("phi_lepton");
  nbins.push_back(50);
  xmin.push_back(-3.2);
  xmax.push_back(3.2); 
  
//  vars.push_back("I_LepCharge");
//  nbins.push_back(3);
//  xmin.push_back(-1.5);
//  xmax.push_back(1.5); 
  
  vars.push_back("pt_jet");
  nbins.push_back(50);
  xmin.push_back(20.);
  xmax.push_back(300.);
  
  vars.push_back("eta_jet");
  nbins.push_back(50);
  xmin.push_back(-2.5);
  xmax.push_back(2.5); 
  
  vars.push_back("phi_jet");
  nbins.push_back(50);
  xmin.push_back(-3.2);
  xmax.push_back(3.2);
  
  vars.push_back("CSVv2");
  nbins.push_back(50);
  xmin.push_back(0.);
  xmax.push_back(1.);
  
  vars.push_back("cMVA");
  nbins.push_back(50);
  xmin.push_back(-1.);
  xmax.push_back(1.);

  vars.push_back("cdiscCvsL_jet");
  nbins.push_back(50);
  xmin.push_back(-1.);
  xmax.push_back(1.);

  vars.push_back("cdiscCvsB_jet");
  nbins.push_back(50);
  xmin.push_back(-1.);
  xmax.push_back(1.);
  
  vars.push_back("jet_matchedMC_pdgID");
  nbins.push_back(60);
  xmin.push_back(-30.5);
  xmax.push_back(29.5);

  vars.push_back("jet_matchedMC_motherpdgID");
  nbins.push_back(60);
  xmin.push_back(-30.5);
  xmax.push_back(29.5);

  vars.push_back("jet_matchedMC_grannypdgID");
  nbins.push_back(60);
  xmin.push_back(-30.5);
  xmax.push_back(29.5);
  
/*
  vars.push_back("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
  vars.push_back("HiggsMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
  vars.push_back("HiggsEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (Higgs)","Events", category);
  vars.push_back("HiggsEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (Higgs)","Events", category);
  vars.push_back("TopLepMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
  vars.push_back("TopLepMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
  vars.push_back("TopLepPt_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepPt_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "Pt(LepTop)","Events", category,"GeV");
  vars.push_back("TopLepPt_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepPt_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "Pt(LepTop)","Events", category,"GeV");
  vars.push_back("TopLepEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (LepTop)","Events", category);
  vars.push_back("TopLepEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (LepTop)","Events", category);
  vars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
  vars.push_back("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
  vars.push_back("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "#Delta R(H,LepTop)","Events", category);
  vars.push_back("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "#Delta R(H,LepTop)","Events", category);
  vars.push_back("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopHadMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopHadMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(HadTop)","Events", category,"GeV");
  vars.push_back("TopLepMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
  vars.push_back("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(HadTop,LepTop)","Events", category);
  vars.push_back("TopLepBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepBJetCSVv2_TOPTOPLEPHAF"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("HiggsMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
  vars.push_back("TopLepMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
  vars.push_back("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
  vars.push_back("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(H,LepTop)","Events", category);
  vars.push_back("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
  vars.push_back("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets, ("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
*/

  TFile *fout = new TFile("NormHistos_Comp2Samples.root","RECREATE");

/*
  for(int i_vars = 0; i_vars < vars.size(); i_vars++)
  {   
      histo1D[(vars[i_vars]+"_sample1").c_str()] = new TH1F();
      histo1D[(vars[i_vars]+"_sample2").c_str()] = new TH1F();
  }
*/
//  gSystem->Load("$ROOTSYS/test/libEvent");

  for(int i_vars = 0; i_vars < vars.size(); i_vars++)
  {   
      histo1D[(vars[i_vars]+"_sample1").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample1").c_str(),vars[i_vars].c_str(),nbins[i_vars],xmin[i_vars],xmax[i_vars]);
      histo1D[(vars[i_vars]+"_sample2").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample2").c_str(),vars[i_vars].c_str(),nbins[i_vars],xmin[i_vars],xmax[i_vars]);
      gStyle->SetOptStat(kFALSE);
      tree_sample1->Draw((vars[i_vars]+">>h_"+vars[i_vars]+"_sample1").c_str());
      tree_sample2->Draw((vars[i_vars]+">>h_"+vars[i_vars]+"_sample2").c_str());


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

      c1->cd();  

      c1->SaveAs((vars[i_vars]+"_Norm.png").c_str());

    delete c1;
  }

  fout->Write();
  

}
