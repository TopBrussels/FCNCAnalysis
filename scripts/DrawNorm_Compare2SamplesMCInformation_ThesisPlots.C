map<string,TH1F*> histo1D;

string intToStr (int number);
string intToStr (int number)
{
  	ostringstream buff;
  	buff<<number;
  	return buff.str();
}


void DrawNorm_Compare2SamplesMCInformation_ThesisPlots()
{
  TFile *sample1_= new TFile("../NtuplerOutput/Ntuples_Mu/Ntuples_2_7_2017/FCNC_1L3B__Run2_TopTree_Study_TTJets_powheg.root","READ");
  TTree* sample1     = (TTree*) sample1_->Get("ObjectVarsTree");
  std::string label_sample1 = "powheg";
  TFile *sample2_= new TFile("../NtuplerOutput/Ntuples_Mu/Ntuples_2_7_2017/FCNC_1L3B__Run2_TopTree_Study_TTJets_FXFX.root","READ");
  TTree* sample2     = (TTree*) sample2_->Get("ObjectVarsTree");
  std::string label_sample2 = "aMC@NLO";

/*  string treename = "eventTree/MCParticles";
  TTree *tree_sample1 = (TTree*)sample1->Get(treename.c_str());
  TTree *tree_sample2 = (TTree*)sample2->Get(treename.c_str());
*/
  string conidition = "";

  vector<std::string> vars;
  std::string xaxis;
  vector<std::string> yaxis;
  int nbins = 50;
  vector<double> xmin;
  vector<double> xmax;


/*
        vars.push_back("mW");
        xmin.push_back(40);
        xmax.push_back(120);
        xaxis = "m(W_{had}) [GeV]";
   */     
        
        vars.push_back("pT_ttbar");
        xmin.push_back(240);
        xmax.push_back(950);
        xaxis = "p_{T}(t #bar{t}) [GeV]";

cout << "DEBUG 1" << endl;
  TFile *fout = new TFile("NormHistos_Comp2Samples_MCinfo_thesis.root","RECREATE");


  for(int i_vars = 0; i_vars < vars.size(); i_vars++)
  {   
      histo1D[(vars[i_vars]+"_sample1").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample1").c_str(),vars[i_vars].c_str(),nbins,xmin[i_vars],xmax[i_vars]);
      histo1D[(vars[i_vars]+"_sample2").c_str()] = new TH1F(("h_"+vars[i_vars]+"_sample2").c_str(),vars[i_vars].c_str(),nbins,xmin[i_vars],xmax[i_vars]);
      gStyle->SetOptStat(kFALSE);
cout << "DEBUG 2" << endl;
/*
      TDirectory* subdir1 = (TDirectory*) sample1->Get("ObjectVarsTree");
      subdir1->cd();
cout << "DEBUG 2a" << endl;
*/
      /*histo1D[(vars[i_vars]+"_sample1").c_str()] = (TH1F*) */sample1->Draw((vars[i_vars]+">>h_temp1").c_str(),(vars[i_vars]+">"+intToStr(xmin[i_vars])+"&&"+vars[i_vars]+"<"+intToStr(xmax[i_vars])).c_str());
      histo1D[(vars[i_vars]+"_sample1").c_str()] = (TH1F*)gDirectory->Get("h_temp1");
      histo1D[(vars[i_vars]+"_sample1").c_str()]->Rebin(2);
/*cout << "DEBUG 2b" << endl;
      TDirectory* subdir2 = (TDirectory*) sample2->Get("ObjectVarsTree");
      subdir2->cd();
*/      /*histo1D[(vars[i_vars]+"_sample2").c_str()] = (TH1F*) */sample2->Draw((vars[i_vars]+">>h_temp2").c_str(),(vars[i_vars]+">"+intToStr(xmin[i_vars])+"&&"+vars[i_vars]+"<"+intToStr(xmax[i_vars])).c_str());
      histo1D[(vars[i_vars]+"_sample2").c_str()] = (TH1F*)gDirectory->Get("h_temp2");
      histo1D[(vars[i_vars]+"_sample2").c_str()]->Rebin(2);
cout << "DEBUG 3" << endl;

      if(!histo1D[(vars[i_vars]+"_sample1").c_str()] || !histo1D[(vars[i_vars]+"_sample2").c_str()])
      {
          cout << "Input histo doesn't exist" << endl;
          continue;
      }
      Double_t norm_sample1 = 1;
      Double_t scale_sample1 = norm_sample1/(histo1D[(vars[i_vars]+"_sample1").c_str()]->Integral());
      histo1D[(vars[i_vars]+"_sample1").c_str()]->Scale(scale_sample1);
      Double_t norm_sample2 = 1;
      Double_t scale_sample2 = norm_sample2/(histo1D[(vars[i_vars]+"_sample2")]->Integral());
      histo1D[(vars[i_vars]+"_sample2")]->Scale(scale_sample2);

      float maxhist= histo1D[(vars[i_vars]+"_sample1").c_str()]->GetMaximum();
      if (histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum() > maxhist) maxhist = histo1D[(vars[i_vars]+"_sample2").c_str()]->GetMaximum();
 cout << "DEBUG 4" << endl;
     
      maxhist = maxhist*1.2;
     
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetMaximum(maxhist);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetMaximum(maxhist);

      double titleoffsety=0.2;
      double titlesizex=0.15;
      double titlesizey=0.1;
      double labelsizex=0.14;
      
      histo1D[(vars[i_vars]+"_sample1").c_str()]->GetXaxis()->SetTitle(xaxis.c_str());
      histo1D[(vars[i_vars]+"_sample1").c_str()]->GetYaxis()->SetTitle("normalized");
      histo1D[(vars[i_vars]+"_sample1").c_str()]->GetYaxis()->SetTitleSize(0.05);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->GetYaxis()->SetTitleOffset(1.);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetTitle("");
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetLineColor(4);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetLineColor(2);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetLineWidth(2);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetLineWidth(2);
/*      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetFillColor(4);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetFillColor(2);
      histo1D[(vars[i_vars]+"_sample1").c_str()]->SetFillStyle(3004);
      histo1D[(vars[i_vars]+"_sample2").c_str()]->SetFillStyle(3005);
*/
cout << "DEBUG 5" << endl;

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
      histo_ratio->GetXaxis()->SetTitle(xaxis.c_str());
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
