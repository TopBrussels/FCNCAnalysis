#include "../interface/MultiSamplePlot.h"

/*MultiSamplePlot::MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float Min, float Max, string XaxisLabel, string YaxisLabel, string Text)
 {
 Nbins_ = Nbins;
 string histoName = "";
 string histoTitle = "";
 plotName_ = PlotName;
 text_ = TString(Text);
 lumi_ = 1.;
 for(unsigned int i=0;i<datasets.size();i++){
 histoName  = PlotName+"_"+datasets[i]->Name();
 histoTitle = datasets[i]->Title();
 TH1F* h = new TH1F(histoName.c_str(),histoTitle.c_str(),Nbins,Min,Max);
 h->SetStats(false);
 h->Sumw2();
 h->GetXaxis()->SetTitle(XaxisLabel.c_str());
 h->GetYaxis()->SetTitle(YaxisLabel.c_str());
 plots_.push_back(pair<TH1F*,Dataset*>(h,datasets[i]));
 if(datasets[i]->Name().find("data") == 0 || datasets[i]->Name().find("Data") == 0 || datasets[i]->Name().find("DATA") == 0 )
 lumi_ = datasets[i]->EquivalentLumi();
 }
 Initialize();
 }*/

MultiSamplePlot::MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float Min, float Max, string XaxisLabel, string YaxisLabel, string Text, string Units)
{
  std::stringstream stream;
  string XaxisLabel_, YaxisLabel_;
  Nbins_ = Nbins;
  string histoName = "";
  string histoTitle = "";
  float binWidth = (Max - Min)/Nbins;
  if(!Units.empty()) stream << std::fixed << std::setprecision(2) << binWidth;
  else stream << std::fixed << std::setprecision(2) << binWidth;
  string sbinWidth = stream.str();
  if(!Units.empty()) XaxisLabel_ = XaxisLabel + " [" + Units + "]";
  else XaxisLabel_ = XaxisLabel;
  if(!Units.empty()) YaxisLabel_ = YaxisLabel + " / " + sbinWidth + " " + Units;
  else YaxisLabel_ = YaxisLabel + " / " + sbinWidth + Units;
  plotName_ = PlotName;
  text_ = TString(Text);
  lumi_ = 1.;
  for(unsigned int i=0;i<datasets.size();i++){
    histoName  = PlotName+"_"+datasets[i]->Name();
    histoTitle = datasets[i]->Title();
    TH1F* h = new TH1F(histoName.c_str(),histoTitle.c_str(),Nbins,Min,Max);
    h->SetStats(false);
    h->Sumw2();
    h->GetXaxis()->SetTitle(XaxisLabel_.c_str());
    h->GetYaxis()->SetTitle(YaxisLabel_.c_str());

    
    plots_.push_back(pair<TH1F*,Dataset*>(h,datasets[i]));
    if(datasets[i]->Name().find("data") == 0 || datasets[i]->Name().find("Data") == 0 || datasets[i]->Name().find("DATA") == 0 )
      lumi_ = datasets[i]->EquivalentLumi();
  }
  Initialize();
}

MultiSamplePlot::MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float* binsX, string XaxisLabel, string YaxisLabel, string Text)
{
  Nbins_ = Nbins;
  string histoName = "";
  string histoTitle = "";
  plotName_ = PlotName;
  text_ = TString(Text);
  lumi_ = 1.;
  for(unsigned int i=0;i<datasets.size();i++){
    histoName = PlotName+"_"+datasets[i]->Name();
    histoTitle = datasets[i]->Title();
    TH1F* h = new TH1F(histoName.c_str(),histoTitle.c_str(),Nbins,binsX);
    h->SetStats(false);
    h->Sumw2();
    h->GetXaxis()->SetTitle(XaxisLabel.c_str());
    h->GetYaxis()->SetTitle(YaxisLabel.c_str());

    plots_.push_back(pair<TH1F*,Dataset*>(h,datasets[i]));
    if(datasets[i]->Name().find("data") == 0 || datasets[i]->Name().find("Data") == 0 || datasets[i]->Name().find("DATA") == 0 )
      lumi_ = datasets[i]->EquivalentLumi();
  }
  Initialize();
}

MultiSamplePlot::MultiSamplePlot(vector<pair<TH1F*,Dataset*> > vec, string PlotName, string XaxisLabel, string YaxisLabel, string Text)
{
  plots_ = vec;
  Nbins_ = (vec[0].first)->GetNbinsX();
  plotName_ = PlotName;
  text_ = TString(Text);
  lumi_ = 1.;
  for(unsigned int i=0;i<vec.size();i++)
  {
    if(vec[i].second->Name().find("data") == 0 || vec[i].second->Name().find("Data") == 0 || vec[i].second->Name().find("DATA") == 0 )
      lumi_ = vec[i].second->EquivalentLumi();
    (vec[i].first)->GetXaxis()->SetTitle(XaxisLabel.c_str());
    (vec[i].first)->GetYaxis()->SetTitle(YaxisLabel.c_str());
  }
  Initialize();
}

void MultiSamplePlot::Initialize()
{
  hData_ = 0;
  hCanvas_ = 0;
  hCanvasStack_ = 0;
  hCanvasStackLogY_ = 0;
  hCanvasStackAreaNorm_ = 0;
  hCanvasStackAreaNormLogY_ = 0;
  hStack_ = 0;
  hStackAreaNorm_ = 0;
  StackErrorGraph_ = 0;
  RatioErrorGraph_ = 0;
  maxY_ = -1;
  minLogY_ = 1.;
  logYMultiplicator_ = 1000;
  maxLogY_ = 1.;
  maxLogYAreaNorm_ = 1.;
  showNumberEntries_ = true;
  //errorbandfile_ = "systematics.root";
  errorbandfile_ = "errorbands.root";
  dosystfile_ = false;
  sqrts_ = 13;
  prelim_=true;
  chan_=false;
  channel_ = " ";
  setBinLabels_ = false;
  vlabel_ = {""};
  doCutFlow_ = false;
}

MultiSamplePlot::~MultiSamplePlot()
{
  if(hStack_)      delete hStack_;
  if(hStackAreaNorm_)      delete hStackAreaNorm_;
  if(hCanvas_)     delete hCanvas_;
  if(hCanvasStack_)     delete hCanvasStack_;
  if(hCanvasStackAreaNorm_)     delete hCanvasStackAreaNorm_;
  if(hCanvasStackLogY_)     delete hCanvasStackLogY_;
  if(hCanvasStackAreaNormLogY_)     delete hCanvasStackAreaNormLogY_;
  if(hData_) delete hData_;
  if(StackErrorGraph_) delete StackErrorGraph_;
  if(RatioErrorGraph_) delete RatioErrorGraph_;
}

void MultiSamplePlot::AddDataHisto(TH1F* histo)
{
  hData_ = (TH1F*) histo->Clone();
  if(doCutFlow_ ){
    for(int iBin = 1; iBin < histo->GetNbinsX()+1; iBin++){
      hData_->SetBinError(iBin, std::sqrt(histo->GetBinContent(iBin)));
      //cout << "bin " << iBin << " error " << std::sqrt(histo->GetBinContent(iBin)) << " content " << histo->GetBinContent(iBin) << endl;
    }
  }
}

void MultiSamplePlot::Fill(float value, Dataset* data, bool scale, float Lumi)
{
  if(scale)
  {
    for(unsigned int i=0;i<plots_.size();i++)
    {
      if(data->Name()==plots_[i].second->Name()) plots_[i].first->Fill(value, data->NormFactor()*Lumi);
      //			if(data->Name()==plots_[i].second->Name()) cout << "NormFactor for " << data->Name() << " = " <<  data->NormFactor()*Lumi << endl;
    }
  }
  else
  {
    for(unsigned int i=0;i<plots_.size();i++)
    {
      if(data->Name()==plots_[i].second->Name()) plots_[i].first->Fill(value);
    }
  }
  
  
}

void MultiSamplePlot::Draw(string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal)
{
  vector<TH1F*> normalizedhistos;
  vector<Dataset*> mergeddatasets;
  
  this->setTDRStyle();
  
  if(RatioType==0 && addRatioErrorBand)
  {
    cout<<"[MultiSamplePlot::Draw] WARNING: addRatioErrorBand set to true but RatioType to 0 (i.e. no ratio) -> resetting addRatioErrorBand to false now."<<endl;
    addRatioErrorBand = false;
  }
  
  TFile* tmpInfile = 0;
  if(addRatioErrorBand || addErrorBand)
  {
    tmpInfile = new TFile(errorbandfile_.c_str(),"READ");
    if(tmpInfile->IsZombie())
    {
      cout<<"[MultiSamplePlot::Draw] WARNING: '"<<errorbandfile_<<"' is no valid error band file -> will not draw error bands. Please set a valid file with the MultiSamplePlot::setErrorBandFile method."<<endl;
      addErrorBand = false;
      addRatioErrorBand = false;
    }
  }
  
  TH1F* hErrorPlus = 0;
  TH1F* hErrorMinus = 0;
  TH1F* hNominal = 0;
  if(addRatioErrorBand || addErrorBand)
  {
    TDirectory* subdir = (TDirectory*) tmpInfile->Get(("MultiSamplePlot_"+plotName_).c_str());
    if(subdir)
    {
      
      
      hErrorPlus = (TH1F*) subdir->Get("Plus")->Clone();
      hErrorMinus = (TH1F*) subdir->Get("Minus")->Clone();
      hNominal = (TH1F*) subdir->Get("Nominal")->Clone();
      //cout << "Integrals for plot " << label << " ->   hErrorPlus: " << hErrorPlus->Integral(1,hErrorPlus->GetNbinsX()) << ",  hNominal: " << hNominal->Integral(1,hNominal->GetNbinsX()) << ",  hErrorMinus: " << hErrorMinus->Integral(1,hErrorMinus->GetNbinsX()) << endl;
    }
    else
    {
      cout << "[MultiSamplePlot::Draw] ERROR: errorband file '"<<errorbandfile_<<"' has unknown structure! (Histograms should be in directories with MultiSamplePlot names) -> will not draw error bands." <<endl;
      addErrorBand = false;
      addRatioErrorBand = false;
    }
  }
  
  
  //loop over the processes and create one new histogram per set of processes to merge. Whether or not processes will be merged, will be based on their title in the xml config; same title -> same merged histogram.
  int dataPlotID=-1;
  unsigned int nDataPlots = 0;
  vector<pair<TH1F*,Dataset*> > SetsOfMergedProcesses; //the dataset will correspond to a set of merged processes (e.g. "TTbar+jets", "W+jets", ...)
  for(unsigned int i=0;i<plots_.size();i++)
  {
    if(plots_[i].second->Name().find("data") == 0 || plots_[i].second->Name().find("Data") == 0 || plots_[i].second->Name().find("DATA") == 0 )
    {
      dataPlotID = i;
      nDataPlots++;
    }
    
    string datasetTitle = plots_[i].second->Title();
    bool histocreated = false;
    //check if the histogram corresponding to this title is already created
    for(unsigned int j=0;j<SetsOfMergedProcesses.size();j++)
    {
      if(SetsOfMergedProcesses[j].second->Title() == datasetTitle) histocreated = true;
    }
    
    (plots_[i].first)->SetLineColor((plots_[i].second)->Color()); //especially relevant for the overlay processes, the stacked histograms get their fill color later on. However, all seperate processes are also saved in MSPlot, so these will get the color specified in the xml too.
    (plots_[i].first)->SetMarkerColor((plots_[i].second)->Color());
    (plots_[i].first)->SetLineStyle((plots_[i].second)->LineStyle());
    (plots_[i].first)->SetLineWidth((plots_[i].second)->LineWidth());
    //making sure that the overflow is transferred to the last 'visible' bin; analogously for underflow...
    (plots_[i].first)->SetBinContent(Nbins_,(plots_[i].first)->GetBinContent(Nbins_)+(plots_[i].first)->GetBinContent(Nbins_+1));
    (plots_[i].first)->SetBinContent(Nbins_+1,0);
    (plots_[i].first)->SetBinContent(1,(plots_[i].first)->GetBinContent(0)+(plots_[i].first)->GetBinContent(1));
    (plots_[i].first)->SetBinContent(0,0);
    
    if(!histocreated)
    {
      TH1F* MergedProcesses = (TH1F*) plots_[i].first->Clone();
      string datasetName = plots_[i].second->Title()+"_Merged";
      if(plots_[i].second->Name().find("NP_overlay_")==0) datasetName = "NP_overlay_"+datasetName;
      else if(plots_[i].second->Name().find("NP_")==0) datasetName = "NP_"+datasetName;
      Dataset* dataSet = new Dataset(datasetName,datasetTitle,true,plots_[i].second->Color(),plots_[i].second->LineStyle(),plots_[i].second->LineWidth(),1.,1.);
      MergedProcesses->GetXaxis()->SetTitle(plots_[i].first->GetXaxis()->GetTitle());
      MergedProcesses->GetYaxis()->SetTitle(plots_[i].first->GetYaxis()->GetTitle());
      SetsOfMergedProcesses.push_back(pair<TH1F*,Dataset*>(MergedProcesses,dataSet));
    }
    else
    {
      //add the process to the histogram corresponding to its title
      for(unsigned int j=0;j<SetsOfMergedProcesses.size();j++)
      {
        if(SetsOfMergedProcesses[j].second->Title() == datasetTitle)
        {
          SetsOfMergedProcesses[j].first->Add(plots_[i].first);
        }
      }
    }
  }
  
  if(nDataPlots>1) cout<<"[MultiSamplePlot::Draw] ERROR: you are not allowed to superimpose two data histograms on the same plot!"<<endl;
  else if(nDataPlots==0)
  {
    //cout<<"[MultiSamplePlot::Draw] WARNING: not running on data, ratio plots will not be drawn
    RatioType=0;
    addRatioErrorBand = false;
  }
  
  //NOTE: in this stage, the SetsOfMergedProcesses vector is basically the same as the plots_ vector, but with the processes merged that have the same title in the xml config.
  for(unsigned int i=0;i<SetsOfMergedProcesses.size();i++)
  {
    TH1F* histo = (TH1F*) SetsOfMergedProcesses[i].first->Clone();
    Dataset* mergeddataSet = SetsOfMergedProcesses[i].second;
    TH1F* normalizedclone = (TH1F*) histo->Clone();
    normalizedclone->Scale(1/normalizedclone->Integral(1,normalizedclone->GetNbinsX()+1));
    //normalizedclone->SetFillColor(10);
    normalizedclone->SetFillStyle(1001);
    normalizedclone->GetYaxis()->SetTitle( ( "Normalized " + (string) histo->GetYaxis()->GetTitle() ).c_str() );
    normalizedhistos.push_back(normalizedclone);
    mergeddatasets.push_back(mergeddataSet);
  }
  
  //Make the first canvas; to plot distributions normalized to unity
  string name = "Canvas_"+label;
  hCanvas_ = TCanvasCreator(normalizedhistos, mergeddatasets, leg_, string("l"), name);
  
  
  //calculating data and MC integrals
  float integralData = 0;
  float integralSM = 0;
  if ( dataPlotID != -1 )
  {
    integralData = plots_[dataPlotID].first->Integral(0,plots_[dataPlotID].first->GetNbinsX()+1); //first argument 0 means you include the underflow bin too
    for(unsigned int i=0;i<plots_.size();i++)
    {
      if ((dataPlotID != (int)i) && (plots_[i].second->Name().find("NP_")!=0))
        integralSM += plots_[i].first->Integral(0,plots_[i].first->GetNbinsX());
    }
    
  }
  
  
  //For histograms to be stacked
  vector<TH1F*> histosForHStack, histosForHStackAreaNorm;
  TH1F* totalSM = 0;
  TH1F* totalSMAreaNorm = 0;
  vector<Dataset*> datasetsForHStack;
  //For histograms (e.g. new-physics signal) to be overlayed
  vector<TH1F*> histosForOverlay, histosForOverlayAreaNorm;
  vector<Dataset*> datasetsForOverlay;
  
  //Loop from back to front because of the way THStack is used later on (related to order as appearing in xml config compared to order in histogram stack)...
  for(int i=SetsOfMergedProcesses.size()-1;i>=0;i--)
  {
    Dataset* mergeddataSet = SetsOfMergedProcesses[i].second;
    if(SetsOfMergedProcesses[i].second->Name().find("data") == 0 || SetsOfMergedProcesses[i].second->Name().find("Data") == 0 || SetsOfMergedProcesses[i].second->Name().find("DATA") == 0 )
    {
      AddDataHisto( SetsOfMergedProcesses[i].first );
    }
    else
    {
      TH1F* clone = (TH1F*) SetsOfMergedProcesses[i].first->Clone();
      if(mergeddataSet->Name().find("NP_")!=0)
      {
        histosForHStack.push_back(clone);
        datasetsForHStack.push_back(mergeddataSet);
        if( ! totalSM ) totalSM = (TH1F*) clone->Clone();
        else totalSM->Add( clone );
      }
      else if (mergeddatasets[i]->Name().find("NP_overlay_")==0)
      {
        histosForOverlay.push_back(clone);
        datasetsForOverlay.push_back(mergeddatasets[i]);
      }
      
      TH1F* cloneAreaNorm = (TH1F*) SetsOfMergedProcesses[i].first->Clone();
      if (dataPlotID != -1  && mergeddataSet->Name().find("NP_")!=0)
      {
//        cloneAreaNorm->Scale(1.);
        cloneAreaNorm->Scale(integralData/integralSM);
      }
      else if (dataPlotID != -1 && mergeddataSet->Name().find("NP_")==0)
      {
//        cloneAreaNorm->Scale(integralSM/cloneAreaNorm->Integral());
        cloneAreaNorm->Scale(integralData/cloneAreaNorm->Integral());
      }
      if(mergeddataSet->Name().find("NP_")!=0)
      {
        histosForHStackAreaNorm.push_back(cloneAreaNorm);
        if( ! totalSMAreaNorm ) totalSMAreaNorm = (TH1F*) cloneAreaNorm->Clone();
        else totalSMAreaNorm->Add( cloneAreaNorm );
      }
      else if (mergeddataSet->Name().find("NP_overlay_")==0)
      {
        histosForOverlayAreaNorm.push_back(cloneAreaNorm);
      }
    }
  }
  
  
  leg_ = new TLegend(0.7,0.55,1.,0.90);
  leg_->SetFillColorAlpha(0,0.0);
  leg_->SetTextFont(42);
  leg_->SetTextSize(0.04);
  leg_->SetLineColor(1);
  leg_->SetLineWidth(1);
  leg_->SetLineStyle(0);
  leg_->SetBorderSize(0);
  if( ! showNumberEntries_ ) leg_->SetX1(0.76);
  
  //a second legend is needed because otherwise it will add the datasets two times to the same legend...
  legAreaNorm_ = new TLegend(0.7,0.55,1.,0.90);
  legAreaNorm_->SetFillColorAlpha(0,0.0);
  legAreaNorm_->SetTextFont(42);
  legAreaNorm_->SetTextSize(0.04);
  legAreaNorm_->SetLineColor(1);
  legAreaNorm_->SetLineWidth(1);
  legAreaNorm_->SetLineStyle(0);
  legAreaNorm_->SetBorderSize(0);
  if( ! showNumberEntries_ ) legAreaNorm_->SetX1(0.76);
  
  //if running over data, add this one first in the legend
  if(hData_)
  {
    char legDataTitle[100];
    //hData_->SetTitle("Observed");
    hData_->SetTitle("Data");
    if (showNumberEntries_)
      sprintf (legDataTitle, "%s (%.1f entries)", hData_->GetTitle(), hData_->Integral(0,hData_->GetNbinsX()+1) );
    else
      sprintf (legDataTitle, "%s", hData_->GetTitle());
    leg_->AddEntry(hData_, legDataTitle ,"L E P");
    legAreaNorm_->AddEntry(hData_, legDataTitle ,"L E P");
  }
  
  
  if( histosForHStack.size() > 0 )
    hStack_ = THStackCreator(histosForHStack, datasetsForHStack, leg_, showNumberEntries_);
  name = "Stack_"+label;
  if(hStack_)
  {
    hStack_->SetName(name.c_str());
    hStack_->SetTitle(""); //hStackAreaNorm_->SetTitle(name.c_str()) gives problem: title and the CMS text get mixed on the canvas...
  }
  
  if( histosForHStackAreaNorm.size() > 0 )
    hStackAreaNorm_ = THStackCreator(histosForHStackAreaNorm, datasetsForHStack, legAreaNorm_, showNumberEntries_);
  name = "StackAreaNorm_"+label;
  if(hStackAreaNorm_)
  {
    hStackAreaNorm_->SetName(name.c_str());
    hStackAreaNorm_->SetTitle(""); //hStackAreaNorm_->SetTitle(name.c_str()) gives problem: title and the CMS text get mixed on the canvas...
  }
  
  
  //just adding the overlay histograms to the legend and calculating if we have to change the y-axis range
  float ymaxoverlay = 0, ymaxoverlayAreaNorm = 0;
  for(unsigned int i=0;i<histosForOverlay.size();i++)
  {
    char legTitle[100];
    if (showNumberEntries_)
    {
      float scalednentries = scaleNPSignal * histosForOverlay[i]->Integral(0,histosForOverlay[i]->GetNbinsX()+1);
      sprintf (legTitle, "%s (%.1f entries) (X %i)", ( datasetsForOverlay[i]->Title() ).c_str(), scalednentries, scaleNPSignal);
    }
    else
      sprintf (legTitle, "%s", ( datasetsForOverlay[i]->Title() ).c_str(), scaleNPSignal);
    leg_->AddEntry(histosForOverlay[i],legTitle,"F");
    legAreaNorm_->AddEntry(histosForOverlay[i],legTitle,"F");
    if(histosForOverlay[i]->GetMaximum()*scaleNPSignal > ymaxoverlay) ymaxoverlay = histosForOverlay[i]->GetMaximum()*scaleNPSignal;
    if(histosForOverlayAreaNorm[i]->GetMaximum() > ymaxoverlayAreaNorm) ymaxoverlayAreaNorm = histosForOverlayAreaNorm[i]->GetMaximum()*scaleNPSignal;
  }
  
  
  //regular 'luminosity normalized'
  name = "CanvasStack_"+label;
  hCanvasStack_ = new TCanvas(name.c_str(), name.c_str(), 1000, 700);
  hCanvasStack_->SetRightMargin(0.25);
  hCanvasStackLogY_ = new TCanvas((name+"_LogY").c_str(), (name+"_LogY").c_str(), 1000, 700);
  hCanvasStackLogY_->SetRightMargin(0.25);
  float ymax = 0;
  if(hStack_)
  {
    if(RatioType>0)
    {
      hCanvasStack_->SetCanvasSize(1000, 800);
      hCanvasStack_->SetBottomMargin(0.3);
      hCanvasStackLogY_->SetCanvasSize(1000, 800);
      hCanvasStackLogY_->SetBottomMargin(0.3);
    }
    
    DrawStackedPlot(hCanvasStack_,hCanvasStackLogY_,hStack_,histosForOverlay,scaleNPSignal,histosForHStack[0]->GetXaxis()->GetTitle(),histosForHStack[0]->GetYaxis()->GetTitle(),RatioType);
    ymax = 1.55*hStack_->GetMaximum();
    if(ymaxoverlay > ymax) ymax = ymaxoverlay;
    maxLogY_ = logYMultiplicator_ * ymax;
  }
  
  //'area normalized'
  name = "CanvasStackAreaNorm_"+label;
  hCanvasStackAreaNorm_ = new TCanvas(name.c_str(), name.c_str(), 1000, 700);
  hCanvasStackAreaNorm_->SetRightMargin(0.25);
  hCanvasStackAreaNormLogY_ = new TCanvas((name+"_LogY").c_str(), (name+"_LogY").c_str(), 1000, 700);
  hCanvasStackAreaNormLogY_->SetRightMargin(0.25);
  float ymaxAreaNorm = 0;
  if(hStackAreaNorm_)
  {
    if(RatioType>0)
    {
      hCanvasStackAreaNorm_->SetCanvasSize(1000, 800);
      hCanvasStackAreaNorm_->SetBottomMargin(0.3);
      hCanvasStackAreaNormLogY_->SetCanvasSize(1000, 800);
      hCanvasStackAreaNormLogY_->SetBottomMargin(0.3);
    }
    
    DrawStackedPlot(hCanvasStackAreaNorm_,hCanvasStackAreaNormLogY_,hStackAreaNorm_,histosForOverlayAreaNorm,scaleNPSignal,histosForHStackAreaNorm[0]->GetXaxis()->GetTitle(),histosForHStackAreaNorm[0]->GetYaxis()->GetTitle(),RatioType);
    ymaxAreaNorm = 1.3*hStackAreaNorm_->GetMaximum();
    if(ymaxoverlayAreaNorm > ymaxAreaNorm) ymaxAreaNorm = ymaxoverlayAreaNorm;
    maxLogYAreaNorm_ = logYMultiplicator_ * ymaxAreaNorm;
  }
  
  
  //Now superimpose data to the plots + draw ratio plots + draw error bands
  TH1F* ratio = 0;
  TPad* pad = 0;
  TPad* padLogY = 0;
  if(hData_)
  {
    hData_->SetLineColor(1);
    hData_->SetMarkerStyle(20);
    
    ratio = (TH1F*) hData_->Clone();
    
    if(hStack_)
    {
      if(RatioType==1)
      {
        ratio->Divide(totalSM); //ratio = data/MC
        ratio->GetYaxis()->SetTitle("Data/MC");
        ratio->SetMaximum(1.6);
        ratio->SetMinimum(0.4);
        for(int iBin = 1; iBin < ratio->GetNbinsX()+1; iBin++){
            double a = std::sqrt(hData_->GetBinContent(iBin));
            double b = totalSM->GetBinContent(iBin);
            ratio->SetBinError(iBin, a/b);
            //cout << "bin " << iBin << " error " << binerror << " content " << ratio->GetBinContent(iBin) << endl;
        }
        
      }
      else if(RatioType==2)
      {
        ratio->Add(totalSM,-1);
        ratio->Divide(totalSM); //ratio = (data-MC)/MC
        ratio->GetYaxis()->SetTitle("(Data-MC)/MC");
        ratio->SetMaximum(0.5);
        ratio->SetMinimum(-0.5);
        for(int iBin = 1; iBin < ratio->GetNbinsX()+1; iBin++){
            double a = std::sqrt(hData_->GetBinContent(iBin));
            double b = totalSM->GetBinContent(iBin);
            ratio->SetBinError(iBin, a/b);
            //cout << "bin " << iBin << " error " << binerror << " content " << ratio->GetBinContent(iBin) << endl;
        }
        
      }
    }
    ratio->SetMarkerStyle(20);
    ratio->GetXaxis()->SetLabelSize(0.04);
    ratio->GetXaxis()->SetTitleSize(0.045);
    ratio->GetYaxis()->SetLabelSize(0.04);
    ratio->GetYaxis()->SetTitleSize(0.045);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->SetMarkerSize(1.);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->SetTitle("");
    if(setBinLabels_ ){
      for(int ibin = 1; ibin < vlabel_.size()+1; ibin++){
        ratio->GetXaxis()->SetBinLabel(ibin,vlabel_[ibin-1].c_str());
      }
    }
    
    pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad->SetTopMargin(0.7);
    pad->SetFillColor(0); // was 1 but could have been related to the problem of black canvases
    pad->SetFillStyle(0); // was 1 but could have been related to the problem of black canvases
    pad->SetGridy(1);
    pad->SetRightMargin(0.25);
    pad->SetLeftMargin(0.12);
    padLogY = (TPad*) pad->Clone();
  }
  
  hCanvasStack_->cd();
  if(hStack_)
  {
    if(addErrorBand)
      DrawErrorBand(totalSM,hErrorPlus,hErrorMinus,hNominal,StackErrorGraph_,ErrorBandAroundTotalInput);
    
    if(hData_)
    {
      hData_->Draw("same E");
     
      
      if(RatioType>0)
      {
        pad->Draw();
        pad->cd(0);
        ratio->Draw("e");
        
        if(addRatioErrorBand)
        {
          if(doCutFlow_) cout << "watch out! uncertainties are not set properly" << endl;
          TH1F* totalSMRatio = (TH1F*) totalSM->Clone();
          TH1F* hErrorPlusRatio = (TH1F*) hErrorPlus->Clone();
          TH1F* hErrorMinusRatio = (TH1F*) hErrorMinus->Clone();
          TH1F* hNominalRatio = (TH1F*) hNominal->Clone();
          if(RatioType==1)
          {
            totalSMRatio->Divide(totalSM);
            hErrorPlusRatio->Divide(totalSM);
            hErrorMinusRatio->Divide(totalSM);
            hNominalRatio->Divide(totalSM);
          }
          else if(RatioType==2)
          {
            totalSMRatio->Add(totalSM,-1);
            totalSMRatio->Divide(totalSM);
            hErrorPlusRatio->Add(totalSM,-1);
            hErrorPlusRatio->Divide(totalSM);
            hErrorMinusRatio->Add(totalSM,-1);
            hErrorMinusRatio->Divide(totalSM);
            hNominalRatio->Add(totalSM,-1);
            hNominalRatio->Divide(totalSM);
          }
          DrawErrorBand(totalSMRatio,hErrorPlusRatio,hErrorMinusRatio,hNominalRatio,RatioErrorGraph_,ErrorBandAroundTotalInput);
        }
      }
    }
  }
  else
  {
    if(hData_) {
      hData_->Draw("E");
    
    }
  }
  
  
  //log-scale plots
  hCanvasStackLogY_->cd();
  hCanvasStackLogY_->SetLogy();
  if(hStack_)
  {
    if(addErrorBand)
      DrawErrorBand(totalSM,hErrorPlus,hErrorMinus,hNominal,StackErrorGraph_,ErrorBandAroundTotalInput);
    
    if(hData_)
    {
      
     hData_->Draw("same E");
      
      
      if(RatioType>0)
      {
        padLogY->Draw();
        padLogY->cd(0);
        ratio->Draw("e");
        if(addRatioErrorBand)
        {
          if(doCutFlow_) cout << "watch out! uncertainties are not set properly" << endl;
          TH1F* totalSMRatio = (TH1F*) totalSM->Clone();
          TH1F* hErrorPlusRatio = (TH1F*) hErrorPlus->Clone();
          TH1F* hErrorMinusRatio = (TH1F*) hErrorMinus->Clone();
          TH1F* hNominalRatio = (TH1F*) hNominal->Clone();
          if(RatioType==1)
          {
            totalSMRatio->Divide(totalSM);
            hErrorPlusRatio->Divide(totalSM);
            hErrorMinusRatio->Divide(totalSM);
            hNominalRatio->Divide(totalSM);
          }
          else if(RatioType==2)
          {
            totalSMRatio->Add(totalSM,-1);
            totalSMRatio->Divide(totalSM);
            hErrorPlusRatio->Add(totalSM,-1);
            hErrorPlusRatio->Divide(totalSM);
            hErrorMinusRatio->Add(totalSM,-1);
            hErrorMinusRatio->Divide(totalSM);
            hNominalRatio->Add(totalSM,-1);
            hNominalRatio->Divide(totalSM);
          }
          DrawErrorBand(totalSMRatio,hErrorPlusRatio,hErrorMinusRatio,hNominalRatio,RatioErrorGraph_,ErrorBandAroundTotalInput);
        }
      }
    }
  }
  else
  {
    if(hData_) {
      hData_->Draw("E");
      
    }
  }
  
  //Now area-normalized plots...
  TH1F* ratioAreaNorm = 0;
  TPad* padAreaNorm = 0;
  TPad* padAreaNormLogY = 0;
  if(hData_)
  {
    ratioAreaNorm = (TH1F*) hData_->Clone();
    double binerror = 0.0;
    if(hStackAreaNorm_)
    {
      if(RatioType==1)
      {
        ratioAreaNorm->Divide(totalSMAreaNorm); //ratio = data/MC
        ratioAreaNorm->GetYaxis()->SetTitle("Data/MC");
        ratioAreaNorm->SetMaximum(1.6);
        ratioAreaNorm->SetMinimum(0.4);
        for(int iBin = 1; iBin < ratioAreaNorm->GetNbinsX()+1; iBin++){
            double a = std::sqrt(hData_->GetBinContent(iBin));
            double b = totalSMAreaNorm->GetBinContent(iBin);
            ratioAreaNorm->SetBinError(iBin, a/b);
            //cout << "bin " << iBin << " error " << binerror << " content " << ratio->GetBinContent(iBin) << endl;
        }
      }
      else if(RatioType==2)
      {
        ratioAreaNorm->Add(totalSMAreaNorm,-1);
        ratioAreaNorm->Divide(totalSMAreaNorm); //ratio = (data-MC)/MC
        ratioAreaNorm->GetYaxis()->SetTitle("(Data-MC)/MC");
        ratioAreaNorm->SetMaximum(0.5);
        ratioAreaNorm->SetMinimum(-0.5);
        for(int iBin = 1; iBin < ratioAreaNorm->GetNbinsX()+1; iBin++){
            double a = std::sqrt(hData_->GetBinContent(iBin));
            double b = totalSMAreaNorm->GetBinContent(iBin);
            ratioAreaNorm->SetBinError(iBin, a/b);
            //cout << "bin " << iBin << " error " << binerror << " content " << ratio->GetBinContent(iBin) << endl;
        }
      }
    }
    ratioAreaNorm->SetMarkerStyle(20);
    ratioAreaNorm->GetXaxis()->SetLabelSize(0.04);
    ratioAreaNorm->GetXaxis()->SetTitleSize(0.045);
    ratioAreaNorm->GetYaxis()->SetLabelSize(0.04);
    //    ratioAreaNorm->GetXaxis()->SetTitle(stack->GetXaxis()->GetTitle());
    ratioAreaNorm->GetYaxis()->SetTitleSize(0.045);
    ratioAreaNorm->GetYaxis()->SetTitleOffset(1.2);
    ratioAreaNorm->SetMarkerSize(1.);
    ratioAreaNorm->GetYaxis()->SetNdivisions(5);
    ratioAreaNorm->SetTitle("");
    
    padAreaNorm = new TPad("padAreaNorm", "padAreaNorm", 0.0, 0.0, 1.0, 1.0);
    padAreaNorm->SetTopMargin(0.7);
    padAreaNorm->SetFillColor(0);
    padAreaNorm->SetFillStyle(0);
    padAreaNorm->SetGridy(1);
    padAreaNorm->SetRightMargin(0.25);
    padAreaNorm->SetLeftMargin(0.12);
    padAreaNormLogY = (TPad*) padAreaNorm->Clone();
  }
  
  hCanvasStackAreaNorm_->cd();
  if(hStackAreaNorm_)
  {
    if(addErrorBand)
      DrawErrorBand(totalSMAreaNorm,hErrorPlus,hErrorMinus,hNominal,StackErrorGraph_,ErrorBandAroundTotalInput); //warning: error band around area-normalized plot may not make sense
    
    if(hData_)
    {
      hData_->Draw("same E");
      
      if(RatioType>0)
      {
        padAreaNorm->Draw();
        padAreaNorm->cd(0);
        ratioAreaNorm->Draw("e");
        if(addRatioErrorBand)
        {
          if(doCutFlow_) cout << "watch out! uncertainties are not set properly" << endl;
          TH1F* totalSMAreaNormRatio = (TH1F*) totalSM->Clone();
          TH1F* hErrorPlusRatio = (TH1F*) hErrorPlus->Clone();
          TH1F* hErrorMinusRatio = (TH1F*) hErrorMinus->Clone();
          TH1F* hNominalRatio = (TH1F*) hNominal->Clone();
          if(RatioType==1)
          {
            totalSMAreaNormRatio->Divide(totalSMAreaNorm);
            hErrorPlusRatio->Divide(totalSMAreaNorm);
            hErrorMinusRatio->Divide(totalSMAreaNorm);
            hNominalRatio->Divide(totalSMAreaNorm);
          }
          else if(RatioType==2)
          {
            totalSMAreaNormRatio->Add(totalSMAreaNorm,-1);
            totalSMAreaNormRatio->Divide(totalSMAreaNorm);
            hErrorPlusRatio->Add(totalSMAreaNorm,-1);
            hErrorPlusRatio->Divide(totalSMAreaNorm);
            hErrorMinusRatio->Add(totalSMAreaNorm,-1);
            hErrorMinusRatio->Divide(totalSMAreaNorm);
            hNominalRatio->Add(totalSMAreaNorm,-1);
            hNominalRatio->Divide(totalSMAreaNorm);
          }
          DrawErrorBand(totalSMAreaNormRatio,hErrorPlusRatio,hErrorMinusRatio,hNominalRatio,RatioErrorGraph_,ErrorBandAroundTotalInput); //warning: error band around area-normalized plot may not make sense
        }
      }
    }
  }
  else
  {
    if(hData_) hData_->Draw("E");
  }
  
  //Now area-normalized log-scale plots...
  hCanvasStackAreaNormLogY_->cd();
  hCanvasStackAreaNormLogY_->SetLogy();
  if(hStackAreaNorm_)
  {
    if(addErrorBand)
      DrawErrorBand(totalSMAreaNorm,hErrorPlus,hErrorMinus,hNominal,StackErrorGraph_,ErrorBandAroundTotalInput); //warning: error band around area-normalized plot may not make sense
    
    if(hData_)
    {
      hData_->Draw("same E");
      
      if(RatioType)
      {
        padAreaNormLogY->Draw();
        padAreaNormLogY->cd(0);
        ratioAreaNorm->Draw("e");
        if(addRatioErrorBand)
        {
          if(doCutFlow_) cout << "watch out! uncertainties are not set properly" << endl;
          TH1F* totalSMAreaNormRatio = (TH1F*) totalSM->Clone();
          TH1F* hErrorPlusRatio = (TH1F*) hErrorPlus->Clone();
          TH1F* hErrorMinusRatio = (TH1F*) hErrorMinus->Clone();
          TH1F* hNominalRatio = (TH1F*) hNominal->Clone();
          if(RatioType==1)
          {
            totalSMAreaNormRatio->Divide(totalSMAreaNorm);
            hErrorPlusRatio->Divide(totalSMAreaNorm);
            hErrorMinusRatio->Divide(totalSMAreaNorm);
            hNominalRatio->Divide(totalSMAreaNorm);
          }
          else if(RatioType==2)
          {
            totalSMAreaNormRatio->Add(totalSMAreaNorm,-1);
            totalSMAreaNormRatio->Divide(totalSMAreaNorm);
            hErrorPlusRatio->Add(totalSMAreaNorm,-1);
            hErrorPlusRatio->Divide(totalSMAreaNorm);
            hErrorMinusRatio->Add(totalSMAreaNorm,-1);
            hErrorMinusRatio->Divide(totalSMAreaNorm);
            hNominalRatio->Add(totalSMAreaNorm,-1);
            hNominalRatio->Divide(totalSMAreaNorm);
          }
          DrawErrorBand(totalSMAreaNormRatio,hErrorPlusRatio,hErrorMinusRatio,hNominalRatio,RatioErrorGraph_,ErrorBandAroundTotalInput); //warning: error band around area-normalized plot may not make sense
        }
      }
    }
  }
  else
  {
    if(hData_) hData_->Draw("E");
  }
  
  if(addErrorBand)
  {
    //make dummy graph for legend; dirty because color/style parameters of the error band now configured in 2 places...
    gStyle->SetHatchesLineWidth(1);
    gStyle->SetHatchesSpacing(0.85);
    TGraphAsymmErrors* ErrorGraph = new TGraphAsymmErrors();
    ErrorGraph->SetFillStyle(3345);//3005 diagonal dashed //3001 ~plain //3013 double diagonal dashed
    ErrorGraph->SetFillColor(kGray);
    ErrorGraph->SetLineColor(kGray);
    ErrorGraph->SetLineWidth(2);
    ErrorGraph->Draw("2SAME");
//    leg_->AddEntry(ErrorGraph,"Uncertainty","F");
  }
  
  if(hData_)
  {
    if(hStack_ && ymax < 1.55*(hData_->GetMaximum())) ymax = 1.55*(hData_->GetMaximum());
    if(hStackAreaNorm_ && ymaxAreaNorm < 1.3*(hData_->GetMaximum())) ymaxAreaNorm = 1.3*(hData_->GetMaximum());
  }
  
  if(maxY_ > 0)
  {
    ymax = maxY_;
    maxLogY_ = logYMultiplicator_ * ymax;
    ymaxAreaNorm = maxY_;
    maxLogYAreaNorm_ = logYMultiplicator_ * ymaxAreaNorm;
  }
  if(hStack_)	hStack_->SetMaximum(ymax);
  if(hStackAreaNorm_) hStackAreaNorm_->SetMaximum(ymaxAreaNorm);
  
  /*	if(pad) pad->cd(0); //otherwise the legend disappears when saving the canvas...
   else hCanvasStack_->cd(0);
   leg_->Draw();
   if(padLogY) padLogY->cd(0); //otherwise the legend disappears when saving the canvas...
   else hCanvasStackLogY_->cd(0);
   leg_->Draw();
   if(padAreaNorm) padAreaNorm->cd(0); //otherwise the legend disappears when saving the canvas...
   else hCanvasStackAreaNorm_->cd(0);
   legAreaNorm_->Draw();
   if(padAreaNormLogY) padAreaNormLogY->cd(0); //otherwise the legend disappears when saving the canvas...
   else hCanvasStackAreaNormLogY_->cd(0);
   legAreaNorm_->Draw();
   */
  hCanvasStack_->cd(0);
  leg_->Draw();
  hCanvasStackLogY_->cd(0);
  leg_->Draw();
  hCanvasStackAreaNorm_->cd(0);
  legAreaNorm_->Draw();
  hCanvasStackAreaNormLogY_->cd(0);
  legAreaNorm_->Draw();
  
  //if(tmpInfile) tmpInfile->Close(); //causes a crash when trying to write hData_ in the Write method...?
}

void MultiSamplePlot::DrawStackedPlot(TCanvas* canvas, TCanvas* canvasLogY, THStack* hstack, vector<TH1F*> histosForOverlay, int scaleNPSignal, const char* xaxistitle, const char* yaxistitle, unsigned int RatioType)
{
  float H = canvas->GetWh();
  float W = canvas->GetWw();
  // references for T, B, L, R
  float T = 0.08*H;
  //float B = 0.12*H;
  float L = 0.12*W;
  float R = 0.25*W;
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  //  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);
  
  canvasLogY->SetFillColor(0);
  canvasLogY->SetBorderMode(0);
  canvasLogY->SetFrameFillStyle(0);
  canvasLogY->SetFrameBorderMode(0);
  canvasLogY->SetLeftMargin( L/W );
  canvasLogY->SetRightMargin( R/W );
  canvasLogY->SetTopMargin( T/H );
  //  canvas->SetBottomMargin( B/H );
  canvasLogY->SetTickx(0);
  canvasLogY->SetTicky(0);
  
  float l = canvas->GetLeftMargin();
  float t = canvas->GetTopMargin();
  float r = canvas->GetRightMargin();
  float b = canvas->GetBottomMargin();
  int iPosX = 11;
  //bool outOfFrame    = false;
  
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;
  
//  stringstream slumi; slumi.precision(2); slumi << lumi_/1000; //luminosity given in picobarns, but will be displayed in femtobarns, because this in not 2010 anymore.
  stringstream slumi; slumi.precision(2); slumi << 36; //luminosity given in picobarns, but will be displayed in femtobarns, because this in not 2010 anymore.
  stringstream ssqrts; ssqrts << sqrts_; //sqrt(s) given in TeV
  TLatex cmstext;
  
  TString lumiText = (slumi.str()+" fb^{-1} ("+ssqrts.str()+" TeV)").c_str();
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;
  if(chan_){
    extraTextSize =  extraTextSize *0.75;
  }
  
  cmstext.SetNDC(true);
  cmstext.SetTextFont(42);
  cmstext.SetTextAlign(31);
  cmstext.SetTextSize(lumiTextSize*t);
  cmstext.SetTextColor(kBlack);
  cmstext.SetTextAngle(0);
  
  float posX_=0;
  if( iPosX%10<=1 )
  {
    posX_ =   l + relPosX*(1-l-r);
  }
  else if( iPosX%10==2 )
  {
    posX_ =  l + 0.5*(1-l-r);
  }
  else if( iPosX%10==3 )
  {
    posX_ =  1-r - relPosX*(1-l-r);
  }
  float posY_ = 1-t - relPosY*(1-t-b);
  
  TLatex text;
  text.SetNDC(true);
  text.SetTextAlign(12);
  text.SetTextFont(42);
  text.SetTextSize(0.04);
  text.SetTextColor(kBlack);
  
  canvas->cd();
  hstack->Draw("HIST");
  hstack->GetXaxis()->SetTitle(xaxistitle);//Make TGaxis to force not too many digits on y-axis
  hstack->GetYaxis()->SetTitle(yaxistitle);
  hstack->GetYaxis()->SetLabelSize(0.04);
  hstack->GetYaxis()->SetTitleSize(0.045);
  hstack->GetYaxis()->SetTitleOffset(1.2);

  TGaxis *myY = (TGaxis*)hstack->GetYaxis();//Make TGaxis to force not too many digits on y-axis
  myY->SetMaxDigits(4);
  if(setBinLabels_ && RatioType == 0){
    for(int ibin = 1; ibin < vlabel_.size()+1; ibin++){
      hstack->GetXaxis()->SetBinLabel(ibin,vlabel_[ibin-1].c_str());
    }
  }
  
  cmstext.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
  
  cmstext.SetTextAlign(align_);
  cmstext.SetTextFont(cmsTextFont);
  cmstext.SetTextSize(cmsTextSize*t);
  
  if(RatioType>0)
  {
    hstack->GetXaxis()->SetLabelSize(0);
    hstack->GetXaxis()->SetTitleSize(0);
  }
  cmstext.DrawLatex(posX_,posY_,"CMS");
  if(chan_)
  {
//      cmstext.DrawLatex(posX_+0.5,posY_,channel_.c_str());
//      extraTextSize =  extraTextSize *0.75;
      cmstext.SetTextSize(0.05);
      cmstext.DrawLatex(posX_+0.5,posY_/*-0.05*/,channel_.c_str());
  }

  if(hData_)
  {
    if(prelim_) {
      cmstext.SetTextFont(extraTextFont);
      cmstext.SetTextAlign(align_);
      cmstext.SetTextSize(extraTextSize*t);
      cmstext.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
    }
  }
  else {
    cmstext.SetTextFont(extraTextFont);
    cmstext.SetTextAlign(align_);
    cmstext.SetTextSize(extraTextSize*t);
    if(chan_) {
      TString text = "Simulation";
      text = text + " - " + channel_ + " channel";
      cmstext.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, text);
    }
    else{
      cmstext.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, "Simulation" );
    }
  }


  
  if(!text_.IsNull()) text.DrawLatex(0.2,0.9,text_);
  
  for(unsigned int i=0;i<histosForOverlay.size();i++)
  {
    histosForOverlay[i]->Scale(scaleNPSignal);
    histosForOverlay[i]->SetFillColor(0);
    histosForOverlay[i]->Draw("HIST SAME");
  }
  
  canvasLogY->cd();
  canvasLogY->SetLogy();
  hstack->Draw("HIST");
  cmstext.SetTextFont(42);
  cmstext.SetTextAlign(31);
  cmstext.SetTextSize(lumiTextSize*t);
  cmstext.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
  
  cmstext.SetTextAlign(align_);
  cmstext.SetTextFont(cmsTextFont);
  cmstext.SetTextSize(cmsTextSize*t);
  cmstext.DrawLatex(posX_,posY_,"CMS");
  
  if(hData_)
  {
    if(prelim_) {
      cmstext.SetTextFont(extraTextFont);
      cmstext.SetTextAlign(align_);
      cmstext.SetTextSize(extraTextSize*t);
      cmstext.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
    }
  }
  else {
    cmstext.SetTextFont(extraTextFont);
    cmstext.SetTextAlign(align_);
    cmstext.SetTextSize(extraTextSize*t);
    cmstext.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, "Simulation");
  }
  
  if(!text_.IsNull()) text.DrawLatex(0.5,0.86,text_);
  
  for(unsigned int i=0;i<histosForOverlay.size();i++)
  {
    histosForOverlay[i]->Draw("HIST SAME");
  }
  
}



void MultiSamplePlot::DrawErrorBand(TH1F* totalSM, TH1F* hErrorPlus, TH1F* hErrorMinus, TH1F* hNominal, TGraphAsymmErrors*  ErrorGraph, bool ErrorBandAroundTotalInput)
{
  //cout <<"In draw error band method..."<< endl;
  
  int nbins = totalSM->GetNbinsX()+1;
  float bins[nbins];
  float bincontents[nbins];
  float dummy[nbins];
  float erroryplus[nbins]; //will be the errors relative to the total histo
  float erroryminus[nbins]; //will be the errors relative to the total histo
  float binwidth = totalSM->GetBinCenter(2) - totalSM->GetBinCenter(1); //assumes a fixed bin width...!
  
  if((totalSM->GetNbinsX() != hErrorPlus->GetNbinsX()) || (totalSM->GetNbinsX() != hErrorMinus->GetNbinsX()))
    cout<<"[MultiSamplePlot::DrawErrorBand] ERROR: error histograms have different binning than total SM histogram!"<<endl;
  
  for(int iBin=0; iBin < totalSM->GetNbinsX()+1; iBin++)
  {
    bins[iBin] = totalSM->GetBinCenter(iBin);
    if(dosystfile_)
    {
      //we still have to take into account if the systematic shifted histogram bin(s) are above the nominal or below, since the input root file contains histograms of just 'raw' systematic shifted distribution(s)
      //			 if(ErrorBandAroundTotalInput){
      bincontents[iBin] = totalSM->GetBinContent(iBin); //meaning an error band around whatever is the total MC SM input (can be 'post-fit' (as in the profile likelihood method), the input is then not the 'nominal' expectation anymore w.r.t which the error band is intrinsically defined)
      
      //First deal with the somewhat rare case where the up/down systematic variation both produce a deviation in the same direction
      if ((  bincontents[iBin] > hErrorPlus->GetBinContent(iBin) )  &&  (  bincontents[iBin] > hErrorMinus->GetBinContent(iBin) ) ){
        erroryplus[iBin] =0;
        if( hErrorMinus->GetBinContent(iBin) <  hErrorPlus->GetBinContent(iBin)){
          erroryminus[iBin] =   bincontents[iBin]  -  hErrorMinus->GetBinContent(iBin);
        }else {
          erroryminus[iBin] =   bincontents[iBin]  -  hErrorPlus->GetBinContent(iBin);
        }
      }
      else if(( bincontents[iBin] < hErrorPlus->GetBinContent(iBin) ) && ( bincontents[iBin] < hErrorMinus->GetBinContent(iBin) ) ){
        erroryminus[iBin] =0;
        
        if( hErrorMinus->GetBinContent(iBin) <  hErrorPlus->GetBinContent(iBin)){
          erroryplus[iBin] =  hErrorPlus->GetBinContent(iBin) - bincontents[iBin] ;
        }else {
          erroryplus[iBin] =  hErrorMinus->GetBinContent(iBin) - bincontents[iBin] ;
        }
        
      }
      //now the normal case where the nominal sits between the up/down variations
      else{
        if ( hErrorPlus->GetBinContent(iBin) > bincontents[iBin]){
          //'standard' expectation up gives up , down gives down
          erroryplus[iBin] = hErrorPlus->GetBinContent(iBin) -  bincontents[iBin];
          erroryminus[iBin] =   bincontents[iBin]  -  hErrorMinus->GetBinContent(iBin);
        }
        else{
          //'flipped'  up gives down , down gives up
          erroryplus[iBin] = hErrorMinus->GetBinContent(iBin) -  bincontents[iBin];
          erroryminus[iBin] =   bincontents[iBin]  -  hErrorPlus->GetBinContent(iBin);
        }
      }
      //			 }
      //			 else {
      //			   bincontents[iBin] = hNominal->GetBinContent(iBin);
      //			   erroryplus[iBin] = hErrorPlus->GetBinContent(iBin) - hNominal->GetBinContent(iBin);
      //			   erroryminus[iBin] = hNominal->GetBinContent(iBin) - hErrorMinus->GetBinContent(iBin);
      //			 }
      dummy[iBin] = binwidth/2;
    }
    else
    {
      //The input systematic shifted histogram has already been produced in the form of the histograms that can directly be used as 'up' and 'down' bands, above and below the nominal, respectively.
      if(ErrorBandAroundTotalInput)
        bincontents[iBin] = totalSM->GetBinContent(iBin); //meaning an error band around whatever is the total MC SM input (can be 'post-fit' (as in the profile likelihood method), the input is then not the 'nominal' expectation anymore w.r.t which the error band is intrinsically defined)
      else
        bincontents[iBin] = hNominal->GetBinContent(iBin);
      erroryplus[iBin] = hErrorPlus->GetBinContent(iBin) - hNominal->GetBinContent(iBin);
      erroryminus[iBin] = hNominal->GetBinContent(iBin) - hErrorMinus->GetBinContent(iBin);
      dummy[iBin] = binwidth/2;
    }
    
    if(iBin == totalSM->GetNbinsX()) cout <<"  bin "<<   iBin << " error plus  " <<  hErrorPlus->GetBinContent(iBin)	  <<  " error minus  " << hErrorMinus->GetBinContent(iBin)  <<" SM total  " <<  bincontents[iBin]   << " band up "  <<  erroryplus[iBin]    << "  band down " <<  erroryminus[iBin]  << endl;
    
  }
  
  ErrorGraph = new TGraphAsymmErrors(nbins,bins,bincontents,dummy,dummy,erroryminus,erroryplus);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.85);
  ErrorGraph->SetFillStyle(3345);//3005 diagonal dashed //3001 ~plain //3013 double diagonal dashed
  ErrorGraph->SetFillColor(kGray);
  ErrorGraph->SetLineColor(kGray);
  ErrorGraph->SetLineWidth(2);
  ErrorGraph->Draw("2SAME");
}

void MultiSamplePlot::Write(TFile* fout, string label, bool savePNG, string pathPNG, string ext)
{
  fout->cd();
  string dirname = "MultiSamplePlot_"+label;
  if(fout->Get(dirname.c_str())==0)
    fout->mkdir(dirname.c_str());
  
  fout->cd(dirname.c_str());
  
  if(hStack_) hStack_->Write();
  if(hStackAreaNorm_) hStackAreaNorm_->Write();
  if(hData_) hData_->Write();
  
  if(hCanvas_)
  {
    if(savePNG)
      hCanvas_->SaveAs( (pathPNG+"/"+label+"_Normalized."+ext).c_str() );
    hCanvas_->Write();
  }
  
  if(hCanvasStack_)
  {
    if(hStack_) hStack_->SetMinimum(0.001);
    hCanvasStack_->Write();
    if(savePNG)
      hCanvasStack_->SaveAs( (pathPNG+"/"+label+"_Stack."+ext).c_str() );
  }
  
  if(hCanvasStackLogY_)
  {
    if(hStack_) hStack_->SetMinimum(minLogY_);
    if(hStack_) hStack_->SetMaximum(maxLogY_);
    hCanvasStackLogY_->Write();
    if(savePNG)
      hCanvasStackLogY_->SaveAs( (pathPNG+"/"+label+"_StackLogY."+ext).c_str() );
  }
  
  if(hCanvasStackAreaNorm_ && hData_)
  {
    if(hStackAreaNorm_) hStackAreaNorm_->SetMinimum(0.001);
    hCanvasStackAreaNorm_->Write();
    if(savePNG)
      hCanvasStackAreaNorm_->SaveAs( (pathPNG+"/"+label+"_StackAreaNorm."+ext).c_str() );
  }
  
  if(hCanvasStackAreaNormLogY_ && hData_)
  {
    if(hStackAreaNorm_) hStackAreaNorm_->SetMinimum(minLogY_);
    if(hStackAreaNorm_) hStackAreaNorm_->SetMaximum(maxLogYAreaNorm_);
    hCanvasStackAreaNormLogY_->Write();
    if(savePNG)
      hCanvasStackAreaNormLogY_->SaveAs( (pathPNG+"/"+label+"_StackAreaNormLogY."+ext).c_str() );
  }
  
  for(unsigned int i=0; i<plots_.size(); i++)
    if(plots_[i].second->Name().find("data") != 0 && plots_[i].second->Name().find("Data") != 0 && plots_[i].second->Name().find("DATA") != 0 )
    {
      plots_[i].first->SetMarkerSize(0.5);

      //making sure that the overflow is transferred to the last 'visible' bin; analogously for underflow...
      Nbins_ = (plots_[i].first)->GetNbinsX();
      (plots_[i].first)->SetBinContent(Nbins_,(plots_[i].first)->GetBinContent(Nbins_)+(plots_[i].first)->GetBinContent(Nbins_+1));
      (plots_[i].first)->SetBinContent(Nbins_+1,0);
      (plots_[i].first)->SetBinContent(1,(plots_[i].first)->GetBinContent(0)+(plots_[i].first)->GetBinContent(1));
      (plots_[i].first)->SetBinContent(0,0);

      plots_[i].first->Write((label+"_"+plots_[i].second->Name()+"_").c_str());
    }
}

void MultiSamplePlot::setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);
  
  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
  
  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);
  
  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);
  
  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);
  
  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);
  
  // For the Global title:
  
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);
  
  // For the axis titles:
  
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
  
  // For the axis labels:
  
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.005, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");
  
  // For the axis:
  
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);
  
  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);
  
  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);
  
  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);
  
  tdrStyle->cd();
  
}

