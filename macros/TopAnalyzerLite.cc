#include "TROOT.h"
#include "TSystem.h"

#include "TCut.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TFile.h"
#include "TChain.h"
#include "TTreePlayer.h"
#include "TFileCollection.h"
#include "THashList.h"
#include "TFileInfo.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TH1F.h"
#include "THStack.h"
#include "TGraph.h"

#include "TMath.h"

#include "TPRegexp.h"

#include "TEntryList.h"

#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace std;

class TopAnalyzerLite
{
public:
  TopAnalyzerLite(const string subDirName = "", const string imageOutDir = "", bool createplots = true, bool printstats = true);
  ~TopAnalyzerLite();

  void addMCSig(const string mcSampleName, const string mcSampleLabel,
                const string fileName, const double xsec, const double nEvents,
                const Color_t color, const bool doStackSignal = true, const TCut cut = "");
  void addMCBkg(const string mcSampleName, const string mcSampleLabel,
                const string fileName, const double xsec, const double nEvents,
                const Color_t color, const TCut cut = "");
  void addDataBkg(const string name, const string label,
                  const string fileName, const double norm,
                  const Color_t color);
  void replaceDataBkgCut(const string name, const string from, const string to);

  void addRealData(const string fileName, const double lumi);

  void addCutStep(const TCut cut, const TString & monitorPlotNamesStr, const double plotScale = 1.0, const string weight = "1", const TString& cutName = "", const TString& postfix = "", const TCut subCut = "1");
  void addMonitorPlot(const string name, const string varexp, const string title,
                      const int nBins, const double xmin, const double xmax,
                      const double ymin = 0, const double ymax = 0, const bool doLogy = true);
  void addMonitorPlot(const string name, const string varexp, const string title,
                      const string xBinsStr,
                      const double ymin = 0, const double ymax = 0, const bool doLogy = true);
  void setEventWeightVar(const string eventWeightVar = "weight");
  void setEventWeight(const string sample, const TString & step, const double & w);
  void setScanVariables(const string scanVariables);

  void applyCutSteps();
  void applySingleCut(const TCut cut, const TString monitirPlotNamesStr, int istep = 0);

  void saveHistograms(TString fileName = "");
  void printCutFlow();

private:
  TObjArray getHistograms();

  struct MCSample
  {
    string name;
    double nEvents;
    double xsec;
    TChain* chain;
    string label;
    Color_t color;
    bool doStack;
    TCut cut;
  };

  struct DataSample
  {
    string name;
    double norm;
    TChain* chain;
    string label;
    Color_t color;
    std::map<std::string, std::string> replaceCuts;
  };

  struct MonitorPlot
  {
    string varexp;
    string title;
    std::vector<double> xBins;
    bool doLogy;
    double ymin, ymax;
  };

  struct CutStep
  {
    TCut cut;
    vector<string> monitorPlotNames;
    double plotScale;
    string weight;
    TString cutName;
    TString postfix;
    TCut subCut;
  };

  struct Stat
  {
    string name, label;
    double nEvents, nEventsErr2;
  };

  double lumi_;
  string subDirName_;
  vector<MCSample> mcSigs_;
  vector<MCSample> mcBkgs_;
  vector<DataSample> dataBkgs_;
  TChain* realDataChain_;

  map<const string, MonitorPlot> monitorPlots_;
  vector<CutStep> cuts_;

  string imageOutDir_;

  void prepareEventList(const TCut &cut, int istep);
  void plot(const string name, TCut cut, MonitorPlot& monitorPlot, const double plotScale = 1.0, const double custStep = 1, const string weight = "1");
  void printStat(const string& name, TCut cut, int cutStep=0);
  void addMC(vector<MCSample>& mcSetup,
             const string name, const string label,
             const string fileName, const double xsec, const double nEvents,
             const Color_t color, bool doStack=true, TCut cut ="");

  TObjArray histograms_;
  ofstream fout_;
  bool writeSummary_;
  string scanVariables_;
  string eventWeightVar_;
  //map<string, vector<double> > wMap_; 
  map<string, map<TString, double> > wMap_; 
  map<TString, vector<Stat> > statsMap_; 
  map<TString, vector<TEntryList *> > entryList_;
  TDirectory* baseRootDir_;

  bool createplots_;
  bool printstats_;
};

TopAnalyzerLite::TopAnalyzerLite(const string subDirName, const string imageOutDir, bool createplots, bool printstats)
{
  subDirName_ = subDirName;
  lumi_ = 0;
  realDataChain_ = 0;
  imageOutDir_ = imageOutDir;

  baseRootDir_ = gROOT->mkdir(subDirName_.c_str());

  if ( imageOutDir != "" )
  {
    gSystem->mkdir(imageOutDir.c_str(), true);
    fout_.open((imageOutDir+"/summary.txt").c_str());
    writeSummary_ = true;
  }
  else writeSummary_ = false;

  scanVariables_ = "";
  eventWeightVar_ = "1";

  createplots_ = createplots;
  printstats_  = printstats;
}

TopAnalyzerLite::~TopAnalyzerLite()
{
  if ( writeSummary_ ) fout_.close();
}

void TopAnalyzerLite::addMC(vector<MCSample>& mcSetup,
                            const string name, const string label,
                            const string fileName, const double xsec, const double nEvents,
                            const Color_t color, bool doStack, TCut cut)
{
  int index = -1;
  for ( unsigned int i=0; i<mcSetup.size(); ++i )
  {
    if ( mcSetup[i].name == name )
    {
      index = i;
      break;
    }
  }

  if ( index == -1 )
  {
    MCSample mc = {name, 0, xsec, 0, label, color, doStack, cut};
    baseRootDir_->cd();
    //mc.chain = new TChain((subDirName_+"/tree").c_str(), (subDirName_+"/tree").c_str());
    mc.chain = new TChain((subDirName_).c_str(), (subDirName_).c_str());
    //mc.chain = new TChain("tree", "tree");
    mcSetup.push_back(mc);
    index = mcSetup.size()-1;
  }
  MCSample& mc = mcSetup[index];
  mc.chain->Add(fileName.c_str());

  if ( nEvents > 0 )
  {
    mc.nEvents += nEvents;
  }
  else
  {
    TFile* f = TFile::Open(fileName.c_str());
    if ( !f || !f->IsOpen() ) cout << "Cannot open file\n";
    else 
    {
      //TH1* hEventSummary = (TH1*)f->Get((subDirName_+"/EventSummary").c_str());
      TH1* hEventSummary = (TH1*)f->Get("EventSummary");
      if ( !hEventSummary ) cout << "Cannot find EventSummary histogram" << endl;
      else mc.nEvents += hEventSummary->GetBinContent(1);
      f->Close();
    }
  }
}

void TopAnalyzerLite::addMCSig(const string name, const string label,
                               const string fileName, const double xsec, const double nEvents,
                               const Color_t color, const bool doStackSignal, const TCut cut)
{
  addMC(mcSigs_, name, label, fileName, xsec, nEvents, color, doStackSignal, cut);
}

void TopAnalyzerLite::addMCBkg(const string name, const string label,
                               const string fileName, const double xsec, const double nEvents,
                               const Color_t color, const TCut cut)
{
  addMC(mcBkgs_, name, label, fileName, xsec, nEvents, color, true, cut); //stack true for background
}

void TopAnalyzerLite::addDataBkg(const string name, const string label,
                                 const string fileName, const double norm,
                                 const Color_t color)
{
  int index = -1;
  for ( unsigned int i=0; i<dataBkgs_.size(); ++i )
  {
    if ( dataBkgs_[i].name == name )
    {
      index = i;
      break;
    }
  }

  if ( index == -1 )
  {
    std::map<std::string, std::string> replaceCuts;
    DataSample data = {name, norm, 0, label, color, replaceCuts};
    baseRootDir_->cd();
    //data.chain = new TChain((subDirName_+"/tree").c_str(), (subDirName_+"/tree").c_str());
    data.chain = new TChain((subDirName_).c_str(), (subDirName_).c_str());
    //data.chain = new TChain("tree", "tree");
    dataBkgs_.push_back(data);
    index = dataBkgs_.size()-1;
  }

  DataSample& data = dataBkgs_[index];
  data.chain->Add(fileName.c_str());
}

void TopAnalyzerLite::prepareEventList(const TCut &cut, int istep)
{
  char entrylistname[100];
  char puttoentrylistname[100];
	
  if (realDataChain_) {
    sprintf(entrylistname, "realdata%d", istep);
    sprintf(puttoentrylistname, ">> %s", entrylistname);
    realDataChain_->Draw(puttoentrylistname, cut, "entrylist");
    TEntryList *list = (TEntryList *)gDirectory->Get(entrylistname);
    if (list==0) cout << "error!" << endl;
    realDataChain_->SetEntryList(list);
    entryList_["realdata"].push_back(list);
    if( entryList_["realdata"].back() == 0 ) cout << "Error!" << endl;
  }

  for ( unsigned int i=0; i<mcSigs_.size(); ++i )
  {
    MCSample& mcSample = mcSigs_[i];
    TCut finalCut = cut + mcSample.cut;
    sprintf(entrylistname, "mcsig%d_%d", i, istep);
    sprintf(puttoentrylistname, ">> %s", entrylistname);
    mcSample.chain->Draw(puttoentrylistname, finalCut, "entrylist");
    TEntryList *list = (TEntryList *)gDirectory->Get(entrylistname);
    mcSample.chain->SetEntryList(list);
    entryList_[Form("mcsig%d",i)].push_back(list);
  }
  for ( unsigned int i=0; i<mcBkgs_.size(); ++i )
  {
    MCSample& mcSample = mcBkgs_[i];
    TCut finalCut = cut + mcSample.cut;
    sprintf(entrylistname, "mcbkg%d_%d", i, istep);
    sprintf(puttoentrylistname, ">> %s", entrylistname);
    mcSample.chain->Draw(puttoentrylistname, finalCut, "entrylist");
    TEntryList *list = (TEntryList *)gDirectory->Get(entrylistname);
    mcSample.chain->SetEntryList(list);
    entryList_[Form("mcbkg%d",i)].push_back(list);
  }
  for ( unsigned int i=0; i<dataBkgs_.size(); ++i )
  {
    DataSample& sample = dataBkgs_[i];
    TString cutStr;
    cutStr = cut;
    map<string, string>::const_iterator cit;
    for(cit=sample.replaceCuts.begin(); cit != sample.replaceCuts.end() ; cit++){
      cutStr.ReplaceAll((*cit).first, (*cit).second);
    }
    sprintf(entrylistname, "databkg%d_%d", i, istep);
    sprintf(puttoentrylistname, ">> %s", entrylistname);
    sample.chain->Draw(puttoentrylistname, cutStr, "entrylist");
    TEntryList *list = (TEntryList *)gDirectory->Get(entrylistname);
    sample.chain->SetEntryList(list);
    entryList_[Form("databkg%d",i)].push_back(list);
  }
}

void TopAnalyzerLite::replaceDataBkgCut(const string name, const string from, const string to)
{
  for ( unsigned int i=0; i<dataBkgs_.size(); ++i )
  {
    DataSample& dataBkg = dataBkgs_[i];
    if ( dataBkg.name != name ) continue;
    dataBkg.replaceCuts[from] = to;
  }
}

void TopAnalyzerLite::addRealData(const string fileName, const double lumi)
{
  lumi_ += lumi;
  if ( fileName == "" ) return;
if ( !realDataChain_ )
  {
    //const string chainName = subDirName_+"/tree";
    const string chainName = subDirName_;
    //const string chainName = "tree";
    baseRootDir_->cd();
    realDataChain_ = new TChain(chainName.c_str(), chainName.c_str());
  }

  realDataChain_->Add(fileName.c_str());
}

void TopAnalyzerLite::addCutStep(const TCut cut, const TString & monitorPlotNamesStr, const double plotScale, const string weight, const TString& cutName, const TString& postfix, TCut subCut)
{
  TObjArray* monitorPlotNames = monitorPlotNamesStr.Tokenize(",");
  const int nPlots = monitorPlotNames->GetSize();

  vector<string> plotNames;
  for ( int i=0; i<nPlots; ++i )
  {
    TObject* obj = monitorPlotNames->At(i);
    if ( !obj ) continue;

    const string plotName = obj->GetName();
  
    plotNames.push_back(plotName);
  }
 
  int nstep = (int)cuts_.size()+1; 
  TString dirName = cutName;
  if( cutName == "" ) dirName = Form("Step_%d", nstep); 
  CutStep cutStep = {cut, plotNames, plotScale, weight, dirName, postfix, subCut};
  cuts_.push_back(cutStep);
}

void TopAnalyzerLite::addMonitorPlot(const string name, const string varexp, const string title,
                                     const int nBins, const double xmin, const double xmax,
                                     const double ymin, const double ymax, const bool doLogy)
{
  std::vector<double> xBins;
  const double dX = (xmax-xmin)/nBins;
  for ( int i=0; i<=nBins; ++i )
  {
    xBins.push_back(xmin+dX*i);
  }

  MonitorPlot monitorPlot = {varexp, title, xBins, doLogy, ymin, ymax};
  monitorPlots_[name] = monitorPlot;
}

void TopAnalyzerLite::addMonitorPlot(const string name, const string varexp, const string title,
                                     const string xBinsStr,
                                     const double ymin, const double ymax, const bool doLogy)
{
  stringstream ss(xBinsStr);
  std::vector<double> xBins;
  double x;
  while(ss >> x ) xBins.push_back(x);

  MonitorPlot monitorPlot = {varexp, title, xBins, doLogy, ymin, ymax};
  monitorPlots_[name] = monitorPlot;
}

void TopAnalyzerLite::applyCutSteps()
{
  cout << "--------------------------------------\n";
  cout << " Cross sections and sample statistics \n";
  if ( writeSummary_ )
  {
    fout_ << "--------------------------------------\n";
    fout_ << " Cross sections and sample statistics \n";
  }
  for ( unsigned int i=0; i<mcSigs_.size(); ++i )
  {
    MCSample& mcSample = mcSigs_[i];
    double effL = mcSample.nEvents/mcSample.xsec;
    cout << " * " << mcSample.name << "\t" << mcSample.xsec << " pb (" << effL << " /pb, " << mcSample.nEvents << ")\n";
    if ( writeSummary_ ) fout_ << " * " << mcSample.name << "\t" << mcSample.xsec << " /pb (" << mcSample.nEvents << ")\n";
  }
  for ( unsigned int i=0; i<mcBkgs_.size(); ++i )
  {
    MCSample& mcSample = mcBkgs_[i];
    double effL = mcSample.nEvents/mcSample.xsec;
    cout << " * " << mcSample.name << "\t" << mcSample.xsec << " pb (" << effL << " /pb, " << mcSample.nEvents << ")\n";
    if ( writeSummary_ ) fout_ << " * " << mcSample.name << "\t" << mcSample.xsec << " /pb (" << mcSample.nEvents << ")\n";
  }
  cout << "--------------------------------------\n";
  if ( writeSummary_ ) fout_ << "--------------------------------------\n";

  TCut cut = "";
  for ( unsigned int i=0; i<cuts_.size(); ++i )
  {
    cut = cut && cuts_[i].cut;
    const vector<string>& monitorPlotNames = cuts_[i].monitorPlotNames;
    const double plotScale = cuts_[i].plotScale;
    const string w = cuts_[i].weight;
    TString cname = cuts_[i].cutName;
    TString postfix = cuts_[i].postfix;
    prepareEventList(cut, i);
    if( printstats_ ) printStat(Form("%s", cname.Data() ), cut, i);
    for ( unsigned int j = 0; j < monitorPlotNames.size(); ++ j)
    {
      const string& plotName = monitorPlotNames[j];

      if ( monitorPlots_.find(plotName) == monitorPlots_.end() ) continue;
      MonitorPlot& monitorPlot = monitorPlots_[plotName];
      plot(Form("%s_%s%s", cname.Data(), plotName.c_str(), postfix.Data() ), cuts_[i].subCut, monitorPlot, lumi_*plotScale, i, w);
    }
  }

  //cout << "Final" << endl;
  //if ( writeSummary_ ) fout_ << "Final" << endl;
  //TCut finalCut = "";
  //for ( unsigned int i=0; i<cuts_.size(); ++i )
  //{
  //  finalCut = finalCut && cuts_[i].cut;
  //}
  //if ( realDataChain_ )
  //{
    //realDataChain_->Scan(scanVariables_.c_str());
    //cout << "Number of entries after final selection = " << entryList_["realdata"].back()->GetN() << endl;
  //}

  if ( writeSummary_ && realDataChain_ && printstats_ )
  {
    printCutFlow();
    const string tmpFileName = imageOutDir_+"/tmp.txt";

    ((TTreePlayer*)(realDataChain_->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(realDataChain_->GetPlayer()))->SetScanFileName(tmpFileName.c_str());
    //realDataChain_->Scan(scanVariables_.c_str(), finalCut);
    realDataChain_->Scan(scanVariables_.c_str());
    ((TTreePlayer*)(realDataChain_->GetPlayer()))->SetScanRedirect(false);

    ifstream tmpFile(tmpFileName.c_str());
    copy(istreambuf_iterator<char>(tmpFile), istreambuf_iterator<char>(), ostreambuf_iterator<char>(fout_));
    //fout_ << "Number of entries after final selection = " << realDataChain_->GetEntries(finalCut) << endl;
    fout_ << "Number of entries after final selection = " << entryList_["realdata"].back()->GetN() << endl;
    gSystem->Exec(("rm -f "+tmpFileName).c_str());
  }
}

void TopAnalyzerLite::plot(const string name, const TCut cut, MonitorPlot& monitorPlot, const double plotScale, const double cutStep, const string weight)
{
  const string& varexp = monitorPlot.varexp;
  const string& title = monitorPlot.title;
  const int nBins = monitorPlot.xBins.size()-1;
  const double* xBins = &(monitorPlot.xBins[0]);
  double ymin = monitorPlot.ymin;
  double ymax = monitorPlot.ymax*plotScale;

  baseRootDir_->cd();
  TLegend* legend = new TLegend(0.73,0.57,0.88,0.88);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->SetLineColor(0);

  //TString dataHistName = Form("hData_%s_%s", subDirName_.c_str(), name.c_str());
  TString dataHistName = Form("hData_%s", name.c_str());
  TH1F* hData = new TH1F(dataHistName, title.c_str(), nBins, xBins);
  histograms_.Add(hData);

  if ( realDataChain_ ) realDataChain_->Project(dataHistName, varexp.c_str(), cut);
  hData->AddBinContent(nBins, hData->GetBinContent(nBins+1));
  hData->Sumw2();
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(1);
  hData->SetTitle(title.c_str());
  hData->SetStats(0);

  legend->AddEntry(hData, "Data", "p");

  TString dataSubHistName = Form("hDataSub_%s", name.c_str());
  TH1F* hDataSub = (TH1F*)hData->Clone(dataSubHistName);
  hDataSub->SetTitle(hData->GetTitle()+TString(" background subtracted"));
  histograms_.Add(hDataSub);

  THStack* hStack = new THStack("hStack", title.c_str());
  typedef vector<pair<string, TH1F*> > LabeledPlots;
  LabeledPlots stackedPlots;
  LabeledPlots sigPlots; // Keep list of signal plots if doStackSignal == false

  TString cutStr;
  cutStr = cut;

  for ( unsigned int i=0; i<mcSigs_.size(); ++i )
  {
    MCSample& mcSample = mcSigs_[i];
    TString mcSigHistName = Form("hMCSig_%s_%s", mcSample.name.c_str(), name.c_str());
    TH1F* hMCSig = new TH1F(mcSigHistName, title.c_str(), nBins, xBins);

    TCut mcWeightStr = Form("(%s)*(%s)*(%s)", eventWeightVar_.c_str(),weight.c_str(),cutStr.Data());
    mcSample.chain->Project(mcSigHistName, varexp.c_str(),mcWeightStr);
    hMCSig->AddBinContent(nBins, hMCSig->GetBinContent(nBins+1));
    hMCSig->Scale(lumi_*mcSample.xsec/mcSample.nEvents);
    
    if ( mcSample.doStack )
    {
      hMCSig->SetFillColor(mcSample.color);
      hMCSig->SetFillStyle(1001);

      stackedPlots.push_back(make_pair(mcSample.label, hMCSig));
      hStack->Add(hMCSig);

      histograms_.Add(hMCSig);
    }
    else
    {
      LabeledPlots::const_iterator matchedPlot = sigPlots.end();
      for ( LabeledPlots::const_iterator plotIter = sigPlots.begin();
            plotIter != sigPlots.end(); ++plotIter )
      {
        if ( plotIter->first == mcSample.label )
        {
          matchedPlot = plotIter;
          break;
        }
      }
      if ( matchedPlot == sigPlots.end() )
      {
        sigPlots.push_back(make_pair(mcSample.label, hMCSig));

//        hMCSig->SetLineWidth(2);
        hMCSig->SetLineStyle(mcSample.color);
        hMCSig->SetLineColor(mcSample.color);
        histograms_.Add(hMCSig);
      }
      else
      {
        TH1F* h = matchedPlot->second;
        if ( h ) h->Add(hMCSig);
        delete hMCSig;
      }

    }
  }

  for ( unsigned int i=0; i<mcBkgs_.size(); ++i )
  {
    MCSample& mcSample = mcBkgs_[i];
    //TString mcHistName = Form("hMC_%s_%s_%s", subDirName_.c_str(), mcSample.name.c_str(), name.c_str());
    TString mcHistName = Form("hMC_%s_%s", mcSample.name.c_str(), name.c_str());
    TH1F* hMC = new TH1F(mcHistName, title.c_str(), nBins, xBins);

    TCut mcWeightStr = Form("(%s)*(%s)*(%s)", eventWeightVar_.c_str(), weight.c_str(), cutStr.Data());
    mcSample.chain->Project(mcHistName, varexp.c_str(),mcWeightStr);
    hMC->AddBinContent(nBins, hMC->GetBinContent(nBins+1));
    hMC->Scale(lumi_*mcSample.xsec/mcSample.nEvents);

    //scale MC
    //map<string, vector<double> >::iterator it;
    map<string, map<TString, double> >::iterator it;
    it = wMap_.find(mcSample.name);
    if( it != wMap_.end() ) {
      //hMC->Scale(it->second[cutStep]);
      hMC->Scale(it->second[cuts_[cutStep].cutName]);
    }

    hMC->SetFillColor(mcSample.color);
    hMC->SetFillStyle(1001);

    // Subtract background from the hDataSub histogram
    hDataSub->Add(hMC, -1);

    // Add to the HStack if there's no duplicated label
    // If duplicated label exists, call TH1::Add
    // First, find if plot with same label already in the THStack
    LabeledPlots::const_iterator matchedPlot = stackedPlots.end();
    for ( LabeledPlots::const_iterator plotIter = stackedPlots.begin();
        plotIter != stackedPlots.end(); ++plotIter )
    {
      if ( plotIter->first == mcSample.label )
      {
        matchedPlot = plotIter;
        break;
      }
    }
    // If the label was not in the stack, insert it
    if ( matchedPlot == stackedPlots.end() )
    {
      stackedPlots.push_back(make_pair(mcSample.label, hMC));
      hStack->Add(hMC);
      histograms_.Add(hMC);
    }
    // If there's plot with same label, sum entries
    else
    {
      TH1F* h = matchedPlot->second;
      if ( h ) h->Add(hMC);
      // In this case, temporary histogram is not needed anymore.
      delete hMC;
    }
  }

  for ( unsigned int i=0; i<dataBkgs_.size(); ++i )
  {
    DataSample& sample = dataBkgs_[i];
    TString histName = Form("hDataBkg_%s_%s", sample.name.c_str(), name.c_str());
    TH1F* hBkg = new TH1F(histName, title.c_str(), nBins, xBins);

    sample.chain->Project(histName, varexp.c_str(), cut);
    hBkg->AddBinContent(nBins, hBkg->GetBinContent(nBins+1));
    hBkg->Scale(sample.norm);

    //scale MC
    //map<string, vector<double> >::iterator it;
    map<string, map<TString, double> >::iterator it;
    it = wMap_.find(sample.name);
    if( it != wMap_.end() ) {
      hBkg->Scale(it->second[cuts_[cutStep].cutName]);
    }

    hBkg->SetFillColor(sample.color);
    hBkg->SetFillStyle(1001);

    // Subtract background from the hDataSub histogram
    hDataSub->Add(hBkg, -1);

    LabeledPlots::const_iterator matchedPlot = stackedPlots.end();
    for ( LabeledPlots::const_iterator plotIter = stackedPlots.begin();
          plotIter != stackedPlots.end(); ++plotIter )
    {
      if ( plotIter->first == sample.label )
      {
        matchedPlot = plotIter;
        break;
      }
    }

    if ( matchedPlot == stackedPlots.end() )
    {
      stackedPlots.push_back(make_pair(sample.label, hBkg));
      hStack->Add(hBkg);
      histograms_.Add(hBkg);
    }
    else
    {
      TH1F* h = matchedPlot->second;
      if ( h ) h->Add(hBkg);
      // In this case, temporary histogram is not needed anymore.
      delete hBkg;
    }
  }

  // Do automatic bin labels
  if ( xBins[0] == 0 and xBins[nBins] == nBins and nBins < 20 )
  {
    const int xmin = xBins[0];
    TList* hList = hStack->GetHists();
    for ( int bin=1; bin<nBins; ++bin )
    {
      hData->GetXaxis()->SetBinLabel(bin, Form("%d", int(xmin+bin-1)));
    }
    hData->GetXaxis()->SetBinLabel(nBins, Form("#geq%d", int(xmin+nBins-1))); 

    for ( int i=0; i<hList->GetSize(); ++i )
    {
      TH1* h = (TH1*)hList->At(i);

      for ( int bin=1; bin<nBins; ++bin )
      {
        h->GetXaxis()->SetBinLabel(bin, Form("%d", int(xmin+bin-1)));
      }
      h->GetXaxis()->SetBinLabel(nBins, Form("#geq%d", int(xmin+nBins-1)));
    }

  }

  // Build legend, legend should be added in reversed order of THStack
  for ( int i=stackedPlots.size()-1; i>=0; --i )
  {
    const char* label = stackedPlots[i].first.c_str();
    TH1F* h = stackedPlots[i].second;
    legend->AddEntry(h, label, "f");
  }

  for ( unsigned int i=0; i<sigPlots.size(); ++i )
  {
    const char* label = sigPlots[i].first.c_str();
    TH1F* h = sigPlots[i].second;
    legend->AddEntry(h, label, "l");
  }

  TCanvas* c = new TCanvas(Form("c_%s", name.c_str()), name.c_str(), 1);
  if ( ymax == 0 )
  {
    const int dataMaxBin = hData->GetMaximumBin();
    const double dataYmax = hData->GetBinContent(dataMaxBin) + hData->GetBinError(dataMaxBin);
    const double mcYmax = hStack->GetMaximum();

    ymax = TMath::Max(dataYmax, mcYmax);
  }

  if ( monitorPlot.doLogy )
  {
    if ( ymin <= 0 ) ymin = 1e-2;
    c->SetLogy();
  }

  hStack->SetMinimum(ymin);
  hStack->SetMaximum(ymax);

  hStack->Draw();
  for ( unsigned int i=0; i<sigPlots.size(); ++i ){ 
    sigPlots[i].second->Draw("same");
  }
  hData->Draw("same");

  legend->Draw();

  if ( createplots_ && imageOutDir_ != "" )
  {
    c->Print((imageOutDir_+"/"+c->GetName()+".png").c_str());
    c->Print((imageOutDir_+"/"+c->GetName()+".eps").c_str());
    //c->Print((imageOutDir_+"/"+c->GetName()+".pdf").c_str());
  }
}

void TopAnalyzerLite::printStat(const string& name, TCut cut, int cutStep)
{
  cout << "-------------------------\n";
  cout << "   " << name << endl;

  if ( writeSummary_ )
  {
    fout_ << "-------------------------\n";
    fout_ << "   " << name << endl;
  }

  const double nData = realDataChain_ ? entryList_["realdata"].at(cutStep)->GetN() : 0;

  double nTotal = 0, nSignal = 0;
  double nTotalErr2 = 0;

  TString cutStr;
  cutStr = cut;

  vector<Stat> stats;
  for ( unsigned int i=0; i<mcSigs_.size(); ++i )
  {
    MCSample& mcSample = mcSigs_[i];

    const double norm = lumi_*mcSample.xsec/mcSample.nEvents;
    const double nEvents = entryList_[Form("mcsig%d",i)].at(cutStep)->GetN()*norm;
    const double nEventsErr2 = nEvents*norm;

    // Merge statistics with same labels
    vector<Stat>::iterator matchedStatObj = stats.end();
    for ( vector<Stat>::iterator statObj = stats.begin();
        statObj != stats.end(); ++statObj )
    {
      if ( statObj->label == mcSample.label )
      {
        matchedStatObj = statObj;
        break;
      }
    }
    if ( matchedStatObj == stats.end() )
    {
      Stat stat = {mcSample.name, mcSample.label, nEvents, nEventsErr2};
      stats.push_back(stat);
    }
  }

  for ( unsigned int i=0; i<mcBkgs_.size(); ++i )
  {
    MCSample& mcSample = mcBkgs_[i];

    //scale MC
    double scale = 1;
    //map<string, vector<double> >::iterator it;
    map<string, map<TString, double> >::iterator it;
    it = wMap_.find(mcSample.name);
    if( it != wMap_.end() ) {
      scale = it->second[cuts_[cutStep].cutName];
    }
    const double norm = lumi_*mcSample.xsec/mcSample.nEvents;
    double rawN = 0;
    if(  entryList_[Form("mcbkg%d",i)].at(cutStep) != NULL ){
      rawN = entryList_[Form("mcbkg%d",i)].at(cutStep)->GetN();
    } 
    const double nEvents = rawN*norm*scale;
    const double nEventsErr2 = nEvents*norm;

    // Merge statistics with same labels
    vector<Stat>::iterator matchedStatObj = stats.end();
    for ( vector<Stat>::iterator statObj = stats.begin();
        statObj != stats.end(); ++statObj )
    {
      if ( statObj->label == mcSample.label )
      {
        matchedStatObj = statObj;
        break;
      }
    }
    if ( matchedStatObj == stats.end() )
    {
      Stat stat = {mcSample.name, mcSample.label, nEvents, nEventsErr2};
      stats.push_back(stat);
    }
    else
    {
      matchedStatObj->nEvents += nEvents;
      matchedStatObj->nEventsErr2 += nEventsErr2;
    }
  }

  for ( unsigned int i=0; i<dataBkgs_.size(); ++i )
  {
    DataSample& sample = dataBkgs_[i];

    map<string, string>::const_iterator cit;
    for(cit=sample.replaceCuts.begin(); cit != sample.replaceCuts.end() ; cit++){
      cutStr.ReplaceAll((*cit).first, (*cit).second);
    }

    //scale MC
    double scale = 1;
    //map<string, vector<double> >::iterator it;
    map<string, map<TString, double> >::iterator it;
    it = wMap_.find(sample.name);
    if( it != wMap_.end() ) {
      scale = it->second[cuts_[cutStep].cutName];
    }

    const double norm = sample.norm;
    const double nEvents = entryList_[Form("databkg%d",i)].at(cutStep)->GetN()*norm*scale;
    const double nEventsErr2 = nEvents*norm;

    vector<Stat>::iterator matchedStatObj = stats.end();
    for ( vector<Stat>::iterator statObj = stats.begin();
          statObj != stats.end(); ++statObj )
    {
      if ( statObj->label == sample.label )
      {
        matchedStatObj = statObj;
        break;
      }
    }
    if ( matchedStatObj == stats.end() )
    {
      Stat stat = {sample.name, sample.label, nEvents, nEventsErr2};
      stats.push_back(stat);
    }
    else
    {
      matchedStatObj->nEvents += nEvents;
      matchedStatObj->nEventsErr2 += nEventsErr2;
    }
  }

  // Get the field width for printing
  int maxFWidth = 0;
  for ( int i=stats.size()-1; i>=0; --i )
  {
    const int fWidth = stats[i].label.size();
    if ( fWidth > maxFWidth ) maxFWidth = fWidth;
  }
  TString form = TString("%-") + Form("%d", maxFWidth) + "s";

  // Print out statistics
  for ( int i=stats.size()-1; i>=0; --i )
  {
    Stat& stat = stats[i];
    const string label = Form(form.Data(), stat.label.c_str());

    cout << label << " = " << stat.nEvents << " +- " << sqrt(stat.nEventsErr2) << endl;
    if ( writeSummary_ ) fout_ << label << " = " << stat.nEvents << " +- " << sqrt(stat.nEventsErr2) << endl;

    bool isSignal = false;
    for ( unsigned int j=0; j<mcSigs_.size(); ++j )
    {
      if ( stat.name == mcSigs_[j].name )
      {
        isSignal = true;
        break;
      }
    }
    if ( isSignal ) nSignal = stat.nEvents;

    nTotal += stat.nEvents;
    nTotalErr2 += stat.nEventsErr2;
  }

  //const double nBkg = nTotal - nSignal;
  const double signif = nSignal/sqrt(nTotal);
  const double nTotalErr = sqrt(nTotalErr2);

  cout << Form(form.Data(), "Total") << " = " << nTotal << " +- " << nTotalErr << endl;
  cout << "-----------------------------------------------" << endl;
  cout << Form(form.Data(), "Data") << " = " << nData << endl;
  cout << Form(form.Data(), "S/sqrt(S+B)") << " = " << signif << endl;
  cout << "-----------------------------------------------" << endl;

  if ( writeSummary_ )
  {
    fout_ << Form(form.Data(), "Total") << " = " << nTotal << " +- " << nTotalErr << endl;
    fout_ << "-----------------------------------------------" << endl;
    fout_ << Form(form.Data(), "Data") << " = " << nData << endl;
    fout_ << Form(form.Data(), "S/sqrt(S+B)") << " = " << signif << endl;
    fout_ << "-----------------------------------------------" << endl;
  }

  statsMap_[name] = stats;
}

void TopAnalyzerLite::applySingleCut(const TCut cut, const TString monitorPlotNamesStr, int istep)
{
  static int singleCutUniqueId = 0;

  TObjArray* monitorPlotNames = monitorPlotNamesStr.Tokenize(",");
  const int nPlots = monitorPlotNames->GetSize();

  cout << "----------------------------\n";
  cout << "Result of single cut" << endl;
  cout << "Cut = " << cut << endl;
  prepareEventList(cut, istep);
  if( printstats_ ) printStat(Form("SingleCut_%d", singleCutUniqueId), cut, istep);
  for ( int i=0; i<nPlots; ++i )
  {
    TObject* obj = monitorPlotNames->At(i);
    if ( !obj ) continue;

    const string plotName = obj->GetName();
    if ( monitorPlots_.find(plotName) == monitorPlots_.end() ) continue;
    MonitorPlot& monitorPlot = monitorPlots_[plotName];
    plot(Form("SingleCut_%d_%s", singleCutUniqueId, plotName.c_str()), cut, monitorPlot);
  }

  monitorPlotNames->Delete();

  ++singleCutUniqueId;
}

TObjArray TopAnalyzerLite::getHistograms()
{
  return histograms_;
}

void TopAnalyzerLite::saveHistograms(TString fileName)
{
  if ( fileName == "" )
  {
    if ( imageOutDir_ != "" ) fileName = imageOutDir_+"/Hist.root";
    else fileName = "Hist.root";
  }

  TFile* f = TFile::Open(fileName, "recreate");

  TParameter<double> lumi("lumi", lumi_);
  lumi.Write();

  TH1F* hScale = new TH1F("hScale", "Scale factors for each samples", mcBkgs_.size()+1, 0, mcBkgs_.size()+1);
  for ( unsigned int i=0; i<mcSigs_.size(); ++i )
  {
    hScale->Fill(mcSigs_[i].name.c_str(), lumi_*mcSigs_[i].xsec/mcSigs_[i].nEvents);
  }
  for ( unsigned int i=0; i<mcBkgs_.size(); ++i )
  {
    hScale->Fill(mcBkgs_[i].name.c_str(), lumi_*mcBkgs_[i].xsec/mcBkgs_[i].nEvents);
  }

  TCut cut;
  TObjArray histograms = getHistograms();

  vector<TString> dirNames;
  for ( unsigned int i=0; i<cuts_.size(); ++i )
  {
    TString dirName = cuts_[i].cutName;
    TDirectory* dir = f->GetDirectory(dirName);
    if ( !dir ) {
      dir = f->mkdir(dirName);
      dir->cd();
      cut += cuts_[i].cut; 
      TNamed cutStr("cut", cut);
      cutStr.Write();
      TNamed subCutStr("subCut", cuts_[i].subCut);
      subCutStr.Write();
      dirNames.push_back(dirName);
    }
  }

  for ( unsigned int i=0; i<dirNames.size(); ++i )
  {
    TString dirName = dirNames[i];
    TDirectory* dir = f->GetDirectory(dirName);
    dir->cd();

    //need to fix : avoid loop  
    //we can make histograms with vector so that only relevant histograms can be taken
    for ( int j=0; j<histograms.GetSize(); ++j ){
      TH1F* h = (TH1F*)histograms.At(j);
      if ( !h ) continue;
      TString hName = h->GetName();
      if( hName.Contains( dirName.Data() )) {
        h->Write();
      }
    }
  }

  f->Write();
  f->Close();
}

void TopAnalyzerLite::setScanVariables(const string scanVariables)
{
  scanVariables_ = scanVariables;
}

void TopAnalyzerLite::setEventWeightVar(const string eventWeightVar)
{
  eventWeightVar_ = eventWeightVar;
}

void TopAnalyzerLite::setEventWeight(const string sample, const TString& step, const double & w)
{
  wMap_[sample.c_str()][step] = w;
}

void TopAnalyzerLite::printCutFlow(){
 
  fout_ << "Cut Flow Table as a latex format" << endl;
  fout_ << "--------------------------------------\n";

  map<TString, vector<Stat> >::iterator it;  
  it = statsMap_.begin();
  int nSample = it->second.size();

  int maxFWidth = 0;
  for ( int i=0; i < nSample; i++ )
  {
    const int fWidth = (*it).second[i].label.size();
    if ( fWidth > maxFWidth ) maxFWidth = fWidth;
  }
  TString form = TString("%-") + Form("%d", maxFWidth) + "s";
 
  map<TString, double> nTotal;
  map<TString, double> nTotalErr2;
  for( int i=0 ; i < nSample ; i++ )
  { 
    it = statsMap_.begin();
    const string label = Form(form.Data(), (*it).second[i].label.c_str());

    fout_ << label << " " ;
    for( int k = 0; k != (int) statsMap_.size() ; k++){
      Stat& stat = (*it).second[i];
      fout_ << " \t&" << setprecision(4) << stat.nEvents ;
      if( k >= (int) statsMap_.size()-1 ) fout_ << " $\\pm$ " << setprecision(4) << sqrt(stat.nEventsErr2) ;
      nTotal[Form("Step_%d", k+1) ] += stat.nEvents;
      nTotalErr2[Form("Step_%d", k+1)] += stat.nEventsErr2;
      it++;
    }
    fout_ << "\\\\ \n" ;
  } 

  map<TString, double >::iterator itTotal;
  map<TString, double >::iterator itTotalErr2;
  itTotal= nTotal.begin();
  itTotalErr2= nTotalErr2.begin();

  fout_ << Form(form.Data(), "Total MC") << " " ;
  for( int k = 0; k != (int) statsMap_.size() ; k++){
    fout_ << " \t&" << setprecision(4) << (*itTotal).second ;
    if( k >= (int) statsMap_.size()-1 ) fout_ << " $\\pm$ " << setprecision(4) << sqrt( (*itTotalErr2).second ) ;
    itTotal++;
    itTotalErr2++;
  }
  fout_ << "\\\\\\hline \n" ;
  fout_ << Form(form.Data(), "Data") << " " ;  
  for( int k = 0; k != (int) statsMap_.size() ; k++){
    const double nData = realDataChain_ ? entryList_["realdata"].at(k)->GetN() : 0;
    fout_ << " \t&" << nData  ;
  }
  fout_ << "\\\\\\hline\\hline \n" ;
}
