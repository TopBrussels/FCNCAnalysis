#define TreeAnalyzer_cxx
//#include "TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include "TStyle.h"
#include "TPaveText.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <TFile.h>
#include <TLeaf.h>
#include <utility>
#include "Style.C"

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"


struct HighestPt{
  bool operator()(TLorentzVector j1, TLorentzVector j2) const
  {
    return j1.Pt() > j2.Pt();
  }
  
  
};

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotCutfl; // TOFIX


map<string,TFile*> tFileMap;

map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;


vector < Dataset* > datasets;

std::vector < int>  decayChannels = {0,1,2,3,-9}; // uuu uue eeu eee all
//std::vector < int>  decayChannels = {-9};

///////////////////////////////////// HOME MADE FUNCTIONS AND VARS  /////////////////////////////////////////


int verbose = 2;


string MakeTimeStamp();
void InitMSPlots(string prefix, vector<int> decayChannels);
void Init1DPlots();
void Init2DPlots();
void InitTree(TTree* tree, bool isData);
// data from global tree
void ClearMetaData();
void GetMetaData(TTree* tree, bool isData,int Entries, bool isAMC);
// put everything to default values
void ClearObjects();
void ClearVars();
void ClearMVAVars();
void ClearTLVs();
void ClearMatchingVars();
void ClearMatchingVarsTLV();
void ClearMatchingSampleVars();
void FillGeneralPlots(int d, string prefix, vector<int>decayChannels, bool isData);
string ConvertIntToString(int nb, bool pad);

int nEntries;
double Xsect;
string pathOutputdate = "";



// Declaration of leaf types
Float_t         MVA_Zboson_pt;
Float_t         MVA_Zboson_eta;
Float_t         MVA_dRWlepb;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_charge_asym;
Float_t         MVA_dRZb;
Float_t         MVA_dRZWlep;
Float_t         MVA_dRZSMtop;
Float_t         MVA_dPhiZb;
Float_t         MVA_dPhiZSMtop;
Float_t         MVA_SMtop_eta;
Float_t         MVA_dPhiZWlep;
Float_t         MVA_region;
Float_t         MVA_weight;
Float_t         MVA_channel;
Float_t         BDT;
Double_t        MVA_EqLumi;

// List of branches
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_Zboson_eta;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_dRZb;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_dRZSMtop;   //!
TBranch        *b_MVA_dPhiZb;   //!
TBranch        *b_MVA_dPhiZSMtop;   //!
TBranch        *b_MVA_SMtop_eta;   //!
TBranch        *b_MVA_dPhiZWlep;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_weight;   //!
TBranch        *b_i_channel;   //!
TBranch        *b_BDT;
TBranch        *b_MVA_EqLumi;


double Luminosity = 32000.;
double EquilumiSF = 1.;
///////////////////////////////////// MAIN CODE /////////////////////////////////////////
int main(int argc, char* argv[]){
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  //setTDRStyle(); // TO FIX
  //setMyStyle(); // TO FIX stat box title
  
  string xmlFileName = "";
  xmlFileName = "config/Run2TriLepton_samples_analy.xml" ;
  const char* xmlFile = xmlFileName.c_str();
  cout << " - Using config file " << xmlFile << endl;
  
  string pathOutput = "OutputPlots/";
  mkdir(pathOutput.c_str(),0777);
  pathOutputdate = pathOutput + dateString + "/"  ;
  mkdir(pathOutputdate.c_str(),0777);
  
  
  
  //  load datasets
  datasets.clear();
  
  
  
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasets, xmlFile);
  cout << "datasets loaded" <<endl;
  bool makePlots = false;
  string placeNtup = "MVAoutput/170214";
  int channel = -999;
  bool datafound = false;
  bool testing = false;
  bool doZct = false;
  bool doZut = false;
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./Analyzer [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   Zct / Zut: make plots for this one" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. MVAoutput/placeNtup " << endl;
      std::cout << "   Test: loop over 100 events" << endl;
      return 0;
    }
    if(string(argv[i]).find("Zct")!=std::string::npos) {
      doZct = true;
    }
    if(string(argv[i]).find("Zut")!=std::string::npos) {
      doZut = true;
    }
    if(string(argv[i]).find("Test")!=std::string::npos) {
      testing = true;
    }
    if(string(argv[i]).find("Ntup")!=std::string::npos) {
      placeNtup = argv[i+1];
      i++;
      
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    
  }
  
  if(makePlots){

    if(doZut) InitMSPlots("singletop_Zut", decayChannels);
    if(doZct) InitMSPlots("singletop_Zct", decayChannels);
    Init1DPlots();
    Init2DPlots();
  }

  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
    {
      Luminosity = datasets[d]->EquivalentLumi();
      datafound = true;
    }
    
    
    
  
    
  }
  
  ////////////////////////////////////
  ///  Open files and get Ntuples  ///
  ////////////////////////////////////
  
  string dataSetName, slumi;
  double timePerDataSet[datasets.size()];
  bool isData = false;
  
  
  /// Loop over datasets
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    clock_t startDataSet = clock();
    
   // ClearMetaData();
    //cout << "meta data cleared" << endl;
    dataSetName = datasets[d]->Name();
    Xsect = datasets[d]->Xsection();
    if (verbose > 1)
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
    }
    
    isData = false;
    if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data")!= std::string::npos || dataSetName.find("DATA")!= std::string::npos )
    {
      isData = true;
    }
    
    string ntupleFileName = "MVAoutput/";
    if(doZut) ntupleFileName += "Zut/"+placeNtup+"/Control_Trees_MVA_regionEq0.root";
    if(doZct) ntupleFileName += "Zct/"+placeNtup+"/Control_Trees_MVA_regionEq0.root";
    tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
    
    string tTreeName = "Control_"+dataSetName;
    
    /// Get data
    if(doZct && dataSetName.find("NP_overlay_ST_FCNC_zut_80X")!=std::string::npos) continue;
    if(doZut && dataSetName.find("NP_overlay_ST_FCNC_zct_80X")!=std::string::npos) continue;
    
    cout << "   treename " << tTreeName << " from " << ntupleFileName <<  endl;
    tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
    nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
    cout << "                nEntries: " << nEntries << endl;
    
    
     // Set branch addresses and branch pointers
     InitTree(tTree[dataSetName.c_str()], isData);
     
     
     int endEvent =nEntries;
     if(testing){
     if(endEvent > 100) endEvent = 100;
     }
     int istartevt = 0;
   
     for (int ievt = 0; ievt < endEvent; ievt++)
     {
   //  ClearObjects(); // put everything to default values
     
       
     if (ievt%1000 == 0)
     std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
     
     
     /// Load event
     tTree[(dataSetName).c_str()]->GetEntry(ievt);
     
       EquilumiSF = MVA_EqLumi;
       
      if(dataSetName.find("NP_overlay")!=std::string::npos) MVA_weight = MVA_weight /10.;
     
     if (makePlots)
     {
     //cout << "ievt " << ievt << endl;
        if(MVA_region == 0 && doZut) FillGeneralPlots(d, "singletop_Zut", decayChannels, isData);
        if(MVA_region == 0 && doZct) FillGeneralPlots(d, "singletop_Zct", decayChannels, isData);
     }
     
     } // events
    
    
    
  } // data
  
  
  
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string rootFileName ="NtuplePlots.root";
  string place =pathOutputdate+"/MSPlotMVA/";
  string placeTH1F = pathOutputdate+"/TH1F/";
  string placeTH2F = pathOutputdate+"/TH2F/";
  vector <string> vlabel_chan = {"uuu", "uue", "eeu", "eee"};
  mkdir(place.c_str(),0777);
  mkdir(placeTH1F.c_str(),0777);
  mkdir(placeTH2F.c_str(),0777);
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
  cout << "   Output file is " << pathOutputdate+rootFileName << endl;
  
  ///Write histograms
  fout->cd();
  
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    if(!datafound) temp->setDataLumi(Luminosity);
    
    if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
    if(name.find("eee")!=std::string::npos) temp->setChannel(true, "eee");
    if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "eeu");
    if(name.find("uue")!=std::string::npos) temp->setChannel(true, "uue");
    if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "uuu");
    temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    cout << "writing to " << pathOutputdate+"MSPlot" << endl;
    cout << "plot " << name << endl;
    cout << "temp " << temp << endl;
    temp->Write(fout, name, true, (pathOutputdate+"/MSPlotMVA").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  

  TDirectory* th1dir = fout->mkdir("1D_histograms");
  th1dir->cd();
  gStyle->SetOptStat(1110);
  for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    TH1F *temp = it->second;
    int N = temp->GetNbinsX();
    temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
    temp->SetBinContent(N+1,0);
    temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
    tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
  }
  
  // 2D
  TDirectory* th2dir = fout->mkdir("2D_histograms");
  th2dir->cd();
  gStyle->SetPalette(55);
  for(std::map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    TH2F *temp = it->second;
    temp->Write();
    TCanvas* tempCanvas = TCanvasCreator(temp, it->first, "colz");
    tempCanvas->SaveAs( (placeTH2F+it->first+".png").c_str() );
  }
  
  fout->Close();
  
  delete fout;
  
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    float secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  
  return 0;
  
}  // end main

///////////////////////////////////// BOOKKEEPING /////////////////////////////////////////
string ConvertIntToString(int Number, int pad){
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}
string MakeTimeStamp(){
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}

///////////////////////////////////// CLEARING /////////////////////////////////////////



///////////////////////////////////// INIT PLOTS /////////////////////////////////////////
void Init1DPlots(){
  TH1::SetDefaultSumw2();
  
  //histo1D["eventID"] = new TH1F("eventID","eventID",5802564204,207137,5802564203);
  
}

void InitMSPlots(string prefix, vector <int> decayChannels){
  clock_t start_sub = clock();
  string decaystring = "";
  // control plots
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
  //  cout << "init " << (prefix+"_BDT_"+decaystring).c_str() << endl;
   MSPlot[(prefix+"_BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), 61, -0.3,1, "BDT");
   
  }
  
  
  
  
  
}

void Init2DPlots(){
  TH2::SetDefaultSumw2();
  //histo2D["CosTheta"]= new TH2F("CosTheta", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);
}

///////////////////////////////////// INIT TREES /////////////////////////////////////////
void InitTree(TTree* tree, bool isData){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("MVA_Zboson_pt", &MVA_Zboson_pt, &b_MVA_Zboson_pt);
  tree->SetBranchAddress("MVA_Zboson_eta", &MVA_Zboson_eta, &b_MVA_Zboson_eta);
  tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
  tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
  tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
  tree->SetBranchAddress("MVA_dRZb", &MVA_dRZb, &b_MVA_dRZb);
  tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
  tree->SetBranchAddress("MVA_dRZSMtop", &MVA_dRZSMtop, &b_MVA_dRZSMtop);
  tree->SetBranchAddress("MVA_dPhiZb", &MVA_dPhiZb, &b_MVA_dPhiZb);
  tree->SetBranchAddress("MVA_dPhiZSMtop", &MVA_dPhiZSMtop, &b_MVA_dPhiZSMtop);
  tree->SetBranchAddress("MVA_SMtop_eta", &MVA_SMtop_eta, &b_MVA_SMtop_eta);
  tree->SetBranchAddress("MVA_dPhiZWlep", &MVA_dPhiZWlep, &b_MVA_dPhiZWlep);
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  tree->SetBranchAddress("MVA_weight", &MVA_weight, &b_weight);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_i_channel);
  tree->SetBranchAddress("BDT", &BDT, &b_BDT);
  tree->SetBranchAddress("MVA_EqLumi",&MVA_EqLumi,&b_MVA_EqLumi);
  
}




void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData){
  //cout << "fill plots" << endl;
  string decaystring = "";
  Double_t eventW = 1.;
  if(!isData) eventW = MVA_weight;
  else eventW = MVA_weight;
  
  
 
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) &&(decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
   // cout << "bdt " << BDT << " in " << (prefix+"_BDT_"+decaystring).c_str()<< endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()]->Fill(BDT , datasets[d], true, eventW);
 
  }

  
  
  //cout << "end plot filling" << endl;
}

