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
void InitMSPlots(string prefix, vector<int> decayChannels, bool istoppair);
void Init1DPlots();
void Init2DPlots();
void InitTree(TTree* tree, bool isData, bool istoppair);
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
void FillGeneralPlots(int d, string prefix, vector<int>decayChannels, bool isData, bool istoppair);
string ConvertIntToString(int nb, bool pad);
vector<double> BDTCUT(string region, string coupling);
int nEntries;
double Xsect;
string pathOutputdate = "";



// Declaration of leaf types
Float_t         MVA_Zboson_pt;
Float_t         MVA_dRWlepb;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_charge_asym;
Float_t         MVA_dRZWlep;
Float_t         MVA_bdiscCSVv2_jet_0;
Float_t         MVA_mlb;
Float_t         MVA_region;
Float_t         MVA_weight;
Float_t         MVA_channel;
Float_t         BDT;
Double_t        MVA_EqLumi;
Float_t         MVA_weight_puSF_up;
Float_t         MVA_weight_puSF_down;
Float_t         MVA_weight_electronSF_up;
Float_t         MVA_weight_electronSF_down;
Float_t         MVA_weight_muonSF_up;
Float_t         MVA_weight_muonSF_down;
Float_t         MVA_weight_btagSF_cferr1_up;
Float_t         MVA_weight_btagSF_cferr1_down;
Float_t         MVA_weight_btagSF_cferr2_up;
Float_t         MVA_weight_btagSF_cferr2_down;
Float_t         MVA_weight_btagSF_hf_up;
Float_t         MVA_weight_btagSF_hf_down;
Float_t         MVA_weight_btagSF_hfstats1_up;
Float_t         MVA_weight_btagSF_hfstats1_down;
Float_t         MVA_weight_btagSF_hfstats2_up;
Float_t         MVA_weight_btagSF_hfstats2_down;
Float_t         MVA_weight_btagSF_lf_up;
Float_t         MVA_weight_btagSF_lf_down;
Float_t         MVA_weight_btagSF_lfstats1_up;
Float_t         MVA_weight_btagSF_lfstats1_down;
Float_t         MVA_weight_btagSF_lfstats2_up;
Float_t         MVA_weight_btagSF_lfstats2_down;
Float_t         MVA_cdiscCvsL_jet_1;
Float_t         MVA_cdiscCvsL_jet_0;
Float_t         MVA_dRZc;


// List of branches
TBranch        *b_MVA_cdiscCvsL_jet_1;
TBranch        *b_MVA_cdiscCvsL_jet_0;
TBranch        *b_MVA_dRZc;
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_0;   //!
TBranch        *b_MVA_mlb;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_MVA_weight;   //!
TBranch        *b_MVA_channel;   //!
TBranch        *b_BDT;   //!
TBranch        *b_MVA_EqLumi;   //!
TBranch        *b_MVA_weight_puSF_up;   //!
TBranch        *b_MVA_weight_puSF_down;   //!
TBranch        *b_MVA_weight_electronSF_up;   //!
TBranch        *b_MVA_weight_electronSF_down;   //!
TBranch        *b_MVA_weight_muonSF_up;   //!
TBranch        *b_MVA_weight_muonSF_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_down;   //!
TBranch        *b_MVA_weight_btagSF_hf_up;   //!
TBranch        *b_MVA_weight_btagSF_hf_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_down;   //!
TBranch        *b_MVA_weight_btagSF_lf_up;   //!
TBranch        *b_MVA_weight_btagSF_lf_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_down;   //!

string placeOutputReading ="";
string template_name = "";
int nbin = 10;
double weight = 1.;
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
  xmlFileName = "config/Run2TriLepton_samples_analyBDT.xml" ;
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
  bool doPseudoData = false;
  bool doSystematics = false;
  string placeNtup = "MVAoutput/170214";
  int channel = -999;
  bool datafound = false;
  bool testing = false;
  bool doZct = false;
  bool doZut = false;
  bool toppair = false;
  bool addData = false;
  bool DetermineCut = false;
  
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./AnalyzerBDT [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   Zct / Zut: make plots for this one" << endl;
      std::cout << "   toppair: make plots for this one" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. MVAoutput/placeNtup " << endl;
      std::cout << "   Test: loop over 100 events" << endl;
      std::cout << "   Data: add CRcut" << endl;
      std::cout << "   PSdata: generate pseudo data" << endl;
      std::cout << "   Systematics: loop over systematics" << endl;
      return 0;
    }
    if(string(argv[i]).find("Systematics")!=std::string::npos) {
      doSystematics= true;
    }
    
    if(string(argv[i]).find("PSdata")!=std::string::npos) {
      doPseudoData = true;
    }
    if(string(argv[i]).find("Data")!=std::string::npos) {
      addData = true;
    }
    if(string(argv[i]).find("toppair")!=std::string::npos) {
      toppair = true;
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
      cout << placeNtup << endl;
      i++;
      
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    if(string(argv[i]).find("DetermineCut")!=string::npos) {
      DetermineCut= true;
    }
    
  }
  
  
  
  
  vector <string> thesystlist;
  thesystlist.clear();
  thesystlist.push_back(""); // nominal
  if(doSystematics){
    cout << "pushing back systematics" << endl;
    thesystlist.push_back("puSF_down");
    thesystlist.push_back("electronSF_down");
    thesystlist.push_back("muonSF_down");
    thesystlist.push_back("btagSF_cferr1_down");
    thesystlist.push_back("btagSF_cferr2_down");
    thesystlist.push_back("btagSF_hf_down");
    thesystlist.push_back("btagSF_hfstats1_down");
    thesystlist.push_back("btagSF_hfstats2_down");
    thesystlist.push_back("btagSF_lf_down");
    thesystlist.push_back("btagSF_lfstats1_down");
    thesystlist.push_back("btagSF_lfstats2_down");
    
    thesystlist.push_back("puSF_up");
    thesystlist.push_back("electronSF_up");
    thesystlist.push_back("muonSF_up");
    thesystlist.push_back("btagSF_cferr1_up");
    thesystlist.push_back("btagSF_cferr2_up");
    thesystlist.push_back("btagSF_hf_up");
    thesystlist.push_back("btagSF_hfstats1_up");
    thesystlist.push_back("btagSF_hfstats2_up");
    thesystlist.push_back("btagSF_lf_up");
    thesystlist.push_back("btagSF_lfstats1_up");
    thesystlist.push_back("btagSF_lfstats2_up");
    
    thesystlist.push_back("JER_up");
    thesystlist.push_back("JER_down");
    thesystlist.push_back("JES_up");
    thesystlist.push_back("JES_down");
  }
  if(makePlots){
    string tempstring = "singletop_Zut";
    TString systematics = "";
    for(int isys = 0; isys < thesystlist.size() ; isys++){
      systematics = thesystlist[isys];
      if( doZut && !toppair) tempstring = "singletop_Zut";
      if( doZct && !toppair) tempstring = "singletop_Zct";;
      if(doZut && toppair)  tempstring = "toppair_Zut";
      if( doZct && toppair) tempstring = "toppair_Zut";
      if(isys != 0 ) tempstring += "_"+ systematics;
      InitMSPlots(tempstring, decayChannels, toppair);
    }
    Init1DPlots();
    Init2DPlots();
  }
  
  string coupling = "Zut";
  if(!doZut && doZct) coupling = "Zct";
  placeOutputReading = "MVAoutput/outputtemplates";
  mkdir(placeOutputReading.c_str(), 0777);
  placeOutputReading += "/" + coupling + "_";
  if(toppair) placeOutputReading += "toppair";
  else placeOutputReading += "singletop";
  mkdir(placeOutputReading.c_str(), 0777);
  template_name = coupling + "_" ;
  if(toppair) template_name += "toppair";
  else template_name += "singletop";
  
  
  if(DetermineCut){
    double cut = -999.;
    string region = "";
    if(toppair) region = "toppair";
    else region = "singletop";
    string coupling = "";
    if(doZut) coupling = "Zut";
    else coupling = "Zct";
    vector<double> v_cut = BDTCUT(region, coupling);
    /* if(MVA_channel == 0) cut = v_cut[0];
     if(MVA_channel == 1) cut = v_cut[1];
     if(MVA_channel == 2) cut = v_cut[2];
     if(MVA_channel == 3) cut = v_cut[3];*/
    return 0;
  }
  
  
  
  
  string combinetemplate_filename = placeOutputReading+"/Reader_" + template_name + ".root";
  TFile* combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "RECREATE" );
  
  
  TH1::SetDefaultSumw2();
  TH1F *hist_BDT(0), *hist_BDTG(0);
  TH1F *hist_uuu = 0, *hist_uue = 0, *hist_eeu = 0, *hist_eee = 0;
  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    string dataSetName = datasets[d]->Name();
    if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  )
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
  string systematic = "";
  string tTreeName = "";
  string postfix = "";
  string ntupleFileName = "MVAoutput/";
  string output_histo_name = "";
  string sample_name = "";
  /// Loop over datasets
  for(int isys = 0; isys < thesystlist.size() ; isys++){
    systematic = thesystlist[isys];
    //cout << "systematic " << systematic << endl;
    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      clock_t startDataSet = clock();
      
      // ClearMetaData();
      //cout << "meta data cleared" << endl;
      dataSetName = datasets[d]->Name();
      Xsect = datasets[d]->Xsection();
      // if(dataSetName.find("NP_overlay_TT_FCNC_T2ZJ_aTleptonic_ZToll_kappa_zut")==std::string::npos) continue;
      
      if (verbose > 1)
      {
        cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
      }
      
      isData = false;
      if ( dataSetName.find("Data") != std::string::npos || dataSetName.find("data")!= std::string::npos || dataSetName.find("DATA")!= std::string::npos || dataSetName.find("fake")!=std::string::npos)
      {
        isData = true;
      }
      
      ntupleFileName = "MVAoutput/";
      ntupleFileName = ntupleFileName + "/" + placeNtup ;
      tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
      
      
      postfix = "";
      if(systematic.find("JES_down")!=std::string::npos) postfix = "_JESdown";
      else if(systematic.find("JES_up")!=std::string::npos) postfix = "_JESup";
      else if(systematic.find("JER_down")!=std::string::npos) postfix = "_JERdown";
      else if(systematic.find("JER_up")!=std::string::npos) postfix = "_JERup";
      else postfix = "";
      //cout << "postfix " << postfix << endl;
      tTreeName = "Control_"+dataSetName + postfix;
      
      /// Get data
      cout << "   treename " << tTreeName << " from " << ntupleFileName <<  endl;
      tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries: " << nEntries << endl;
      
      
      // Set branch addresses and branch pointers
      InitTree(tTree[dataSetName.c_str()], isData, toppair);
      
      std::cout << "--- Select "<<datasets[d]->Name()<<" sample" << std::endl;
      // prepare outpout histograms
      hist_uuu     = new TH1F( (template_name+"_uuu").c_str(),           (template_name+"_uuu").c_str(),           nbin, -1, 1 );
      hist_uue     = new TH1F( (template_name+"_uue").c_str(),           (template_name+"_uue").c_str(),           nbin, -1, 1 );
      hist_eeu     = new TH1F( (template_name+"_eeu").c_str(),           (template_name+"_eeu").c_str(),           nbin, -1, 1 );
      hist_eee     = new TH1F( (template_name+"_eee").c_str(),           (template_name+"_eee").c_str(),           nbin, -1, 1 );
      
      
      if(isData && isys != 0) continue;
      if(!doSystematics && isys != 0) continue;
      int endEvent =nEntries;
      if(testing){
        if(endEvent > 100) endEvent = 100;
      }
      int istartevt = 0;
      
      for (int ievt = 0; ievt < endEvent; ievt++)
      {
        //  ClearObjects(); // put everything to default values
        weight = 1.;
        
        if (ievt%1000 == 0)
          std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
        
        
        /// Load event
        tTree[(dataSetName).c_str()]->GetEntry(ievt);
        
        EquilumiSF = MVA_EqLumi;
        weight = MVA_weight;
        if(!isData){
          if(systematic.find("puSF_up")) weight = MVA_weight_puSF_up;
          if(systematic.find("puSF_down")) weight = MVA_weight_puSF_down;
          if(systematic.find("electronSF_up")) weight = MVA_weight_electronSF_up;
          if(systematic.find("electronSF_down")) weight = MVA_weight_electronSF_down;
          if(systematic.find("muonSF_up")) weight = MVA_weight_muonSF_up;
          if(systematic.find("muonSF_down")) weight = MVA_weight_muonSF_down;
          if(systematic.find("btagSF_cferr1_up")) weight = MVA_weight_btagSF_cferr1_up;
          if(systematic.find("btagSF_cferr1_down")) weight = MVA_weight_btagSF_cferr1_down;
          if(systematic.find("btagSF_cferr2_up")) weight = MVA_weight_btagSF_cferr2_up;
          if(systematic.find("btagSF_cferr2_down")) weight = MVA_weight_btagSF_cferr2_down;
          if(systematic.find("btagSF_lf_up")) weight = MVA_weight_btagSF_lf_up;
          if(systematic.find("btagSF_lf_down")) weight = MVA_weight_btagSF_lf_down;
          if(systematic.find("btagSF_hf_up")) weight = MVA_weight_btagSF_hf_up;
          if(systematic.find("btagSF_hf_down")) weight = MVA_weight_btagSF_hf_down;
          if(systematic.find("btagSF_hfstats1_up")) weight = MVA_weight_btagSF_hfstats1_up;
          if(systematic.find("btagSF_hfstats1_down")) weight = MVA_weight_btagSF_hfstats1_down;
          if(systematic.find("btagSF_hfstats2_up")) weight = MVA_weight_btagSF_hfstats2_up;
          if(systematic.find("btagSF_hfstats2_down")) weight = MVA_weight_btagSF_hfstats2_down;
          if(systematic.find("btagSF_lfstats2_up")) weight = MVA_weight_btagSF_lfstats2_up;
          if(systematic.find("btagSF_lfstats2_down")) weight = MVA_weight_btagSF_lfstats2_down;
          
        }
        
        
        if(MVA_channel== 0) 		{hist_uuu->Fill( BDT, weight);}
        else if(MVA_channel== 1) {hist_uue->Fill( BDT, weight);}
        else if(MVA_channel== 2) {hist_eeu->Fill( BDT, weight);}
        else if(MVA_channel == 3) {hist_eee->Fill( BDT, weight);}
        else if(MVA_channel== 9 || weight == 0) {cout<<__LINE__<<" : problem  weight" << weight <<endl;}
        
        
        //if(dataSetName.find("NP_overlay")!=std::string::npos) MVA_weight = MVA_weight /10.;
        
        /*  if(addData){
         string region = "";
         if(toppair) region = "toppair";
         else region = "singletop";
         string coupling = "";
         if(doZut) coupling = "Zut";
         else coupling = "Zct";
         v_cut = BDTCUT(region, coupling);
         if(MVA_channel == 0) cut = v_cut[0];
         if(MVA_channel == 1) cut = v_cut[1];
         if(MVA_channel == 2) cut = v_cut[2];
         if(MVA_channel == 3) cut = v_cut[3];
         
         }
         */
        
        // if(addData && BDT < cut) continue;
        if (makePlots)
        {
          //cout << "ievt " << ievt << endl;
          string tempstring = "singletop_Zut";
          if( doZut && !toppair) tempstring = "singletop_Zut";
          if( doZct && !toppair) tempstring = "singletop_Zct";;
          if(doZut && toppair) tempstring = "toppair_Zut";
          if( doZct && toppair)tempstring = "toppair_Zut";
          if(isys != 0) tempstring += "_"+ systematic;
          FillGeneralPlots(d, tempstring, decayChannels, isData, toppair);
        }
        
      } // events
      
      // --- Write histograms
      combinetemplate_file->cd();
      
      //NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)] FIX ME
      output_histo_name = "";
      sample_name = datasets[d]->Name();
      if (sample_name.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
      {
        if(isys!=0) output_histo_name = template_name+"_uuu_FakeMu_"  + systematic ;
        else output_histo_name = template_name+"_uuu_FakeMu"  ;
        hist_uuu->SetTitle(output_histo_name.c_str());
        hist_uuu->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_uue_FakeEl_"  + systematic ;
        else output_histo_name = template_name+"_uue_FakeEl"  ;
        hist_uue->SetTitle(output_histo_name.c_str());
        hist_uue->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_eeu_FakeMu_"  + systematic ;
        else output_histo_name = template_name+"_eeu_FakeMu"  ;
        hist_eeu->SetTitle(output_histo_name.c_str());
        hist_eeu->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_eee_FakeEl_"  + systematic ;
        else output_histo_name = template_name+"_eee_FakeEl"  ;
        hist_eee->SetTitle(output_histo_name.c_str());
        hist_eee->Write(output_histo_name.c_str());
      }
      else //If fakes are not considered, or if sample is not fake --> write directly !
      {
        if(isys!=0) output_histo_name = template_name+"_uuu_"  + sample_name + "_" + systematic ;
        else output_histo_name = template_name+"_uuu_"  + sample_name ;
        hist_uuu->SetTitle(output_histo_name.c_str());
        hist_uuu->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_uue_"  + sample_name + "_" + systematic ;
        else output_histo_name = template_name+"_uue_"  + sample_name ;
        hist_uue->SetTitle(output_histo_name.c_str());
        hist_uue->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_eeu_"  + sample_name + "_" + systematic ;
        else output_histo_name = template_name+"_eeu_"  + sample_name ;
        hist_eeu->SetTitle(output_histo_name.c_str());
        hist_eeu->Write(output_histo_name.c_str());
        if(isys!=0) output_histo_name = template_name+"_eee_"  + sample_name + "_" + systematic ;
        else output_histo_name = template_name+"_eee_"  + sample_name ;
        hist_eee->SetTitle(output_histo_name.c_str());
        hist_eee->Write(output_histo_name.c_str());
      }
      
      cout<<"Done with "<<datasets[d]->Name()<<" sample"<<endl;
      
    } // data
    
    //don't delete if processing MC fake templates (unless all the loops have reached their ends)
    if(isys == (thesystlist.size() - 1))
    {
      //cout<<"deleting dynamic histograms"<<endl;
      delete hist_uuu; delete hist_uue; delete hist_eeu; delete hist_eee;
    }
    
    if(isys != 0) cout<<"Done with "<< systematic <<" systematic"<<endl;
    else cout<<"Done with nominal sample"<<endl;
    
    
    
  } // systlist
  
  
  combinetemplate_file->Close();
  std::cout << "--- Created root file: \""<<combinetemplate_file->GetName()<<"\" containing the output histograms" << std::endl;
  std::cout << "==> Reader() is done!" << std::endl << std::endl;
  
  delete combinetemplate_file;
  
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  if(makePlots){
    string rootFileName ="NtupleMVAPlots.root";
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
      if(name.find("eee")!=std::string::npos) temp->setChannel(true, "3e");
      if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "2e1#mu");
      if(name.find("uue")!=std::string::npos) temp->setChannel(true, "1e2#mu");
      if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "3#mu");
      temp->Draw(name, 1, false, false, false, 50);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
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
  }
  
  if(doPseudoData && !addData){
    cout << "generating pseudo data" << endl;
    TRandom3 therand(0); //Randomization
    
    string pseudodata_input_name = placeOutputReading+"/Reader_" + template_name + ".root";
    TFile* pseudodata_file = TFile::Open( pseudodata_input_name.c_str(), "UPDATE" );
    
    cout << pseudodata_input_name << endl;
    
    
    
    TH1F *h_sum = 0, *h_tmp = 0;
    string histo_name = "";
    string template_fake_name = "";
    
    
    std::vector<string> channel_list;
    channel_list.push_back("eee");
    channel_list.push_back("uue");
    channel_list.push_back("eeu");
    channel_list.push_back("uuu");
    for(int ichan=0; ichan<channel_list.size(); ichan++)
    {
      h_sum = 0;
      histo_name = "";
      if(channel_list[ichan] == "uuu" || channel_list[ichan] == "eeu") {template_fake_name = "FakeMu";}
      else {template_fake_name = "FakeEl";}
      
      
      
      
      for(int isample = 0; isample < datasets.size(); isample++)
      {
        string dataSetName = datasets[isample]->Name();
        // cout << dataSetName << endl;
        if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {continue; } // no signal in data
        if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // no signal in data
        if(datasets[isample]->Name().find("fake")==std::string::npos) {
          // cout << "  -- sample " << datasets[isample]->Name() << endl;
          h_tmp = 0;
          histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
          //histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
          //  cout << "  --- histo " << histo_name << endl;
          if(!pseudodata_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
          h_tmp = (TH1F*) pseudodata_file->Get(histo_name.c_str());
          if(h_sum == 0) {h_sum = (TH1F*) h_tmp->Clone();}
          else {h_sum->Add(h_tmp);}
        }
        else{
          
          histo_name = template_name + "_" + channel_list[ichan] + "_" + template_fake_name;
          cout << "  --- histo " << histo_name << endl;
          if(!pseudodata_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<histo_name<<" : not found"<<endl;}
          else
          {
            h_tmp = (TH1F*) pseudodata_file->Get(histo_name.c_str());
            if(h_sum == 0) {h_sum = (TH1F*) h_tmp->Clone();}
            else {h_sum->Add(h_tmp);}
          }
          
        }
      }
      
      if(h_sum == 0) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; }
      int nofbins = h_sum->GetNbinsX();
      
      for(int i=0; i<nofbins; i++)
      {
        double bin_content = h_sum->GetBinContent(i+1); //cout<<"bin "<<i+1<<endl; cout<<"initial content = "<<bin_content<<endl;
        double new_bin_content = therand.Poisson(bin_content);// cout<<"new content = "<<new_bin_content<<endl;
        h_sum->SetBinContent(i+1, new_bin_content);
        h_sum->SetBinError(i+1, sqrt(new_bin_content)); //Poissonian error
      }
      
      pseudodata_file->cd();
      string output_histo_name = template_name + "_" + channel_list[ichan] + "_data_obs"; // TO FIX
      h_sum->Write(output_histo_name.c_str(), TObject::kOverwrite);
      
    } // chan
    
    pseudodata_file->Close();
    
    cout<<"--- Done with generation of pseudo-data"<<endl<<endl;
    
    delete pseudodata_file;
    
  }
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

///////////////////////////////////// BDT cyt /////////////////////////////////////////
vector<double> BDTCUT(string region, string coupling){
  cout << "Determine BDT cut " << endl;
  
  
  string bdtinput_name = placeOutputReading+"/Reader_" + template_name + ".root";
  TFile* bdt_file = TFile::Open( bdtinput_name.c_str(), "READ" );
  
  cout << bdtinput_name << endl;
  
  
  
  
  string histo_name = "";
  string template_fake_name = "";
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  
  double cut_eee, cut_eeu, cut_uuu, cut_uue;
  cut_eee = cut_uue = cut_eeu = cut_uuu = -1.;
  for(int ichan=0; ichan<channel_list.size(); ichan++)
  {
    string channel = channel_list[ichan];
    TH1F *h_sum_bkg(0), *h_sum_sig(0), *h_tmp(0);
    histo_name = "";
    if(channel.find("uuu") != std::string::npos || channel.find("eeu") != std::string::npos) {template_fake_name = "FakeMu";}
    else {template_fake_name = "FakeEl";}
    
    
    for(int isample = 0; isample < datasets.size(); isample++)
    {
      string dataSetName = datasets[isample]->Name();
      //cout << dataSetName << endl;
      if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // no signal in data
      if(datasets[isample]->Name().find("fake")==std::string::npos && datasets[isample]->Name().find("FCNC")==std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << endl;
        h_tmp = 0;
        histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
        if(h_tmp->GetEntries() != 0){
          if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
          else {h_sum_bkg->Add(h_tmp);}}
      }
      else if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << " is signal "<< endl;
        h_tmp = 0;
        histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = template_name + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        if(h_tmp->GetEntries() != 0){
          if(h_sum_sig == 0) {h_sum_sig = (TH1F*) h_tmp->Clone();}
          else {h_sum_sig->Add(h_tmp);}
        }
      }
      else{
        h_tmp = 0;
        histo_name = template_name + "_" + channel_list[ichan] + "_" + template_fake_name;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<histo_name<<" : not found"<<endl;}
        else
        {
          h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
          //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
          if(h_tmp->GetEntries() != 0){
            if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
            else {h_sum_bkg->Add(h_tmp);}
          }
        }
        
      }
    }
    //cout << "h_sum_bkg " << h_sum_bkg << " h_sum_sig "<< h_sum_sig << endl;
    if(h_sum_bkg == 0 || h_sum_sig == 0) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; }
    //S+B histogram
    TH1F* h_total = (TH1F*) h_sum_bkg->Clone();
    h_total->Add(h_sum_sig);
    
    //Normalization
  h_total->Scale(1/h_total->Integral()); //Integral() : Return integral of bin contents in range [binx1,binx2] (inclusive !)
    h_sum_bkg->Scale(1/h_sum_bkg->Integral());
    h_sum_sig->Scale(1/h_sum_sig->Integral());
    
    double sig_over_total = 100; //initialize to unreasonable value
    int bin_cut = -1; //initialize to false value
    int nofbins = h_total->GetNbinsX();
    
    for(int ibin=nofbins; ibin>0; ibin--)
    {
      //Search the bin w/ lowest sig/total, while keeping enough bkg events (criterion needs to be optimized/tuned)
      if( (h_sum_sig->Integral(1, ibin) / h_total->Integral(1, ibin)) < sig_over_total && (h_sum_bkg->Integral(1, ibin) / h_sum_bkg->Integral()) >= 0.6 )
      {
        bin_cut = ibin;
        sig_over_total = h_sum_sig->Integral(1, bin_cut) / h_total->Integral(1,bin_cut);
      }
    }
    
    double cut = h_total->GetBinLowEdge(bin_cut+1); //Get the BDT cut value to apply to create a BDT CR control tree
    double min = TMath::Min(h_sum_bkg->GetMinimum(), h_sum_sig->GetMinimum());
    double max = TMath::Max(h_sum_bkg->GetMaximum(), h_sum_sig->GetMaximum());
    //Create plot to represent the cut on BDT
    TCanvas* c = new TCanvas("c", "Signal vs Background");
    gStyle->SetOptStat(0);
    h_sum_bkg->GetYaxis()->SetRange(min,1.049*max);
    h_sum_bkg->GetXaxis()->SetTitle("BDT Discriminant");
    h_sum_bkg->SetTitle("Signal vs Background");
    h_sum_bkg->SetLineColor(kBlue);
    h_sum_sig->SetLineColor(kRed-3);
    h_sum_bkg->SetLineWidth(3);
    h_sum_sig->SetLineWidth(3);
    //h_sum_sig->Scale(100);
    h_sum_bkg->Draw("HIST");
    h_sum_sig->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.7,0.75,0.88,0.85);
    leg->SetHeader("");
    leg->AddEntry(h_sum_sig,"Signal","L");
    leg->AddEntry(h_sum_bkg,"Background","L");
    leg->Draw();
    //Draw vertical line at cut value
    
    TLine* l = new TLine(cut,min, cut, 1.049*max);
    l->SetLineWidth(3);
    l->SetLineColor(921);
    l->Draw("");
    TString outputsaving = "CutPlots";
    mkdir(outputsaving,0777);
    outputsaving += "/"+region + "_" + coupling + "_";
    outputsaving += "Signal_Background_BDT_"+channel+".png";
    c->SaveAs(outputsaving.Data());
    
    //Cout some results
    cout<<"---------------------------------------"<<endl;
    cout<<"* Cut Value = "<<cut<<endl;
    cout<<"-> BDT_CR defined w/ all events inside bins [1 ; "<<bin_cut<<"] of the BDT distribution!"<<endl<<endl;
    cout<<"* Signal integral = "<<h_sum_sig->Integral(1, bin_cut)<<" / Total integral "<<h_total->Integral(1, bin_cut)<<endl;
    cout<<"Signal contamination in CR --> Sig/Total = "<<sig_over_total<<endl;
    cout<<"Bkg(CR) / Bkg(Total) = "<<h_sum_bkg->Integral(1,bin_cut) / h_sum_bkg->Integral()<<endl;
    cout<<"---------------------------------------"<<endl<<endl;
    
    //for(int i=0; i<h_sig->GetNbinsX(); i++) {cout<<"bin content "<<i+1<<" = "<<h_sig->GetBinContent(i+1)<<endl;} //If want to verify that the signal is computed correctly
     delete c; delete leg; delete l;
    if(channel.find("uuu") != std::string::npos) cut_uuu = cut;
    else if(channel.find("eeu") != std::string::npos) cut_eeu = cut;
    else if(channel.find("eee") != std::string::npos) cut_eee = cut;
    else if(channel.find("uue") != std::string::npos) cut_uue = cut;
    
    delete h_sum_sig;
    delete h_sum_bkg;
    delete h_tmp;
  }
  bdt_file->Close();
  vector<double> v_return;
  v_return.push_back(cut_uuu);
  v_return.push_back(cut_uue);
  v_return.push_back(cut_eeu);
  v_return.push_back(cut_eee);
  return v_return;
  
  
}



///////////////////////////////////// INIT PLOTS /////////////////////////////////////////
void Init1DPlots(){
  TH1::SetDefaultSumw2();
  
  // histo1D["eventID"] = new TH1F("eventID","eventID",5802564204,207137,5802564203);
  
}

void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair){
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
    
    
    //cout << "init " << (prefix+"_BDT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), 10, -1,1, "BDT");
    
    MSPlot[(prefix+"_MVA_dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRWlepb_"+decaystring).c_str(),30,0, 6, "dR(l_{W}b)");
    MSPlot[(prefix+"_MVA_charge_asym_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_charge_asym_"+decaystring).c_str(),40, -4, 4, "Q(l_{W})|#eta(W)|");
    MSPlot[(prefix+"_MVA_dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZWlep_"+decaystring).c_str(),30,0, 6, "dR(Z,l_{W})");
    
    
    if(istoppair){
      MSPlot[(prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str(),50,-1, 1, "CvsL 2nd Highest pt jet");
      MSPlot[(prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str(),50, -1, 1, "CvsL Highest pt jet");
      MSPlot[(prefix+"_MVA_dRZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dRZc_"+decaystring).c_str(),30,0, 6, "dR(Z,Ljet)");
      MSPlot[(prefix+"_MVA_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mlb_"+decaystring).c_str(),100, 0, 500, "M(l_{W}b)");
      
    }
    else{
      MSPlot[(prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str(),20, 0.4, 1, "CSVv2 Highest pt jet");
      MSPlot[(prefix+"_MVA_Zboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"_MVA_Zboson_pt_"+decaystring).c_str(), 25,0, 500, "p_{T} (Z)[GeV]");
      MSPlot[(prefix+"_MVA_dPhiWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_dPhiWlepb_"+decaystring).c_str(),40,-4, 4, "dphi(l_{W}b)");
      MSPlot[(prefix+"_MVA_mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_MVA_mlb_"+decaystring).c_str(),50, 0, 500, "M(l_{W}b)");
    }
    
    
    
  }
  
  
  
  
  
}

void Init2DPlots(){
  TH2::SetDefaultSumw2();
  //histo2D["CosTheta"]= new TH2F("CosTheta", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);
}

///////////////////////////////////// INIT TREES /////////////////////////////////////////
void InitTree(TTree* tree, bool isData, bool istoppair){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  if(istoppair){
    tree->SetBranchAddress("MVA_cdiscCvsL_jet_1", &MVA_cdiscCvsL_jet_1, &b_MVA_cdiscCvsL_jet_1);
    tree->SetBranchAddress("MVA_cdiscCvsL_jet_0", &MVA_cdiscCvsL_jet_0, &b_MVA_cdiscCvsL_jet_0);
    tree->SetBranchAddress("MVA_dRZc", &MVA_dRZc, &b_MVA_dRZc);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
    
  }
  else if(!istoppair){
    tree->SetBranchAddress("MVA_Zboson_pt", &MVA_Zboson_pt, &b_MVA_Zboson_pt);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    
    
  }
  
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  tree->SetBranchAddress("MVA_weight", &MVA_weight, &b_MVA_weight);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress("BDT", &BDT, &b_BDT);
  tree->SetBranchAddress("MVA_EqLumi", &MVA_EqLumi, &b_MVA_EqLumi);
  tree->SetBranchAddress("MVA_weight", &MVA_weight, &b_MVA_weight);
  if(!isData){
    tree->SetBranchAddress("MVA_weight_puSF_up", &MVA_weight_puSF_up, &b_MVA_weight_puSF_up);
    tree->SetBranchAddress("MVA_weight_puSF_down", &MVA_weight_puSF_down, &b_MVA_weight_puSF_down);
    tree->SetBranchAddress("MVA_weight_electronSF_up", &MVA_weight_electronSF_up, &b_MVA_weight_electronSF_up);
    tree->SetBranchAddress("MVA_weight_electronSF_down", &MVA_weight_electronSF_down, &b_MVA_weight_electronSF_down);
    tree->SetBranchAddress("MVA_weight_muonSF_up", &MVA_weight_muonSF_up, &b_MVA_weight_muonSF_up);
    tree->SetBranchAddress("MVA_weight_muonSF_down", &MVA_weight_muonSF_down, &b_MVA_weight_muonSF_down);
    tree->SetBranchAddress("MVA_weight_btagSF_cferr1_up", &MVA_weight_btagSF_cferr1_up, &b_MVA_weight_btagSF_cferr1_up);
    tree->SetBranchAddress("MVA_weight_btagSF_cferr1_down", &MVA_weight_btagSF_cferr1_down, &b_MVA_weight_btagSF_cferr1_down);
    tree->SetBranchAddress("MVA_weight_btagSF_cferr2_up", &MVA_weight_btagSF_cferr2_up, &b_MVA_weight_btagSF_cferr2_up);
    tree->SetBranchAddress("MVA_weight_btagSF_cferr2_down", &MVA_weight_btagSF_cferr2_down, &b_MVA_weight_btagSF_cferr2_down);
    tree->SetBranchAddress("MVA_weight_btagSF_hf_up", &MVA_weight_btagSF_hf_up, &b_MVA_weight_btagSF_hf_up);
    tree->SetBranchAddress("MVA_weight_btagSF_hf_down", &MVA_weight_btagSF_hf_down, &b_MVA_weight_btagSF_hf_down);
    tree->SetBranchAddress("MVA_weight_btagSF_hfstats1_up", &MVA_weight_btagSF_hfstats1_up, &b_MVA_weight_btagSF_hfstats1_up);
    tree->SetBranchAddress("MVA_weight_btagSF_hfstats1_down", &MVA_weight_btagSF_hfstats1_down, &b_MVA_weight_btagSF_hfstats1_down);
    tree->SetBranchAddress("MVA_weight_btagSF_hfstats2_up", &MVA_weight_btagSF_hfstats2_up, &b_MVA_weight_btagSF_hfstats2_up);
    tree->SetBranchAddress("MVA_weight_btagSF_hfstats2_down", &MVA_weight_btagSF_hfstats2_down, &b_MVA_weight_btagSF_hfstats2_down);
    tree->SetBranchAddress("MVA_weight_btagSF_lf_up", &MVA_weight_btagSF_lf_up, &b_MVA_weight_btagSF_lf_up);
    tree->SetBranchAddress("MVA_weight_btagSF_lf_down", &MVA_weight_btagSF_lf_down, &b_MVA_weight_btagSF_lf_down);
    tree->SetBranchAddress("MVA_weight_btagSF_lfstats1_up", &MVA_weight_btagSF_lfstats1_up, &b_MVA_weight_btagSF_lfstats1_up);
    tree->SetBranchAddress("MVA_weight_btagSF_lfstats1_down", &MVA_weight_btagSF_lfstats1_down, &b_MVA_weight_btagSF_lfstats1_down);
    tree->SetBranchAddress("MVA_weight_btagSF_lfstats2_up", &MVA_weight_btagSF_lfstats2_up, &b_MVA_weight_btagSF_lfstats2_up);
    tree->SetBranchAddress("MVA_weight_btagSF_lfstats2_down", &MVA_weight_btagSF_lfstats2_down, &b_MVA_weight_btagSF_lfstats2_down);
    
  }
}

void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData, bool toppair){
  //cout << "fill plots" << endl;
  string decaystring = "";
  Double_t eventW = 1.;
  //if(!isData) eventW = MVA_weight;
  eventW = weight;
  
  if(datasets[d]->Name().find("fake")!=std::string::npos) eventW *= 0.0001;
  
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
    string sregion = prefix;
    MSPlot[(sregion+"_MVA_dRWlepb_"+decaystring).c_str()]->Fill(MVA_dRWlepb , datasets[d], true,eventW);
    MSPlot[(sregion+"_MVA_mlb_"+decaystring).c_str()]->Fill(MVA_mlb , datasets[d], true,eventW);
    MSPlot[(sregion+"_MVA_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym , datasets[d], true,eventW);
    MSPlot[(sregion+"_MVA_dRZWlep_"+decaystring).c_str()]->Fill(MVA_dRZWlep , datasets[d], true, eventW);
    
    if(toppair){
      MSPlot[(sregion+"_MVA_cdiscCvsL_jet_1_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_1 , datasets[d], true, eventW);
      MSPlot[(sregion+"_MVA_cdiscCvsL_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_0 , datasets[d], true, eventW);
      MSPlot[(sregion+"_MVA_dRZc_"+decaystring).c_str()]->Fill(MVA_dRZc , datasets[d], true, eventW);
    }
    else{
      
      MSPlot[(sregion+"_MVA_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0 , datasets[d], true, eventW);
      MSPlot[(sregion+"_MVA_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, eventW);
      MSPlot[(sregion+"_MVA_dPhiWlepb_"+decaystring).c_str()]->Fill(MVA_dPhiWlepb , datasets[d], true,eventW);
      
      
      
    }
    
    
    
    
    
    
    
  }
  
  
  
  //cout << "end plot filling" << endl;
}






