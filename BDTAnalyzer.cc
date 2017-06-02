
#include <TH1.h>
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
#include <map>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"

//#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/LHAPDF.h"
//#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/Reweighting.h"

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH1F*> histo1DPDF;
map<string,TH1F*> histo1DSys;
map<string,TH1F*> histo1DSysMTW;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlotMTW;
////////////////////////////////// mapping ///////////////////////////////
map<string,TFile*> tFileMap;
TFile *fin;
map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;
vector < Dataset* > datasets;
vector < Dataset* > datasetsbf;
////////////////////////////////// functions ////////////////////////////////////////////
// bookkeeping
std::vector<std::string> split(const std::string &text,  char sep) ;
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
string intToStr (int number);
double maximumValue(vector<double> array);
double minimumValue(vector<double> array);
// initialisations

void InitMSPlotsMTW(string prefix, vector <int> decayChannels);
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair, bool isZut);
void InitCalculatePDFWeightHisto(string dataSetName);
void InitMTWShapeHisto(string dataSetName, string systematic, int isys,  vector <int> decayChannels);
void InitSystematicHisto(string dataSetName, string systematic, int isys, bool doMTWtemplate);
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut);
void InitAnalyzerTree(TTree* tree);
void Init1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels);
// functions
vector<double> BDTCUT(string region, string coupling);
void CalculatePDFWeight(string dataSetName, double BDT, double MVA_weight_nom, int MVA_channel);
void FillMTWPlots(int d, string postfix, vector <int> decayChannels, double weight_, int MVA_channel);
void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData, bool toppair, double weight_, int MVA_channel);
void GetPDFEnvelope(string dataSetName);
void FillMTWShapeHisto(string dataSetName, string systematic, double weight_,int isys, int MVA_channel, vector <int> decayChannels);
void FillSystematicHisto(string dataSetName, string systematic, double weight_, int isys, bool doMTWtemplate);
void Fill1DHisto(string dataSetName,string systematic, bool istoppair, bool isZut, vector <int> decayChannels, double weight_, int MVA_channel);

//////////////////////////////// settings ////////////////////////////////
bool makePlots = false;
bool doMTWtemplate = false;
bool PlotSystematics = false;
bool doPseudoData = false;
bool doSystematics = false;
int channel = -999;
bool datafound = false;
bool testing = false;
bool doZut = false;
bool toppair = false;
bool addData = false;
bool DetermineCut = false;
bool CalculateSign = false;
bool doPDFunc  = false;
bool PlotMVAvars = false;
string placeNtup = "";
string tempstring = "";
string systematic = "";
string decaystring = "";
string coupling = "Zct";
string region = "singletop";
string placeOutputReading = "";
string combinetemplate_filename = "";
Double_t Luminosity = 36000.;
string dataSetName = "";
bool isData = false;
string tTreeName = "";
string postfix = "";
string output_histo_name = "";
string ntupleFileName ="";
int nbin = 50;
int nbinMTW = 20;
double endMTW = 300.;
int nEntries = -1;
int scaleNP = 10;

//////////////////////////// branches //////////////////////////////
// Declaration of leaf types
Float_t         MVA_Zboson_pt;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_charge_asym;
Float_t         MVA_bdiscCSVv2_jet_0;

Float_t         MVA_cdiscCvsB_jet_1;
Float_t         MVA_cdiscCvsB_jet_0;
Float_t         MVA_cdiscCvsL_jet_1;
Float_t         MVA_cdiscCvsL_jet_0;
Float_t         MVA_dRZc;
Float_t         MVA_dRWlepb;
Float_t         MVA_dRZWlep;
Float_t         MVA_mlb;
Float_t         MVA_mWt2;
Float_t         MVA_mWt;
Float_t         MVA_FCNCtop_M;
Float_t           MVA_nJets_CharmL;
Float_t           MVA_NJets_CSVv2M;
Float_t         MVA_region;
Double_t       MVA_x1;
Double_t       MVA_x2;
Int_t       MVA_id1;
Int_t       MVA_id2;
Double_t       MVA_q;
//Float_t         MVA_weight;
Double_t       MVA_weight_nom;
Int_t         MVA_channel;
Double_t       MVA_BDT;
Double_t        MVA_EqLumi;
Double_t        MVA_Luminosity;
Double_t       MVA_weight_puSF_up;
Double_t       MVA_weight_puSF_down;
Double_t       MVA_weight_electronSF_up;
Double_t       MVA_weight_electronSF_down;
Double_t       MVA_weight_muonSF_up;
Double_t       MVA_weight_muonSF_down;
Double_t       MVA_weight_btagSF_cferr1_up;
Double_t       MVA_weight_btagSF_cferr1_down;
Double_t       MVA_weight_btagSF_cferr2_up;
Double_t       MVA_weight_btagSF_cferr2_down;
Double_t       MVA_weight_btagSF_hf_up;
Double_t       MVA_weight_btagSF_hf_down;
Double_t       MVA_weight_btagSF_hfstats1_up;
Double_t       MVA_weight_btagSF_hfstats1_down;
Double_t       MVA_weight_btagSF_hfstats2_up;
Double_t       MVA_weight_btagSF_hfstats2_down;
Double_t       MVA_weight_btagSF_lf_up;
Double_t       MVA_weight_btagSF_lf_down;
Double_t       MVA_weight_btagSF_lfstats1_up;
Double_t       MVA_weight_btagSF_lfstats1_down;
Double_t       MVA_weight_btagSF_lfstats2_up;
Double_t       MVA_weight_btagSF_lfstats2_down;

// List of branches
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_0;   //!

TBranch        *b_MVA_cdiscCvsB_jet_1;   //!
TBranch        *b_MVA_cdiscCvsB_jet_0;   //!
TBranch        *b_MVA_cdiscCvsL_jet_1;   //!
TBranch        *b_MVA_cdiscCvsL_jet_0;   //!
TBranch        *b_MVA_dRZc;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_mlb;   //!
TBranch        *b_MVA_mWt2;   //!
TBranch        *b_MVA_mWt;   //!
TBranch        *b_MVA_FCNCtop_M;   //!
TBranch        *b_MVA_nJets_CharmL;   //!
TBranch        *b_MVA_NJets_CSVv2M;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_MVA_x1;   //!
TBranch        *b_MVA_x2;   //!
TBranch        *b_MVA_id1;   //!
TBranch        *b_MVA_id2;   //!
TBranch        *b_MVA_q;   //!
TBranch        *b_MVA_weight;   //!
TBranch        *b_MVA_channel;   //!
TBranch        *b_MVA_BDT;   //!
TBranch        *b_MVA_EqLumi;   //!
TBranch        *b_MVA_Luminosity;   //!
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
TBranch        *b_MVA_weight_nom; //


map<string,TH1F*> histo1DMTW;
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
  
  
  int testnr = -1;
  //////////// Settings of the analysis //////////////////
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./AnalyzerBDT [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   Zct / Zut: make plots for this one" << endl;
      std::cout << "   toppair: make plots for this one" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. MVAoutput/placeNtup " << endl;
      std::cout << "   Test nb: loop over nb events" << endl;
      std::cout << "   Data: add CRcut" << endl;
      std::cout << "   PSdata: generate pseudo data" << endl;
      std::cout << "   Systematics: loop over systematics" << endl;
      std::cout << "   doPDFunc: calculate PDF unc" << endl;
      std::cout << "   PlotSystematics: make sys plots fo WZ" << endl;
      std::cout << "   PlotMVAvars: plot mva vars" << endl;
      std::cout << "  doMTWtemplate make mtw templates" << endl;
      return 0;
    }
    if(string(argv[i]).find("doMTWtemplate")!=std::string::npos){
      doMTWtemplate = true;
      cout << "******* making MTW templates *********" << endl;
    }
    if(string(argv[i]).find("doPDFunc")!=std::string::npos){
      doPDFunc = true;
      cout << "******* calculating pdf uncertainties *********" << endl;
    }
    if(string(argv[i]).find("Systematics")!=std::string::npos) {
      doSystematics= true;
    }
    if(string(argv[i]).find("PlotSystematics")!=std::string::npos) {
      PlotSystematics= true;
      //cout << "plotting sys" << endl;
    }
    if(string(argv[i]).find("PSdata")!=std::string::npos) {
      doPseudoData = true;
    }
    if(string(argv[i]).find("Data")!=std::string::npos) {
      addData = true;
    }
    if(string(argv[i]).find("toppair")!=std::string::npos) {
      toppair = true;
      region = "toppair";
    }
    if(string(argv[i]).find("singletop")!=std::string::npos) {
      toppair = false;
      region = "singletop";
    }
    if(string(argv[i]).find("Zut")!=std::string::npos) {
      doZut = true;
      coupling = "Zut";
    }
    if(string(argv[i]).find("Zct")!=std::string::npos) {
      doZut = false;
      coupling = "Zct";
    }
    if(string(argv[i]).find("Test")!=std::string::npos) {
      testing = true;
      testnr = strtol(argv[i+1], NULL, 10);
      i++;
    }
    if(string(argv[i]).find("Ntup")!=std::string::npos) {
      placeNtup = argv[i+1];
      i++;
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    if(string(argv[i]).find("DetermineCut")!=string::npos) {
      DetermineCut= true;
    }
    if(string(argv[i]).find("Significane")!=string::npos) {
      CalculateSign = true;
    }
    if(string(argv[i]).find("PlotMVAvars")!=string::npos) {
      PlotMVAvars= true;
    }
  }
  string xmlFileName = "";
  xmlFileName = "config/Run2TriLepton_samples_analy.xml" ;
  const char* xmlFile = xmlFileName.c_str();
  cout << " - Using config file " << xmlFile << endl;
  cout << " - Using mvatrees of " << placeNtup << endl;
  placeOutputReading = "MVAoutput/";
  mkdir(placeOutputReading.c_str(), 0777);
  placeOutputReading += "outputtemplates";
  mkdir(placeOutputReading.c_str(), 0777);
  if(!doMTWtemplate) placeOutputReading += "/" + coupling + "_" + region;
  else placeOutputReading += "/MTWtemplate" ;
  mkdir(placeOutputReading.c_str(), 0777);
  combinetemplate_filename = placeOutputReading+"/Reader_" + coupling +"_" + region + ".root";
  if(doMTWtemplate) combinetemplate_filename = placeOutputReading+"/Reader_"+coupling+"_MTW.root";
  cout <<" - Combine templates stored at " << combinetemplate_filename.c_str() << endl;
  
  std::vector < int>  decayChannels = {0,1,2,3,-9}; // uuu uue eeu eee all
  vector <string> thesystlist;
  vector <string> thesystlistnames;
  thesystlistnames.clear();
  thesystlist.clear();
  thesystlist.push_back(""); // nominal
  if(doSystematics){
    cout << "pushing back systematics" << endl;
    thesystlist.push_back("puSFDown");
    thesystlist.push_back("electronSFDown");
    thesystlist.push_back("muonSFDown");
    thesystlist.push_back("btagSF_cferr1Down");
    thesystlist.push_back("btagSF_cferr2Down");
    thesystlist.push_back("btagSF_hfDown");
    thesystlist.push_back("btagSF_hfstats1Down");
    thesystlist.push_back("btagSF_hfstats2Down");
    thesystlist.push_back("btagSF_lfDown");
    thesystlist.push_back("btagSF_lfstats1Down");
    thesystlist.push_back("btagSF_lfstats2Down");
    thesystlist.push_back("JERUp");
    thesystlist.push_back("JESUp");
    
    thesystlist.push_back("puSFUp");
    thesystlist.push_back("electronSFUp");
    thesystlist.push_back("muonSFUp");
    thesystlist.push_back("btagSF_cferr1Up");
    thesystlist.push_back("btagSF_cferr2Up");
    thesystlist.push_back("btagSF_hfUp");
    thesystlist.push_back("btagSF_hfstats1Up");
    thesystlist.push_back("btagSF_hfstats2Up");
    thesystlist.push_back("btagSF_lfUp");
    thesystlist.push_back("btagSF_lfstats1Up");
    thesystlist.push_back("btagSF_lfstats2Up");
    
    
    thesystlist.push_back("JERDown");
    
    thesystlist.push_back("JESDown");
  }
  cout << "Number of systematics " << thesystlist.size() << endl;
  //for plotting
  thesystlistnames.push_back("puSF");
  thesystlistnames.push_back("electronSF");
  thesystlistnames.push_back("muonSF");
  thesystlistnames.push_back("btagSF_cferr1");
  thesystlistnames.push_back("btagSF_cferr2");
  thesystlistnames.push_back("btagSF_hf");
  thesystlistnames.push_back("btagSF_hfstats1");
  thesystlistnames.push_back("btagSF_hfstats2");
  thesystlistnames.push_back("btagSF_lf");
  thesystlistnames.push_back("btagSF_lfstats1");
  thesystlistnames.push_back("btagSF_lfstats2");
  thesystlistnames.push_back("JES");
  thesystlistnames.push_back("JER");
  
  
  ///////////////:  load datasets
  datasets.clear();
  datasetsbf.clear();
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasetsbf, xmlFile);
  cout << "datasets loaded " << datasetsbf.size() << " samples" <<  endl;
  for (int d = 0; d < datasetsbf.size(); d++){   //Loop through datasets to get lumi setting
    
    dataSetName = datasetsbf[d]->Name();
   // cout << "sample " << dataSetName << endl;
    if( dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
      Luminosity = datasetsbf[d]->EquivalentLumi();
      cout << " - lumi set to " << Luminosity << endl;
    }
    
    if((dataSetName.find("Zct")!=std::string::npos || dataSetName.find("zct")!=std::string::npos) && doZut && !doMTWtemplate){
      cout << " - removing " << dataSetName << " from samples" << endl;
      continue;
    }
    else if ((dataSetName.find("Zut")!=std::string::npos || dataSetName.find("zut")!=std::string::npos) && !doZut && !doMTWtemplate){
      cout << " - removing " << dataSetName << " from samples" << endl;
      continue;
    }
    else{
      
      if( (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ) && doPseudoData ) {
        cout << " - removing " << dataSetName << " from samples" << endl;
      }
      else if(dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
        datafound = true;
        //cout << "pushing back " << dataSetName << endl;
        datasets.push_back(datasetsbf[d]);
      }
      else {
        //cout << "pushing back " << dataSetName << endl;
        datasets.push_back(datasetsbf[d]);
      }
    }
  }
  
  cout << datasets.size() << " samples will be used " << endl;
  ///////////////// Initialisation ////////////////////

  //if(toppair) nbin = 50;
  if(makePlots && !doMTWtemplate){
    for(int isys = 0; isys < thesystlist.size() ; isys++){
      systematic = thesystlist[isys];
      tempstring = region + "_"+coupling;
      if(isys != 0 ) tempstring += "_"+ systematic;
      InitMSPlots(tempstring, decayChannels, toppair, doZut);
    }
  }
  
  if(makePlots && doMTWtemplate){
    for(int isys = 0; isys < thesystlist.size() ; isys++){
      systematic = thesystlist[isys];
      if(isys != 0 ) tempstring = "_" + systematic;
      InitMSPlotsMTW(tempstring, decayChannels);
    }
  }
  
  
  ///////////// General function //////////
  if(DetermineCut && !doMTWtemplate){
    vector<double> v_cut = BDTCUT(region, coupling);
    return 0;
  }
  
  //////////// START LOOPING ON SYS - DATASETS - EVENTS //////////
  ///*****************///
  ///   MAIN CODE  ///
  ///*****************///
  ntupleFileName = placeNtup;
  if(!doMTWtemplate) fin = new TFile((ntupleFileName).c_str(),"READ");
  bool onlynomforsys = false;
  int WZregionEntries = 0;
 
  
  TH1::SetDefaultSumw2();
  TH1F* hist_WZ = new TH1F("MTW_WZ","trans. mass W boson in WZ region: WZ (GeV)",           nbinMTW, 0., endMTW);
  TH1F* hist_TT_FCNC = new TH1F("MTW_TT_FCNC","transv. mass W boson in TTSR (GeV)",           nbinMTW, 0., endMTW);
  TH1F* hist_fakes = new TH1F("MTW_fakes","transv. mass W boson WZ region: fakes (GeV)",           nbinMTW, 0., endMTW);
  histo1DMTW["MTW_WZ"] = hist_WZ ;
  histo1DMTW["MTW_TT_FCNC"] = hist_TT_FCNC ;
  histo1DMTW["MTW_fakes"] = hist_fakes ;
  
  for(int isys = 0; isys < thesystlist.size() ; isys++){
    systematic = thesystlist[isys];
    cout << "looking at " << systematic << " systematics " << endl;

    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
      // settings
      isData = false;
      onlynomforsys = false;
      dataSetName = datasets[d]->Name();
      if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
        isData = true;
      }
      
      if ((dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos || dataSetName.find("fake")!=std::string::npos   ) && isys != 0){
        onlynomforsys = true;
      }
      
      
      postfix = "";
      if(systematic.find("JESDown")!=std::string::npos) postfix = "_JESdown";
      else if(systematic.find("JESUp")!=std::string::npos) postfix = "_JESup";
      else if(systematic.find("JERDown")!=std::string::npos) postfix = "_JERdown";
      else if(systematic.find("JERUp")!=std::string::npos) postfix = "_JERup";
      else postfix = "";
      if(onlynomforsys) postfix = "";
      
      if(doMTWtemplate) {
        ntupleFileName = placeNtup + "MVA_tree_" + dataSetName + postfix + ".root";
        tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
      }
      
      if(!doMTWtemplate) tTreeName = "Control_"+dataSetName + postfix;
      else tTreeName = "mvatree" + postfix;
      /// Get data
      cout << "   treename " << tTreeName << " from " << ntupleFileName <<  endl;
      if(doMTWtemplate){
        tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      }
      else if(!doMTWtemplate) tTree[dataSetName.c_str()] = (TTree*)fin->Get(tTreeName.c_str());
      nEntries = -1;
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries: " << nEntries << endl;
      
      // Initialise tree
      if(!doMTWtemplate) InitTree(tTree[dataSetName.c_str()], isData, toppair, doZut);
      else if(doMTWtemplate) InitAnalyzerTree(tTree[dataSetName.c_str()]);
      
      
      // Initialise plots
      if(!doMTWtemplate && PlotMVAvars && isys == 0){
        if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("TT_FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && toppair ) Init1DHisto(dataSetName, systematic, toppair, doZut, decayChannels);
        else if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("ST_FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && !toppair ) Init1DHisto(dataSetName, systematic, toppair, doZut, decayChannels);
      }
      // initialise combine output histograms
      TH1::SetDefaultSumw2();
      //cout << "create template histo" << endl;
      TH1F *hist_uuu(0);
      if(!doMTWtemplate) hist_uuu = new TH1F( (coupling + "_" + region+"_uuu").c_str(),           (coupling + "_" + region+"_uuu").c_str(),           nbin, -1, 1 );
      else hist_uuu = new TH1F( (coupling + "_mTW_uuu").c_str(),           (coupling + "_mTW_uuu").c_str(),           nbinMTW,0, endMTW );
      TH1F *hist_uue(0);
      if(!doMTWtemplate) hist_uue = new TH1F( (coupling + "_" + region+"_uue").c_str(),           (coupling + "_" + region+"_uue").c_str(),           nbin, -1, 1 );
      else hist_uue = new TH1F( (coupling + "_mTW_uue").c_str(),           (coupling + "_mTW_uue").c_str(),           nbinMTW,0, endMTW );
      TH1F *hist_eeu(0);
      if(!doMTWtemplate) hist_eeu = new TH1F( (coupling + "_" + region+"_eeu").c_str(),           (coupling + "_" + region+"_eeu").c_str(),           nbin, -1, 1 );
      else hist_eeu = new TH1F( (coupling + "_mTW_eeu").c_str(),           (coupling + "_mTW_eeu").c_str(),           nbinMTW,0, endMTW );
      TH1F *hist_eee(0);
      if(!doMTWtemplate) hist_eee = new TH1F( (coupling + "_" + region+"_eee").c_str(),           (coupling + "_" + region+"_eee").c_str(),           nbin, -1, 1 );
      else hist_eee = new TH1F( (coupling + "_mTW_eee").c_str(),           (coupling + "_mTW_eee").c_str(),           nbinMTW,0, endMTW );
      //cout << "created template histo" << endl;
      /// Initialise WZ plots
      if((dataSetName.find("WZTo3LNu_3Jets_MLL50_80X")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos) && doPDFunc && !doMTWtemplate){
        InitCalculatePDFWeightHisto(dataSetName);
      }
      if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos) && PlotSystematics ){
        InitSystematicHisto(dataSetName, systematic, isys, doMTWtemplate);
      }
      
     
    
      if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && doMTWtemplate){
        // InitMTWShapeHisto(dataSetName, systematic, isys, decayChannels);
       
      }
      
      // safeties
      if(!doSystematics && isys != 0) continue;
      int endEvent =nEntries;
      if(testing){
        if(endEvent > testnr) endEvent = testnr;
      }
      /// loop on events
      double weight = 1.;
      WZregionEntries = 0;
      int WZregionEntriesuuu = 0;
      for (int ievt = 0; ievt < endEvent; ievt++)
      {
        if (ievt%100 == 0)
          std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
        
        
        /// Load event
        tTree[(dataSetName).c_str()]->GetEntry(ievt);
        if(datafound && MVA_BDT > -0.68 && !doMTWtemplate && !toppair && doZut){ continue;}
        else if(datafound && MVA_BDT > -0.6 && !doMTWtemplate && !toppair && !doZut){ continue;}
        else if(datafound && MVA_BDT > -0.12 && !doMTWtemplate && toppair && !doZut){ continue;}
        else if(datafound && MVA_BDT > -0.2 && !doMTWtemplate && toppair && doZut){ continue;}
        //if(isData) cout << "region " << MVA_region << endl;
        
        if(doMTWtemplate && MVA_region != 2){ continue ;} // only in WZ control region}
        else if(doMTWtemplate) { WZregionEntries++; }
        
        if(doMTWtemplate && MVA_region == 2 && MVA_channel == 0) { WZregionEntriesuuu++;}
        
        
        if(!isData && !onlynomforsys && isys != 0){
          //cout <<  "getting weight for " << systematic << endl;
          if(systematic.find("puSFUp")) weight = MVA_weight_puSF_up;
          else if(systematic.find("puSFDown")) weight = MVA_weight_puSF_down;
          else if(systematic.find("electronSFUp")) weight = MVA_weight_electronSF_up;
          else if(systematic.find("electronSFDown")) weight = MVA_weight_electronSF_down;
          else if(systematic.find("muonSFUp")) weight = MVA_weight_muonSF_up;
          else if(systematic.find("muonSFDown")) weight = MVA_weight_muonSF_down;
          else if(systematic.find("btagSF_cferr1Up")) weight = MVA_weight_btagSF_cferr1_up;
          else if(systematic.find("btagSF_cferr1Down")) weight = MVA_weight_btagSF_cferr1_down;
          else if(systematic.find("btagSF_cferr2Up")) weight = MVA_weight_btagSF_cferr2_up;
          else if(systematic.find("btagSF_cferr2Down")) weight = MVA_weight_btagSF_cferr2_down;
          else if(systematic.find("btagSF_lfUp")) weight = MVA_weight_btagSF_lf_up;
          else if(systematic.find("btagSF_lfDown")) weight = MVA_weight_btagSF_lf_down;
          else if(systematic.find("btagSF_hfUp")) weight = MVA_weight_btagSF_hf_up;
          else if(systematic.find("btagSF_hfDown")) weight = MVA_weight_btagSF_hf_down;
          else if(systematic.find("btagSF_hfstats1Up")) weight = MVA_weight_btagSF_hfstats1_up;
          else if(systematic.find("btagSF_hfstats1Down")) weight = MVA_weight_btagSF_hfstats1_down;
          else if(systematic.find("btagSF_hfstats2Up")) weight = MVA_weight_btagSF_hfstats2_up;
          else if(systematic.find("btagSF_hfstats2Down")) weight = MVA_weight_btagSF_hfstats2_down;
          else if(systematic.find("btagSF_lfstats2Up")) weight = MVA_weight_btagSF_lfstats2_up;
          else if(systematic.find("btagSF_lfstats2Down")) weight = MVA_weight_btagSF_lfstats2_down;
          
        }
        else if(!isData && !onlynomforsys ){ weight = MVA_weight_nom; }
        //if(Luminosity/MVA_Luminosity != 1. ) cout << "lumi "  << Luminosity << " while tuples are made with " << MVA_Luminosity << endl;
        
       // if(MVA_Luminosity != 0 && !isData) weight = (weight * Luminosity)/ MVA_Luminosity;
       // if(!datafound){ Luminosity = MVA_Luminosity; cout << "lumi set to " << Luminosity << endl; }
         if(dataSetName.find("fake")!=std::string::npos && (MVA_channel == 0 || MVA_channel == 2)){ weight *= 0.0802 * 0.0001 ;}
         if(dataSetName.find("fake")!=std::string::npos && (MVA_channel == 1 || MVA_channel == 3)){ weight *= 0.238 * 0.0001;}
        
        if( isData){ weight = 1.;}
        if(!doMTWtemplate){
          if(MVA_channel== 0) 		{hist_uuu->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 1) {hist_uue->Fill( MVA_BDT, weight);}
          else if(MVA_channel== 2) {hist_eeu->Fill( MVA_BDT, weight);}
          else if(MVA_channel == 3) {hist_eee->Fill( MVA_BDT, weight);}
        }
        else if(doMTWtemplate){
          if(MVA_channel== 0) 		{hist_uuu->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 1) {hist_uue->Fill( MVA_mWt, weight);}
          else if(MVA_channel== 2) {hist_eeu->Fill( MVA_mWt, weight);}
          else if(MVA_channel == 3) {hist_eee->Fill( MVA_mWt, weight);}
        }
        
        // for MS plots
        double weightMSPlot = weight;
        if(isData) weightMSPlot = Luminosity;
        /// Fill plots
        if(doPDFunc && !doMTWtemplate){
          if(dataSetName.find("WZTo3LNu_3Jets_MLL50_80X")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos) CalculatePDFWeight(dataSetName, MVA_BDT, MVA_weight_nom, MVA_channel);
        }
        if(PlotMVAvars  && isys == 0 && !doMTWtemplate){
          if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("TT_FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos )&& toppair) Fill1DHisto(dataSetName, systematic, toppair, doZut, decayChannels, weight, MVA_channel);
          else  if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("ST_FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && !toppair) Fill1DHisto(dataSetName, systematic, toppair, doZut, decayChannels, weight, MVA_channel);
        }
        
        
        if (makePlots && !doMTWtemplate)
        {
          //cout << "ievt " << ievt << endl;
          tempstring = region + "_"+coupling;
          if(isys != 0) tempstring += "_"+ systematic;
          FillGeneralPlots(d, tempstring, decayChannels, doZut, toppair, weightMSPlot, MVA_channel);
        }
        if (makePlots && doMTWtemplate)
        {
          if(isys != 0) tempstring = "_" + systematic;
         // if(isData) cout << "fill data " << endl;
          FillMTWPlots(d, tempstring, decayChannels, weightMSPlot, MVA_channel);
        }
        if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos) && PlotSystematics ){
          FillSystematicHisto(dataSetName, systematic, weight, isys, doMTWtemplate);
        }
       
        if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos)&& isys == 0 && doMTWtemplate) histo1DMTW["MTW_WZ"]->Fill(MVA_mWt2, weight);
        if(dataSetName.find("fake")!=std::string::npos && isys == 0 && doMTWtemplate) histo1DMTW["MTW_fakes"]->Fill(MVA_mWt2, weight);
        if(dataSetName.find("TT_FCNC")!=std::string::npos && isys == 0 && doMTWtemplate) histo1DMTW["MTW_TT_FCNC"]->Fill(MVA_mWt2, weight);
        
        
        if((dataSetName.find("WZTo3LNu")!=std::string::npos || dataSetName.find("WZJTo3LNu")!=std::string::npos || dataSetName.find("FCNC")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && doMTWtemplate){
          //cout << "filling FillMTWShapeHisto" << endl;
          //FillMTWShapeHisto(dataSetName, systematic, weight, isys, MVA_channel,decayChannels);
        }
        
      } // events
      
      cout << endl;
     if(doMTWtemplate) cout << "                WZ entries " << WZregionEntries << " uuu " << WZregionEntriesuuu << endl;
      /// Write combine histograms
      // --- Write histograms
      //cout << "DATASET " << dataSetName << " ISYS " << isys << endl;
      if((dataSetName.find("data")!=std::string::npos || dataSetName.find("fake")!=std::string::npos) && isys != 0) {
        delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
        continue;
      }
      if(doMTWtemplate && WZregionEntries == 0) {
        delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
        continue;
      }
      //cout << "making template" << endl;
      TFile* combinetemplate_file(0);
      if(d == 0 && isys == 0) combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "RECREATE" );
      else combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "UPDATE" );
      combinetemplate_file->cd();
      //cout << "opened " << combinetemplate_filename.c_str() << endl;
      //NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)] FIX ME
      output_histo_name = "";
      if(!doMTWtemplate){
        if (dataSetName.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uuu_FakeMu_80X"  ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uue_FakeEl_80X"  ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eeu_FakeMu_80X"  ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eee_FakeEl_80X"  ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uuu_"  + dataSetName ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_uue_"  + dataSetName ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eeu_"  + dataSetName ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = coupling + "_BDT_" + region+"_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = coupling + "_BDT_" + region+"_eee_"  + dataSetName ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
      }
      else if(doMTWtemplate){
        if (dataSetName.find("fake")!=std::string::npos ) //Last fake MC sample or data-driven fakes -> write fake histo w/ special name (for THETA)
        {
          if(isys!=0) output_histo_name = "MTW_uuu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "MTW_uuu_FakeMu_80X"  ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_uue_FakeEl_80X_"  + systematic ;
          else output_histo_name = "MTW_uue_FakeEl_80X"  ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eeu_FakeMu_80X_"  + systematic ;
          else output_histo_name = "MTW_eeu_FakeMu_80X"  ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eee_FakeEl_80X_"  + systematic ;
          else output_histo_name = "MTW_eee_FakeEl_80X"  ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
        else //If fakes are not considered, or if sample is not fake --> write directly !
        {
          if(isys!=0) output_histo_name = "MTW_uuu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_uuu_"  + dataSetName ;
          hist_uuu->SetTitle(output_histo_name.c_str());
          hist_uuu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_uue_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_uue_"  + dataSetName ;
          hist_uue->SetTitle(output_histo_name.c_str());
          hist_uue->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eeu_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_eeu_"  + dataSetName ;
          hist_eeu->SetTitle(output_histo_name.c_str());
          hist_eeu->Write(output_histo_name.c_str());
          if(isys!=0) output_histo_name = "MTW_eee_"  + dataSetName + "_" + systematic ;
          else output_histo_name = "MTW_eee_"  + dataSetName ;
          hist_eee->SetTitle(output_histo_name.c_str());
          hist_eee->Write(output_histo_name.c_str());
        }
      }
      combinetemplate_file->Close();
      //cout << "closed " << combinetemplate_filename.c_str() << endl;
      delete combinetemplate_file;
      delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
      
      if(doMTWtemplate) tFileMap[dataSetName.c_str()]->Close();
    }// datasets
    
    if(isys != 0) cout<<"Done with "<< systematic <<" systematic"<<endl;
    else cout<<"Done with nominal sample"<<endl;
    systematic = "";
  } // systematics
  
 //cout << "ENTRIES " << hist_WZ->GetEntries() << endl;
  
  cout << "ENTRIES " << histo1DMTW["MTW_WZ"]->GetEntries() << endl;
  if(!doMTWtemplate){
    fin->Close();
    delete fin;
  }
  ///*****************///
  ///   PDF envelope   ///
  ///*****************///
  
  
  if(doPDFunc && !doMTWtemplate)  GetPDFEnvelope("WZTo3LNu_3Jets_MLL50_80X");
  
  ///*****************///
  ///   Pseudodata   ///
  ///*****************///
  if(doPseudoData && !addData){
    cout << "generating pseudo data" << endl;
    TRandom3 therand(0); //Randomization
    
    string pseudodata_input_name = placeOutputReading+"/Reader_" + coupling + "_" + region + ".root";
    
    if(doMTWtemplate) pseudodata_input_name = placeOutputReading+"/Reader_"+coupling+"_MTW.root";
    TFile* pseudodata_file = TFile::Open( pseudodata_input_name.c_str(), "UPDATE" );
    
    cout << "Generating pseudo data from " << pseudodata_input_name << endl;
    
    
    
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
      if(channel_list[ichan] == "uuu" || channel_list[ichan] == "eeu") {template_fake_name = "FakeMu_80X";}
      else {template_fake_name = "FakeEl_80X";}
      
      
      
      
      for(int isample = 0; isample < datasets.size(); isample++)
      {
        string dataSetName = datasets[isample]->Name();
        // cout << dataSetName << endl;
        if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {continue; } // no signal in data
        else if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // safety
        else if(datasets[isample]->Name().find("fake")==std::string::npos){
          // cout << "  -- sample " << datasets[isample]->Name() << endl;
          h_tmp = 0;
          if(!doMTWtemplate) histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
          else histo_name = "MTW_" + channel_list[ichan] + "_" + datasets[isample]->Name();
          //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
          //  cout << "  --- histo " << histo_name << endl;
          if(!pseudodata_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
          h_tmp = (TH1F*) pseudodata_file->Get(histo_name.c_str());
          if(h_sum == 0) {h_sum = (TH1F*) h_tmp->Clone();}
          else {h_sum->Add(h_tmp);}
        }
        else{
          
          if(!doMTWtemplate) histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + template_fake_name;
          else if(doMTWtemplate) histo_name = "MTW_" + channel_list[ichan] + "_" + template_fake_name;
          // cout << "  --- histo " << histo_name << endl;
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
      string output_histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_data";
      if(doMTWtemplate) output_histo_name = "MTW_" + channel_list[ichan] + "_data";
      h_sum->SetTitle(output_histo_name.c_str());
      h_sum->SetName(output_histo_name.c_str());
      h_sum->Write(output_histo_name.c_str(), TObject::kOverwrite);
      
    } // chan
    
    pseudodata_file->Close();
    
    cout<<"--- Done with generation of pseudo-data"<<endl<<endl;
    
    delete pseudodata_file;
    
  }
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  if((makePlots || doPDFunc || PlotMVAvars || PlotSystematics) ){
    string pathOutput = "OutputPlots/";
    mkdir(pathOutput.c_str(),0777);
    string pathOutputdate = pathOutput + dateString + "/"  ;
    mkdir(pathOutputdate.c_str(),0777);
    string place =pathOutputdate+"MSPlot/";
    if(doMTWtemplate) place = pathOutputdate + "MSPlotMTW/";
    string placeTH1F = pathOutputdate+"TH1F/";
    string placeTH2F = pathOutputdate+"TH2F/";
    vector <string> vlabel_chan = {"3#mu", "1e2#mu", "2e1#mu", "3e"};
    mkdir(place.c_str(),0777);
    mkdir(placeTH1F.c_str(),0777);
    mkdir(placeTH2F.c_str(),0777);
    string rootFileName ="NtuplePlotsMVA.root";
    if(doMTWtemplate) rootFileName = "NtuplePlotsMTW.root";
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutputdate+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    
   if(makePlots && doMTWtemplate){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlotMTW.begin(); it != MSPlotMTW.end(); it++)
      {
        //cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        if(!datafound){
          cout << "no data found, setting lumi as " << Luminosity << endl;
          temp->setDataLumi(Luminosity);
        }
        //temp->SetPreliminary(false);
        temp->setDataLumi(Luminosity);
        if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
        if(name.find("eee")!=std::string::npos) temp->setChannel(true, "3e");
        if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "2e1#mu");
        if(name.find("uue")!=std::string::npos) temp->setChannel(true, "1e2#mu");
        if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "3#mu");
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        cout << "writing to " << pathOutputdate+"MSPlotMTW" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        temp->Write(fout, name, true, (pathOutputdate+"MSPlotMTW").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
    }
    
   if(makePlots && !doMTWtemplate){
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
        if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);
        //temp->SetPreliminary(false);
        temp->Draw(name, 1, false, false, false, 10);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        cout << "writing to " << pathOutputdate+"MSPlot" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        temp->Write(fout, name, true, (pathOutputdate+"MSPlot").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
    }
    if(doPDFunc && !doMTWtemplate){
      TDirectory* th1dir = fout->mkdir("1D_PDF_histograms");
      th1dir->cd();
      gStyle->SetOptStat(1110);
      for (std::map<std::string,TH1F*>::const_iterator it = histo1DPDF.begin(); it != histo1DPDF.end(); it++)
      {
        TH1F *temp = it->second;
        int N = temp->GetNbinsX();
        temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
        temp->SetBinContent(N+1,0);
        temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
        temp->Write();
        if(it->first.find("nominal")!=std::string::npos || it->first.find("Up")!=std::string::npos || it->first.find("Down")!=std::string::npos) {
          TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
          tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
        }
      }
    }
    if(PlotMVAvars && !doMTWtemplate){
      // cout << "plot mva vars " << endl;
      TDirectory* th1dirmva = fout->mkdir("1D_MVA_histograms");
      th1dirmva->cd();
      gStyle->SetOptStat(0);
      
      std::vector<string> channellist;
      channellist.push_back("all");
      channellist.push_back("eee");
      channellist.push_back("uue");
      channellist.push_back("eeu");
      channellist.push_back("uuu");
      // find all variables
      std::string splitname = "";
      char seperator = '_';
      std::vector < std::string > variables;
      for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
      {
        
        if(toppair && it->first.find("TT_FCNC")!=std::string::npos && it->first.find("all")!=std::string::npos){}
        else if(!toppair && it->first.find("ST_FCNC")!=std::string::npos && it->first.find("all")!=std::string::npos){
          //cout << it->first << endl;
          
        }
        else continue;
        
        splitname = (split(it->first, seperator))[0];
        
        variables.push_back(splitname);
      }
      for(int i = 0 ; i < variables.size(); i++){
        splitname = variables[i];
        // cout << "name " << splitname << endl;
        for(int iChan = 0; iChan < channellist.size(); iChan++){
          TH1F *tempBKG(0);
          TH1F *tempSignal(0);
          TH1F *tempfake(0);
          
          for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
          {
            if(it->first.find(channellist[iChan].c_str())==std::string::npos) continue;
            if(it->first.find(splitname.c_str())==std::string::npos) continue;
            
            //cout << "looking at " << it->first << endl;
            TH1F *temp = it->second;
            if(it->first.find("WZTo3LNu")!=std::string::npos || it->first.find("WZJTo3LNu")!=std::string::npos) {
              if(tempBKG == 0) tempBKG = (TH1F*) temp->Clone();
              else tempBKG->Add(temp);
            }
            if(it->first.find("fake")!=std::string::npos) {
              cout << "filling fake " << tempfake << " " << temp << endl;
              if(tempfake == 0) tempfake = (TH1F*) it->second->Clone();
              else tempfake->Add(temp);
              cout << "filling fake " << tempfake << endl;
            }
            if(it->first.find("T_FCNC")!=std::string::npos) {
              if(tempSignal == 0) tempSignal = (TH1F*) temp->Clone();
              else tempSignal->Add(temp);
            }
            //delete temp;
          }
          //cout << "filled histos " << endl;
          if(tempBKG == 0) cout << "ERROR tempBKG is null" << endl ;
          if(tempfake == 0) cout << "ERROR tempfake is null" << endl ;
          if(tempSignal == 0) cout << "ERROR tempSignal is null" << endl ;
          
          tempBKG->SetLineColor(kBlue);
          tempSignal->SetLineColor(kRed);
          // tempfake->SetLineColor(kGreen);
          tempBKG->SetName(splitname.c_str());
          tempBKG->SetTitle(("Shape comparison - " + channellist[iChan]).c_str());
          
          Double_t scaleBKG = 1./tempBKG->Integral();
          tempBKG->Scale(scaleBKG);
          Double_t scaleSig = 1./tempSignal->Integral();
          tempSignal->Scale(scaleSig);
          // Double_t scalefake = 1./tempfake->Integral();
          // tempfake->Scale(scalefake);
          double max = TMath::Max(tempSignal->GetMaximum(), tempBKG->GetMaximum());
          // double max = TMath::Max(max0, tempfake->GetMaximum());
          tempBKG->SetMaximum(max*1.2);
          tempBKG->GetXaxis()->SetTitle(splitname.c_str());
          tempBKG->GetYaxis()->SetTitle("Nb. Events");
          
          
          Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
          TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
          leg->AddEntry(tempSignal,"Signal","L");   // h1 and h2 are histogram pointers
          leg->AddEntry(tempBKG,"WZ background","L");
          //leg->AddEntry(tempfake,"DD non prompt","L");
          
          TCanvas* tempCanvas = TCanvasCreator(tempBKG,"Normalised MVA distribution" );
          tempBKG->Draw("h");
          tempSignal->Draw("h,Sames");
          //tempfake->Draw("L,Sames");
          leg->Draw("Same");
          tempCanvas->SaveAs( (placeTH1F+splitname+"_"+channellist[iChan]+".png").c_str() );
          
          tempBKG->Write();
          tempSignal->Write();
          // tempfake->Write();
          delete tempCanvas;
          delete tempSignal;
          //delete tempfake;
          delete tempBKG;
        } // channellist
      }
    }
    if(PlotSystematics && !doMTWtemplate){
      //cout << "plot systematics" << endl;
      TDirectory* th1dirsys = fout->mkdir("1D_Systematic_histograms");
      th1dirsys->cd();
      gStyle->SetOptStat(1110);
      string systematic = "";
      TH1F *temp_nom(0);
      for (std::map<std::string,TH1F*>::const_iterator it = histo1DSys.begin(); it != histo1DSys.end(); it++)
      {
        if(it->first.find("nominal")!=std::string::npos){
          if(temp_nom == 0) temp_nom = (TH1F*) it->second->Clone();
          else temp_nom->Add(it->second);
          //cout << "found nominal " << it->first << endl;
        }
      }
      
      for(int isys = 1; isys < thesystlistnames.size(); isys++){ // first value was nominal
        systematic = thesystlistnames[isys];
        //cout << "looking at " << systematic.c_str() << endl;
        TH1F *temp_up(0);
        TH1F *temp_down(0);
        TCanvas* Canvas = 0;
        string nameplot = systematic;
        for (std::map<std::string,TH1F*>::const_iterator it = histo1DSys.begin(); it != histo1DSys.end(); it++)
        {
          if(it->first.find(systematic.c_str())!=std::string::npos){
            if(it->first.find("Up")!=std::string::npos){
              if(temp_up ==0) temp_up = (TH1F*) it->second->Clone();
              else temp_up->Add(it->second);
            }//cout << "found " << it->first << endl; }
            else if(it->first.find("Down")!=std::string::npos){
              if(temp_down == 0) temp_down = (TH1F*) it->second->Clone();
              else temp_down->Add(it->second);
            }// cout << "found " << it->first << endl;}
            
            
          }
        }
        
        if(temp_up == 0 || temp_down == 0 || temp_up == 0){ cout << "Error someting went wrong with " << systematic.c_str() << " ! Skipping..." << endl; continue;}
        
        temp_down->Write();
        temp_nom->Write();
        temp_up->Write();
        temp_nom->SetLineColor(kRed);
        temp_up->SetLineColor(kBlue);
        temp_down->SetLineColor(kViolet);
        temp_up->SetTitle(("Influence of "+systematic).c_str());
        Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
        TLegend *legend = new TLegend(xl1,yl1,xl2,yl2);
        //TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);//(0.55,0.65,0.76,0.82);
        legend->AddEntry(temp_nom,"nominal","L");
        legend->AddEntry(temp_down,(systematic+" down").c_str(),"L");
        legend->AddEntry(temp_up,(systematic+" up").c_str(),"L");
        
        
        Canvas =  TCanvasCreator(temp_up, systematic.c_str() );//new TCanvas("Canvas_PU","Canvas_PU");
        Canvas->cd();
        gStyle->SetOptStat(0);
        temp_up->GetXaxis()->SetTitle("BDT");
        
        temp_up->Draw("h");
        temp_nom->Draw("SAME,h");
        temp_down->Draw("SAME,h");
        legend->Draw("SAME");
        Canvas->SaveAs( (placeTH1F+nameplot+".png").c_str() );
        Canvas->SetLogy();
        Canvas->Update();
        Canvas->SaveAs( (placeTH1F+nameplot+"_LogY.png").c_str() );
      }
    }
    if(PlotSystematics && doMTWtemplate){
      cout << "plot systematics" << endl;
      TDirectory* th1dirsys = fout->mkdir("1D_Systematic_histograms");
      th1dirsys->cd();
      gStyle->SetOptStat(1110);
      string systematic = "";
      TH1F *temp_nom(0);
      for (std::map<std::string,TH1F*>::const_iterator it = histo1DSysMTW.begin(); it != histo1DSysMTW.end(); it++)
      {
        cout << "looking at " << it->first<<endl;
        if(it->first.find("nominal")!=std::string::npos){
          cout << "found nominal " << it->first << " " << it->second << " " << temp_nom <<  endl;
          if(temp_nom == 0){ temp_nom = (TH1F*) (it->second); cout << "cloned" << endl;}
          else temp_nom->Add(it->second);
          cout << "added nominal " << it->first << endl;
        }
      }
      
      for(int isys = 1; isys < thesystlistnames.size(); isys++){ // first value was nominal
        systematic = thesystlistnames[isys];
        cout << "looking at " << systematic.c_str() << endl;
        TH1F *temp_up(0);
        TH1F *temp_down(0);
        TCanvas* Canvas = 0;
        string nameplot = systematic;
        for (std::map<std::string,TH1F*>::const_iterator it = histo1DSysMTW.begin(); it != histo1DSysMTW.end(); it++)
        {
          if(it->first.find(systematic.c_str())!=std::string::npos){
            if(it->first.find("Up")!=std::string::npos){
              if(temp_up ==0){ temp_up = (TH1F*) (it->second);}
              else temp_up->Add(it->second);
              cout << "found " << it->first << endl; }
            else if(it->first.find("Down")!=std::string::npos){
              if(temp_down == 0){ temp_down = (TH1F*) (it->second);}
              else temp_down->Add(it->second);
              cout << "found " << it->first << endl;}
            
            
          }
        }
        
        if(temp_up == 0 || temp_down == 0 || temp_up == 0){ cout << "Error someting went wrong with " << systematic.c_str() << " ! Skipping..." << endl; continue;}
        cout << "writing" << endl;
        temp_down->Write();
        temp_nom->Write();
        temp_up->Write();
        cout << "written" << endl;
        temp_nom->SetLineColor(kRed);
        temp_up->SetLineColor(kBlue);
        temp_down->SetLineColor(kViolet);
        temp_up->SetTitle(("Influence of "+systematic).c_str());
        Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
        TLegend *legend = new TLegend(xl1,yl1,xl2,yl2);
        //TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);//(0.55,0.65,0.76,0.82);
        legend->AddEntry(temp_nom,"nominal","L");
        legend->AddEntry(temp_down,(systematic+" down").c_str(),"L");
        legend->AddEntry(temp_up,(systematic+" up").c_str(),"L");
        
        
        Canvas =  TCanvasCreator(temp_up, systematic.c_str() );//new TCanvas("Canvas_PU","Canvas_PU");
        Canvas->cd();
        gStyle->SetOptStat(0);
        temp_up->GetXaxis()->SetTitle("M_T(W)");
        
        temp_up->Draw("h");
        temp_nom->Draw("SAME,h");
        temp_down->Draw("SAME,h");
        legend->Draw("SAME");
        Canvas->SaveAs( (placeTH1F+nameplot+".png").c_str() );
        Canvas->SetLogy();
        Canvas->Update();
        Canvas->SaveAs( (placeTH1F+nameplot+"_LogY.png").c_str() );
      }
    }
 // if(false){ // TO FIX
    if(makePlots && doMTWtemplate){
      hist_WZ->SetLineColor(kBlue);
      hist_fakes->SetLineColor(kGreen);
     hist_TT_FCNC->SetLineColor(kRed-2);

      
      
      hist_WZ->SetName("M_T(W)");
      hist_WZ->SetTitle("Shape comparison");
      
      Double_t scaleBKG_nom = 1./hist_WZ->Integral();
      hist_WZ->Scale(scaleBKG_nom);
      Double_t scaleSigTT_nom= 1./hist_TT_FCNC->Integral();
      hist_TT_FCNC->Scale(scaleSigTT_nom);
      Double_t scalefake_nom = 1./hist_fakes->Integral();
      cout << "scale fakes " << scalefake_nom << endl;
      hist_fakes->Scale(scalefake_nom);
      double max0 = TMath::Max(hist_fakes->GetMaximum(), hist_TT_FCNC->GetMaximum());
      double max1 = TMath::Max(hist_TT_FCNC->GetMaximum(), hist_WZ->GetMaximum());
      double maximum = TMath::Max(max0, max1);
      hist_WZ->SetMaximum(maximum*1.2);
      hist_WZ->GetXaxis()->SetTitle("M_T(W)");
      hist_WZ->GetYaxis()->SetTitle("Nb. Norm. Events");
      
      
      Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
      TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
      leg->AddEntry(hist_TT_FCNC,"Signal TT","L");
      leg->AddEntry(hist_WZ,"WZ background","L");
      leg->AddEntry(hist_fakes,"DD non prompt","L");
      
      TCanvas* tempCanvas = TCanvasCreator(hist_WZ,"Normalised distribution" );
      hist_WZ->Draw("h");
      hist_TT_FCNC->Draw("h,sames");
      hist_fakes->Draw("h,sames");
      leg->Draw("Same");
      tempCanvas->SaveAs( (placeTH1F+"MWT.png").c_str() );
      
      hist_WZ->Write();
      hist_TT_FCNC->Write();
      hist_fakes->Write();
      
      
      /*
      cout << "plotting mtW shapes " << endl;
      std::vector<string> channellist;
      channellist.push_back("all");
      channellist.push_back("eee");
      channellist.push_back("uue");
      channellist.push_back("eeu");
      channellist.push_back("uuu");

      for(int iChan = 0; iChan < channellist.size(); iChan++){
        TH1F *tempBKG_nom(0);
        TH1F *tempSignalST_nom(0);
        TH1F *tempfake_nom(0);
        TH1F *tempBKG_up(0);
        TH1F *tempSignalST_up(0);
        TH1F *tempBKG_down(0);
        TH1F *tempSignalST_down(0);
        TH1F *tempSignalTT_down(0);
        TH1F *tempSignalTT_up(0);
        TH1F *tempSignalTT_nom(0);
        
             //cout << "histo mtw size " << histo1DMTW.size() << endl;
       for (map<string,TH1F*>::const_iterator it = histo1DMTW.begin(); it != histo1DMTW.end(); it++)
        {
         // cout << "looking at " << it->first << " and the channel to keep " << channellist[iChan].c_str() << endl;
        
         if(it->first.find(channellist[iChan].c_str())==std::string::npos){
            //cout << "continuing " << endl;
            continue;
          }

          TH1F* temp = it->second;
          if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
            if(tempBKG_nom == 0) tempBKG_nom = (TH1F*) temp->Clone();
            else tempBKG_nom->Add(temp);
          }
          else if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
            if(tempBKG_up == 0) tempBKG_up = (TH1F*) temp->Clone();
            else tempBKG_up->Add(temp);
          }
          else if(it->first.find("WZTo3LNu")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
            if(tempBKG_down == 0) tempBKG_down= (TH1F*) temp->Clone();
            else tempBKG_down->Add(temp);
          }
          if(it->first.find("fake")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
            if(tempfake_nom == 0) tempfake_nom = (TH1F*) temp->Clone();
            else tempfake_nom->Add(temp);
          }
          if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
            cout << "st nom" << endl;
            if(tempSignalST_nom == 0) tempSignalST_nom = (TH1F*) temp->Clone();
            else tempSignalST_nom->Add(temp);
          }
          else if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
            if(tempSignalST_up == 0) tempSignalST_up = (TH1F*) temp->Clone();
            else tempSignalST_up->Add(temp);
          }
          else if(it->first.find("ST_FCNC")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
            if(tempSignalST_down == 0) tempSignalST_down= (TH1F*) temp->Clone();
            else tempSignalST_down->Add(temp);
          }
          if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("nominal")!= std::string::npos) {
             cout << "tt nom" << endl;
            if(tempSignalTT_nom == 0) tempSignalTT_nom = (TH1F*) temp->Clone();
            else tempSignalTT_nom->Add(temp);
          }
          else if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("Up")!= std::string::npos) {
            if(tempSignalTT_up == 0) tempSignalTT_up = (TH1F*) temp->Clone();
            else tempSignalTT_up->Add(temp);
          }
          else if(it->first.find("TT_FCNC")!=std::string::npos && it->first.find("Down")!= std::string::npos) {
            if(tempSignalTT_down == 0) tempSignalTT_down= (TH1F*) temp->Clone();
            else tempSignalTT_down->Add(temp);
          }
          delete temp;
        }
        //cout << "filled histos " << endl;
        if(tempBKG_nom == 0) cout << "ERROR tempBKG is null" << endl ;
        if(tempfake_nom == 0) cout << "ERROR tempfake is null" << endl ;
        if(tempSignalST_nom == 0) cout << "ERROR tempSignalST is null" << endl ;
        if(tempSignalTT_nom == 0) cout << "ERROR tempSignalTT is null" << endl ;
        
        tempBKG_nom->SetLineColor(kBlue);
        tempBKG_up->SetLineColor(kBlue);
        tempBKG_down->SetLineColor(kBlue);
        tempBKG_up->SetLineStyle(8);
        tempBKG_down->SetLineStyle(3);
        tempSignalST_nom->SetLineColor(kRed);
        tempSignalST_up->SetLineColor(kRed);
        tempSignalST_down->SetLineColor(kRed);
        tempSignalST_up->SetLineStyle(8);
        tempSignalST_down->SetLineStyle(3);
        tempfake_nom->SetLineColor(kGreen);
        tempSignalTT_nom->SetLineColor(kRed-2);
        tempSignalTT_up->SetLineColor(kRed-2);
        tempSignalTT_down->SetLineColor(kRed-2);
        tempSignalTT_up->SetLineStyle(8);
        tempSignalTT_down->SetLineStyle(3);
        
        
        tempBKG_nom->SetName("M_T(W)");
        tempBKG_nom->SetTitle(("Shape comparison - " + channellist[iChan]).c_str());
        
        Double_t scaleBKG_nom = 1./tempBKG_nom->Integral();
        tempBKG_nom->Scale(scaleBKG_nom);
        Double_t scaleBKG_up= 1./tempBKG_up->Integral();
        tempBKG_up->Scale(scaleBKG_up);
        Double_t scaleBKG_down = 1./tempBKG_down->Integral();
        tempBKG_down->Scale(scaleBKG_down);
        Double_t scaleSigST_nom = 1./tempSignalST_nom->Integral();
        tempSignalST_nom->Scale(scaleSigST_nom);
        Double_t scaleSigST_down = 1./tempSignalST_down->Integral();
        tempSignalST_down->Scale(scaleSigST_down);
        Double_t scaleSigST_up = 1./tempSignalST_up->Integral();
        tempSignalST_up->Scale(scaleSigST_up);
        Double_t scaleSigTT_nom = 1./tempSignalTT_nom->Integral();
        tempSignalST_nom->Scale(scaleSigST_nom);
        Double_t scaleSigTT_down = 1./tempSignalTT_down->Integral();
        tempSignalTT_down->Scale(scaleSigTT_down);
        Double_t scaleSigTT_up = 1./tempSignalTT_up->Integral();
        tempSignalTT_up->Scale(scaleSigTT_up);
        Double_t scalefake_nom = 1./tempfake_nom->Integral();
        tempfake_nom->Scale(scalefake_nom);
        double max0 = TMath::Max(tempSignalST_nom->GetMaximum(), tempBKG_nom->GetMaximum());
        double max1 = TMath::Max(tempSignalST_nom->GetMaximum(), tempfake_nom->GetMaximum());
        double max2 = TMath::Max(tempSignalST_nom->GetMaximum(), tempSignalTT_nom->GetMaximum());
        double max00 = TMath::Max(tempSignalST_up->GetMaximum(), tempBKG_up->GetMaximum());
        double max10 = TMath::Max(tempSignalST_up->GetMaximum(), tempfake_nom->GetMaximum());
        double max20 = TMath::Max(tempSignalST_up->GetMaximum(), tempSignalTT_up->GetMaximum());
        double max01 = TMath::Max(tempSignalST_down->GetMaximum(), tempBKG_down->GetMaximum());
        double max11 = TMath::Max(tempSignalST_down->GetMaximum(), tempfake_nom->GetMaximum());
        double max21 = TMath::Max(tempSignalST_down->GetMaximum(), tempSignalTT_down->GetMaximum());
        double maxA = TMath::Max(max0, max1);
        double maxB = TMath::Max(max2,max00);
        double maxC = TMath::Max(max10,max20);
        double maxD = TMath::Max(max01,max11);
        double maxX = TMath::Max(max21,maxA);
        double maxY = TMath::Max(maxB,maxC);
        double maxi = TMath::Max(maxD,maxX);
        double maximum = TMath::Max(maxi, maxY);
        tempBKG_nom->SetMaximum(maximum*1.2);
        tempBKG_nom->GetXaxis()->SetTitle("M_T(W)");
        tempBKG_nom->GetYaxis()->SetTitle("Nb. Events");
        
        
        Double_t xl1=0.7, yl1=.7, xl2=xl1+.2, yl2=yl1+.2;
        TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
        leg->AddEntry(tempSignalST_nom,"Signal ST","L");   // h1 and h2 are histogram pointers
        leg->AddEntry(tempSignalTT_nom,"Signal TT","L");
        leg->AddEntry(tempBKG_nom,"WZ background","L");
        leg->AddEntry(tempfake_nom,"DD non prompt","L");
        
        TCanvas* tempCanvas = TCanvasCreator(tempBKG_nom,"Normalised distribution" );
        tempBKG_nom->Draw("h");
        tempBKG_up->Draw("h,sames");
        tempBKG_down->Draw("h,sames");
        tempSignalST_nom->Draw("h,Sames");
        tempSignalST_up->Draw("h,Sames");
        tempSignalST_down->Draw("h,Sames");
        tempSignalTT_nom->Draw("h,Sames");
        tempSignalTT_up->Draw("h,Sames");
        tempSignalTT_down->Draw("h,Sames");
        tempfake_nom->Draw("L,Sames");
        leg->Draw("Same");
        tempCanvas->SaveAs( (placeTH1F+"MWT_"+channellist[iChan]+".png").c_str() );
        
        tempBKG_nom->Write();
        tempSignalST_nom->Write();
        tempSignalTT_nom->Write();
        tempfake_nom->Write();
        tempBKG_up->Write();
        tempSignalST_up->Write();
        tempSignalTT_up->Write();
        tempBKG_down->Write();
        tempSignalST_down->Write();
        tempSignalTT_down->Write();
       
        delete tempCanvas;
        delete tempSignalST_nom;
        delete tempSignalTT_nom;
        delete tempfake_nom;
        delete tempBKG_nom;
        delete tempBKG_up;
        delete tempSignalST_up;
        delete tempSignalTT_up;
        delete tempBKG_down;
        delete tempSignalST_down;
        delete tempSignalTT_down;
        
      } // channellist
      */
    }
    fout->Write();
    fout->Close();
    
    delete fout;
    
  }
  
 /* delete hist_WZ;
  delete hist_TT_FCNC;
  delete hist_fakes;*/
  ///************************************///
  ///   ADD PDF UNC TO COMBINE TEMPLATE  ///
  ///************************************///
  if(doPDFunc && doMTWtemplate){
    TFile* combinetemplate_file = TFile::Open( combinetemplate_filename.c_str(), "UPDATE" );
    combinetemplate_file->cd();
    
    TH1::SetDefaultSumw2();
    TH1F *hist_uuu     = new TH1F( (coupling + "_" + region+"_uuu").c_str(),           (coupling + "_" + region+"_uuu").c_str(),           nbin, -1, 1 );
    TH1F *hist_uue     = new TH1F( (coupling + "_" + region+"_uue").c_str(),           (coupling + "_" + region+"_uue").c_str(),           nbin, -1, 1 );
    TH1F *hist_eeu     = new TH1F( (coupling + "_" + region+"_eeu").c_str(),           (coupling + "_" + region+"_eeu").c_str(),           nbin, -1, 1 );
    TH1F *hist_eee     = new TH1F( (coupling + "_" + region+"_eee").c_str(),           (coupling + "_" + region+"_eee").c_str(),           nbin, -1, 1 );
    
    //NB : theta name convention = <observable>__<process>[__<uncertainty>__(plus,minus)] FIX ME
    output_histo_name = "";
    vector<string> v_sys = {"PDFup", "PDFdown"};
    dataSetName  = "WZTo3LNu_3Jets_MLL50_80X";
    for(int isys = 0; isys < v_sys.size() ; isys++)
    {
      systematic = v_sys[isys];
      output_histo_name = coupling + "_" + region+"_uuu_"  + dataSetName + "_" + systematic ;
      hist_uuu->SetTitle(output_histo_name.c_str());
      hist_uuu->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_uue_"  + dataSetName + "_" + systematic ;
      hist_uue->SetTitle(output_histo_name.c_str());
      hist_uue->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_eeu_"  + dataSetName + "_" + systematic ;
      hist_eeu->SetTitle(output_histo_name.c_str());
      hist_eeu->Write(output_histo_name.c_str());
      output_histo_name = coupling + "_" + region+"_eee_"  + dataSetName + "_" + systematic ;
      hist_eee->SetTitle(output_histo_name.c_str());
      hist_eee->Write(output_histo_name.c_str());
    }
    combinetemplate_file->Close();
    delete combinetemplate_file;
    //delete  hist_eee; delete hist_uuu; delete hist_uue; delete hist_eeu;
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
}



///// BOOK KEEPING FUNCTIONS
std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    if (end != start) {
      tokens.push_back(text.substr(start, end - start));
    }
    start = end + 1;
  }
  if (end != start) {
    tokens.push_back(text.substr(start));
  }
  return tokens;
}
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
string intToStr (int number){
  ostringstream buff;
  buff<<number;
  return buff.str();
}
double maximumValue(vector<double> array){
  int length = array.size();  // establish size of array
  double max = array[0];       // start with max = first element
  
  for(int i = 1; i<length; i++)
  {
    if(array[i] > max)
      max = array[i];
  }
  return max;                // return highest value in array
}
double minimumValue(vector<double> array){
  int length = array.size();  // establish size of array
  double max = array[0];       // start with max = first element
  
  for(int i = 1; i<length; i++)
  {
    if(array[i] < max)
      max = array[i];
  }
  return max;                // return highest value in array
}

//// INITIALISATIONS
void InitMSPlotsMTW(string prefix, vector <int> decayChannels){
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    decaystring += prefix;
    
    //cout << "init msplots " << endl;
    MSPlotMTW[("MTW_"+decaystring).c_str()] = new MultiSamplePlot(datasets, ("MTW_"+decaystring).c_str(), nbinMTW, 0,endMTW, "transv. mass W boson (GeV)","GeV");
  }
  decaystring = "";
}
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair, bool isZut){
  clock_t start_sub = clock();
  
  prefix = prefix + "_";
  
  // control plots
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    //cout << "init " << (prefix+"_BDT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), nbin, -1.,1., "BDT");
    MSPlot[ (prefix+"channel_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"channel_"+decaystring).c_str(), 5,-0.5, 4.5, "decay");
    MSPlot[ (prefix+"weight_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"weight_"+decaystring).c_str(), 100,0, 0.3, "eventweight");
    
    if(!istoppair){
      MSPlot[(prefix+"mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"mlb_"+decaystring).c_str(),10, 0, 500, "inv. mass l_{W}b (GeV)","GeV");
      MSPlot[(prefix+"dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRWlepb_"+decaystring).c_str(),10,0, 5, "#Delta R(l_{W},b)");
      MSPlot[(prefix+"dPhiWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dPhiWlepb_"+decaystring).c_str(),10,-4, 4, "#Delta #Phi (l_{W},b)");
      MSPlot[(prefix+"Zboson_pt_"+decaystring).c_str()]= new MultiSamplePlot(datasets, (prefix+"Zboson_pt_"+decaystring).c_str(), 20,0, 500, "Z boson p_{T} (GeV)", "GeV");
      MSPlot[(prefix+"dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZWlep_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,l_{W})");
      MSPlot[(prefix+"bdiscCSVv2_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"bdiscCSVv2_jet_0_"+decaystring).c_str(),20, 0.5, 1, "leading jet CSVv2");
      MSPlot[(prefix+"cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_0_"+decaystring).c_str(),10, 0, 0.8, "leading jet charm vs b disc.");
      
      if(isZut){
        MSPlot[(prefix+"charge_asym_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"charge_asym_"+decaystring).c_str(),10, -4, 4, "charge l_{W} x|W boson #eta|");
        
      }
      else{
        MSPlot[(prefix+"cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_0_"+decaystring).c_str(),25, 0, 1, "leading jet charm vs light disc.");
        
      }
      
    }
    else if(istoppair ){
      MSPlot[(prefix+"mlb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"mlb_"+decaystring).c_str(),25, 0, 500, "inv mass l_{W}b (GeV)","GeV");
      MSPlot[(prefix+"FCNCtop_M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"FCNCtop_M_"+decaystring).c_str(),20, 100, 500, "inv. mass Zq (GeV)","GeV");
      MSPlot[(prefix+"dRWlepb_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRWlepb_"+decaystring).c_str(),20,0, 6, "#Delta R(b,l_{W})");
      MSPlot[(prefix+"nJets_CSVv2M_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"nJets_CSVv2M_"+decaystring).c_str(),10,-0.5, 9.5, "# CSVv2M jets");
      MSPlot[(prefix+"nJets_CharmL_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"nJets_CharmL_"+decaystring).c_str(),10,-0.5, 9.5, "# charm loose jets");
      MSPlot[(prefix+"dRZWlep_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZWlep_"+decaystring).c_str(),20,0, 6, "#Delta R(Z,l_{W})");
      MSPlot[(prefix+"dRZc_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"dRZc_"+decaystring).c_str(),30,0, 6, "#Delta R(Z,q)");
      
      if(isZut){
        MSPlot[(prefix+"cdiscCvsB_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_0_"+decaystring).c_str(),20, 0, 0.8, "leading jet vharm vs b disc.");
        MSPlot[(prefix+"cdiscCvsB_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsB_jet_1_"+decaystring).c_str(),20, 0, 0.85, "2nd leading charm vs b disc. ");
        
      }
      else{
        MSPlot[(prefix+"cdiscCvsL_jet_0_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_0_"+decaystring).c_str(),20, 0, 1, "leading jet charm vs light disc.");
        MSPlot[(prefix+"cdiscCvsL_jet_1_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"cdiscCvsL_jet_1_"+decaystring).c_str(),20, 0, 1, "2nd leading jet charm vs light disc.");
        
        
      }
    }
    
  }
  decaystring = "";
  
  
  
  
}
void InitCalculatePDFWeightHisto(string dataSetName){
  TH1::SetDefaultSumw2();
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  for(int iChan = 0; iChan < channel_list.size(); iChan++){
    channel = channel_list[iChan];
    for ( int i=0; i<101; i++)
    {
      output_histo_name = dataSetName+"_BDT_"+channel+"_"+intToStr(i);
      histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    }
    output_histo_name = dataSetName+"_BDT_"+channel+"_nominal";
    histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    output_histo_name = dataSetName+"_BDT_"+channel+"_PDFEnvelopeUp";
    histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    output_histo_name = dataSetName+"_BDT_"+channel+"_PDFEnvelopeDown";
    histo1DPDF[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    
    output_histo_name = "";
  }
}
void InitMTWShapeHisto(string dataSetName, string systematic, int isys,  vector <int> decayChannels){
  TH1::SetDefaultSumw2();
  //histo1DMTW.clear();
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    

    if(isys == 0) output_histo_name = dataSetName+"_MTW_nominal_"+decaystring;
    else output_histo_name = dataSetName+"_MTW_"+systematic + "_" + decaystring;
   //cout << "init " << output_histo_name.c_str() << endl;
    histo1DMTW[output_histo_name.c_str()] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbinMTW,0.,endMTW);
    
    output_histo_name = "";
  }
  decaystring = "";
}
void InitSystematicHisto(string dataSetName, string systematic, int isys, bool doMTWtemplate){
  TH1::SetDefaultSumw2();
  
  if(isys == 0 && !doMTWtemplate) output_histo_name = dataSetName+"_BDT_nominal";
  else if(isys == 0 && doMTWtemplate) output_histo_name = dataSetName+"_MTW_nominal";
  else if(!doMTWtemplate) output_histo_name = dataSetName+"_BDT_"+systematic;
  else if(doMTWtemplate) output_histo_name = dataSetName+"_MTW_"+systematic;
  
  if(!doMTWtemplate) histo1DSys[output_histo_name] = new TH1F(output_histo_name.c_str(), dataSetName.c_str(), nbin,-1.,1.);
  else if(doMTWtemplate) histo1DSysMTW[output_histo_name] = new TH1F(output_histo_name.c_str(), dataSetName.c_str(), nbinMTW,0.,endMTW);
  output_histo_name = "";
  
}
void InitAnalyzerTree(TTree* tree){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  tree->SetBranchAddress("MVA_weight_nom", &MVA_weight_nom, &b_MVA_weight_nom);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress("MVA_EqLumi", &MVA_EqLumi, &b_MVA_EqLumi);
  tree->SetBranchAddress("MVA_Luminosity", &MVA_Luminosity, &b_MVA_Luminosity);
  tree->SetBranchAddress("MVA_mWt2", &MVA_mWt2, &b_MVA_mWt2);
    tree->SetBranchAddress("MVA_mWt", &MVA_mWt, &b_MVA_mWt);
  
  tree->SetBranchAddress("MVA_x1", &MVA_x1, &b_MVA_x1);
  tree->SetBranchAddress("MVA_x2", &MVA_x2, &b_MVA_x2);
  tree->SetBranchAddress("MVA_id1", &MVA_id1, &b_MVA_id1);
  tree->SetBranchAddress("MVA_id2", &MVA_id2, &b_MVA_id2);
  tree->SetBranchAddress("MVA_q", &MVA_q, &b_MVA_q);
  
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
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  if(istoppair && doZut){
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_1", &MVA_cdiscCvsB_jet_1, &b_MVA_cdiscCvsB_jet_1);
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
    tree->SetBranchAddress("MVA_dRZc", &MVA_dRZc, &b_MVA_dRZc);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_FCNCtop_M", &MVA_FCNCtop_M, &b_MVA_FCNCtop_M);
    tree->SetBranchAddress("MVA_nJets_CharmL", &MVA_nJets_CharmL, &b_MVA_nJets_CharmL);
    tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
  }
  else  if(istoppair && !doZut){
    //tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
    tree->SetBranchAddress("MVA_cdiscCvsL_jet_1", &MVA_cdiscCvsL_jet_1, &b_MVA_cdiscCvsL_jet_1);
    tree->SetBranchAddress("MVA_cdiscCvsL_jet_0", &MVA_cdiscCvsL_jet_0, &b_MVA_cdiscCvsL_jet_0);
    tree->SetBranchAddress("MVA_dRZc", &MVA_dRZc, &b_MVA_dRZc);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_FCNCtop_M", &MVA_FCNCtop_M, &b_MVA_FCNCtop_M);
    tree->SetBranchAddress("MVA_nJets_CharmL", &MVA_nJets_CharmL, &b_MVA_nJets_CharmL);
    tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
  }
  else if(!istoppair && doZut){
    
    tree->SetBranchAddress("MVA_Zboson_pt", &MVA_Zboson_pt, &b_MVA_Zboson_pt);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    
  }
  else if(!istoppair && !doZut){
    
    tree->SetBranchAddress("MVA_Zboson_pt", &MVA_Zboson_pt, &b_MVA_Zboson_pt);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_cdiscCvsL_jet_0", &MVA_cdiscCvsL_jet_0, &b_MVA_cdiscCvsL_jet_0);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
  }
  
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  tree->SetBranchAddress("MVA_weight_nom", &MVA_weight_nom, &b_MVA_weight_nom);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress("MVA_BDT", &MVA_BDT, &b_MVA_BDT);
  tree->SetBranchAddress("MVA_EqLumi", &MVA_EqLumi, &b_MVA_EqLumi);
  tree->SetBranchAddress("MVA_Luminosity", &MVA_Luminosity, &b_MVA_Luminosity);
  
  tree->SetBranchAddress("MVA_x1", &MVA_x1, &b_MVA_x1);
  tree->SetBranchAddress("MVA_x2", &MVA_x2, &b_MVA_x2);
  tree->SetBranchAddress("MVA_id1", &MVA_id1, &b_MVA_id1);
  tree->SetBranchAddress("MVA_id2", &MVA_id2, &b_MVA_id2);
  tree->SetBranchAddress("MVA_q", &MVA_q, &b_MVA_q);
  
  
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
void Init1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels){
  TH1::SetDefaultSumw2();
  cout << "initialising MVA var histo" << endl;
  // control plots
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    output_histo_name = "BDT_"+dataSetName +"_"+decaystring+"_"+systematic;
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    
    if(!istoppair){
      output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,0,500);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,0,5);
      output_histo_name = "dPhiWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,-4,4);
      output_histo_name = "ZbosonPt_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
      output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),30,0,6);
      output_histo_name = "bdiscCSVv2jet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0.5,1);
      output_histo_name = "cdiscCvsBjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
      if(isZut){
        output_histo_name = "chargeAsym_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,-4,4);
        
      }
      else
      {
        output_histo_name = "cdiscCvsLjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),25,0,1);
        
      }
      
      
    }
    else if(istoppair ){
      output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),25,0,500);
      output_histo_name = "FCNCtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,6);
      output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,6);
      output_histo_name = "dRZc_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),30,0,6);
      output_histo_name = "nJetsCSVv2M_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,-0.5,9.5);
      output_histo_name = "nJetsCharmL_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),10,-0.5,9.5);
      
      if(isZut){
        output_histo_name = "cdiscCvsBjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
        output_histo_name = "cdiscCvsBjet1_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
        
      }
      else{
        output_histo_name = "cdiscCvsLjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
        output_histo_name = "cdiscCvsLjet1_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(),20,0,500);
      }
    }
    
  }
  
  output_histo_name = "";
  decaystring = "";
}
////////// FUNCTIONS
vector<double> BDTCUT(string region, string coupling){
  
  cout << "Determine BDT cut " << endl;
  
  
  string bdtinput_name = placeOutputReading+"/Reader_" + coupling + "_" + region + ".root";
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
        histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        cout << "  --- histo " << histo_name << endl;
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
        histo_name = coupling + "_BDT_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        if(h_tmp->GetEntries() != 0){
          if(h_sum_sig == 0) {h_sum_sig = (TH1F*) h_tmp->Clone();}
          else {h_sum_sig->Add(h_tmp);}
        }
      }
      else{
        h_tmp = 0;
        histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + template_fake_name;
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
void CalculatePDFWeight(string dataSetName, double BDT, double MVA_weight_nom, int MVA_channel){
  // cout << "calculate pdf" << endl;
  //std::vector<double> pdfweights;
  //cout << "MVA channel " << MVA_channel << endl;
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  /*
   //PDF weights calculation
   LHAPDF::setVerbosity(0);
   string sbase = "";
   if(dataSetName.find("WZTo3LNu")!=std::string::npos) sbase = "NNPDF30_nlo_as_0118"; //kevin (ttbar)
   //cout << "base set " << sbase << endl;
   LHAPDF::PDFSet basepdfSet(sbase.c_str());
   //LHAPDF::PDFSet basepdfSet("NNPDF30_lo_as_0130"); // base from main MC // WZ
   //LHAPDF::PDFSet basepdfSet("NNPDF23_lo_as_0130_qed");
   LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_100"); // give the correct name see https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#PDF_uncertainties
   //LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_30");
   const LHAPDF::PDF* basepdf = basepdfSet.mkPDF(0);
   channel = channel_list[MVA_channel];
   // cout << "MVA_x1 "<< MVA_x1 <<" MVA_x2 "<< MVA_x2 <<" MVA_id1 "<< MVA_id1 <<" MVA_id2 "<< MVA_id2 <<" MVA_q "<< MVA_q << endl;
   
   
   for ( size_t i=0; i<newpdfSet.size(); i++)
   {
   const LHAPDF::PDF* newpdf = newpdfSet.mkPDF(i);
   
   double weightpdf = LHAPDF::weightxxQ(MVA_id1, MVA_id2, MVA_x1, MVA_x2, MVA_q, *basepdf, *newpdf);
   output_histo_name = dataSetName+"_BDT_"+channel+"_"+intToStr(i);
   histo1DPDF[output_histo_name]->Fill(MVA_BDT, MVA_weight_nom*weightpdf);
   // cout << "fill " << (dataSetName+"_BDT"+"_"+intToStr(i)).c_str() << " with " << BDT << " " <<  weightpdf << endl;
   //cout << "pdf weight " << weight << endl;
   //pdfweights.push_back(weight);
   
   delete newpdf;
   
   }
   output_histo_name = dataSetName+"_BDT_"+channel+"_nominal";
   histo1DPDF[output_histo_name]->Fill(MVA_BDT, MVA_weight_nom);
   output_histo_name = "";
   delete basepdf;
   */
  
}
void FillMTWPlots(int d, string postfix, vector <int> decayChannels, double weight_, int MVA_channel){
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  

  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    decaystring += postfix;
    
    //if(datasets[d]->Name().find("data")!=std::string::npos) cout << "filling " << ("MTW_"+decaystring).c_str() << " with " << MVA_mWt2 << " " << weight_ << endl;
    MSPlotMTW[("MTW_"+decaystring).c_str()]->Fill(MVA_mWt2 , datasets[d], true, weight_);
  }
  decaystring = "";
}
void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isZut , bool istoppair, double weight_, int MVA_channel){
  
  //cout << "fill plots" << endl;
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  
  
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    
    //cout << "bdt " << MVA_BDT << " in " << (prefix+"_BDT_"+decaystring).c_str()<< endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()]->Fill(MVA_BDT , datasets[d], true, weight_);
    MSPlot[ (prefix+"_channel_"+decaystring).c_str()]->Fill(MVA_channel, datasets[d], true, weight_);
    MSPlot[ (prefix+"_weight_"+decaystring).c_str()]->Fill(weight_, datasets[d], true, 1.);
    
    if(!istoppair){
      MSPlot[(prefix+"_mlb_"+decaystring).c_str()] ->Fill(MVA_mlb, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRWlepb_"+decaystring).c_str()] ->Fill(MVA_dRWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_dPhiWlepb_"+decaystring).c_str()] ->Fill(MVA_dPhiWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_Zboson_pt_"+decaystring).c_str()]->Fill(MVA_Zboson_pt, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZWlep_"+decaystring).c_str()] ->Fill(MVA_dRZWlep, datasets[d], true, weight_);
      MSPlot[(prefix+"_bdiscCSVv2_jet_0_"+decaystring).c_str()]->Fill(MVA_bdiscCSVv2_jet_0, datasets[d], true, weight_);
      MSPlot[(prefix+"_cdiscCvsB_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsB_jet_0, datasets[d], true, weight_);
      
      if(isZut){
        MSPlot[(prefix+"_charge_asym_"+decaystring).c_str()]->Fill(MVA_charge_asym, datasets[d], true, weight_);
        
      }
      else{
        MSPlot[(prefix+"_cdiscCvsL_jet_0_"+decaystring).c_str()]->Fill(MVA_cdiscCvsL_jet_0, datasets[d], true, weight_);
        
      }
      
    }
    else if(istoppair ){
      MSPlot[(prefix+"_mlb_"+decaystring).c_str()] ->Fill(MVA_mlb, datasets[d], true, weight_);
      MSPlot[(prefix+"_FCNCtop_M_"+decaystring).c_str()] ->Fill(MVA_FCNCtop_M, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRWlepb_"+decaystring).c_str()] ->Fill(MVA_dRWlepb, datasets[d], true, weight_);
      MSPlot[(prefix+"_nJets_CSVv2M_"+decaystring).c_str()]->Fill(MVA_NJets_CSVv2M,  datasets[d], true, weight_);
      MSPlot[(prefix+"_nJets_CharmL_"+decaystring).c_str()] ->Fill(MVA_nJets_CharmL, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZWlep_"+decaystring).c_str()] ->Fill(MVA_dRZWlep, datasets[d], true, weight_);
      MSPlot[(prefix+"_dRZc_"+decaystring).c_str()] ->Fill(MVA_dRZc, datasets[d], true, weight_);
      
      if(isZut){
        MSPlot[(prefix+"_cdiscCvsB_jet_0_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsB_jet_0, datasets[d], true, weight_);
        MSPlot[(prefix+"_cdiscCvsB_jet_1_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsB_jet_1, datasets[d], true, weight_);
        
      }
      else{
        MSPlot[(prefix+"_cdiscCvsL_jet_0_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsL_jet_0, datasets[d], true, weight_);
        MSPlot[(prefix+"_cdiscCvsL_jet_1_"+decaystring).c_str()] ->Fill(MVA_cdiscCvsL_jet_1, datasets[d], true, weight_);
        
        
      }
    }
  }
  decaystring = "";
}

void GetPDFEnvelope(string dataSetName){
  double binContentMax = -1.;
  double binContentMin = 1000000000000000.;
  vector<double> bincontents;
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  // loop over channels
  for(int iChan = 0; iChan < channel_list.size(); iChan++){
    channel = channel_list[iChan];
    output_histo_name = dataSetName+"_BDT_" + channel + "_nominal";
    // get nominal th1F
    TH1F* histo_nom = (TH1F*) histo1DPDF[output_histo_name]->Clone();
    
    // loop over bins
    for( int ibin = 1; ibin <histo_nom->GetNbinsX(); ibin++)
    {
      binContentMax = -1.;
      binContentMin = 1000000000000000.;
      bincontents.clear();
      // get bincontents of each histo for this bin
      for(int iCount = 0; iCount < 101; iCount++)
      {
        output_histo_name = dataSetName+"_BDT_"+channel + "_" +intToStr(iCount);
        bincontents.push_back(histo1DPDF[output_histo_name]->GetBinContent(ibin));
      }
      if(binContentMin > minimumValue(bincontents)) binContentMin = minimumValue(bincontents);
      else binContentMin = histo_nom->GetBinContent(ibin);
      if(binContentMax < maximumValue(bincontents)) binContentMax = maximumValue(bincontents);
      else binContentMax = histo_nom->GetBinContent(ibin);
      
      output_histo_name = dataSetName+"_BDT_"+channel + "_PDFEnvelopeUp";
      histo1DPDF[output_histo_name]->SetBinContent(ibin, binContentMax);
      output_histo_name = dataSetName+"_BDT_"+channel + "_PDFEnvelopeDown";
      histo1DPDF[output_histo_name]->SetBinContent(ibin, binContentMin);;
      output_histo_name = "";
    }// bins
  }//channels
  
}
void Fill1DHisto(string dataSetName, string systematic, bool istoppair, bool isZut, vector <int> decayChannels, double weight_, int MVA_channel){
  
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    output_histo_name = "BDT_"+dataSetName +"_"+decaystring+"_"+systematic;
    histo1D[output_histo_name] ->Fill(  MVA_BDT        ,weight_);
    
    if(!istoppair){
      output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill( MVA_mlb         ,weight_);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_dRWlepb       ,weight_);
      output_histo_name = "dPhiWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(  MVA_dPhiWlepb        ,weight_);
      output_histo_name = "ZbosonPt_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_Zboson_pt       ,weight_);
      output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_dRZWlep        ,weight_);
      output_histo_name = "bdiscCSVv2jet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(    MVA_bdiscCSVv2_jet_0      ,weight_);
      output_histo_name = "cdiscCvsBjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_cdiscCvsB_jet_0       ,weight_);
      if(isZut){
        output_histo_name = "chargeAsym_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill( MVA_charge_asym         ,weight_);
        
      }
      else
      {
        output_histo_name = "cdiscCvsLjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(   MVA_cdiscCvsL_jet_0       ,weight_);
        
      }
      
      
    }
    else if(istoppair ){
      output_histo_name = "mlb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_mlb       ,weight_);
      output_histo_name = "FCNCtopM_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_FCNCtop_M       ,weight_);
      output_histo_name = "dRWlepb_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(    MVA_dRWlepb      ,weight_);
      output_histo_name = "dRZWlep_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_dRZWlep       ,weight_);
      output_histo_name = "dRZc_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_dRZc       ,weight_);
      output_histo_name = "nJetsCSVv2M_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill(   MVA_NJets_CSVv2M       ,weight_);
      output_histo_name = "nJetsCharmL_"+dataSetName + "_" +decaystring+"_"+systematic;
      histo1D[output_histo_name] ->Fill( MVA_nJets_CharmL        ,weight_);
      
      if(isZut){
        output_histo_name = "cdiscCvsBjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(  MVA_cdiscCvsB_jet_0        ,weight_);
        output_histo_name = "cdiscCvsBjet1_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(   MVA_cdiscCvsB_jet_1       ,weight_);
        
      }
      else{
        output_histo_name = "cdiscCvsLjet0_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(   MVA_cdiscCvsL_jet_0       ,weight_);
        output_histo_name = "cdiscCvsLjet1_"+dataSetName + "_" +decaystring+"_"+systematic;
        histo1D[output_histo_name] ->Fill(   MVA_cdiscCvsL_jet_1       ,weight_);
      }
    }
    
  }
  
  output_histo_name = "";
  decaystring = "";
  
}
void FillMTWShapeHisto(string dataSetName, string systematic, double weight_,int isys, int MVA_channel, vector <int> decayChannels){
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";

    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    if(isys == 0) output_histo_name = dataSetName+"_MTW_nominal_"+decaystring;
    else output_histo_name = dataSetName+"_MTW_"+systematic + "_" + decaystring;
    
   // cout << "fill " << output_histo_name << " " << MVA_mWt2 << " " << weight_ <<  endl;
   histo1DMTW[output_histo_name.c_str()]->Fill(MVA_mWt2, 1.);
    // histo1DMTW[output_histo_name.c_str()]->Fill(1, 1);
    output_histo_name = "";
  }
  decaystring = "";
}

void FillSystematicHisto(string dataSetName, string systematic, double weight_, int isys, bool doMTWtemplate ){

  if(isys == 0 && !doMTWtemplate) output_histo_name = dataSetName+"_BDT_nominal";
  else if(isys == 0 && doMTWtemplate) output_histo_name = dataSetName+"_MTW_nominal";
  else if(!doMTWtemplate) output_histo_name = dataSetName+"_BDT_"+systematic;
  else if(doMTWtemplate) output_histo_name = dataSetName+"_MTW_"+systematic;
  
  
  if(!doMTWtemplate) histo1DSys[output_histo_name]->Fill(MVA_BDT, weight_);
  else histo1DSysMTW[output_histo_name]->Fill(MVA_mWt2, weight_);
  
  
  output_histo_name = "";
  
}

