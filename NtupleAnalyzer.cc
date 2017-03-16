#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <sstream>
#include <TFile.h>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"




using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TFile*> globalFileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,TNtuple*> ntree;
map<string,TNtuple*> otree;
map<string,TTree*> globalttree;
map<string,MultiSamplePlot*> MSPlot;


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG);
void MSPCreator (string pathPNG);
void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();





// CONFIGURATION
Bool_t debug = false;
bool mumumu  = false;
bool eee = false;
double DataLumi = -1;
bool elecPlot = false;
bool muPlot = false;
//applying all appropriate scale factors for individual objects if the bool is set to true
Bool_t applyElectronSF = false;
Bool_t applyMuonSF = false;
Bool_t applyPUSF = false;
Bool_t applyGlobalSF = false;
Bool_t applyAMC = false;
Bool_t applyBtagSF = false;
Bool_t NewPhys = false;
Bool_t applyMET = false;
Bool_t makeTree = false;
string dateString;
int dChan = -1;
string decayChan = "all";
/*if(selectedMuons.size() == 3) {channelInt = 0; i_channel = 0;}
else if(selectedElectrons.size() == 3) {channelInt = 3; i_channel = 3;}
else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {channelInt = 2; i_channel = 2; }
else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){channelInt = 1; i_channel = 1; }*/


int main(int argc, char* argv[])
{
  if (debug){
    cout << "argc = " << argc << endl;
    for(int i = 0; i < argc; i++)
    {
      cout << "argv[" << i << "] = " << argv[i] << endl;
    }
  }

  
  
  //Placing arguments in properly typed variables
  string channel = "eee";
/*  debug = false;
  applyElectronSF = false;
  applyMuonSF = false;
  applyPUSF = false;
  applyGlobalSF = false;
  debug = strtol(argv[2],NULL,10);
  applyElectronSF = strtol(argv[3],NULL,10);
  applyMuonSF = strtol(argv[4],NULL,10);
  applyPUSF = strtol(argv[5],NULL,10);
  applyBtagSF = strtol(argv[6],NULL,10);
  applyGlobalSF = strtol(argv[7],NULL,10);
  applyAMC = strtol(argv[8],NULL,10);
  applyMET = strtol(argv[9],NULL,10);
  string decayChan = argv[10];
  */
  
  
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("channel")!=std::string::npos) {
      channel = argv[i+1];
      i++;
    }
    if(string(argv[i]).find("debug")!=std::string::npos) {
      debug =true;
    }
    if(string(argv[i]).find("applyPU")!=std::string::npos) {
      applyPUSF =true;
    }
    if(string(argv[i]).find("applyGlobalSF")!=string::npos) {
       applyGlobalSF=true;
    }
    if(string(argv[i]).find("applyElectronSF")!=string::npos) {
      applyElectronSF =true;
    }
    if(string(argv[i]).find("applyMuonSF")!=string::npos) {
      applyMuonSF =true;
    }
    if(string(argv[i]).find("applyMET")!=string::npos) {
      applyMET =true;
    }
    if(string(argv[i]).find("applyBtagSF")!=string::npos) {
      applyBtagSF =true;
    }
    if(string(argv[i]).find("applyAMC")!=string::npos) {
      applyAMC =true;
    }
    if(string(argv[i]).find("decay")!=string::npos) {
      decayChan = argv[i+1];
      i++;
    }
    if(string(argv[i]).find("MakeTree")!=string::npos) {
      makeTree = true;
    }
    
  }
  
  
  
  if(decayChan.find("eee")!=string::npos) dChan = 3;
  if(decayChan.find("eeu")!=string::npos) dChan = 2;
  if(decayChan.find("uue")!=string::npos) dChan = 1;
  if(decayChan.find("uuu")!=string::npos) dChan = 0;
  if(decayChan.find("all")!=string::npos) dChan = 9;
  
  
  
  string xmlFileName;
  string CraneenPath;
  cout << " --> Using the all channel..." << endl;
  xmlFileName = "config/Run2TriLepton_samples_analy.xml";
  mumumu = false;
  eee = false;
  DataLumi  =   36417.624743959; //  21016; // 36417624743.959/1000000 ;// ub * 1000000 =   pb-1
  cout << "datalumi " << DataLumi << endl;
  
  dateString = MakeTimeStamp();
  
  // where the tuples are
  CraneenPath = "NtupleMakerOutput/MergedTuples/";
  CraneenPath += "singletop/170211/";
  
  
  
  
  
  // where the plots go
  string pathPNG = "myOutput";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/" + dateString + "/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MSPlots/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;
  
  
  // info
  cout << "xmlFileName is " << xmlFileName << endl;
  cout << "NtupleFolder is " << CraneenPath << endl;
  cout << "decaychannel is " << decayChan << endl;
  
  if(debug) cout << "debugging on" << endl;
  if(applyElectronSF) cout << "Electron SF on " << endl;
  if(applyMuonSF) cout << "Muon SF on " << endl;
  if(applyPUSF) cout << "PU SF on" <<endl;
  if(applyAMC) cout << "nlo reweight on " << endl;
  if(applyBtagSF) cout << "BtagSF on" << endl;
  if(applyMET) cout << "MET filter on "<< endl;
  if(!applyElectronSF && !applyMuonSF && !applyPUSF && !applyAMC) applyGlobalSF = false;
  if(applyGlobalSF) cout << "applying SF " << endl;
  else cout << "not applying SF " << endl;
  
  
  // calling datasetPlotter to create MSPplots
  //  DatasetPlotter(25, -1, 1, "BDTscore", xmlFileName,CraneenPath,pathPNG);
  
  
  // event plots
  // DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath,pathPNG);
  //
 /*
  if(dChan != 0){
    elecPlot = true;
    muPlot = false;
    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath,pathPNG);
    DatasetPlotter(50, 0, 500, "pt_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
    
    DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
    DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
    //DatasetPlotter(100, -0.1, 0.1, "d0_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
    //DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
    DatasetPlotter(100, 0, 1000, "E_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
  }
  
  if(dChan != 3){
    // muon plots
    muPlot = true;
    elecPlot = false;
    DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath,pathPNG);
    DatasetPlotter(50, 0, 500, "pt_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
    
    DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
    DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
    //DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
    //DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
    //DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
  }
  elecPlot = false;
  muPlot = false;
  DatasetPlotter(11, -0.5, 10.5, "nJets", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(70, 0, 700, "pt_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50, -3.15, 3.15, "eta_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(30, -3.15, 3.15, "phi_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(25, 0, 1, "bdisc_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  // DatasetPlotter(25,-1, 1, "cdiscCvsL_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  //  DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(11, -0.5, 10.5, "nJets_CSVL", xmlFileName,CraneenPath,pathPNG);
  
  
  DatasetPlotter(11, -0.5, 10.5, "nJets_CSVM", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(11, -0.5, 10.5, "nJets_CSVT", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(11, -0.5, 10.5, "nJets_nonCSVL", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(11, -0.5, 10.5, "nJets_nonCSVM", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(11, -0.5, 10.5, "nJets_nonCSVT", xmlFileName,CraneenPath,pathPNG);
*/  /*
   
   // MVA plots
   DatasetPlotter(25,-1, 1, "CosTheta", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(25,-1, 1, "CosTheta_alt", xmlFileName,CraneenPath,pathPNG);
   
   DatasetPlotter(25,-6, 6, "SMtop_Eta", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(40,80, 400, "SMtop_M", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(25,-3.15, 3.15, "SMtop_Phi", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50,0, 500, "SMtop_Pt", xmlFileName,CraneenPath,pathPNG);
   
   DatasetPlotter(100,0, 2000, "TotalHt", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(175,0, 35000, "TotalInvMass", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50,0, 500, "TotalPt", xmlFileName,CraneenPath,pathPNG);
  
  DatasetPlotter(4,-2, 2, "Wlep_Charge", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50,-2.5, 2.5, "Wlep_Eta", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(25,-3.15, 3.15, "Wlep_Phi", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50,0, 500, "Wlep_Pt", xmlFileName,CraneenPath,pathPNG);
  
  DatasetPlotter(30, 60, 120, "Zboson_M", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50, 0, 1000, "Zboson_Energy", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50,-6, 6, "Zboson_Eta", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(25,-3.15, 3.15, "Zboson_Phi", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(50,0, 500, "Zboson_Pt", xmlFileName,CraneenPath,pathPNG);
  
  DatasetPlotter(25,0, 1, "bdiscCSVv2_jet_1", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(25,0, 1, "bdiscCSVv2_jet_2", xmlFileName,CraneenPath,pathPNG);
  
 // DatasetPlotter(25,-1, 1, "cdiscCvsL_jet_1", xmlFileName,CraneenPath,pathPNG);
 // DatasetPlotter(25,-1, 1, "cdiscCvsB_jet_1", xmlFileName,CraneenPath,pathPNG);
 // DatasetPlotter(25,-1, 1, "cdiscCvsL_jet_2", xmlFileName,CraneenPath,pathPNG);
 // DatasetPlotter(25,-1, 1, "cdiscCvsB_jet_2", xmlFileName,CraneenPath,pathPNG);
  
   DatasetPlotter(50,-2.5, 2.5, "charge_asym", xmlFileName,CraneenPath,pathPNG);
  
  DatasetPlotter(20,-3.15, 3.15, "dPhiWlepb", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(20,-3.15, 3.15, "dPhiZMET", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(20,-3.15, 3.15, "dPhiZSMtop", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(20,-3.15, 3.15, "dPhiZWlep", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(20,-3.15, 3.15, "dPhiZb", xmlFileName,CraneenPath,pathPNG);
  
  DatasetPlotter(10,0, 6, "dRWlepb", xmlFileName,CraneenPath,pathPNG);
 
  DatasetPlotter(10,0, 6, "dRZSMtop", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(10,0, 6, "dRZWlep", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(10,0,6, "dRZb", xmlFileName,CraneenPath,pathPNG);
  
   DatasetPlotter(20, 0, 400, "mWt", xmlFileName,CraneenPath,pathPNG);
  
   DatasetPlotter(70, 0, 700, "met_Pz", xmlFileName,CraneenPath,pathPNG);
  DatasetPlotter(70, 0, 700, "met_Pt", xmlFileName,CraneenPath,pathPNG);

   DatasetPlotter(200, 0, 200, "mlb", xmlFileName,CraneenPath,pathPNG);
 */
 //

  
  MSPCreator (pathPNG);
  
}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string PathPNG)
{
  cout<<""<<endl;
  cout<<"RUNNING NOMINAL DATASETS: "<< sVarofinterest <<endl;
  cout<<""<<endl;
  
  const char *xmlfile = xmlNom.c_str();
  if(debug) cout << "used config file: " << xmlfile << endl;
  
  
  
  string histo_dir = "NtupleMakerOutput/MyAnalysis/";
  string histo_dirdecay = histo_dir +"singletop";
  string histo_dir_date = histo_dirdecay+"/ControlHisto_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
  mkdir(histo_dirdecay.c_str(),0777);
  mkdir(histo_dir_date.c_str(),0777);
  
  
  
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  if (debug) cout << "will start loading from xml file ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  if (debug) cout << "finished loading from xml file ..." << endl;
  
  
  //***************************************************CREATING PLOTS****************************************************
  //  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  //  outfile->cd();
  string plotname = sVarofinterest;   ///// Non Jet Split plot
  // make for loop here!!!
  MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
  MSPlot[plotname.c_str()]->setChannel(true, decayChan); 
  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath , slumi;
  
  int nEntries;
  float ScaleFactor, NormFactor;
  Int_t  varofInterest;
  float varofInterest_double[20];
  
 
  
  
  vector<string> v;
  // to avoid modifying original string
  // first duplicate the original string and return a char pointer then free the memory
  if(debug) cout << "LOOKING at " << sVarofinterest.c_str() << endl;
  char delim[] = " []";
  char * dup = strdup(sVarofinterest.c_str());
  char * token = strtok(dup, delim);
  while(token != NULL){
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
  
  //  if (debug) cout << v[0] << "  " << v[1] << endl;
  double weightv2 = 0. ;
  double weightv3 = 0.;
  double Xsect = 0.;
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    dataSetName = datasets[d]->Name();
  //  if(dataSetName.find("zut_80X")==string::npos && dataSetName.find("tZq_80X")==string::npos && dataSetName.find("Zjets50_80X")==string::npos ) continue;
    
        Xsect = datasets[d]->Xsection();
    cout<<"Dataset:  :"<<dataSetName<< " with xsec: " << Xsect << endl;
    
    string rootFileName (histo_dir_date+"/FCNC_3L_"+dataSetName+".root");
    cout << "Histofile: " << rootFileName << endl;
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
    
    
    histo1D["muonSF"]                                  = new TH1F("muonSF", "Muon ScaleFactors", 100, 0, 1);
    histo1D["electronSF"]                                  = new TH1F("electronSF", "Electron ScaleFactors", 100, 0, 1);
    histo1D["btagSF"]  = new TH1F("btagSF", "Btag ScaleFactors", 100, 0, 1);
    histo1D["puSF"]  = new TH1F("puSF", "PUScaleFactors", 100, 0, 1);
    
    
    
    NewPhys = false;
    if(dataSetName.find("NP")!=string::npos ) NewPhys = true;
    // get the tree corresponding to the final state of interest
    string stree = "tree";
    string sglobaltree = "globaltree";
    string sbaselinetree = "baselinetree";
    filepath = TreePath + dataSetName + ".root";
/
		  
    if(makeTree) FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"update"); //create TFile for each dataset
    else FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ");
//    string TTreename = sbaselinetree;
    string TTreename = stree;
    
    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
    nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
    cout<<"                 nEntries: "<<nEntries<<endl;
    
    
    globalFileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
    string globalTTreename = sglobaltree;
    
    if(debug) cout << "globalTTreename: " << globalTTreename.c_str() << endl;
    if(debug) cout << "globalFileObj " << globalFileObj[dataSetName.c_str()] << endl;
    globalttree[dataSetName.c_str()] = (TTree*)globalFileObj[dataSetName.c_str()]->Get(globalTTreename.c_str()); //get ttre for each dataset
    if(debug) cout << "globalttree " << globalttree[dataSetName.c_str()]<< endl;
    int globalnEntries = (int)globalttree[dataSetName.c_str()]->GetEntries();
    cout<<"                 nEntries gt: "<<globalnEntries<<endl;
    
    
    
    // bo logic to set the right branch address depending on the string given as argument of the datasetplotter  (int or double[n] )
    if (v.size() == 2){
      ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&varofInterest);
      ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),varofInterest_double);
    }
    
    else if (v.size() == 1){
      if (debug)	cout << "v.size is to 1" << " and v[0] is " << v[0] << endl ;
      ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);//&varofInterest // faco To be fixed!
      
    }
    else {
      cout << "Vector of string does not have the good size!!!" << endl;
    }
    // eo logic to set the right branch address depending on the string given as argument of the datasetplotter
    
    bool isData= false;
    bool isAMC = false;
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
    if(debug) cout << "isData? " << isData << endl;
    if(dataSetName.find("amc")!=string::npos) isAMC =true;
    cout << "                 isAMC? " << isAMC << endl;
    ///////////////////////////////////
    // determine event scalefactor ///
    //////////////////////////////////
    
    if(applyGlobalSF) cout << "                 Applying scale factors (not for data)" << endl;
    
    // get the SF from the corresponding branch
    
    // new branch
    Float_t Weight = 1;
    TBranch *bWeight= ttree[dataSetName.c_str()]->Branch("Weight",&Weight, "Weight/F");
    
    
    
    Int_t DecayChannel = 0;
    ttree[dataSetName.c_str()]->SetBranchAddress("channelInt",&DecayChannel);
    
    
    Int_t PassedMET = 0;
    ttree[dataSetName.c_str()]->SetBranchAddress("PassedMETFilter",&PassedMET);
    
    Double_t puSF = 1. ;
    ttree[dataSetName.c_str()]->SetBranchAddress("puSF",&puSF);
    
    Double_t nloW;
    ttree[dataSetName.c_str()]->SetBranchAddress("nloWeight",&nloW);
    
    Double_t electronSF[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("ElectronSF",&electronSF);
    
    Double_t muonID[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonIDSF", &muonID);
    
    Double_t muonIso[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonIsoSF", &muonIso);
    
    
    Int_t nEl;
    ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons",&nEl);
    
    Int_t nMu;
    ttree[dataSetName.c_str()]->SetBranchAddress("nMuons",&nMu);
    
    Int_t nPosW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofPosWeights",&nPosW);
    
    Int_t nNegW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofNegWeights",&nNegW);
    
    Int_t nEvents;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nEv",&nEvents);
    
    Int_t SumW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("sumW",&SumW);
    
    Int_t nbHLTv2;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofEventsHLTv2",&nbHLTv2);
    
    Int_t nbHLTv3;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofEventsHLTv3", &nbHLTv3);

    Double_t BSF;
    ttree[dataSetName.c_str()]->SetBranchAddress("btagSF",&BSF);
    

    
    
    
    
    Int_t NbCuts;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nCuts", &NbCuts);
    
    double  CutSteps[10];
    //      globalttree[dataSetName.c_str()]->SetBranchAddress("cutstep[nCuts]", &CutSteps);
    
    if(debug) cout << "done setting SF addresses " << endl;
    
    // -----------
    // eo of event SF
    
    
    
    double globalScaleFactor= 1.;
    double muonSF = 1.;
    double BtagSFtemp = 1.;
    double puSFtemp= 1.;
    double electronSFtemp = 1.;
    double EquilumiSF = 1.;
    double nloSF = 1.;
    int nPos = 0;
    int nNeg = 0;
    int Ev = 0;
    int Weights = 0;
    
    if(!isData){
      int eventCounter = 0;
      for (int k = 0; k<globalnEntries; k++)
      {
        globalttree[(dataSetName).c_str()]->GetEntry(k);
        eventCounter += nEvents;
      
      }
      EquilumiSF = eventCounter / Xsect;
      cout << "                 equilumi = " <<  eventCounter <<" / " << Xsect <<" = " << EquilumiSF << endl;
    }
    
    if(!isData)
    {
      
      for (int k = 0; k<globalnEntries; k++)
      {
        globalttree[(dataSetName).c_str()]->GetEntry(k);
        if(debug) cout << "get globaltree" << endl;
        if(applyAMC && isAMC && !isData) nPos += nPosW;
        if(applyAMC && isAMC && !isData) nNeg += nNegW;
        if(applyAMC && isAMC && !isData) Ev += nEvents;
        if(applyAMC && isAMC && !isData) Weights += SumW;
        //        if(dataSetName.find("NP_overlay_FCNC_TT")!= string::npos)  Matched += nMatch;
        //        if(dataSetName.find("NP_overlay_FCNC_TT")!= string::npos) nonMatched += nNonMatch;
        // cout << "nPos " << nPos << " vs " << nPosW << " nNeg " << nNeg << " vs " << nNegW << " + " << nPos + nNeg << " - " << nPos - nNeg  << endl;
        // cout << "nEvents " << nEvents << " vs " << Ev << " sumWeights " << SumW << " vs " << Weights << endl;
        if(NewPhys){
          NmatchCharm += nMatched_c;
          NonmatchCharm += nNonMatched_c;
          NmatchBottom += nMatched_b;
          NonmatchBottom += nNonMatched_b;
          NmatchZelec += nMatched_Zel;
          NonmatchZelec += nNonMatched_Zel;
          NmatchZmu += nMatched_Zm;
          NonmatchZmu += nNonMatched_Zm;
          NmatchWelec += nMatched_Wel;
          NonmatchWelec += nNonMatched_Wel;
          NmatchWmu += nMatched_Wm;
          NonmatchWmu += nNonMatched_Wm;
          TagEqM += nMatched_tm;
          TagNeqM += nNonMatched_tm;
          NmatchTcharm += nMatched_ct;
          NonmatchTcharm += nNonMatched_ct;
        }
      }
      //          if(!isData) nloSF *= (double) Weights/(double) Ev; //
      if(applyAMC && isAMC && !isData) nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
      if(applyAMC && isAMC && !isData) cout << "                 nloSF: " << nloSF << endl;
      
    }
    for (int j = 0; j<nEntries; j++)
    {
      ttree[(dataSetName).c_str()]->GetEntry(j);
      //          cout << "nEl " << nEl << " nMu " << nMu << endl;
      globalScaleFactor = 1.;
      muonSF = 1.;
      electronSFtemp = 1.;
      puSFtemp = 1.; 
      BtagSFtemp = 1.;
      
      if(v.size() == 1 && sVarofinterest.find("nElectrons")!=string::npos) {varofInterest = nEl;}
      if(v.size() == 1 && sVarofinterest.find("nMuons")!=string::npos) {varofInterest = nMu;}
      
      if(applyMET && PassedMET == 0 ){continue; }
      if( dChan != 9 && DecayChannel != dChan){
        //cout << dChan << " - " << DecayChannel<< endl;
        continue;}
      
      if(applyGlobalSF && !isData) // sf on and not data
      {
        // Electron scale factors
        if(applyElectronSF)
        {
          for(unsigned int i = 0; i < nEl ; i ++)
          {
            // if(debug) cout << "lepton sf at index " << i << " is " << electronSF[i] << endl;
            globalScaleFactor *= electronSF[i];
            //if(debug) cout << "the globalScaleFactor is " << globalScaleFactor << endl;
            electronSFtemp *= electronSF[i];
          }
           histo1D["electronSF"]->Fill(electronSFtemp);
        }
        if(applyMuonSF)
        {
          for(unsigned int i = 0; i < nMu ; i ++)
          {
            //	   if(debug) cout << "Muon ID sf at index " << i << " is " << muonID[i] << endl;
            //	   if(debug) cout << "Muon Iso sf at index " << i << " is " << muonIso[i] << endl;
            //	   if(debug) cout << "Muon trig v2 sf at index " << i << " is " << muonTrigv2[i] << endl;
            //	   if(debug) cout << "Muon trig v3 sf at index " << i << " is " << muonTrigv3[i] << endl;
            //		   if(isData) weightv2 = (double) nbHLTv2 / (double) (nbHLTv2 + nbHLTv3);
            //           if(isData) weightv3 = (double) nbHLTv3 / (double) (nbHLTv2 + nbHLTv3);
            //		   cout << "weightv2 " << weightv2 << " weightv3 " << weightv3 << endl;
            globalScaleFactor *= muonID[i] *  muonIso[i]  ;
            muonSF *= muonID[i] *  muonIso[i];
            //	   if(debug) cout << "the globalScaleFactor is " << globalScaleFactor << endl;
          }
          histo1D["muonSF"]->Fill(muonSF);
        }
        if(applyPUSF)
        {
          globalScaleFactor *= puSF;
          puSFtemp *= puSF;
          histo1D["puSF"]->Fill(puSFtemp);
          if (debug){
            //      	cout << "puSF is " << puSF << endl;
            //	        cout << "the globalScaleFactor is " << globalScaleFactor << endl;
          }
          
        }
        if(applyBtagSF)
        {
          globalScaleFactor *= BSF;
          BtagSFtemp *= BSF;
          histo1D["btagSF"]->Fill(BtagSFtemp);
          
        }
        
        
        
      }
      
      if(applyAMC && !isData) globalScaleFactor =globalScaleFactor * nloSF * nloW;
      if(NewPhys) globalScaleFactor = 1.;
      
      // ----------------
      // eo event SF
      // make MS plot for single value
      if (v.size() == 1){
       /* if (isData)
        {
          // for data, fill once per event, weighted with the event scale factor only ???? what??
          MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, 1);
        }
        else
        {
          // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file (set to 1 everywhere)
          MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*DataLumi);
        }*/
        if(isData) globalScaleFactor = 1.;
        MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*DataLumi/EquilumiSF);
      }
      // make MS plot for vector
      if (v.size() == 2){
        
        // bo of loop over the number of object per entry
        //if(elecPlot) varofInterest = nEl;
       // if(muPlot) varofInterest = nMu;
        for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
        {
          if (debug) cout << "varofInterest is " << varofInterest_double[i_object] << endl;
         /* if (isData)
          {
            // for data, fill once per event, weighted with the event scale factor
            MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], false,1);
          }
          else
          {
            // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi = 1 is gotten from the xml file
            MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, globalScaleFactor*DataLumi);
            
          }
          */
          if(isData) globalScaleFactor = 1.;
          MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, globalScaleFactor*DataLumi/EquilumiSF);
        }
        
      }
      Weight = globalScaleFactor*EquilumiSF;
      if(makeTree) bWeight->Fill();
      
    } // nentries
    cout<<"                 event SF: "<<globalScaleFactor << " / " << EquilumiSF << " = " << globalScaleFactor/EquilumiSF <<endl;
    TCanvas *canv = new TCanvas(("canv_"+v[0]+dataSetName).c_str(),("canv_"+v[0]+dataSetName).c_str());
    string writename = "";
    if(isData)
    {
      writename = "data_nominal";
    }
    else
    {
      writename = dataSetName +"_nominal";
    }
    if(makeTree){
      
      FileObj[dataSetName.c_str()]->cd();
      ttree[dataSetName.c_str()]->Write();
      FileObj[dataSetName.c_str()]->Close();
    }
    
    
    
    
    
    fout-> cd();
    for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
      cout << "1D Plot: " << it->first << endl;
      TCanvas *ctemp = new TCanvas();
      ctemp->cd();
      TH1F *temp = it->second;
      temp->Draw();
      delete ctemp;
    }
    for (map<string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {
      cout << "2D Plot: " << it->first << endl;
      TCanvas *ctemp = new TCanvas();
      ctemp->cd();
      TH2F *temp = it->second;
      temp->Draw();
      delete ctemp;
    }
    fout->Write();
    fout->Close();
    delete fout;
    
  } // datasets
  
  
  
  if (debug){
    cout << "before cleaning" << endl;
    if (v.size() == 2){
      cout << " v[0] is " << v[0] << " and v[1] is " << v[1] << endl;
    }
    else if (v.size() == 1){
      cout << " v[0] is " << v[0] << endl;
    }
  }
  
  
  // clearing vector
  v.clear();
  if (debug){
    cout << "after cleaning" << endl ;
    cout << "v.size() is " << v.size() << endl;
  }
}; // datasetplotter


// function that writes all the MSPlots created in a root file
void MSPCreator (string pathPNG)
{
  // cout << pathPNG.c_str() << endl;
  
  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  outfile->cd();
  cout << "created " << (pathPNG+"/Output.root").c_str() << endl;
  
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    
    string name = "MyMSP_" + it->first;
    cout << " name " << name << endl;
    MultiSamplePlot *temp = it->second;
    if (debug){
      cout << "Saving the MSP" << endl;
      cout << " and it->first is " << it->first << endl;
    }
    // void MultiSamplePlot::Draw(string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal)
    temp->Draw("MyMSP"+it->first, 1, false, false, false, 1);// 0 = no ratio 1 = ratio
    name += "_"+decayChan;
    if(!applyGlobalSF) name += "_noSF";
    if(!applyPUSF) name += "_noPUSF";
    if(!applyElectronSF) name += "_noElSF";
    if(!applyMuonSF) name+= "_noMuSF";
    if(!applyAMC) name+= "_noAMCcor";
    if(!applyBtagSF) name+= "_noBtagSF";
    cout << "name " << name << endl;
    temp->Write(outfile, name, true,pathPNG.c_str() , "png");
    //      vector<string> temp_histo = it->GetTH1FNames();
    //      for (int i_hist=0; i_hist < temp_histo.size();i_hist++  ){
    //	cout << "hist is" << temp_histo[i_hist] << endl;
    //	cout << "integral is " << it->GetTH1F(temp_histo[i_hist].GetSum()) << endl;
    //      }
  }
  
  outfile->Write("kOverwrite");
}



string ConvertIntToString(int Number, bool pad)
{
  ostringstream convert;
  convert.clear();
  if ( pad && Number < 10 ) { convert << std::setw(2) << std::setfill('0');}
  convert << Number;
  return convert.str();
};


string MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, true);
  string month_str = ConvertIntToString(month, true);
  string day_str = ConvertIntToString(day, true);
  string hour_str = ConvertIntToString(hour, true);
  string min_str = ConvertIntToString(min, true);
  //string sec_str = ConvertIntToString(sec, true);
  
  string date_str = year_str + month_str + day_str; //+ "_" + hour_str + min_str;
  return date_str;
};
