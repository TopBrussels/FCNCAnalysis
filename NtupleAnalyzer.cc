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


//inlcludes for TMVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"


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
void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel);
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();





// CONFIGURATION
Bool_t debug = true;
bool mumumu  = false;
bool eee = false;
string channelpostfix = "";
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
Bool_t train_mva = false;
string dateString;

int main(int argc, char* argv[])
{
  if (debug){
    cout << "argc = " << argc << endl;
    for(int i = 0; i < argc; i++)
    {
      cout << "argv[" << i << "] = " << argv[i] << endl;
    }
  }
  if(argc < 7) cout << " ERROR: 6 arguments expected" << endl;
  
  
  //Placing arguments in properly typed variables
  const string channel = argv[1];
  debug = false;
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
  
  
  
  string xmlFileName;
  string CraneenPath;
  CraneenPath = "NtupleMakerOutput/MergedTuples/";
  if(channel=="MuMuMu")
  {
    cout << " --> Using the TriMuon channel..." << endl;
    channelpostfix = "_mumumu";
    xmlFileName = "config/Run2TriLepton_samples_analyzer_mumumu.xml";
    mumumu = true;
    eee = false;
    DataLumi  = 2100 ;//2612.180735004;//  pb-1
    CraneenPath += "mumumu/";
  }
  else if(channel=="ElElEl")
  {
    cout << " --> Using the TriElectron channel..." << endl;
    channelpostfix = "_eee";
    xmlFileName = "config/Run2TriLepton_samples_analyzer_eee.xml";
    mumumu = false;
    eee = true;
    DataLumi  = 2612.180735004;//  pb-1
    CraneenPath += "eee/";
  }
  else if(channel=="All")
  {
    cout << " --> Using the all channel..." << endl;
    channelpostfix = "_all";
    xmlFileName = "config/Run2TriLepton_samples_analy.xml";
    mumumu = false;
    eee = false;
    DataLumi  = 2612.180735004;//  pb-1
    CraneenPath += "all/";
  }
  else
  {
    cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
    exit(1);
  }
  dateString = MakeTimeStamp();
  //    CraneenPath += dateString + "/";
  CraneenPath += "160602_trial/";
  string pathPNG = "myOutput";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/" + dateString + "/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MSPlots"+channelpostfix+"/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;
  
  cout << "xmlFileName is " << xmlFileName << endl;
  cout << "NtupleFolder is " << CraneenPath << endl;
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
  
  // event plots
  //    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath,pathPNG);
  //    DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath,pathPNG);
  //    DatasetPlotter(10, -0.5, 9.5, "nLeptons", xmlFileName,CraneenPath,pathPNG);
  // DatasetPlotter(20, 0, 200, "Zboson_M", xmlFileName,CraneenPath,pathPNG);
  //DatasetPlotter(20, 0, 400, "mWt", xmlFileName,CraneenPath,pathPNG);
  //  DatasetPlotter(30, 0, 300, "met_Pt", xmlFileName,CraneenPath,pathPNG);
  // DatasetPlotter(20,-1, 1, "cdiscCvsL_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  // DatasetPlotter(20,-1, 1, "cdiscCvsB_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
  /*
   
   
   elecPlot = true;
   muPlot = false;
   DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50, 0, 500, "pt_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50, -3.15, 3.15, "eta_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(30, -3.15, 3.15, "phi_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, -0.1, 0.1, "d0_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, 0.0, 0.2, "pfIso_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, 0, 1000, "E_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG);
   
   
   // muon plots
   muPlot = true;
   elecPlot = false;
   DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50, 0, 500, "pt_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50, -3.15, 3.15, "eta_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(30, -3.15, 3.15, "phi_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(100, 0.0, 0.2, "pfIso_muon[nMuons]", xmlFileName,CraneenPath,pathPNG);
   
   elecPlot = false;
   muPlot = false;
   DatasetPlotter(11, -0.5, 10.5, "nJets", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(70, 0, 700, "pt_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50, -3.15, 3.15, "eta_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(30, -3.15, 3.15, "phi_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(25, 0, 1, "bdisc_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(25,-1, 1, "cdiscCvsL_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[nJets]", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(11, -0.5, 10.5, "nJets_CSVL", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(11, -0.5, 10.5, "nJets_CSVM", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(11, -0.5, 10.5, "nJets_CSVT", xmlFileName,CraneenPath,pathPNG);
   
   
   DatasetPlotter(70, 0, 700, "met_Pt", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50,0, 500, "pt_electron_1", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(20, 0, 200, "Zboson_M", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(20, 0, 400, "mWt", xmlFileName,CraneenPath,pathPNG);
   DatasetPlotter(50,0, 500, "pt_electron_2", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_electron_3", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_muon_1", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_muon_2", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_muon_3", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_jet_1", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_jet_2", xmlFileName, CraneenPath, pathPNG);
   DatasetPlotter(50,0, 500, "pt_jet_3", xmlFileName, CraneenPath, pathPNG);
   //    DatasetPlotter(10,-0.5, 9.5, "cutstep[nCuts], xmlFileName,, CraneenPath, pathPNG);
   */
  // calling the function that writtes all the MSPlots in a root file
  //void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel)
  MVAanalysis(train_mva, "BDT", 2, "tZq_amc", xmlFileName, xmlFileName, CraneenPath, channel); // divide sample in 2
  if(!train_mva) MSPCreator (pathPNG);
  
}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string PathPNG)
{
  cout<<""<<endl;
  cout<<"RUNNING NOMINAL DATASETS: "<< sVarofinterest <<endl;
  cout<<""<<endl;
  
  const char *xmlfile = xmlNom.c_str();
  if(debug) cout << "used config file: " << xmlfile << endl;
  
  
  
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
  
  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath , slumi;
  
  int nEntries;
  float ScaleFactor, NormFactor;
  int  varofInterest;
  double varofInterest_double [20];
  
  
  
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
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    dataSetName = datasets[d]->Name();
    cout<<"Dataset:  :"<<dataSetName<<endl;
    NewPhys = false;
    if(dataSetName.find("NP")!=string::npos) NewPhys = true;
    // get the tree corresponding to the final state of interest
    string stree = "tree";
    string sglobaltree = "globaltree";
    string sbaselinetree = "baselinetree";
    filepath = TreePath + dataSetName + ".root";
    
		  
    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
    string TTreename = stree;
    //     string TTreename = stree;
    
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
    cout << "isAMC? " << isAMC << endl;
    ///////////////////////////////////
    // determine event scalefactor ///
    //////////////////////////////////
    
    if(applyGlobalSF) cout << "                 Applying scale factors (not for data)" << endl;
    
    // get the SF from the corresponding branch
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
    
    Double_t muonTrigv2[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonTrigSFv2", &muonTrigv2);
    
    Double_t muonTrigv3[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonTrigSFv3", &muonTrigv3);
    
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
    double nloSF = 1.;
    int nPos = 0;
    int nNeg = 0;
    int Ev = 0;
    int Weights = 0;
    if(applyAMC && isAMC && !isData)
    {
      
      for (int k = 0; k<globalnEntries; k++)
      {
        globalttree[(dataSetName).c_str()]->GetEntry(k);
        nPos += nPosW;
        nNeg += nNegW;
        Ev += nEvents;
        Weights += SumW;
        // cout << "nPos " << nPos << " vs " << nPosW << " nNeg " << nNeg << " vs " << nNegW << " + " << nPos + nNeg << " - " << nPos - nNeg  << endl;
        // cout << "nEvents " << nEvents << " vs " << Ev << " sumWeights " << SumW << " vs " << Weights << endl;
      }
      //          if(!isData) nloSF *= (double) Weights/(double) Ev; //
      if(!isData) nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
      cout << "                nloSF: " << nloSF << endl;
    }
    for (int j = 0; j<nEntries; j++)
    {
      ttree[(dataSetName).c_str()]->GetEntry(j);
      //          cout << "nEl " << nEl << " nMu " << nMu << endl;
      globalScaleFactor = 1.;
      if(v.size() == 1 && sVarofinterest.find("nElectrons")!=string::npos) {varofInterest = nEl;}
      if(v.size() == 1 && sVarofinterest.find("nMuons")!=string::npos) {varofInterest = nMu;}
      
      if(applyMET && PassedMET == 0 ){continue; }
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
          }
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
            //	   if(debug) cout << "the globalScaleFactor is " << globalScaleFactor << endl;
          }
        }
        if(applyPUSF)
        {
          globalScaleFactor *= puSF;
          if (debug){
            //      	cout << "puSF is " << puSF << endl;
            //	        cout << "the globalScaleFactor is " << globalScaleFactor << endl;
          }
          
        }
        if(applyBtagSF)
        {
          globalScaleFactor *= BSF;
        }
        
        
        
      }
      
      if(applyAMC && !isData) globalScaleFactor =globalScaleFactor * nloSF * nloW;
      if(NewPhys) globalScaleFactor = 1.;
      
      // ----------------
      // eo event SF
      // make MS plot for single value
      if (v.size() == 1){
        if (isData)
        {
          // for data, fill once per event, weighted with the event scale factor only ???? what??
          MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, 1);
        }
        else
        {
          // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
          MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*DataLumi);
        }
      }
      // make MS plot for vector
      if (v.size() == 2){
        
        // bo of loop over the number of object per entry
        if(elecPlot) varofInterest = nEl;
        if(muPlot) varofInterest = nMu;
        for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
        {
          if (debug) cout << "varofInterest is " << varofInterest_double[i_object] << endl;
          if (isData)
          {
            // for data, fill once per event, weighted with the event scale factor
            MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], false,1);
          }
          else
          {
            // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
            MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, globalScaleFactor*DataLumi);
            
          }
          
        }
        
      }
      
    } // nentries
    cout<<"                 event SF: "<<globalScaleFactor <<endl;
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
    temp->Draw("MyMSP"+it->first, 1, false, false, false, 10);// 0 = no ratio
    //    name += "_3L";
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


void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel)
{
  
  MVAComputer* Eventcomputer_ =0;
  MVATrainer* Eventtrainer_ = 0;
  if(doTraining) Eventtrainer_ = new MVATrainer(MVAmethod,"TrainedEventMVA"+channel, "TrainedEventMVA"+channel+".root");
  vector<std::string> MVAvars;
  if(debug) cout << "event trainier initialised " << endl;
  
  // the name of the variables to be used
  MVAvars.push_back("NumberOfElectrons");
  MVAvars.push_back("NumberOfMuons");
  //   MVAvars.push_back("Zboson_M");
  
  
  if(doTraining){
    for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
    {
      Eventtrainer_->bookInputVar(MVAvars[N_var]);
    }
    if(debug) cout << "input variables booked" << endl;
  }
  if(!doTraining) Eventcomputer_ = new MVAComputer(MVAmethod,"TrainedEventMVA"+channel+".root", "TrainedEventMVA"+channel,MVAvars, "test");
  
  
  
  cout<<""<<endl;
  cout<<"RUNNING MVA: ";
  if(doTraining) cout << " training " <<endl;
  else cout << "testing " << endl;
  cout<<""<<endl;
  
  const char *xmlfile_train = xmlNom_train.c_str();
  const char *xmlfile_eval = xmlNom_evaluate.c_str();
  if(doTraining) cout << "used config file: " << xmlfile_train << endl;
  else cout << "used config file: " << xmlfile_eval << endl;
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  if(doTraining) treeLoader.LoadDatasets (datasets, xmlfile_train);	//cout<<"datasets loaded"<<endl;
  else treeLoader.LoadDatasets(datasets, xmlfile_eval);
  
  //***************************************************CREATING PLOT****************************************************
  //  MSPlot[MVAmethod.c_str()] = new MultiSamplePlot(datasets, MVAmethod.c_str() , 50, (float) -1, (float) 1, MVAmethod.c_str());
  //  MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
  
  
  // ******************************************* outputfiles ************************************************************
  
//  string branchName = "BDT";
//  string branchDescrip = channel + " additional " + branchName;
//  TNtuple* output = new TNtuple(branchDescrip.c_str(), branchDescrip.c_str(), branchName.c_str());
  
  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;
  float ScaleFactor, NormFactor;
  
  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    dataSetName = datasets[d]->Name();
    cout<<"Dataset:  :"<<dataSetName<<endl;
    filepath = CraneenPath+ dataSetName + ".root";
    if (debug) cout<<"filepath: "<<filepath<<endl;
    
    
   FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"update"); //create TFile for each dataset
    
   string TTreename = "tree"; // tree containing the variablesi
  // ttree[dataSetName.c_str()] = (TNtuple*) FileObj[dataSetName.c_str()]->Get(TTreename.c_str());
   //string tupname = "MVAoutput/"+channel+"/"+dateString+"/"+dataSetName+"BDT.root"; 
   //TFile *tupfile = new TFile(tupname.c_str(),"RECREATE");
  //  TTree* BDTtree = new TTree("bdttree","bdttree");
   // if(debug) cout << "got ntuple "<< endl;  
    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
    nEntries = ttree[dataSetName.c_str()]->GetEntries();
    cout<<"                 nEntries: "<<nEntries<<endl;
      
    /////////////////////////////////////////
    // Define variables relevant for MVA
    ////////////////////////////////////////
    int NumberOfElectrons;
    int NumberOfMuons;
    //    double Zmass;
    ttree[(dataSetName).c_str()]->SetBranchAddress("nElectrons",&NumberOfElectrons);
    ttree[(dataSetName).c_str()]->SetBranchAddress("nMuons",&NumberOfMuons);
    //    ttree[(dataSetName).c_str()]->SetBranchAddress("Zboson_M", &Zmass);
    
    bool isData= false;
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
    {
      if(debug) cout << "Data found" << endl;
      isData =true;
    }
    
    ////////////////////////////////////////////////////////////
    // Tree for reweighting
    ////////////////////////////////////////////////////////////
    // put in reweighing
    
    
//    int Train_nEntries = int(nEntries/skipEvents);
//    if(isData && !doTraining) Train_nEntries = int(nEntries/skipEvents);
    
    cout << "Number of entries: " << nEntries << endl; //", number of train Entries: " << Train_nEntries << endl;
    
    //////////////////////////////////////////////////////////
    // Running on events
    //////////////////////////////////////////////////////////
    
    if(doTraining)
    {
      
      if(isData && SignalName != "Data") continue;
      
      for (int j = 0; j<nEntries; j++)
      {
        ScaleFactor = 1.; // event scale factor
        //   ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
        //   ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
        //  if(ScaleFactor < 0) ScaleFactor = 0;
        ttree[(dataSetName).c_str()]->GetEntry(j);
        
        if(dataSetName == SignalName)
        {
          Eventtrainer_->Fill("S","NumberOfElectrons", NumberOfElectrons);
          Eventtrainer_->Fill("S","NumberOfMuons", NumberOfMuons);
          //          Eventtrainer_->Fill("S","Zmass", Zmass);
          
        }
        else
        {

          Eventtrainer_->Fill("B","NumberOfElectrons", NumberOfElectrons);
          Eventtrainer_->Fill("B","NumberOfMuons", NumberOfMuons);
          //           Eventtrainer_->Fill("B","Zmass", Zmass);
        }
      }//for-loop events
    }//If-statement doTraining
    else if(!doTraining) //not training but computing
    {
      double BDTscore;
      TBranch *newb =  ttree[(dataSetName).c_str()]->Branch("BDTscore",&BDTscore,"BDTscore/D");;
      for (int j = 0; j<nEntries; j++)
      {
        ScaleFactor = 1.; // event scale factor
        //ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
        //ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
        //if(ScaleFactor < 0) ScaleFactor = 0;
        ttree[(dataSetName).c_str()]->GetEntry(j);
        
        if (Eventcomputer_ == 0) cout <<"null computer...." <<endl;
        Eventcomputer_->FillVar("NumberOfElectrons", NumberOfElectrons);
        Eventcomputer_->FillVar("NumberOfMuons", NumberOfMuons);
        
        
        
        
       std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
          
       for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
       {
          if(fabs(it->second > 1)) cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
          BDTscore = it->second;
       }
    //   cout << "BDT score " << BDTscore << endl; 
       newb->Fill();
//        if(isData) MSPlot[MVAmethod.c_str()]->Fill(BDTscore, datasets[d], true, 1.);
//        else  MSPlot[MVAmethod.c_str()]->Fill(BDTscore, datasets[d], true, ScaleFactor);
      }
    
    }
    ttree[(dataSetName).c_str()]->Write(); 
  }//for-loop datasets
  
  if(doTraining) Eventtrainer_->TrainMVA("Random","",0,0,"",0,0,"test",false);
  
  delete Eventtrainer_;
  delete Eventcomputer_;
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
