////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <TLorentzVector.h>

//used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"




using namespace std;
using namespace reweight;
using namespace TopTree;



/// Normal Plots (TH1F* and TH2F*)
//map<string,TH1F*> histo1D;
//map<string,TH2F*> histo2D;
//map<string,TProfile*> histoProfile;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
//map<string,MultiSamplePlot*> MultiPadPlot;


// Homemade functions
std::vector <int> OSSFLeptonPairCalculator(std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 
TLorentzVector CreateZboson(std::vector<int> Lep, std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 


/// Some variables from POG/PAG
float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.


int main (int argc, char *argv[])
{
  
  string rootFileName = "ControlPlots.root";
  
  clock_t start = clock();
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool eventSelected = false;
  int nofSelectedEvents = 0;
  bool ee = false; 
  bool emu = false; 
  bool mumu = false; 
  bool runHLT = true; 
  std::string sWPMuon = "Tight"; 
  std::string sWPElectron = "Tight";
  /// xml file
  string xmlFileName ="config/Run2TriLepton_samples.xml";
  
  for(int iarg = 0; iarg < argc && argc>1 ; iarg++){
       std::string argval=argv[iarg];
       if(argval=="--help" || argval =="--h"){
   	       cout << "--ee:	Di-Electron" << endl;
   	       cout << "--emu:  Electron-Muon" << endl;
   	       cout << "--mumu: Di-Muon" << endl;
	       cout << "--noHLT: disable triggering" << endl; 
	       cout << "--WPmu WPpoint: working point muon" << endl; 
	       cout << "--WPel WPpoint: working point electron" << endl; 
// 		 cout << "--JESplus: JES sys +1 sigma MET included" << endl;
// 		  cout << "--JESminus: JES sys -1 sigma MET included" << endl;
// 		  cout << "--JERplus: JER +" << endl;
// 		  cout << "--JERminus: JER -" << endl;
// 		  cout << "--SFplus: SF up +10% syst for b quarks" << endl;
// 		  cout << "--SFminus: SF down -10% syst for b quarks" << endl;
//             cout << "--SFplus_c: SF up +10% systfor c quarks" << endl;
// 		  cout << "--SFminus_c: SF down -10% syst for c quarks" << endl;
//             cout << "--SFplus_l: SF up +10% syst for light quarks" << endl;
// 		  cout << "--SFminus_l: SF down -10% syst for light quarks" << endl;
// 		  cout << "--PUup: PU reweghting scaled up " << endl;
// 		  cout << "--PUdown: PU reweghting scaled down " << endl;
// 		  cout << "--NoPU: Do not apply pileup re-weighting" << endl;
// 		  cout << "--NoSF: Do not apply b-tag scale factor" << endl;
// 		  cout << "--RAW: Do not apply pileup re-weighting or b-tag scale factor" << endl;
               cout << "--xml dir/myxml.xml Xml file" << endl;  
   	       return 0;
       }
       if(argval=="--noHLT"){
               runHLT = false; 
       }
       if (argval=="--ee"){
   	       ee = true ;
       }
       if (argval=="--emu"){
   	       emu = true;
       }
       if (argval=="--mumu"){
   	       mumu = true ;
       }
// 	  if (argval=="--PUup" ){
// 		  PUsysUp = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--PUdown" ){
// 		  PUsysDown = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--JESplus") {
// 		  JESPlus = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--JESminus") {
// 		  JESMinus = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--JERplus") {
// 		  JERPlus = true;
// 	  }
// 	  if (argval=="--JERminus") {
// 		  JERMinus = true;
// 	  }
// 	  if (argval=="--SFplus") {
// 		  SFplus = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--SFminus"){
// 		  SFminus = true;
//             pdf = false; 
// 	  }	   
// 	  if (argval=="--NoPU") {
// 		  reweightPU = false;
//             pdf = false; 
// 	  }
// 	  if (argval=="--NoSF") {
// 		  scaleFactor = false;
//             pdf = false; 
// 	  }
// 	  if (argval=="--RAW") {
// 		  reweightPU = false; 
// 		  scaleFactor = false; 
// 		  isRAW = true;
//             pdf = false; 
// 	  }
// 	  if (argval=="--3D") {
// 		  Pu3D = true;
//             pdf = false; 
// 	  }
      if(argval=="--WPel")
      {
      	iarg++; 
	sWPElectron =argv[iarg];
      }
      if(argval=="--WPmu")
      {
      	iarg++; 
	sWPMuon = argv[iarg];
      }
      if (argval=="--xml") {
   	       iarg++;
               xmlFileName = argv[iarg];
       }
   }
     
 

  const char *xmlfile = xmlFileName.c_str();
  
  
  //Setting Lepton Channels 

  if(emu)
  {
      cout << " --> Using the Muon-Electron channel..." << endl;
      //channelpostfix = "_MuEl";
  }
  else if(ee)
  {
      cout << " --> Using the Electron-Electron channel..." << endl;
      //channelpostfix = "_ElEl";
  }
  else if(mumu)
  {
      cout << " --> Using the Muon-Muon channel..." << endl;
      //channelpostfix = "_MuMu";
  }
  else
  {
      cerr<<"ERROR: Correct Di-lepton Channel not selected."<<endl;
      exit(1);
  }  
  
  
  cout << " - Using config file " << xmlfile << endl;
  
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  
  
  ////////////////////////////////////
  ///  AnalysisEnvironment
  ////////////////////////////////////
  
  AnalysisEnvironment anaEnv;
  cout << " - Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  
  cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
  
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  
  ////////////
  /// object selection and identification
  //////////////////////
  int PVertexNdofCut = 4; // anaEnv.PVertexNdofCut; 
  int PVertexZCut =24;// anaEnv.PVertexZCut; 
  int PVertexRhoCut = 2; // anaEnv.PVertexRhoCut; 
  int MuonPtCut = 20;  //anaEnv.MuonPtCutSR; 
  int MuonEtaCut = 2.4;  //anaEnv.MuonEtaCutSR; 
  int MuonRelIsoCut = 0.15; //anaEnv.MuonRelIsoCutSR;
  std::string WPMuon = sWPMuon; // https://indico.cern.ch/event/450085/contribution/4/attachments/1169767/1688138/20151013_MOC.pdf
  std::string CampaignMuon = "Spring15"; 
  int ElectronPtCut = 20.; //anaEnv.ElectronPtCut; 
  int ElectronEtaCut = 2.4; //anaEnv.ElectronEtaCut; 
  std::string CampaignElectron = "Spring15_25ns"; 
  std::string WPElectron = sWPElectron; // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  int cutsBasedElectron = 1; 
  int JetsPtCut = 40; //anaEnv.JetsPtCutSR; 
  int applyJetID = anaEnv.applyJetID; 
  int JetsEtaCut = 2.4; 
  std::string WPJet = "Tight"; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

  /////////
  /// lumi
  /////////
  float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
  cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
  
 
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;

  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); 
  cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  
  /// //////////
  /// determine lumi
  ////////////////////////
  float Luminosity = oldLuminosity;


  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    Luminosity = datasets[d]->EquivalentLumi();
    if( dataSetName.find("WZ") == 0 ){      datasets[d]->SetColor(kBlue);   }
  }
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
   
   
   // set rootfile to store controlplots 
   TFile *fout = new TFile(rootFileName.c_str(), "RECREATE"); 
    

  
  
  //Global variable
  TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  
  ///////////////////////////////////
  /// MultisamplePlots            ///
  //////////////////////////////////
   MSPlot["cutFlow"]              = new MultiSamplePlot(datasets, "cutFlow", 15, -0.5, 14.5, "cutFlow");
   MSPlot["BasecutFlow"]              = new MultiSamplePlot(datasets, "BasecutFlow", 15, -0.5, 14.5, "Baseline cutFlow");   

   MSPlot["init_NbOfVertices"]    = new MultiSamplePlot(datasets, "init_NbOfVertices", 60, 0, 60, "initial nb of vertices");
   MSPlot["init_NbOfJets"]	  = new MultiSamplePlot(datasets, "init_NbOfJets", 16, -0.5, 15, "initial nb of jets"); 
   MSPlot["init_NbOfCSVLJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVLJets", 16, -0.5, 15, "initial nb of CSV loose jets");
   MSPlot["init_NbOfCSVMJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVMJets", 16, -0.5, 15, "initial nb of CSV medium jets"); 
   MSPlot["init_NbOfCSVTJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVTJets", 16, -0.5, 15, "initial nb of CSV tight jets");
   MSPlot["init_NbOfMuons"]       = new MultiSamplePlot(datasets, "init_NbOfMuons", 16, -0.5, 15,"initial nb of muons"); 
   MSPlot["init_NbOfElectrons"]   = new MultiSamplePlot(datasets, "init_NbOfElectrons", 16, -0.5, 15, "initial nb of electrons"); 
   MSPlot["init_NbOfLeptons"]     = new MultiSamplePlot(datasets, "init_NbOfLeptons", 16, -0.5, 15, "initial nb of leptons"); 
   
   MSPlot["3L_NbOfVertices"]      = new MultiSamplePlot(datasets, "3L_NbOfVertices", 60, 0, 60, "After 3 lepton req.:  nb of vertices");
   MSPlot["3L_NbOfJets"]          = new MultiSamplePlot(datasets, "3L_NbOfJets", 16, -0.5, 15, "After 3 lepton req.:  nb of jets");
   MSPlot["3L_NbOfCSVLJets"]      = new MultiSamplePlot(datasets, "3L_NbOfCSVLJets", 16, -0.5, 15, "After 3 lepton req.:  nb of CSV loose jets");
   MSPlot["3L_NbOfCSVMJets"]      = new MultiSamplePlot(datasets, "3L_NbOfCSVMJets", 16, -0.5, 15, "After 3 lepton req.:  nb of CSV medium jets");
   MSPlot["3L_NbOfCSVTJets"]      = new MultiSamplePlot(datasets, "3L_NbOfCSVTJets", 16, -0.5, 15, "After 3 lepton req.:  nb of CSV tight jets");
   MSPlot["3L_NbOfMuons"]         = new MultiSamplePlot(datasets, "3L_NbOfMuons", 16, -0.5, 15,"After 3 lepton req.:  nb of muons");
   MSPlot["3L_NbOfElectrons"]     = new MultiSamplePlot(datasets, "3L_NbOfElectrons", 16, -0.5, 15, "After 3 lepton req.:  nb of electrons");
   MSPlot["3L_NbOfLeptons"]       = new MultiSamplePlot(datasets, "3L_NbOfLeptons", 16, -0.5, 15, "After 3 lepton req.:  nb of leptons");

   MSPlot["2J_NbOfVertices"]      = new MultiSamplePlot(datasets, "2J_NbOfVertices", 60, 0, 60, "After 2 jet req.:  nb of vertices");
   MSPlot["2J_NbOfJets"]          = new MultiSamplePlot(datasets, "2J_NbOfJets", 16, -0.5, 15, "After 2 jet req.:  nb of jets");
   MSPlot["2J_NbOfCSVLJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVLJets", 16, -0.5, 15, "After 2 jet req.:  nb of CSV loose jets");
   MSPlot["2J_NbOfCSVMJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVMJets", 16, -0.5, 15, "After 2 jet req.:  nb of CSV medium jets");
   MSPlot["2J_NbOfCSVTJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVTJets", 16, -0.5, 15, "After 2 jet req.:  nb of CSV tight jets");
   MSPlot["2J_NbOfMuons"]         = new MultiSamplePlot(datasets, "2J_NbOfMuons", 16, -0.5, 15,"After 2 jet req.:  nb of muons");
   MSPlot["2J_NbOfElectrons"]     = new MultiSamplePlot(datasets, "2J_NbOfElectrons", 16, -0.5, 15, "After 2 jet req.:  nb of electrons");
   MSPlot["2J_NbOfLeptons"]       = new MultiSamplePlot(datasets, "2J_NbOfLeptons", 16, -0.5, 15, "After 2 jet req.:  nb of leptons");

   MSPlot["1BJ_NbOfVertices"]     = new MultiSamplePlot(datasets, "1BJ_NbOfVertices", 60, 0, 60, "After 1 Bjet req.:  nb of vertices");
   MSPlot["1BJ_NbOfJets"]         = new MultiSamplePlot(datasets, "1BJ_NbOfJets", 16, -0.5, 15, "After 1 Bjet req.:  nb of jets");
   MSPlot["1BJ_NbOfCSVLJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVLJets", 16, -0.5, 15, "After 1 Bjet req.:  nb of CSV loose jets");
   MSPlot["1BJ_NbOfCSVMJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVMJets", 16, -0.5, 15, "After 1 Bjet req.:  nb of CSV medium jets");
   MSPlot["1BJ_NbOfCSVTJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVTJets", 16, -0.5, 15, "After 1 Bjet req.:  nb of CSV tight jets");
   MSPlot["1BJ_NbOfMuons"]        = new MultiSamplePlot(datasets, "1BJ_NbOfMuons", 16, -0.5, 15,"After 1 Bjet req.:  nb of muons");
   MSPlot["1BJ_NbOfElectrons"]    = new MultiSamplePlot(datasets, "1BJ_NbOfElectrons", 16, -0.5, 15, "After 1 Bjet req.:  nb of electrons");
   MSPlot["1BJ_NbOfLeptons"]      = new MultiSamplePlot(datasets, "1BJ_NbOfLeptons", 16, -0.5, 15, "After 1 Bjet req.:  nb of leptons");

   MSPlot["OSSF_NbOfVertices"]    = new MultiSamplePlot(datasets, "OSSF_NbOfVertices", 60, 0, 60, "After OSSF req.:  nb of vertices");
   MSPlot["OSSF_NbOfJets"]        = new MultiSamplePlot(datasets, "OSSF_NbOfJets", 16, -0.5, 15, "After OSSF req.:  nb of jets");
   MSPlot["OSSF_NbOfCSVLJets"]    = new MultiSamplePlot(datasets, "OSSF_NbOfCSVLJets", 16, -0.5, 15, "After OSSF req.:  nb of CSV loose jets");
   MSPlot["OSSF_NbOfCSVMJets"]    = new MultiSamplePlot(datasets, "OSSF_NbOfCSVMJets", 16, -0.5, 15, "After OSSF req.:  nb of CSV medium jets");
   MSPlot["OSSF_NbOfCSVTJets"]    = new MultiSamplePlot(datasets, "OSSF_NbOfCSVTJets", 16, -0.5, 15, "After OSSF req.:  nb of CSV tight jets");
   MSPlot["OSSF_NbOfMuons"]       = new MultiSamplePlot(datasets, "OSSF_NbOfMuons", 16, -0.5, 15,"After OSSF req.:  nb of muons");
   MSPlot["OSSF_NbOfElectrons"]   = new MultiSamplePlot(datasets, "OSSF_NbOfElectrons", 16, -0.5, 15, "After OSSF req.:  nb of electrons");
   MSPlot["OSSF_NbOfLeptons"]     = new MultiSamplePlot(datasets, "OSSF_NbOfLeptons", 16, -0.5, 15, "After OSSF req.:  nb of leptons");

   MSPlot["ZMASS_NbOfVertices"]   = new MultiSamplePlot(datasets, "ZMASS_NbOfVertices", 60, 0, 60, "After Zmass window req.:  nb of vertices");
   MSPlot["ZMASS_NbOfJets"]       = new MultiSamplePlot(datasets, "ZMASS_NbOfJets", 16, -0.5, 15, "After Zmass window req.:  nb of jets");
   MSPlot["ZMASS_NbOfCSVLJets"]   = new MultiSamplePlot(datasets, "ZMASS_NbOfCSVLJets", 16, -0.5, 15, "After Zmass window req.:  nb of CSV loose jets");
   MSPlot["ZMASS_NbOfCSVMJets"]   = new MultiSamplePlot(datasets, "ZMASS_NbOfCSVMJets", 16, -0.5, 15, "After Zmass window req.:  nb of CSV medium jets");
   MSPlot["ZMASS_NbOfCSVTJets"]   = new MultiSamplePlot(datasets, "ZMASS_NbOfCSVTJets", 16, -0.5, 15, "After Zmass window req.:  nb of CSV tight jets");
   MSPlot["ZMASS_NbOfMuons"]      = new MultiSamplePlot(datasets, "ZMASS_NbOfMuons", 16, -0.5, 15,"After Zmass window req.:  nb of muons");
   MSPlot["ZMASS_NbOfElectrons"]  = new MultiSamplePlot(datasets, "ZMASS_NbOfElectrons", 16, -0.5, 15, "After Zmass window req.:  nb of electrons");
   MSPlot["ZMASS_NbOfLeptons"]    = new MultiSamplePlot(datasets, "ZMASS_NbOfLeptons", 16, -0.5, 15, "After Zmass window req.:  nb of leptons");

  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("preselected"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("3 selected leptons"));
  CutsSelecTable.push_back(string("At least 2 jets")); 
  CutsSelecTable.push_back(string("At least 1 CSV loose jet")); 
  CutsSelecTable.push_back(string("At least 1 OSSF lepton pair"));
  CutsSelecTable.push_back(string("Z mass window of 15 GeV")); 
   
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
    
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    int isdata = 0; 
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
    if (verbose > 1)
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    
    // make root tree file name
    string roottreename = "Ntuples/";
    roottreename+= datasets[d]->Name();
    roottreename+="_tree.root";

    cout << "  - Recreate outputfile for ntuples ... " << roottreename.c_str() << endl; 
    // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
    TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
    fileout->cd();
    
     //////////////////////////////
     // My tree - variables //
     //////////////////////////////
     // various weights
     Double_t pu_weight;
     Int_t run_num;
     Int_t evt_num;
     Int_t lumi_num;
     Int_t nvtx;
     Int_t npu;
    
     //doubles
     double ptZboson;
     double pxZboson;
     double pyZboson;
     double pzZboson;
     double etaZboson;
     double eZboson; 
     double mZboson; 

     double ptWboson_lep;
     double pxWboson_lep;
     double pyWboson_lep;
     double pzWboson_lep;
     double etaWboson_lep;
     double eWboson_lep; 
 
     //vectors
     std::vector<double> *ptMuon; 
     std::vector<double> *pxMuon; 
     std::vector<double> *pyMuon; 
     std::vector<double> *pzMuon; 
     std::vector<double> *etaMuon;
     std::vector<double> *eMuon; 
     std::vector<double> *qMuon;
     
     std::vector<double> *ptElectron; 
     std::vector<double> *pxElectron; 
     std::vector<double> *pyElectron; 
     std::vector<double> *pzElectron; 
     std::vector<double> *etaElectron; 
     std::vector<double> *eElectron; 
     std::vector<double> *qElectron;
     
     std::vector<double> *ptJet; 
     std::vector<double> *pxJet; 
     std::vector<double> *pyJet; 
     std::vector<double> *pzJet; 
     std::vector<double> *eJet; 
     std::vector<double> *etaJet;
     std::vector<double> *qJet;  
     std::vector<double> *BtagCSVjet;
     std::vector<bool> *BtagCSVL;

     std::vector<double> *ptCSVLJet;
     std::vector<double> *pxCSVLJet;
     std::vector<double> *pyCSVLJet;
     std::vector<double> *pzCSVLJet;
     std::vector<double> *eCSVLJet;
     std::vector<double> *etaCSVLJet;
     std::vector<double> *qCSVLJet;

     std::vector<double> *ptCSVMJet;
     std::vector<double> *pxCSVMJet;
     std::vector<double> *pyCSVMJet;
     std::vector<double> *pzCSVMJet;
     std::vector<double> *eCSVMJet;
     std::vector<double> *etaCSVMJet;
     std::vector<double> *qCSVMJet;

     std::vector<double> *ptCSVTJet;
     std::vector<double> *pxCSVTJet;
     std::vector<double> *pyCSVTJet;
     std::vector<double> *pzCSVTJet;
     std::vector<double> *eCSVTJet;
     std::vector<double> *etaCSVTJet;
     std::vector<double> *qCSVTJet; 

    
    
     ///////////////////////////////
     // My trees                  //
     ///////////////////////////////
     TTree *bookkeeping = new TTree("startevents","startevents");
     bookkeeping->Branch("run_num",&run_num,"run_num/I");
     bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
     bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
     bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
     bookkeeping->Branch("npu",&npu,"npu/I");


     // define the output tree
     TTree* myTree = new TTree("tree","tree");
     myTree->Branch("isdata",&isdata,"isdata/I");
     myTree->Branch("run_num",&run_num,"run_num/I");
     myTree->Branch("evt_num",&evt_num,"evt_num/I");
     myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
     myTree->Branch("nvtx",&nvtx,"nvtx/I");
     myTree->Branch("npu",&npu,"npu/I");
     
     //Set branches for doubles 
 ///   myTree -> Branch("metPt", &metPt, "metPt/D");
 //    myTree -> Branch("metPx", &metPx, "metPx/D");
 //    myTree -> Branch("metPy", &metPy, "metPy/D");
     
     myTree->Branch("ptZboson", &ptZboson,"ptZboson/D");
     myTree->Branch("pxZboson", &pxZboson,"pxZboson/D");
     myTree->Branch("pyZboson", &pyZboson,"pyZboson/D");
     myTree->Branch("pzZboson", &pzZboson,"pzZboson/D");  
     myTree->Branch("etaZboson", &etaZboson,"etaZboson/D");
     myTree->Branch("eZboson", &eZboson,"eZboson/D");
     myTree->Branch("mZboson", &mZboson,"mZboson/D");
     
     myTree->Branch("ptWboson_lep", &ptWboson_lep,"ptWboson_lep/D");
     myTree->Branch("pxWboson_lep", &pxWboson_lep,"pxWboson_lep/D");
     myTree->Branch("pyWboson_lep", &pyWboson_lep,"pyWboson_lep/D");
     myTree->Branch("pzWboson_lep", &pzWboson_lep,"pzWboson_lep/D");
     myTree->Branch("etaWboson_lep", &etaWboson_lep,"etaWboson_lep/D");
     myTree->Branch("eWboson_lep", &eWboson_lep,"eWboson_lep/D");
    
     // Set branches for vectors
     myTree->Branch("ptMuon","std::vector<double>",&ptMuon);   // Make a branch with name ptMuon, from type vector(double), on loaction ptMuon
     myTree->Branch("pxMuon","std::vector<double>",&pxMuon);
     myTree->Branch("pyMuon","std::vector<double>",&pyMuon);
     myTree->Branch("pzMuon","std::vector<double>",&pzMuon);
     myTree->Branch("eMuon","std::vector<double>",&eMuon);
     myTree->Branch("etaMuon","std::vector<double>",&etaMuon);
     myTree->Branch("qMuon","std::vector<double>",&qMuon);
     
     myTree->Branch("ptElectron","std::vector<double>",&ptElectron);   // Make a branch with name ptElectron, from type vector(double), on loaction ptElectron
     myTree->Branch("pxElectron","std::vector<double>",&pxElectron);
     myTree->Branch("pyElectron","std::vector<double>",&pyElectron);
     myTree->Branch("pzElectron","std::vector<double>",&pzElectron);
     myTree->Branch("eElectron","std::vector<double>",&eElectron);
     myTree->Branch("etaElectron","std::vector<double>",&etaElectron);
     myTree->Branch("qElectron","std::vector<double>",&qElectron);     
     
     myTree->Branch("ptJet","std::vector<double>",&ptJet);
     myTree->Branch("pxJet","std::vector<double>",&pxJet);
     myTree->Branch("pyJet","std::vector<double>",&pyJet);
     myTree->Branch("pzJet","std::vector<double>",&pzJet);
     myTree->Branch("eJet","std::vector<double>",&eJet);
     myTree->Branch("etaJet","std::vector<double>",&etaJet);
     myTree->Branch("qJet","std::vector<double>",&qJet);
     myTree->Branch("BtagCSVjet", "std::vector<double>",&BtagCSVjet);
     myTree->Branch("BtagCSVL","std::vector<bool>",&BtagCSVL);
    

     myTree->Branch("ptCSVLJet","std::vector<double>",&ptCSVLJet);
     myTree->Branch("pxCSVLJet","std::vector<double>",&pxCSVLJet);
     myTree->Branch("pyCSVLJet","std::vector<double>",&pyCSVLJet);
     myTree->Branch("pzCSVLJet","std::vector<double>",&pzCSVLJet);
     myTree->Branch("eCSVLJet","std::vector<double>",&eCSVLJet);
     myTree->Branch("etaCSVLJet","std::vector<double>",&etaCSVLJet);
     myTree->Branch("qCSVLJet","std::vector<double>",&qCSVLJet);


     myTree->Branch("ptCSVMJet","std::vector<double>",&ptCSVMJet);
     myTree->Branch("pxCSVMJet","std::vector<double>",&pxCSVMJet);
     myTree->Branch("pyCSVMJet","std::vector<double>",&pyCSVMJet);
     myTree->Branch("pzCSVMJet","std::vector<double>",&pzCSVMJet);
     myTree->Branch("eCSVMJet","std::vector<double>",&eCSVMJet);
     myTree->Branch("etaCSVMJet","std::vector<double>",&etaCSVMJet);
     myTree->Branch("qCSVMJet","std::vector<double>",&qCSVMJet);


     myTree->Branch("ptCSVTJet","std::vector<double>",&ptCSVTJet);
     myTree->Branch("pxCSVTJet","std::vector<double>",&pxCSVTJet);
     myTree->Branch("pyCSVTJet","std::vector<double>",&pyCSVTJet);
     myTree->Branch("pzCSVTJet","std::vector<double>",&pzCSVTJet);
     myTree->Branch("eCSVTJet","std::vector<double>",&eCSVTJet);
     myTree->Branch("etaCSVTJet","std::vector<double>",&etaCSVTJet);
     myTree->Branch("qCSVTJet","std::vector<double>",&qCSVTJet);


    
    
//     TTree *noselTree = new TTree("noselTree","noselTree");
/*      noselTree->Branch("isdata",&isdata,"isdata/I");
     noselTree->Branch("run_num",&run_num,"run_num/I");
     noselTree->Branch("evt_num",&evt_num,"evt_num/I");
     noselTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
     noselTree->Branch("nvtx",&nvtx,"nvtx/I");
 */    // noselTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
    
    
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    
    nofSelectedEvents = 0; 
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    //some bookkeeping variables
    nEvents[d] = 0;
    int previousRun = -1;
    int itrigger = -1; 
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 2000; ievt++)
    {
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      //vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%500 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);


      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if (previousFilename != currentFilename)
      {
        previousFilename = currentFilename;
        iFile++;
        cout << "WARNING: File changed => iFile = " << iFile << endl;
      }
      
      
      int currentRun = event->runId(); 
 
      run_num=event->runId();
      evt_num=event->eventId();
      lumi_num=event->lumiBlockId();
      nvtx=vertex.size();
      npu=(int)event->nTruePU();
      if( run_num > 10000){//data
         isdata=1;
      }
      bookkeeping->Fill(); 
      
      /////////////////////////////
      /// Trigger
      ///////////////////////////
      
      bool trigged = false; 
      //If the HLT is applied 
      if(runHLT && previousRun != currentRun){
        //The HLT is only used for data
        if(isdata == 1){
          //The HLT path is dependent of the mode, these paths are the several steps or software modules. Each module performs a well defined task 
          // such as reconstruction of physics objects, making intermediate decisions, triggering more refined reconstructions in subsequent modules, 
          // or calculating the final decision for that trigger path.
	  if(emu) itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"), currentRun, iFile);
          else if(mumu) itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2"), currentRun, iFile);
          else if(ee) itrigger = treeLoader.iTrigger (string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2"), currentRun, iFile);

          if(itrigger == 9999) 
	  {
	      cout << "ERROR: no valid trigger found for this event/data in run " << event->runId() << endl; 
	  }   
	  
        } // closing the HLT for data loop
        //For the MC, there is no triggerpath
        else itrigger = true;     	
      } // closing the HLT run loop
      
      ////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      
      // PU reweighting
      
      // old method
      //cout << "scalefactor " << scaleFactor << endl;
      double lumiWeight = 1; //LumiWeights.ITweight( (int)event->nTruePU() ); // currently no pile-up reweighting applied
      
      if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
      
      scaleFactor = scaleFactor*lumiWeight;
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool isGoodPV = selection.isPVSelected(vertex, PVertexNdofCut, PVertexZCut,PVertexRhoCut); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(JetsPtCut, JetsEtaCut, applyJetID, WPJet);  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(MuonPtCut, MuonEtaCut, MuonRelIsoCut,WPMuon,CampaignMuon);  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(ElectronPtCut, ElectronEtaCut, WPElectron, CampaignElectron, cutsBasedElectron);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      
      sort(selectedJets.begin(), selectedJets.end(),HighestPt()); 
      sort(selectedMuons.begin(), selectedMuons.end(), HighestPt()); 
      sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt()); 
      
      vector<bool> BtagBooleans; 
      BtagBooleans.clear(); 
      vector<TRootJet*> selectedBCSVLJets; 
      vector<TRootJet*> selectedBCSVMJets;
      vector<TRootJet*> selectedBCSVTJets;
      for(unsigned int i = 0; i < selectedJets.size() ; i++)
      {
         bool Btagged = false;
         TRootJet* tempJet = (TRootJet*) selectedJets[i];
         if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose)//loose WP
         {
           Btagged = true;
           selectedBCSVLJets.push_back(tempJet); 
         }
         BtagBooleans.push_back(Btagged);
         if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium)//medium WP
         {
           selectedBCSVMJets.push_back(tempJet);
         }
         if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight)//tight WP
         {
           selectedBCSVTJets.push_back(tempJet);
         }
       }
       if(verbose > 3) cout << "btagging done" << endl; 
      
      // Start analysis selection
      eventSelected = false;
      TLorentzVector Zboson;
      TLorentzVector Wlep; 
      
      
      /// Initial nbrs
      
      selecTable.Fill(d,0,scaleFactor);
      MSPlot["cutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor );
      MSPlot["init_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor); 
      MSPlot["init_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["init_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor); 
      MSPlot["init_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor); 
      MSPlot["init_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor); 
      
      
      /// Trigger
      if(runHLT) trigged = treeLoader.EventTrigged(itrigger);
      else trigged = true; 
       
      if(trigged)
      { 
       selecTable.Fill(d,1,scaleFactor);
       MSPlot["cutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor );
      
       if (isGoodPV)
       {
        if(verbose>3) cout << "GoodPV" << endl; 
        selecTable.Fill(d,2,scaleFactor);
	MSPlot["cutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor );

        if (selectedMuons.size() + selectedElectrons.size()== 3)
        {
            if(verbose>3) cout << "3 electrons "<< endl; 
 	    selecTable.Fill(d,3,scaleFactor);
	    MSPlot["cutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor );
            MSPlot["BasecutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor );
            MSPlot["3L_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["3L_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);	    
	    
	    if(selectedJets.size() > 1)
	    {
	       if(verbose>3) cout << " at least 2 jets " << endl; 
	       selecTable.Fill(d,4,scaleFactor); 
	       MSPlot["cutFlow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor );
	       MSPlot["BasecutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor );
               MSPlot["2J_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
               MSPlot["2J_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);

	       if(selectedBCSVLJets.size()>1)
	       {
	         if(verbose>3) cout << " at least 1 bjet " << endl; 
	         selecTable.Fill(d,5,scaleFactor); 
		 MSPlot["cutFlow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor );
 		 MSPlot["BasecutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor );
                 MSPlot["1BJ_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                 MSPlot["1BJ_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
     		 MSPlot["1BJ_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
	         MSPlot["1BJ_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
		 MSPlot["1BJ_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
                 MSPlot["1BJ_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                 MSPlot["1BJ_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
                 MSPlot["1BJ_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);	 
		 
		 std::vector<int> Leptons; 
		 Leptons.clear(); 
		 Leptons = OSSFLeptonPairCalculator(selectedElectrons, selectedMuons, verbose);
		 if(verbose>3) cout <<   Leptons[0]<< " , " <<  Leptons[1]<< " , " <<  Leptons[2]<< " , " <<  Leptons[3]<< " , " <<  Leptons[4]<< " , " <<  Leptons[5]   << endl; 
		 
		 bool OSSFpair = false; 
		 if( (Leptons[0] != -5 && Leptons[1] != -5) | (Leptons[3] != -5 && Leptons[4] != -5) ) OSSFpair = true; 
		 if(OSSFpair)
	         {
		    if(verbose>3) cout << " OSSF "<< endl; 
		    selecTable.Fill(d,6,scaleFactor);
		    MSPlot["cutFlow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor );
		    MSPlot["BasecutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor );
	            MSPlot["OSSF_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
     		    MSPlot["OSSF_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
      		    MSPlot["OSSF_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
      		    MSPlot["OSSF_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
      		    MSPlot["OSSF_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
    		    MSPlot["OSSF_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
     		    MSPlot["OSSF_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["OSSF_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);   
		    
		    
		    Zboson.Clear(); 
		    Zboson = CreateZboson(Leptons, selectedElectrons, selectedMuons, verbose); 
		    
                    Wlep.Clear(); 
                    if(fabs(Leptons[2]) != 5)
		    {
		      if(verbose>3) cout << " the W lepton is an electron " << endl; 
		      Wlep.SetPxPyPzE(selectedElectrons[Leptons[2]]->Px(), selectedElectrons[Leptons[2]]->Py(), selectedElectrons[Leptons[2]]->Pz(), selectedElectrons[Leptons[2]]->Energy()); 
                    }
		    else if(fabs(Leptons[5]) != 5)
		    {
		      if(verbose>3) cout << " the W lepton is a muon " << endl; 
		      Wlep.SetPxPyPzE(selectedMuons[Leptons[5]]->Px(), selectedMuons[Leptons[5]]->Py(), selectedMuons[Leptons[5]]->Pz(), selectedMuons[Leptons[5]]->Energy()); 		    
		    }
		    bool ZmassWindow = false; 
		    if(fabs(Zboson.M()-90.0) < 15.0) ZmassWindow = true; 
		    if(ZmassWindow)
		    {
		      if(verbose>3) cout << " Zmass window " << endl; 
		      selecTable.Fill(d,7,scaleFactor);
		      MSPlot["cutFlow"]->Fill(7, datasets[d], true, Luminosity*scaleFactor );
		      MSPlot["BasecutFlow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor );
                      MSPlot["ZMASS_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
  		      MSPlot["ZMASS_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
		      MSPlot["ZMASS_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);

		      eventSelected = true;
	            }
		 }
	       } 
	    }
         
        }  // 3 leptons
       }  // good PV
      } // trigger
      
      
      if (! eventSelected )
      {
        continue;
      }
      
      if(verbose>3) cout << "filling the tree" << endl; 
      ptMuon = new std::vector<double>; 
      pxMuon = new std::vector<double>; 
      pyMuon = new std::vector<double>; 
      pzMuon = new std::vector<double>; 
      etaMuon = new std::vector<double>;
      eMuon = new std::vector<double>; 
      qMuon = new std::vector<double>;
      
      ptElectron = new std::vector<double>; 
      pxElectron = new std::vector<double>; 
      pyElectron = new std::vector<double>; 
      pzElectron = new std::vector<double>; 
      etaElectron = new std::vector<double>; 
      eElectron = new std::vector<double>; 
      qElectron = new std::vector<double>;
      
 
      ptJet = new std::vector<double>; 
      pxJet = new std::vector<double>; 
      pyJet = new std::vector<double>; 
      pzJet = new std::vector<double>; 
      eJet = new std::vector<double>; 
      etaJet = new std::vector<double>; 
      qJet = new std::vector<double>; 
      BtagCSVjet = new std::vector<double>;
      BtagCSVL = new std::vector<bool>; 
      
      ptCSVLJet = new std::vector<double>;
      pxCSVLJet = new std::vector<double>;
      pyCSVLJet = new std::vector<double>;
      pzCSVLJet = new std::vector<double>;
      eCSVLJet = new std::vector<double>;
      etaCSVLJet = new std::vector<double>;
      qCSVLJet = new std::vector<double>;

      ptCSVMJet = new std::vector<double>;
      pxCSVMJet = new std::vector<double>;
      pyCSVMJet = new std::vector<double>;
      pzCSVMJet = new std::vector<double>;
      eCSVMJet = new std::vector<double>;
      etaCSVMJet = new std::vector<double>;
      qCSVMJet = new std::vector<double>;

      ptCSVTJet = new std::vector<double>;
      pxCSVTJet = new std::vector<double>;
      pyCSVTJet = new std::vector<double>;
      pzCSVTJet = new std::vector<double>;
      eCSVTJet = new std::vector<double>;
      etaCSVTJet = new std::vector<double>;
      qCSVTJet = new std::vector<double>;
      
      ptZboson = Zboson.Pt();
      pxZboson = Zboson.Px();
      pyZboson = Zboson.Py();
      pzZboson = Zboson.Pt();
      etaZboson = Zboson.Eta();
      eZboson = Zboson.Energy();
      mZboson = Zboson.M();  

      ptWboson_lep = Wlep.Pt();
      pxWboson_lep = Wlep.Px();
      pyWboson_lep = Wlep.Py();
      pzWboson_lep = Wlep.Pz();
      etaWboson_lep = Wlep.Eta();
      eWboson_lep = Wlep.Energy();
      
      for (unsigned int i = 0; i < selectedElectrons.size(); i++) 
      {
      	TRootElectron* tempElectron = (TRootElectron*) selectedElectrons[i];
        ptElectron->push_back(tempElectron->Pt()); 
	pxElectron->push_back(tempElectron->Px()); 
	pyElectron->push_back(tempElectron->Py()); 
	pzElectron->push_back(tempElectron->Pz()); 
	eElectron->push_back(tempElectron->Energy());
	etaElectron->push_back(tempElectron->Eta()); 
	qElectron->push_back(tempElectron->charge());
      }
      
      for (unsigned int i = 0; i < selectedMuons.size(); i++) 
      {
      	TRootMuon* tempMuon = (TRootMuon*) selectedMuons[i];
        ptMuon->push_back(tempMuon->Pt()); 
	pxMuon->push_back(tempMuon->Px()); 
	pyMuon->push_back(tempMuon->Py()); 
	pzMuon->push_back(tempMuon->Pz()); 
	eMuon->push_back(tempMuon->Energy());
	etaMuon->push_back(tempMuon->Eta()); 
	qMuon->push_back(tempMuon->charge());
      }
            
      for (unsigned int i =0; i < selectedJets.size(); i ++)
      {
        TRootJet* tempJet = (TRootJet*) selectedJets[i];
        ptJet->push_back(tempJet->Pt());
        pxJet->push_back(tempJet->Px());
        pyJet->push_back(tempJet->Py());
        pzJet->push_back(tempJet->Pz());
        eJet->push_back(tempJet->Energy());
        etaJet->push_back(tempJet->Eta());
        qJet->push_back(tempJet->charge());
        BtagCSVjet->push_back(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags()); 
        BtagCSVL->push_back(BtagBooleans[i]); 
      }
      for (unsigned int i = 0; i< selectedBCSVLJets.size(); i++)
      {
        TRootJet* tempJet = (TRootJet*) selectedBCSVLJets[i];
	ptCSVLJet->push_back(tempJet->Pt());
        pxCSVLJet->push_back(tempJet->Px());
        pyCSVLJet->push_back(tempJet->Py());
        pzCSVLJet->push_back(tempJet->Pz());
        eCSVLJet->push_back(tempJet->Energy());
        etaCSVLJet->push_back(tempJet->Eta());
        qCSVLJet->push_back(tempJet->charge());
      }
      for (unsigned int i = 0; i< selectedBCSVMJets.size(); i++)
      {
        TRootJet* tempJet = (TRootJet*) selectedBCSVMJets[i];
        ptCSVMJet->push_back(tempJet->Pt());
        pxCSVMJet->push_back(tempJet->Px());
        pyCSVMJet->push_back(tempJet->Py());
        pzCSVMJet->push_back(tempJet->Pz());
        eCSVMJet->push_back(tempJet->Energy());
        etaCSVMJet->push_back(tempJet->Eta());
        qCSVMJet->push_back(tempJet->charge());
      }
      for (unsigned int i = 0; i< selectedBCSVTJets.size(); i++)
      {
        TRootJet* tempJet = (TRootJet*) selectedBCSVTJets[i];
        ptCSVTJet->push_back(tempJet->Pt());
        pxCSVTJet->push_back(tempJet->Px());
        pyCSVTJet->push_back(tempJet->Py());
        pzCSVTJet->push_back(tempJet->Pz());
        eCSVTJet->push_back(tempJet->Energy());
        etaCSVTJet->push_back(tempJet->Eta());
        qCSVTJet->push_back(tempJet->charge());
      }


      nofSelectedEvents++; 
      myTree->Fill();
      
      delete pxElectron;
      delete pyElectron;
      delete pzElectron;
      delete etaElectron; 
      delete eElectron;
      delete qElectron;
      
      delete ptMuon;
      delete pxMuon;
      delete pyMuon;
      delete pzMuon;
      delete etaMuon; 
      delete eMuon;
      delete qMuon;
            
      delete ptJet;
      delete pxJet;
      delete pyJet;
      delete pzJet;
      delete eJet;
      delete etaJet;
      delete qJet;
      delete BtagCSVjet;
      delete BtagCSVL; 
      
      delete ptCSVLJet;
      delete pxCSVLJet;
      delete pyCSVLJet;
      delete pzCSVLJet;
      delete eCSVLJet;
      delete etaCSVLJet;
      delete qCSVLJet;

      delete ptCSVMJet;
      delete pxCSVMJet;
      delete pyCSVMJet;
      delete pzCSVMJet;
      delete eCSVMJet;
      delete etaCSVMJet;
      delete qCSVMJet;

      delete ptCSVTJet;
      delete pxCSVTJet;
      delete pyCSVTJet;
      delete pzCSVTJet;
      delete eCSVTJet;
      delete etaCSVTJet;
      delete qCSVTJet;
      
      if (verbose > 2)
        cout << "  Event " << ievt << " is selected" << endl;
      

      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  //loop on events
   
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    
    selecTable.TableCalculator(false, true, true, true, true);
    string selectionTable = "SelectionTables/SelectionTable_"+datasets[d]->Name() +".tex";
    selecTable.Write(selectionTable.c_str(), true, true, true, true, true, true, false);
    
    myTree->Write();
    fileout->Write();
    fileout->Close();
    delete fileout; 

   
    
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // Loop on datasets
  
  
  
  
    ///////////////
    /// TABLES
    //////////////////
   selecTable.TableCalculator(false, true, true, true, true); 
   string selectionTableAll = "SelectionTables/SelectionTable_allSamples.tex";
   selecTable.Write(selectionTableAll.c_str(), true, true, true, true, true, true, false);
   
   
   ///*****************///
  ///   Write plots   ///
  ///*****************///
  
  string pathPNG = "ControlPlots/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"MSPlot/").c_str(),0777); // 0777 if it doesn't exist already, make it
  
  ///Write histograms
  fout->cd();
  for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    if(verbose>3) cout << "MSPlot: " << it->first << endl;
    MultiSamplePlot *temp = it->second;
    string name = it->first;
    temp->Draw(name, 1, false, false, false, 1);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    temp->Write(fout, name, true, pathPNG+"MSPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
   
   
   
   ///////////////////
   /// CLEANING
   /////////////////
  
  fout->Close();
  
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
 
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
  
}


std::vector<int> OSSFLeptonPairCalculator(vector<TRootElectron*> electrons, vector<TRootMuon*> muons, int verbose)
{
  if(verbose > 3) cout << "In OSSFLeptonPairCalculator " << endl; 
  std::vector<int> leptons;
  leptons.clear();  
  int ElecZ0 = -5; 
  int ElecZ1 = -5; 
  int ElecW = -5; 
  int MuZ0 = -5; 
  int MuZ1 = -5; 
  int MuW = -5;
  
  if(electrons.size() == 3)
  {
    if(electrons[0]->charge() != electrons[1]->charge())
    {
      ElecZ0 = 0; 
      ElecZ1 = 1;
      ElecW = 2; 
      if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
    }    
    else if (electrons[0]->charge() != electrons[2]->charge())
    {
      ElecZ0 = 0; 
      ElecZ1 = 2; 
      ElecW = 1; 
      if(verbose>3) cout << " the Zboson consistes of electrons " << endl; 
    }
    else if (electrons[1]->charge() != electrons[2]->charge())
    {
      ElecZ0 = 1; 
      ElecZ1 = 2; 
      ElecW = 0;
      if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
    }

  }
  else if ( (electrons.size() == 2) && (electrons[0]->charge() != electrons[1]->charge()) ) 
  {
    ElecZ0 = 0; 
    ElecZ1 = 1; 
    MuW = 0; 
    if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
  }
  else if ( (muons.size() == 2) && (muons[0]->charge() != muons[1]->charge()) )
  {
    ElecW = 0; 
    MuZ0 = 0; 
    MuZ1 = 1;
    if(verbose>3) cout << " the Zboson consists of muons " << endl; 
  }
  else if (muons.size() == 3)
  {
    if(muons[0]->charge() != muons[1]->charge())
    {
      MuZ0 = 0; 
      MuZ1 = 1;
      MuW = 2; 
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }    
    else if (muons[0]->charge() != muons[2]->charge())
    {
      MuZ0 = 0; 
      MuZ1 = 2; 
      MuW = 1; 
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }
    else if (muons[1]->charge() != muons[2]->charge())
    {
      MuZ0 = 1; 
      MuZ1 = 2; 
      MuW = 0;
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }
     
  }
  leptons.push_back(ElecZ0); 
  leptons.push_back(ElecZ1); 
  leptons.push_back(ElecW); 
  leptons.push_back(MuZ0); 
  leptons.push_back(MuZ1); 
  leptons.push_back(MuW); 
  if(verbose>3) cout << " out OSSF.. " << endl; 
  if(verbose>3) cout <<   leptons[0]<< " , " <<  leptons[1]<< " , " <<  leptons[2]<< " , " <<  leptons[3]<< " , " <<  leptons[4]<< " , " <<  leptons[5]   << endl; 
  
  return leptons; 
}


TLorentzVector CreateZboson(std::vector<int> leptons, std::vector<TRootElectron*> electrons, std::vector<TRootMuon*> muons, int verbose)
{
  if(verbose>3) cout << " in Zboson creator " << endl; 
  TLorentzVector Zbos;
  Zbos.Clear(); 
  TLorentzVector Zbos_lep0; 
  Zbos_lep0.Clear(); 
  TLorentzVector Zbos_lep1; 
  Zbos_lep1.Clear(); 
  
  if(verbose>3) cout <<   leptons[0]<< " , " <<  leptons[1]<< " , " <<  leptons[2]<< " , " <<  leptons[3]<< " , " <<  leptons[4]<< " , " <<  leptons[5]   << endl; 
  
  if(fabs(leptons[0]) < 3 && fabs(leptons[1]) < 3) 
  {
      Zbos_lep0.SetPxPyPzE(electrons[leptons[0]]->Px(), electrons[leptons[0]]->Py(), electrons[leptons[0]]->Pz(), electrons[leptons[0]]->Energy()); 
      Zbos_lep1.SetPxPyPzE(electrons[leptons[1]]->Px(), electrons[leptons[1]]->Py(),electrons[leptons[1]]->Pz(), electrons[leptons[1]]->Energy()); 
  }
  else if(fabs(leptons[3]) < 3 && fabs(leptons[4]) < 3)
  {
     Zbos_lep0.SetPxPyPzE(muons[leptons[3]]->Px(),muons[leptons[3]]->Py(), muons[leptons[3]]->Pz(),muons[leptons[3]]->Energy()); 
     Zbos_lep1.SetPxPyPzE(muons[leptons[4]]->Px(),muons[leptons[4]]->Py(), muons[leptons[4]]->Pz(),muons[leptons[4]]->Energy());     
  }

  Zbos = Zbos_lep0 + Zbos_lep1; 
  if(verbose>3) cout << " out Zboson creator " <<endl; 
  return Zbos; 
}
