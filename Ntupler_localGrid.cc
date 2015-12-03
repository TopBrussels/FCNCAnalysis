
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
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
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
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
//This header file is taken directly from the BTV wiki. It contains
//
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"


using namespace std;
using namespace reweight;
using namespace TopTree;




// Homemade functions
std::vector <int> OSSFLeptonPairCalculator(std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 
TLorentzVector CreateZboson(std::vector<int> Lep, std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 


/// Some variables from POG/PAG
float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.


int main (int argc, char *argv[])
{
  
  clock_t start = clock();


  

  
    
  /////////////////////
  ///  Configuration
  /////////////////////

  bool eventSelected = false;
  int nofSelectedEvents = 0;
  bool ee = false; 
  bool emu = false; 
  bool mumu = true; 
  bool runHLT = false; 
  std::string sWPMuon = "Tight"; 
  std::string sWPElectron = "Medium";
  /// xml file
  string xmlFileName ="config/Run2TriLepton_samples.xml";
  float Luminosity = 1263.885980236;  ; //  rereco run D + prompt v4 
  const char *xmlfile = xmlFileName.c_str();
  std::string channelpostfix = ""; 
  //UNCERTAINTIES
  bool doJESup = false; 
  bool doJESdown= false; 
  bool doJERup = false; 
  bool doJERdown =false;  
  //Setting Lepton Channels 

  if(emu)
  {
      cout << " --> Using the Muon-Electron channel..." << endl;
      channelpostfix = "_MuEl_";
  }
  else if(ee)
  {
      cout << " --> Using the Electron-Electron channel..." << endl;
      channelpostfix = "_ElEl_";
  }
  else if(mumu)
  {
      cout << " --> Using the Muon-Muon channel..." << endl;
      channelpostfix = "_MuMu_";
  }
  else
  {
      cerr<<"ERROR: Correct Di-lepton Channel not selected."<<endl;
      exit(1);
  }    
  cout << " - Using config file " << xmlfile << endl;

  /// setting all arguments to be able to run on local grid
  // check if all arguments are there
  if(argc < 15)
  {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
    for (int n_arg=1; n_arg<argc; n_arg++)
    {
       std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl; 
    }
    return 1;
  }
  // placing argument in variables
  const string dName              = argv[1];
  const string dTitle             = argv[2];
  const int color                 = strtol(argv[4], NULL, 10);
  const int ls                    = strtol(argv[5], NULL, 10);
  const int lw                    = strtol(argv[6], NULL, 10);
  const float normf               = strtod(argv[7], NULL);
  const float EqLumi              = strtod(argv[8], NULL);
  const float xSect               = strtod(argv[9], NULL);
  const float PreselEff           = strtod(argv[10], NULL);
  string fileName                 = argv[11];
  const int startEvent            = strtol(argv[argc-3], NULL, 10);
  const int endEvent              = strtol(argv[argc-2], NULL, 10);
  const int JobNum                = strtol(argv[argc-1], NULL, 10);
  


    
  vector<string> vecfileNames;
  for(int args = 11; args < argc-3; args++)
  {
     vecfileNames.push_back(argv[args]);
  }


  //info
    cout << "---Dataset accepted from command line---" << endl;
    cout << "Dataset Name: " << dName << endl;
    cout << "Dataset Title: " << dTitle << endl;
    cout << "Dataset color: " << color << endl;
    cout << "Dataset ls: " << ls << endl;
    cout << "Dataset lw: " << lw << endl;
    cout << "Dataset normf: " << normf << endl;
    cout << "Dataset EqLumi: " << EqLumi << endl;
    cout << "Dataset xSect: " << xSect << endl;
    cout << "Dataset File Name: " << vecfileNames[0] << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "JobNum: " << JobNum << endl;
    

  
  //Configuration output format
//  TTree *configTree = new TTree("configTree","configuration Tree");
//  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
//  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
//  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  
  
  ////////////////////////////////////
  ///  AnalysisEnvironment
  ////////////////////////////////////
  
  AnalysisEnvironment anaEnv;
  cout << " - Loading environment ..." << endl;
  // AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile); doesn't work on localgird
  anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = true;
    anaEnv.loadGenJetCollection = true;// changed on 31okt
    anaEnv.loadGenEventCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.loadTrackMETCollection = false;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;  



  cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
  
  //new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = 2; // anaEnv.Verbose;
  
  ////////////
  /// object selection and identification
  //////////////////////
/*  int PVertexNdofCut = 4; // anaEnv.PVertexNdofCut; 
  int PVertexZCut =24;// anaEnv.PVertexZCut; 
  int PVertexRhoCut = 2; // anaEnv.PVertexRhoCut; 
  int MuonPtCut = 20;  //anaEnv.MuonPtCutSR; 
  int MuonEtaCut = 2.4;  //anaEnv.MuonEtaCutSR; 
  int MuonRelIsoCut = 0.15; //anaEnv.MuonRelIsoCutSR;
  std::string WPMuon = sWPMuon; // https://indico.cern.ch/event/450085/contribution/4/attachments/1169767/1688138/20151013_MOC.pdf
  std::string CampaignMuon = "Spring15"; 
  int ElectronPtCut = 20.; //anaEnv.ElectronPtCut; 
  int ElectronEtaCut = 2.5; //anaEnv.ElectronEtaCut; 
  std::string CampaignElectron = "Spring15_25ns"; 
  std::string WPElectron = sWPElectron; // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  int cutsBasedElectron = 1; 
  int JetsPtCut = 30; //anaEnv.JetsPtCutSR; 
  int applyJetID = true; //anaEnv.applyJetID; 
  int JetsEtaCut = 2.4; 
  std::string WPJet = "Tight"; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
*/
  /////////
  /// lumi
  /////////
  float oldLuminosity =  1263.885980236;    //anaEnv.Luminosity;  // in 1/pb
  cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
  
 
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  cout << " - Creating Dataset ..." << endl;
  Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
  theDataset->SetEquivalentLuminosity(EqLumi*normf);
  datasets.push_back(theDataset);
  cout << "Number of datasets: " << datasets.size() << endl;

  
  /// //////////
  /// determine lumi
  ////////////////////////
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    string dataSetName = datasets[d]->Name();
    if(dataSetName.find("Data")==0 || dataSetName.find("data")==0 || dataSetName.find("DATA")==0)
    {
        cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
    }   
  }
//  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
    cout << "lumi is " << Luminosity << endl;  
   stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str(); 
  
  
  //Global variable
  TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  

  ///////////////////////
  //  Event SF
  ///////////////
  Double_t scaleFactor = 1.;

  Double_t muonScaleFactor, electronScaleFactor, puScaleFactor, btagScaleFactor;
  muonScaleFactor = electronScaleFactor = puScaleFactor = btagScaleFactor = 1.0;

  Bool_t applyMuonSF , applyElectronSF, applyPUSF, applyGlobalSF, applyBtagSF, fillingbTagHistos, applyJetCleaning;
  applyMuonSF = false;
  applyElectronSF = false;
  applyPUSF = false;
  applyGlobalSF = false;
  applyBtagSF = false; // doesn't work in 74X
  fillingbTagHistos = false;
  applyJetCleaning = false; 
  string pathToCaliDir = "/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/TopTreeAnalysisBase/Calibrations/";

  //Muon SF 
  string muonFile= "Muon_SF_TopEA.root";
  MuonSFWeight *muonSFWeight = new MuonSFWeight (pathToCaliDir+"LeptonSF/"+muonFile,"SF_totErr", false, false); // (... , ... , debug, print warning)
  
  //Electron SF
  string electronFile= "Elec_SF_TopEA.root";
  ElectronSFWeight *electronSFWeight = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electronFile,"GlobalSF", false, false); // (... , ... , debug, print warning)

  //PU SF 
  LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root","pileup","pileup");
 
  // Btag SF 
  BTagCalibration * bTagCalib;   
  BTagCalibrationReader * bTagReader;
  BTagWeightTools *btwt;
  
 
  
 
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    cout << "luminosity of dataset " << Luminosity << endl; 
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
 
    int isdata = 0; 
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
    if (verbose > 1)
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    double lumiWeight = -99.;
    if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0){
        lumiWeight=1;
	isdata = 1;
    }
    else{
        lumiWeight = Luminosity/ datasets[d]->EquivalentLumi();
        cout << "the weight to apply for each event of this data set is " << "Lumi * (xs/NSample) = Lumi / EqLumi = "  << Luminosity  /  datasets[d]->EquivalentLumi()  <<  endl;
      }
    
     
    
    
    // Controlplots
   // set rootfile to store controlplots
    string pathRoot = "ControlPlots/";
    mkdir(pathRoot.c_str(),0777); 
    string pathRootSample = channelpostfix + datasets[d]->Name();
    if(doJERup) pathRootSample += "_JERup";
    if(doJERdown) pathRootSample += "_JERdown";
    if(doJESup) pathRootSample += "_JERup";
    if(doJESdown) pathRootSample += "_JERdown";
   
    mkdir((pathRoot+pathRootSample+"/").c_str(),0777);
    string rootFileName = pathRoot + pathRootSample +"/ControlPlots"+channelpostfix + datasets[d]->Name() + "_" + strJobNum+".root";
    TFile *fout = new TFile(rootFileName.c_str(), "RECREATE");
    
   
       ////////////////////////////
  /// HISTO
  ////////////////////////////:
    map <string,TH1F*> histo1D;
    map <string,TH2F*> histo2D;
 
    std::string  titlePlot = ""; 
    titlePlot = "initial_Nb_Jets"+channelpostfix; 	
    histo1D["h_initial_Nb_Jets"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets",  16, - 0.5, 15.5 ); 
    titlePlot = "3L_Nb_Jets"+channelpostfix;
    histo1D["h_3L_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After 3L cut: nb. of jets",  16, - 0.5, 15.5 );	
    titlePlot = "2J_Nb_Jets"+channelpostfix;
    histo1D["h_2J_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After at least 2 jets cut: nb. of jets",  16, - 0.5, 15.5 );     
    titlePlot = "1BJ_Nb_Jets"+channelpostfix;
    histo1D["h_1BJ_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After at least 1 CSVL jet cut: nb. of jets",  16, - 0.5, 15.5 );     
    titlePlot = "OSSF_Nb_Jets"+channelpostfix;
    histo1D["h_OSSF_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After OSSF lepton pair req: nb. of jets",  16, - 0.5, 15.5 );     
    titlePlot = "ZMASS_Nb_Jets"+channelpostfix;
    histo1D["h_ZMASS_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After Zmass window: nb. of jets",  16, - 0.5, 15.5 );     

    titlePlot = "cutFlow"+channelpostfix; 
    histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 13,-0.5,12.5);
    titlePlot = "raw_cutFlow"+channelpostfix;
    histo1D["h_raw_cutFlow"] = new TH1F(titlePlot.c_str(), "Raw cutflow", 13,-0.5,12.5);

   // some kinetic variables
    titlePlot = "initial_met"+channelpostfix;
    histo1D["h_initial_met"] = new TH1F(titlePlot.c_str(), "missing E_{T}; E_{T}^{mis} [GeV]", 200, 0,200); 
    titlePlot = "Ht"+channelpostfix;
    histo1D["h_Ht"] = new TH1F(titlePlot.c_str(), "Scalar sum of transverse momenta of the jets; H_{T} [GeV]", 120, 0,1200);
    titlePlot = "Pt_first_Muon"+channelpostfix;
    histo1D["h_Pt_first_Muon"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading muon", 200,0,200); 
    titlePlot = "Pt_2nd_Muon"+channelpostfix;
    histo1D["h_Pt_2nd_Muon"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading muon", 200,0,200);
    titlePlot = "Pt_3d_Muon"+channelpostfix;
    histo1D["h_Pt_3d_Muon"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 3d leading muon", 200,0,200);
    titlePlot = "Pt_first_Electron"+channelpostfix;
    histo1D["h_Pt_first_Electron"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading electron", 200,0,200);
    titlePlot = "Pt_2nd_Electron"+channelpostfix;
    histo1D["h_Pt_2nd_Electron"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading electron", 200,0,200);
    titlePlot = "Pt_3d_Electron"+channelpostfix;
    histo1D["h_Pt_3d_Electron"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 3d leading electron", 200,0,200);
    titlePlot = "Pt_first_Jet"+channelpostfix;
    histo1D["h_Pt_first_Jet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading jet", 200,0,200);
    titlePlot = "Pt_2nd_Jet"+channelpostfix;
    histo1D["h_Pt_2nd_Jet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading jet", 200,0,200);
    titlePlot = "Pt_first_BCSVLJet"+channelpostfix;
    histo1D["h_Pt_first_BCSVLJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading CSVL jet", 200,0,200);
    titlePlot = "Pt_2nd_BCSVLJet"+channelpostfix;
    histo1D["h_Pt_2nd_BCSVLJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading CSVL jet", 200,0,200);
    titlePlot = "Pt_first_BCSVMJet"+channelpostfix;
    histo1D["h_Pt_first_BCSVMJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading CSVM jet", 200,0,200);
    titlePlot = "Pt_2nd_BCSVMJet"+channelpostfix;
    histo1D["h_Pt_2nd_BCSVMJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading CSVM jet", 200,0,200);
    titlePlot = "Pt_first_BCSVTJet"+channelpostfix;
    histo1D["h_Pt_first_BCSVTJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the leading CSVT jet", 200,0,200);
    titlePlot = "Pt_2nd_BCSVTJet"+channelpostfix;
    histo1D["h_Pt_2nd_BCSVTJet"] = new TH1F(titlePlot.c_str(), "Transverse momentum of the 2nd leading CSVT jet", 200,0,200);

    titlePlot = "Zmass"+channelpostfix;
    histo1D["h_Zmass"] = new TH1F(titlePlot.c_str(), "Invariant mass of the Z boson", 200, 0, 200);
    titlePlot = "Zmass_bf"+channelpostfix;
    histo1D["h_Zmass_bf"] = new TH1F(titlePlot.c_str(), "Invariant mass of the Z boson before cuts", 200, 0, 200);
    titlePlot = "mWT"+channelpostfix;
    histo1D["h_mWT"] = new TH1F(titlePlot.c_str(), "Transverse  mass of the W boson", 200, 0, 200);
    titlePlot = "topmass"+channelpostfix;
    histo1D["h_topmass"] = new TH1F(titlePlot.c_str(), "Invariant mass of the SM top (b+l)", 200, 0, 200);

    // plots to check reweighting
    titlePlot = "initial_Nb_CSVLJets_beforeBtagSF"+channelpostfix;
    histo1D["h_initial_Nb_CSVLJets_beforeBtagSF"] = new TH1F(titlePlot.c_str(), "Initial nb. of CSVL jets (no b tag SF) ",  16, - 0.5, 15.5 );
    titlePlot = "initial_Nb_CSVLJets"+channelpostfix;
    histo1D["h_initial_Nb_CSVLJets"] = new TH1F(titlePlot.c_str(), "Initial nb. of CSVL jets ",  16, - 0.5, 15.5 );
    titlePlot = "initial_Nb_Jets_bfCleaning"+channelpostfix;
    histo1D["h_initial_Nb_Jets_bfCleaning"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets (no lepton cleaning) ",  16, - 0.5, 15.5 );     
    titlePlot = "initial_Nb_Jets_unCORJER"+channelpostfix;
    histo1D["h_initial_Nb_Jets_unCORJER"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets (no JER) ",  16, - 0.5, 15.5 );
    titlePlot = "initial_Jet_unCORJER_Pt"+channelpostfix; 
    histo1D["h_initial_Jet_unCORJER_Pt"]  = new TH1F(titlePlot.c_str(), "Initial Pt jet (no JER)",  200, 0., 400 );
    titlePlot = "initial_Nb_Jets_unCORJES"+channelpostfix;
    histo1D["h_initial_Nb_Jets_unCORJES"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets (no JES) ",  16, - 0.5, 15.5 );
    titlePlot = "initial_Jet_unCORJES_Pt"+channelpostfix;
    histo1D["h_initial_Jet_unCORJES_Pt"]  = new TH1F(titlePlot.c_str(), "Initial Pt jet (no JES)",  200, 0., 400 );
    titlePlot = "initial_Jet_Pt"+channelpostfix;
    histo1D["h_initial_Jet_Pt"]  = new TH1F(titlePlot.c_str(), "Initial Pt jet ",  200, 0., 400 );
    titlePlot = "initial_nPV_beforePUSF" + channelpostfix;
    histo1D["h_initial_nPV_beforePUSF"] = new TH1F(titlePlot.c_str(),"Initial nb of vertices before PU reweighting", 13,-0.5,12.5 );
    titlePlot = "initial_nPV_afterPUSF" + channelpostfix;
    histo1D["h_initial_nPV_afterPUSF"] = new TH1F(titlePlot.c_str(),"Initial nb of vertices after PU reweighting", 13,-0.5,12.5 );
    titlePlot = "muonSF"+channelpostfix;
    histo1D["h_muonSF"] = new TH1F(titlePlot.c_str(), "Muon scale factors", 100, 0, 1);
    titlePlot = "electronSF"+channelpostfix;
    histo1D["h_electronSF"] = new TH1F(titlePlot.c_str(), "Electron scale factors", 100, 0, 1);
    titlePlot = "puSF"+channelpostfix;
    histo1D["h_puSF"] = new TH1F(titlePlot.c_str(), "PU scale factors", 100, 0, 1);
    titlePlot = "btagSF"+channelpostfix;
    histo1D["h_btagSF"] = new TH1F(titlePlot.c_str(), "Btag scale factors", 100, 0, 1);
        

 
    titlePlot = "GmuonSF"+channelpostfix;
    histo2D["h2_muonSF"]= new TH2F(titlePlot.c_str(), "Muon scale factors in function of p_{T} and #eta; p_{T} [GeV]; #eta", 60, 0, 600, 21, 0, 2.1);
    titlePlot = "GelectronSF"+channelpostfix;
    histo2D["h2_electronSF"]= new TH2F(titlePlot.c_str(), "Electron scale factors in function of p_{T} and #eta; p_{T} [GeV]; #eta", 60, 0, 600, 21, 0, 2.1); 
   
    
    
    // make root tree file name
    string roottreename = "/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/FCNCAnalysis/";
    roottreename+= "Ntuples/";

    mkdir((roottreename).c_str(),0777);
    roottreename+= datasets[d]->Name();
    roottreename+="_";
    roottreename+= strJobNum;
    if(doJESup) roottreename+= "_JESup";
    if(doJESdown) roottreename+= "_JESdown";
    if(doJERup) roottreename+="_JERup";
    if(doJERdown) roottreename+="_JERdown"; 
    roottreename+="_tree.root";

    cout << "  - Recreate outputfile for ntuples ... " << roottreename.c_str() << endl; 
    // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
    TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
    fileout->cd();
   
    ////////////////////
    //  SF
    //  ////////////////
    if(applyBtagSF && !isdata){
       // http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v2.pdf
       // removed in the first line of the csv file "CSVv2;"  
       bTagCalib = new BTagCalibration("CSVv2",pathToCaliDir+"BTagging/CSVv2_13TeV_25ns.csv");
       bTagReader = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_LOOSE,"comb","central");
       if(fillingbTagHistos) btwt = new BTagWeightTools(bTagReader,30,999,2.4,"forBtagSF/HistosPtEta_"+dataSetName+".root");
       else btwt = new BTagWeightTools(bTagReader,30,999,2.4,"forBtagSF/HistosPtEta_tmp.root");
    }
    
    

   
    

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

     double metPt; 
     double metPx; 
     double metPy; 
 
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
    myTree -> Branch("metPt", &metPt, "metPt/D");
    myTree -> Branch("metPx", &metPx, "metPx/D");
    myTree -> Branch("metPy", &metPy, "metPy/D");
     
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

    //////////////////////////
    // Initialize JEC  Factor
    ///////////////(order matters! )
     vector<JetCorrectorParameters> vCorrParam;

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
   {
     JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
       vCorrParam.push_back(*L3JetCorPar);
       JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
       vCorrParam.push_back(*L2L3ResJetCorPar);
   }
   else
   {
       JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
       vCorrParam.push_back(*L1JetCorPar);
       JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt");
       vCorrParam.push_back(*L2JetCorPar);
       JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
        vCorrParam.push_back(*L3JetCorPar);
   }
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV2_MC_Uncertainty_AK4PFchs.txt");

 
   JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
     
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
    bool trigged = false; 
    int trigEMU, trigMUMU,trigEE; 
    trigEMU = trigMUMU = trigEE = -1; 
    if (verbose > 1)
      cout << "	Loop over events " << endl;

    vector < TRootVertex* > vertex;
    vector < TRootMuon* > init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* > init_jets;
    vector < TRootJet* > init_jets_unCORJER;
    vector < TRootJet* > init_jets_unCORJES;
    vector < TRootMET* > mets;
    int currentRun; 
    vector<TRootGenJet*> genjets;

    bool isGoodPV;
    vector<TRootMuon*> selectedMuons;
    vector<TRootMuon*> selectedVetoMuons;

    vector<TRootPFJet*> selectedJets_bfCleaning;
    vector<TRootPFJet*> selectedJets ;
    vector<TRootPFJet*> selectedJets_unCORJER ;
    vector<TRootPFJet*> selectedJets_unCORJES ;
    vector<TRootElectron*> selectedElectrons ;
    vector<TRootElectron*> selectedVetoElectrons ;

    vector<bool> BtagBooleans; 
    vector<TRootJet*> selectedBCSVLJets; 
    vector<TRootJet*> selectedBCSVMJets;
    vector<TRootJet*> selectedBCSVTJets;

    std::vector<int> Leptons;
    TRootJet* tempJet;
    
    TLorentzVector Zboson;
    TLorentzVector Wlep;
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    {
      
      vertex.clear();
      init_muons.clear();
      init_electrons.clear(); 
      init_jets.clear(); 
      init_jets_unCORJER.clear(); 
      init_jets_unCORJES.clear();
      mets.clear();
      
      currentRun = -99999;
      isGoodPV=false;
      selectedMuons.clear();
      selectedVetoMuons.clear();
      selectedJets_bfCleaning.clear();
      selectedJets.clear() ;
      selectedJets_unCORJER.clear() ;
      selectedJets_unCORJES.clear() ; 
      selectedElectrons.clear() ;
      selectedVetoElectrons.clear() ;

      BtagBooleans.clear();
      selectedBCSVLJets.clear();
      selectedBCSVMJets.clear();
      selectedBCSVTJets.clear();
      
      Zboson.Clear();
      Wlep.Clear(); 

      nEvents[d]++;
      
      if (ievt%500 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);


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
      
      
      currentRun = event->runId();  
      run_num=event->runId();
      evt_num=event->eventId();
      lumi_num=event->lumiBlockId();
      nvtx=vertex.size();
      npu=(int)event->nTruePU();
      if( run_num > 10000){//data
         isdata=1;
      }
      bookkeeping->Fill();

      //////////////////////
      // Load genjets for JER smearing 
      /////////////////////////////////
      genjets.clear();
      if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
       {
          // loading GenJets as I need them for JER
          genjets = treeLoader.LoadGenJet(ievt);
       } 
      //JER
     
     if(doJERdown){
         init_jets_unCORJER = init_jets; 
        jetTools->correctJetJER(init_jets,genjets,mets[0], "minus", false); // false means don't use old numbers, but new ones
      }
      else if(doJERup){
         init_jets_unCORJER= init_jets;
        jetTools->correctJetJER(init_jets,genjets,mets[0], "plus", false); // false means don't use old numbers, but new ones
      }
      else{
        init_jets_unCORJER= init_jets;
	jetTools->correctJetJER(init_jets,genjets,mets[0], "nominal", false); // false means don't use old numbers, but new ones
     }

      // JES
      if(doJESup){
         init_jets_unCORJES = init_jets; 
         jetTools->correctJetJESUnc(init_jets, "plus", 1);


     }
     else if(doJESdown){
         init_jets_unCORJES = init_jets;
         jetTools->correctJetJESUnc(init_jets, "minus", 1);


     }
     else init_jets_unCORJES = init_jets;
     
      /////////////////////////////
      /// Trigger
      ///////////////////////////
      
       trigEMU = trigMUMU = trigEE = -1; 
       itrigger = -1; 

      //If the HLT is applied 
      if(runHLT && previousRun != currentRun){
        //The HLT is only used for data
        if(isdata == 1){
          //The HLT path is dependent of the mode, these paths are the several steps or software modules. Each module performs a well defined task 
          // such as reconstruction of physics objects, making intermediate decisions, triggering more refined reconstructions in subsequent modules, 
          // or calculating the final decision for that trigger path.
          trigEMU = treeLoader.iTrigger (string ("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"), currentRun, iFile);
          trigEE = treeLoader.iTrigger (string("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2"), currentRun, iFile);
          trigMUMU =treeLoader.iTrigger (string ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2"), currentRun, iFile); 
	  if(emu) itrigger = trigEMU;
          else if(mumu && !trigEMU) itrigger = trigMUMU;
          else if(ee && !trigEMU) itrigger = trigEE;

          if(itrigger == 9999) 
	  {
	      cout << "ERROR: no valid trigger found for this event/data in run " << event->runId() << endl; 
	  }   
	  
        } // closing the HLT for data loop
        //For the MC, there is no triggerpath
        else
        {
          trigEMU =  itrigger = treeLoader.iTrigger ("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1", currentRun, iFile);
          trigMUMU =  itrigger = treeLoader.iTrigger (string ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1"), currentRun, iFile);
          trigEE = itrigger = treeLoader.iTrigger (string ("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1"), currentRun, iFile);
          if(emu) itrigger = trigEMU;
          else if(mumu && !trigEMU) itrigger = trigMUMU;
          else if(ee && !trigEMU) itrigger = trigEE;

        }
      } // closing the HLT run loop
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets, init_muons, init_electrons, mets);
      Run2Selection selection_unCORJER( init_jets_unCORJER, init_muons, init_electrons, mets);    
      Run2Selection selection_unCORJES( init_jets_unCORJES, init_muons, init_electrons, mets);

      bool isGoodPV = selection.isPVSelected(vertex, 4, 24,2); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
      selectedMuons = selection.GetSelectedMuons(20, 2.4, 0.15,"Tight","Spring15");  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
      selectedVetoMuons = selection.GetSelectedMuons(8, 2.4, 0.15,"Loose","Spring15");  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)

      selectedJets_bfCleaning = selection.GetSelectedJets(30, 2.4, true, "Tight"); 
      selectedJets = selection.GetSelectedJets(30, 2.4, true, "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
      selectedJets_unCORJER = selection_unCORJER.GetSelectedJets(30, 2.4, true, "Tight"); 
      selectedJets_unCORJES = selection_unCORJES.GetSelectedJets(30, 2.4, true, "Tight");
      selectedElectrons = selection.GetSelectedElectrons(20, 2.5, "Tight", "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      selectedVetoElectrons = selection.GetSelectedElectrons(12, 2.5, "Veto", "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)

      sort(selectedJets.begin(), selectedJets.end(),HighestPt());
      sort(selectedJets_bfCleaning.begin(), selectedJets_bfCleaning.end(),HighestPt()); 
      sort(selectedJets_unCORJER.begin(), selectedJets_unCORJER.end(),HighestPt());
      sort(selectedJets_unCORJES.begin(), selectedJets_unCORJES.end(),HighestPt());
      sort(selectedMuons.begin(), selectedMuons.end(), HighestPt()); 
      sort(selectedVetoMuons.begin(), selectedVetoMuons.end(), HighestPt());
      sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt()); 
      sort(selectedVetoElectrons.begin(), selectedVetoElectrons.end(), HighestPt());
      
      //jetcleaning
      if(applyJetCleaning){
       for (int origJets=0; origJets<selectedJets.size(); origJets++){
         bool erased = false;
         if(selectedMuons.size()>0){
             if(selectedJets[origJets]->DeltaR(*selectedMuons[0])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }
         if(selectedMuons.size()>1 && !erased){
             if(selectedJets[origJets]->DeltaR(*selectedMuons[1])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }            
         if(selectedMuons.size()>2 && !erased){
             if(selectedJets[origJets]->DeltaR(*selectedMuons[2])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }
         if(selectedElectrons.size()>0 && !erased){        
             if(selectedJets[origJets]->DeltaR(*selectedElectrons[0])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }                       
         if(selectedElectrons.size()>1 && !erased){
             if(selectedJets[origJets]->DeltaR(*selectedElectrons[1])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }
         if(selectedElectrons.size()>2 && !erased){
             if(selectedJets[origJets]->DeltaR(*selectedElectrons[2])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
         }
       }
       if(verbose>3) if( selectedJets_bfCleaning.size() != selectedJets.size()) cout << "original = " << selectedJets_bfCleaning.size() << " after cleaning = " << selectedJets.size() << endl; 
      }

      
      for(unsigned int i = 0; i < selectedJets.size() ; i++)
      {
         bool Btagged = false;
         tempJet = (TRootJet*) selectedJets[i];
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
	
       double met_px = mets[0]->Px();
       double met_py = mets[0]->Py();
       double met_pt = sqrt(met_px*met_px + met_py*met_py);

      ////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      // PU SF
       if (applyPUSF && !isdata ){
          double puWeight = LumiWeights.ITweight( nvtx ); // simplest reweighting, just use reconstructed number of PV. faco
          puScaleFactor=puWeight;
          if (verbose>3) cout << "puScaleFactor is " << puScaleFactor << endl;
      }

      // Lepton SF
      float muon1SF, muon2SF,muon3SF, electron1SF, electron2SF, electron3SF; 
      muon1SF =  muon2SF = muon3SF = electron1SF = electron2SF = electron3SF = 0.;
       if(applyMuonSF && !isdata){
         if(selectedMuons.size()>0) {muon1SF = muonSFWeight->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0); muonScaleFactor = muon1SF; 
             histo2D["h2_muonSF"]->Fill(selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), muon1SF);}
         if(selectedMuons.size()>1) {muon2SF = muonSFWeight->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0); muonScaleFactor *= muon2SF; 
             histo2D["h2_muonSF"]->Fill(selectedMuons[1]->Pt(), selectedMuons[1]->Eta(), muon2SF);}
         if(selectedMuons.size()>2) {muon3SF = muonSFWeight->at(selectedMuons[2]->Eta(), selectedMuons[2]->Pt(), 0); muonScaleFactor *= muon3SF; 
             histo2D["h2_muonSF"]->Fill(selectedMuons[2]->Pt(), selectedMuons[2]->Eta(), muon3SF);}
       }
       if(applyElectronSF && !isdata){
         if(selectedElectrons.size()>0)  {electron1SF =  electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0); electronScaleFactor = electron1SF;  
             histo2D["h2_electronSF"]->Fill(selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), electron1SF);}
         if(selectedElectrons.size()>1)  {electron2SF =  electronSFWeight->at(selectedElectrons[1]->Eta(),selectedElectrons[1]->Pt(),0); electronScaleFactor *= electron2SF;  
             histo2D["h2_electronSF"]->Fill(selectedElectrons[1]->Pt(), selectedElectrons[1]->Eta(), electron2SF); }
         if(selectedElectrons.size()>2)  {electron3SF =  electronSFWeight->at(selectedElectrons[2]->Eta(),selectedElectrons[2]->Pt(),0); electronScaleFactor *= electron3SF;  
             histo2D["h2_electronSF"]->Fill(selectedElectrons[2]->Pt(), selectedElectrons[2]->Eta(), electron3SF); }
       }
      if(fillingbTagHistos){
         if(!isdata && applyBtagSF) btwt->FillMCEfficiencyHistos(selectedJets); 
                
      }
      if (verbose>3) cout<<"getMCEventWeight for btag"<<endl;
      if(applyBtagSF && !isdata && !fillingbTagHistos){
           btagScaleFactor =  btwt->getMCEventWeight(selectedJets,(TFile*) "HistosPtEta.root", false);
           // cout<<"btag weight "<<btagWeight<<endl;
       }     
         
       
      histo1D["h_muonSF"]->Fill(muonScaleFactor); 
      histo1D["h_electronSF"]->Fill(electronScaleFactor);
      histo1D["h_puSF"]->Fill(puScaleFactor);
      histo1D["h_btagSF"]->Fill(btagScaleFactor);

      if(applyMuonSF) scaleFactor *= muonScaleFactor;
      if(applyElectronSF) scaleFactor *= electronScaleFactor;
      if(applyPUSF) scaleFactor *= puScaleFactor;
      if(applyBtagSF) scaleFactor *= btagScaleFactor;
      if(!applyGlobalSF || isdata) scaleFactor = 1.;

      
      if(applyPUSF) { histo1D["h_initial_nPV_beforePUSF"]->Fill(vertex.size(), scaleFactor*lumiWeight/puScaleFactor);}
      else { histo1D["h_initial_nPV_beforePUSF"]->Fill(vertex.size(), scaleFactor*lumiWeight);}
      histo1D["h_initial_nPV_afterPUSF"]->Fill(vertex.size(), scaleFactor*lumiWeight);
      if(applyBtagSF) {histo1D["h_initial_Nb_CSVLJets_beforeBtagSF"]->Fill(vertex.size(), scaleFactor*lumiWeight/btagScaleFactor);}
      else { histo1D["h_initial_Nb_CSVLJets_beforeBtagSF"]->Fill(vertex.size(), scaleFactor*lumiWeight);}
      histo1D["h_initial_Nb_CSVLJets"]->Fill(vertex.size(), scaleFactor*lumiWeight);
      
      // Start analysis selection
      eventSelected = false;

      
      
      /// Initial nbrs
      
      histo1D["h_cutFlow"]->Fill(0., scaleFactor*lumiWeight);
      histo1D["h_raw_cutFlow"]->Fill(0.);

      /// Trigger
      if(runHLT) trigged = treeLoader.EventTrigged(itrigger);
      else trigged = true; 
      //SELECTION 
      if(trigged)
      { 
       histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
       histo1D["h_raw_cutFlow"]->Fill(1.);

       if (isGoodPV)
       {
        if(verbose>3) cout << "GoodPV" << endl; 
        histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
        histo1D["h_raw_cutFlow"]->Fill(2.);

        histo1D["h_initial_Nb_Jets_bfCleaning"]->Fill(selectedJets_bfCleaning.size(), scaleFactor*lumiWeight); 
        histo1D["h_initial_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);
        histo1D["h_initial_Nb_Jets_unCORJER"]->Fill(selectedJets_unCORJER.size(), scaleFactor*lumiWeight);
        histo1D["h_initial_Nb_Jets_unCORJES"]->Fill(selectedJets_unCORJES.size(), scaleFactor*lumiWeight);     
        histo1D["h_initial_met"]->Fill(met_pt, scaleFactor*lumiWeight);
	

        float tempPt = 0.;           
        for(unsigned int i = 0 ; i < selectedJets.size(); i++){
            TRootJet* jet = (TRootJet*) selectedJets[i];
            histo1D["h_initial_Jet_Pt"]->Fill(jet->Pt(),  scaleFactor*lumiWeight);
            tempPt += jet->Pt();
        }
        histo1D["h_Ht"]->Fill(tempPt, scaleFactor*lumiWeight);
        
        for(unsigned int i = 0 ; i < selectedJets_unCORJER.size(); i++){
          TRootJet* jet_unCORJER = (TRootJet*) selectedJets_unCORJER[i];
          histo1D["h_initial_Jet_unCORJER_Pt"]->Fill(jet_unCORJER->Pt(),  scaleFactor*lumiWeight);
        
        }
        for(unsigned int i = 0 ; i < selectedJets_unCORJES.size(); i++){
          TRootJet* jet_unCORJES = (TRootJet*) selectedJets_unCORJES[i];
          histo1D["h_initial_Jet_unCORJES_Pt"]->Fill(jet_unCORJES->Pt(),  scaleFactor*lumiWeight);

        }
        
        if(selectedMuons.size() > 0) histo1D["h_Pt_first_Muon"]->Fill(selectedMuons[0]->Pt(),	scaleFactor*lumiWeight);
        if(selectedMuons.size() > 1) histo1D["h_Pt_2nd_Muon"]->Fill(selectedMuons[1]->Pt(),   scaleFactor*lumiWeight);
        if(selectedMuons.size() > 2) histo1D["h_Pt_3d_Muon"]->Fill(selectedMuons[2]->Pt(),   scaleFactor*lumiWeight);
        if(selectedElectrons.size() > 0) histo1D["h_Pt_first_Electron"]->Fill(selectedElectrons[0]->Pt(),   scaleFactor*lumiWeight);
        if(selectedElectrons.size() > 1) histo1D["h_Pt_2nd_Electron"]->Fill(selectedElectrons[1]->Pt(),   scaleFactor*lumiWeight);
        if(selectedElectrons.size() > 2) histo1D["h_Pt_first_Electron"]->Fill(selectedElectrons[2]->Pt(),  scaleFactor*lumiWeight);
        if(selectedJets.size()>0) histo1D["h_Pt_first_Jet"]->Fill(selectedJets[0]->Pt(), scaleFactor*lumiWeight);
        if(selectedJets.size()>1) histo1D["h_Pt_2nd_Jet"]->Fill(selectedJets[1]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVLJets.size()>0) histo1D["h_Pt_first_BCSVLJet"]->Fill(selectedBCSVLJets[0]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVLJets.size()>1) histo1D["h_Pt_2nd_BCSVLJet"]->Fill(selectedBCSVLJets[1]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVMJets.size()>0) histo1D["h_Pt_first_BCSVMJet"]->Fill(selectedBCSVMJets[0]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVMJets.size()>1) histo1D["h_Pt_2nd_BCSVMJet"]->Fill(selectedBCSVMJets[1]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVTJets.size()>0) histo1D["h_Pt_first_BCSVTJet"]->Fill(selectedBCSVTJets[0]->Pt(), scaleFactor*lumiWeight);
        if(selectedBCSVTJets.size()>1) histo1D["h_Pt_2nd_BCSVTJet"]->Fill(selectedBCSVTJets[1]->Pt(), scaleFactor*lumiWeight);
	  


        if (selectedMuons.size() + selectedElectrons.size()== 3 )
        {
	    
            if(verbose>3) cout << "3 electrons "<< endl; 
            histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
            histo1D["h_raw_cutFlow"]->Fill(3.);
            histo1D["h_3L_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);
	    
	    
	    if (selectedVetoMuons.size()+ selectedVetoElectrons.size() == 3)
	    {
               histo1D["h_cutFlow"]->Fill(4., scaleFactor*lumiWeight);
               histo1D["h_raw_cutFlow"]->Fill(4.);
	    
	       Leptons.clear(); 
	       Leptons = OSSFLeptonPairCalculator(selectedElectrons, selectedMuons, verbose);
	       if(verbose>3) cout <<   Leptons[0]<< " , " <<  Leptons[1]<< " , " <<  Leptons[2]<< " , " <<  Leptons[3]<< " , " <<  Leptons[4]<< " , " <<  Leptons[5]   << endl; 
	
	       bool OSSFpair = false; 
	       if( (Leptons[0] != -5 && Leptons[1] != -5) | (Leptons[3] != -5 && Leptons[4] != -5) ) OSSFpair = true; 
	       if(OSSFpair)
	       {
	          if(verbose>3) cout << " OSSF "<< endl; 
 	          histo1D["h_OSSF_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);
	          
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
                  histo1D["h_Zmass"]->Fill(Zboson.M(), scaleFactor*lumiWeight);
 
	          if(Zboson.M() < 106 && Zboson.M() > 76) ZmassWindow = true; 
	          if(ZmassWindow)
	          {
	            if(verbose>3) cout << " Zmass window " << endl; 
	            histo1D["h_cutFlow"]->Fill(5, scaleFactor*lumiWeight);
                    histo1D["h_raw_cutFlow"]->Fill(5);
	            histo1D["h_ZMASS_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);

		    if(selectedJets.size() > 1)
		    {
		       if(verbose>3) cout << " at least 2 jets " << endl; 
                       histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
                       histo1D["h_raw_cutFlow"]->Fill(6.);
                       histo1D["h_2J_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);

                       if(selectedBCSVLJets.size()>0)
	               {
	       		  if(verbose>3) cout << " at least 1 bjet " << endl; 
			  histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                          histo1D["h_raw_cutFlow"]->Fill(7.);
               		  histo1D["h_1BJ_Nb_Jets"]->Fill(selectedJets.size(), scaleFactor*lumiWeight);

                          float Phi_Wlep_MET = mets[0]->DeltaPhi(Wlep);
                 	  float CosPhi_Wlep_MET = cos(Phi_Wlep_MET);
                 	  float mWT = TMath::Sqrt(2*met_pt*Wlep.Pt()*(1-CosPhi_Wlep_MET));
                          histo1D["h_mWT"]->Fill(mWT, scaleFactor*lumiWeight);

                          if(mWT > 20)
			  {
                             histo1D["h_cutFlow"]->Fill(8., scaleFactor*lumiWeight);
                             histo1D["h_raw_cutFlow"]->Fill(8.);

                             TLorentzVector Bjet;
                             Bjet.Clear();
                             Bjet.SetPxPyPzE(selectedBCSVLJets[0]->Px(),selectedBCSVLJets[0]->Py(),selectedBCSVLJets[0]->Pz(),selectedBCSVLJets[0]->Energy());

			     TLorentzVector SMtop; 
			     SMtop.Clear(); 
			     SMtop = Bjet + Wlep;     
			     float topmass = SMtop.M();
                             histo1D["h_topmass"]->Fill(topmass, scaleFactor*lumiWeight);

			     if( topmass < 155 && topmass > 95)
                             {
                                histo1D["h_cutFlow"]->Fill(9., scaleFactor*lumiWeight);
                                histo1D["h_raw_cutFlow"]->Fill(9.);
 
			        eventSelected = true;
                            } // topmass
			  } //mWt
		       } // >0 bjets

	            }// > 1 jet
	          }//Zmass
	       } //OSSFpair	  
	    
	    } //lepton veto
	   
         
        }  // 3 leptons
       }  // good PV
      } // trigger
      
      
      if (! eventSelected)
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


       metPt = met_pt; 
       metPx = met_px; 
       metPy = met_py;      
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
        tempJet = (TRootJet*) selectedJets[i];
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
        tempJet = (TRootJet*) selectedBCSVLJets[i];
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
        tempJet = (TRootJet*) selectedBCSVMJets[i];
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
        tempJet = (TRootJet*) selectedBCSVTJets[i];
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
    
    
    myTree->Write();

    fileout->Write();
    fileout->Close();
    delete fileout; 

   
   
   
   ///*****************///
  ///   Write plots   ///
  ///*****************///
  string pathPNG = "ControlPlots/";
  mkdir(pathPNG.c_str(),0777);
  mkdir((pathPNG+"1DPlot/").c_str(),0777); // 0777 if it doesn't exist already, make it
  if(doJERup) mkdir((pathPNG+"1DPlot/JERup").c_str(),0777);
  if(doJERdown) mkdir((pathPNG+"1DPlot/JERdown").c_str(),0777);
  if(doJESup) mkdir((pathPNG+"1DPlot/JESup").c_str(),0777);
  if(doJESdown) mkdir((pathPNG+"1DPlot/JESdown").c_str(),0777);

  mkdir((pathPNG+"2DPlot/").c_str(),0777); // 0777 if it doesn't exist already, make it
  if(doJERup) mkdir((pathPNG+"2DPlot/JERup").c_str(),0777);
  if(doJERdown) mkdir((pathPNG+"2DPlot/JERdown").c_str(),0777);
  if(doJESup) mkdir((pathPNG+"2DPlot/JESup").c_str(),0777);
  if(doJESdown) mkdir((pathPNG+"2DPlot/JESdown").c_str(),0777);

  ///Write histograms
  fout->cd();
  for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    cout << "1D Plot: " << it->first << endl;
   // TCanvas ctemp = 
    
    TH1F *temp = it->second;
    string name = it->first;
    cout << name << endl; 
    temp->Draw();  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
    //c1->SaveAs(pathPNG+"1DPlot/"+name modeString[mode] + "_" + cutLabel + ".png");
    //temp->Write(fout, name, true, pathPNG+"1DPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
  }
  for (map<string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
    if(verbose>3) cout << "2D Plot: " << it->first << endl;
   
     TH2F *temp = it->second;
     string name = it->first;
     temp->Draw();
  }
   
   
   ///////////////////
   /// CLEANING
   /////////////////
  fout->Write();   
  fout->Close();
  
  delete fout;

   treeLoader.UnLoadDataset();
 } // end datasetloop


 
// delete tcdatasets;
  delete tcAnaEnv;
 // delete configTree;
 
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
