//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for FCNC tZq                     ////
//////////////////////////////////////////////////////////////////////////////


#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <iostream>
#include <map>
#include <cstdlib>

//user code

#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"

#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Tools/interface/SourceDate.h"
#include "TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h

//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"

#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"




using namespace std;
using namespace TopTree;
using namespace reweight;


/// TH1F
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
//Initializing CSVv2 b-tag WP
float workingpointvalue_Loose = -1;
float workingpointvalue_Medium = -1;
float workingpointvalue_Tight = -1;


//What you want to do
bool synchex = false;
bool Assigned = false;


// home made functions
int FCNCjetCalculator(std::vector<TRootPFJet*> Jets, TLorentzVector recoZ ,int index, int verb);
int FCNCjetCalculatorTagger(std::vector<TRootPFJet*> Jets, int index, int verb);
int SMjetCalculator(std::vector<TRootPFJet*> Jets,int verb);
double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy);
float EffectiveAreaRho(TRootElectron *el, float _rho) ;
float EffectiveArea(TRootElectron *el) ;
float relPfIsoEl(TRootElectron *el, float _rho);
float IsoDBeta(TRootMuon *mu);
vector<TLorentzVector> LeptonAssigner(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons);
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met);
vector< pair<unsigned int, unsigned int> > JetPartonPair_charm;
vector< pair<unsigned int, unsigned int> > JetPartonPair_bottom;
vector< pair<unsigned int, unsigned int> > JetPartonPair_electron;
vector< pair<unsigned int, unsigned int> > JetPartonPair_muon;
vector< pair<unsigned int, unsigned int> > JetPartonPair_Welectron;
vector< pair<unsigned int, unsigned int> > JetPartonPair_Wmuon;
// administration functions
string ConvertIntToString(int Number, bool pad);
string MakeTimeStamp();


// members
//   bool stop_program;
double M_W  = 80.4;
double M_mu =  0.10566; // 105.66 MeV/c^2
double M_el = 0.000510999; // 0.510998910 Mev/c^2
int nMatched_charm = 0;
int nNonMatched_charm = 0;
int nMatched_charm_tag = 0;
int nNonMatched_charm_tag = 0;
int nMatched_bottom = 0;
int nNonMatched_bottom = 0;
int nMatched_Zelec = 0;
int nNonMatched_Zelec = 0;
int nMatched_Zmu = 0;
int nNonMatched_Zmu = 0;
int nMatched_Welec = 0;
int nNonMatched_Welec = 0;
int nMatched_Wmu = 0;
int nNonMatched_Wmu = 0;
int nTagEqMass = 0;
int nTagNotEqMass = 0;
bool matching = false;
bool elecbool = false;
bool mubool = false;
vector <int> muIndices;
vector <int> elecIndices;
vector <int> WmuIndices;
vector <int> WelecIndices;


int main (int argc, char *argv[])
{
  
  string dateString = MakeTimeStamp();
  cout << "***********************************" << endl;
  cout << "***   Beginning of program: tZq FCNC      ***" << endl;
  cout << "***********************************" << endl;
  cout << "Current time: " << dateString << endl;
  
  clock_t start = clock();
  
  
  
  ///////////////////////////
  /// Configuration      ///
  //////////////////////////
  int verbose = 1; // 0 = cout alll
  bool eventSelected = false;
  bool baseSelected = false;
  int nbTrig = 0;
  int nbBaseline = 0;
  int nbGPV = 0;
  int nbSelectedEvents = 0;
  int nbEvents = 0;
  double dataLumi = 0; //pb
  bool runHLT = true;
  bool applyJetLeptonCleaning = true;
  bool fillBtagHisto = false;
  bool printTrigger = false;
  bool printLeptonSF = false;
  bool applyJER = false;
  bool applyJES = false;
  bool applyNegWeightCorrection = false;
  bool applyPU = false;
  bool applyLeptonSF = false;
  bool btagShape = true;
  string xmlFileName = "";
  int maxMCParticles = -1;
  
  //////////////////////////////////////////////
  /// Set up everything for local submission ////
  ///////////////////////////////////////////////
  // check the arguments passed
  if(verbose == 0)
  {
    cout << " The list of arguments are: " << endl;
    for (int n_arg=1; n_arg<argc; n_arg++)
    {
      std:: cerr << "  - arg number " << n_arg << " is " << argv[n_arg] << std::endl;
    }
  }
  if(argc < 19)
  {
    std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
    for (int n_arg=1; n_arg<argc; n_arg++)
    {
      std:: cerr << "  - arg number " << n_arg << " is " << argv[n_arg] << std::endl;
    }
    exit(2);
    
  }
  
  // Put the argument in a format we can use
  const string dName             = argv[1];
  const string dTitle		 = argv[2];
  const int color		  = strtol(argv[4], NULL, 10);
  const int ls  		  = strtol(argv[5], NULL, 10);
  const int lw  		  = strtol(argv[6], NULL, 10);
  const float normf		   = strtod(argv[7], NULL);
  const float EqLumi		  = strtod(argv[8], NULL);
  const float xSect		  = strtod(argv[9], NULL);
  const float PreselEff 	  = strtod(argv[10], NULL);
  string fileName		  = argv[11];
  // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list)
  const int JES                 =  strtol(argv[argc-6], NULL,10);
  const int JER                 =  strtol(argv[argc-5], NULL,10);
  const int FillBtagHisto	 =  strtol(argv[argc-4], NULL,10);
  const int JobNum		  = strtol(argv[argc-3], NULL, 10);
  const int startEvent  	  = strtol(argv[argc-2], NULL, 10);
  const int endEvent		  = strtol(argv[argc-1], NULL, 10);
  
  applyJES = JES;
  applyJER = JER;
  fillBtagHisto = FillBtagHisto;
  // all the files are stored from arg 11 to argc-2
  vector<string> vecfileNames;
  for(int args = 11; args < argc-7; args++)
  {
    vecfileNames.push_back(argv[args]);
  }
  
  if (verbose==0)
  {
    cout << "The list of file to run over will be printed..." << endl;
    for ( int nfiles = 0; nfiles < vecfileNames.size(); nfiles++)
    {
      cout << "file number " << nfiles << " is " << vecfileNames[nfiles] << endl;
    }
  }
  /// define channels
  //
  
    cout << " --> Using the all  channel <-- " << endl;
    xmlFileName = "config/Run2TriLepton.xml" ;
    dataLumi = 2700; //pb
    
  
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
  bool isData= false;
  if(dName.find("Data")!=string::npos || dName.find("data")!=string::npos || dName.find("DATA")!=string::npos){
    isData = true;
    cout << "running on data !!!!" << endl;
    cout << "luminosity is " << dataLumi << endl;
  }
  cout << "----------------------------------------" << endl;
  
  
  /////////////////////////////////
  //  Set up AnalysisEnvironment
  /////////////////////////////////
  
  AnalysisEnvironment anaEnv;
  cout<<" - Creating environment ..."<<endl;
  anaEnv.PrimaryVertexCollection = "PrimaryVertex";
  anaEnv.JetCollection = "PFJets_slimmedJets";
  anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
  anaEnv.METCollection = "PFMET_slimmedMETs";
  anaEnv.MuonCollection = "Muons_slimmedMuons";
  anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
  anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
  //  anaEnv.TrackMETCollection = "";
  //  anaEnv.GenEventCollection = "GenEvent";
  anaEnv.NPGenEventCollection = "NPGenEvent";
  anaEnv.MCParticlesCollection = "MCParticles";
  anaEnv.loadFatJetCollection = false;
  anaEnv.loadGenJetCollection = true;
  //  anaEnv.loadGenEventCollection = false;
  anaEnv.loadNPGenEventCollection = false;
  anaEnv.loadMCParticles = true;
  //  anaEnv.loadTrackMETCollection = false;
  anaEnv.JetType = 2;
  anaEnv.METType = 2;
  
  ////////////////////////////////
  //  Load datasets
  ////////////////////////////////
  
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
  theDataset->SetEquivalentLuminosity(EqLumi);
  datasets.push_back(theDataset);
  int ndatasets = datasets.size() - 1 ;
  
  ////////////////////////////
  ///  Initialise trigger  ///
  ////////////////////////////
  
  if(verbose == 0) cout << "Initializing trigger" << endl;
  Trigger* trigger_mumu  = new Trigger(1, 0, 0, 1,0); // mu , el, single, double, tri
  Trigger* trigger_ee  = new Trigger(0, 1, 0, 1,0);
  Trigger* trigger_emu  = new Trigger(1, 1, 0, 1,0) ;
  Trigger* trigger_mumumu  = new Trigger(1, 0, 0, 0,1);
  Trigger* trigger_eee  = new Trigger(0, 1, 0, 0,1);
  Trigger* trigger_emumu_mumue  = new Trigger(1, 1, 0, 0,1) ;
  Trigger* trigger_mu  = new Trigger(1, 0, 1, 0,0);
  Trigger* trigger_e  = new Trigger(0, 1, 1, 0,0);
  
  
  ////////////////////////
  // intialize  Calibrations      //
  ///////////////////////
  BTagCalibration *btagcalib;
  BTagCalibrationReader *btagreader;
  BTagWeightTools *btwt;
  BTagCalibrationReader * reader_csvv2;
  // for pu
  LumiReWeighting LumiWeights;
  
  // JER / JEC
  vector<JetCorrectorParameters> vCorrParam;
  string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV6/";
  
  
  ///////////////////////////////
  //  Set up Output ROOT file  ///
  //////////////////////////////
  stringstream ss;
  ss << JobNum;
  string strJobNum = ss.str();
  string histo_dir = "NtupleMakerOutput/TriLepton_histos";
  string histo_dir_date = histo_dir+"/TriLepton_histos_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
  mkdir(histo_dir_date.c_str(),0777);
  
  string rootFileName (histo_dir_date+"/FCNC_3L_"+dName+".root");
  if (strJobNum != "0")
  {
    if(verbose == 0) cout << "strJobNum is " << strJobNum << endl;
    rootFileName = histo_dir_date+"/FCNC_3L_"+dName + "_"+strJobNum+".root";
  }
  cout << "Histofile: " << rootFileName << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  ///////////////////////////
  /// Global variables ////
  //////////////////////////
  TRootEvent* event = 0;
  // TRootRun *runInfos = new TRootRun();
  
  /////////////////////////////
  /// Object ID              ///
  /////////////////////////////
  // electron
  float el_pt_cut =20.; // 42
  float el_eta_cut = 2.5;
  float el_iso_cone  = 0.3;
  // reliso cut fabs(eta supercluster) <= 1.479 --> 0.107587 // (fabs(eta supercluster) > 1.479 && fabs(eta supercluster) < 2.5) --> 0.113254
  // muon
  float mu_pt_cut = 20.; // 40
  float mu_eta_cut = 2.4;
  float mu_iso_cut = 0.15;
  //jets
  float jet_pt_cut = 30.;
  float jet_eta_cut = 2.4;
  
  // convert into string
  
  std::ostringstream el_pt_cut_strs, el_eta_cut_strs, mu_pt_cut_strs, mu_eta_cut_strs, mu_iso_cut_strs, jet_pt_cut_strs, jet_eta_cut_strs;
  std::string el_pt_cut_str, el_eta_cut_str, mu_pt_cut_str, mu_eta_cut_str, mu_iso_cut_str, jet_pt_cut_str, jet_eta_cut_str;
  el_pt_cut_strs << el_pt_cut;
  el_eta_cut_strs << el_eta_cut;
  mu_pt_cut_strs << mu_pt_cut;
  mu_eta_cut_strs << mu_eta_cut;
  mu_iso_cut_strs << mu_iso_cut;
  jet_pt_cut_strs << jet_pt_cut;
  jet_eta_cut_strs << jet_eta_cut;
  el_pt_cut_str = el_pt_cut_strs.str();
  el_eta_cut_str = el_eta_cut_strs.str();
  mu_pt_cut_str = mu_pt_cut_strs.str();
  mu_eta_cut_str = mu_eta_cut_strs.str();
  mu_iso_cut_str = mu_iso_cut_strs.str();
  jet_pt_cut_str = jet_pt_cut_strs.str();
  jet_eta_cut_str = jet_eta_cut_strs.str();
  
  
  
  
  ////////////////////////////////////////////////////////////////////
  //////////////////  1D plots  //////////////////////////////
  ////////////////////////////////////////////////////////////////////
  histo1D["NbOfVertices"]                                  = new TH1F("NbOfVertices", "Nb. of vertices", 60, 0, 60);
  histo1D["cutFlow"]                                  = new TH1F( "cutFlow", "cutFlow", 15, -0.5, 14.5);
  histo1D["weightIndex"]				= new TH1F("weightIndex", "weightIndex", 5, -2.5,2.5); // 0: None; 1: scale_variation 1; 2: Central scale variation 1
  histo1D["nloweight"]				= new TH1F("nloweight", "nloweight", 200, -2.0, 2.0);
  histo1D["init_nPVs_before"]	                       = new TH1F("init_nPVs_before", "init_nPVs_before", 41,-0.5,40.5);
  histo1D["init_nPVs_after"]                        = new TH1F("init_nPVs_after", "init_nPVs_after", 41,-0.5,40.5);
  
  histo1D["nbMuons"]					= new TH1F("nbMuons","nbMuons",10,-0.5,9.5);
  histo1D["nbElectrons"]                                  = new TH1F("nbElectrons","nbElectrons",10,-0.5,9.5);
  histo1D["nbJets"]                                  = new TH1F("nbJets","nbJets",10,-0.5,9.5);
  
  
  
  /////////////////////////////////
  /// Matching
  ///////////////////////////////
  vector<TLorentzVector> mcParticlesTLV_charm, mcParticlesTLV_bottom, selectedJetsTLV, selectedMuonsTLV, selectedElectronsTLV, mcParticlesTLV_electrons, mcParticlesTLV_muons, mcParticlesTLV_Welectrons, mcParticlesTLV_Wmuons;
  // TLorentzVector cQuark, anticQuark;
  
  
  
  
  /////////////////////////////////
  //       Loop on datasets      //
  /////////////////////////////////
  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
  bool nlo = false;
  for (unsigned int d = 0; d < datasets.size(); d++)
  {
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
    
    
    double nloSF = 1;
    double sumWeights = 0;
    nlo = false;
    /////////
    string daName = datasets[d]->Name();
    float normfactor = datasets[d]->NormFactor();
    cout <<"found sample " << daName.c_str() << " with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    if(daName.find("Data")!=string::npos || daName.find("data")!=string::npos || daName.find("DATA")!=string::npos){
      isData = true;
    }
    if(daName.find("amc")!=string::npos) nlo = true;
    /////////////////////////////////////////
    ///    Calibrations                  ///
    ////////////////////////////////////////
    string CaliPath = "../TopTreeAnalysisBase/Calibrations/";
    string BCaliPath = CaliPath + "BTagging/CSVv2_13TeV_25ns_combToMujets.csv";
    if(!isData && !btagShape)
    {
      // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
      //	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_com@
      btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_76X_combToMujets.csv");
      btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets","central");
      if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
      {
        btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_mujets_central.root",false,30,999,2.4);
      }
      else
      {
        btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+"_mujets_central.root",false,30,999,2.4);
        //btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
      }
      
      
    }
    else if(!isData)
    {
      BTagCalibration calib_csvv2("csvv2", "../TopTreeAnalysisBase/Calibrations/BTagging/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv");
      reader_csvv2 = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                               BTagEntry::OP_RESHAPING, // operating point
                                               "iterativefit", // measurement type
                                               "central"); // systematics type  --> depending on JES up/Down andother reader is needed
      
      
    }
    
    LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring16MiniAODv2-Asympt.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting//pileup_2016Data80X_Run273158-276811Cert.root", "pileup", "pileup");
    
    //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)
    
    MuonSFWeight* muonSFWeightID_T = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);
    MuonSFWeight* muonSFWeightID_M = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true,  printLeptonSF, printLeptonSF);
    MuonSFWeight* muonSFWeightID_L = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, printLeptonSF, printLeptonSF);
    MuonSFWeight* muonSFWeightIso_TT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);  // Tight RelIso, Tight ID
    MuonSFWeight* muonSFWeightIso_TM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Tight RelIso, Medium ID
    MuonSFWeight* muonSFWeightIso_LT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Tight ID
    MuonSFWeight* muonSFWeightIso_LM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Medium ID
    
    
    
    
    string electronFile= "CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
    string electronRecoFile = "eleRECO.txt.egamma_SF2D.root";
    string elecHistName = "EGamma_SF2D";
    ElectronSFWeight* electronSFWeight = new ElectronSFWeight (CaliPath+"LeptonSF/"+electronFile,elecHistName, true,printLeptonSF, printLeptonSF); // (... , ... , debug, print warning)  i
    ElectronSFWeight* electronSFWeightReco = new ElectronSFWeight(CaliPath+"LeptonSF/"+electronRecoFile,elecHistName, true,printLeptonSF, printLeptonSF);
    
    vCorrParam.clear();
    if (isData)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    
    ////////////////////////////////////////////////////////////
    // Setup Date string and nTuple for output
    ///////////////////////////////////////////////////////////
    
    string channel_dir = "NtupleMakerOutput/Ntuples";
    string date_dir = channel_dir+"/Ntuples_" + dateString +"/";
    mkdir(channel_dir.c_str(),0777);
    mkdir(date_dir.c_str(),0777);
    
    
    string Ntupname = date_dir +"FCNC_3L_" + dName + "_"+  strJobNum + ".root";
    
    TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
    tupfile->cd();
    TTree* myTree = new TTree("tree","tree");
    TTree* baselineTree = new TTree("baselinetree","baselinetree");
    TTree* globalTree = new TTree("globaltree","globaltree");
    ///////////////////////////
    /// output tree
    ///////////////////////////
    // event related variables
    Int_t run_num;
    float i_channel; 
    Long64_t evt_num;
    Int_t lumi_num;
    Int_t nvtx;
    Int_t npu;
    Int_t PassedMETFilter;
    Int_t PassedGoodPV;
    Double_t cutstep[10];
    Int_t nCuts;
    Double_t puSF;
    Double_t btagSF;
    Double_t MuonIDSF[10];
    Double_t MuonIsoSF[10];
    Double_t MuonTrigSFv2[10];
    Double_t MuonTrigSFv3[10];
    Double_t ElectronSF[10];
    Int_t nofPosWeights;
    Int_t nofNegWeights;
    Int_t sumW;
    Int_t nEv;
    Double_t nloWeight; // for amc@nlo samples
    Int_t JERon;
    Int_t JESon;
    Double_t WPb_L;
    Double_t WPb_M;
    Double_t WPb_T;
    Int_t PassedMET;
    Int_t channelInt;
    
    Double_t pt_electron_1;
    Double_t pt_electron_2;
    Double_t pt_electron_3;
    Double_t pt_muon_1;
    Double_t pt_muon_2;
    Double_t pt_muon_3;
    Double_t pt_jet_1;
    Double_t pt_jet_2;
    Double_t pt_jet_3;
    
    
    Int_t nLeptons;
    // variables for electrons
    Int_t nElectrons;
    Double_t pt_electron[10];
    Double_t phi_electron[10];
    Double_t eta_electron[10];
    Double_t eta_superCluster_electron[10];
    Double_t E_electron[10];
    Double_t d0_electron[10];
    Double_t d0BeamSpot_electron[10];
    Double_t chargedHadronIso_electron[10];
    Double_t neutralHadronIso_electron[10];
    Double_t photonIso_electron[10];
    Double_t pfIso_electron[10];
    Int_t charge_electron[10];
    
    Double_t sigmaIEtaIEta_electron[10];
    Double_t deltaEtaIn_electron[10];
    Double_t deltaPhiIn_electron[10];
    Double_t hadronicOverEm_electron[10];
    Int_t missingHits_electron[10];
    Bool_t passConversion_electron[10];
    Bool_t isId_electron[10];
    Bool_t isIso_electron[10];
    
    Bool_t isEBEEGap[10];
    
    //variable for muons
    Int_t nMuons;
    Double_t pt_muon[10];
    Double_t phi_muon[10];
    Double_t eta_muon[10];
    Double_t E_muon[10];
    Double_t d0_muon[10];
    Double_t d0BeamSpot_muon[10];
    Double_t chargedHadronIso_muon[10];
    Double_t neutralHadronIso_muon[10];
    Double_t photonIso_muon[10];
    Double_t relIso_muon[10];
    Bool_t isId_muon[10];
    Bool_t isIso_muon[10];
    Double_t pfIso_muon[10];
    Int_t charge_muon[10];
    
    //variable for jets
    Int_t nJets;
    Int_t nJets_CSVL;
    Int_t nJets_CSVM;
    Int_t nJets_CSVT;
    Double_t pt_jet[20];
    Double_t px_jet[20];
    Double_t py_jet[20];
    Double_t pz_jet[20];
    Double_t phi_jet[20];
    Double_t eta_jet[20];
    Double_t E_jet[20];
    Int_t charge_jet[20];
    Double_t bdisc_jet[20];
    Double_t cdiscCvsL_jet[20];
    Double_t cdiscCvsL_jet_1;
    Double_t cdiscCvsB_jet_1;
    Double_t cdiscCvsB_jet[20];
    
    
    // variables for Zboson
    Double_t Zboson_M;
    Double_t Zboson_Px;
    Double_t Zboson_Py;
    Double_t Zboson_Pz;
    Double_t Zboson_Energy;
    
    // met
    Double_t met_Pt;
    Double_t met_Ptbf;
    Double_t met_Px;
    Double_t met_Py;
    Double_t met_Pz;
    Double_t met_Phi;
    Double_t met_Eta;
    
    Double_t mWt;
    Double_t FCNCtop_M;
    Double_t SMtop_M;
    Double_t cjet_Pt;
    Double_t mlb;
    Double_t dRWlepc;
    Double_t dRZb;
    Double_t dRZc;
    Double_t dRWlepb;
    Double_t dRSMFCNCtop;
    Double_t dPhiSMFCNCtop;
    Double_t dPhiWlepb;
    Double_t dPhiWlepc;
    Double_t dPhiZb;
    Double_t dPhiZc;
    
    Double_t FCNCtop_M_tagger;
    Double_t cjet_Pt_tagger;
    Double_t dRWlepc_tagger;
   
    Double_t dRZc_tagger;
  
    Double_t dRSMFCNCtop_tagger;
    Double_t dPhiSMFCNCtop_tagger;
    
    Double_t dPhiWlepc_tagger;
    
    Double_t dPhiZc_tagger;
    
    // mcparicles
    Int_t nMCParticles;
    Int_t mc_status[200];
    Int_t mc_pdgId[200];
    Int_t mc_mother[200];
    Int_t mc_granny[200];
    Double_t mc_pt[200];
    Double_t mc_phi[200];
    Double_t mc_eta[200];
    Double_t mc_E[200];
    Double_t mc_M[200];

    
    if(dName.find("NP_overlay_FCNC_TT")!=string::npos || dName.find("tZq")!=string::npos )
    {
      globalTree->Branch("nMatched_charm",&nMatched_charm,"nMatched_charm/I");
      globalTree->Branch("nNonMatched_charm",&nNonMatched_charm,"nNonMatched_charm/I");
      globalTree->Branch("nMatched_Zelec",&nMatched_Zelec,"nMatched_Zelec/I");
      globalTree->Branch("nNonMatched_Zelec",&nNonMatched_Zelec,"nNonMatched_Zelec/I");
      globalTree->Branch("nMatched_Zmu",&nMatched_Zmu,"nMatched_Zmu/I");
      globalTree->Branch("nNonMatched_Zmu",&nNonMatched_Zmu,"nNonMatched_Zmu/I");
      globalTree->Branch("nMatched_Welec",&nMatched_Welec,"nMatched_Welec/I");
      globalTree->Branch("nNonMatched_Welec",&nNonMatched_Welec,"nNonMatched_Welec/I");
      globalTree->Branch("nMatched_Wmu",&nMatched_Wmu,"nMatched_Wmu/I");
      globalTree->Branch("nNonMatched_Wmu",&nNonMatched_Wmu,"nNonMatched_Wmu/I");
      globalTree->Branch("nMatched_charm_tag",&nMatched_charm_tag,"nMatched_charm_tag/I");
      globalTree->Branch("nNonMatched_charm_tag",&nNonMatched_charm_tag,"nNonMatched_charm_tag/I");
      globalTree->Branch("nMatched_bottom",&nMatched_bottom,"nMatched_bottom/I");
      globalTree->Branch("nNonMatched_bottom",&nNonMatched_bottom,"nNonMatched_bottom/I");
      globalTree->Branch("nTagEqMass",&nTagEqMass,"nTagEqMass/I");
      globalTree->Branch("nTagNotEqMass", &nTagNotEqMass, "nTagNotEqMass/I");
      
      myTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
      myTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
      myTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
      myTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
      myTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
      myTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
      myTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
      myTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
      myTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
      myTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
      
      baselineTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
      baselineTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
      baselineTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
      baselineTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
      baselineTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
      baselineTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
      baselineTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
      baselineTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
      baselineTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
      baselineTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
    }
    
    
    
    
    // global data set variables
    Int_t nofEventsHLTv2;
    Int_t nofEventsHLTv3;
    int nTrigg;
    int n3lep;
    int nVetoMu;
    int nVetoEl;
    int nOS;
    int nZmass;
    int nJet;
    int nBJet;
    int nMWT;
    int nSMtop;
    int nMET;
    globalTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
    globalTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
    globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
    globalTree->Branch("nEv" , &nEv, "nEv/I");
    globalTree->Branch("sumW", &sumW, "sumW/I");
    globalTree->Branch("nCuts",&nCuts, "nCuts/I");
    globalTree->Branch("cutstep",&cutstep,"cutstep[nCuts]/D");
    globalTree->Branch("JERon",&JERon,"JERon/I");
    globalTree->Branch("JESon", &JESon, "JESon/I");
    globalTree->Branch("WPb_L", &WPb_L, "WPb_L/D");
    globalTree->Branch("WPb_M", &WPb_M, "WPb_M/D");
    globalTree->Branch("WPb_T", &WPb_T, "WPb_T/D");
    
    globalTree->Branch("nTrigg", &nTrigg, "nTrigg/I");
    globalTree->Branch("n3lep", &n3lep, "n3lep/I");
    globalTree->Branch("nVetoMu", &nVetoMu, "nVetoMu/I");
    globalTree->Branch("nVetoEl", &nVetoEl, "nVetoEl/I");
    globalTree->Branch("nOS", &nOS, "nOS/I");
    globalTree->Branch("nZmass",&nZmass, "nZmass/I");
    globalTree->Branch("nJet", &nJet, "nJet/I");
    globalTree->Branch("nBJet",&nBJet, "nBJet/I");
    globalTree->Branch("nMWT", &nMWT, "nMWT/I");
    globalTree->Branch("nSMtop",&nSMtop, "nSMtop/I");
    globalTree->Branch("nMET",&nMET, "nMET/I");
    
    
    
    
    
    // event related variables
    myTree->Branch("channelInt", &channelInt, "channelInt/I");
    baselineTree->Branch("channelInt", &channelInt, "channelInt/I");
    myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/L");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("puSF",&puSF,"puSF/D");
    myTree->Branch("btagSF",&btagSF,"btagSF/D");
    myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
    myTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I");
    myTree->Branch("PassedGoodPV", &PassedGoodPV,"PassedGoodPV/I");
    
    baselineTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I");
    baselineTree->Branch("PassedGoodPV", &PassedGoodPV,"PassedGoodPV/I");
    baselineTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    baselineTree->Branch("run_num",&run_num,"run_num/I");
    baselineTree->Branch("evt_num",&evt_num,"evt_num/I");
    baselineTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    baselineTree->Branch("nvtx",&nvtx,"nvtx/I");
    baselineTree->Branch("npu",&npu,"npu/I");
    baselineTree->Branch("puSF",&puSF,"puSF/D");
    baselineTree->Branch("btagSF",&btagSF,"btagSF/D");
    baselineTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
    // electrons
    myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//
    myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
    myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
    myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
    myTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
    myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
    myTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
    myTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
    myTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
    myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
    myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
    myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
    myTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
    myTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
    myTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
    myTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
    myTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
    myTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
    myTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
    myTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
    myTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
    myTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
    myTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/D");
    myTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/D");
    myTree->Branch("pt_electron_3",&pt_electron_3,"pt_electron_3/D");
    
    
    baselineTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//
    baselineTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    baselineTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
    baselineTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
    baselineTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
    baselineTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
    baselineTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
    baselineTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
    baselineTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
    baselineTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
    baselineTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
    baselineTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
    baselineTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
    baselineTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
    baselineTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
    baselineTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
    baselineTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
    baselineTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
    baselineTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
    baselineTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
    baselineTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
    baselineTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
    baselineTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
    baselineTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/D");
    baselineTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/D");
    baselineTree->Branch("pt_electron_3",&pt_electron_3,"pt_electron_3/D");
    
    // muons
    myTree->Branch("nMuons",&nMuons, "nMuons/I");
    myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
    myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
    myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
    myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
    myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
    myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
    myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
    myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
    myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
    myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
    myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
    myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
    myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
    myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
    myTree->Branch("pt_muon_1",&pt_muon_1,"pt_muon_1/D");
    myTree->Branch("pt_muon_2",&pt_muon_2,"pt_muon_2/D");
    myTree->Branch("pt_muon_3",&pt_muon_3,"pt_muon_3/D");
    
    baselineTree->Branch("nMuons",&nMuons, "nMuons/I");
    baselineTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    baselineTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    baselineTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    baselineTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
    baselineTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
    baselineTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
    baselineTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
    baselineTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
    baselineTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
    baselineTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
    baselineTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
    baselineTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
    baselineTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
    baselineTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
    baselineTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
    baselineTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
    baselineTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
    baselineTree->Branch("pt_muon_1",&pt_muon_1,"pt_muon_1/D");
    baselineTree->Branch("pt_muon_2",&pt_muon_2,"pt_muon_2/D");
    baselineTree->Branch("pt_muon_3",&pt_muon_3,"pt_muon_3/D");
    
    // jets
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
    myTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
    myTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
    myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
    myTree->Branch("px_jet",px_jet,"px_jet[nJets]/D");
    myTree->Branch("py_jet",py_jet,"py_jet[nJets]/D");
    myTree->Branch("pz_jet",pz_jet,"pz_jet[nJets]/D");
    myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
    myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
    myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
    myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
    myTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
    myTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
    myTree->Branch("cdiscCvsL_jet_1",&cdiscCvsL_jet_1,"cdiscCvsL_jet_1/D");
    myTree->Branch("cdiscCvsB_jet_1",&cdiscCvsB_jet_1,"cdiscCvsB_jet_1/D");
    myTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/D");
    myTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/D");
    myTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/D");
    
    baselineTree->Branch("nJets",&nJets,"nJets/I");
    baselineTree->Branch("cdiscCvsL_jet_1",&cdiscCvsL_jet_1,"cdiscCvsL_jet_1/D");
    baselineTree->Branch("cdiscCvsB_jet_1",&cdiscCvsB_jet_1,"cdiscCvsB_jet_1/D");
    baselineTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
    baselineTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
    baselineTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
    baselineTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
    baselineTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
    baselineTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
    baselineTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
    baselineTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    baselineTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
    baselineTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
    baselineTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
    baselineTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/D");
    baselineTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/D");
    baselineTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/D");
    
    // Zboson
    myTree->Branch("Zboson_M",&Zboson_M,"Zboson_M/D");
    baselineTree->Branch("Zboson_M",&Zboson_M,"Zboson_M/D");
    myTree->Branch("mWt",&mWt,"mWt/D");
    baselineTree->Branch("mWt",&mWt,"mWt/D");
    myTree->Branch("FCNCtop_M",&FCNCtop_M,"FCNCtop_M/D");
    myTree->Branch("FCNCtop_tagger",&FCNCtop_M_tagger,"FCNCtop_M_tagger/D");
    baselineTree->Branch("FCNCtop_tagger",&FCNCtop_M_tagger,"FCNCtop_M_tagger/D");
    myTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/D");
    baselineTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/D");
    myTree->Branch("Zboson_Px",&Zboson_Px,"Zboson_Px/D");
    myTree->Branch("Zboson_Py",&Zboson_Py,"Zboson_Py/D");
    myTree->Branch("Zboson_Pz",&Zboson_Pz,"Zboson_Pz/D");
    myTree->Branch("Zboson_Energy",&Zboson_Energy,"Zboson_Energy/D");
    myTree->Branch("cjet_Pt",&cjet_Pt,"cjet_Pt/D");
    baselineTree->Branch("cjet_Pt",&cjet_Pt,"cjet_Pt/D");
    myTree->Branch("cjet_Pt_tagger",&cjet_Pt_tagger,"cjet_Pt_tagger/D");
    baselineTree->Branch("cjet_Pt_tagger",&cjet_Pt_tagger,"cjet_Pt_tagger/D");
    myTree->Branch("mlb",&mlb,"mlb/D");
    baselineTree->Branch("mlb",&mlb,"mlb/D");
    myTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/D");
     myTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/D");
    baselineTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/D");
    myTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/D");
    myTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/D");
    myTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/D");
    baselineTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/D");
    myTree->Branch("dRZb",&dRZb,"dRZb/D");
    myTree->Branch("dRZc",&dRZc,"dRZc/D");
    myTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/D");
    baselineTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/D");
    myTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/D");
    myTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/D");
    baselineTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/D");
    myTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/D");
    myTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/D");
    baselineTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/D");
    myTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/D");
    myTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/D");
    myTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/D");
    myTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/D");
    baselineTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/D");
    baselineTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/D");
    baselineTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/D");
    baselineTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/D");
    baselineTree->Branch("dRZb",&dRZb,"dRZb/D");
    baselineTree->Branch("dRZc",&dRZc,"dRZc/D");
    baselineTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/D");
    baselineTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/D");
    baselineTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/D");
    baselineTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/D");
    baselineTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/D");
    
    // met
    myTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
    myTree->Branch("met_Ptbf", &met_Ptbf, "met_Ptbf/D");
    myTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
    myTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
    myTree->Branch("met_Px", &met_Px, "met_Px/D");
    myTree->Branch("met_Py", &met_Py, "met_Py/D");
    myTree->Branch("met_Pz", &met_Pz, "met_Pz/D");
    
    baselineTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
    baselineTree->Branch("met_Ptbf", &met_Ptbf, "met_Ptbf/D");
    baselineTree->Branch("met_Px", &met_Px, "met_Px/D");
    baselineTree->Branch("met_Py", &met_Py, "met_Py/D");
    baselineTree->Branch("met_Pz", &met_Pz, "met_Pz/D");
    baselineTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
    baselineTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
    
    
    
    /////////////////////////
    //// Corrections/trigger ///
    ///////////////////////////
    
    /// book triggers
    trigger_mumu->bookTriggers(isData);
    trigger_ee->bookTriggers(isData);
    trigger_emu->bookTriggers(isData);
    trigger_mumumu->bookTriggers(isData);
    trigger_eee->bookTriggers(isData);
    trigger_emumu_mumue->bookTriggers(isData);
    trigger_mu->bookTriggers(isData);
    trigger_e->bookTriggers(isData);
    
    
    
    //////////////////////////////////////////////////
    // Pre-event loop definitions
    /////////////////////////////////////////////////
    
    int itrigger = -1, previousRun = -1, start = 0;
    int currentRun;
    int iFile = -1;
    unsigned int ending = datasets[d]->NofEvtsToRunOver();
    cout <<"Number of events = "<<  ending  <<endl;
    
    string previousFilename = "";
    int event_start = startEvent;
    
    double currentfrac =0.;
    double end_d;
    if(endEvent > ending)
      end_d = ending;
    else
      end_d = endEvent;
    
    int nEvents = end_d - event_start;
    cout <<"Will run over "<<  (end_d - event_start) << " events..."<<endl;
    cout <<"Starting event = = = = "<< event_start  << endl;
    if(end_d < startEvent)
    {
      cout << "Starting event larger than number of events.  Exiting." << endl;
      exit(3);
    }
    
    if (verbose == 0) cout << " - Loop over events " << endl;
    
    //define object containers
    
    // initial variables
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_jets_corrected;
    vector < TRootGenJet* >   genjets;
    vector < TRootMET* >      mets;
    vector<TRootElectron*> selectedElectrons;
    vector<TRootElectron*> selectedLooseElectrons;
    vector<TRootPFJet*>    selectedJets;
    vector<TRootPFJet*>    PreselectedJets;
    vector<TRootMuon*>     selectedMuons;
    vector<TRootMuon*>     selectedLooseMuons;
    vector<TRootPFJet*>      selectedCSVLBJets;
    vector<TRootPFJet*>      selectedCSVMBJets;
    vector<TRootPFJet*>      selectedCSVTBJets;
    vector<TRootPFJet*>      selectedCSVLLJets;
    vector<TRootPFJet*>      selectedCSVMLJets;
    vector<TRootPFJet*>      selectedCSVTLJets;
    vector<TRootMCParticle*> mcParticles;
    vector <TRootPFJet*>     selectednonCSVLJets;
    
    TLorentzVector Zboson;
    TLorentzVector Zlep0;
    TLorentzVector Zlep1;
    TLorentzVector Wlep;
    TLorentzVector SMbjet;
    TLorentzVector cjet;
    TLorentzVector cjet_tagger;
    TLorentzVector SMtop;
    TLorentzVector FCNCtop;
    TLorentzVector FCNCtop_tagger;
    vector<TLorentzVector> AssignedLeptons;
    //////////////////////////////////////
    // Begin Event Loop
    //////////////////////////////////////
    nbEvents = 0;
    nofEventsHLTv2 = 0;
    nofEventsHLTv3 = 0;
    nofPosWeights = 0;
    nofNegWeights = 0;
    float eventweight = 1;
    bool continueFlow ;
    nbSelectedEvents = 0;
    int nbEvents_0 = 0;
    int nbEvents_test = 0;
    int nbEvents_1 = 0;
    int nbEvents_1m = 0;
    int nbEvents_2m = 0;
    int nbEvents_2 = 0;
    int nbEvents_3 = 0;
    int nbEvents_4 = 0;
    int nbEvents_5 = 0;
    int nbEvents_6 = 0;
    int nbEvents_7 = 0;
    int nbEvents_8 = 0;
    int nbEvents_9 = 0;
    bool debug = false;
    vector <int> selections;
    std::ostringstream  selectionsnb;
    bool   passedMET = false;
    bool   HBHEnoise = false;
    bool   HBHEIso = false;
    bool   CSCTight = false;
    bool   EcalDead = false;
    bool    eeBad = false;
    bool   lep3 = false;
    TLorentzVector metTLV;
    TLorentzVector metTLVbf;
    string TriggBits;
    float  pt_lept1;
    float  pt_lept2;
    float  pt_lept3;
    float  iso_lept1;
    float  iso_lept2;
    float  iso_lept3;
    bool id_lept1 = 1;
    bool id_lept2  = 1;
    bool id_lept3 = 1;
    float leading_jet_btagDiscr;
    float leading_jetPt;
    float met;
    nMatched_bottom = 0;
    nNonMatched_bottom = 0;
    nMatched_charm = 0;
    nNonMatched_charm = 0;
    nMatched_charm_tag = 0;
    nNonMatched_charm_tag = 0;
    nNonMatched_Zmu = 0;
    nNonMatched_Zelec = 0;
    nMatched_Zmu = 0;
    nMatched_Zelec = 0;
    nNonMatched_Wmu = 0;
    nNonMatched_Welec = 0;
    nMatched_Wmu = 0;
    nMatched_Welec = 0;
    nTagEqMass = 0;
    nTagNotEqMass = 0;
    for (unsigned int ievt = event_start; ievt < end_d; ievt++)
    {
      elecbool = false;
      mubool = false;
      muIndices.clear();
      elecIndices.clear();
      WmuIndices.clear();
      WelecIndices.clear();
      eventSelected = false;
      baseSelected = false;
      continueFlow = true;
      lep3 = false;
      AssignedLeptons.clear();
      leading_jetPt = 0.;
      met = 0.;
      leading_jet_btagDiscr = 0.;
      TriggBits = "";
      pt_lept1 = pt_lept2 = pt_lept3 = 0. ;
      metTLV.Clear();
      metTLVbf.Clear();
      metTLV.SetPxPyPzE(0,0,0,0);
      selections.clear();
      bool lepsel = false;
      selectionsnb.clear();
      selectionsnb.str(std::string());
      mcParticles.clear();
      /// mcparticles
      nMCParticles = -1;
      for (Int_t i = 0; i < 200; i++)
      {
        mc_status[i] = -1;
        mc_pdgId[i] = 0;
        mc_mother[i] = 0;
        mc_granny[i] = 0;
        mc_pt[i] = 0.;
        mc_phi[i] = 0.;
        mc_eta[i] = 0.;
        mc_E[i] = 0.;
        mc_M[i] = 0.;
      }
      nCuts = 0;
      passedMET = false;
      HBHEnoise = false;
      HBHEIso = false;
      CSCTight = false;
      EcalDead = false;
      eeBad = false;
      eventweight = 1;
      if(verbose == 0 ) cout << "new event " << ievt << endl;
      double ievt_d = ievt;
      debug = false;
      if (verbose == 0 ) debug = true;
      currentfrac = ievt_d/end_d;
      if (debug)cout << endl << endl << "Starting a new event loop!"<<endl;
      
      if(ievt%10000 == 0)
      {
        std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC
        << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
      }
      
      float scaleFactor = 1.;  // scale factor for the event
      
      
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
      if(debug) cout << "event loaded" << endl;
      genjets.clear();
      if(!isData){
        genjets = treeLoader.LoadGenJet(ievt,false);  //needed for JER
      }
      init_jets_corrected = init_jets;
      
      if(verbose==0)
      {
        cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
        cout <<"Number of Muons Loaded: " << init_muons.size() <<endl;
        cout << "Number of Jets  Loaded: " << init_jets.size() << endl;
      }
      
      //  take the event
      datasets[d]->eventTree()->LoadTree(ievt);
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      int currentRun = event->runId();
      run_num = event->runId();
      evt_num = event->eventId();
      HBHEnoise = event->getHBHENoiseFilter();
      HBHEIso = event->getHBHENoiseIsoFilter();
      CSCTight = event->getCSCTightHalo2015Filter();
      EcalDead = event->getEcalDeadCellTriggerPrimitiveFilter();
      eeBad = event->getEEBadScFilter();
      
      lumi_num=event->lumiBlockId();
      nvtx = vertex.size();
      npu = (int) event->nTruePU();
      
      /////////////////////////////////////
      //  fix negative weights for amc@nlo///
      /////////////////////////////////////
      if(debug) cout << "amc fixing" << endl;
      double hasNegWeight = false;
      double mc_baseweight = 1;
      if((!isData && dName.find("NP")==string::npos) && (event->getWeight(1001) != -9999.))
      {
        mc_baseweight =  event->getWeight(1001)/abs(event->originalXWGTUP());
        //mc_scaleupweight = event->getWeight(1005)/abs(event->originalXWGTUP());
        //mc_scaledownweight = event->getWeight(1009)/abs(event->originalXWGTUP());
        if(mc_baseweight >= 0)
        {
          nofPosWeights++;
          histo1D["weightIndex"]->Fill(1.,1.);
          
        }
        else
        {
          if(nlo) hasNegWeight = true;
          nofNegWeights++;
          histo1D["weightIndex"]->Fill(-1.,1.);
        }
      }
      if( (!isData && dName.find("NP")==string::npos) && (event->getWeight(1) != -9999. ))
      {
        mc_baseweight =  event->getWeight(1)/abs(event->originalXWGTUP());
        //mc_scaleupweight = event->getWeight(5)/abs(event->originalXWGTUP());
        //mc_scaledownweight = event->getWeight(9)/abs(event->originalXWGTUP());
        if(mc_baseweight >= 0)
        {
          nofPosWeights++;
          histo1D["weightIndex"]->Fill(2.,1.);
          
        }
        else
        {
          if(nlo) hasNegWeight = true;
          nofNegWeights++;
          histo1D["weightIndex"]->Fill(-2.,1.);
        }
        
        
      }
      if((!isData && dName.find("FCNC")==string::npos))
      {
        if ( event->getWeight(1001) == -9999. && event->getWeight(1) == -9999. )
        {
          cout << "WARNING: No weight found for event " << ievt << " in dataset " << dName << endl;
          cout << "         Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
          cout << "         Weight type is different from 'scale_variation' (1001) or 'Central scale variation' (1)." << endl;
        }
        if ( event->getWeight(1001) != -9999. && event->getWeight(1) != -9999. )
        {
          cout << "WARNING: Two weight types found for event " << ievt << " in dataset " << dName << endl;
          cout << "         Event Id: " << event->eventId() << "  Run Id: " << event->runId() << "  Lumi block Id: " << event->lumiBlockId() << endl;
          cout << "         Check which weight type should be used when." << endl;
        }
        
        nloWeight = mc_baseweight;
        histo1D["nloweight"]->Fill(mc_baseweight, 1.);
        sumWeights += mc_baseweight;
        
        
      }
      ///////////////////////////////////////////
      //  Trigger
      ///////////////////////////////////////////
      
      bool trigged = false;
      bool trigged_mumu = false;
      bool trigged_ee = false;
      bool trigged_emu = false;
      bool trigged_mumumu = false;
      bool trigged_eee = false;
      bool trigged_emumu_mumue = false;
      bool trigged_mu = false;
      bool trigged_e = false;
      bool filechanged = false;
      bool runchanged = false;
      
      if(runHLT)
      {
        trigger_mumu->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_mumu =  trigger_mumu->checkIfFired();
        trigger_ee->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_ee =  trigger_ee->checkIfFired();
        trigger_emu->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_emu =  trigger_emu->checkIfFired();
        trigger_mumumu->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_mumumu =  trigger_mumumu->checkIfFired();
        trigger_eee->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_eee =  trigger_eee->checkIfFired();
        trigger_emumu_mumue->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_emumu_mumue =  trigger_emumu_mumue->checkIfFired();
        trigger_mu->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_mu =  trigger_mu->checkIfFired();
        trigger_e->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_e =  trigger_e->checkIfFired();
        
        
        bool emdataset = dName.find("MuonEG")!=string::npos;
        bool mmdataset = dName.find("DoubleM")!=string::npos;
        bool eedataset = dName.find("DoubleE")!=string::npos;
        bool mdataset = dName.find("SingleM")!=string::npos;
        bool edataset = dName.find("SingleE")!=string::npos;
        
        bool EM = false;
        bool MM = false;
        bool EE = false;
        bool E = false;
        bool M = false;
        int result_trigger = 0;
        
        
        if(isData){
          EM = (trigged_emumu_mumue|| trigged_emu);
          MM = (trigged_mumu || trigged_mumumu ) ;
          EE = (trigged_ee || trigged_eee );
          M  = ( trigged_mu );
          E  = (trigged_e);
        }
        else{
          EM = (trigged_emumu_mumue|| trigged_emu);
          MM = (trigged_mumu || trigged_mumumu ) ;
          EE = (trigged_ee || trigged_eee );
          M  = ( trigged_mu );
          E  = (trigged_e);
        }
        if ( EM &&                               (emdataset) ) result_trigger = 1;
        if ( MM && !EM &&                        (mmdataset) ) result_trigger = 1;
        if ( EE && !EM && !MM &&                 (eedataset) ) result_trigger = 1;
        if ( M  && !EM && !MM && !EE &&          (mdataset ) ) result_trigger = 1;
        if ( E  && !EM && !MM && !EE && !M &&    (edataset ) ) result_trigger = 1;
        
        if ( EM &&                               !isData ) result_trigger = 1;
        if ( MM && !EM &&                        !isData ) result_trigger = 1;
        if ( EE && !EM && !MM &&                 !isData ) result_trigger = 1;
        if ( M  && !EM && !MM && !EE &&          !isData ) result_trigger = 1;
        if ( E  && !EM && !MM && !EE && !M &&    !isData ) result_trigger = 1;
        
        
        
        
        trigged = result_trigger;
        if(dName.find("NP")!=string::npos) trigged = true;
        
        
      }
      else if(!runHLT && previousFilename != currentFilename)
      {
        filechanged = true;
        previousFilename = currentFilename;
        iFile++;
        cout << "File changed!!! => iFile = " << iFile << endl;
        trigged = true;
        
      }
      else if(!runHLT)
      {
        trigged = true;
      }
      
      if(verbose == 0) cout << "Apply trigger? " << runHLT << " trigged? " << trigged << endl;
      
      ////////////////////////////
      ///// JES - JER smearing     ////
      //////////////////////////
      JERon = 0;
      if(applyJER && !isData)
      {
        jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        JERon = 1;
      }
      JESon = 0;
      if(applyJES && !isData)
      {
        jetTools->correctJets(init_jets_corrected,event->fixedGridRhoFastjetAll() ,false);
        JESon = 1;
      }
      
      ///////////////////////////////////////////////////////////
      // Event selection
      ///////////////////////////////////////////////////////////
      
      // Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets,event->fixedGridRhoFastjetAll());
      PreselectedJets.clear();
      PreselectedJets  = selection.GetSelectedJets(jet_pt_cut,jet_eta_cut, true, "Loose");
      selectedMuons.clear();
      selectedLooseMuons.clear();
      selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Tight", "Spring15");
      selectedLooseMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut,0.2, "Loose", "Spring15");
      // pt, eta, iso // run normally
      selectedElectrons.clear();
      selectedLooseElectrons.clear();
      selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Tight","Spring15_25ns",true);// pt, eta
      selectedLooseElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Veto","Spring15_25ns",true);// pt, eta
      /// For MC Information
      mcParticles.clear();
      if(!isData) treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
      if(!isData) sort(mcParticles.begin(),mcParticles.end(),HighestPt());
      if (dName.find("NP_overlay_FCNC_TT")!=string::npos || dName.find("tZq")!=string::npos )
      {
        nMCParticles = mcParticles.size();
        if (nMCParticles > maxMCParticles) maxMCParticles = nMCParticles;
        for (Int_t iMC = 0; iMC < nMCParticles; iMC++)
        {
          mc_status[iMC] = mcParticles[iMC]->status();
          mc_pdgId[iMC] = mcParticles[iMC]->type();
          mc_mother[iMC] = mcParticles[iMC]->motherType();
          mc_granny[iMC] = mcParticles[iMC]->grannyType();
          mc_pt[iMC] = mcParticles[iMC]->Pt();
          mc_phi[iMC] = mcParticles[iMC]->Phi();
          mc_eta[iMC] = mcParticles[iMC]->Eta();
          mc_E[iMC] = mcParticles[iMC]->E();
          mc_M[iMC] = mcParticles[iMC]->M();
        }
      }
      // void TTreeLoader::LoadMCEvent(int, TopTree::TRootNPGenEvent*, std::vector<TopTree::TRootMCParticle*>&, bool)
      if (verbose==0) cout <<"Number of Muons, Electrons, Jets  ===>  " << endl << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< PreselectedJets.size()   << endl;
      selectedJets.clear();
      if(applyJetLeptonCleaning){
        bool PushBack = true;
        for(int iJ = 0; iJ < PreselectedJets.size() ; iJ++)
        {
          PushBack = true;
          for(int iM = 0; iM < selectedMuons.size(); iM++){
            if( PreselectedJets[iJ]->DeltaR(*selectedMuons[iM]) < 0.4) {
              PushBack = false;
              break;
            }
          }
          if(!PushBack) continue;
          for(int iE = 0; iE < selectedElectrons.size(); iE++){
            if( PreselectedJets[iJ]->DeltaR(*selectedElectrons[iE]) < 0.3) {
              PushBack = false;
              break;
            }
          }
          if(PushBack) selectedJets.push_back(PreselectedJets[iJ]);
        }
      }
      else if(!applyJetLeptonCleaning)   selectedJets = PreselectedJets;
      if(debug) cout << evt_num << " init " << init_jets_corrected.size() << " sel "  << selectedJets.size() << " bf cleaning " << PreselectedJets.size() << endl;
      
      ////////////////////////////////////////////////
      // Pre cut operations
      ////////////////////////////////////////////////
      // Apply primary vertex selection
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      // Met filters
      if(HBHEnoise && HBHEIso && CSCTight && EcalDead && eeBad && isGoodPV) passedMET = true;
      PassedMETFilter = passedMET;
      
      //////////////////////////////////////
      //   B jet selection	       ////
      ///////////////////////////////////////
      
      selectedCSVLBJets.clear();
      selectedCSVMBJets.clear();
      selectedCSVTBJets.clear();
      selectedCSVLLJets.clear();
      selectedCSVMLJets.clear();
      selectedCSVTLJets.clear();
      selectednonCSVLJets.clear();
      for(unsigned int iJ = 0; iJ < selectedJets.size(); iJ++)
      {
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose) selectedCSVLBJets.push_back(selectedJets[iJ]);
        else selectedCSVLLJets.push_back(selectedJets[iJ]);
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium) selectedCSVMBJets.push_back(selectedJets[iJ]);
        else selectedCSVMLJets.push_back(selectedJets[iJ]);
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight) selectedCSVTBJets.push_back(selectedJets[iJ]);
        else selectedCSVTLJets.push_back(selectedJets[iJ]);
        
      }
      WPb_L =  workingpointvalue_Loose;
      WPb_M =  workingpointvalue_Medium;
      WPb_T =  workingpointvalue_Tight;
      
      ////////////////////////////////////
      //   Event Weights               ///
      ///////////////////////////////////
      float btagWeight  =  1.;
      float bTagEff = 1.;
      if( fillBtagHisto && !isData && !btagShape)
      {
        btwt->FillMCEfficiencyHistos(selectedJets);
        
      }
      else if( !fillBtagHisto && !isData && !btagShape)
      {
        btagWeight =  btwt->getMCEventWeight(selectedJets,false);
        
      }
      else if( !isData && btagShape)
      {
        for(int intJet = 0; intJet < selectedJets.size(); intJet++)
        {
          float jetpt = selectedJets[intJet]->Pt();
          if(jetpt > 1000.) jetpt = 999.;
          float jeteta = selectedJets[intJet]->Eta();
          float jetdisc = selectedJets[intJet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
          BTagEntry::JetFlavor jflav;
          int jetpartonflav = std::abs(selectedJets[intJet]->partonFlavour());
          if(debug) cout<<"parton flavour: "<<jetpartonflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
          if(jetpartonflav == 5){
            jflav = BTagEntry::FLAV_B;
          }
          else if(jetpartonflav == 4){
            jflav = BTagEntry::FLAV_C;
          }
          else{
            jflav = BTagEntry::FLAV_UDSG;
          }
          bTagEff = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          btagWeight *= bTagEff;
          
        }
        
      }
      float PUweight = 1;
      if(!isData)
      {
        PUweight = LumiWeights.ITweight((int)event->nTruePU());
        
        
      }
      
      
      ///////////////////////////////
      //// Matching
      //////////////////////////////
      if(dName.find("NP_overlay_FCNC_TT")!=string::npos || dName.find("tZq")!=string::npos) matching = true;
      //cout << "matching " << matching << endl;
      
      if(matching){
        //cout << "in matching" << endl;
        int pdgID_charm = 4;
        int pdgID_bottom = 5;
        int pdgId_Z = 23;
        int pdgId_top = 6;
        int pdgId_W = 24;
        int pdgID_electron = 11;
        int pdgID_muon = 13;
        vector<TRootMCParticle*> mcParticlesMatching_;
        
        
        if(dName.find("NP_overlay_FCNC_TT")!=string::npos || dName.find("tZq")!=string::npos){
          mcParticlesTLV_charm.clear(); selectedJetsTLV.clear();
          mcParticlesTLV_bottom.clear();
          selectedElectronsTLV.clear();
          selectedMuonsTLV.clear();
          mcParticlesTLV_electrons.clear();
          mcParticlesTLV_muons.clear();
          mcParticlesTLV_Welectrons.clear();
          mcParticlesTLV_Wmuons.clear();
          
          for (unsigned int i = 0; i < mcParticles.size(); i++)
          {
            if(verbose>3)  cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
            
            if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
            

            
            if ( fabs(mcParticles[i]->type()) == pdgID_charm && fabs(mcParticles[i]->motherType()) == pdgId_top){
              mcParticlesTLV_charm.push_back(*mcParticles[i]);
              //mcParticlesMatching_.push_back(mcParticles[i]);
              
            }
            if ( fabs(mcParticles[i]->type()) == pdgID_bottom && fabs(mcParticles[i]->motherType()) == pdgId_top){
              mcParticlesTLV_bottom.push_back(*mcParticles[i]);
              //mcParticlesMatching_.push_back(mcParticles[i]);
              
            }
            
            if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_Z){
              mcParticlesTLV_electrons.push_back(*mcParticles[i]);
            }
            if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_Z){
              mcParticlesTLV_muons.push_back(*mcParticles[i]);
            }
            if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_W){
              mcParticlesTLV_Welectrons.push_back(*mcParticles[i]);
            }
            if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_W){
              mcParticlesTLV_Wmuons.push_back(*mcParticles[i]);
            }
            
            
          }
        }
        
        
        
        // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
        for (unsigned int i = 0; i < selectedJets.size(); i++)
        {
          selectedJetsTLV.push_back(*selectedJets[i]);
        }
        for (unsigned int i = 0; i < selectedMuons.size(); i++)
        {
          selectedMuonsTLV.push_back(*selectedMuons[i]);
        }
        for (unsigned int i = 0; i < selectedElectrons.size(); i++)
        {
          selectedElectronsTLV.push_back(*selectedElectrons[i]);
        }
        
        //cout << "selectedJetsTLV.size() = " << selectedJetsTLV.size() << endl;
        
        JetPartonMatching matching_charm = JetPartonMatching(mcParticlesTLV_charm, selectedJetsTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
        JetPartonMatching matching_bottom = JetPartonMatching(mcParticlesTLV_bottom, selectedJetsTLV, 2, true, true, 0.3);
        JetPartonMatching matching_elec = JetPartonMatching(mcParticlesTLV_electrons, selectedElectronsTLV, 2, true, true, 0.3);
        JetPartonMatching matching_muon = JetPartonMatching(mcParticlesTLV_muons, selectedMuonsTLV,2, true, true, 0.3);
        JetPartonMatching matching_Welec = JetPartonMatching(mcParticlesTLV_Welectrons, selectedElectronsTLV, 2, true, true, 0.3);
        JetPartonMatching matching_Wmuon = JetPartonMatching(mcParticlesTLV_Wmuons, selectedMuonsTLV,2, true, true, 0.3);
        
        if (matching_charm.getNumberOfAvailableCombinations() != 1)
          cerr << "matching_charm.getNumberOfAvailableCombinations() = " << matching_charm.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
        if (matching_bottom.getNumberOfAvailableCombinations() != 1)
          cerr << "matching_bottom.getNumberOfAvailableCombinations() = " << matching_bottom.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
        
        
        JetPartonPair_bottom.clear(); // First one is jet number, second one is mcParticle number
        JetPartonPair_charm.clear();
        JetPartonPair_electron.clear();
        JetPartonPair_muon.clear();
        JetPartonPair_Welectron.clear();
        JetPartonPair_Wmuon.clear();
        //cout << "mcParticlesTLV.size() " << mcParticlesTLV.size() << endl;
        
        for (unsigned int i = 0; i < mcParticlesTLV_charm.size(); i++)
        {
          
          int matchedJetNumber_charm = matching_charm.getMatchForParton(i, 0);
          if (matchedJetNumber_charm > -1){
            JetPartonPair_charm.push_back( pair<unsigned int, unsigned int> (matchedJetNumber_charm, i) );
            //cout << "Matched Jet number " << matchedJetNumber << endl;
          }
        }
        for (unsigned int i = 0; i < mcParticlesTLV_bottom.size(); i++)
        {
          
          int matchedJetNumber_bottom = matching_bottom.getMatchForParton(i, 0);
          if (matchedJetNumber_bottom > -1){
            JetPartonPair_bottom.push_back( pair<unsigned int, unsigned int> (matchedJetNumber_bottom, i) );
            //cout << "Matched Jet number " << matchedJetNumber << endl;
          }
        }
        for (unsigned int i = 0; i < mcParticlesTLV_electrons.size(); i++)
        {
          
          int matchedNumber = matching_elec.getMatchForParton(i, 0);
          if (matchedNumber > -1){
            JetPartonPair_electron.push_back( pair<unsigned int, unsigned int> (matchedNumber, i) );
            
          }
        }
        
        for (unsigned int i = 0; i < mcParticlesTLV_muons.size(); i++)
        {
          
          int matchedNumber = matching_muon.getMatchForParton(i, 0);
          if (matchedNumber > -1){
            JetPartonPair_muon.push_back( pair<unsigned int, unsigned int> (matchedNumber, i) );
            
          }
        }
        
        for (unsigned int i = 0; i < mcParticlesTLV_Welectrons.size(); i++)
        {
          
          int matchedNumber = matching_Welec.getMatchForParton(i, 0);
          if (matchedNumber > -1){
            JetPartonPair_Welectron.push_back( pair<unsigned int, unsigned int> (matchedNumber, i) );
            
          }
        }
        
        for (unsigned int i = 0; i < mcParticlesTLV_Wmuons.size(); i++)
        {
          
          int matchedNumber = matching_Wmuon.getMatchForParton(i, 0);
          if (matchedNumber > -1){
            JetPartonPair_Wmuon.push_back( pair<unsigned int, unsigned int> (matchedNumber, i) );
            
          }
        }
        
        
      } // end matching
      
      ////////////////////////////////////
      //   Determine eventweight        ///
      /////////////////////////////////
      
      histo1D["init_nPVs_before"]->Fill(vertex.size(), eventweight);
      if(applyPU && !isData)  eventweight *= PUweight;
      histo1D["init_nPVs_after"]->Fill(vertex.size(), eventweight);
      
      //////////////////////////////////////////////////////
      // Applying baseline selection
      //////////////////////////////////////////////////////
      continueFlow = true;
      nbEvents++;
      eventweight = 1.;
      if(trigged){
        selections.push_back(1);
        if(continueFlow){
          histo1D["cutFlow"]->Fill(0., eventweight);
          nCuts++;
          nbEvents_0++;
        }
      }
      else{
        selections.push_back(0);
        continueFlow = false;
      }
      
      // to be ok with triggers
      if(dName.find("DoubleEG")!=string::npos && selectedElectrons.size() < 2) { continueFlow = false; }
      else if(dName.find("DoubleEG")!=string::npos) { nbEvents_test++ ;}
      if(dName.find("DoubleMu")!=string::npos && selectedMuons.size() < 2) { continueFlow = false; }
      else if(dName.find("DoubleMu")!=string::npos) { nbEvents_test++ ;}
      if(dName.find("MuonEG")!=string::npos && (selectedElectrons.size() < 1 || selectedMuons.size() < 1)) { continueFlow = false; }
      else if(dName.find("MuonEG")!=string::npos){ nbEvents_test++ ;}
      
      if(((selectedMuons.size() + selectedElectrons.size()) != 3)){
        selections.push_back(0);
        continueFlow = false;
      }
      else if( ((selectedMuons.size() + selectedElectrons.size()) == 3)){
        selections.push_back(1);
	       if(continueFlow){
           histo1D["cutFlow"]->Fill(1., eventweight);
           nCuts++;
           nbEvents_1++;
         }
        lep3 = true;
        if(selectedMuons.size() == 3) {channelInt = 0; i_channel = 0;} 
        else if(selectedElectrons.size() == 3) {channelInt = 3; i_channel = 3;}
        else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {channelInt = 2; i_channel = 2; }
        else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){channelInt = 1; i_channel = 1; }
        else cout << "ERROR no channel selected" << endl;
      }
      
      //cout << "LOOKING AT CHANNEL " << channelInt << endl;
      
      if(selectedMuons.size() == selectedLooseMuons.size() && continueFlow) nbEvents_1m++;
      else continueFlow = false;
      if(selectedLooseElectrons.size() == selectedElectrons.size() && continueFlow)  nbEvents_2m++;
      else continueFlow = false;
      
      if((selectedMuons.size() != selectedLooseMuons.size()) || (selectedLooseElectrons.size() != selectedElectrons.size())){
        selections.push_back(0);
        continueFlow = false;
      }
      else {
        selections.push_back(1);
        if(continueFlow) {
          lepsel = true;
          histo1D["cutFlow"]->Fill(2., eventweight);
        }
      }
      
      double met_px = mets[0]->Px();
      double met_py = mets[0]->Py();
      met_Pt = sqrt(met_px*met_px + met_py*met_py);
      met = met_Pt;
      met_Phi = mets[0]->Phi();
      met_Eta = mets[0]->Eta();
      
      puSF = PUweight;
      if(!isData) btagSF = btagWeight;
      if(isData) btagSF = 1.;
      
      Zlep0.Clear();
      Zlep1.Clear();
      Wlep.Clear();
      Wlep.SetPxPyPzE(0,0,0,0);
      
      // check sign
      
      if(lep3){
        //cout << "assigning leptons " << endl;
        AssignedLeptons = LeptonAssigner(selectedElectrons, selectedMuons);
        
        if(Assigned){
          
          Zlep0.SetPxPyPzE(AssignedLeptons[0].Px(), AssignedLeptons[0].Py(), AssignedLeptons[0].Pz(), AssignedLeptons[0].Energy());
          Zlep1.SetPxPyPzE(AssignedLeptons[1].Px(), AssignedLeptons[1].Py(), AssignedLeptons[1].Pz(), AssignedLeptons[1].Energy());
          Wlep.SetPxPyPzE(AssignedLeptons[2].Px(), AssignedLeptons[2].Py(), AssignedLeptons[2].Pz(), AssignedLeptons[2].Energy());
          
          //double phis = Wlep.Phi() - mets[0]->Phi();
          //double cosphis = TMath::Cos(phis);
          mWt = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_px)*(Wlep.Px() + met_px) - (Wlep.Py() + met_py)* (Wlep.Py() + met_py));
          //mWtsecond = TMath::Sqrt(2*Wlep.Pt() * met_Pt*(1-cosphis));
          
          
          nCuts++;
          nbEvents_2++;
          
          Zboson.Clear();
          
          Zboson.SetPxPyPzE(( Zlep0 + Zlep1).Px() ,( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Pz(),( Zlep0 + Zlep1).Energy()) ;
          Zboson_M = (Zlep0+Zlep1).M();
          Zboson_Px = ( Zlep0 + Zlep1).Px();
          Zboson_Py = ( Zlep0 + Zlep1).Py();
          Zboson_Pz = ( Zlep0 + Zlep1).Pz();
          Zboson_Energy = ( Zlep0 + Zlep1).Energy();
          
          
        }
        else{
          continueFlow = false;
          Zboson_M = -5;
          Zboson_Px = -5;
          Zboson_Py = -5;
          Zboson_Pz = -5;
          Zboson_Energy = -5;
          mWt = -5;
        }
      }
      else{
        continueFlow = false;
        Zboson_M = -5;
        Zboson_Px = -5;
        Zboson_Py = -5;
        Zboson_Pz = -5;
        Zboson_Energy = -5;
        mWt = -5;
      }
      
      
      
      if(Zboson_M < 76 || Zboson_M > 106)
      {
        selections.push_back(0);
        continueFlow = false;
      }
      else{
        selections.push_back(1);
        if(continueFlow){
          nCuts++;
          nbEvents_3++;
          histo1D["cutFlow"]->Fill(3., eventweight);
          baseSelected = true;
        }
      }
      
      
      if(selectedJets.size() < 2){
        selections.push_back(0);
        continueFlow = false;
      }
      else{
        selections.push_back(1);
        if(continueFlow){
          histo1D["cutFlow"]->Fill(4., eventweight);
          nCuts++;
          nbEvents_4++;
        }
      }
      if(selectedCSVLBJets.size()  < 1){
        selections.push_back(0);
	       continueFlow = false;
      }
      else{
        selections.push_back(1);
        if(continueFlow){
          histo1D["cutFlow"]->Fill(5., eventweight);
          nCuts++;
          nbEvents_5++;
        }
      }
      

      
      if(false){ // no mWT cut
        selections.push_back(0);
	       continueFlow = false;
      }
      else{
        selections.push_back(1);
        if(continueFlow){
          //			 mWtFile << evt_num << endl;
          histo1D["cutFlow"]->Fill(6., eventweight);
          nCuts++;
          nbEvents_6++;
        }
      }
      
      //            double met_pz =  MEtz(Wmu, Wel, Wlep, met_px, met_py);
    //  cout << "MET reconstruc" << endl;
      double met_pz = 0.; // has to be adapted !!!
      metTLVbf.SetPxPyPzE(met_px,met_py,met_pz,TMath::Sqrt(met_px*met_px+met_py*met_py+met_pz*met_pz));
      met_Ptbf = metTLVbf.Pt();
      metTLV = MetzCalculator(Wlep, metTLVbf);
      
    //  cout << "Met reconstructed" << endl;
      met_Px = metTLV.Px();
      met_Py = metTLV.Py();
      met_Pz = metTLV.Pz();
      SMbjet.Clear();
      SMtop.Clear();
      int SMbjetindex = -5;
      if(selectedCSVLBJets.size() > 0 ){
        // cout << "bjets " << selectedCSVLBJets.size() << " jets " << selectedJets.size() << endl;
        SMbjetindex = SMjetCalculator(selectedJets, verbose);
        // cout << "SMbjetindex " << SMbjetindex << endl;
        SMbjet.SetPxPyPzE(selectedJets[SMbjetindex]->Px(),selectedJets[SMbjetindex]->Py(),selectedJets[SMbjetindex]->Pz(),selectedJets[SMbjetindex]->Energy());
        
        if(Assigned && continueFlow)  {
          SMtop_M = (Wlep+SMbjet+metTLV).M();
          SMtop.SetPxPyPzE((SMbjet.Px()+Wlep.Px()+metTLV.Px()),(SMbjet.Py()+Wlep.Py()+metTLV.Py()),(SMbjet.Pz()+Wlep.Pz()+metTLV.Pz()),(SMbjet.Energy()+Wlep.Energy()+metTLV.Energy()));
          mlb = (Wlep+SMbjet).M();
          dRWlepb = Wlep.DeltaR(SMbjet);
          dRZb = Zboson.DeltaR(SMbjet);
          dPhiWlepb = Wlep.DeltaPhi(SMbjet);
          dPhiZb = Zboson.DeltaPhi(SMbjet);
        }
        else {
          SMtop_M = -5.;
          mlb = -5.;
          dRWlepb = -5;
          dRZb = -5;
        }
      }
      else  SMtop_M = -5. ;
      
      cjet.Clear();
      cjet_tagger.Clear();
      FCNCtop.Clear();
      FCNCtop_tagger.Clear();
      int cjetindex = -5;
      int cjetindex_tagger = -5;
      if(Assigned && continueFlow && selectedJets.size()>1) {
        cjetindex = FCNCjetCalculator(selectedJets,Zboson ,SMbjetindex, 3);
        //cout << "bjet_index " << SMbjetindex << endl;
        cjetindex_tagger = FCNCjetCalculatorTagger(selectedJets,SMbjetindex, 3);
        // cout << "cjet index " << cjetindex << endl;
        cjet.SetPxPyPzE(selectedJets[cjetindex]->Px(),selectedJets[cjetindex]->Py(),selectedJets[cjetindex]->Pz(),selectedJets[cjetindex]->Energy());
        //cout << "cjetindex_tagger " << cjetindex_tagger << endl;
        cjet_tagger.SetPxPyPzE(selectedJets[cjetindex_tagger]->Px(),selectedJets[cjetindex_tagger]->Py(),selectedJets[cjetindex_tagger]->Pz(),selectedJets[cjetindex_tagger]->Energy());
        
        FCNCtop.SetPxPyPzE((cjet+Zboson).Px(), (cjet+Zboson).Py(), (cjet+Zboson).Pz(), (cjet+Zboson).Energy());
        FCNCtop_M = (Zlep0+Zlep1+cjet).M();
        cjet_Pt = TMath::Sqrt(cjet.Px()*cjet.Px()+cjet.Py()*cjet.Py());
        dRZc = Zboson.DeltaR(cjet);
        dRWlepc = Wlep.DeltaR(cjet);
        dPhiZc = Zboson.DeltaPhi(cjet);
        dPhiWlepc = Wlep.DeltaPhi(cjet);
        dRSMFCNCtop = SMtop.DeltaR(FCNCtop);
        dPhiSMFCNCtop = SMtop.DeltaPhi(FCNCtop);
        
        
        FCNCtop_tagger.SetPxPyPzE((cjet_tagger+Zboson).Px(), (cjet_tagger+Zboson).Py(), (cjet_tagger+Zboson).Pz(), (cjet_tagger+Zboson).Energy());
        FCNCtop_M_tagger = (Zlep0+Zlep1+cjet_tagger).M();
        cjet_Pt_tagger = TMath::Sqrt(cjet_tagger.Px()*cjet_tagger.Px()+cjet_tagger.Py()*cjet_tagger.Py());
        dRZc_tagger = Zboson.DeltaR(cjet_tagger);
        dRWlepc_tagger = Wlep.DeltaR(cjet_tagger);
        dPhiZc_tagger = Zboson.DeltaPhi(cjet_tagger);
        dPhiWlepc_tagger = Wlep.DeltaPhi(cjet_tagger);
        dRSMFCNCtop_tagger = SMtop.DeltaR(FCNCtop_tagger);
        dPhiSMFCNCtop_tagger = SMtop.DeltaPhi(FCNCtop_tagger);
        
      }
      else {
        FCNCtop_M = -5.;
        // cout << "event: " << evt_num << " - Zboson.M()= " << Zboson.M() << " - cjet.M()= " << cjet.M() << " - top.M()= " << (Zboson+cjet).M() << endl;
        dRZc = -5;
        dRWlepc = -5;
        dPhiWlepc = -5;
        dPhiZc = -5;
        dRSMFCNCtop = -5 ;
        dPhiSMFCNCtop = -5;
        cjet_Pt = -5;
        FCNCtop_M_tagger = -5.;
        // cout << "event: " << evt_num << " - Zboson.M()= " << Zboson.M() << " - cjet.M()= " << cjet.M() << " - top.M()= " << (Zboson+cjet).M() << endl;
        dRZc_tagger = -5;
        dRWlepc_tagger = -5;
        dPhiWlepc_tagger = -5;
        dPhiZc_tagger = -5;
        dRSMFCNCtop_tagger = -5 ;
        dPhiSMFCNCtop_tagger = -5;
        cjet_Pt_tagger = -5;
      }
      
      
      
      if(SMtop_M < 95 || SMtop_M > 200 ){
        selections.push_back(0);
        continueFlow = false;
        //		 continue;
      }
      else{
        selections.push_back(1);
        if(continueFlow){
          histo1D["cutFlow"]->Fill(7., eventweight);
          nCuts++;
          nbEvents_7++;
          
        }
      }
      
      //	   if(continueFlow)  eventSelected = true;
      //	   else eventSelected = false;
      if(passedMET && continueFlow){
        histo1D["cutFlow"]->Fill(8., eventweight);
        nCuts++;
        nbEvents_8++;
        eventSelected = true;
      }
      //////////////////////////////////////
      //  DO STUFF WITH SELECTED EVENTS ////
      //////////////////////////////////////
      // fill the tree
      if(eventSelected || baseSelected){
        
        nJets = 0;
        for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
        {
          
          pt_jet[nJets]=selectedJets[seljet]->Pt();
          px_jet[nJets]=selectedJets[seljet]->Px();
          py_jet[nJets]=selectedJets[seljet]->Py();
          pz_jet[nJets]=selectedJets[seljet]->Pz();
          phi_jet[nJets]=selectedJets[seljet]->Phi();
          eta_jet[nJets]=selectedJets[seljet]->Eta();
          E_jet[nJets]=selectedJets[seljet]->E();
          charge_jet[nJets]=selectedJets[seljet]->charge();
          bdisc_jet[nJets]=selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() ;
          cdiscCvsB_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsBJetTags() ;
          cdiscCvsL_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsLJetTags() ;
          nJets++;
        }
        if(selectedJets.size()>0) cdiscCvsB_jet_1 = selectedJets[0]->ctag_pfCombinedCvsBJetTags();
        if(selectedJets.size()>0) cdiscCvsL_jet_1 = selectedJets[0]->ctag_pfCombinedCvsLJetTags();
        if(selectedJets.size()>0) pt_jet_1 = selectedJets[0]->Pt();
        if(selectedJets.size()>1) pt_jet_2 = selectedJets[1]->Pt();
        if(selectedJets.size()>2) pt_jet_3 = selectedJets[2]->Pt();
        nJets_CSVT =  selectedCSVTBJets.size();
        nJets_CSVM =  selectedCSVMBJets.size();
        nJets_CSVL =  selectedCSVLBJets.size();
        nMuons = 0;
        for (Int_t selmu =0; selmu < selectedMuons.size() ; selmu++ )
        {
          
        	 pt_muon[nMuons]=selectedMuons[selmu]->Pt();
        	 phi_muon[nMuons]=selectedMuons[selmu]->Phi();
        	 eta_muon[nMuons]=selectedMuons[selmu]->Eta();
        	 E_muon[nMuons]=selectedMuons[selmu]->E();
          
          pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
          if(!isData)
        	 {
             MuonIDSF[nMuons] = muonSFWeightID_T->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
             MuonIsoSF[nMuons] =  muonSFWeightIso_TT->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
           }
        	 else
           {
             MuonIDSF[nMuons] = 1.;
             MuonIsoSF[nMuons] = 1.;
           }
        	 charge_muon[nMuons]=selectedMuons[selmu]->charge();
        	 nMuons++;
        }
        if(selectedMuons.size()>0) pt_muon_1 = selectedMuons[0]->Pt();
        if(selectedMuons.size()>1) pt_muon_2 = selectedMuons[1]->Pt();
        if(selectedMuons.size()>2) pt_muon_3 = selectedMuons[2]->Pt();
        nElectrons=0;
        for (Int_t selel =0; selel < selectedElectrons.size() ; selel++ )
        {
          
        	 pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
        	 phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
        	 eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
        	 eta_superCluster_electron[nElectrons]=selectedElectrons[selel]->superClusterEta();
        	 E_electron[nElectrons]=selectedElectrons[selel]->E();
        	 pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
        	 charge_electron[nElectrons]=selectedElectrons[selel]->charge();
        	 if(!isData) ElectronSF[nElectrons] = electronSFWeight->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
           else ElectronSF[nElectrons] = 1.;
          
        	 nElectrons++;
        }
        if(selectedElectrons.size()>0) pt_electron_1 = selectedElectrons[0]->Pt();
        if(selectedElectrons.size()>1) pt_electron_2 = selectedElectrons[1]->Pt();
        if(selectedElectrons.size()>2) pt_electron_3 = selectedElectrons[2]->Pt();
        
        
        nLeptons = nMuons + nElectrons;
        
      }
      
      if(eventSelected){
        nbSelectedEvents++;
        if(matching && JetPartonPair_bottom.size()>0){
          
          if(JetPartonPair_bottom[0].first == SMbjetindex)  nMatched_bottom++;
          else nNonMatched_bottom++;
          
        }
        if(matching && JetPartonPair_charm.size()>0){
          
          if(JetPartonPair_charm[0].first == cjetindex)  nMatched_charm++;
          else nNonMatched_charm++;
        }
        if(matching && JetPartonPair_charm.size()>0){
          
          if(JetPartonPair_charm[0].first == cjetindex_tagger)  nMatched_charm_tag++;
          else nNonMatched_charm_tag++;
        }
        if(matching)
        {
          if(cjetindex_tagger == cjetindex) nTagEqMass++;
          else nTagNotEqMass++;
          
        }
        if(matching && elecIndices.size() == 2){
          if( (elecIndices[0] == JetPartonPair_electron[0].first) && (elecIndices[1] == JetPartonPair_electron[1].first)) nMatched_Zelec++;
          else if((elecIndices[1] == JetPartonPair_electron[0].first) && (elecIndices[0] == JetPartonPair_electron[1].first)) nMatched_Zelec++;
          else nNonMatched_Zelec++;
        
        }
        if(matching && muIndices.size() == 2){
          if( (muIndices[0] == JetPartonPair_muon[0].first) && (muIndices[1] == JetPartonPair_muon[1].first)) nMatched_Zmu++;
          else if((muIndices[1] == JetPartonPair_muon[0].first) && (muIndices[0] == JetPartonPair_muon[1].first)) nMatched_Zmu++;
          else nNonMatched_Zmu++;
          
        }
        if(matching && WelecIndices.size() == 1){
          if( (WelecIndices[0] == JetPartonPair_Welectron[0].first)) nMatched_Welec++;
          else nNonMatched_Welec++;
          
        }
        if(matching && WmuIndices.size() == 1){
          if( (WmuIndices[0] == JetPartonPair_Wmuon[0].first) ) nMatched_Wmu++;
          else nNonMatched_Wmu++;
          
        }
        
        //cout << "SIZE mu " << JetPartonPair_muon.size() << " elec " << JetPartonPair_electron.size() << endl;  ;
        myTree->Fill();
 	    }
      if(baseSelected){ baselineTree->Fill(); }
      if(selections.size() != 8) cout << "ERROR SOMETHING WENT WRONG WITH THE SELECTIONS " << endl;
      for(int inb = 0; inb <selections.size(); inb++)
      {
        selectionsnb << selections[inb];
      }
      
      
      
      
      
      
    } // end eventloop
    cutstep[0] = nbEvents_0;
    cutstep[1] = nbEvents_1;
    cutstep[2] = nbEvents_2;
    cutstep[3] = nbEvents_3;
    cutstep[4] = nbEvents_4;
    cutstep[5] = nbEvents_5;
    cutstep[6] = nbEvents_6;
    cutstep[7] = nbEvents_7;
    cutstep[8] = nbEvents_8;
    cutstep[9] = nbEvents_9;
    
    cout << "nbEvents_0 trigg: " << nbEvents_0 << endl;
    cout << "trigger req for data " << nbEvents_test << endl;
    cout << "nbEvents_1 3 lep: " << nbEvents_1 << endl;
    cout << "nbEvents_1m  veto mu: " << nbEvents_1m << endl;
    cout << "nbEvents_2m veto el: " << nbEvents_2m << endl;
    cout << "nbEvents_2 OS: " << nbEvents_2 << endl;
    cout << "nbEvents_3 Z mass: " << nbEvents_3 << endl;
    cout << "nbEvents_4 >0 jet: " << nbEvents_4 << endl;
    cout << "nbEvents_5 1 bjet: " << nbEvents_5 << endl;
    cout << "nbEvents_6 mWt: " << nbEvents_6 << endl;
    cout << "nbEvents_7 SMtop: " << nbEvents_7 << endl;
    cout << "nbEvents_8 MET: " << nbEvents_8 << endl;
    
    
    nTrigg = nbEvents_0;
    n3lep = nbEvents_1;
    nVetoMu = nbEvents_1m;
    nVetoEl = nbEvents_2m;
    nOS = nbEvents_2;
    nZmass = nbEvents_3;
    nJet = nbEvents_4;
    nBJet = nbEvents_5;
    nMWT = nbEvents_6;
    nSMtop = nbEvents_7;
    nMET= nbEvents_8;
    
    
    
    
    
    if(matching) {
      cout << "Percentage matched charm: " << (double) nMatched_charm / (nMatched_charm + nNonMatched_charm) << endl;
      cout << "Percentage matched charm tag: " << (double) nMatched_charm_tag / (nMatched_charm_tag + nNonMatched_charm_tag) << endl;
      cout << "Percentage tag equal mass charm: " << (double) nTagEqMass / (nTagEqMass+nTagNotEqMass) << endl;
      cout << "Percentage matched bottom: " << (double) nMatched_bottom / (nMatched_bottom + nNonMatched_bottom) << endl;
      cout << "Percentage matched Z elec: " << (double) nMatched_Zelec/ (nMatched_Zelec + nNonMatched_Zelec) << endl;
      cout << "Percentage matched Z mu: " << (double) nMatched_Zmu/ (nMatched_Zmu + nNonMatched_Zmu) << endl;
      cout << "Percentage matched W elec: " << (double) nMatched_Welec/ (nMatched_Welec + nNonMatched_Welec) << endl;
      cout << "Percentage matched W mu: " << (double) nMatched_Wmu/ (nMatched_Wmu + nNonMatched_Wmu) << endl;
      cout << "Percentage matched Z: " << (double) (nMatched_Zmu+nMatched_Zelec)/ (nMatched_Zmu +nMatched_Zelec+nNonMatched_Zelec+ nNonMatched_Zmu) << endl;
      cout << "Percentage matched W: " << (double) (nMatched_Welec+nMatched_Wmu)/ (nMatched_Welec + nMatched_Wmu + nNonMatched_Wmu+ nNonMatched_Welec) << endl;
    }
    //	for(int j = 0; j < 9; j++){       cout << cutstep[j] << endl; }
    sumW = (int) sumWeights;
    nEv = (int) nEvents;
    
    globalTree->Fill();
    if(verbose == 0) cout << "end eventloop" << endl;
    
    cout << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl;
    cout << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
    cout << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
    
    if (! isData  )
    {
      cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
      cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
      cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents) << endl;
      
      // Determine scale factor due to negative weights
      nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
      cout << "This corresponds to an event scale factor of " << nloSF  << endl;
    }
    tupfile->cd();
    myTree->Write();
    globalTree->Write();
    baselineTree->Write();
    tupfile->Close();
    delete tupfile;
    if(!isData && !btagShape) delete btwt;
    //    if(!isData && fillBtagHisto) delete btwt;
    treeLoader.UnLoadDataset();
  } //End Loop on Datasets
  
  
  
  /////////////
  // Writing //
  /////////////
  
  cout << " - Writing outputs to the files ..." << endl;
  
  
  
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
  
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  
  return 0;
};


/////////////////////////////////////// FUNCTIONS
//



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
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, true);
  string month_str = ConvertIntToString(month, true);
  string day_str = ConvertIntToString(day, true);
  string hour_str = ConvertIntToString(hour, true);
  string min_str = ConvertIntToString(min, true);
  //string sec_str = ConvertIntToString(sec, true);
  
  string date_str = year_str + month_str + day_str; //+ "_" + hour_str + min_str;
  return date_str;
};


double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy)
{
  
  double emu = Wlep.E();
  double pxmu = Wlep.Px();
  double pymu = Wlep.Py();
  double pzmu = Wlep.Pz();
  double pxnu = MetPx;
  double pynu = MetPy;
  double pznu = 0.;
  if(el && ! mu) M_mu = M_el;
  
  double a = M_W*M_W - M_mu*M_mu + 2.0*pxmu*pxnu + 2.0*pymu*pynu;
  double A = 4.0*(emu*emu - pzmu*pzmu);
  double B = -4.0*a*pzmu;
  double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
  
  
  bool isComplex_ = false;
  double tmproot = B*B - 4.0*A*C;
  
  if (tmproot<0) {
    isComplex_= true;
    pznu = - B/(2*A); // take real part of complex roots
  }
  else {
    isComplex_ = false;
    double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
    double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
    
    if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
    else pznu = tmpsol1;
    
    
  }
  return pznu;
  
}
;

int FCNCjetCalculator(std::vector<TRootPFJet*> Jets, TLorentzVector recoZ ,int index, int verb)
{
  
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -1;
  if(Jets.size() > 1){
    //cout << " non bjets: " << nonBJets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size(); iJ++)
    {
      if(iJ == index) continue;
      TLorentzVector Jet;
      Jet.SetPxPyPzE(Jets[iJ]->Px(),Jets[iJ]->Py(),Jets[iJ]->Pz(),Jets[iJ]->Energy());
      //cout << iJ << " tempMinM " << TempMinMass << " newmass " << (recoZ+Jet).M() ;
      if(fabs((recoZ+Jet).M() - TopMass) < TempMinMass)
      {
        TempMinMass = fabs((recoZ+Jet).M() - TopMass);
        Jetcandidate.SetPxPyPzE(Jet.Px(), Jet.Py(), Jet.Pz(), Jet.E());
        NbInColl = iJ;
        
      }
      //cout << " NbInColl is " << iJ << endl;
      //  269297.249181
    }
  }
  else{
    NbInColl = -5,
    cout << "no cjets available" << endl;
  }
  return NbInColl;
};

int FCNCjetCalculatorTagger(std::vector<TRootPFJet*> Jets,int index, int verb)
{
  
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -5;
  if(Jets.size() > 2){
   // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      for(int kJ = 1; kJ < Jets.size(); kJ++){
        if(kJ == index) continue;
        if((Jets[iJ]->ctag_pfCombinedCvsLJetTags()+Jets[iJ]->ctag_pfCombinedCvsBJetTags())>=(Jets[kJ]->ctag_pfCombinedCvsLJetTags()+Jets[kJ]->ctag_pfCombinedCvsBJetTags())) NbInColl = iJ;
        else NbInColl = kJ;
      }
    }
  }
  else if(Jets.size() == 2){
    if(index == 0) NbInColl = 1;
    if(index == 1) NbInColl = 0;
  }
  else{
    NbInColl = -5,
    cout << "no cjets available" << endl;
  }
  return NbInColl;
};


float EffectiveAreaRho(TRootElectron *el, float rho_)
{
  double EffectiveArea = 0.;
  // Updated to Spring 2015 EA from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_14/RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt#L8
  if (fabs(el->superClusterEta()) >= 0.0   && fabs(el->superClusterEta()) < 1.0   ) EffectiveArea = 0.1752;
  if (fabs(el->superClusterEta()) >= 1.0   && fabs(el->superClusterEta()) < 1.479 ) EffectiveArea = 0.1862;
  if (fabs(el->superClusterEta()) >= 1.479 && fabs(el->superClusterEta()) < 2.0   ) EffectiveArea = 0.1411;
  if (fabs(el->superClusterEta()) >= 2.0   && fabs(el->superClusterEta()) < 2.2   ) EffectiveArea = 0.1534;
  if (fabs(el->superClusterEta()) >= 2.2   && fabs(el->superClusterEta()) < 2.3   ) EffectiveArea = 0.1903;
  if (fabs(el->superClusterEta()) >= 2.3   && fabs(el->superClusterEta()) < 2.4   ) EffectiveArea = 0.2243;
  if (fabs(el->superClusterEta()) >= 2.4   && fabs(el->superClusterEta()) < 5.0   ) EffectiveArea = 0.2687;
  if (fabs(el->superClusterEta()) >= 5.0) EffectiveArea = -9999;
  
  double isocorr = 0;
  
  isocorr = rho_*EffectiveArea;
  
  return isocorr;
};

float EffectiveArea(TRootElectron *el)
{
  double EffectiveArea = 0.;
  
  if (fabs(el->superClusterEta()) >= 0.0   && fabs(el->superClusterEta()) < 1.0   ) EffectiveArea = 0.1752;
  if (fabs(el->superClusterEta()) >= 1.0   && fabs(el->superClusterEta()) < 1.479 ) EffectiveArea = 0.1862;
  if (fabs(el->superClusterEta()) >= 1.479 && fabs(el->superClusterEta()) < 2.0   ) EffectiveArea = 0.1411;
  if (fabs(el->superClusterEta()) >= 2.0   && fabs(el->superClusterEta()) < 2.2   ) EffectiveArea = 0.1534;
  if (fabs(el->superClusterEta()) >= 2.2   && fabs(el->superClusterEta()) < 2.3   ) EffectiveArea = 0.1903;
  if (fabs(el->superClusterEta()) >= 2.3   && fabs(el->superClusterEta()) < 2.4   ) EffectiveArea = 0.2243;
  if (fabs(el->superClusterEta()) >= 2.4   && fabs(el->superClusterEta()) < 5.0   ) EffectiveArea = 0.2687;
  if (fabs(el->superClusterEta()) >= 5.0) EffectiveArea = -9999;
  
  
  return EffectiveArea;
};


float relPfIsoEl(TRootElectron *el, float _rho)
{
  float isoCorr = (el->neutralHadronIso(3) + el->photonIso(3) - EffectiveAreaRho(el,_rho));
  //    float isolation = (el->chargedHadronIso(3) + (isoCorr > 0.0 ? isoCorr : 0.0))/(el->Pt());
  float isolation = (el->chargedHadronIso(3) + std::max(el->neutralHadronIso(3)+el->photonIso(3)-EffectiveAreaRho(el,_rho),float(0.)))/(el->Pt());
  return isolation;
  
};


float IsoDBeta(TRootMuon *mu)
{
  float iso = (mu->chargedHadronIso(4) + std::max(0.0, mu->neutralHadronIso(4) + mu->photonIso(4) - 0.5*mu->puChargedHadronIso(4)))/mu->Pt();
  
  return iso;
  
}

vector <TLorentzVector> LeptonAssigner(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons)
{
  //  cout << " in assigner " << endl;
  vector<TLorentzVector> ReturnColl;
  Assigned = false;
  
  if(electrons.size() + muons.size() != 3){
    cout << " WARNING: not 3 leptons " << endl;
    cout << "muons " << muons.size() << " electrons " << electrons.size() << endl;
    return ReturnColl;
  }
  elecbool = false;
  mubool = false;
  elecIndices.clear();
  muIndices.clear();
  WelecIndices.clear();
  WmuIndices.clear();
  //  cout << " in 3 lep " << endl;
  
  TLorentzVector Zlepcan0;
  Zlepcan0.SetPxPyPzE(0.,0.,0.,0.);
  TLorentzVector Zlepcan1;
  Zlepcan1.SetPxPyPzE(0.,0.,0.,0.);
  TLorentzVector Wlepcan;
  Wlepcan.SetPxPyPzE(0.,0.,0.,0.);
  
  if(electrons.size() == 2){
    //cout << "2 electr " << electrons[0]->charge() << " " << electrons[1]->charge() << endl;
    if(electrons[0]->charge() != electrons[1]->charge()){
      Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Wlepcan.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Assigned = true;
      elecbool = true;
      elecIndices.push_back(0);
      elecIndices.push_back(1);
      WmuIndices.push_back(0);
    }
  }
  else if(muons.size() == 2){
    //    cout << "2 muons" << endl;
    if(muons[0]->charge() != muons[1]->charge()){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Wlepcan.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Assigned = true;
      mubool = true;
      muIndices.push_back(0);
      muIndices.push_back(1);
      WelecIndices.push_back(0);
    }
  }
  else if(electrons.size() ==3){
    //    cout << " 3 electrons " << endl;
    bool can01 = false;
    bool can02= false;
    bool can12 = false;
    elecbool = true;
    if(electrons[0]->charge() != electrons[1]->charge()) can01 = true;
    if(electrons[0]->charge() != electrons[2]->charge()) can02 = true;
    if(electrons[2]->charge() != electrons[1]->charge()) can12 = true;
    
    double mass01 = 9999.;
    double mass02 = 9999.;
    double mass12 = 9999.;
    TLorentzVector temp0;
    temp0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
    TLorentzVector temp1;
    temp1.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
    TLorentzVector temp2;
    temp2.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
    if(can01) mass01 = fabs(91.1-(temp1+temp0).M());
    if(can02) mass02 = fabs(91.1-(temp2+temp0).M());
    if(can12) mass12 = fabs(91.1-(temp1+temp2).M());
    if(mass01 <= mass02 && mass01 <= mass12){
      Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Wlepcan.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Assigned = true;
      elecIndices.push_back(0); elecIndices.push_back(1); WelecIndices.push_back(2);
    }
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Wlepcan.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Assigned = true;
      elecIndices.push_back(0); elecIndices.push_back(2); WelecIndices.push_back(1);
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Wlepcan.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Assigned = true;
      elecIndices.push_back(1); elecIndices.push_back(2); WelecIndices.push_back(0);
    }
  }
  else if(muons.size() == 3){
    bool can01 = false;
    bool can02= false;
    bool can12 = false;
    mubool = true;
    if(muons[0]->charge() != muons[1]->charge()) can01 = true;
    if(muons[0]->charge() != muons[2]->charge()) can02 = true;
    if(muons[2]->charge() != muons[1]->charge()) can12 = true;
    
    double mass01 = 9999.;
    double mass02 = 9999.;
    double mass12 = 9999.;
    TLorentzVector temp0;
    temp0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
    TLorentzVector temp1;
    temp1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
    TLorentzVector temp2;
    temp2.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
    if(can01) mass01 = fabs(91.1-(temp1+temp0).M());
    if(can02) mass02 = fabs(91.1-(temp2+temp0).M());
    if(can12) mass12 = fabs(91.1-(temp1+temp2).M());
    if(mass01 <= mass02 && mass01 <= mass12){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Wlepcan.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Assigned = true;
      muIndices.push_back(0); muIndices.push_back(1); WmuIndices.push_back(2);
    }
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Wlepcan.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Assigned = true;
      muIndices.push_back(0); muIndices.push_back(2); WmuIndices.push_back(1);
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Wlepcan.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Assigned = true;
      muIndices.push_back(1); muIndices.push_back(2); WmuIndices.push_back(0);
    }
  }
  if(Assigned){
    ReturnColl.push_back(Zlepcan0);
    ReturnColl.push_back(Zlepcan1);
    ReturnColl.push_back(Wlepcan);
  }
  if(!Assigned){
    //    cout << " WARNING: leptons not set for assignment " << endl;
    return ReturnColl;
  }
  
  
  return ReturnColl;
}

TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met)
{
  
  double term1 = leptW.Pz() * ( leptW.Px()* v_met.Px() + leptW.Py()*v_met.Py() + pow(80.399, 2)/2.);
  
  double det = pow(leptW.Px() * v_met.Px() + leptW.Py() * v_met.Py() + pow(80.399, 2)/2., 2) - v_met.Pt()*v_met.Pt() * (leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz() );
  
  if(det<0) det=0;
  
  double term2 = leptW.E() * pow(det, 0.5);
  double denom = leptW.E()*leptW.E() - leptW.Pz()*leptW.Pz();
  double sol1 = (term1 - term2) / denom;
  //double sol2 = (term1 + term2) / denom;
  double nu_E = 0;
  
  TLorentzVector neutrino;
  
  nu_E = pow( pow(v_met.Px(),2) + pow(v_met.Py(),2) + pow(sol1,2), 0.5);//neglecting neutrino mass
  neutrino.SetPxPyPzE( v_met.Px(), v_met.Py(), sol1, nu_E);
  
  return neutrino;
  
  
}



int SMjetCalculator(std::vector<TRootPFJet*> Jets,int verb){
  int index_ = -5 ;
  if(Jets.size()>1){
    for(int iJ = 0; iJ < Jets.size()-1 ; iJ++){
      for(int kJ = 1; kJ < Jets.size(); kJ++){
        if(Jets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= Jets[kJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags())
        {
          index_ = iJ;
          // cout << "index is " << iJ << endl ;
        }
        else {
          index_ = kJ;
          // cout << "index is " << kJ << endl;
        }
      }
      
    }
    // cout << "index is " << index_ << endl;
  }
  else if(Jets.size() == 1){
    index_ = 0;
    //  cout << "index is " << index_ << endl;
    
  }
  
  
  return index_;
};






