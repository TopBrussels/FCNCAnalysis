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



//Initializing b-tag WP
// taken from https://indico.cern.ch/event/600194/contributions/2423599/attachments/1396766/2129590/kskovpenTOP20170117.pdf
Double_t workingpointvalue_Loose = 0.5426;//working points updated to 2016 ReReco BTV-POG recommendations.
Double_t workingpointvalue_Medium = 0.8484;//working points updated to 2016 BTV-POG recommendations.
Double_t workingpointvalue_Tight = 0.9535;//working points updated to 2016 BTV-POG recommendations.

std::pair <Double_t,Double_t> c_workingpointvalue_Loose(-0.48, -0.17); // reduces b -jets (cvsln cvsb)
std::pair < Double_t, Double_t > c_workingpointvalue_Medium(-0.1, -0.08); // reduces light and b jets
std::pair <Double_t, Double_t> c_workingpointvalue_Tight(0.69, -0.45); // reduces light




//What you want to do
bool synchex = false;


float lum_RunsBCDEF = 15.658183109;// /fb
float lum_RunsGH = 15.199167277;// /fb



// home made functions
std::pair <Float_t,Float_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo);
int FCNCjetCalculator(std::vector<TRootPFJet*> Jets, TLorentzVector recoZ ,int index, int verb);
int FCNCjetCalculatorCvsBTagger(std::vector<TRootPFJet*> Jets, int index, int verb);
int FCNCjetCalculatorCvsLTagger(std::vector<TRootPFJet*> Jets, int index, int verb);
int FCNCjetCalculatorCwp(std::vector<TRootPFJet*> Jets, std::vector <int> cjetindex, int index, int verb);
int SMjetCalculator(std::vector<TRootPFJet*> Jets,int verb);
double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy);
pair< vector <TLorentzVector> , vector < pair < string , int > > >  LeptonAssigner(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons);
//vector<TLorentzVector> LeptonAssignerv2(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons);

bool isVetoElectronSpring2016(TRootElectron electron);
bool isTightElectronSpring2016(TRootElectron electron);
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met);
// administration functions
string ConvertIntToString(int Number, bool pad);
string MakeTimeStamp();
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  ObjectMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedobjects);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatcherST(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatchertZq(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatcherST(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets);
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatchertZq(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets);
// members
//   bool stop_program;
double M_W  = 80.4;
double M_mu =  0.10566; // 105.66 MeV/c^2
double M_el = 0.000510999; // 0.510998910 Mev/c^2



pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   LeptonMatcherPair;
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   JetMatcherPair;
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   ObjectMatcherPair;

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
  int verbose = 0; // 0 = cout alll
  bool eventSelected = false;
  bool baseSelected = false;
  int nbTrig = 0;
  int nbBaseline = 0;
  int nbGPV = 0;
  int nbSelectedEvents = 0;
  int nbEvents = 0;
  double dataLumi = 0; //pb
  
  // to put  on here
  bool runHLT = true;
  bool applyJetLeptonCleaning = true;
  bool btagShape = false;
  bool printTrigger = false;
  bool printLeptonSF = false;
  bool applyPU = true;
  
  int maxMCParticles = -1;
  
  // to put on with agrs
  bool applyJER = false;
  bool applyJES = false;
  bool fillBtagHisto = false;
  
  
  
  
  //////////////////////////////////////////////
  /// Set up everything for local submission ////
  ///////////////////////////////////////////////
  // check the arguments passed
  if(verbose>3)
  {
    cout << " The list of arguments are: " << endl;
    for (int n_arg=1; n_arg<argc; n_arg++)
    {
      std:: cerr << "  - arg number " << n_arg << " is " << argv[n_arg] << std::endl;
    }
  }
  if(argc < 18)
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
  const int JES                 =  strtol(argv[argc-7], NULL,10);
  const int JER                 =  strtol(argv[argc-6], NULL,10);
  const int FillBtagHisto	 =  strtol(argv[argc-5], NULL,10);
  const int Usettbar = strtol(argv[argc-4], NULL,10);
  const int JobNum		  = strtol(argv[argc-3], NULL, 10);
  const int startEvent  	  = strtol(argv[argc-2], NULL, 10);
  const int endEvent		  = strtol(argv[argc-1], NULL, 10);
  
  applyJES = JES;
  applyJER = JER;
  fillBtagHisto = FillBtagHisto;
  
  string btag_dir = "BTagHistosPtEta";
  if(fillBtagHisto)
  {
    
    mkdir(btag_dir.c_str(),0777);
  }
  
  // all the files are stored from arg 11 to argc-2
  vector<string> vecfileNames;
  
  for(int args = 11; args < argc-7; args++)
  {
    if(verbose > 3){cout << "pushing back " << argv[args] << endl;}
    vecfileNames.push_back(argv[args]);
    
  }
  
  if (verbose>0)
  {
    cout << "The list of file to run over will be printed..." << endl;
    for ( int nfiles = 0; nfiles < vecfileNames.size(); nfiles++)
    {
      cout << "file number " << nfiles << " is " << vecfileNames[nfiles] << endl;
    }
  }
  /// define decays
  //
  
  string sdecay = "";
  if(Usettbar) sdecay = "ttbar";
  else sdecay = "singletop";
  cout << " --> Using the " << sdecay << "  decay <-- " << endl;
  string xmlFileName = "config/Run2TriLepton_samples.xml" ;
  dataLumi = 36000; //pb
  
  
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
  bool matching = false;
  bool istZq = false;
  if(dName.find("Data")!=string::npos || dName.find("data")!=string::npos || dName.find("DATA")!=string::npos){
    isData = true;
    cout << "running on data !!!!" << endl;
    cout << "luminosity is " << dataLumi << endl;
  }
  if( dName.find("NP_overlay_TT_FCNC")!=string::npos || dName.find("tZq")!=string::npos || dName.find("NP_overlay_ST_FCNC")!=string::npos  )
  {
    matching = true;
    cout << "WARNING: looking at mcParticles !! " << endl;
  }
  if(dName.find("tZq")!=string::npos )
  {
    istZq = true;
    
  }
  cout << "----------------------------------------" << endl;
  matching = false;
  
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
  anaEnv.ElectronCollection = "Electrons_calibratedPatElectrons";
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
  
  if(verbose > 0) cout << "Initializing trigger" << endl;
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
  BTagWeightTools *btwt = 0;
  BTagCalibrationReader * reader_csvv2;
 
  BTagCalibration *btagcalib_c;
  BTagCalibrationReader *btagreader_c;
  BTagWeightTools *btwt_c = 0;
  BTagCalibrationReader * reader_charm;
  //TFile *histoFileHandle = 0;
  // for pu
  LumiReWeighting LumiWeights;
  
  // JER / JEC
  vector<JetCorrectorParameters> vCorrParam;
  
  
  ///////////////////////////////
  //  Set up Output ROOT file  ///
  //////////////////////////////
  stringstream ss;
  ss << JobNum;
  string strJobNum = ss.str();
  string histo_dir = "NtupleMakerOutput/TriLepton_histos/";
  string histo_dirdecay = histo_dir +sdecay;
  string histo_dir_date = histo_dirdecay+"/TriLepton_histos_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
  mkdir(histo_dirdecay.c_str(),0777);
  mkdir(histo_dir_date.c_str(),0777);
  
  string rootFileName (histo_dir_date+"/FCNC_3L_"+dName+".root");
  if (strJobNum != "0")
  {
    if(verbose > 0) cout << "strJobNum is " << strJobNum << endl;
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
  float el_pt_cut =35.; // 42
  float el_eta_cut = 2.1;
  float el_iso_cone  = 0.3;
  // reliso cut fabs(eta supercluster) <= 1.479 --> 0.107587 // (fabs(eta supercluster) > 1.479 && fabs(eta supercluster) < 2.5) --> 0.113254
  // muon
  float mu_pt_cut = 25.; // 40
  float mu_eta_cut = 2.1;
  float mu_iso_cut = 0.15;
  float mu_iso_cut_loose = 0.25;
  //jets
  float jet_pt_cut = 30.;
  float jet_eta_cut = 2.4;
  
  // convert into string
  
  std::ostringstream el_pt_cut_strs, el_eta_cut_strs, mu_pt_cut_strs, mu_eta_cut_strs, mu_iso_cut_strs, jet_pt_cut_strs, jet_eta_cut_strs, mu_iso_cut_loose_strs;
  std::string el_pt_cut_str, el_eta_cut_str, mu_pt_cut_str, mu_eta_cut_str, mu_iso_cut_str, jet_pt_cut_str, jet_eta_cut_str, mu_iso_cut_loose_str;
  el_pt_cut_strs << el_pt_cut;
  el_eta_cut_strs << el_eta_cut;
  mu_pt_cut_strs << mu_pt_cut;
  mu_eta_cut_strs << mu_eta_cut;
  mu_iso_cut_strs << mu_iso_cut;
  jet_pt_cut_strs << jet_pt_cut;
  jet_eta_cut_strs << jet_eta_cut;
  mu_iso_cut_loose_strs << mu_iso_cut_loose;
  el_pt_cut_str = el_pt_cut_strs.str();
  el_eta_cut_str = el_eta_cut_strs.str();
  mu_pt_cut_str = mu_pt_cut_strs.str();
  mu_eta_cut_str = mu_eta_cut_strs.str();
  mu_iso_cut_str = mu_iso_cut_strs.str();
  jet_pt_cut_str = jet_pt_cut_strs.str();
  jet_eta_cut_str = jet_eta_cut_strs.str();
  mu_iso_cut_loose_str = mu_iso_cut_loose_strs.str();
  
  
  
  
  
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
  
  histo1D["recoZmass"] = new TH1F("recoZmass","recoZmass",500,0,500);
  histo1D["recomWt"] = new TH1F("recomWt","recomWt",150,0,150);
  if(Usettbar && !istZq) histo1D["recoFCNCTopmass"] = new TH1F("recoFCNCTopmass","recoFCNCTopmass",500,0,500);
  histo1D["recoSMTopmass"] = new TH1F("recoSMTopmass","recoSMTopmass",500,0,500);
  
  histo1D["LP"]= new TH1F("LP", "alt cos theta", 200, -1,1);
  histo1D["CosThetaWRF"]= new TH1F("CosThetaWRF", "CosTheta* in the W restframe reconstructed", 200, -1,1);
  histo1D["CosThetaWRFTRF"]=new TH1F("CosThetaWRTRF", "CosTheta* in the W and top RF reconstructed", 200, -1,1);
  histo2D["CosTheta"]= new TH2F("CosTheta", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);
  
  if(matching){
    histo1D["matchedZmass"] = new TH1F("matchedZmass","matchedZmass",500,0,500);
    if(Usettbar && !istZq) histo1D["matchedFCNCTopmass"] = new TH1F("matchedFCNCTopmass","matchedFCNCTopmass",500,0,500);
    histo1D["matchedSMTopmass"] = new TH1F("matchedSMTopmass","matchedSMTopmass",500,0,500);
    
    histo1D["CosThetaWRF_gen"]= new TH1F("CosThetaWRF_gen", "CosTheta* in the W restframe reconstructed", 200, -1,1);
    histo1D["CosThetaWRFTRF_gen"]=new TH1F("CosThetaWRTRF_gen", "CosTheta* in the W and top RF reconstructed", 200, -1,1);
    histo2D["CosTheta_gen"]= new TH2F("CosTheta_gen", "CosTheta* in the W RF vs W RF en Top RF", 200, -1,1, 200, -1,1);

    
    histo1D["Topmass_Wb"] = new TH1F("Topmass_Wb","Topmass_Wb",500,0,500);
    histo1D["pt_Wb"] = new TH1F("pt_Wb","pt_Wb",500,0,500);
    histo1D["eta_Wb"]= new TH1F("eta_Wb","eta_Wb",30,-3,3);
    histo1D["phi_Wb"]= new TH1F("phi_Wb","phi_Wb",32,-3.2,3.2);
    
    histo1D["Topmass_lvb"] = new TH1F("Topmass_lvb","Topmass_lvb",500,0,500);
    histo1D["pt_lvb"] = new TH1F("pt_lvb","pt_lvb",500,0,500);
    histo1D["eta_lvb"]= new TH1F("eta_lvb","eta_lvb",30,-3,3);
    histo1D["phi_lvb"]= new TH1F("phi_lvb","phi_lvb",32,-3.2,3.2);
    
    
    histo1D["Topmass_tq"] = new TH1F("Topmass_tq","Topmass_tq",500,0,500);
    histo1D["pt_tq"] = new TH1F("pt_tq","pt_tq",500,0,500);
    histo1D["eta_tq"]= new TH1F("eta_tq","eta_tq",30,-3,3);
    histo1D["phi_tq"]= new TH1F("phi_tq","phi_tq",32,-3.2,3.2);
    
    histo1D["dR_lvb"] =  new TH1F("DRlvb", "dR lvb", 500,0, 5);
    histo1D["dR_Wb"] =  new TH1F("DRWb", "dR Wb", 500,0, 5);
    histo1D["dR_Wbtop"]  = new TH1F("DRWbtop", "dR Wb-top", 500,0, 5);
    histo1D["dPhi_lvb"]  = new TH1F("DPhiWb", "dPhi Wb", 140,-7, 7);
    histo1D["dPhi_Wb"]  = new TH1F("DPhilvb", "dPhi lvb", 140,-7, 7);
    histo1D["dPhi_Wbtop"]  = new TH1F("DPhiWbtop", "dPhi Wb-top", 140,-7, 7);
    
    histo2D["Topmass_top_Wb"] = new TH2F("Topmass_top_Wb","Topmass:Wb:t",500,0,500,500,0,500);
    histo2D["pt_top_Wb"] = new TH2F("pt_top_Wb","pt:Wb:t",500,0,500,500,0,500);
    histo2D["eta_top_Wb"]= new TH2F("eta_top_Wb","eta:Wb:t",30,-3,3,30,-3,3);
    histo2D["phi_top_Wb"]= new TH2F("phi_top_Wb","phi:Wb:t",32,-3.2,3.2,32,-3.2,3.2);
    
    
    
    
    if(Usettbar && !istZq){
      histo1D["Topmass_Zq"] = new TH1F("Topmass_Zq","Topmass_Zq",500,0,500);
      histo1D["pt_Zq"] = new TH1F("pt_Zq","pt_Zq",500,0,500);
      histo1D["eta_Zq"]= new TH1F("eta_Zq","eta_Zq",30,-3,3);
      histo1D["phi_Zq"]= new TH1F("phi_Zq","phi_Zq",32,-3.2,3.2);
      
      histo1D["Topmass_llq"] = new TH1F("Topmass_llq","Topmass_llq",500,0,500);
      histo1D["pt_llq"] = new TH1F("pt_llq","pt_llq",500,0,500);
      histo1D["eta_llq"]= new TH1F("eta_llq","eta_llq",30,-3,3);
      histo1D["phi_llq"]= new TH1F("phi_llq","phi_llq",32,-3.2,3.2);
      
      histo1D["Topmass_fcnctq"] = new TH1F("Topmass_fcnctq","Topmass_fcnctq",500,0,500);
      histo1D["pt_fcnctq"] = new TH1F("pt_fcnctq","pt_fcnctq",500,0,500);
      histo1D["eta_fcnctq"]= new TH1F("eta_fcnctq","eta_fcnctq",30,-3,3);
      histo1D["phi_fcnctq"]= new TH1F("phi_fcnctq","phi_fcnctq",32,-3.2,3.2);
      
      histo1D["mass_FCNCq"] = new TH1F("mass_FCNCq","mass_fcnc q",100,0,10);
      histo1D["pt_FCNCq"] = new TH1F("pt_FCNCq","pt_fcnc q",500,0,500);
      histo1D["eta_FCNCq"]= new TH1F("eta_FCNCq","eta_fcnc q",30,-3,3);
      histo1D["phi_FCNCq"]= new TH1F("phi_FCNCq","phi_fcnc q",32,-3.2,3.2);
      
      histo1D["dR_llq"] =  new TH1F("DRllq", "dR llq", 500,0, 5);
      histo1D["dR_Zq"] =  new TH1F("DRZq", "dR Zq", 500,0, 5);
      histo1D["dR_Zqtop"]  = new TH1F("DRZqtop", "dR Zq-top", 500,0, 5);
      histo1D["dPhi_llq"]  = new TH1F("DPhiZq", "dPhi Zq", 140,-7, 7);
      histo1D["dPhi_Zq"]  = new TH1F("DPhillq", "dPhi llq", 140,-7, 7);
      histo1D["dPhi_Zqtop"]  = new TH1F("DPhiZqtop", "dPhi Zq-top", 140,-7, 7);
      
      histo2D["Topmass_top_Zq"] = new TH2F("Topmass_top_Zq","Topmass:Zq:t",500,0,500,500,0,500);
      histo2D["pt_top_Zq"] = new TH2F("pt_top_Zq","pt:Zq:t",500,0,500,500,0,500);
      histo2D["eta_top_Zq"]= new TH2F("eta_top_Zq","eta:Zq:t",30,-3,3,30,-3,3);
      histo2D["phi_top_Zq"]= new TH2F("phi_top_Zq","phi:Zq:t",32,-3.2,3.2,32,-3.2,3.2);
    }
    histo1D["mass_lep1"]                                  = new TH1F("mass_lep1","mass_lep1",250,0,0.5);
    histo1D["mass_lep2"]                                  = new TH1F("mass_lep2","mass_lep2",250,0,0.5);
    histo1D["Zmass_Zlep"]                                  = new TH1F("Zmass_Zlep","Zmass_lep",200,0,200);
    histo1D["Zmass_Zbos"]                                  = new TH1F("Zmass_Zbos","Zmass_Zbos",200,0,200);
    histo2D["mass_lep"]            = new TH2F("mass_lep", "mass lep;mass lep1;mass lep2", 250,0,0.5,250,0,0.5 );
    histo2D["Zmass_Zbos_Zlep"]            = new TH2F("Zmass", "Zmass;Zmass_lep;Zmass_Zboson", 200,0,200,200,0,200 );
    
    histo1D["pt_lep1"]          = new TH1F("ptlep1", "Pt lep 1", 200,0,400);
    histo1D["pt_lep2"]          = new TH1F("ptlep2", "Pt lep 2", 200,0,400);
    histo1D["pt_Zbos"]          = new TH1F("ptZbos", "Pt Zboson", 200,0,400);
    histo1D["pt_Zlep"]          = new TH1F("ptZlep", "Pt Zlep", 200,0,400);
    histo2D["pt_lep"]          = new TH2F("ptlep", "Pt leptons;Pt lep 1;Pt lep 2", 200,0,400, 200,0,400);
    histo2D["pt_Z"]          = new TH2F("ptZ", "Pt Z;Pt Zlep;Pt Zbos", 200,0,400, 200,0,400);
    
    histo1D["phi_lep1"]          = new TH1F("philep1", "phi lep 1",32,-3.2,3.2);
    histo1D["phi_lep2"]          = new TH1F("philep2", "phi lep 2",32,-3.2,3.2);
    histo1D["phi_Zbos"]          = new TH1F("phiZbos", "phi Zboson",32,-3.2,3.2);
    histo1D["phi_Zlep"]          = new TH1F("phiZlep", "phi Zlep",32,-3.2,3.2);
    histo2D["phi_lep"]          = new TH2F("philep", "phi lep;phi lep 1;phi lep 2",32,-3.2,3.2,32,-3.2,3.2);
    histo2D["phi_Z"]          = new TH2F("phiZ", "phi Z;phi Zlep;phi Zbos",32,-3.2,3.2,32,-3.2,3.2);
    
    histo1D["eta_lep1"]          = new TH1F("etalep1", "eta lep 1", 30,-3,3);
    histo1D["eta_lep2"]          = new TH1F("etalep2", "eta lep 2", 30,-3,3);
    histo1D["eta_Zbos"]          = new TH1F("etaZbos", "eta Zboson", 30,-3,3);
    histo1D["eta_Zlep"]          = new TH1F("etaZlep", "eta Zlep", 30,-3,3);
    histo2D["eta_lep"]          = new TH2F("etalep", "eta lep;eta lep 1;eta lep 2", 30,-3,3, 30,-3,3);
    histo2D["eta_Z"]          = new TH2F("etaZ", "eta Z;eta Zlep;eta Zbos", 30,-3,3, 30,-3,3);
    
    histo1D["dR_lep"]          = new TH1F("DRlep", "dR lep", 500,0, 5);
    histo1D["dPhi_lep"]          = new TH1F("DPhilep", "dPhi lep", 140,-7, 7);
  }
  histo1D["nRecoLeptons"] = new TH1F("nRecoLeptons","nRecoLeptons", 10,-0.5,9.5);
  histo1D["nIniRecoLeptons"] = new TH1F("nIniRecoLeptons","nIniRecoLeptons", 10,-0.5,9.5);
  histo1D["nRecoElectrons"] = new TH1F("nRecoElectrons","nRecoElectrons", 10,-0.5,9.5);
  histo1D["nIniRecoElectrons"] = new TH1F("nIniRecoElectrons","nIniRecoElectrons", 10,-0.5,9.5);
  histo1D["nRecoMuons"] = new TH1F("nRecoMuons","nRecoMuons", 10,-0.5,9.5);
  histo1D["nIniRecoMuons"] = new TH1F("nIniRecoMuons","nIniRecoMuons", 10,-0.5,9.5);
  histo2D["nIniRecoMuonsnRecoMuons"] =  new TH2F("nIniRecoMuonsnRecoMuons", "nRecoMuons;initial;selected", 10,-0.5,9.5, 10,-0.5,9.5);
  histo2D["nIniRecoElectronsnRecoElectrons"] =  new TH2F("nIniRecoElectronsnRecoElectrons", "nRecoElectrons;initial;selected", 10,-0.5,9.5, 10,-0.5,9.5);
  histo2D["nIniRecoLeptonsnRecoLeptons"] =  new TH2F("nIniRecoLeptonsnRecoLeptons", "nRecoLeptons;initial;selected", 10,-0.5,9.5, 10,-0.5,9.5);
  
  //  histo1D["mc_nZ"] = new TH1F("mc_nZ","mc_nZ", 10,-0.5,9.5);
  histo1D["mc_nZEl"]= new TH1F("mc_nZEl","mc_nZEl", 10,-0.5,9.5);
  histo1D["mc_nZMu"]= new TH1F("mc_nZMu","mc_nZMu", 10,-0.5,9.5);
  histo1D["mc_nZLep"]= new TH1F("mc_nZLep","mc_nZLep", 10,-0.5,9.5);
  histo1D["mc_nWEl"]= new TH1F("mc_nWEl","mc_nWEl", 10,-0.5,9.5);
  histo1D["mc_nWMu"]= new TH1F("mc_nWMu","mc_nWMu", 10,-0.5,9.5);
  histo1D["mc_nWLep"]= new TH1F("mc_nWLep","mc_nWLep", 10,-0.5,9.5);
  /* histo1D["mc_nW"]= new TH1F("mc_nW","mc_nW", 10,-0.5,9.5);
   histo1D["mc_nTW"]= new TH1F("mc_nTW","mc_nTW", 10,-0.5,9.5);
   histo1D["mc_nTZ"]= new TH1F("mc_nTZ","mc_nTZ", 10,-0.5,9.5);
   histo1D["mc_nTWelectrons"]= new TH1F("mc_nTWelectrons","mc_nTW elec", 10,-0.5,9.5);
   histo1D["mc_nTZelectrons"]= new TH1F("mc_nTZelectrins","mc_nTZ elec", 10,-0.5,9.5);
   histo1D["mc_nTWmuons"]= new TH1F("mc_nTWmuons","mc_nTW mu", 10,-0.5,9.5);
   histo1D["mc_nTZmuons"]= new TH1F("mc_nTZmuons","mc_nTZ mu", 10,-0.5,9.5);
   
   histo2D["nWboson"] =  new TH2F("nWboson", "Wboson;nW;nLep", 10,-0.5,9.5, 10,-0.5,9.5);
   histo2D["nZboson"] =  new TH2F("nZboson", "Zboson;nZ;nLep", 10,-0.5,9.5, 10,-0.5,9.5);
   */  histo2D["nZbosonnWboson"] =  new TH2F("nZbosonnWboson", "Boson;nZlep;nWlep", 10,-0.5,9.5, 10,-0.5,9.5);
  
  
  
  
  /// LUMIREWEIGHING
  
  //LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring16MiniAODv2-Asympt.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting//pileup_2016Data80X_Run273158-276811Cert.root", "pileup", "pileup");
  
  LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
  
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
    string histfile = "BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_comb_central_" + sdecay + ".root";
    string histreadfile = "BTagHistosPtEta/Merged/"+daName+ "_comb_central_"+ sdecay + ".root";
    if(!isData && !btagShape)
    {
      // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
      //	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_com@
      //     btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_80X_ichep_incl_ChangedTo_mujets.csv");
      btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
      btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_LOOSE, "comb","central");
      
      
      
      if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
      {
        cout << "filling btag histo's" << endl;
        //histoFileHandle = TFile::Open(histfile.c_str(), "UPDATE");
        btwt = new BTagWeightTools(btagreader,histfile.c_str(),false,30,999,2.4);
      }
      else if(!fillBtagHisto)
      {
        cout << "reading btag histo's from " << histreadfile.c_str() << endl;
       // histoFileHandle = TFile::Open("BTagHistosPtEta/HistosPtEta_"+daName+ "_comb_central.root", "READ");
        //histoFileHandle = TFile::Open(histreadfile.c_str());
        
        btwt = new BTagWeightTools(btagreader,histreadfile.c_str(),true,30,999,2.4);
        //btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
      }
      
      
    }
    else if(!isData) // NEEDS TO BE CHECKED FOR 80X
    {
      cout << " WARNING: btag shape is used" << endl;
      BTagCalibration calib_csvv2("csvv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
      reader_csvv2 = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                               BTagEntry::OP_RESHAPING, // operating point
                                               "iterativefit", // measurement type
                                               "central"); // systematics type  --> depending on JES up/Down andother reader is needed
      
      
    }
    
    if(verbose>1) cout << "btag done" << endl;
    
    
    MuonSFWeight* muonSFWeightID_BCDEF;
    MuonSFWeight* muonSFWeightID_GH;
    MuonSFWeight* muonSFWeightIso_BCDEF;
    MuonSFWeight* muonSFWeightIso_GH;
    MuonSFWeight* muonSFWeightTrig_BCDEF;
    MuonSFWeight* muonSFWeightTrig_GH;
    
    
    ElectronSFWeight* electronSFWeightID;
    ElectronSFWeight* electronSFWeightReco;
    
    muonSFWeightID_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
    muonSFWeightID_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonID_EfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
    muonSFWeightIso_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonIso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
    muonSFWeightIso_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonIso_EfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
    // muonSFWeightTrig_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsBCDEF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
    //  muonSFWeightTrig_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsGH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
    
    
    if(verbose>1) cout << "muon SF loaded" << endl;
    
    electronSFWeightID = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/Moriond17/egammaEffi.txt_EGM2D_CutBasedTightID.root","EGamma_SF2D",true,false,false);
    electronSFWeightReco = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/Moriond17/egammaEffi.txt_EGM2D_RecoEff.root","EGamma_SF2D",true,false,false);
    
    if(verbose >1) cout << "electron SF loaded " << endl;
    
    vCorrParam.clear();
    JetCorrectionUncertainty *jecUnc;
    
    if(dName.find("Data_Run2016B")!=string::npos || dName.find("Data_Run2016C")!=string::npos || dName.find("Data_Run2016D")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data_Run2016E")!=string::npos || dName.find("Data_Run2016F")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data_Run2016G")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data_Run2016H")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data")==string::npos ) // for Summer16 MC
    {
      //cout << "loading MC JEC / JER" << endl;
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt");
    }
    if(verbose>1) cout << "jec and jer loaded"<< endl;
    
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    if(verbose>1) cout << "jet tools initialised " << endl;
    
    ////////////////////////////////////////////////////////////
    // Setup Date string and nTuple for output
    ///////////////////////////////////////////////////////////
    
    string channel_dir = "NtupleMakerOutput/Ntuples/"+sdecay+"/" ;
    string date_dir = channel_dir+ "/" + dateString +"/";
    mkdir(channel_dir.c_str(),0777);
    mkdir(date_dir.c_str(),0777);
    
    
    string Ntupname = date_dir +"FCNC_3L_" + dName + "_"+  strJobNum + ".root";
    cout << "Ntuple " << Ntupname << " created " << endl;
    
    TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
    tupfile->cd();
    cout << "outputfile " << Ntupname.c_str()<< " made and accessed " <<endl ;
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
    int PassedTrigger;
    vector <int> cutstep;
    vector <int> cutstep_eee;
    vector <int> cutstep_eeu;
    vector <int> cutstep_uuu;
    vector <int> cutstep_uue;
    vector <string> cutstep_string;
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
    Int_t METon;
    Double_t WPb_L;
    Double_t WPb_M;
    Double_t WPb_T;
    Double_t WPc_CvsB_Loose;
    Double_t WPc_CvsB_Medium;
    Double_t WPc_CvsB_Tight;
    Double_t WPc_CvsL_Loose;
    Double_t WPc_CvsL_Medium;
    Double_t WPc_CvsL_Tight;
    Int_t PassedMET;
    Int_t channelInt;
    Int_t signalInt;
    
    Float_t pt_electron_1;
    Float_t pt_electron_2;
    Float_t pt_electron_3;
    Float_t pt_muon_1;
    Float_t pt_muon_2;
    Float_t pt_muon_3;
    Float_t pt_jet_1;
    Float_t pt_jet_2;
    Float_t pt_jet_3;
    
    
    Int_t nLeptons;
    // variables for electrons
    Int_t nElectrons;
    Float_t pt_electron[10];
    Float_t phi_electron[10];
    Float_t eta_electron[10];
    Float_t eta_superCluster_electron[10];
    Float_t E_electron[10];
    Float_t d0_electron[10];
    Float_t d0BeamSpot_electron[10];
    Float_t chargedHadronIso_electron[10];
    Float_t neutralHadronIso_electron[10];
    Float_t photonIso_electron[10];
    Float_t pfIso_electron[10];
    Int_t charge_electron[10];
    
    Float_t sigmaIEtaIEta_electron[10];
    Float_t deltaEtaIn_electron[10];
    Float_t deltaPhiIn_electron[10];
    Float_t hadronicOverEm_electron[10];
    Int_t missingHits_electron[10];
    Bool_t passConversion_electron[10];
    Bool_t isId_electron[10];
    Bool_t isIso_electron[10];
    
    Bool_t isEBEEGap[10];
    
    //variable for muons
    Int_t nMuons;
    Float_t pt_muon[10];
    Float_t phi_muon[10];
    Float_t eta_muon[10];
    Float_t E_muon[10];
    Float_t d0_muon[10];
    Float_t d0BeamSpot_muon[10];
    Float_t chargedHadronIso_muon[10];
    Float_t neutralHadronIso_muon[10];
    Float_t photonIso_muon[10];
    Float_t relIso_muon[10];
    Bool_t isId_muon[10];
    Bool_t isIso_muon[10];
    Float_t pfIso_muon[10];
    Int_t charge_muon[10];
    
    //variable for jets
    Int_t nJets;
    Int_t nJets_CSVL;
    Int_t nJets_CSVM;
    Int_t nJets_CSVT;
    Int_t nJets_CharmL;
    Int_t nJets_CharmM;
    Int_t nJets_CharmT;
    Int_t nJets_nonCSVL;
    Int_t nJets_nonCSVM;
    Int_t nJets_nonCSVT;
    Int_t nJets_nonCharmL;
    Int_t nJets_nonCharmM;
    Int_t nJets_nonCharmT;
    Int_t nJets_nonCharmLCSVL;
    Int_t nJets_nonCharmLCSVM;
    Int_t nJets_nonCharmLCSVT;
    Int_t nJets_nonCharmMCSVL;
    Int_t nJets_nonCharmMCSVM;
    Int_t nJets_nonCharmMCSVT;
    Int_t nJets_nonCharmTCSVL;
    Int_t nJets_nonCharmTCSVM;
    Int_t nJets_nonCharmTCSVT;
    Int_t nJets_nonCharmLnonCSVL;
    Int_t nJets_nonCharmLnonCSVM;
    Int_t nJets_nonCharmLnonCSVT;
    Int_t nJets_nonCharmMnonCSVL;
    Int_t nJets_nonCharmMnonCSVM;
    Int_t nJets_nonCharmMnonCSVT;
    Int_t nJets_nonCharmTnonCSVL;
    Int_t nJets_nonCharmTnonCSVM;
    Int_t nJets_nonCharmTnonCSVT;
    Int_t nJets_CharmLnonCSVL;
    Int_t nJets_CharmLnonCSVM;
    Int_t nJets_CharmLnonCSVT;
    Int_t nJets_CharmMnonCSVL;
    Int_t nJets_CharmMnonCSVM;
    Int_t nJets_CharmMnonCSVT;
    Int_t nJets_CharmTnonCSVL;
    Int_t nJets_CharmTnonCSVM;
    Int_t nJets_CharmTnonCSVT;
    Int_t nJets_CharmLCSVL;
    Int_t nJets_CharmLCSVM;
    Int_t nJets_CharmLCSVT;
    Int_t nJets_CharmMCSVL;
    Int_t nJets_CharmMCSVM;
    Int_t nJets_CharmMCSVT;
    Int_t nJets_CharmTCSVL;
    Int_t nJets_CharmTCSVM;
    Int_t nJets_CharmTCSVT;
    Int_t charmL_jet[20];
    Int_t charmM_jet[20];
    Int_t charmT_jet[20];
    Int_t btagL_jet[20];
     Int_t btagM_jet[20];
     Int_t btagT_jet[20];
    Float_t pt_jet[20];
    Float_t px_jet[20];
    Float_t py_jet[20];
    Float_t pz_jet[20];
    Float_t phi_jet[20];
    Float_t eta_jet[20];
    Float_t E_jet[20];
    Int_t charge_jet[20];
    Float_t bdisc_jet[20];
    Float_t cdiscCvsL_jet[20];
    Float_t cdiscCvsL_jet_1;
    Float_t cdiscCvsB_jet_1;
    Float_t cdiscCvsB_jet[20];
    double orig_jet_px;
    double orig_jet_py ;
    double orig_jet_pt;
    double corrected_jet_px ;
    double corrected_jet_py ;
    double corrected_jet_pt ;
    double jet_pt_check;
    
    // variables for Zboson
    Float_t Zboson_M;
    Float_t Zboson_Px;
    Float_t Zboson_Py;
    Float_t Zboson_Pz;
    Float_t Zboson_Energy;
    
    Float_t Zboson2_M;
    Float_t Zboson2_Px;
    Float_t Zboson2_Py;
    Float_t Zboson2_Pz;
    Float_t Zboson2_Energy;
    
    // met
    Float_t met_Pt;
    Float_t met_Ptbf;
    Float_t met_Px;
    Float_t met_Py;
    Float_t met_Pz;
    Float_t met_Phi;
    Float_t met_Eta;
    
    
    
    double orig_met_px;
    double orig_met_py ;
    double orig_met_pt;
    double corrected_met_px ;
    double corrected_met_py ;
    double corrected_met_pt ;
    
    Float_t mWt;
    Float_t FCNCtop_M;
    Float_t SMtop_M;
    Float_t cjet_Pt;
    Float_t mlb;
    Float_t dRWlepc;
    Float_t dRZb;
    Float_t dRZc;
    Float_t dRWlepb;
    Float_t dRSMFCNCtop;
    Float_t dPhiSMFCNCtop;
    Float_t dPhiWlepb;
    Float_t dPhiWlepc;
    Float_t dPhiZb;
    Float_t dPhiZc;
    
    Float_t FCNCtop_M_tagger;
    Float_t cjet_Pt_tagger;
    Float_t dRWlepc_tagger;
    
    Float_t dRZc_tagger;
    
    Float_t dRSMFCNCtop_tagger;
    Float_t dPhiSMFCNCtop_tagger;
    
    Float_t dPhiWlepc_tagger;
    
    Float_t dPhiZc_tagger;
    
    
    Float_t Wlep_Charge;
    Float_t Wlep_Pt;
    Float_t Wlep_Phi;
    Float_t Wlep_Eta;
    Float_t Zboson_Eta;
    Float_t Zboson_Phi;
    Float_t Zboson_Pt;
    Float_t charge_asym;
    Float_t dPhiZWlep;
    Float_t dRZWlep;
    Float_t SMtop_Eta;
    Float_t SMtop_Phi;
    Float_t SMtop_Pt;
    Float_t dPhiZMET;
    Float_t dPhiZSMtop;
    Float_t dRZSMtop;
    Float_t TotalPt;
    Float_t TotalHt;
    Float_t TotalInvMass;
    Float_t cdiscCvsB_jet_2;
    Float_t cdiscCvsL_jet_2;
    Float_t bdiscCSVv2_jet_2;
    Float_t bdiscCSVv2_jet_1;
    Float_t CosTheta;
    Float_t CosTheta_alt;
    Float_t LP;
    
    
    
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
    Double_t mc_Zmass;
    Double_t reco_Zmass;
    int nRecoLeptons;
    int nRecoElectrons;
    int nRecoMuons;
    int nIniRecoLeptons;
    int nIniRecoElectrons;
    int nIniRecoMuons;
    
    
    
    
    
    
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
    int nTrigg_eee;
    int n3lep_eee;
    int nVetoMu_eee;
    int nVetoEl_eee;
    int nOS_eee;
    int nZmass_eee;
    int nJet_eee;
    int nBJet_eee;
    int nMWT_eee;
    int nSMtop_eee;
    int nMET_eee;
    int nTrigg_eeu;
    int n3lep_eeu;
    int nVetoMu_eeu;
    int nVetoEl_eeu;
    int nOS_eeu;
    int nZmass_eeu;
    int nJet_eeu;
    int nBJet_eeu;
    int nMWT_eeu;
    int nSMtop_eeu;
    int nMET_eeu;
    int nTrigg_uuu;
    int n3lep_uuu;
    int nVetoMu_uuu;
    int nVetoEl_uuu;
    int nOS_uuu;
    int nZmass_uuu;
    int nJet_uuu;
    int nBJet_uuu;
    int nMWT_uuu;
    int nSMtop_uuu;
    int nMET_uuu;
    int nTrigg_uue;
    int n3lep_uue;
    int nVetoMu_uue;
    int nVetoEl_uue;
    int nOS_uue;
    int nZmass_uue;
    int nJet_uue;
    int nBJet_uue;
    int nMWT_uue;
    int nSMtop_uue;
    int nMET_uue;
    int nEvPassed;
    double xsec;
    int matchedZlep_1;
    int matchedZlep_0;
    int matchedWlep;
    int matchedBjet;
    int matchedEvents_Bjet;
    int matchedEvents_Wlep;
    int matchedEvents_Zlep;
    int matchedCjet_CvsBtagger;
    int matchedEvents_Cjet_CvsBtagger;
    int int_eventForCjetmatching_CvsBtagger;
    int int_eventForCjetmatchingmatched_CvsBtagger;
    int matchedCjet_CvsLtagger;
    int matchedEvents_Cjet_CvsLtagger;
    int int_eventForCjetmatching_CvsLtagger;
    int int_eventForCjetmatchingmatched_CvsLtagger;
    int matchedCjet;
    int matchedEvents_Cjet;
    int int_eventForCjetmatching;
    int int_eventForCjetmatchingmatched;
    
    int matchedCjet_Cloose;
    int matchedEvents_Cjet_Cloose;
    int int_eventForCjetmatching_Cloose;
    int int_eventForCjetmatchingmatched_Cloose;
    int matchedCjet_Cmedium;
    int matchedEvents_Cjet_Cmedium;
    int int_eventForCjetmatching_Cmedium;
    int int_eventForCjetmatchingmatched_Cmedium;
    int matchedCjet_Ctight;
    int matchedEvents_Cjet_Ctight;
    int int_eventForCjetmatching_Ctight;
    int int_eventForCjetmatchingmatched_Ctight;
    bool eventForCjetmatchingmatched_Cloose;
    bool eventForCjetmatching_Cloose;
    bool eventForCjetmatchingmatched_Cmedium;
    bool eventForCjetmatching_Cmedium;
    bool eventForCjetmatchingmatched_Ctight;
    bool eventForCjetmatching_Ctight;
    
    int int_eventForWlepmatching;
    int int_eventForBjetmatching;
    int int_eventForBjetmatchingmatched;
    int int_eventForWlepmatchingmatched;
    int int_eventForZlepmatching;
    int int_eventForZlepmatchingmatched0;
    int int_eventForZlepmatchingmatched1;
    bool eventForWlepmatchingmatched;
    bool eventForBjetmatchingmatched;
    bool eventForBjetmatching;
    bool eventForZlepmatching;
    bool eventForZlepmatchingmatched0;
    bool eventForZlepmatchingmatched1;
    bool eventForWlepmatching;
    bool eventForCjetmatchingmatched;
    bool eventForCjetmatching;
    bool eventForCjetmatchingmatched_CvsBtagger;
    bool eventForCjetmatching_CvsBtagger;
    bool eventForCjetmatchingmatched_CvsLtagger;
    bool eventForCjetmatching_CvsLtagger;
    globalTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
    globalTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
    globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
    globalTree->Branch("nEv" , &nEv, "nEv/I");
    globalTree->Branch("matchedZlep_1" , &matchedZlep_1, "matchedZlep_1/I");
    globalTree->Branch("matchedZlep_0" , &matchedZlep_0, "matchedZlep_0/I");
    globalTree->Branch("matchedWlep" , &matchedWlep, "matchedWlep/I");
    globalTree->Branch("matchedEvents_Wlep" , &matchedEvents_Wlep, "matchedEvents_Wlep/I");
    globalTree->Branch("matchedBjet" , &matchedBjet, "matchedBjet/I");
    globalTree->Branch("matchedCjet" , &matchedCjet, "matchedCjet/I");
    globalTree->Branch("matchedCjet_Cloose" , &matchedCjet_Cloose, "matchedCjet_Cloose/I");
    globalTree->Branch("matchedCjet_Cmedium" , &matchedCjet_Cmedium, "matchedCjet_Cmedium/I");
    globalTree->Branch("matchedCjet_Ctight" , &matchedCjet_Ctight, "matchedCjet_Ctight/I");
    globalTree->Branch("matchedCjet_CvsBtagger" , &matchedCjet_CvsBtagger, "matchedCjet_CvsBtagger/I");
    globalTree->Branch("matchedCjet_CvsLtagger" , &matchedCjet_CvsLtagger, "matchedCjet_CvsLtagger/I");
    globalTree->Branch("matchedEvents_Bjet" , &matchedEvents_Bjet, "matchedEvents_Bjet/I");
    globalTree->Branch("matchedEvents_Cjet" , &matchedEvents_Cjet, "matchedEvents_Cjet/I");
    globalTree->Branch("matchedEvents_Cjet_Cloose" , &matchedEvents_Cjet_Cloose, "matchedEvents_Cjet_Cloose/I");
    globalTree->Branch("matchedEvents_Cjet_Cmedium" , &matchedEvents_Cjet_Cmedium, "matchedEvents_Cjet_Cmedium/I");
    globalTree->Branch("matchedEvents_Cjet_Ctight" , &matchedEvents_Cjet_Ctight, "matchedEvents_Cjet_Ctight/I");
    globalTree->Branch("matchedEvents_Cjet_CvsBtagger" , &matchedEvents_Cjet_CvsBtagger, "matchedEvents_Cjet_CvsBtagger/I");
    globalTree->Branch("matchedEvents_Cjet_CvsLtagger" , &matchedEvents_Cjet_CvsLtagger, "matchedEvents_Cjet_CvsLtagger/I");
    globalTree->Branch("matchedEvents_Zlep" , &matchedEvents_Zlep, "matchedEvents_Zlep/I");
    globalTree->Branch("sumW", &sumW, "sumW/I");
    globalTree->Branch("nCuts",&nCuts, "nCuts/I");
    globalTree->Branch("cutstep_string", &cutstep_string, "cutstep_string");
    globalTree->Branch("cutstep",&cutstep,"cutstep");
    globalTree->Branch("cutstep_eee",&cutstep_eee,"cutstep_eee");
    globalTree->Branch("cutstep_eeu",&cutstep_eeu,"cutstep_eeu");
    globalTree->Branch("cutstep_uuu",&cutstep_uuu,"cutstep_uuu");
    globalTree->Branch("cutstep_uue",&cutstep_uue,"cutstep_uue");
    globalTree->Branch("JERon",&JERon,"JERon/I");
    globalTree->Branch("JESon", &JESon, "JESon/I");
    globalTree->Branch("METon", &METon, "METon/I");
    globalTree->Branch("WPb_L", &WPb_L, "WPb_L/D");
    globalTree->Branch("WPb_M", &WPb_M, "WPb_M/D");
    globalTree->Branch("WPb_T", &WPb_T, "WPb_T/D");
    if(Usettbar){
      globalTree->Branch("WPc_CvsB_Loose", &WPc_CvsB_Loose, "WPc_CvsB_Loose/D");
      globalTree->Branch("WPc_CvsB_Medium", &WPc_CvsB_Medium, "WPc_CvsB_Medium/D");
      globalTree->Branch("WPc_CvsB_Tight", &WPc_CvsB_Tight, "WPc_CvsB_Tight/D");
      globalTree->Branch("WPc_CvsL_Loose", &WPc_CvsL_Loose, "WPc_CvsL_Loose/D");
      globalTree->Branch("WPc_CvsL_Medium", &WPc_CvsL_Medium, "WPc_CvsL_Medium/D");
      globalTree->Branch("WPc_CvsL_Tight", &WPc_CvsL_Tight, "WPc_CvsL_Tight/D");
    }
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
    globalTree->Branch("nEvPassed",&nEvPassed, "nEvPassed/I");
    globalTree->Branch("xsec",&xsec,"xsec/D");
    globalTree->Branch("nTrigg_eee",&nTrigg_eee, "nTrigg_eee/I");
    globalTree->Branch("n3lep_eee",&n3lep_eee, "n3lep_eee/I");
    globalTree->Branch("nVetoMu_eee",&nVetoMu_eee, "nVetoMu_eee/I");
    globalTree->Branch("nVetoEl_eee",&nVetoEl_eee, "nVetoEl_eee/I");
    globalTree->Branch("nOS_eee",&nOS_eee, "nOS_eee/I");
    globalTree->Branch("nZmass_eee",&nZmass_eee, "nZmass_eee/I");
    globalTree->Branch("nJet_eee",&nJet_eee, "nJet_eee/I");
    globalTree->Branch("nBJet_eee",&nBJet_eee, "nBJet_eee/I");
    globalTree->Branch("nMWT_eee",&nMWT_eee, "nMWT_eee/I");
    globalTree->Branch("nSMtop_eee",&nSMtop_eee, "nSMtop_eee/I");
    globalTree->Branch("nMET_eee",&nMET_eee, "nMET_eee/I");
    globalTree->Branch("nTrigg_eeu",&nTrigg_eeu, "nTrigg_eeu/I");
    globalTree->Branch("n3lep_eeu",&n3lep_eeu, "n3lep_eeu/I");
    globalTree->Branch("nVetoMu_eeu",&nVetoMu_eeu, "nVetoMu_eeu/I");
    globalTree->Branch("nVetoEl_eeu",&nVetoEl_eeu, "nVetoEl_eeu/I");
    globalTree->Branch("nOS_eeu",&nOS_eeu, "nOS_eeu/I");
    globalTree->Branch("nZmass_eeu",&nZmass_eeu, "nZmass_eeu/I");
    globalTree->Branch("nJet_eeu",&nJet_eeu, "nJet_eeu/I");
    globalTree->Branch("nBJet_eeu",&nBJet_eeu, "nBJet_eeu/I");
    globalTree->Branch("nMWT_eeu",&nMWT_eeu, "nMWT_eeu/I");
    globalTree->Branch("nSMtop_eeu",&nSMtop_eeu, "nSMtop_eeu/I");
    globalTree->Branch("nMET_eeu",&nMET_eeu, "nMET_eeu/I");
    globalTree->Branch("nTrigg_uuu",&nTrigg_uuu, "nTrigg_uuu/I");
    globalTree->Branch("n3lep_uuu",&n3lep_uuu, "n3lep_uuu/I");
    globalTree->Branch("nVetoMu_uuu",&nVetoMu_uuu, "nVetoMu_uuu/I");
    globalTree->Branch("nVetoEl_uuu",&nVetoEl_uuu, "nVetoEl_uuu/I");
    globalTree->Branch("nOS_uuu",&nOS_uuu, "nOS_uuu/I");
    globalTree->Branch("nZmass_uuu",&nZmass_uuu, "nZmass_uuu/I");
    globalTree->Branch("nJet_uuu",&nJet_uuu, "nJet_uuu/I");
    globalTree->Branch("nBJet_uuu",&nBJet_uuu, "nBJet_uuu/I");
    globalTree->Branch("nMWT_uuu",&nMWT_uuu, "nMWT_uuu/I");
    globalTree->Branch("nSMtop_uuu",&nSMtop_uuu, "nSMtop_uuu/I");
    globalTree->Branch("nMET_uuu",&nMET_uuu, "nMET_uuu/I");
    globalTree->Branch("nTrigg_uue",&nTrigg_uue, "nTrigg_uue/I");
    globalTree->Branch("n3lep_uue",&n3lep_uue, "n3lep_uue/I");
    globalTree->Branch("nVetoMu_uue",&nVetoMu_uue, "nVetoMu_uue/I");
    globalTree->Branch("nVetoEl_uue",&nVetoEl_uue, "nVetoEl_uue/I");
    globalTree->Branch("nOS_uue",&nOS_uue, "nOS_uue/I");
    globalTree->Branch("nZmass_uue",&nZmass_uue, "nZmass_uue/I");
    globalTree->Branch("nJet_uue",&nJet_uue, "nJet_uue/I");
    globalTree->Branch("nBJet_uue",&nBJet_uue, "nBJet_uue/I");
    globalTree->Branch("nMWT_uue",&nMWT_uue, "nMWT_uue/I");
    globalTree->Branch("nSMtop_uue", &nSMtop_uue, "nSMtop_uue/I");
    globalTree->Branch("nMET_uue",&nMET_uue, "nMET_uue/I");
    
    
    
    
    // event related variables
    baselineTree->Branch("signalInt", &signalInt, "signalInt/I");

    baselineTree->Branch("channelInt", &channelInt, "channelInt/I");
    
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
    
    
    baselineTree->Branch("LP", &LP, "LP/F");
    baselineTree->Branch("CosTheta",&CosTheta,"CosTheta/F");
    baselineTree->Branch("CosTheta_alt",&CosTheta_alt,"CosTheta_alt/F");
    baselineTree->Branch("Wlep_Charge",&Wlep_Charge,"Wlep_Charge/F");
    baselineTree->Branch("Wlep_Phi", &Wlep_Phi, "Wlep_Phi/F");
    baselineTree->Branch("Wlep_Eta", &Wlep_Eta, "Wlep_Eta/F");
    baselineTree->Branch("Wlep_Pt", &Wlep_Pt, "Wlep_Pt/F");
    baselineTree->Branch("Zboson_Phi", &Zboson_Phi, "Zboson_Phi/F");
    baselineTree->Branch("Zboson_Pt", &Zboson_Pt, "Zboson_Pt/F");
    baselineTree->Branch("Zboson_Eta", &Zboson_Eta, "Zboson_Eta/F");
    baselineTree->Branch("charge_asym", &charge_asym, "charge_asym/F");
    baselineTree->Branch("dPhiZWlep", &dPhiZWlep, "dPhiZWlep/F");
    baselineTree->Branch("dRZWlep", &dRZWlep, "dRZWlep/F");
    baselineTree->Branch("SMtop_Pt", &SMtop_Pt, "SMtop_Pt/F");
    baselineTree->Branch("SMtop_Phi", &SMtop_Phi, "SMtop_Phi/F");
    baselineTree->Branch("SMtop_Eta", &SMtop_Eta,"SMtop_Eta/F");
    baselineTree->Branch("dPhiZMET", &dPhiZMET, "dPhiZMET/F");
    baselineTree->Branch("dPhiZSMtop", &dPhiZSMtop, "dPhiZSMtop/F");
    baselineTree->Branch("dRZSMtop", &dRZSMtop,"dRZSMtop/F");
    baselineTree->Branch("TotalInvMass", &TotalInvMass, "TotalInvMass/F");
    baselineTree->Branch("TotalPt", &TotalPt, "TotalPt/F");
    baselineTree->Branch("TotalHt", &TotalHt, "TotalHt/F");
    baselineTree->Branch("cdiscCvsB_jet_2", &cdiscCvsB_jet_2, "cdiscCvsB_jet_2/F");
    baselineTree->Branch("cdiscCvsL_jet_2", &cdiscCvsL_jet_2, "cdiscCvsL_jet_2/F");
    baselineTree->Branch("bdiscCSVv2_jet_1", &bdiscCSVv2_jet_1, "bdiscCSVv2_jet_1/F");
    baselineTree->Branch("bdiscCSVv2_jet_2", &bdiscCSVv2_jet_2, "bdiscCSVv2_jet_2/F");
    baselineTree->Branch("PassedTrigger", &PassedTrigger, "PassedTrigger/I");
    
    
    
    baselineTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//
    baselineTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/F");
    baselineTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/F");
    baselineTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/F");
    baselineTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/F");
    baselineTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/F");
    baselineTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/F");
    baselineTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/F");
    baselineTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/F");
    baselineTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/F");
    baselineTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/F");
    baselineTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
    baselineTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/F");
    baselineTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/F");
    baselineTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/F");
    baselineTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/F");
    baselineTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/F");
    baselineTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/F");
    baselineTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
    baselineTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
    baselineTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
    baselineTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
    baselineTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
    baselineTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/F");
    baselineTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/F");
    baselineTree->Branch("pt_electron_3",&pt_electron_3,"pt_electron_3/F");
    
    
    baselineTree->Branch("nMuons",&nMuons, "nMuons/I");
    baselineTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    baselineTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    baselineTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    baselineTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
    baselineTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/F");
    baselineTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/F");
    baselineTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/F");
    baselineTree->Branch("E_muon",E_muon,"E_muon[nMuons]/F");
    baselineTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/F");
    baselineTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/F");
    baselineTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/F");
    baselineTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
    baselineTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
    baselineTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/F");
    baselineTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
    baselineTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/F");
    baselineTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/F");
    baselineTree->Branch("pt_muon_1",&pt_muon_1,"pt_muon_1/F");
    baselineTree->Branch("pt_muon_2",&pt_muon_2,"pt_muon_2/F");
    baselineTree->Branch("pt_muon_3",&pt_muon_3,"pt_muon_3/F");
    
       baselineTree->Branch("nJets",&nJets,"nJets/I");
    baselineTree->Branch("cdiscCvsL_jet_1",&cdiscCvsL_jet_1,"cdiscCvsL_jet_1/F");
    baselineTree->Branch("cdiscCvsB_jet_1",&cdiscCvsB_jet_1,"cdiscCvsB_jet_1/F");
    baselineTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
    baselineTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
    baselineTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
    if(Usettbar){
      baselineTree->Branch("nJets_CharmL",&nJets_CharmL,"nJets_CharmL/I");
      baselineTree->Branch("nJets_CharmM",&nJets_CharmM,"nJets_CharmM/I");
      baselineTree->Branch("nJets_CharmT",&nJets_CharmT,"nJets_CharmT/I");
      baselineTree->Branch("nJets_nonCSVL",&nJets_nonCSVL,"nJets_nonCSVL/I");
      baselineTree->Branch("nJets_nonCSVM",&nJets_nonCSVM,"nJets_nonCSVM/I");
      baselineTree->Branch("nJets_nonCSVT",&nJets_nonCSVT,"nJets_nonCSVT/I");
      baselineTree->Branch("nJets_nonCharmLCSVL",&nJets_nonCharmLCSVL,"nJets_nonCharmLCSVL/I");
      baselineTree->Branch("nJets_nonCharmLCSVM",&nJets_nonCharmLCSVM,"nJets_nonCharmLCSVM/I");
      baselineTree->Branch("nJets_nonCharmLCSVT",&nJets_nonCharmLCSVT,"nJets_nonCharmLCSVT/I");
      baselineTree->Branch("nJets_nonCharmMCSVL",&nJets_nonCharmMCSVL,"nJets_nonCharmMCSVL/I");
      baselineTree->Branch("nJets_nonCharmMCSVM",&nJets_nonCharmMCSVM,"nJets_nonCharmMCSVM/I");
      baselineTree->Branch("nJets_nonCharmMCSVT",&nJets_nonCharmMCSVT,"nJets_nonCharmMCSVT/I");
      baselineTree->Branch("nJets_nonCharmTCSVL",&nJets_nonCharmTCSVL,"nJets_nonCharmTCSVL/I");
      baselineTree->Branch("nJets_nonCharmTCSVM",&nJets_nonCharmTCSVM,"nJets_nonCharmTCSVM/I");
      baselineTree->Branch("nJets_nonCharmTCSVT",&nJets_nonCharmTCSVT,"nJets_nonCharmTCSVT/I");
      baselineTree->Branch("nJets_nonCharmLnonCSVL",&nJets_nonCharmLnonCSVL,"nJets_nonCharmLnonCSVL/I");
      baselineTree->Branch("nJets_nonCharmLnonCSVM",&nJets_nonCharmLnonCSVM,"nJets_nonCharmLnonCSVM/I");
      baselineTree->Branch("nJets_nonCharmLnonCSVT",&nJets_nonCharmLnonCSVT,"nJets_nonCharmLnonCSVT/I");
      baselineTree->Branch("nJets_nonCharmMnonCSVL",&nJets_nonCharmMnonCSVL,"nJets_nonCharmMnonCSVL/I");
      baselineTree->Branch("nJets_nonCharmMnonCSVM",&nJets_nonCharmMnonCSVM,"nJets_nonCharmMnonCSVM/I");
      baselineTree->Branch("nJets_nonCharmMnonCSVT",&nJets_nonCharmMnonCSVT,"nJets_nonCharmMnonCSVT/I");
      baselineTree->Branch("nJets_nonCharmTnonCSVL",&nJets_nonCharmTnonCSVL,"nJets_nonCharmTnonCSVL/I");
      baselineTree->Branch("nJets_nonCharmTnonCSVM",&nJets_nonCharmTnonCSVM,"nJets_nonCharmTnonCSVM/I");
      baselineTree->Branch("nJets_nonCharmTnonCSVT",&nJets_nonCharmTnonCSVT,"nJets_nonCharmTnonCSVT/I");
      
      baselineTree->Branch("nJets_CharmLnonCSVL",&nJets_CharmLnonCSVL,"nJets_CharmLnonCSVL/I");
      baselineTree->Branch("nJets_CharmLnonCSVM",&nJets_CharmLnonCSVM,"nJets_CharmLnonCSVM/I");
      baselineTree->Branch("nJets_CharmLnonCSVT",&nJets_CharmLnonCSVT,"nJets_CharmLnonCSVT/I");
      baselineTree->Branch("nJets_CharmMnonCSVL",&nJets_CharmMnonCSVL,"nJets_CharmMnonCSVL/I");
      baselineTree->Branch("nJets_CharmMnonCSVM",&nJets_CharmMnonCSVM,"nJets_CharmMnonCSVM/I");
      baselineTree->Branch("nJets_CharmMnonCSVT",&nJets_CharmMnonCSVT,"nJets_CharmMnonCSVT/I");
      baselineTree->Branch("nJets_CharmTnonCSVL",&nJets_CharmTnonCSVL,"nJets_CharmTnonCSVL/I");
      baselineTree->Branch("nJets_CharmTnonCSVM",&nJets_CharmTnonCSVM,"nJets_CharmTnonCSVM/I");
      baselineTree->Branch("nJets_CharmTnonCSVT",&nJets_CharmTnonCSVT,"nJets_CharmTnonCSVT/I");
      
      
      baselineTree->Branch("nJets_CharmL",&nJets_CharmL,"nJets_CharmL/I");
      baselineTree->Branch("nJets_CharmM",&nJets_CharmM,"nJets_CharmM/I");
      baselineTree->Branch("nJets_CharmT",&nJets_CharmT,"nJets_CharmT/I");
      baselineTree->Branch("nJets_nonCharmL",&nJets_nonCharmL,"nJets_nonCharmL/I");
      baselineTree->Branch("nJets_nonCharmM",&nJets_nonCharmM,"nJets_nonCharmM/I");
      baselineTree->Branch("nJets_nonCharmT",&nJets_nonCharmT,"nJets_nonCharmT/I");
      
      baselineTree->Branch("nJets_CharmLCSVL",&nJets_CharmLCSVL,"nJets_CharmLCSVL/I");
      baselineTree->Branch("nJets_CharmLCSVM",&nJets_CharmLCSVM,"nJets_CharmLCSVM/I");
      baselineTree->Branch("nJets_CharmLCSVT",&nJets_CharmLCSVT,"nJets_CharmLCSVT/I");
      baselineTree->Branch("nJets_CharmMCSVL",&nJets_CharmMCSVL,"nJets_CharmMCSVL/I");
      baselineTree->Branch("nJets_CharmMCSVM",&nJets_CharmMCSVM,"nJets_CharmMCSVM/I");
      baselineTree->Branch("nJets_CharmMCSVT",&nJets_CharmMCSVT,"nJets_CharmMCSVT/I");
      baselineTree->Branch("nJets_CharmTCSVL",&nJets_CharmTCSVL,"nJets_CharmTCSVL/I");
      baselineTree->Branch("nJets_CharmTCSVM",&nJets_CharmTCSVM,"nJets_CharmTCSVM/I");
      baselineTree->Branch("nJets_CharmTCSVT",&nJets_CharmTCSVT,"nJets_CharmTCSVT/I");
      
    }
    
    baselineTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/F");
    baselineTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/F");
    baselineTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/F");
    baselineTree->Branch("E_jet",E_jet,"E_jet[nJets]/F");
    baselineTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    baselineTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/F");
    if(Usettbar){
      baselineTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/F");
      baselineTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/F");
    }
    baselineTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/F");
    baselineTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/F");
    baselineTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/F");
    
    // Zboson
       baselineTree->Branch("Zboson_M",&Zboson_M,"Zboson_M/F");
  
    baselineTree->Branch("Zboson2_M",&Zboson2_M,"Zboson2_M/F");
    myTree->Branch("mWt",&mWt,"mWt/F");
    baselineTree->Branch("mWt",&mWt,"mWt/F");
    if(Usettbar){
    baselineTree->Branch("FCNCtop_tagger",&FCNCtop_M_tagger,"FCNCtop_M_tagger/F");
    }
       baselineTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/F");
    if(Usettbar){
         baselineTree->Branch("cjet_Pt",&cjet_Pt,"cjet_Pt/F");
      baselineTree->Branch("cjet_Pt_tagger",&cjet_Pt_tagger,"cjet_Pt_tagger/F");
    }
    myTree->Branch("mlb",&mlb,"mlb/F");
    baselineTree->Branch("mlb",&mlb,"mlb/F");
    myTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/F");
    myTree->Branch("dRZb",&dRZb,"dRZb/F");
    if(Usettbar){
      myTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/F");
      myTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/F");
      baselineTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/F");
      myTree->Branch("dRZc",&dRZc,"dRZc/F");
      myTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/F");
      baselineTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/F");
      myTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/F");
      myTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/F");
      baselineTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/F");
      myTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/F");
      baselineTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/F");
      myTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/F");
      myTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/F");
      myTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/F");
      baselineTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/F");
      baselineTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/F");
      baselineTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/F");
      baselineTree->Branch("dRZc",&dRZc,"dRZc/F");
      baselineTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/F");
      baselineTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/F");
      baselineTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/F");
      myTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/F");
      myTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/F");
      baselineTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/F");
    }
    myTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/F");
    myTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/F");
    baselineTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/F");
    baselineTree->Branch("dRZb",&dRZb,"dRZb/F");
    baselineTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/F");
    baselineTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/F");
    
    
    // met
    myTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
    myTree->Branch("met_Ptbf", &met_Ptbf, "met_Ptbf/F");
    myTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
    myTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
    myTree->Branch("met_Px", &met_Px, "met_Px/F");
    myTree->Branch("met_Py", &met_Py, "met_Py/F");
    myTree->Branch("met_Pz", &met_Pz, "met_Pz/F");
    
    myTree->Branch("orig_met_pt", &orig_met_pt, "orig_met_pt/F");
    myTree->Branch("orig_met_px", &orig_met_px, "orig_met_px/F");
    myTree->Branch("orig_met_py", &orig_met_py, "orig_met_py/F");
    myTree->Branch("corrected_met_pt", &corrected_met_pt, "corrected_met_pt/F");
    myTree->Branch("corrected_met_px", &corrected_met_px, "corrected_met_px/F");
    myTree->Branch("corrected_met_py", &corrected_met_py, "corrected_met_py/F");
    
    myTree->Branch("orig_jet_pt", &orig_jet_pt, "orig_jet_pt/F");
    myTree->Branch("orig_jet_px", &orig_jet_px, "orig_jet_px/F");
    myTree->Branch("orig_jet_py", &orig_jet_py, "orig_jet_py/F");
    myTree->Branch("corrected_jet_pt", &corrected_jet_pt, "corrected_jet_pt/F");
    myTree->Branch("corrected_jet_px", &corrected_jet_px, "corrected_jet_px/F");
    myTree->Branch("corrected_jet_py", &corrected_jet_py, "corrected_jet_py/F");
    myTree->Branch("jet_pt_check", &jet_pt_check, "jet_pt_check/F");
    
    
    baselineTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
    baselineTree->Branch("met_Ptbf", &met_Ptbf, "met_Ptbf/F");
    baselineTree->Branch("met_Px", &met_Px, "met_Px/F");
    baselineTree->Branch("met_Py", &met_Py, "met_Py/F");
    baselineTree->Branch("met_Pz", &met_Pz, "met_Pz/F");
    baselineTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
    baselineTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
    
    
    baselineTree->Branch("orig_met_pt", &orig_met_pt, "orig_met_pt/F");
    baselineTree->Branch("orig_met_px", &orig_met_px, "orig_met_px/F");
    baselineTree->Branch("orig_met_py", &orig_met_py, "orig_met_py/F");
    baselineTree->Branch("corrected_met_pt", &corrected_met_pt, "corrected_met_pt/F");
    baselineTree->Branch("corrected_met_px", &corrected_met_px, "corrected_met_px/F");
    baselineTree->Branch("corrected_met_py", &corrected_met_py, "corrected_met_py/F");
    
    baselineTree->Branch("orig_jet_pt", &orig_jet_pt, "orig_jet_pt/F");
    baselineTree->Branch("orig_jet_px", &orig_jet_px, "orig_jet_px/F");
    baselineTree->Branch("orig_jet_py", &orig_jet_py, "orig_jet_py/F");
    baselineTree->Branch("corrected_jet_pt", &corrected_jet_pt, "corrected_jet_pt/F");
    baselineTree->Branch("corrected_jet_px", &corrected_jet_px, "corrected_jet_px/F");
    baselineTree->Branch("corrected_jet_py", &corrected_jet_py, "corrected_jet_py/F");
    baselineTree->Branch("jet_pt_check", &jet_pt_check, "jet_pt_check/F");
    
    
    int ZmuIndiceF_0 = -999;
    int ZmuIndiceF_1 = -999;
    int ZelecIndiceF_0 = -999;
    int ZelecIndiceF_1= -999;
    int WmuIndiceF = -999;
    int WelecIndiceF = -999;
    baselineTree->Branch("ZmuIndiceF_0", &ZmuIndiceF_0, "ZmuIndiceF_0/I");
    baselineTree->Branch("ZelecIndiceF_0", &ZelecIndiceF_0,"ZelecIndiceF_0/I");
    baselineTree->Branch("ZmuIndiceF_1", &ZmuIndiceF_1, "ZmuIndiceF_1/I");
    baselineTree->Branch("ZelecIndiceF_1", &ZelecIndiceF_1,"ZelecIndiceF_1/I");
    baselineTree->Branch("WmuIndiceF", &WmuIndiceF, "WmuIndiceF/I");
    baselineTree->Branch("WelecIndiceF", &WelecIndiceF, "WelecIndiceF/I");
    
    
    if(verbose>1) cout << "trees created " << endl;
    
    /////////////////////////
    //// Corrections/trigger ///
    ///////////////////////////
    
    /// book triggers
    trigger_mumu->bookTriggers(isData,dName);
    trigger_ee->bookTriggers(isData,dName);
    trigger_emu->bookTriggers(isData,dName);
    trigger_mumumu->bookTriggers(isData, dName);
    trigger_eee->bookTriggers(isData,dName);
    trigger_emumu_mumue->bookTriggers(isData,dName);
    trigger_mu->bookTriggers(isData,dName);
    trigger_e->bookTriggers(isData, dName);
    
    
    
    if(verbose>1) cout << "triggers booked " << endl;
    
    //////////////////////////////////////////////////
    // Pre-event loop definitions
    /////////////////////////////////////////////////
    
    int itrigger = -1, previousRun = -1, start = 0;
    int currentRun;
    int iFile = -1;
    //cout << "before" << endl; if you get a error for the line below, check vecFileNames
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
    
    if (verbose > 0) cout << " - Loop over events " << endl;
    
    //define object containers
    
    // initial variables
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_jets_corrected;
    vector < TRootGenJet* >   genjets;
    vector < TRootMET* >      mets;
    vector<TRootElectron*>    selectedElectrons;
    vector<TRootElectron*>    selectedVetoElectrons;
    vector<TRootPFJet*>       selectedJets;
    vector<TRootPFJet*>       PreselectedJets;
    vector<TRootMuon*>        selectedMuons;
    vector<TRootMuon*>        selectedLooseMuons;
    vector<TRootPFJet*>      selectedCSVLBJets;
    vector<TRootPFJet*>      selectedCSVMBJets;
    vector<TRootPFJet*>      selectedCSVTBJets;
    
    
    vector<TRootMCParticle*> mcParticles;
    vector <TRootPFJet*>     selectednonCSVLJets;
    vector<TRootPFJet*>       selectedCharmLJets;
    vector<TRootPFJet*>       selectedCharmMJets;
    vector<TRootPFJet*>       selectedCharmTJets;
    vector <int>            selectedCharmLJetsindex;
    vector <int>            selectedCharmMJetsindex;
    vector <int>            selectedCharmTJetsindex;
    vector <TRootPFJet*>     selectednonCSVMJets;
    vector <TRootPFJet*>     selectednonCSVTJets;
    vector<TRootPFJet*>      selectednonCharmLJets;
    vector<TRootPFJet*>     selectednonCharmMJets;
    vector<TRootPFJet*>     selectednonCharmTJets;
    vector<TRootPFJet*>      selectednonCLBLJets;
    vector<TRootPFJet*>     selectednonCLBMJets;
    vector<TRootPFJet*>     selectednonCLBTJets;
    vector<TRootPFJet*>      selectednonCMBLJets;
    vector<TRootPFJet*>     selectednonCMBMJets;
    vector<TRootPFJet*>     selectednonCMBTJets;
    vector<TRootPFJet*>      selectednonCTBLJets;
    vector<TRootPFJet*>     selectednonCTBMJets;
    vector<TRootPFJet*>     selectednonCTBTJets;
    vector<TRootPFJet*>      selectednonCLnonBLJets;
    vector<TRootPFJet*>     selectednonCLnonBMJets;
    vector<TRootPFJet*>     selectednonCLnonBTJets;
    vector<TRootPFJet*>      selectednonCMnonBLJets;
    vector<TRootPFJet*>     selectednonCMnonBMJets;
    vector<TRootPFJet*>     selectednonCMnonBTJets;
    vector<TRootPFJet*>      selectednonCTnonBLJets;
    vector<TRootPFJet*>     selectednonCTnonBMJets;
    vector<TRootPFJet*>     selectednonCTnonBTJets;
    
    vector<TRootPFJet*>      selectedCLnonBLJets;
    vector<TRootPFJet*>     selectedCLnonBMJets;
    vector<TRootPFJet*>     selectedCLnonBTJets;
    vector<TRootPFJet*>      selectedCMnonBLJets;
    vector<TRootPFJet*>     selectedCMnonBMJets;
    vector<TRootPFJet*>     selectedCMnonBTJets;
    vector<TRootPFJet*>      selectedCTnonBLJets;
    vector<TRootPFJet*>     selectedCTnonBMJets;
    vector<TRootPFJet*>     selectedCTnonBTJets;
    vector<TRootPFJet*> selectedCLBLJets;
    vector<TRootPFJet*> selectedCLBMJets;
    vector<TRootPFJet*> selectedCLBTJets;
    vector<TRootPFJet*> selectedCMBLJets;
    vector<TRootPFJet*> selectedCMBMJets;
    vector<TRootPFJet*> selectedCMBTJets;
    vector<TRootPFJet*> selectedCTBLJets;
    vector<TRootPFJet*> selectedCTBMJets;
    vector<TRootPFJet*> selectedCTBTJets;
    
    
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
    pair< vector <TLorentzVector> , vector < pair < string , int > > >  AssignedLeptons;
    //////////////////////////////////////
    // Begin Event Loop
    //////////////////////////////////////
    nbEvents = 0;
    nofEventsHLTv2 = 0;
    nofEventsHLTv3 = 0;
    nofPosWeights = 0;
    nofNegWeights = 0;
    float eventweight = 1.;
    bool continueFlow ;
    nbSelectedEvents = 0;
    
    matchedWlep = 0;
    matchedZlep_0 = 0;
    matchedZlep_1 = 0;
    matchedEvents_Wlep = 0;
    matchedEvents_Zlep = 0;
    matchedEvents_Bjet = 0;
    matchedBjet = 0;
    matchedEvents_Cjet_Cmedium = 0;
    matchedEvents_Cjet_Ctight = 0;
    matchedEvents_Cjet_Cloose = 0;
    matchedCjet_Cloose = 0;
    matchedCjet_Cmedium = 0;
    matchedCjet_Ctight = 0;
    matchedEvents_Cjet = 0;
    matchedCjet = 0;
    matchedEvents_Cjet_CvsBtagger = 0;
    matchedCjet_CvsBtagger = 0;
    matchedEvents_Cjet_CvsLtagger = 0;
    matchedCjet_CvsLtagger = 0;
    int_eventForCjetmatchingmatched = 0;
    int_eventForCjetmatching = 0;
    int_eventForCjetmatchingmatched_CvsBtagger = 0;
    int_eventForCjetmatching_CvsBtagger = 0;
    int_eventForCjetmatchingmatched_CvsLtagger = 0;
    int_eventForCjetmatching_CvsLtagger = 0;
    int_eventForZlepmatchingmatched1 = 0;
    int_eventForWlepmatching = 0;
    int_eventForWlepmatchingmatched = 0;
    int_eventForZlepmatchingmatched0 = 0;
    int_eventForZlepmatching = 0;
    int_eventForBjetmatching = 0;
    int_eventForBjetmatchingmatched = 0;
    int_eventForCjetmatching_Cloose = 0;
    int_eventForCjetmatching_Cmedium = 0;
    int_eventForCjetmatching_Ctight = 0;
    int_eventForCjetmatchingmatched_Cloose = 0;
    int_eventForCjetmatchingmatched_Cmedium = 0;
    int_eventForCjetmatchingmatched_Ctight = 0;
    
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
    
    int nbEvents_eee_0 = 0;
    int nbEvents_eee_test = 0;
    int nbEvents_eee_1 = 0;
    int nbEvents_eee_1m = 0;
    int nbEvents_eee_2m = 0;
    int nbEvents_eee_2 = 0;
    int nbEvents_eee_3 = 0;
    int nbEvents_eee_4 = 0;
    int nbEvents_eee_5 = 0;
    int nbEvents_eee_6 = 0;
    int nbEvents_eee_7 = 0;
    int nbEvents_eee_8 = 0;
    int nbEvents_eee_9 = 0;
    
    int nbEvents_eeu_0 = 0;
    int nbEvents_eeu_test = 0;
    int nbEvents_eeu_1 = 0;
    int nbEvents_eeu_1m = 0;
    int nbEvents_eeu_2m = 0;
    int nbEvents_eeu_2 = 0;
    int nbEvents_eeu_3 = 0;
    int nbEvents_eeu_4 = 0;
    int nbEvents_eeu_5 = 0;
    int nbEvents_eeu_6 = 0;
    int nbEvents_eeu_7 = 0;
    int nbEvents_eeu_8 = 0;
    int nbEvents_eeu_9 = 0;
    
    int nbEvents_uuu_0 = 0;
    int nbEvents_uuu_test = 0;
    int nbEvents_uuu_1 = 0;
    int nbEvents_uuu_1m = 0;
    int nbEvents_uuu_2m = 0;
    int nbEvents_uuu_2 = 0;
    int nbEvents_uuu_3 = 0;
    int nbEvents_uuu_4 = 0;
    int nbEvents_uuu_5 = 0;
    int nbEvents_uuu_6 = 0;
    int nbEvents_uuu_7 = 0;
    int nbEvents_uuu_8 = 0;
    int nbEvents_uuu_9 = 0;
    
    int nbEvents_uue_0 = 0;
    int nbEvents_uue_test = 0;
    int nbEvents_uue_1 = 0;
    int nbEvents_uue_1m = 0;
    int nbEvents_uue_2m = 0;
    int nbEvents_uue_2 = 0;
    int nbEvents_uue_3 = 0;
    int nbEvents_uue_4 = 0;
    int nbEvents_uue_5 = 0;
    int nbEvents_uue_6 = 0;
    int nbEvents_uue_7 = 0;
    int nbEvents_uue_8 = 0;
    int nbEvents_uue_9 = 0;
    
    
    bool debug = false;
    vector <int> selections;
    std::ostringstream  selectionsnb;
    bool   passedMET = false;
    PassedGoodPV = false;
    bool   HBHEnoise = false;
    bool   HBHEIso = false;
    bool   CSCTight = false;
    bool badchan = false;
    bool badmu = false;
    bool   EcalDead = false;
    //bool    eeBad = false; not recommended
    bool   lep3 = false;
    bool lep2 = false;
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
    cutstep_string.clear();
    cutstep.clear();
    cutstep_eee.clear();
    cutstep_eeu.clear();
    cutstep_uue.clear();
    cutstep_uuu.clear();
    bool leptonsAssigned ;
    
    
    
    TLorentzVector tempObj;
    vector <TLorentzVector> selectedleptons_;
    
    vector <TLorentzVector> selectedobjects_;
    vector <TLorentzVector> selectedjets_;
    TLorentzVector totalOfObjects;
    Float_t pttotal_x;
    Float_t pttotal_y;
    Float_t httemp;
    
    for (unsigned int ievt = event_start; ievt < end_d; ievt++)
    {
      pttotal_x = 0.;
      pttotal_y = 0.;
      httemp = 0.;
      totalOfObjects.Clear();
      eventForCjetmatching_Ctight = false;
      eventForCjetmatching_Cmedium = false;
      eventForCjetmatching_Cloose = false;
      eventForCjetmatchingmatched_Cloose = false;
      eventForCjetmatchingmatched_Cmedium = false;
      eventForCjetmatchingmatched_Ctight = false;
      eventForCjetmatching = false;
      eventForCjetmatchingmatched = false;
      eventForCjetmatching_CvsBtagger = false;
      eventForCjetmatchingmatched_CvsBtagger = false;
      eventForCjetmatching_CvsLtagger = false;
      eventForCjetmatchingmatched_CvsLtagger = false;
      eventForBjetmatching = false;
      eventForBjetmatchingmatched = false;
      eventForWlepmatchingmatched = false;
      eventForZlepmatching = false;
      eventForZlepmatchingmatched0 = false;
      eventForZlepmatchingmatched1 = false;
      eventForWlepmatching = false;
      leptonsAssigned = false;
      //elecbool = false;
      // mubool = false;
      ZmuIndiceF_0 = -999;
      ZmuIndiceF_1 = -999;
      ZelecIndiceF_0 = -999;
      ZelecIndiceF_1= -999;
      WmuIndiceF= -999;
      WelecIndiceF = -999;
      eventSelected = false;
      baseSelected = false;
      continueFlow = true;
      lep3 = false;
      lep2 = false;
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
      mc_Zmass = -5.;
      reco_Zmass = -5.;
      /* for (Int_t i = 0; i < 200; i++)
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
       
       }*/
      nCuts = 0;
      
      passedMET = false;
      PassedGoodPV = false;
      HBHEnoise = false;
      HBHEIso = false;
      CSCTight = false;
      badchan = false;
      badmu = false;
      EcalDead = false;
      PassedTrigger = false;
      //eeBad = false;
      eventweight = 1.;
      if(verbose > 0 ) cout << "new event " << ievt << endl;
      double ievt_d = ievt;
      debug = false;
      if (verbose > 4 ) debug = true;
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
      
      if(verbose>1)
      {
        cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
        cout <<"Number of Muons Loaded: " << init_muons.size() <<endl;
        cout << "Number of Jets  Loaded: " << init_jets.size() << endl;
        cout << "Met px / py loaded: "<< mets[0]->Px() << " / " << mets[0]->Py() << endl;
      }
      
      
      nIniRecoLeptons=0;
      nIniRecoElectrons = 0;
      nIniRecoMuons = 0;
      for(int iMu = 0 ; iMu < init_muons.size(); iMu++){
        if(init_muons[iMu]->Pt() > 10.0) nIniRecoMuons++;
      }
      for(int iEl= 0 ; iEl < init_electrons.size(); iEl++){
        if(init_electrons[iEl]->Pt() > 10.0) nIniRecoElectrons++;
      }
      nIniRecoLeptons = nIniRecoMuons + nIniRecoElectrons;
      histo1D["nIniRecoLeptons"]->Fill(nIniRecoLeptons);
      histo1D["nIniRecoElectrons"]->Fill(nIniRecoElectrons);
      histo1D["nIniRecoMuons"]->Fill(nIniRecoMuons);
      
      
      //  take the event
      datasets[d]->eventTree()->LoadTree(ievt);
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      int currentRun = event->runId();
      run_num = event->runId();
      evt_num = event->eventId();
      HBHEnoise = event->getHBHENoiseFilter();
      HBHEIso = event->getHBHENoiseIsoFilter();
      CSCTight = event->getglobalTightHalo2016Filter();
      EcalDead = event->getEcalDeadCellTriggerPrimitiveFilter();
      //eeBad = event->getEEBadScFilter();
      badchan   = event-> getBadChCandFilter();
      badmu	    = event-> getBadPFMuonFilter();
      
      
      
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
      
      if(runHLT && dTitle.find("noTrig")==string::npos )  // FIXME old samples (not withHLT not reHLT) don't contain trigger info
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
        
        
        EM = (trigged_emumu_mumue|| trigged_emu);
        MM = (trigged_mumu || trigged_mumumu ) ;
        EE = (trigged_ee || trigged_eee );
        M  = ( trigged_mu );
        E  = (trigged_e);
        
        //for data
        if ( EM &&                               (emdataset) ) result_trigger = 1;
        if ( MM && !EM &&                        (mmdataset) ) result_trigger = 1;
        if ( EE && !EM && !MM &&                 (eedataset) ) result_trigger = 1;
        if ( M  && !EM && !MM && !EE &&          (mdataset ) ) result_trigger = 1;
        if ( E  && !EM && !MM && !EE && !M &&    (edataset ) ) result_trigger = 1;
        // for MC
        if ( EM &&                               !isData ) result_trigger = 1;
        if ( MM && !EM &&                        !isData ) result_trigger = 1;
        if ( EE && !EM && !MM &&                 !isData ) result_trigger = 1;
        if ( M  && !EM && !MM && !EE &&          !isData ) result_trigger = 1;
        if ( E  && !EM && !MM && !EE && !M &&    !isData ) result_trigger = 1;
        
        trigged = result_trigger;
        if(dName.find("NP")!=string::npos) trigged = true; // needs to be fixed with the new MC
        
        
      }
      else if(runHLT && dTitle.find("noTrig")!=string::npos ){
        trigged = true;
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
      
      if(verbose > 1) cout << "Apply trigger? " << runHLT << " trigged? " << trigged << endl;
      if(trigged) PassedTrigger = true;
      
      ////////////////////////////
      ///// JES - JER smearing     ////
      //////////////////////////
      JERon = 0;
      
      
      
      if(applyJER && !isData)
      {
       // cout << "applying JER" << endl;
        jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        JERon = 1;
      }
      JESon = 0;
      if(applyJES && !isData)
      {
       // cout << "applying JES" << endl;
        jetTools->correctJets(init_jets_corrected,event->fixedGridRhoFastjetAll() ,false);
        JESon = 1;
      }
      
      orig_jet_pt = 0.;
      orig_jet_px = 0.;
      orig_jet_py = 0.;
      for(int intJet = 0; intJet < init_jets.size(); intJet++)
      {
        orig_jet_pt = orig_jet_pt + init_jets[intJet]->Pt();
        orig_jet_px = orig_jet_px + init_jets[intJet]->Px();
        orig_jet_py = orig_jet_py + init_jets[intJet]->Py();
      }
      
      corrected_jet_pt = 0.;
      corrected_jet_px = 0.;
      corrected_jet_py = 0.;
      for(int intJet = 0; intJet < init_jets_corrected.size(); intJet++)
      {
        corrected_jet_pt = corrected_jet_pt + init_jets_corrected[intJet]->Pt();
        corrected_jet_px = corrected_jet_px + init_jets_corrected[intJet]->Px();
        corrected_jet_py = corrected_jet_py + init_jets_corrected[intJet]->Py();
      }
      for(int intJet = 0; intJet < init_jets_corrected.size(); intJet++)
      {
        jet_pt_check = init_jets_corrected[intJet]->Pt() - init_jets[intJet]->Pt();
      }
      
      /// propagate JEC to MET

      orig_met_px = mets[0]->Px();
      orig_met_py = mets[0]->Py();
      orig_met_pt = sqrt(orig_met_px*orig_met_px + orig_met_py*orig_met_py);
      if(applyJES ) // jer doesn't need to be applied ||  applyJER)) --> smeared type-1 corrected MET,  NOW only yes --> Type 1 corrected MET
      {
        jetTools->correctMETTypeOne(init_jets_corrected, mets[0], isData);
        METon = 1;
        //  if JES applied: replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which JEC is applied
        //  if JER applied:  replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which smearing is applied.
        // type 1 correction / sleard pmet correction
        
      }
      corrected_met_px = mets[0]->Px();
      corrected_met_py = mets[0]->Py();
      corrected_met_pt = sqrt(corrected_met_px*corrected_met_px + corrected_met_py*corrected_met_py);
      
      ///////////////////////////////////////////////////////////
      // Event selection
      ///////////////////////////////////////////////////////////
      
      // Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets,event->fixedGridRhoFastjetAll());
      PreselectedJets.clear();
      PreselectedJets  = selection.GetSelectedJets(jet_pt_cut,jet_eta_cut, true, "Loose");
      selectedMuons.clear();
      selectedLooseMuons.clear();
      selectedMuons = selection.GetSelectedMuons(25., 2.1, mu_iso_cut, "Tight", "Spring15");   // spring 15 still counts for 2016
      selectedLooseMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut_loose, "Loose", "Spring15"); // spring 15 still counts for 2016
      
      // pt, eta, iso // run normally
      selectedElectrons.clear();
      selectedVetoElectrons.clear();
      selectedElectrons = selection.GetSelectedElectrons(35., 2.1, "Tight","Spring16_80X",true,true);// pt, eta, WP point, campaign, cutbased, VID EA
      selectedVetoElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Veto","Spring16_80X",true,true);// pt, eta
      /// For MC Information
      mcParticles.clear();
      if(!isData) treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
      if(!isData) sort(mcParticles.begin(),mcParticles.end(),HighestPt());
      
      
      
      /*  bool tightID = false;
       bool vetoID =false;
       for(int iEl = 0; iEl < init_electrons.size() ; iEl++){
       
       if(init_electrons[iEl]->isCB_TightID()) tightID = true;
       if(init_electrons[iEl]->isCB_VetoID()) vetoID = true;
       
       cout << "|" << evt_num << "|" << init_electrons[iEl]->Pt() << "|" << vetoID << "|" << tightID << "|" <<  endl;
       
       }
       */
      
      
      
      /*
       
       3 parts: unc related to jets: JEC - JER
       unc related to unclustered E (see below)
       unc related to leptons --> can be accounted for in uncl energy
        SEE https://twiki.cern.ch/twiki/bin/view/CMS/TopJME#MET
       
       NOW different!!! see https://indico.cern.ch/event/510570/contributions/1190890/attachments/1246661/1836294/JME_METUnc_210316.pdf
       // Set up the unclustered MET systematic
       double uncmet_px = mets[0]->Px();
       double uncmet_py = mets[0]->Py();
       for(unsigned int i=0; i<init_jets_corrected.size(); i++){
       uncmet_px += init_jets[i]->Px();
       uncmet_py += init_jets[i]->Py();
       }
       for(unsigned int i=0; i<init_muons.size(); i++){
       uncmet_px += init_muons[i]->Px();
       uncmet_py += init_muons[i]->Py();
       }
       for(unsigned int i=0; i<init_electrons.size(); i++){
       uncmet_px += init_electrons[i]->Px();
       uncmet_py += init_electrons[i]->Py();
       }
       
       double met_px = mets[0]->Px();
       double met_py = mets[0]->Py();
       
       if(unclusteredUp){  --> now it should be each within their resolution
       met_px += uncmet_px*0.1;
       met_py += uncmet_py*0.1;
       } if(unclusteredDown){
       met_px -= uncmet_px*0.1;
       met_py -= uncmet_py*0.1;
       }
       
       double met_pt = sqrt(met_px*met_px + met_py*met_py);
       
       */
      
      
      // void TTreeLoader::LoadMCEvent(int, TopTree::TRootNPGenEvent*, std::vector<TopTree::TRootMCParticle*>&, bool)
      if (verbose>1) cout <<"Number of Muons, Electrons, Jets  ===>  " << endl << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< PreselectedJets.size()   << endl;
      

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
      
      nRecoLeptons=0;
      nRecoElectrons = 0;
      nRecoMuons = 0;
      for(int iMu = 0 ; iMu < selectedMuons.size(); iMu++){
        if(selectedMuons[iMu]->Pt() > 10.0) nRecoMuons++;
      }
      for(int iEl= 0 ; iEl < selectedElectrons.size(); iEl++){
        if(selectedElectrons[iEl]->Pt() > 10.0) nRecoElectrons++;
      }
      nRecoLeptons = nRecoMuons + nRecoElectrons;
      histo1D["nRecoLeptons"]->Fill(nRecoLeptons);
      histo1D["nRecoElectrons"]->Fill(nRecoElectrons);
      histo1D["nRecoMuons"]->Fill(nRecoMuons);
      histo2D["nIniRecoMuonsnRecoMuons"]->Fill(nIniRecoMuons, nRecoMuons);
      histo2D["nIniRecoElectronsnRecoElectrons"]->Fill(nIniRecoElectrons, nRecoElectrons);
      histo2D["nIniRecoLeptonsnRecoLeptons"]->Fill(nIniRecoLeptons, nRecoLeptons);
      
      
      
      tempObj.Clear();
      selectedleptons_.clear();
      selectedobjects_.clear();
      selectedjets_.clear();
      for(unsigned int iLep = 0 ; iLep < selectedElectrons.size(); iLep++)
      {
        tempObj.Clear();
        tempObj.SetPtEtaPhiE(selectedElectrons[iLep]->Pt(), selectedElectrons[iLep]->Eta(), selectedElectrons[iLep]->Phi(), selectedElectrons[iLep]->E());
        selectedleptons_.push_back(tempObj);
        selectedobjects_.push_back(tempObj);
        
      }
      for(unsigned int iLep = 0 ; iLep < selectedMuons.size(); iLep++)
      {
        tempObj.Clear();
        tempObj.SetPtEtaPhiE(selectedMuons[iLep]->Pt(), selectedMuons[iLep]->Eta(), selectedMuons[iLep]->Phi(), selectedMuons[iLep]->E());
        selectedleptons_.push_back(tempObj);
        selectedobjects_.push_back(tempObj);
        
      }
      for(unsigned int iLep = 0 ; iLep < selectedJets.size(); iLep++)
      {
        tempObj.Clear();
        tempObj.SetPtEtaPhiE(selectedJets[iLep]->Pt(), selectedJets[iLep]->Eta(), selectedJets[iLep]->Phi(), selectedJets[iLep]->E());
        selectedjets_.push_back(tempObj);
        selectedobjects_.push_back(tempObj);
        
      }
      
      bool foundAllObjects = true;
      bool foundAllJets = true;
      bool foundAllLeptons = true;
      if(matching && !istZq && Usettbar) {
        LeptonMatcherPair =  LeptonMatcher(mcParticles, evt_num, selectedleptons_);
        if((LeptonMatcherPair.second).size() < 3) foundAllLeptons = false; // only when all partons have found a match
        //cout << "matching (MatcherPair.second).size() " << (MatcherPair.second).size() << endl;
        //cout << "matching (MatcherPair.first).size() " << (MatcherPair.first).size() << endl;
        
        JetMatcherPair =  JetMatcher(mcParticles, evt_num, selectedjets_);
        if((JetMatcherPair.second).size() < 2) foundAllJets = false;
        
        ObjectMatcherPair =  ObjectMatcher(mcParticles, evt_num, selectedjets_);
        if((ObjectMatcherPair.second).size() < 2) foundAllObjects = false;
        
      }
      if(matching && !istZq && !Usettbar) {
        LeptonMatcherPair =  LeptonMatcherST(mcParticles, evt_num, selectedleptons_);
        if((LeptonMatcherPair.second).size() < 3) foundAllLeptons = false; // only when all partons have found a match
        //cout << "matching (MatcherPair.second).size() " << (MatcherPair.second).size() << endl;
        //cout << "matching (MatcherPair.first).size() " << (MatcherPair.first).size() << endl;
        
        JetMatcherPair =  JetMatcherST(mcParticles, evt_num, selectedjets_);
        if((JetMatcherPair.second).size() < 1) foundAllJets = false;
      }
      if(matching && istZq ) {
        LeptonMatcherPair =  LeptonMatchertZq(mcParticles, evt_num, selectedleptons_);
        if((LeptonMatcherPair.second).size() <3) foundAllLeptons = false; // only when all partons have found a match
        //cout << "matching (MatcherPair.second).size() " << (MatcherPair.second).size() << endl;
        //cout << "matching (MatcherPair.first).size() " << (MatcherPair.first).size() << endl;
        
        JetMatcherPair =  JetMatchertZq(mcParticles, evt_num, selectedjets_);
        if((JetMatcherPair.second).size() < 1) foundAllJets = false;
        
      }
      
      ////////////////////////////////////////////////
      // Pre cut operations
      ////////////////////////////////////////////////
      // Apply primary vertex selection
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      // Met filters
      //if(HBHEnoise && HBHEIso && CSCTight && EcalDead && eeBad && isGoodPV && badchan && badmu) passedMET = true;
      if(HBHEnoise && HBHEIso && CSCTight && EcalDead  && isGoodPV && badchan && badmu) passedMET = true;
      PassedMETFilter = passedMET;
      PassedGoodPV = isGoodPV;
      
      //////////////////////////////////////
      //   B jet selection	       ////
      ///////////////////////////////////////
      
      selectednonCLnonBLJets.clear();
      selectednonCLnonBMJets.clear();
      selectednonCLnonBTJets.clear();
      selectednonCMnonBLJets.clear();
      selectednonCMnonBMJets.clear();
      selectednonCMnonBTJets.clear();
      selectednonCTnonBLJets.clear();
      selectednonCTnonBMJets.clear();
      selectednonCTnonBTJets.clear();
      
      selectedCLnonBLJets.clear();
      selectedCLnonBMJets.clear();
      selectedCLnonBTJets.clear();
      selectedCMnonBLJets.clear();
      selectedCMnonBMJets.clear();
      selectedCMnonBTJets.clear();
      selectedCTnonBLJets.clear();
      selectedCTnonBMJets.clear();
      selectedCTnonBTJets.clear();
      
      
      
      
      selectedCSVLBJets.clear();
      selectedCSVMBJets.clear();
      selectedCSVTBJets.clear();
      selectednonCSVLJets.clear();
      selectedCharmLJets.clear();
      selectedCharmMJets.clear();
      selectedCharmTJets.clear();
      selectednonCSVMJets.clear();
      selectednonCSVTJets.clear();
      selectednonCharmLJets.clear();
      selectednonCharmMJets.clear();
      selectednonCharmTJets.clear();
      selectednonCLBLJets.clear();
      selectednonCLBMJets.clear();
      selectednonCLBTJets.clear();
      selectednonCMBLJets.clear();
      selectednonCMBMJets.clear();
      selectednonCMBTJets.clear();
      selectednonCTBLJets.clear();
      selectednonCTBMJets.clear();
      selectednonCTBTJets.clear();
      selectedCLBLJets.clear();
      selectedCLBMJets.clear();
      selectedCLBTJets.clear();
      selectedCMBLJets.clear();
      selectedCMBMJets.clear();
      selectedCMBTJets.clear();
      selectedCTBLJets.clear();
      selectedCTBMJets.clear();
      selectedCTBTJets.clear();
      selectedCharmLJetsindex.clear();
      selectedCharmMJetsindex.clear();
      selectedCharmTJetsindex.clear();
      for(unsigned int iJ = 0; iJ < selectedJets.size(); iJ++)
      {
        //bjets
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose) { selectedCSVLBJets.push_back(selectedJets[iJ]);    }
        else{          selectednonCSVLJets.push_back(selectedJets[iJ]);        }
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium){ selectedCSVMBJets.push_back(selectedJets[iJ]);}
        else { selectednonCSVMJets.push_back(selectedJets[iJ]);}
        if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight){ selectedCSVTBJets.push_back(selectedJets[iJ]);}
        else {selectednonCSVTJets.push_back(selectedJets[iJ]);}
        

          //cjets
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first){   selectedCharmLJets.push_back(selectedJets[iJ]);  selectedCharmLJetsindex.push_back(iJ); }
          else{   selectednonCharmLJets.push_back(selectedJets[iJ]);}
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first){   selectedCharmMJets.push_back(selectedJets[iJ]);  selectedCharmMJetsindex.push_back(iJ); }
          else{   selectednonCharmMJets.push_back(selectedJets[iJ]);    }
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first){   selectedCharmTJets.push_back(selectedJets[iJ]); selectedCharmTJetsindex.push_back(iJ);  }
          else{   selectednonCharmTJets.push_back(selectedJets[iJ]);    }
      }
      
                WPb_L =  workingpointvalue_Loose;
      WPb_M =  workingpointvalue_Medium;
      WPb_T =  workingpointvalue_Tight;
   
        WPc_CvsL_Loose = c_workingpointvalue_Loose.first;
        WPc_CvsB_Loose = c_workingpointvalue_Loose.second;
        WPc_CvsL_Medium = c_workingpointvalue_Medium.first;
        WPc_CvsB_Medium = c_workingpointvalue_Medium.second;
        WPc_CvsL_Tight = c_workingpointvalue_Tight.first;
        WPc_CvsB_Tight = c_workingpointvalue_Tight.second;
      
      
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
        btagWeight =  btwt->getMCEventWeight(selectedJets,false);  // use parton flavour = true or hadron flavour = false
        
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
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_0++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_0++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_0++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_0++; }
          
        }
      }
      else{
        selections.push_back(0);
        continueFlow = false;
      }
      
      
      if(dName.find("DoubleEG")!=string::npos && selectedElectrons.size() < 2) { continueFlow = false; }
      else if(dName.find("DoubleEG")!=string::npos) { nbEvents_test++ ;}
      if(dName.find("DoubleMu")!=string::npos && selectedMuons.size() < 2) { continueFlow = false; }
      else if(dName.find("DoubleMu")!=string::npos) { nbEvents_test++ ;}
      if(dName.find("MuonEG")!=string::npos && (selectedElectrons.size() < 1 || selectedMuons.size() < 1)) { continueFlow = false; }
      else if(dName.find("MuonEG")!=string::npos){ nbEvents_test++ ;}
      if((dName.find("SingleMu") != string::npos) && selectedMuons.size() < 1) {continueFlow= false; }
      else if(dName.find("SingleMu")!= string::npos){ nbEvents_test++;}
      if((dName.find("SingleEl") != string::npos) && selectedElectrons.size() < 1) {continueFlow= false; }
      else if(dName.find("SingleEl") != string::npos){ nbEvents_test++;}
      /*if(trigged && (selectedElectrons.size()!=0 || selectedMuons.size() !=0)){
       cout << "Name, Trigged, Flow" << dName << " " << trigged << " " << continueFlow << endl;
       }*/
      
      
      if(((selectedMuons.size() + selectedElectrons.size()) != 3)){
        selections.push_back(0);
        if((selectedMuons.size() + selectedElectrons.size()) >1){
          
          if(continueFlow){
            
            //baseSelected = true;
          }
          lep2 = true;
          if(selectedMuons.size() == 3) {channelInt = 0; i_channel = 0;}
          else if(selectedElectrons.size() == 3) {channelInt = 3; i_channel = 3;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {channelInt = 2; i_channel = 2; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){channelInt = 1; i_channel = 1; }
          
        }
        
        continueFlow = false;
      }
      else if((selectedMuons.size() + selectedElectrons.size()) == 3){
        selections.push_back(1);
        if(continueFlow){
          histo1D["cutFlow"]->Fill(1., eventweight);
          nCuts++;
          nbEvents_1++;
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_1++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_1++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_1++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_1++; }
          
        }
        lep3 = true;
        if(selectedMuons.size() == 3) {channelInt = 0; i_channel = 0;}
        else if(selectedElectrons.size() == 3) {channelInt = 3; i_channel = 3;}
        else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {channelInt = 2; i_channel = 2; }
        else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){channelInt = 1; i_channel = 1; }
        else cout << "ERROR no channel selected" << endl;
      }
      
      //cout << "LOOKING AT CHANNEL " << channelInt << endl;
      
      if(selectedMuons.size() == selectedLooseMuons.size() && continueFlow) {
        nbEvents_1m++;
        nCuts++;
        
        
        if(selectedMuons.size() == 3) {nbEvents_uuu_1m++;}
        else if(selectedElectrons.size() == 3) {nbEvents_eee_1m++;}
        else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_1m++; }
        else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_1m++; }
      }
      else continueFlow = false;
      if(selectedVetoElectrons.size() == selectedElectrons.size() && continueFlow){
        nbEvents_2m++;
        nCuts++;
        
        
        if(selectedMuons.size() == 3) {nbEvents_uuu_2m++;}
        else if(selectedElectrons.size() == 3) {nbEvents_eee_2m++;}
        else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_2m++; }
        else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_2m++; }
      }
      else continueFlow = false;
      
      if((selectedMuons.size() != selectedLooseMuons.size()) || (selectedVetoElectrons.size() != selectedElectrons.size())){
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
      if(lep2)
      {
        // cout << "in two leptons " << selectedElectrons.size() << " " << selectedMuons.size() << endl;
        // vector <TLorentzVector> leptons = LeptonAssignerv2(selectedElectrons, selectedMuons);
        /* if(Assigned){
         TLorentzVector lep0;
         TLorentzVector lep1;
         lep0.SetPxPyPzE(leptons[0].Px(), leptons[0].Py(), leptons[0].Pz(), leptons[0].Energy());
         lep1.SetPxPyPzE(leptons[1].Px(), leptons[1].Py(), leptons[1].Pz(), leptons[1].Energy());
         TLorentzVector Zboson2;
         Zboson2.SetPxPyPzE(( lep0 + lep1).Px() ,( lep0 + lep1).Py(),( lep0 + lep1).Pz(),( lep0 + lep1).Energy()) ;
         Zboson2_M = (lep0+lep1).M();
         Zboson2_Px = ( lep0 + lep1).Px();
         Zboson2_Py = ( lep0 + lep1).Py();
         Zboson2_Pz = ( lep0 + lep1).Pz();
         Zboson2_Energy = ( lep0 + lep1).Energy();
         }*/
        
      }
      
      if(lep3){
        //cout << "assigning leptons " << endl;
        
        AssignedLeptons = LeptonAssigner(selectedElectrons, selectedMuons);
        
        
        //vector < TLorentzVector > FoundLeptons;
        //FoundLeptons = AssignedLeptons.first;
        //vector< pair < string, int > > FoundLeptonsIndices;
        //FoundLeptonsIndices = AssignedLeptons.second;
       
        for(unsigned int iF = 0; iF < (AssignedLeptons.second).size(); iF++){
          if((((AssignedLeptons.second)[iF]).first).find("Wmu")!=string::npos){ WmuIndiceF = ((AssignedLeptons.second)[iF]).second;  }
          if((((AssignedLeptons.second)[iF]).first).find("Wel")!=string::npos){ WelecIndiceF = ((AssignedLeptons.second)[iF]).second; }
          if((((AssignedLeptons.second)[iF]).first).find("Zmu_0")!=string::npos) ZmuIndiceF_0 = ((AssignedLeptons.second)[iF]).second;
          if((((AssignedLeptons.second)[iF]).first).find("Zmu_1")!=string::npos) ZmuIndiceF_1 = ((AssignedLeptons.second)[iF]).second;
          if((((AssignedLeptons.second)[iF]).first).find("Zel_0")!=string::npos) ZelecIndiceF_0 = ((AssignedLeptons.second)[iF]).second;
          if((((AssignedLeptons.second)[iF]).first).find("Zel_1")!=string::npos) ZelecIndiceF_1 = ((AssignedLeptons.second)[iF]).second;
        }
        
        int WlepIndice = -999;
        int ZlepIndice_0 = -999;
        int ZlepIndice_1 = -999;
        
        if(WelecIndiceF != -999){ WlepIndice = WelecIndiceF; }
        if(WmuIndiceF != -999){ WlepIndice = WmuIndiceF; }
        if(ZmuIndiceF_0 != -999) ZlepIndice_0 = ZmuIndiceF_0;
        if(ZelecIndiceF_0 != -999) ZlepIndice_0 = ZelecIndiceF_0;
        if(ZmuIndiceF_1 != -999) ZlepIndice_1 = ZmuIndiceF_1;
        if(ZelecIndiceF_1 != -999) ZlepIndice_1 = ZelecIndiceF_1;
        if(WlepIndice != -999 && ZlepIndice_0 != -999 && ZlepIndice_1 != -999){ leptonsAssigned = true; }
         //cout << "evt " << evt_num << " assigned " << leptonsAssigned <<  " found all objects " << foundAllLeptons <<  endl;
         //cout << "WmuIndice " << WmuIndiceF << " WelecIndice "<< WelecIndiceF << " ZmuIndice_0 "<< ZmuIndiceF_0 << " ZmuIndice_1 "<< ZmuIndiceF_1 <<" ZelecIndice_0 "<< ZelecIndiceF_0 <<" ZelecIndice_1 "<< ZelecIndiceF_1 << endl;
        // cout << "WlepIndice " << WlepIndice << " ZlepIndice_0 "<< ZlepIndice_0 << " ZlepIndice_1 "<< ZlepIndice_1 << endl;
        
        
        if(leptonsAssigned){
          if(matching && foundAllLeptons && !istZq && Usettbar ){
            
            
            //NPair = LeptonMatcherPair.second
            //PPair = vector< pair<unsigned int, unsigned int>> (LeptonMatcherPair.first)
            
            // cout << "(LeptonMatcherPair.second).size() " << (LeptonMatcherPair.second).size() << endl;
            // cout << "(LeptonMatcherPair.first).size() " << (LeptonMatcherPair.first).size() << endl;
            int WmuIndiceM = -999;
            int WelecIndiceM = -999;
            int ZelecIndiceM_0 = -999;
            int ZelecIndiceM_1 = -999;
            int ZmuIndiceM_0 = -999;
            int ZmuIndiceM_1 = -999;
            
            
            
            for(unsigned int iPart = 0 ; iPart < (LeptonMatcherPair.second).size(); iPart++){
              if((LeptonMatcherPair.second)[iPart].find("SMmu")!=string::npos){ WmuIndiceM = (LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ; }
              if((LeptonMatcherPair.second)[iPart].find("SMel")!=string::npos){ WelecIndiceM = (LeptonMatcherPair.first)[iPart].second ; }
              if((LeptonMatcherPair.second)[iPart].find("FCNCmumin")!=string::npos){ ZmuIndiceM_0 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("FCNCelmin")!=string::npos){ ZelecIndiceM_0 = (LeptonMatcherPair.first)[iPart].second; }
              if((LeptonMatcherPair.second)[iPart].find("FCNCmuplus")!=string::npos){ ZmuIndiceM_1 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("FCNCelplus")!=string::npos){ZelecIndiceM_1 = (LeptonMatcherPair.first)[iPart].second; }
            }
            
            // cout << "WmuIndiceM " << WmuIndiceM << " WelecIndiceM "<< WelecIndiceM << " ZmuIndiceM_0 "<< ZmuIndiceM_0 << " ZmuIndiceM_1 "<< ZmuIndiceM_1 <<" ZelecIndiceM_0 "<< ZelecIndiceM_0 <<" ZelecIndiceM_1 "<< ZelecIndiceM_1 << endl;
            
            if( WmuIndiceM != -999 && WmuIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WmuIndiceM == WmuIndiceF){ matchedWlep++;    eventForWlepmatchingmatched = true; }
            }
            if( WelecIndiceM != -999 && WelecIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WelecIndiceM == WelecIndiceF){ matchedWlep++;  eventForWlepmatchingmatched = true;}
            }
            if(ZmuIndiceM_0 != -999 && ZmuIndiceM_1 != -999 && ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZmuIndiceM_0 == ZmuIndiceF_0) { matchedZlep_0++;  eventForZlepmatchingmatched0 = true;  }
              else if(ZmuIndiceM_0 == ZmuIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true; }
              if(ZmuIndiceM_1 == ZmuIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
              else if(ZmuIndiceM_1 == ZmuIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
            }
            if(ZelecIndiceM_0 != -999 && ZelecIndiceM_1 != -999 && ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZelecIndiceM_0 == ZelecIndiceF_0) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;    }
              else if(ZelecIndiceM_0 == ZelecIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;   }
              if(ZelecIndiceM_1 == ZelecIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
              else if(ZelecIndiceM_1 == ZelecIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
            }
            
            
          }
          else if(matching && foundAllLeptons && !istZq && !Usettbar ){
            
            
            //NPair = LeptonMatcherPair.second
            //PPair = vector< pair<unsigned int, unsigned int>> (LeptonMatcherPair.first)
            
            // cout << "(LeptonMatcherPair.second).size() " << (LeptonMatcherPair.second).size() << endl;
            // cout << "(LeptonMatcherPair.first).size() " << (LeptonMatcherPair.first).size() << endl;
            int WmuIndiceM = -999;
            int WelecIndiceM = -999;
            int ZelecIndiceM_0 = -999;
            int ZelecIndiceM_1 = -999;
            int ZmuIndiceM_0 = -999;
            int ZmuIndiceM_1 = -999;
            
            
            
            for(unsigned int iPart = 0 ; iPart < (LeptonMatcherPair.second).size(); iPart++){
              if((LeptonMatcherPair.second)[iPart].find("SMmu")!=string::npos){ WmuIndiceM = (LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ; }
              if((LeptonMatcherPair.second)[iPart].find("SMel")!=string::npos){ WelecIndiceM = (LeptonMatcherPair.first)[iPart].second ; }
              if((LeptonMatcherPair.second)[iPart].find("FCNCmumin")!=string::npos){ ZmuIndiceM_0 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("FCNCelmin")!=string::npos){ ZelecIndiceM_0 = (LeptonMatcherPair.first)[iPart].second; }
              if((LeptonMatcherPair.second)[iPart].find("FCNCmuplus")!=string::npos){ ZmuIndiceM_1 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("FCNCelplus")!=string::npos){ZelecIndiceM_1 = (LeptonMatcherPair.first)[iPart].second; }
            }
            
            //cout << "WmuIndiceM " << WmuIndiceM << " WelecIndiceM "<< WelecIndiceM << " ZmuIndiceM_0 "<< ZmuIndiceM_0 << " ZmuIndiceM_1 "<< ZmuIndiceM_1 <<" ZelecIndiceM_0 "<< ZelecIndiceM_0 <<" ZelecIndiceM_1 "<< ZelecIndiceM_1 << endl;
            
            
            if( WmuIndiceM != -999 && WmuIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WmuIndiceM == WmuIndiceF){ matchedWlep++;    eventForWlepmatchingmatched = true; }
            }
            if( WelecIndiceM != -999 && WelecIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WelecIndiceM == WelecIndiceF){ matchedWlep++;  eventForWlepmatchingmatched = true;}
            }
            if(ZmuIndiceM_0 != -999 && ZmuIndiceM_1 != -999 && ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZmuIndiceM_0 == ZmuIndiceF_0) { matchedZlep_0++;  eventForZlepmatchingmatched0 = true;  }
              else if(ZmuIndiceM_0 == ZmuIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true; }
              if(ZmuIndiceM_1 == ZmuIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
              else if(ZmuIndiceM_1 == ZmuIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
            }
            if(ZelecIndiceM_0 != -999 && ZelecIndiceM_1 != -999 && ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZelecIndiceM_0 == ZelecIndiceF_0) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;    }
              else if(ZelecIndiceM_0 == ZelecIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;   }
              if(ZelecIndiceM_1 == ZelecIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
              else if(ZelecIndiceM_1 == ZelecIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
            }

            
            
          } // ST matching
          else if(matching && foundAllLeptons && istZq  ){
            
            
            //NPair = LeptonMatcherPair.second
            //PPair = vector< pair<unsigned int, unsigned int>> (LeptonMatcherPair.first)
            
            // cout << "(LeptonMatcherPair.second).size() " << (LeptonMatcherPair.second).size() << endl;
            // cout << "(LeptonMatcherPair.first).size() " << (LeptonMatcherPair.first).size() << endl;
            int WmuIndiceM = -999;
            int WelecIndiceM = -999;
            int ZelecIndiceM_0 = -999;
            int ZelecIndiceM_1 = -999;
            int ZmuIndiceM_0 = -999;
            int ZmuIndiceM_1 = -999;
            
            
            
            for(unsigned int iPart = 0 ; iPart < (LeptonMatcherPair.second).size(); iPart++){
              if((LeptonMatcherPair.second)[iPart].find("SMmu")!=string::npos){ WmuIndiceM = (LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ; }
              if((LeptonMatcherPair.second)[iPart].find("SMel")!=string::npos){ WelecIndiceM = (LeptonMatcherPair.first)[iPart].second ; }
              if((LeptonMatcherPair.second)[iPart].find("Radmumin")!=string::npos){ ZmuIndiceM_0 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("Radelmin")!=string::npos){ ZelecIndiceM_0 = (LeptonMatcherPair.first)[iPart].second; }
              if((LeptonMatcherPair.second)[iPart].find("Radmuplus")!=string::npos){ ZmuIndiceM_1 = ((LeptonMatcherPair.first)[iPart].second -  selectedElectrons.size() ); }
              if((LeptonMatcherPair.second)[iPart].find("Radelplus")!=string::npos){ZelecIndiceM_1 = (LeptonMatcherPair.first)[iPart].second; }
            }
            
           // cout << "WmuIndiceM " << WmuIndiceM << " WelecIndiceM "<< WelecIndiceM << " ZmuIndiceM_0 "<< ZmuIndiceM_0 << " ZmuIndiceM_1 "<< ZmuIndiceM_1 <<" ZelecIndiceM_0 "<< ZelecIndiceM_0 <<" ZelecIndiceM_1 "<< ZelecIndiceM_1 << endl;
            
            //cout << "WmuIndiceF " << WmuIndiceF << " WelecIndiceF "<< WelecIndiceF << " ZmuIndiceF_0 "<< ZmuIndiceF_0 << " ZmuIndiceF_1 "<< ZmuIndiceF_1 <<" ZelecIndiceF_0 "<< ZelecIndiceF_0 <<" ZelecIndiceF_1 "<< ZelecIndiceF_1 << endl;
            
            
            if( WmuIndiceM != -999 && WmuIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WmuIndiceM == WmuIndiceF){ matchedWlep++;    eventForWlepmatchingmatched = true; }
            }
            if( WelecIndiceM != -999 && WelecIndiceF != -999){
              matchedEvents_Wlep++;
              eventForWlepmatching = true;
              if(WelecIndiceM == WelecIndiceF){ matchedWlep++;  eventForWlepmatchingmatched = true;}
            }
            if(ZmuIndiceM_0 != -999 && ZmuIndiceM_1 != -999 && ZmuIndiceF_0 != -999 && ZmuIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZmuIndiceM_0 == ZmuIndiceF_0) { matchedZlep_0++;  eventForZlepmatchingmatched0 = true;  }
              else if(ZmuIndiceM_0 == ZmuIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true; }
              if(ZmuIndiceM_1 == ZmuIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
              else if(ZmuIndiceM_1 == ZmuIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;  }
            }
            if(ZelecIndiceM_0 != -999 && ZelecIndiceM_1 != -999 && ZelecIndiceF_0 != -999 && ZelecIndiceF_1 != -999){
              matchedEvents_Zlep++;
              eventForZlepmatching = true;
              if(ZelecIndiceM_0 == ZelecIndiceF_0) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;    }
              else if(ZelecIndiceM_0 == ZelecIndiceF_1) { matchedZlep_0++; eventForZlepmatchingmatched0 = true;   }
              if(ZelecIndiceM_1 == ZelecIndiceF_0) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
              else if(ZelecIndiceM_1 == ZelecIndiceF_1) { matchedZlep_1++; eventForZlepmatchingmatched1 = true;   }
            }
            
           // cout <<" matchedEvents_Wlep "<< matchedEvents_Wlep <<" matchedWlep "<< matchedWlep <<" matchedEvents_Zlep "<< matchedEvents_Zlep << " matchedZlep_0 "<< matchedZlep_0 <<" matchedZlep_1 "<< matchedZlep_1 << endl;
            
          } // tZq matching
          

          if(ZmuIndiceF_0 != -999) Zlep0.SetPxPyPzE(selectedMuons[ZmuIndiceF_0]->Px(), selectedMuons[ZmuIndiceF_0]->Py(), selectedMuons[ZmuIndiceF_0]->Pz(), selectedMuons[ZmuIndiceF_0]->Energy());
          if(ZmuIndiceF_1 != -999) Zlep1.SetPxPyPzE(selectedMuons[ZmuIndiceF_1]->Px(), selectedMuons[ZmuIndiceF_1]->Py(), selectedMuons[ZmuIndiceF_1]->Pz(), selectedMuons[ZmuIndiceF_1]->Energy());
          if(WmuIndiceF != -999 ) Wlep.SetPxPyPzE(selectedMuons[WmuIndiceF]->Px(), selectedMuons[WmuIndiceF]->Py(), selectedMuons[WmuIndiceF]->Pz(),selectedMuons[WmuIndiceF]->Energy());
          
          if(ZelecIndiceF_0 != -999) Zlep0.SetPxPyPzE(selectedElectrons[ZelecIndiceF_0]->Px(), selectedElectrons[ZelecIndiceF_0]->Py(), selectedElectrons[ZelecIndiceF_0]->Pz(), selectedElectrons[ZelecIndiceF_0]->Energy());
          if(ZelecIndiceF_1 != -999) Zlep1.SetPxPyPzE(selectedElectrons[ZelecIndiceF_1]->Px(), selectedElectrons[ZelecIndiceF_1]->Py(), selectedElectrons[ZelecIndiceF_1]->Pz(), selectedElectrons[ZelecIndiceF_1]->Energy());
          if(WelecIndiceF != -999 ) Wlep.SetPxPyPzE(selectedElectrons[WelecIndiceF]->Px(), selectedElectrons[WelecIndiceF]->Py(), selectedElectrons[WelecIndiceF]->Pz(),selectedElectrons[WelecIndiceF]->Energy());
          
          //double phis = Wlep.Phi() - mets[0]->Phi();
          //double cosphis = TMath::Cos(phis);
          mWt = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_px)*(Wlep.Px() + met_px) - (Wlep.Py() + met_py)* (Wlep.Py() + met_py));
          //mWtsecond = TMath::Sqrt(2*Wlep.Pt() * met_Pt*(1-cosphis));
          
          
          nCuts++;
          nbEvents_2++;
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_2++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_2++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_2++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_2++; }
          
          Zboson.Clear();
          
          Zboson.SetPxPyPzE(( Zlep0 + Zlep1).Px() ,( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Pz(),( Zlep0 + Zlep1).Energy()) ;
          Zboson_M = (Zlep0+Zlep1).M();
          Zboson_Px = ( Zlep0 + Zlep1).Px();
          Zboson_Py = ( Zlep0 + Zlep1).Py();
          Zboson_Pz = ( Zlep0 + Zlep1).Pz();
          Zboson_Energy = ( Zlep0 + Zlep1).Energy();
          Wlep_Pt = Wlep.Pt();
          Wlep_Eta = Wlep.Eta();
          Wlep_Phi = Wlep.Phi();
          if(WelecIndiceF != -999){ if(selectedElectrons[WelecIndiceF]->charge() > 0){ Wlep_Charge = 1.;}else{Wlep_Charge = -1.; }}
          if(WmuIndiceF != -999){ if(selectedMuons[WmuIndiceF]->charge() > 0){ Wlep_Charge = 1.;}else{Wlep_Charge = -1; }}
          //cout << "Wlep_charge " << Wlep_Charge << endl;
          Zboson_Pt = Zboson.Pt();
          Zboson_Eta = Zboson.Eta();
          Zboson_Phi = Zboson.Phi();
          charge_asym =Wlep_Charge*fabs(Wlep.Eta());
          
          dPhiZWlep = Zboson.DeltaPhi(Wlep);
          dRZWlep = Zboson.DeltaR(Wlep);
          
          
          histo1D["recoZmass"]->Fill((Zlep0+Zlep1).M());
          histo1D["recomWt"]->Fill(mWt);
          
          
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
          
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_3++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_3++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_3++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_3++; }
          histo1D["cutFlow"]->Fill(3., eventweight);
          // baseSelected = true;
        }
      }
      
      if( selectedJets.size() >0){
        if(continueFlow){ baseSelected = true};
      }
      
      
      
      //////////////////////////////////////
      //  DO STUFF WITH SELECTED EVENTS ////
      //////////////////////////////////////
      // fill the tree
      if(baseSelected){
        
        nJets = 0;
        TLorentzVector tempObject;
        for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
        {
          tempObject.Clear();
          tempObject.SetPxPyPzE(selectedJets[seljet]->Px(),selectedJets[seljet]->Py(),selectedJets[seljet]->Py(), selectedJets[seljet]->Energy());
          totalOfObjects =  totalOfObjects+tempObject;
          pttotal_x = pttotal_x + selectedJets[seljet]->Px();
          pttotal_y = pttotal_y + selectedJets[seljet]->Py();
          httemp = httemp +selectedJets[seljet]->Pt();
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
          
          if(selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose) {
            btagL_jet[nJets] = 1; }
          else{          btagL_jet[nJets] = 0;      }
          if(selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium) {
            btagM_jet[nJets] = 1; }
          else{          btagM_jet[nJets] = 0;      }
          if(selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Ti) {
            btagL_jet[nJets] = 1; }
          else{          btagL_jet[nJets] = 0;      }
          if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium){ selectedCSVMBJets.push_back(selectedJets[iJ]);}
          else { selectednonCSVMJets.push_back(selectedJets[iJ]);}
          if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight){ selectedCSVTBJets.push_back(selectedJets[iJ]);}
          else {selectednonCSVTJets.push_back(selectedJets[iJ]);}
          
          
          //cjets
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first){   selectedCharmLJets.push_back(selectedJets[iJ]);  selectedCharmLJetsindex.push_back(iJ); }
          else{   selectednonCharmLJets.push_back(selectedJets[iJ]);}
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first){   selectedCharmMJets.push_back(selectedJets[iJ]);  selectedCharmMJetsindex.push_back(iJ); }
          else{   selectednonCharmMJets.push_back(selectedJets[iJ]);    }
          if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first){   selectedCharmTJets.push_back(selectedJets[iJ]); selectedCharmTJetsindex.push_back(iJ);  }
          else{   selectednonCharmTJets.push_back(selectedJets[iJ]);    }
          
          nJets++;
        }
        if(selectedJets.size()>0 && Usettbar) cdiscCvsB_jet_1 = selectedJets[0]->ctag_pfCombinedCvsBJetTags();
        if(selectedJets.size()>0 && Usettbar) cdiscCvsL_jet_1 = selectedJets[0]->ctag_pfCombinedCvsLJetTags();
        if(selectedJets.size()>1 && Usettbar) cdiscCvsB_jet_2 = selectedJets[1]->ctag_pfCombinedCvsBJetTags();
        if(selectedJets.size()>1 && Usettbar) cdiscCvsL_jet_2 = selectedJets[1]->ctag_pfCombinedCvsLJetTags();
        if(selectedJets.size()>0 ) bdiscCSVv2_jet_1 = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
        if(selectedJets.size()>1 ) bdiscCSVv2_jet_2 = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
        if(selectedJets.size()>0) pt_jet_1 = selectedJets[0]->Pt();
        if(selectedJets.size()>1) pt_jet_2 = selectedJets[1]->Pt();
        if(selectedJets.size()>2) pt_jet_3 = selectedJets[2]->Pt();
        // bjets
        nJets_CSVT =  selectedCSVTBJets.size();
        nJets_CSVM =  selectedCSVMBJets.size();
        nJets_CSVL =  selectedCSVLBJets.size();
        nJets_nonCSVL = selectednonCSVLJets.size();
        nJets_nonCSVM = selectednonCSVMJets.size();
        nJets_nonCSVT = selectednonCSVTJets.size();
        if(Usettbar){
          // charm jets
          nJets_CharmL = selectedCharmLJets.size();
          nJets_CharmM = selectedCharmMJets.size();
          nJets_CharmT = selectedCharmTJets.size();
          nJets_nonCharmL = selectednonCharmLJets.size();
          nJets_nonCharmM = selectednonCharmMJets.size();
          nJets_nonCharmT = selectednonCharmTJets.size();
          
          
          // charm b jets
          nJets_CharmLCSVL = selectedCLBLJets.size();
          nJets_CharmLCSVM = selectedCLBMJets.size();
          nJets_CharmLCSVT = selectedCLBTJets.size();
          nJets_CharmMCSVL = selectedCMBLJets.size();
          nJets_CharmMCSVM = selectedCMBMJets.size();
          nJets_CharmMCSVT = selectedCMBTJets.size();
          nJets_CharmTCSVL = selectedCTBLJets.size();
          nJets_CharmTCSVM = selectedCTBMJets.size();
          nJets_CharmTCSVT = selectedCTBTJets.size();
          
          // non charm b jets
          nJets_nonCharmLCSVL = selectednonCLBLJets.size();
          nJets_nonCharmLCSVM = selectednonCLBMJets.size();
          nJets_nonCharmLCSVT = selectednonCLBTJets.size();
          nJets_nonCharmMCSVL = selectednonCMBLJets.size();
          nJets_nonCharmMCSVM = selectednonCMBMJets.size();
          nJets_nonCharmMCSVT = selectednonCMBTJets.size();
          nJets_nonCharmTCSVL = selectednonCTBLJets.size();
          nJets_nonCharmTCSVM = selectednonCTBMJets.size();
          nJets_nonCharmTCSVT = selectednonCTBTJets.size();
          
          //  charm non b jets
          nJets_CharmLnonCSVL = selectedCLnonBLJets.size();
          nJets_CharmLnonCSVM = selectedCLnonBMJets.size();
          nJets_CharmLnonCSVT = selectedCLnonBTJets.size();
          nJets_CharmMnonCSVL = selectedCMnonBLJets.size();
          nJets_CharmMnonCSVM = selectedCMnonBMJets.size();
          nJets_CharmMnonCSVT = selectedCMnonBTJets.size();
          nJets_CharmTnonCSVL = selectedCTnonBLJets.size();
          nJets_CharmTnonCSVM = selectedCTnonBMJets.size();
          nJets_CharmTnonCSVT = selectedCTnonBTJets.size();
          
          // non charm non b jets
          nJets_nonCharmLnonCSVL = selectednonCLnonBLJets.size();
          nJets_nonCharmLnonCSVM = selectednonCLnonBMJets.size();
          nJets_nonCharmLnonCSVT = selectednonCLnonBTJets.size();
          nJets_nonCharmMnonCSVL = selectednonCMnonBLJets.size();
          nJets_nonCharmMnonCSVM = selectednonCMnonBMJets.size();
          nJets_nonCharmMnonCSVT = selectednonCMnonBTJets.size();
          nJets_nonCharmTnonCSVL = selectednonCTnonBLJets.size();
          nJets_nonCharmTnonCSVM = selectednonCTnonBMJets.size();
          nJets_nonCharmTnonCSVT = selectednonCTnonBTJets.size();
          
        }
        
        
        
        nMuons = 0;
        double muonSFtemp = 1.;
        for (Int_t selmu =0; selmu < selectedMuons.size() ; selmu++ )
        {
          tempObject.Clear();
          tempObject.SetPxPyPzE(selectedMuons[selmu]->Px(),selectedMuons[selmu]->Py(),selectedMuons[selmu]->Py(), selectedMuons[selmu]->Energy());
          totalOfObjects =  totalOfObjects+tempObject;
          pttotal_x = pttotal_x + selectedMuons[selmu]->Px();
          pttotal_y = pttotal_y + selectedMuons[selmu]->Py();
          httemp = httemp +selectedMuons[selmu]->Pt();
          pt_muon[nMuons]=selectedMuons[selmu]->Pt();
          phi_muon[nMuons]=selectedMuons[selmu]->Phi();
          eta_muon[nMuons]=selectedMuons[selmu]->Eta();
          E_muon[nMuons]=selectedMuons[selmu]->E();
          
          pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
          if(!isData)
          {
            
            
            MuonIDSF[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            MuonIsoSF[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            
            
            
          }
          else
          {
            MuonIDSF[nMuons] = 1.;
            MuonIsoSF[nMuons] = 1.;
          }
          if(MuonIDSF[nMuons]*MuonIsoSF[nMuons] == 0 ) cout << "  MuonIDSF[nMuons] " <<  MuonIDSF[nMuons] << " MuonIsoSF[nMuons] " << MuonIsoSF[nMuons] << "  MuonIDSF[nMuons]*MuonIsoSF[nMuons] " <<    MuonIDSF[nMuons]*MuonIsoSF[nMuons]     << endl;
          if(muonSFtemp == 0) cout << " muon SF " << muonSFtemp * MuonIDSF[nMuons]*MuonIsoSF[nMuons] << endl;
          charge_muon[nMuons]=selectedMuons[selmu]->charge();
          nMuons++;
        }
        
        
        if(selectedMuons.size()>0) pt_muon_1 = selectedMuons[0]->Pt();
        if(selectedMuons.size()>1) pt_muon_2 = selectedMuons[1]->Pt();
        if(selectedMuons.size()>2) pt_muon_3 = selectedMuons[2]->Pt();
        nElectrons=0;
        for (Int_t selel =0; selel < selectedElectrons.size() ; selel++ )
        {
          tempObject.Clear();
          tempObject.SetPxPyPzE(selectedElectrons[selel]->Px(),selectedElectrons[selel]->Py(),selectedElectrons[selel]->Py(), selectedElectrons[selel]->Energy());
          totalOfObjects =  totalOfObjects+tempObject;
          pttotal_x = pttotal_x + selectedElectrons[selel]->Px();
          pttotal_y = pttotal_y + selectedElectrons[selel]->Py();
          httemp = httemp +selectedElectrons[selel]->Pt();
          pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
          phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
          eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
          eta_superCluster_electron[nElectrons]=selectedElectrons[selel]->superClusterEta();
          E_electron[nElectrons]=selectedElectrons[selel]->E();
          pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
          charge_electron[nElectrons]=selectedElectrons[selel]->charge();
          if(!isData){
            ElectronSF[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
            
          }
          else ElectronSF[nElectrons] = 1.;
          
          nElectrons++;
        }
        if(selectedElectrons.size()>0) pt_electron_1 = selectedElectrons[0]->Pt();
        if(selectedElectrons.size()>1) pt_electron_2 = selectedElectrons[1]->Pt();
        if(selectedElectrons.size()>2) pt_electron_3 = selectedElectrons[2]->Pt();
        
        TotalInvMass = totalOfObjects.M();
        TotalPt = sqrt(pttotal_x*pttotal_x+pttotal_y*pttotal_y); //totalOfObjects.Pt();
        TotalHt = httemp;
        nLeptons = nMuons + nElectrons;
        
      }
      
      
      
    if(baseSelected){ baselineTree->Fill(); }
      //if(selections.size() != 8) cout << "ERROR SOMETHING WENT WRONG WITH THE SELECTIONS " << endl;
      for(int inb = 0; inb <selections.size(); inb++)
      {
        selectionsnb << selections[inb];
      }
      
      
      
      
      
      
    } // end eventloop
    
    
    
    cutstep_string.push_back("trigger");
    cutstep_string.push_back("3lep");
    cutstep_string.push_back("VetoMu");
    cutstep_string.push_back("VetoEl");
    cutstep_string.push_back("OSSF");
    cutstep_string.push_back("Zmass");
    cutstep_string.push_back(">1jet");
    cutstep_string.push_back(">0CSVL");
    cutstep_string.push_back("mWt");
    cutstep_string.push_back("SMtop");
    cutstep_string.push_back("METfilter");
    
    
    for(int iC = 0; iC < cutstep_string.size() ; iC++){
      cout << "cut " << iC << " label " << cutstep_string[iC] << endl;
    }
    
    cutstep.push_back(nbEvents_0);
    cutstep.push_back(nbEvents_1);
    cutstep.push_back(nbEvents_2m);
    cutstep.push_back(nbEvents_3);
    cutstep.push_back(nbEvents_4);
    cutstep.push_back(nbEvents_5);
    cutstep.push_back(nbEvents_6);
    cutstep.push_back(nbEvents_7);
    cutstep.push_back(nbEvents_8);
    cutstep.push_back(nbEvents_9);
    
    cutstep_eee.push_back(nbEvents_eee_0);
    cutstep_eee.push_back(nbEvents_eee_1);
    cutstep_eee.push_back(nbEvents_eee_2m);
    cutstep_eee.push_back(nbEvents_eee_3);
    cutstep_eee.push_back(nbEvents_eee_4);
    cutstep_eee.push_back(nbEvents_eee_5);
    cutstep_eee.push_back(nbEvents_eee_6);
    cutstep_eee.push_back(nbEvents_eee_7);
    cutstep_eee.push_back(nbEvents_eee_8);
    cutstep_eee.push_back(nbEvents_eee_9);
    
    cutstep_eeu.push_back( nbEvents_eeu_0);
    cutstep_eeu.push_back(nbEvents_eeu_1);
    cutstep_eeu.push_back(nbEvents_eeu_2m);
    cutstep_eeu.push_back(nbEvents_eeu_3);
    cutstep_eeu.push_back(nbEvents_eeu_4);
    cutstep_eeu.push_back(nbEvents_eeu_5);
    cutstep_eeu.push_back(nbEvents_eeu_6);
    cutstep_eeu.push_back(nbEvents_eeu_7);
    cutstep_eeu.push_back(nbEvents_eeu_8);
    cutstep_eeu.push_back(nbEvents_eeu_9);
    
    
    cutstep_uuu.push_back( nbEvents_uuu_0);
    cutstep_uuu.push_back(nbEvents_uuu_1);
    cutstep_uuu.push_back(nbEvents_uuu_2m);
    cutstep_uuu.push_back(nbEvents_uuu_3);
    cutstep_uuu.push_back(nbEvents_uuu_4);
    cutstep_uuu.push_back(nbEvents_uuu_5);
    cutstep_uuu.push_back(nbEvents_uuu_6);
    cutstep_uuu.push_back(nbEvents_uuu_7);
    cutstep_uuu.push_back(nbEvents_uuu_8);
    cutstep_uuu.push_back(nbEvents_uuu_9);
    
    
    cutstep_uue.push_back( nbEvents_uue_0);
    cutstep_uue.push_back(nbEvents_uue_1);
    cutstep_uue.push_back(nbEvents_uue_2m);
    cutstep_uue.push_back(nbEvents_uue_3);
    cutstep_uue.push_back(nbEvents_uue_4);
    cutstep_uue.push_back(nbEvents_uue_5);
    cutstep_uue.push_back(nbEvents_uue_6);
    cutstep_uue.push_back(nbEvents_uue_7);
    cutstep_uue.push_back(nbEvents_uue_8);
    cutstep_uue.push_back(nbEvents_uue_9);
    
    
    for( int i =0 ; i < 10; i++){
      cout << "cutstep " << i << " has " << cutstep[i] << " events" << endl;
    }
    
    
    
    cout << "nbEvents_0 trigg: " <<  nbEvents_0 << endl;
    cout << "trigger req check for data " << nbEvents_test << endl;
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
    
    nTrigg_eee = nbEvents_eee_0;
    n3lep_eee = nbEvents_eee_1;
    nVetoMu_eee = nbEvents_eee_1m;
    nVetoEl_eee = nbEvents_eee_2m;
    nOS_eee = nbEvents_eee_2;
    nZmass_eee = nbEvents_eee_3;
    nJet_eee = nbEvents_eee_4;
    nBJet_eee = nbEvents_eee_5;
    nMWT_eee = nbEvents_eee_6;
    nSMtop_eee = nbEvents_eee_7;
    nMET_eee = nbEvents_eee_8;
    
    nTrigg_eeu = nbEvents_eeu_0;
    n3lep_eeu = nbEvents_eeu_1;
    nVetoMu_eeu = nbEvents_eeu_1m;
    nVetoEl_eeu = nbEvents_eeu_2m;
    nOS_eeu = nbEvents_eeu_2;
    nZmass_eeu = nbEvents_eeu_3;
    nJet_eeu = nbEvents_eeu_4;
    nBJet_eeu = nbEvents_eeu_5;
    nMWT_eeu = nbEvents_eeu_6;
    nSMtop_eeu = nbEvents_eeu_7;
    nMET_eeu = nbEvents_eeu_8;
    
    nTrigg_uuu = nbEvents_uuu_0;
    n3lep_uuu = nbEvents_uuu_1;
    nVetoMu_uuu = nbEvents_uuu_1m;
    nVetoEl_uuu = nbEvents_uuu_2m;
    nOS_uuu = nbEvents_uuu_2;
    nZmass_uuu = nbEvents_uuu_3;
    nJet_uuu = nbEvents_uuu_4;
    nBJet_uuu = nbEvents_uuu_5;
    nMWT_uuu = nbEvents_uuu_6;
    nSMtop_uuu = nbEvents_uuu_7;
    nMET_uuu = nbEvents_uuu_8;
    
    nTrigg_uue = nbEvents_uue_0;
    n3lep_uue = nbEvents_uue_1;
    nVetoMu_uue = nbEvents_uue_1m;
    nVetoEl_uue = nbEvents_uue_2m;
    nOS_uue = nbEvents_uue_2;
    nZmass_uue = nbEvents_uue_3;
    nJet_uue = nbEvents_uue_4;
    nBJet_uue = nbEvents_uue_5;
    nMWT_uue = nbEvents_uue_6;
    nSMtop_uue = nbEvents_uue_7;
    nMET_uue = nbEvents_uue_8;
    
    
    
    
    
    //	for(int j_eeu = 0; j < 9; j++){       cout << cutstep[j] << endl; }
    sumW = (int) sumWeights;
    nEv = (int) nEvents;
    xsec = xSect;
    nEvPassed  = (int) nbSelectedEvents;
    nbTrig = nbEvents_0;
    globalTree->Fill();
    if(verbose > 0) cout << "end eventloop" << endl;
    
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
    if(matching){
      cout << " ******************** MATCHING INFO **************************" << endl;
      cout << " W lepton " << matchedWlep << " from " << matchedEvents_Wlep << " or " <<((double) matchedWlep / (double)matchedEvents_Wlep)*100 << " % matched" << endl;
      cout << " Z lepton " << (double) (matchedZlep_1+matchedZlep_0)/2 << " from " << matchedEvents_Zlep << " or " << ((double)(matchedZlep_1 + (double)matchedZlep_0) / (2*(double)matchedEvents_Zlep))*100 << " % matched" << endl;
       cout << " B jet " << matchedBjet << " from " << matchedEvents_Bjet<< " or " <<((double) matchedBjet / (double)matchedEvents_Bjet)*100 << " % matched" << endl;
      cout << "** for selected lepton matching events " << endl;
      cout << int_eventForZlepmatching << " events out of " << nbSelectedEvents << " could be used for matching or " << ((double) int_eventForZlepmatching / (double) nbSelectedEvents)*100 << " %" << endl;
      cout << " W lepton " << int_eventForWlepmatchingmatched << " from " << int_eventForWlepmatching << " or " <<((double) int_eventForWlepmatchingmatched / (double)int_eventForWlepmatching)*100 << " % matched" << endl;
      cout << " Z lepton " << (double) (int_eventForZlepmatchingmatched0 + int_eventForZlepmatchingmatched1)/2 << " from " << int_eventForZlepmatching << " or " << ((double)(int_eventForZlepmatchingmatched0 + (double)int_eventForZlepmatchingmatched1) / (2*(double)int_eventForZlepmatching))*100 << " % matched" << endl;
      cout << "** for selected jet matching events " << endl;
      cout << int_eventForBjetmatching << " events out of " << nbSelectedEvents << " could be used for matching or " << ((double) int_eventForBjetmatching / (double) nbSelectedEvents)*100 << " %" << endl;
      cout << " B jet " << int_eventForBjetmatchingmatched << " from " << int_eventForBjetmatching << " or " <<((double) int_eventForBjetmatchingmatched / (double)int_eventForBjetmatching)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet " << int_eventForCjetmatchingmatched << " from " << int_eventForCjetmatching << " or " <<((double) int_eventForCjetmatchingmatched / (double)int_eventForCjetmatching)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet CvsL " << int_eventForCjetmatchingmatched_CvsLtagger << " from " << int_eventForCjetmatching_CvsLtagger << " or " <<((double) int_eventForCjetmatchingmatched_CvsLtagger / (double)int_eventForCjetmatching_CvsLtagger)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet CvsB " << int_eventForCjetmatchingmatched_CvsBtagger << " from " << int_eventForCjetmatching_CvsBtagger << " or " <<((double) int_eventForCjetmatchingmatched_CvsBtagger / (double)int_eventForCjetmatching_CvsBtagger)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet C loose  " << int_eventForCjetmatchingmatched_Cloose << " from " << int_eventForCjetmatching_Cloose << " or " <<((double) int_eventForCjetmatchingmatched_Cloose / (double)int_eventForCjetmatching_Cloose)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet C medium " << int_eventForCjetmatchingmatched_Cmedium << " from " << int_eventForCjetmatching_Cmedium << " or " <<((double) int_eventForCjetmatchingmatched_Cmedium/ (double)int_eventForCjetmatching_Cmedium)*100 << " % matched" << endl;
      if(Usettbar && !istZq) cout << " FCNC jet C tight  " << int_eventForCjetmatchingmatched_Ctight << " from " << int_eventForCjetmatching_Ctight<< " or " <<((double) int_eventForCjetmatchingmatched_Ctight/ (double)int_eventForCjetmatching_Ctight)*100 << " % matched" << endl;
      cout << " ******************************************************************" << endl;
    }
    tupfile->cd();
    myTree->Write();
    globalTree->Write();
    baselineTree->Write();
    tupfile->Close();

    delete tupfile;
    if(!isData && !btagShape){
      delete btwt;
      //histoFileHandle->Close();
      //delete histoFileHandle;
    }
    
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


/////////////////////////////////////// FUNCTIONS  //////////////////////////////////////////////



string ConvertIntToString(int Number, bool pad){
  ostringstream convert;
  convert.clear();
  if ( pad && Number < 10 ) { convert << std::setw(2) << std::setfill('0');}
  convert << Number;
  return convert.str();
};
string MakeTimeStamp(){
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


double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy){
  
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
  
};


///////////////// FCNC JET
int FCNCjetCalculator(std::vector<TRootPFJet*> Jets, TLorentzVector recoZ ,int index, int verb){
  
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -5;
  if(Jets.size() > 0){
    //cout << " non bjets: " << nonBJets.size() << " possibilities " <<endl;
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
    }
  }
  else{
    NbInColl = -5,
    cout << "no cjets available" << endl;
  }
  return NbInColl;
};
int FCNCjetCalculatorCvsBTagger(std::vector<TRootPFJet*> Jets,int index, int verb){
  
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
        //if((Jets[iJ]->ctag_pfCombinedCvsLJetTags()+Jets[iJ]->ctag_pfCombinedCvsBJetTags())>=(Jets[kJ]->ctag_pfCombinedCvsLJetTags()+Jets[kJ]->ctag_pfCombinedCvsBJetTags())) NbInColl = iJ;
        //else NbInColl = kJ;
        if((Jets[iJ]->ctag_pfCombinedCvsBJetTags())>=(Jets[kJ]->ctag_pfCombinedCvsBJetTags())) NbInColl = iJ;
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
int FCNCjetCalculatorCvsLTagger(std::vector<TRootPFJet*> Jets,int index, int verb){
  
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
        //if((Jets[iJ]->ctag_pfCombinedCvsLJetTags()+Jets[iJ]->ctag_pfCombinedCvsBJetTags())>=(Jets[kJ]->ctag_pfCombinedCvsLJetTags()+Jets[kJ]->ctag_pfCombinedCvsBJetTags())) NbInColl = iJ;
        //else NbInColl = kJ;
        if((Jets[iJ]->ctag_pfCombinedCvsLJetTags())>=(Jets[kJ]->ctag_pfCombinedCvsLJetTags())) NbInColl = iJ;
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
int FCNCjetCalculatorCwp(std::vector<TRootPFJet*> Jets, std::vector <int> cjetindex, int index, int verb){
  
  double TempMinMass = 100000.00;
  double TopMass = 172.9;
  TLorentzVector Jetcandidate;
  int NbInColl = -5;
  bool isCjet = false;
  if(Jets.size() > 2){
    // cout << " jets: " << Jets.size() << " possibilities " <<endl;  ;
    for( int iJ = 0; iJ < Jets.size()-1; iJ++)
    {
      if(iJ == index) continue;
      for(int iC = 0; iC < cjetindex.size() ; iC++){
        if(iJ == cjetindex[iJ]) isCjet = true;
      }
      if(!isCjet) continue;
      for(int kJ = 1; kJ < Jets.size(); kJ++){
        if(kJ == index) continue;
        for(int iC = 0; iC < cjetindex.size() ; iC++){
          if(kJ == cjetindex[kJ]) isCjet = true;
        }
        if(!isCjet) continue;
        if(Jets[iJ]->Pt()>=Jets[kJ]->Pt()) NbInColl = iJ;
        else NbInColl = kJ;
      }
    }
  }
  else if(Jets.size() == 2 && cjetindex.size() > 0){
    if(index == 0) NbInColl = 1;
    if(index == 1) NbInColl = 0;
  }
  else{
    NbInColl = -5;
   // cout << "no cjets available" << endl;
  }
  return NbInColl;
};

//////////////////// Z leptons
pair< vector <TLorentzVector> , vector < pair < string , int > > > LeptonAssigner(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons){
  //  cout << " in assigner " << endl;
  pair< vector <TLorentzVector> , vector < pair < string , int > > >  Returner;
  vector < pair < string , int > > Indices;
  Indices.clear();
  vector<TLorentzVector> ReturnColl;
  ReturnColl.clear();
  bool Assigned = false;
  
  if(electrons.size() + muons.size() != 3){
    cout << " WARNING: not 3 leptons " << endl;
    cout << "muons " << muons.size() << " electrons " << electrons.size() << endl;
    return Returner;
  }
  //elecbool = false;
  //mubool = false;
  int ZelecIndice_0 = -999;
  int ZelecIndice_1 = -999;
  int ZmuIndice_0 = -999;
  int ZmuIndice_1 = -999;
  int WelecIndice = -999;
  int WmuIndice = -999;
  //cout << " in 3 lep " << endl;
  
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
      ZelecIndice_0 = 0;
      ZelecIndice_1 = 1;
      WmuIndice = 0;
    }
    else Assigned = false;
  }
  else if(muons.size() == 2){
    //    cout << "2 muons" << endl;
    if(muons[0]->charge() != muons[1]->charge()){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Wlepcan.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Assigned = true;
      
      ZmuIndice_0 = 0;
      ZmuIndice_1 = 1;
      WelecIndice = 0;
    }
    else Assigned = false;
  }
  else if(electrons.size() ==3){
    //cout << " 3 electrons " << endl;
    bool can01 = false;
    bool can02= false;
    bool can12 = false;
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
      ZelecIndice_0 = 0;ZelecIndice_1 = 1;WelecIndice = 2;
    }
    
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Wlepcan.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Assigned = true;
      ZelecIndice_0 = 0; ZelecIndice_1=2;WelecIndice = 1;
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Wlepcan.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Assigned = true;
      ZelecIndice_0 = 1; ZelecIndice_1=2;WelecIndice = 0;
    }
    else Assigned = false;
  }
  else if(muons.size() == 3){
    bool can01 = false;
    bool can02= false;
    bool can12 = false;
    
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
      ZmuIndice_0 = 0; ZmuIndice_1=1; WmuIndice = 2;
    }
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Wlepcan.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Assigned = true;
      ZmuIndice_0 = 0; ZmuIndice_1=2;WmuIndice = 1;
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Wlepcan.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Assigned = true;
      ZmuIndice_0 = 1; ZmuIndice_1=2;WmuIndice = 0;
    }
    else Assigned = false;
  }
  if(Assigned){
    ReturnColl.push_back(Zlepcan0);
    ReturnColl.push_back(Zlepcan1);
    ReturnColl.push_back(Wlepcan);
    
    Indices.push_back( pair< string, int > ("Wmu", WmuIndice));
    Indices.push_back( pair< string, int > ("Wel", WelecIndice));
    Indices.push_back( pair< string, int > ("Zel_1", ZelecIndice_1));
    Indices.push_back( pair< string, int > ("Zel_0", ZelecIndice_0));
    Indices.push_back( pair< string, int > ("Zmu_1", ZmuIndice_1));
    Indices.push_back( pair< string, int > ("Zmu_0", ZmuIndice_0));
    
  }
  
  Returner = pair < vector < TLorentzVector > , vector < pair < string, int > > > (ReturnColl, Indices);
  return Returner;
};
/*vector <TLorentzVector> LeptonAssignerv2(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons)
 {
 // cout << " in assigner " << endl;
 vector<TLorentzVector> ReturnColl;
 Assigned = false;
 
 
 
 // elecbool = false;
 // mubool = false;
 ZelecIndice_0 = -999;
 ZelecIndice_1 = -999;
 ZmuIndice_1 = -999;
 ZmuIndice_0 = -999;
 
 TLorentzVector Zlepcan0;
 Zlepcan0.SetPxPyPzE(0.,0.,0.,0.);
 TLorentzVector Zlepcan1;
 Zlepcan1.SetPxPyPzE(0.,0.,0.,0.);
 
 
 if(electrons.size() == 2){
 //  cout << "2 electr " << electrons[0]->charge() << " " << electrons[1]->charge() << endl;
 if(electrons[0]->charge() != electrons[1]->charge()){
 Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
 Zlepcan1.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
 
 Assigned = true;
 // elecbool = true;
 ZelecIndice_0 = 0;
 ZelecIndice_1 = 1;
 //cout << "assigned " <<endl;
 }
 }
 else if(muons.size() == 2){
 //    cout << "2 muons" << endl;
 if(muons[0]->charge() != muons[1]->charge()){
 Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
 Zlepcan1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
 
 Assigned = true;
 // mubool = true;
 ZmuIndice_0 = 0;
 ZmuIndice_1 = 1;
 
 }
 }
 else if(electrons.size() ==3){
 //    cout << " 3 electrons " << endl;
 bool can01 = false;
 bool can02= false;
 bool can12 = false;
 // elecbool = true;
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
 Assigned = true;
 ZelecIndice_0 = 0;ZelecIndice_1 = 1;
 }
 else if(mass02 <= mass12 && mass02 < mass01){
 Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
 Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
 Assigned = true;
 ZelecIndice_0 = 0; ZelecIndice_1=2;
 }
 else if(mass12 < mass01 && mass12 < mass02){
 Zlepcan0.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
 Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
 Assigned = true;
 ZelecIndice_0 = 1; ZelecIndice_1=2;
 }
 }
 else if(muons.size() == 3){
 bool can01 = false;
 bool can02= false;
 bool can12 = false;
 // mubool = true;
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
 Assigned = true;
 ZmuIndice_0 = 0; ZmuIndice_1=1;
 }
 else if(mass02 <= mass12 && mass02 < mass01){
 Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
 Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
 Assigned = true;
 ZmuIndice_0 = 0; ZmuIndice_1=2;
 }
 else if(mass12 < mass01 && mass12 < mass02){
 Zlepcan0.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
 Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
 Assigned = true;
 ZmuIndice_0 = 1; ZmuIndice_1=2;
 }
 }
 if(Assigned){
 ReturnColl.push_back(Zlepcan0);
 ReturnColl.push_back(Zlepcan1);
 //   cout << "filled" << endl;
 
 }
 if(!Assigned){
 //    cout << " WARNING: leptons not set for assignment " << endl;
 return ReturnColl;
 }
 
 //cout << "returned" << endl;
 return ReturnColl;
 }
 */



///////////////// MET
TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met){
  
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
  
  
};

/////////////////SM B
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

//////////////////// GEN INFO
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int SMtau = -999;
  int FCNCmuPlus = -999;
  int FCNCmuMin = -999;
  int FCNCelPlus = -999;
  int FCNCelMin = -999;
  int FCNCtauMin = -999;
  int FCNCtauPlus = -999;
  int FCNCZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  int FCNCq = -999;
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMtauATop = false;
  bool SMtauTop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  bool FCNCZATop = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool FCNCmuPlusFound = false;
  bool FCNCelPlusFound = false;
  bool FCNCelMinFound = false;
  bool FCNCmuMinFound = false;
  bool FCNCtau = false;
  bool FCNCZATopEl = false;
  bool FCNCZTopEl = false;
  bool FCNCZTopMu = false;
  bool FCNCZATopMu = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  bool FCNCqATop = false;
  bool FCNCqTop = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  // plots for mcParticles
  
  //cout << "event " << evt_num_ << endl;
  // search for the right events
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << "found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    else if( abs(mcParticles_[iMC]->type()) ==  5  && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    
    
    //FCNC
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqATop = true;
      
      //cout << "q from tbar found  " << endl;
    } // q  from tbar
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqTop = true;
      
      // cout << "q from t found  " << endl;
    } // q  from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      FCNCZATop = true;
      FCNCmuMinFound = true;
      
      // cout << "mu from Z from tbar found  " << endl;
    } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus= iMC;}
      FCNCZATop = true;
      FCNCmuPlusFound = true;
      
      //cout << "mu from Z from tbar found  " << endl;
    } // mu + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus = iMC;}
      FCNCmuMinFound = true;
      
      //cout << "mu from Z from t found  " << endl;
    } // mu+ from Z from top
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      
      FCNCmuPlusFound = true;
      
      // cout << "mu from Z from t found  " << endl;
    } // mu- from Z from top
    
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      FCNCZATop = true;
      FCNCelMinFound = true;
      
      //cout << "el from Z from tbar found  " << endl;
    } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      FCNCZATop = true;
      FCNCelPlusFound = true;
      
      // cout << "el from Z from tbar found  " << endl;
    } // el + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      
      FCNCelPlusFound = true;
      
      //cout << "el from Z from t found  " << endl;
    } // el+ from Z from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      
      FCNCelMinFound = true;
      
      // cout << "el from Z from t found  " << endl;
    } // el- from Z from top
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopMu = true;
    } //Z from top with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopEl = true;
    } //Z from top with el decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC;nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopMu = true;
    } //Z from atop with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopEl = true;
    } //Z from atop with el decay
    
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool FCNCZfound = false;
  if(FCNCZATopEl || FCNCZTopEl || FCNCZATopMu || FCNCZTopMu) {FCNCZfound = true; }
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool FCNCmufound = false;
  if(FCNCmuPlusFound && FCNCmuMinFound){ FCNCmufound = true; }
  bool FCNCelfound = false;
  if(FCNCelPlusFound && FCNCelMinFound){ FCNCelfound = true; }
  bool FCNCqfound = false;
  if(FCNCqATop || FCNCqTop ){FCNCqfound = true; }
  
  bool SMTopFCNCATop = false;
  if((SMmuTop|| SMelTop) && (FCNCZATopMu || FCNCZATopEl)){ SMTopFCNCATop = true; }
  bool SMATopFCNCTop = false;
  if((SMmuATop|| SMelATop) && (FCNCZTopMu || FCNCZTopEl) && SMbfound && FCNCqfound ){ SMATopFCNCTop = true; }
  bool foundDecay = false;
  // if( SMTopFCNCATop || SMATopFCNCTop ){ foundDecay = true;}
  
  
  if(( FCNCmufound || FCNCelfound ) && (SMmufound || SMelfound) && SMbfound && FCNCqfound){
    foundDecay = true;
    //cout << "found decay" << endl;
  }
  if(!foundDecay) {
    
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  
  // begin matching
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  //cout << "event " << evt_num_ << endl;
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    
    if( abs(mcParticles_[iMC]->type()) <=  5 || abs(mcParticles_[iMC]->type()) == 21 ){
     partons.push_back(mcpartTLV[iMC]);
     partonID.push_back(iMC);
     
    }
    
  }
  
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedjets,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedjets collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 2){
    // cout << "ERROR " << JetPartonPair.size() <<   endl;
    // cout << "partons " << partons.size() << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  // find objects to be matched
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMb");
    } // b from t
    //FCNC
    if( (abs(mcParticles_[partonID[partonIDnb]]->type()) ==  2 || abs(mcParticles_[partonID[partonIDnb]]->type()) ==  4 )&& abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCq");
    } // q from t
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempQ;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempQ.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
    if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedjets[PPair[iPart].first].Px(), selectedjets[PPair[iPart].first].Py(), selectedjets[PPair[iPart].first].Pz(), selectedjets[PPair[iPart].first].E());}
    if(NPair[iPart].find("FCNCq")!=string::npos){ tempQ.SetPxPyPzE(selectedjets[PPair[iPart].first].Px(), selectedjets[PPair[iPart].first].Py(), selectedjets[PPair[iPart].first].Pz(), selectedjets[PPair[iPart].first].E());  }
  }
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int SMtau = -999;
  int FCNCmuPlus = -999;
  int FCNCmuMin = -999;
  int FCNCelPlus = -999;
  int FCNCelMin = -999;
  int FCNCtauMin = -999;
  int FCNCtauPlus = -999;
  int FCNCZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  int FCNCq = -999;
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMtauATop = false;
  bool SMtauTop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  bool FCNCZATop = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool FCNCmuPlusFound = false;
  bool FCNCelPlusFound = false;
  bool FCNCelMinFound = false;
  bool FCNCmuMinFound = false;
  bool FCNCtau = false;
  bool FCNCZATopEl = false;
  bool FCNCZTopEl = false;
  bool FCNCZTopMu = false;
  bool FCNCZATopMu = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  bool FCNCqATop = false;
  bool FCNCqTop = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  

  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  // plots for mcParticles
  
  //cout << "event " << evt_num_ << endl;
  // search for the right events
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << "found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    else if( abs(mcParticles_[iMC]->type()) ==  5  && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    
    
    //FCNC
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqATop = true;
      
      //cout << "q from tbar found  " << endl;
    } // q  from tbar
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqTop = true;
      
      // cout << "q from t found  " << endl;
    } // q  from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      FCNCZATop = true;
      FCNCmuMinFound = true;
      
      // cout << "mu from Z from tbar found  " << endl;
    } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus= iMC;}
      FCNCZATop = true;
      FCNCmuPlusFound = true;
      
      //cout << "mu from Z from tbar found  " << endl;
    } // mu + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus = iMC;}
      FCNCmuMinFound = true;
      
      //cout << "mu from Z from t found  " << endl;
    } // mu+ from Z from top
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      
      FCNCmuPlusFound = true;
      
      // cout << "mu from Z from t found  " << endl;
    } // mu- from Z from top
    
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      FCNCZATop = true;
      FCNCelMinFound = true;
      
      //cout << "el from Z from tbar found  " << endl;
    } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      FCNCZATop = true;
      FCNCelPlusFound = true;
      
      // cout << "el from Z from tbar found  " << endl;
    } // el + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      
      FCNCelPlusFound = true;
      
      //cout << "el from Z from t found  " << endl;
    } // el+ from Z from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      
      FCNCelMinFound = true;
      
      // cout << "el from Z from t found  " << endl;
    } // el- from Z from top
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopMu = true;
    } //Z from top with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopEl = true;
    } //Z from top with el decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC;nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopMu = true;
    } //Z from atop with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopEl = true;
    } //Z from atop with el decay
    
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool FCNCZfound = false;
  if(FCNCZATopEl || FCNCZTopEl || FCNCZATopMu || FCNCZTopMu) {FCNCZfound = true; }
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool FCNCmufound = false;
  if(FCNCmuPlusFound && FCNCmuMinFound){ FCNCmufound = true; }
  bool FCNCelfound = false;
  if(FCNCelPlusFound && FCNCelMinFound){ FCNCelfound = true; }
  bool FCNCqfound = false;
  if(FCNCqATop || FCNCqTop ){FCNCqfound = true; }
  
  bool SMTopFCNCATop = false;
  if((SMmuTop|| SMelTop) && (FCNCZATopMu || FCNCZATopEl)){ SMTopFCNCATop = true; }
  bool SMATopFCNCTop = false;
  if((SMmuATop|| SMelATop) && (FCNCZTopMu || FCNCZTopEl) && SMbfound && FCNCqfound ){ SMATopFCNCTop = true; }
  bool foundDecay = false;
  // if( SMTopFCNCATop || SMATopFCNCTop ){ foundDecay = true;}
  
  
  if(( FCNCmufound || FCNCelfound ) && (SMmufound || SMelfound) && SMbfound && FCNCqfound){
    foundDecay = true;
    //cout << "found decay" << endl;
  }
  if(!foundDecay) {
    
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  // theta
  // cos theta check
  if(SMbfound && (SMmufound || SMelfound) && (SMnuelfound || SMnumufound)){
    TLorentzVector tempLep;
    //cout << "SMmu " << SMmu << " SMel " << SMel << " SMnumu " << SMnumu << " SMnuel " << SMnuel << " SMb " << SMb << endl;
    
    if(SMmufound ){   tempLep.SetPxPyPzE(mcpartTLV[SMmu].Px(), mcpartTLV[SMmu].Py(), mcpartTLV[SMmu].Pz(), mcpartTLV[SMmu].Energy());}
    if(SMelfound ){   tempLep.SetPxPyPzE(mcpartTLV[SMel].Px(), mcpartTLV[SMel].Py(), mcpartTLV[SMel].Pz(), mcpartTLV[SMel].Energy());}
    TLorentzVector tempB;
    tempB.SetPxPyPzE(mcpartTLV[SMb].Px(), mcpartTLV[SMb].Py(), mcpartTLV[SMb].Pz(), mcpartTLV[SMb].Energy());
    
    TLorentzVector tempNu;
    if(SMnumufound ){  tempNu.SetPxPyPzE(mcpartTLV[SMnumu].Px(), mcpartTLV[SMnumu].Py(), mcpartTLV[SMnumu].Pz(), mcpartTLV[SMnumu].Energy());}
    if(SMnuelfound ){  tempNu.SetPxPyPzE(mcpartTLV[SMnuel].Px(), mcpartTLV[SMnuel].Py(), mcpartTLV[SMnuel].Pz(), mcpartTLV[SMnuel].Energy());}
    
    //cout << "templep " << tempLep.Pt() << " x " << tempLep.Px() << " y " << tempLep.Py() << " z " << tempLep.Pz() << endl;
    //cout << " tempB " << tempB.Pt() << " x " << tempB.Px() <<" y " << tempB.Py() << " z " << tempB.Pz() << endl;
    //cout <<  " tempNu " << tempNu.Pt() << " x " << tempNu.Px() << " y " << tempNu.Py() << " z " << tempNu.Pz() << " E " << tempNu.Energy() << " Phi " << tempNu.Phi() << " Eta " << tempNu.Eta()<< endl;
    pair <float, float> tempcosPair = CosThetaCalculation(tempLep, tempNu, tempB, true);
    
  }
  
  
  
  //SM TOP
  
  if( SMbfound && SMWfound){
    histo1D["Topmass_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).M());
    histo1D["pt_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Pt());
    histo1D["eta_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Eta());
    histo1D["phi_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Phi());
    histo1D["dPhi_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb],mcpartTLV[SMW]));
    histo1D["dR_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb],mcpartTLV[SMW]));
    
  }
  if( SMbATop && SMWATop && antitopfound){
    //cout << "in tbar" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[antitopQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[antitopQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[antitopQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[antitopQ].Phi());
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Eta() );
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    
  }
  if( SMbTop && SMWTop && topfound){
    //cout << "in top" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[topQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[topQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[topQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[topQ].Phi());
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Eta() );
    
  }
  if( SMbfound && SMnumufound && SMmufound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    
    
  }
  if( SMbfound && SMnuelfound && SMelfound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    
    
  }
  
  
  //FCNC TOP
  if( FCNCqfound && FCNCZfound ){
    histo1D["Topmass_Zq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCZ]).M());
    histo1D["pt_Zq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCZ]).Pt());
    histo1D["eta_Zq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCZ]).Eta());
    histo1D["phi_Zq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCZ]).Phi());
    
    histo1D["dPhi_Zq"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCq],mcpartTLV[FCNCZ]));
    histo1D["dR_Zq"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCq],mcpartTLV[FCNCZ]));
    
  }
  if( FCNCqATop && (FCNCZATopMu || FCNCZATopEl) && antitopfound){
    //cout << "in tbar" << endl;
    histo1D["Topmass_fcnctq"]->Fill(mcpartTLV[antitopQ].M());
    histo1D["pt_fcnctq"]->Fill(mcpartTLV[antitopQ].Pt());
    histo1D["eta_fcnctq"]->Fill(mcpartTLV[antitopQ].Eta());
    histo1D["phi_fcnctq"]->Fill(mcpartTLV[antitopQ].Phi());
    
    histo2D["Topmass_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[antitopQ].M() );
    histo2D["pt_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[antitopQ].Pt() );
    histo2D["phi_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[antitopQ].Phi() );
    histo2D["eta_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[antitopQ].Eta() );
    
    histo1D["dPhi_Zqtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCq]+mcpartTLV[FCNCZ],mcpartTLV[antitopQ]));
    histo1D["dR_Zqtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCq]+mcpartTLV[FCNCZ],mcpartTLV[antitopQ]));
    
  }
  
  if( FCNCqTop && (FCNCZTopMu || FCNCZTopEl) && topfound){
    // cout << "in top" << endl;
    // cout << "topQ " << topQ << " FCNCq " << FCNCq << " FCNCZ " << FCNCZ << endl;
    histo1D["Topmass_fcnctq"]->Fill(mcpartTLV[topQ].M());
    histo1D["pt_fcnctq"]->Fill(mcpartTLV[topQ].Pt());
    histo1D["eta_fcnctq"]->Fill(mcpartTLV[topQ].Eta());
    histo1D["phi_fcnctq"]->Fill(mcpartTLV[topQ].Phi());
    
    histo1D["dPhi_Zqtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCq]+mcpartTLV[FCNCZ],mcpartTLV[topQ]));
    histo1D["dR_Zqtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCq]+mcpartTLV[FCNCZ],mcpartTLV[topQ]));
    
    
    histo2D["Topmass_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[topQ].M() );
    histo2D["pt_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[topQ].Pt() );
    histo2D["phi_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[topQ].Phi() );
    histo2D["eta_top_Zq"]->Fill((mcpartTLV[FCNCZ] + mcpartTLV[FCNCq]).M(),mcpartTLV[topQ].Eta() );
    
  }
  if( FCNCqfound &&  FCNCmufound){
    histo1D["Topmass_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCmuPlus]+ mcpartTLV[FCNCmuMin]).M());
    histo1D["pt_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCmuPlus]+ mcpartTLV[FCNCmuMin]).Pt());
    histo1D["eta_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCmuPlus]+ mcpartTLV[FCNCmuMin]).Eta());
    histo1D["phi_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCmuPlus]+ mcpartTLV[FCNCmuMin]).Phi());
    
    histo1D["dPhi_llq"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCq] , mcpartTLV[FCNCmuPlus]+ mcpartTLV[ FCNCmuMin]));
    histo1D["dR_llq"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCq] , mcpartTLV[FCNCmuPlus]+ mcpartTLV[ FCNCmuMin]));
    
    
  }
  if( FCNCqfound &&  FCNCelfound){
    histo1D["Topmass_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[ FCNCelPlus]+ mcpartTLV[ FCNCelMin]).M());
    histo1D["pt_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCelPlus]+ mcpartTLV[FCNCelMin]).Pt());
    histo1D["eta_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCelPlus]+ mcpartTLV[FCNCelMin]).Eta());
    histo1D["phi_llq"]->Fill((mcpartTLV[FCNCq] + mcpartTLV[FCNCelPlus]+ mcpartTLV[FCNCelMin]).Phi());
    
    histo1D["dPhi_llq"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCq] , mcpartTLV[ FCNCelPlus]+ mcpartTLV[FCNCelMin]));
    histo1D["dR_llq"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCq] , mcpartTLV[FCNCelPlus]+ mcpartTLV[FCNCelMin]));
    
    
  }
  
  if( FCNCmufound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[FCNCmuMin].Pt(), mcpartTLV[FCNCmuPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[FCNCmuMin].Phi(), mcpartTLV[FCNCmuPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[FCNCmuMin].Eta(), mcpartTLV[FCNCmuPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCmuMin],mcpartTLV[FCNCmuPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCmuMin],mcpartTLV[FCNCmuPlus]));
    histo1D["mass_lep1"]->Fill(mcpartTLV[FCNCmuMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[FCNCmuPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[FCNCmuMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[FCNCmuMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[FCNCmuMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[FCNCmuMin].M(), mcpartTLV[FCNCmuPlus].M());
    
  }
  if( FCNCelfound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[FCNCelMin].Pt(), mcpartTLV[FCNCelPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[FCNCelMin].Phi(), mcpartTLV[FCNCelPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[FCNCelMin].Eta(), mcpartTLV[FCNCelPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCelMin],mcpartTLV[FCNCelPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCelMin],mcpartTLV[FCNCelPlus]));
    
    histo1D["mass_lep1"]->Fill(mcpartTLV[FCNCelMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[FCNCelPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[FCNCelMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[FCNCelPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[FCNCelMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[FCNCelPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[FCNCelMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[FCNCelPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[FCNCelMin].M(), mcpartTLV[FCNCelPlus].M());
    
  }
  
  if( FCNCZfound){
    histo1D["Zmass_Zbos"]->Fill(mcpartTLV[FCNCZ].M());
    histo1D["pt_Zbos"]->Fill(mcpartTLV[FCNCZ].Pt());
    histo1D["phi_Zbos"]->Fill(mcpartTLV[FCNCZ].Phi());
    histo1D["eta_Zbos"]->Fill(mcpartTLV[FCNCZ].Eta());
    histo1D["mc_nZLep"]->Fill(nbZDaughters);
    histo1D["mc_nZMu"]->Fill(nbZDaughters);
  }
  if( FCNCZfound && FCNCmufound){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).M(),mcpartTLV[FCNCZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Eta() );
  }
  
  if( FCNCelfound && FCNCZfound ){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).M(),mcpartTLV[FCNCZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Eta() );
    
  }
  
  if( SMWfound && SMelfound ){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWEl"]->Fill(nbWDaughters);
    
  }
  if( SMWfound && SMmufound){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWMu"]->Fill(nbWDaughters);
  }
  if( SMWfound && FCNCZfound){
    histo2D["nZbosonnWboson"]->Fill(nbZDaughters, nbWDaughters);
  }
  
  if( FCNCqfound)
  {
    histo1D["mass_FCNCq"]->Fill(mcpartTLV[FCNCq].M());
    histo1D["pt_FCNCq"]->Fill(mcpartTLV[FCNCq].Pt());
    histo1D["phi_FCNCq"]->Fill(mcpartTLV[FCNCq].Phi());
    histo1D["eta_FCNCq"]->Fill(mcpartTLV[FCNCq].Eta());
    
  }
  
  
  // end plots for mc particles
  // begin matching
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  //cout << "event " << evt_num_ << endl;
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    
/*if( abs(mcParticles_[iMC]->type()) <=  5 || abs(mcParticles_[iMC]->type()) == 21 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    }
    else*/
    if( abs(mcParticles_[iMC]->type()) ==  13 ||  abs(mcParticles_[iMC]->type()) ==  11 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    }
  }
  
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedleptons,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedleptons collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 3){
   // cout << "ERROR " << JetPartonPair.size() <<   endl;
   // cout << "partons " << partons.size() << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  // find objects to be matched
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
   /* if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMb");
    } // b from t */
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMmu");
    } // mu from W from t
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMel");
    } // el from W from t
    
    
    //FCNC
    /*if( (abs(mcParticles_[partonID[partonIDnb]]->type()) ==  2 || abs(mcParticles_[partonID[partonIDnb]]->type()) ==  4 )&& abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCq");
    } // q from t*/
    if( mcParticles_[partonID[partonIDnb]]->type() ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmumin");
    } // mu- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelmin");
    } // el- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmuplus");
    } // mu+ from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelplus");
    } // el+ from Z from t
    
    
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempQ;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempQ.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
   /* if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());}*/
    if(NPair[iPart].find("SMmu")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("SMel")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());}
    
    if(NPair[iPart].find("FCNCmumin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());   }
    if(NPair[iPart].find("FCNCelmin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());  }
    if(NPair[iPart].find("FCNCmuplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());   }
    if(NPair[iPart].find("FCNCelplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
   /* if(NPair[iPart].find("FCNCq")!=string::npos){ tempQ.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());  }*/
  }
  histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
 // histo1D["matchedFCNCTopmass"]->Fill((tempZlepm+tempZlepp+tempQ).M());
 // histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  ObjectMatcher(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedobjects){ // to be checked
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int SMtau = -999;
  int FCNCmuPlus = -999;
  int FCNCmuMin = -999;
  int FCNCelPlus = -999;
  int FCNCelMin = -999;
  int FCNCtauMin = -999;
  int FCNCtauPlus = -999;
  int FCNCZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  int FCNCq = -999;
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMtauATop = false;
  bool SMtauTop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  bool FCNCZATop = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool FCNCmuPlusFound = false;
  bool FCNCelPlusFound = false;
  bool FCNCelMinFound = false;
  bool FCNCmuMinFound = false;
  bool FCNCtau = false;
  bool FCNCZATopEl = false;
  bool FCNCZTopEl = false;
  bool FCNCZTopMu = false;
  bool FCNCZATopMu = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  bool FCNCqATop = false;
  bool FCNCqTop = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  // plots for mcParticles
  
  //cout << "event " << evt_num_ << endl;
  // search for the right events
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << "found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    else if( abs(mcParticles_[iMC]->type()) ==  5  && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    
    
    //FCNC
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqATop = true;
      
      //cout << "q from tbar found  " << endl;
    } // q  from tbar
    else if(( abs(mcParticles_[iMC]->type()) ==  2 ||  abs(mcParticles_[iMC]->type()) ==  4) && mcParticles_[iMC]->motherType()  ==  6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCq== -999) {FCNCq = iMC;}
      FCNCqTop = true;
      
      // cout << "q from t found  " << endl;
    } // q  from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      FCNCZATop = true;
      FCNCmuMinFound = true;
      
      // cout << "mu from Z from tbar found  " << endl;
    } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus= iMC;}
      FCNCZATop = true;
      FCNCmuPlusFound = true;
      
      //cout << "mu from Z from tbar found  " << endl;
    } // mu + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus = iMC;}
      FCNCmuMinFound = true;
      
      //cout << "mu from Z from t found  " << endl;
    } // mu+ from Z from top
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      
      FCNCmuPlusFound = true;
      
      // cout << "mu from Z from t found  " << endl;
    } // mu- from Z from top
    
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      FCNCZATop = true;
      FCNCelMinFound = true;
      
      //cout << "el from Z from tbar found  " << endl;
    } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      FCNCZATop = true;
      FCNCelPlusFound = true;
      
      // cout << "el from Z from tbar found  " << endl;
    } // el + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      
      FCNCelPlusFound = true;
      
      //cout << "el from Z from t found  " << endl;
    } // el+ from Z from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      
      FCNCelMinFound = true;
      
      // cout << "el from Z from t found  " << endl;
    } // el- from Z from top
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopMu = true;
    } //Z from top with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopEl = true;
    } //Z from top with el decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC;nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopMu = true;
    } //Z from atop with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == -6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZATopEl = true;
    } //Z from atop with el decay
    
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool FCNCZfound = false;
  if(FCNCZATopEl || FCNCZTopEl || FCNCZATopMu || FCNCZTopMu) {FCNCZfound = true; }
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool FCNCmufound = false;
  if(FCNCmuPlusFound && FCNCmuMinFound){ FCNCmufound = true; }
  bool FCNCelfound = false;
  if(FCNCelPlusFound && FCNCelMinFound){ FCNCelfound = true; }
  bool FCNCqfound = false;
  if(FCNCqATop || FCNCqTop ){FCNCqfound = true; }
  
  bool SMTopFCNCATop = false;
  if((SMmuTop|| SMelTop) && (FCNCZATopMu || FCNCZATopEl)){ SMTopFCNCATop = true; }
  bool SMATopFCNCTop = false;
  if((SMmuATop|| SMelATop) && (FCNCZTopMu || FCNCZTopEl) && SMbfound && FCNCqfound ){ SMATopFCNCTop = true; }
  bool foundDecay = false;
  // if( SMTopFCNCATop || SMATopFCNCTop ){ foundDecay = true;}
  
  
  if(( FCNCmufound || FCNCelfound ) && (SMmufound || SMelfound) && SMbfound && FCNCqfound){
    foundDecay = true;
    //cout << "found decay" << endl;
  }
  if(!foundDecay) {
    
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  

  
  
  // end plots for mc particles
  // begin matching
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  //cout << "event " << evt_num_ << endl;
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    
    if( abs(mcParticles_[iMC]->type()) <=  5 || abs(mcParticles_[iMC]->type()) == 21 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    }
    else if( abs(mcParticles_[iMC]->type()) ==  13 ||  abs(mcParticles_[iMC]->type()) ==  11 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    }
  }
  
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedobjects,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedobjects collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 5){
    // cout << "ERROR " << JetPartonPair.size() <<   endl;
    // cout << "partons " << partons.size() << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  // find objects to be matched
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
     PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
     NPair.push_back("SMb");
     } // b from t
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMmu");
    } // mu from W from t
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMel");
    } // el from W from t
    
    
    //FCNC
    if( (abs(mcParticles_[partonID[partonIDnb]]->type()) ==  2 || abs(mcParticles_[partonID[partonIDnb]]->type()) ==  4 )&& abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
     PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
     NPair.push_back("FCNCq");
     } // q from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmumin");
    } // mu- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelmin");
    } // el- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmuplus");
    } // mu+ from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelplus");
    } // el+ from Z from t
    
    
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempQ;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempQ.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
    if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());}
    if(NPair[iPart].find("SMmu")!=string::npos){ tempWlep.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E()); }
    if(NPair[iPart].find("SMel")!=string::npos){ tempWlep.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());}
    
    if(NPair[iPart].find("FCNCmumin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());   }
    if(NPair[iPart].find("FCNCelmin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());  }
    if(NPair[iPart].find("FCNCmuplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());   }
    if(NPair[iPart].find("FCNCelplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E()); }
     if(NPair[iPart].find("FCNCq")!=string::npos){ tempQ.SetPxPyPzE(selectedobjects[PPair[iPart].first].Px(), selectedobjects[PPair[iPart].first].Py(), selectedobjects[PPair[iPart].first].Pz(), selectedobjects[PPair[iPart].first].E());  }
  }
  //histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
  histo1D["matchedFCNCTopmass"]->Fill((tempZlepm+tempZlepp+tempQ).M());
  histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatcherST(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int FCNCmuPlus = -999;
  int FCNCmuMin = -999;
  int FCNCelPlus = -999;
  int FCNCelMin = -999;
  
  int FCNCZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMtauATop = false;
  bool SMtauTop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool FCNCmuPlusFound = false;
  bool FCNCelPlusFound = false;
  bool FCNCelMinFound = false;
  bool FCNCmuMinFound = false;
  
  bool FCNCZTopEl = false;
  bool FCNCZTopMu = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  
  // Plots of mc particles
  
  //cout << "event " << evt_num_ << endl;
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << "found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      //cout << "found mu from tbar" << endl;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    
    
    //FCNC
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      FCNCmuMinFound = true;
      
      // cout << "mu from Z from tbar found  " << endl;
    } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus= iMC;}
      FCNCmuPlusFound = true;
      
      //cout << "mu from Z from tbar found  " << endl;
    } // mu + from Z  from tbar
    
    
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      FCNCelMinFound = true;
      //cout << "el from Z from tbar found  " << endl;
    } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      FCNCelPlusFound = true;
      
      // cout << "el from Z from tbar found  " << endl;
    } // el + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  23 &&  abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopMu = true;
    } //Z from top with mu decay
    else if( mcParticles_[iMC]->type() ==  23 &&  abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopEl = true;
    } //Z from top with el decay
    
    
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool FCNCZfound = false;
  if(FCNCZTopEl  || FCNCZTopMu) {FCNCZfound = true; }
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool FCNCmufound = false;
  if(FCNCmuPlusFound && FCNCmuMinFound){ FCNCmufound = true; }
  bool FCNCelfound = false;
  if(FCNCelPlusFound && FCNCelMinFound){ FCNCelfound = true; }
  
  bool foundDecay = false;
  if(  (SMmufound || SMelfound) && (FCNCelfound || FCNCmufound) && SMbfound ){ foundDecay = true;}
  if(!foundDecay) {
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }

  
  
  // end plots
  // begin matching
  
  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( abs(mcParticles_[iMC]->type()) <= 5 || abs(mcParticles_[iMC]->type()) == 21){
     partons.push_back(mcpartTLV[iMC]);
     partonID.push_back(iMC);
     } // jets
    
  }
  
 // cout << "partons " << partons.size() << endl;
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedjets,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedjets collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 1){
    //cout << "ERROR" << endl;
    //cout << JetPartonPair.size() << " " << partons.size() << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMb");
    } // b from t
    
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
      if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedjets[PPair[iPart].first].Px(), selectedjets[PPair[iPart].first].Py(), selectedjets[PPair[iPart].first].Pz(), selectedjets[PPair[iPart].first].E()); }
    
  }
  
 // histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
  //histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatcherST(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int FCNCmuPlus = -999;
  int FCNCmuMin = -999;
  int FCNCelPlus = -999;
  int FCNCelMin = -999;
  
  int FCNCZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMtauATop = false;
  bool SMtauTop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool FCNCmuPlusFound = false;
  bool FCNCelPlusFound = false;
  bool FCNCelMinFound = false;
  bool FCNCmuMinFound = false;
  
  bool FCNCZTopEl = false;
  bool FCNCZTopMu = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  
  // Plots of mc particles
  
  //cout << "event " << evt_num_ << endl;
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << "found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      //cout << "found mu from tbar" << endl;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      //cout << "found el from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    
    
    //FCNC
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23 ){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuMin== -999) {FCNCmuMin = iMC;}
      FCNCmuMinFound = true;
      
      // cout << "mu from Z from tbar found  " << endl;
    } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23){
      if(mcParticles_[iMC]->status() == 23) {FCNCmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCmuPlus== -999) {FCNCmuPlus= iMC;}
      FCNCmuPlusFound = true;
      
      //cout << "mu from Z from tbar found  " << endl;
    } // mu + from Z  from tbar
    
    
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelMin== -999) {FCNCelMin = iMC;}
      FCNCelMinFound = true;
      //cout << "el from Z from tbar found  " << endl;
    } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCelPlus= iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && FCNCelPlus== -999) {FCNCelPlus = iMC;}
      FCNCelPlusFound = true;
      
      // cout << "el from Z from tbar found  " << endl;
    } // el + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  23 &&  abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopMu = true;
    } //Z from top with mu decay
    else if( mcParticles_[iMC]->type() ==  23 &&  abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      
      if(mcParticles_[iMC]->status() == 23) {FCNCZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && FCNCZ== -999) {FCNCZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      FCNCZTopEl = true;
    } //Z from top with el decay
    
    
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool FCNCZfound = false;
  if(FCNCZTopEl  || FCNCZTopMu) {FCNCZfound = true; }
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool FCNCmufound = false;
  if(FCNCmuPlusFound && FCNCmuMinFound){ FCNCmufound = true; }
  bool FCNCelfound = false;
  if(FCNCelPlusFound && FCNCelMinFound){ FCNCelfound = true; }
  
  bool foundDecay = false;
  if(  (SMmufound || SMelfound) && (FCNCelfound || FCNCmufound) && SMbfound ){ foundDecay = true;}
  if(!foundDecay) {
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  // cos theta check
  if(SMbfound && (SMmufound || SMelfound) && (SMnuelfound || SMnumufound)){
    TLorentzVector tempLep;
    //cout << "SMmu " << SMmu << " SMel " << SMel << " SMnumu " << SMnumu << " SMnuel " << SMnuel << " SMb " << SMb << endl;
    
    if(SMmufound ){   tempLep.SetPxPyPzE(mcpartTLV[SMmu].Px(), mcpartTLV[SMmu].Py(), mcpartTLV[SMmu].Pz(), mcpartTLV[SMmu].Energy());}
    if(SMelfound ){  tempLep.SetPxPyPzE(mcpartTLV[SMel].Px(), mcpartTLV[SMel].Py(), mcpartTLV[SMel].Pz(), mcpartTLV[SMel].Energy());}
    TLorentzVector tempB;
    tempB.SetPxPyPzE(mcpartTLV[SMb].Px(), mcpartTLV[SMb].Py(), mcpartTLV[SMb].Pz(), mcpartTLV[SMb].Energy());
        TLorentzVector tempNu;
    if(SMnumufound ){ tempNu.SetPxPyPzE(mcpartTLV[SMnumu].Px(), mcpartTLV[SMnumu].Py(), mcpartTLV[SMnumu].Pz(), mcpartTLV[SMnumu].Energy());}
    if(SMnuelfound ){tempNu.SetPxPyPzE(mcpartTLV[SMnuel].Px(), mcpartTLV[SMnuel].Py(), mcpartTLV[SMnuel].Pz(), mcpartTLV[SMnuel].Energy());}
    
    //cout << "templep " << tempLep.Pt() << " x " << tempLep.Px() << " y " << tempLep.Py() << " z " << tempLep.Pz() << endl;
    //cout << " tempB " << tempB.Pt() << " x " << tempB.Px() <<" y " << tempB.Py() << " z " << tempB.Pz() << endl;
    //cout <<  " tempNu " << tempNu.Pt() << " x " << tempNu.Px() << " y " << tempNu.Py() << " z " << tempNu.Pz() << " E " << tempNu.Energy() << " Phi " << tempNu.Phi() << " Eta " << tempNu.Eta()<< endl;
    pair <float, float> tempcosPair = CosThetaCalculation(tempLep, tempNu, tempB, true);
    
  }
  
  
  
  //SM TOP
  
  if( SMbfound && SMWfound){
    histo1D["Topmass_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).M());
    histo1D["pt_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Pt());
    histo1D["eta_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Eta());
    histo1D["phi_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Phi());
    histo1D["dPhi_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb],mcpartTLV[SMW]));
    histo1D["dR_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb],mcpartTLV[SMW]));
    
  }
  if( SMbATop && SMWATop && antitopfound){
    //cout << "in tbar" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[antitopQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[antitopQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[antitopQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[antitopQ].Phi());
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Eta() );
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    
  }
  if( SMbTop && SMWTop && topfound){
    //cout << "in top" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[topQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[topQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[topQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[topQ].Phi());
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Eta() );
    
  }
  if( SMbfound && SMnumufound && SMmufound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    
    
  }
  if( SMbfound && SMnuelfound && SMelfound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    
    
  }
  
  
  //FCNC TOP
  if( FCNCmufound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[FCNCmuMin].Pt(), mcpartTLV[FCNCmuPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[FCNCmuMin].Phi(), mcpartTLV[FCNCmuPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[FCNCmuMin].Eta(), mcpartTLV[FCNCmuPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCmuMin],mcpartTLV[FCNCmuPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCmuMin],mcpartTLV[FCNCmuPlus]));
    histo1D["mass_lep1"]->Fill(mcpartTLV[FCNCmuMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[FCNCmuPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[FCNCmuMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[FCNCmuMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[FCNCmuMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[FCNCmuPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[FCNCmuMin].M(), mcpartTLV[FCNCmuPlus].M());
    
  }
  if( FCNCelfound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[FCNCelMin].Pt(), mcpartTLV[FCNCelPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[FCNCelMin].Phi(), mcpartTLV[FCNCelPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[FCNCelMin].Eta(), mcpartTLV[FCNCelPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[FCNCelMin],mcpartTLV[FCNCelPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[FCNCelMin],mcpartTLV[FCNCelPlus]));
    
    histo1D["mass_lep1"]->Fill(mcpartTLV[FCNCelMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[FCNCelPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[FCNCelMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[FCNCelPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[FCNCelMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[FCNCelPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[FCNCelMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[FCNCelPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[FCNCelMin].M(), mcpartTLV[FCNCelPlus].M());
    
  }
  
  if( FCNCZfound){
    histo1D["Zmass_Zbos"]->Fill(mcpartTLV[FCNCZ].M());
    histo1D["pt_Zbos"]->Fill(mcpartTLV[FCNCZ].Pt());
    histo1D["phi_Zbos"]->Fill(mcpartTLV[FCNCZ].Phi());
    histo1D["eta_Zbos"]->Fill(mcpartTLV[FCNCZ].Eta());
    histo1D["mc_nZLep"]->Fill(nbZDaughters);
    histo1D["mc_nZMu"]->Fill(nbZDaughters);
  }
  if( FCNCZfound && FCNCmufound){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).M(),mcpartTLV[FCNCZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[FCNCmuMin] + mcpartTLV[FCNCmuPlus]).Pt(),mcpartTLV[FCNCZ].Eta() );
  }
  
  if( FCNCelfound && FCNCZfound ){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).M(),mcpartTLV[FCNCZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[FCNCelMin] + mcpartTLV[FCNCelPlus]).Pt(),mcpartTLV[FCNCZ].Eta() );
    
  }
  
  if( SMWfound && SMelfound ){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWEl"]->Fill(nbWDaughters);
    
  }
  if( SMWfound && SMmufound){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWMu"]->Fill(nbWDaughters);
  }
  if( SMWfound && FCNCZfound){
    histo2D["nZbosonnWboson"]->Fill(nbZDaughters, nbWDaughters);
  }
  
  // end plots
  // begin matching
  
  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
   /* if( abs(mcParticles_[iMC]->type()) <= 5 || abs(mcParticles_[iMC]->type()) == 21){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    } // jets
    else */if( abs(mcParticles_[iMC]->type()) ==  13 || abs(mcParticles_[iMC]->type()) ==  11 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    } // leptons
  }
  
  
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedleptons,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedleptons collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 3){
    //cout << "ERROR" << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    /*if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMb");
    } // b from t */
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMmu");
    } // mu from W from t
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMel");
    } // el from W from t
    
    
    //FCNC
    if( mcParticles_[partonID[partonIDnb]]->type() ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmumin");
    } // mu- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelmin");
    } // el- from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCmuplus");
    } // mu+ from Z from t
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("FCNCelplus");
    } // el+ from Z from t
    
    
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
  /*  if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); } */
    if(NPair[iPart].find("SMmu")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("SMel")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    
    if(NPair[iPart].find("FCNCmumin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("FCNCelmin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("FCNCmuplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("FCNCelplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    
  }
  
  histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
  //histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  JetMatchertZq(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedjets){
  //cout << "in jet matcher for tZq" << endl;
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int radZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  int SMq = -999;
  int RadmuPlus = -999;
  int RadmuMin = -999;
  int RadelMin = -999;
  int RadelPlus = -999;
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  bool RadelMinFound = false;
  bool RadelPlusFound = false;
  bool RadmuMinFound = false;
  bool RadmuPlusFound = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  bool SMqfound = false;
  bool radZfound = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  //cout << "event " << evt_num_ << endl;
  
  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << iMC << " found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      //cout << iMC << " found b from t" << endl;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      //cout << iMC << " found mu- from tbar" << endl;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      //cout << iMC << " found m+ from t" << endl;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
      // cout << iMC << " found e- from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      //cout << iMC << " found e+ from t" << endl;
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    else if( abs(mcParticles_[iMC]->type()) >= 1 &&  abs(mcParticles_[iMC]->type())<  5){
      if(mcParticles_[iMC]->status() == 23) {SMq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMq== -999) {SMq = iMC;}
      //cout << iMC << " found q " << endl;
      SMqfound = true;
    } // q
    
    //radiated Z
    
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadmuMin== -999) {RadmuMin = iMC;}
      RadmuMinFound = true;
      
    } // mu - from Z
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadmuPlus== -999) {RadmuPlus= iMC;}
      RadmuPlusFound = true;
      
    } // mu + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadelMin== -999) {RadelMin = iMC;}
      RadelMinFound = true;
      
    } // e - from Z
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadelPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadelPlus== -999) {RadelPlus= iMC;}
      RadelPlusFound = true;
      
    } // e + from Z
    
    
    else if( mcParticles_[iMC]->type() ==  23 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {radZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && radZ== -999) {radZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      radZfound = true;
    } //Z  with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      if(mcParticles_[iMC]->status() == 23) {radZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && radZ== -999) {radZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      radZfound = true;
    } //Z  with el decay
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool Radmufound = false;
  if(RadmuPlusFound && RadmuMinFound){ Radmufound = true; }
  bool Radelfound = false;
  if(RadelPlusFound && RadelMinFound){ Radelfound = true; }
  //bool Radqfound = false;
  
  bool foundDecay = false;
  if((Radmufound || Radelfound) && (SMmufound || SMelfound) && SMbfound   ){
    foundDecay = true;
  }
  
  if(!foundDecay) {
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }

  //cout << "found decay " << endl;
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( abs(mcParticles_[iMC]->type()) <=  5 || abs(mcParticles_[iMC]->type()) ==  21 ){
      
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    } //jets and gluons
    
    
  }
  //cout << "selected jets " << selectedjets.size() <<endl;
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedjets,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedjets collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 1){
    //cout << "NOT FOUND IT" << endl;
    //cout << " JetPartonPair.size() " << JetPartonPair.size() << endl;
    //cout << " partons.size() " << partons.size() << endl;
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
     PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
     NPair.push_back("SMb");
     } // b from t
   
}
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempQ;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempQ.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
    if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedjets[PPair[iPart].first].Px(), selectedjets[PPair[iPart].first].Py(), selectedjets[PPair[iPart].first].Pz(), selectedjets[PPair[iPart].first].E());}
  }
 // histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
  //histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};
pair< vector< pair<unsigned int, unsigned int>>, vector <string> >  LeptonMatchertZq(vector<TRootMCParticle*> mcParticles_ , Long64_t evt_num_, vector<TLorentzVector> selectedleptons){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;
  
  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  int radZ = -999;
  int SMnuel = -999;
  int SMnumu = -999;
  int SMW = -999;
  int SMb = -999;
  int SMq = -999;
  int RadmuPlus = -999;
  int RadmuMin = -999;
  int RadelMin = -999;
  int RadelPlus = -999;
  
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool SMelTop = false;
  bool SMelATop = false;
  bool RadelMinFound = false;
  bool RadelPlusFound = false;
  bool RadmuMinFound = false;
  bool RadmuPlusFound = false;
  
  bool SMWATop = false;
  bool SMWTop = false;
  bool SMbATop = false;
  bool SMbTop = false;
  bool SMnuelfound = false;
  bool SMnumufound = false;
  bool topfound = false;
  bool antitopfound = false;
  bool SMqfound = false;
  bool radZfound = false;
  int nbZDaughters = 0;
  int nbWDaughters = 0;
  
  
  
  
  
  vector <TLorentzVector> partons;
  partons.clear();
  vector <int> partonID;
  partonID.clear();
  
  vector <int> storedMCParticles;
  storedMCParticles.clear();
  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "ERROR mcP not filled correctly" << endl;  }
  
  //cout << "event " << evt_num_ << endl;

  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; topfound = true;  }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; antitopfound = true; }
    
    //SM
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W- from tbar" << endl;
    } // nu el  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  12 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnuel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnuel == -999) {SMnuel = iMC;}
      SMnuelfound = true;
      //cout << "found nu from W+ from t" << endl;
    } // nu el  from W from t
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==-24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu  from W from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  14 && mcParticles_[iMC]->motherType() ==24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMnumu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMnumu == -999) {SMnumu = iMC;}
      SMnumufound = true;
    } // nu mu from W from t
    
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbATop = true;
      //cout << iMC << " found b from tbar" << endl;
      
    } // b from tbar
    else if( abs(mcParticles_[iMC]->type()) ==  5 && mcParticles_[iMC]->motherType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMb = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMb == -999) {SMb = iMC;}
      SMbTop = true;
      //cout << iMC << " found b from t" << endl;
      
    } // b from t
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuATop = true;
      //cout << iMC << " found mu- from tbar" << endl;
      
    } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMmu = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMmu == -999) {SMmu = iMC;}
      SMmuTop = true;
      //cout << iMC << " found m+ from t" << endl;
      
    } // mu+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == -24 && mcParticles_[iMC]->grannyType()  == -6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelATop = true;
     // cout << iMC << " found e- from tbar" << endl;
      
    } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 24 && mcParticles_[iMC]->grannyType()  == 6 ){
      if(mcParticles_[iMC]->status() == 23) {SMel = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMel== -999) {SMel = iMC;}
      SMelTop = true;
      //cout << iMC << " found e+ from t" << endl;
    } // el+ from W+ from top
    else if( mcParticles_[iMC]->type() ==  24 && mcParticles_[iMC]->motherType()  == 6   ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC; nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWTop = true;
      //cout << "W from t  found " << endl;
    } //W from top
    else if( mcParticles_[iMC]->type() ==  -24 && mcParticles_[iMC]->motherType()  == -6  ){
      if(mcParticles_[iMC]->status() == 23) {SMW= iMC; nbWDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && SMW== -999) {SMW = iMC;nbWDaughters = mcParticles_[iMC]->nDau();}
      SMWATop = true;
      //cout << "W from tbar found  " << endl;
    } //W from atop
    else if( abs(mcParticles_[iMC]->type()) >= 1 &&  abs(mcParticles_[iMC]->type())<  5){
      if(mcParticles_[iMC]->status() == 23) {SMq = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && SMq== -999) {SMq = iMC;}
      //cout << iMC << " found q " << endl;
      SMqfound = true;
    } // q
    
    //radiated Z
    
    else if( mcParticles_[iMC]->type() ==  13 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadmuMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadmuMin== -999) {RadmuMin = iMC;}
      RadmuMinFound = true;
      
    } // mu - from Z
    else if( mcParticles_[iMC]->type() ==  -13 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadmuPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadmuPlus== -999) {RadmuPlus= iMC;}
      RadmuPlusFound = true;
      
    } // mu + from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadelMin = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadelMin== -999) {RadelMin = iMC;}
      RadelMinFound = true;
      
    } // e - from Z
    else if( mcParticles_[iMC]->type() ==  -11 && mcParticles_[iMC]->motherType()  == 23  ){
      if(mcParticles_[iMC]->status() == 23) {RadelPlus = iMC;  }
      else if( mcParticles_[iMC]->status() != 23 && RadelPlus== -999) {RadelPlus= iMC;}
      RadelPlusFound = true;
      
    } // e + from Z
    
    
    else if( mcParticles_[iMC]->type() ==  23 && abs(mcParticles_[iMC]->dauOneId())  == 13  ){
      if(mcParticles_[iMC]->status() == 23) {radZ= iMC; nbZDaughters = mcParticles_[iMC]->nDau(); }
      else if( mcParticles_[iMC]->status() != 23 && radZ== -999) {radZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      radZfound = true;
    } //Z  with mu decay
    else if( mcParticles_[iMC]->type() ==  23 && mcParticles_[iMC]->motherType()  == 6 && abs(mcParticles_[iMC]->dauOneId())  == 11  ){
      if(mcParticles_[iMC]->status() == 23) {radZ= iMC;  nbZDaughters = mcParticles_[iMC]->nDau();}
      else if( mcParticles_[iMC]->status() != 23 && radZ== -999) {radZ = iMC; nbZDaughters = mcParticles_[iMC]->nDau();}
      radZfound = true;
    } //Z  with el decay
  }
  
  bool SMWfound = false;
  if(SMWATop || SMWTop){ SMWfound = true;}
  bool SMbfound = false;
  if(SMbATop || SMbTop){SMbfound = true;}
  bool SMmufound = false;
  if(SMmuTop || SMmuATop){ SMmufound = true;}
  bool SMelfound = false;
  if(SMelTop || SMelATop){ SMelfound = true; }
  bool Radmufound = false;
  if(RadmuPlusFound && RadmuMinFound){ Radmufound = true; }
  bool Radelfound = false;
  if(RadelPlusFound && RadelMinFound){ Radelfound = true; }
  //bool Radqfound = false;
  
  bool foundDecay = false;
  if((Radmufound || Radelfound) && (SMmufound || SMelfound) && SMbfound   ){
    foundDecay = true;
  }
  
  if(!foundDecay) {
    pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  // theta calcul
  // cos theta check
  if(SMbfound && (SMmufound || SMelfound) && (SMnuelfound || SMnumufound)){
    TLorentzVector tempLep;
    //cout << "SMmu " << SMmu << " SMel " << SMel << " SMnumu " << SMnumu << " SMnuel " << SMnuel << " SMb " << SMb << endl;
    
    if(SMmufound ){   tempLep.SetPxPyPzE(mcpartTLV[SMmu].Px(), mcpartTLV[SMmu].Py(), mcpartTLV[SMmu].Pz(), mcpartTLV[SMmu].Energy());}
    if(SMelfound ){ tempLep.SetPxPyPzE(mcpartTLV[SMel].Px(), mcpartTLV[SMel].Py(), mcpartTLV[SMel].Pz(), mcpartTLV[SMel].Energy());}
    TLorentzVector tempB;
    tempB.SetPxPyPzE(mcpartTLV[SMb].Px(), mcpartTLV[SMb].Py(), mcpartTLV[SMb].Pz(), mcpartTLV[SMb].Energy());
    
    TLorentzVector tempNu;
    if(SMnumufound ){ tempNu.SetPxPyPzE(mcpartTLV[SMnumu].Px(), mcpartTLV[SMnumu].Py(), mcpartTLV[SMnumu].Pz(), mcpartTLV[SMnumu].Energy());}
    if(SMnuelfound ){tempNu.SetPxPyPzE(mcpartTLV[SMnuel].Px(), mcpartTLV[SMnuel].Py(), mcpartTLV[SMnuel].Pz(), mcpartTLV[SMnuel].Energy());}
    
    //cout << "templep " << tempLep.Pt() << " x " << tempLep.Px() << " y " << tempLep.Py() << " z " << tempLep.Pz() << endl;
    //cout << " tempB " << tempB.Pt() << " x " << tempB.Px() <<" y " << tempB.Py() << " z " << tempB.Pz() << endl;
    //cout <<  " tempNu " << tempNu.Pt() << " x " << tempNu.Px() << " y " << tempNu.Py() << " z " << tempNu.Pz() << " E " << tempNu.Energy() << " Phi " << tempNu.Phi() << " Eta " << tempNu.Eta()<< endl;
     pair <float, float> tempcosPair = CosThetaCalculation(tempLep, tempNu, tempB, true);
    
  }
  
  
  
  
  //SM TOP
  
  if( SMbfound && SMWfound){
    histo1D["Topmass_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).M());
    histo1D["pt_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Pt());
    histo1D["eta_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Eta());
    histo1D["phi_Wb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMW]).Phi());
    histo1D["dPhi_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb],mcpartTLV[SMW]));
    histo1D["dR_Wb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb],mcpartTLV[SMW]));
    
  }
  if( SMbATop && SMWATop && antitopfound){
    //cout << "in tbar" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[antitopQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[antitopQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[antitopQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[antitopQ].Phi());
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[antitopQ].Eta() );
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[antitopQ]));
    
  }
  if( SMbTop && SMWTop && topfound){
    //cout << "in top" << endl;
    histo1D["Topmass_tq"]->Fill(mcpartTLV[topQ].M());
    histo1D["pt_tq"]->Fill(mcpartTLV[topQ].Pt());
    histo1D["eta_tq"]->Fill(mcpartTLV[topQ].Eta());
    histo1D["phi_tq"]->Fill(mcpartTLV[topQ].Phi());
    
    histo1D["dPhi_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    histo1D["dR_Wbtop"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb]+mcpartTLV[SMW],mcpartTLV[topQ]));
    
    
    histo2D["Topmass_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].M() );
    histo2D["pt_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Pt() );
    histo2D["phi_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Phi() );
    histo2D["eta_top_Wb"]->Fill((mcpartTLV[SMW] + mcpartTLV[SMb]).M(),mcpartTLV[topQ].Eta() );
    
  }
  if( SMbfound && SMnumufound && SMmufound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMmu]+ mcpartTLV[SMnumu]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMmu]+ mcpartTLV[SMnumu]));
    
    
  }
  if( SMbfound && SMnuelfound && SMelfound){
    histo1D["Topmass_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).M());
    histo1D["pt_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Pt());
    histo1D["eta_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Eta());
    histo1D["phi_lvb"]->Fill((mcpartTLV[SMb] + mcpartTLV[SMel]+ mcpartTLV[SMnuel]).Phi());
    
    histo1D["dPhi_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    histo1D["dR_lvb"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[SMb] , mcpartTLV[SMel]+ mcpartTLV[SMnuel]));
    
    
  }
  
  if( Radmufound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[RadmuMin].Pt(), mcpartTLV[RadmuPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[RadmuMin].Phi(), mcpartTLV[RadmuPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[RadmuMin].Eta(), mcpartTLV[RadmuPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[RadmuMin],mcpartTLV[RadmuPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[RadmuMin],mcpartTLV[RadmuPlus]));
    histo1D["mass_lep1"]->Fill(mcpartTLV[RadmuMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[RadmuPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[RadmuMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[RadmuPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[RadmuMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[RadmuPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[RadmuMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[RadmuPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[RadmuMin].M(), mcpartTLV[RadmuPlus].M());
    
  }
  if( Radelfound){
    histo1D["Zmass_Zlep"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).M());
    histo1D["pt_Zlep"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Pt());
    histo1D["eta_Zlep"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Eta());
    histo1D["phi_Zlep"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Phi());
    histo2D["pt_lep"]->Fill(mcpartTLV[RadelMin].Pt(), mcpartTLV[RadelPlus].Pt());
    histo2D["phi_lep"]->Fill(mcpartTLV[RadelMin].Phi(), mcpartTLV[RadelPlus].Phi());
    histo2D["eta_lep"]->Fill(mcpartTLV[RadelMin].Eta(), mcpartTLV[RadelPlus].Eta());
    histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcpartTLV[RadelMin],mcpartTLV[RadelPlus]));
    histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcpartTLV[RadelMin],mcpartTLV[RadelPlus]));
    
    histo1D["mass_lep1"]->Fill(mcpartTLV[RadelMin].M());
    histo1D["mass_lep2"]->Fill(mcpartTLV[RadelPlus].M());
    histo1D["pt_lep1"]->Fill(mcpartTLV[RadelMin].Pt());
    histo1D["pt_lep2"]->Fill(mcpartTLV[RadelPlus].Pt());
    histo1D["eta_lep1"]->Fill(mcpartTLV[RadelMin].Eta());
    histo1D["eta_lep2"]->Fill(mcpartTLV[RadelPlus].Eta());
    histo1D["phi_lep1"]->Fill(mcpartTLV[RadelMin].Phi());
    histo1D["phi_lep2"]->Fill(mcpartTLV[RadelPlus].Phi());
    histo2D["mass_lep"]->Fill(mcpartTLV[RadelMin].M(), mcpartTLV[RadelPlus].M());
    
  }
  
  if( radZfound){
    histo1D["Zmass_Zbos"]->Fill(mcpartTLV[radZ].M());
    histo1D["pt_Zbos"]->Fill(mcpartTLV[radZ].Pt());
    histo1D["phi_Zbos"]->Fill(mcpartTLV[radZ].Phi());
    histo1D["eta_Zbos"]->Fill(mcpartTLV[radZ].Eta());
    histo1D["mc_nZLep"]->Fill(nbZDaughters);
    histo1D["mc_nZMu"]->Fill(nbZDaughters);
  }
  if( radZfound && Radmufound){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).M(),mcpartTLV[radZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Pt(),mcpartTLV[radZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Pt(),mcpartTLV[radZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[RadmuMin] + mcpartTLV[RadmuPlus]).Pt(),mcpartTLV[radZ].Eta() );
  }
  
  if( Radelfound && radZfound ){
    histo2D["Zmass_Zbos_Zlep"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).M(),mcpartTLV[radZ].M() );
    histo2D["pt_Z"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Pt(),mcpartTLV[radZ].Pt() );
    histo2D["phi_Z"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Pt(),mcpartTLV[radZ].Phi() );
    histo2D["eta_Z"]->Fill((mcpartTLV[RadelMin] + mcpartTLV[RadelPlus]).Pt(),mcpartTLV[radZ].Eta() );
    
  }
  
  if( SMWfound && SMelfound ){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWEl"]->Fill(nbWDaughters);
    
  }
  if( SMWfound && SMmufound){
    histo1D["mc_nWLep"]->Fill(nbWDaughters);
    histo1D["mc_nWMu"]->Fill(nbWDaughters);
  }
  if( SMWfound && radZfound){
    histo2D["nZbosonnWboson"]->Fill(nbZDaughters, nbWDaughters);
  }
  
  
  
  vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
  PPair.clear();
  vector<string > NPair; // First one is jet number, second one is mcParticle number
  NPair.clear();
  
  
  for (unsigned int iMC = 0; iMC < mcpartTLV.size(); iMC++)
  {
    if (false)
      cout << setw(3) << right << iMC << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    /*if( abs(mcParticles_[iMC]->type()) <=  5 || abs(mcParticles_[iMC]->type()) ==  21 ){
     
     partons.push_back(mcpartTLV[iMC]);
     partonID.push_back(iMC);
     } //jets
    else*/
    if( abs(mcParticles_[iMC]->type() ) ==  13 || abs(mcParticles_[iMC]->type() ) ==  11 ){
      partons.push_back(mcpartTLV[iMC]);
      partonID.push_back(iMC);
    } // leptons
    
  }
  
  JetPartonMatching matchingTool = JetPartonMatching(partons, selectedleptons,2,true,true,0.3 );
  
  if (matchingTool.getNumberOfAvailableCombinations() != 1)
    cerr << "matching.getNumberOfAvailableCombinations() = " << matchingTool.getNumberOfAvailableCombinations() << " .  This should be equal to 1 !!!" << endl;
  
  
  /// Fill match in JetPartonPair;
  vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
  JetPartonPair.clear();
  
  
  for (unsigned int i = 0; i < partons.size(); i++)
  {
    int matchedJetNumber = matchingTool.getMatchForParton(i, 0);
    if (matchedJetNumber > -1)
      JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
    // matched jet number is nb in selectedleptons collection
    // i is nb in partons
    // partonID contains place in mcParticles_ vector
  }
  
  
  
  if(JetPartonPair.size() < 3){
    ///cout << "NOT FOUND IT" << endl;
    //cout << " JetPartonPair.size() " << JetPartonPair.size() << endl;
    //cout << " partons.size() " << partons.size() << endl;
   pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVectorEr;
    return returnVectorEr;
  }
  
  for (unsigned int i = 0; i < JetPartonPair.size(); i++)
  {
    unsigned int partonIDnb = JetPartonPair[i].second; // place in mcParticles_ vector
    unsigned int particlenb = JetPartonPair[i].first;  // place in selectedLeptons vector
    
    //SM
    
   /* if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  5 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMb");
    } // b from t*/
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMmu");
    } // mu from W from t
    if( abs(mcParticles_[partonID[partonIDnb]]->type()) ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 24 && abs(mcParticles_[partonID[partonIDnb]]->grannyType())  == 6 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMel");
    } // el from W from t
   /* if( abs(mcParticles_[partonID[partonIDnb]]->type()) >=  1 && abs(mcParticles_[partonID[partonIDnb]]->type()) < 5 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("SMq");
    } // q*/
    if( mcParticles_[partonID[partonIDnb]]->type() ==  13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("Radmumin");
    } // mu- from Z
    if( mcParticles_[partonID[partonIDnb]]->type() ==  11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("Radelmin");
    } // el- from Z
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -13 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23 ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("Radmuplus");
    } // mu+ from Z
    if( mcParticles_[partonID[partonIDnb]]->type() ==  -11 && abs(mcParticles_[partonID[partonIDnb]]->motherType())  == 23  ){
      PPair.push_back(pair<unsigned int,unsigned int> (JetPartonPair[i].first,JetPartonPair[i].second));
      NPair.push_back("Radelplus");
    } // el+ from Z
    
    
  }
  
  //cout << "PPair.size() " << PPair.size() << endl;
  //cout << "NPair.size() "  << NPair.size() << endl;
  
  
  TLorentzVector tempB;
  TLorentzVector tempQ;
  TLorentzVector tempWlep;
  TLorentzVector tempZlepm;
  TLorentzVector tempZlepp;
  tempB.Clear();
  tempQ.Clear();
  tempZlepm.Clear();
  tempZlepp.Clear();
  tempWlep.Clear();
  
  
  for(unsigned int iPart = 0 ; iPart < PPair.size(); iPart++){
    
    // cout << " iPart " << iPart << endl;
   /* if(NPair[iPart].find("SMb")!=string::npos){ tempB.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());}*/
    if(NPair[iPart].find("SMmu")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
    if(NPair[iPart].find("SMel")!=string::npos){ tempWlep.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());}
    
    if(NPair[iPart].find("Radmumin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());   }
    if(NPair[iPart].find("Radelmin")!=string::npos){ tempZlepm.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());  }
    if(NPair[iPart].find("Radmuplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());   }
    if(NPair[iPart].find("Radelplus")!=string::npos){ tempZlepp.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E()); }
   /* if(NPair[iPart].find("SMq")!=string::npos){ tempQ.SetPxPyPzE(selectedleptons[PPair[iPart].first].Px(), selectedleptons[PPair[iPart].first].Py(), selectedleptons[PPair[iPart].first].Pz(), selectedleptons[PPair[iPart].first].E());  }*/
  }
  histo1D["matchedZmass"]->Fill((tempZlepm+tempZlepp).M());
  //histo1D["matchedSMTopmass"]->Fill((tempWlep+tempB).M());
  
  
  
  
  pair< vector< pair<unsigned int, unsigned int>>, vector <string> >   returnVector;
  returnVector = pair< vector< pair<unsigned int, unsigned int>>, vector <string> >( PPair , NPair ) ;
  ////  cout << " (returnVector.second).size() " << (returnVector.second).size() << endl;
  //  cout << " (returnVector.first).size() " << (returnVector.first).size() << endl;
  //  cout << "Matcher PPair.size() " << PPair.size() << endl;
  //  cout << "Matcher NPair.size() "  << NPair.size() << endl;
  return returnVector;
};


std::pair <Float_t,Float_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo){
  // see https://github.com/TopBrussels/TopTreeAnalysis/blob/CMSSW_53X/WHelicities/src/BTagCosThetaCalculation.cc
  // ttbar http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2012/157
  // single top : https://cds.cern.ch/record/1601800
  
 
  float CosTheta = 999;
  
  //----------------------------------------------
  //  Calculating cos theta value
  //----------------------------------------------
  
  TRootMCParticle WLeptonic = (Neutrino+lepton);
  TRootMCParticle TopLeptonic = (Neutrino+lepton+leptonicBJet);
  TLorentzVector TopWRF = (Neutrino+lepton+leptonicBJet);
  TLorentzVector leptWRF = lepton;
  TLorentzVector WTRF;
  WTRF.SetPxPyPzE(WLeptonic.Px(), WLeptonic.Py(), WLeptonic.Pz(), WLeptonic.Energy());
 // TLorentzVector leptTRF = lepton;
  
  //Angle between Top in WRF and lepton in WRF
  TopWRF.Boost(-WLeptonic.BoostVector());
  leptWRF.Boost(-WLeptonic.BoostVector());
  
  // boost to top RF
  WTRF.Boost(-TopLeptonic.BoostVector());
  //leptTRF.Boost(-TopLeptonic.BoostVector());
  
  //Calculating cos:
  float ThetaTevatron = ROOT::Math::VectorUtil::Angle( TopWRF, leptWRF );
  CosTheta = -(TMath::Cos(ThetaTevatron));
  //Cos theta is defined as the angle between the lepton and the reversed direction of the top quark, both boosted to the W-boson rest frame.
  //Still reversed direction doesn't need to be calculated since angle between lepton and top and between lepton and reversed top is proportional to theta and Pi-theta.
  //For these the angles the following relation holds: cos(theta) = - cos(Pi-theta)
  // --> Hence the need of the minus sign in the CosTheta definition!!
  
  float ThetaSTAR = ROOT::Math::VectorUtil::Angle( WTRF, leptWRF );
  float CosThetaSTAR= TMath::Cos(ThetaSTAR);
  
  
  
  if(WLeptonic.E() < 0.){
    cout << " Event with negative WLept energy!!! (BTagCosThetaCalculation class) Cos theta = " << CosTheta << endl;
  }
  
  //cout << "cos " << CosThetaSTAR << " cos " << CosTheta <<  endl;
  /*
  
  As the W-boson helicity fractions are very sensitive to the Wtb vertex, their measurements can be used to investigate contribution from non-standard models.
  In this study, the t ̄t signal events are reconstructed where both W-bosons decay leptonically (electron or muon).
  In this case, the relevant angle θ∗ is defined in top quark rest frame, as the angle between the 3-momentum of the charged lepton in the rest
  frame of the decaying W-boson and the momentum of the W-boson in the rest frame of the decaying top quark.
  */
  
  if(!geninfo){
    histo1D["CosThetaWRF"]->Fill(CosTheta);
    histo1D["CosThetaWRFTRF"]->Fill(CosThetaSTAR);
    histo2D["CosTheta"]->Fill(CosTheta, CosThetaSTAR);
  }
  else if(geninfo){
    //cout << "filling gen info" << endl;
    histo1D["CosThetaWRF_gen"]->Fill(CosTheta);
    histo1D["CosThetaWRFTRF_gen"]->Fill(CosThetaSTAR);
    histo2D["CosTheta_gen"]->Fill(CosTheta, CosThetaSTAR);
  }
  
  std::pair <Float_t,Float_t> CosThetaPair(CosTheta, CosThetaSTAR);
  
  return CosThetaPair;
}



