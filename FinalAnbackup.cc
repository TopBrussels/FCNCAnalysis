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
Double_t workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations.
Double_t workingpointvalue_Medium = 0.800;//working points updated to 2016 BTV-POG recommendations.
Double_t workingpointvalue_Tight = 0.935;//working points updated to 2016 BTV-POG recommendations.

std::pair <Double_t,Double_t> c_workingpointvalue_Loose(-0.48, -0.17); // reduces b -jets (cvsln cvsb)
std::pair < Double_t, Double_t > c_workingpointvalue_Medium(-0.1, -0.08); // reduces light and b jets
std::pair <Double_t, Double_t> c_workingpointvalue_Tight(0.69, -0.45); // reduces light




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
vector<TLorentzVector> LeptonAssignerv2(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons);

bool isVetoElectronSpring2016(TRootElectron electron);
bool isTightElectronSpring2016(TRootElectron electron);
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
int Matcher(vector<TRootMCParticle*> mcParticles_)


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
  string xmlFileName = "";
  int maxMCParticles = -1;
  
  // to put on with agrs
  bool applyJER = false;
  bool applyJES = false;
  bool fillBtagHisto = false;
  
  
  
  
  //////////////////////////////////////////////
  /// Set up everything for local submission ////
  ///////////////////////////////////////////////
  // check the arguments passed
  if(true)
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
  const int JES                 =  strtol(argv[argc-6], NULL,10);
  const int JER                 =  strtol(argv[argc-5], NULL,10);
  const int FillBtagHisto	 =  strtol(argv[argc-4], NULL,10);
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
  for(int args = 11; args < argc-6; args++)
  {
    cout << "pushing back " << argv[args] << endl;
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
  /// define channels
  //
  
  cout << " --> Using the all  channel <-- " << endl;
  xmlFileName = "config/Run2TriLepton.xml" ;
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
  if( dName.find("NP_overlay_TT_FCNC")!=string::npos || dName.find("tZq")!=string::npos )
  {
    matching = true;
    cout << " looking at mcParticles !! " << endl;
  }
  if(dName.find("tZq")!=string::npos )
  {
    istZq = true;
   
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
  string histo_dir = "NtupleMakerOutput/TriLepton_histos";
  string histo_dir_date = histo_dir+"/TriLepton_histos_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
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
  
  
  histo1D["mass_lep1"]                                  = new TH1F("mass_lep1","mass_lep1",200,0,1);
  histo1D["mass_lep2"]                                  = new TH1F("mass_lep2","mass_lep2",200,0,200);
  histo1D["Zmass_Zlep"]                                  = new TH1F("Zmass_Zlep","Zmass_lep",200,0,200);
  histo1D["Zmass_Zbos"]                                  = new TH1F("Zmass_Zbos","Zmass_Zbos",200,0,200);
  histo2D["mass_lep"]            = new TH2F("mass_lep", "mass lep;mass lep1;mass lep2", 200,0,200,200,0,200 );
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
  histo1D["dR_elec"]          = new TH1F("DRelec", "dR electrons", 500,0, 5);
  histo1D["dR_mu"]          = new TH1F("DRmu", "dR muons", 500,0, 5);
  
  histo1D["dPhi_lep"]          = new TH1F("DPhilep", "dPhi lep", 140,-7, 7);
  histo1D["dPhi_elec"]          = new TH1F("DPhielec", "dPhi electrons", 140,-7, 7);
  histo1D["dPhi_mu"]          = new TH1F("DPhimu", "dPhi muons", 140,-7, 7);
  
  histo1D["mass_elec1"]       = new TH1F("mass_elec1","mass_elec1",200,0,200);
  histo1D["mass_elec2"]       = new TH1F("mass_elec2","mass_elec2",200,0,200);
  histo1D["Zmass_Zelec"]      = new TH1F("Zmass_Zelec","Zmass_elec",200,0,200);
  histo1D["Zmass_Zbos_elec"]  = new TH1F("Zmass_Zbos_elec","Zmass_Zbos_elec",200,0,200);
  histo2D["mass_elec"]        = new TH2F("mass_elec", "mass elec;mass elec1;mass elec2", 200,0,200,200,0,200 );
  histo2D["Zmass_Zbos_Zelec"] = new TH2F("Zmass_Zbos_Zelec", "Zmass;Zmass_elec;Zmass_Zboson", 200,0,200,200,0,200 );
  
  histo1D["pt_elec1"]          = new TH1F("ptelec1", "Pt elec 1", 200,0,400);
  histo1D["pt_elec2"]          = new TH1F("ptelec2", "Pt elec 2", 200,0,400);
  histo1D["pt_Zbos_elec"]          = new TH1F("ptZbos_elec", "Pt Zboson", 200,0,400);
  histo1D["pt_Zelec"]          = new TH1F("ptZelec", "Pt Zelec", 200,0,400);
  histo2D["pt_elec"]          = new TH2F("ptelec", "Pt electons;Pt elec 1;Pt elec 2", 200,0,400, 200,0,400);
  histo2D["pt_Zbos_Zelec"]          = new TH2F("ptZbos_Zelec", "Pt Z;Pt Zelec;Pt Zbos_elec", 200,0,400, 200,0,400);
  
  histo1D["phi_elec1"]          = new TH1F("phielec1", "phi elec 1",32,-3.2,3.2);
  histo1D["phi_elec2"]          = new TH1F("phielec2", "phi elec 2",32,-3.2,3.2);
  histo1D["phi_Zbos_elec"]          = new TH1F("phiZbos_elec", "phi Zboson",32,-3.2,3.2);
  histo1D["phi_Zelec"]          = new TH1F("phiZelec", "phi Zelec",32,-3.2,3.2);
  histo2D["phi_elec"]          = new TH2F("phielec", "phi elec;phi elec 1;phi elec 2",32,-3.2,3.2,32,-3.2,3.2);
  histo2D["phi_Zbos_Zelec"]          = new TH2F("phiZ_Zbos_Zelec", "phi Z;phi Zelec;phi Zbos_elec",32,-3.2,3.2,32,-3.2,3.2);
  
  histo1D["eta_elec1"]          = new TH1F("etaelec1", "eta elec 1", 30,-3,3);
  histo1D["eta_elec2"]          = new TH1F("etaelec2", "eta elec 2", 30,-3,3);
  histo1D["eta_Zbos_elec"]          = new TH1F("etaZbos_elec", "eta Zboson", 30,-3,3);
  histo1D["eta_Zelec"]          = new TH1F("etaZelec", "eta Zelec", 30,-3,3);
  histo2D["eta_elec"]          = new TH2F("etaelec", "eta elec;eta elec 1;eta elec 2", 30,-3,3, 30,-3,3);
  histo2D["eta_Zbos_Zelec"]          = new TH2F("etaZ_Zbos_Zelec", "eta Z;eta Zelec;eta Zbos", 30,-3,3, 30,-3,3);
  
  histo1D["mass_mu1"]                                  = new TH1F("mass_mu1","mass_mu1",200,0,200);
  histo1D["mass_mu2"]                                  = new TH1F("mass_mu2","mass_mu2",200,0,200);
  histo1D["Zmass_Zmu"]                                  = new TH1F("Zmass_Zmu","Zmass_mu",200,0,200);
  histo1D["Zmass_Zbos_mu"]                                  = new TH1F("Zmass_Zbos_mu","Zmass_Zbos_mu",200,0,200);
  histo2D["mass_mu"]            = new TH2F("mass_mu", "mass mu;mass mu1;mass mu2", 200,0,200,200,0,200 );
  histo2D["Zmass_Zbos_Zmu"]            = new TH2F("Zmass_Zbos_Zmu", "Zmass;Zmass_mu;Zmass_Zboson", 200,0,200,200,0,200 );
  
  histo1D["pt_mu1"]          = new TH1F("ptmu1", "Pt mu 1", 200,0,400);
  histo1D["pt_mu2"]          = new TH1F("ptmu2", "Pt mu 2", 200,0,400);
  histo1D["pt_Zbos_mu"]          = new TH1F("ptZbos_mu", "Pt Zboson", 200,0,400);
  histo1D["pt_Zmu"]          = new TH1F("ptZmu", "Pt Zmu", 200,0,400);
  histo2D["pt_mu"]          = new TH2F("ptmu", "Pt mutons;Pt mu 1;Pt mu 2", 200,0,400, 200,0,400);
  histo2D["pt_Zbos_Zmu"]          = new TH2F("ptZbos_Zmu", "Pt Z;Pt Zmu;Pt Zbos_mu", 200,0,400, 200,0,400);
  
  histo1D["phi_mu1"]          = new TH1F("phimu1", "phi mu 1",32,-3.2,3.2);
  histo1D["phi_mu2"]          = new TH1F("phimu2", "phi mu 2",32,-3.2,3.2);
  histo1D["phi_Zbos_mu"]          = new TH1F("phiZbos_mu", "phi Zboson",32,-3.2,3.2);
  histo1D["phi_Zmu"]          = new TH1F("phiZmu", "phi Zmu",32,-3.2,3.2);
  histo2D["phi_mu"]          = new TH2F("phimu", "phi mu;phi mu 1;phi mu 2",32,-3.2,3.2,32,-3.2,3.2);
  histo2D["phi_Zbos_Zmu"]          = new TH2F("phiZbos_Zmu", "phi Z;phi Zmu;phi Zbos_mu",32,-3.2,3.2,32,-3.2,3.2);
  
  histo1D["eta_mu1"]          = new TH1F("etamu1", "eta mu 1", 30,-3,3);
  histo1D["eta_mu2"]          = new TH1F("etamu2", "eta mu 2", 30,-3,3);
  histo1D["eta_Zbos_mu"]          = new TH1F("etaZbos_mu", "eta Zboson", 30,-3,3);
  histo1D["eta_Zmu"]          = new TH1F("etaZmu", "eta Zmu", 30,-3,3);
  histo2D["eta_mu"]          = new TH2F("etamu", "eta mu;eta mu 1;eta mu 2", 30,-3,3, 30,-3,3);
  histo2D["eta_Zbos_Zmu"]          = new TH2F("etaZbos_Zmu", "eta Z;eta Zmu;eta Zbos", 30,-3,3, 30,-3,3);
  
  histo1D["nRecoLeptons"] = new TH1F("nRecoLeptons","nRecoLeptons", 10,-0.5,9.5);
  histo1D["nIniRecoLeptons"] = new TH1F("nIniRecoLeptons","nIniRecoLeptons", 10,-0.5,9.5);
  histo1D["nRecoElectrons"] = new TH1F("nRecoElectrons","nRecoElectrons", 10,-0.5,9.5);
  histo1D["nIniRecoElectrons"] = new TH1F("nIniRecoElectrons","nIniRecoElectrons", 10,-0.5,9.5);
  histo1D["nRecoMuons"] = new TH1F("nRecoMuons","nRecoMuons", 10,-0.5,9.5);
  histo1D["nIniRecoMuons"] = new TH1F("nIniRecoMuons","nIniRecoMuons", 10,-0.5,9.5);
  
  histo1D["mc_nZ"] = new TH1F("mc_nZ","mc_nZ", 10,-0.5,9.5);
  histo1D["mc_nZEl"]= new TH1F("mc_nZEl","mc_nZEl", 10,-0.5,9.5);
  histo1D["mc_nZMu"]= new TH1F("mc_nZMu","mc_nZMu", 10,-0.5,9.5);
  histo1D["mc_nZLep"]= new TH1F("mc_nZLep","mc_nZLep", 10,-0.5,9.5);
  histo1D["mc_nWEl"]= new TH1F("mc_nWEl","mc_nWEl", 10,-0.5,9.5);
  histo1D["mc_nWMu"]= new TH1F("mc_nWMu","mc_nWMu", 10,-0.5,9.5);
  histo1D["mc_nWLep"]= new TH1F("mc_nWLep","mc_nWLep", 10,-0.5,9.5);
  histo1D["mc_nW"]= new TH1F("mc_nW","mc_nW", 10,-0.5,9.5);
  histo1D["mc_nTW"]= new TH1F("mc_nTW","mc_nTW", 10,-0.5,9.5);
  histo1D["mc_nTZ"]= new TH1F("mc_nTZ","mc_nTZ", 10,-0.5,9.5);
  histo1D["mc_nTWelectrons"]= new TH1F("mc_nTWelectrons","mc_nTW elec", 10,-0.5,9.5);
  histo1D["mc_nTZelectrons"]= new TH1F("mc_nTZelectrins","mc_nTZ elec", 10,-0.5,9.5);
  histo1D["mc_nTWmuons"]= new TH1F("mc_nTWmuons","mc_nTW mu", 10,-0.5,9.5);
  histo1D["mc_nTZmuons"]= new TH1F("mc_nTZmuons","mc_nTZ mu", 10,-0.5,9.5);
  
  histo2D["nWboson"] =  new TH2F("nWboson", "Wboson;nW;nLep", 10,-0.5,9.5, 10,-0.5,9.5);
  histo2D["nZboson"] =  new TH2F("nZboson", "Zboson;nZ;nLep", 10,-0.5,9.5, 10,-0.5,9.5);
  histo2D["nZbosonnWboson"] =  new TH2F("nZbosonnWboson", "Boson;nZ;nW", 10,-0.5,9.5, 10,-0.5,9.5);
  
  /////////////////////////////////
  /// Matching
  ///////////////////////////////
  //vector<TLorentzVector> mcParticlesTLV_charm, mcParticlesTLV_bottom, selectedJetsTLV, selectedMuonsTLV, selectedElectronsTLV, mcParticlesTLV_electrons, mcParticlesTLV_muons, mcParticlesTLV_Welectrons, mcParticlesTLV_Wmuons, mcParticlesTLV_Zboson, mcParticlesTLV_Wboson, mcParticlesTLV_Wtauelectrons, mcParticlesTLV_Wtaumuons, mcParticlesTLV_Ztauelectrons, mcParticlesTLV_Ztaumuons, mcParticlesTLV_Ztaus, mcParticlesTLV_Wtaus;
  // TLorentzVector cQuark, anticQuark;
  
  LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring16MiniAODv2-Asympt.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting//pileup_2016Data80X_Run273158-276811Cert.root", "pileup", "pileup");
  
  
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
      btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_80X_ichep_incl_ChangedTo_mujets.csv");
      btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_LOOSE, "mujets","central");
      if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
      {
        cout << "filling btag histo's" << endl;
        btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_mujets_central.root",false,30,999,2.4);
      }
      else
      {
        btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/Merged/"+daName+".root",false,30,999,2.4);
        //btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
      }
      
      
    }
    else if(!isData) // NEEDS TO BE CHECKED FOR 80X
    {
      BTagCalibration calib_csvv2("csvv2", "../TopTreeAnalysisBase/Calibrations/BTagging/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv");
      reader_csvv2 = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                               BTagEntry::OP_RESHAPING, // operating point
                                               "iterativefit", // measurement type
                                               "central"); // systematics type  --> depending on JES up/Down andother reader is needed
      
      
    }
    
    if(verbose>1) cout << "btag done" << endl;
    
    
    //MuonSFWeight(const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)
    
    MuonSFWeight* muonSFWeightID_T = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonID_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);
    // MuonSFWeight* muonSFWeightID_M = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonID_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true,  printLeptonSF, printLeptonSF);
    //MuonSFWeight* muonSFWeightID_L = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonID_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, printLeptonSF, printLeptonSF);
    MuonSFWeight* muonSFWeightIso_TT = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonIso_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);  // Tight RelIso, Tight ID
    //  MuonSFWeight* muonSFWeightIso_TM = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonIso_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Tight RelIso, Medium ID
    //   MuonSFWeight* muonSFWeightIso_LT = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonIso_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Tight ID
    // MuonSFWeight* muonSFWeightIso_LM = new MuonSFWeight(CaliPath+"LeptonSF/MuonSF/"+"MuonIso_Z_RunBCD_prompt80X_7p65.root", "MC_NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Medium ID
    
    
    if(verbose>1) cout << "muon SF loaded" << endl;
    
    
    ElectronSFWeight* electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/egammaEffi.txt_SF2D_CutBasedTightID.root","EGamma_SF2D",true,false,false);
    ElectronSFWeight* electronSFWeightReco = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/egammaEffi.txt_SF2D_GsfTrackingEff.root","EGamma_SF2D",true,false,false);
    
    if(verbose >1) cout << "electron SF loaded " << endl;
    
    vCorrParam.clear();
    JetCorrectionUncertainty *jecUnc;
    
    if(dName.find("Data_Run2016B")!=string::npos || dName.find("Data_Run2016C")!=string::npos || dName.find("Data_Run2016D")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10BCD_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10BCD_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10BCD_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10BCD_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10BCD_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data_Run2016E")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10E_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10E_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10E_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10E_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10E_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data_Run2016F")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10F_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10F_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10F_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10F_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10F_DATA_Uncertainty_AK4PFchs.txt");
    }
    else if(dName.find("Data")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
      isData = true;
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_Uncertainty_AK4PFchs.txt");
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10p2_DATA_Uncertainty_AK4PFchs.txt");
    }
    
    if(verbose>1) cout << "jec and jer loaded"<< endl;
    
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    if(verbose>1) cout << "jet tools initialised " << endl;
    
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
    
    
    if(matching )
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
      
      
      myTree->Branch("nRecoLeptons", &nRecoLeptons, "nRecoLeptons/I");
      myTree->Branch("nRecoElectrons", &nRecoElectrons, "nRecoElectrons/I");
      myTree->Branch("nRecoMuons", &nRecoMuons, "nRecoMuons/I");
      myTree->Branch("nIniRecoLeptons", &nIniRecoLeptons, "nIniRecoLeptons/I");
      myTree->Branch("nIniRecoElectrons", &nIniRecoElectrons, "nIniRecoElectrons/I");
      myTree->Branch("nIniRecoMuons", &nIniRecoMuons, "nIniRecoMuons/I");
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
      myTree->Branch("mc_Zmass", &mc_Zmass, "mc_Zmass/D");
      myTree->Branch("reco_Zmass", &reco_Zmass, "reco_Zmass/D");
      
      baselineTree->Branch("nRecoLeptons", &nRecoLeptons, "nRecoLeptons/I");
      baselineTree->Branch("nRecoElectrons", &nRecoElectrons, "nRecoElectrons/I");
      baselineTree->Branch("nRecoMuons", &nRecoMuons, "nRecoMuons/I");
      baselineTree->Branch("nIniRecoLeptons", &nIniRecoLeptons, "nIniRecoLeptons/I");
      baselineTree->Branch("nIniRecoElectrons", &nIniRecoElectrons, "nIniRecoElectrons/I");
      baselineTree->Branch("nIniRecoMuons", &nIniRecoMuons, "nIniRecoMuons/I");
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
      baselineTree->Branch("mc_Zmass", &mc_Zmass, "mc_Zmass/D");
      baselineTree->Branch("reco_Zmass", &reco_Zmass, "reco_Zmass/D");
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
    globalTree->Branch("nofEventsHLTv2",&nofEventsHLTv2,"nofEventsHLTv2/I");
    globalTree->Branch("nofEventsHLTv3",&nofEventsHLTv3,"nofEventsHLTv3/I");
    globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
    globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
    globalTree->Branch("nEv" , &nEv, "nEv/I");
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
    globalTree->Branch("WPc_CvsB_Loose", &WPc_CvsB_Loose, "WPc_CvsB_Loose/D");
    globalTree->Branch("WPc_CvsB_Medium", &WPc_CvsB_Medium, "WPc_CvsB_Medium/D");
    globalTree->Branch("WPc_CvsB_Tight", &WPc_CvsB_Tight, "WPc_CvsB_Tight/D");
    globalTree->Branch("WPc_CvsL_Loose", &WPc_CvsL_Loose, "WPc_CvsL_Loose/D");
    globalTree->Branch("WPc_CvsL_Medium", &WPc_CvsL_Medium, "WPc_CvsL_Medium/D");
    globalTree->Branch("WPc_CvsL_Tight", &WPc_CvsL_Tight, "WPc_CvsL_Tight/D");
    
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
    myTree->Branch("signalInt", &signalInt, "signalInt/I");
    baselineTree->Branch("signalInt", &signalInt, "signalInt/I");
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
    myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/F");
    myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/F");
    myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/F");
    myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/F");
    myTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/F");
    myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/F");
    myTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/F");
    myTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/F");
    myTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/F");
    myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/F");
    myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
    myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/F");
    myTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/F");
    myTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/F");
    myTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/F");
    myTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/F");
    myTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/F");
    myTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
    myTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
    myTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
    myTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
    myTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
    myTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/F");
    myTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/F");
    myTree->Branch("pt_electron_3",&pt_electron_3,"pt_electron_3/F");
    
    
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
    
    // muons
    myTree->Branch("nMuons",&nMuons, "nMuons/I");
    myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/F");
    myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/F");
    myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/F");
    myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/F");
    myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/F");
    myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/F");
    myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/F");
    myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/F");
    myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/F");
    myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/F");
    myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/F");
    myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
    myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
    myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/F");
    myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
    myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/F");
    myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/F");
    myTree->Branch("pt_muon_1",&pt_muon_1,"pt_muon_1/F");
    myTree->Branch("pt_muon_2",&pt_muon_2,"pt_muon_2/F");
    myTree->Branch("pt_muon_3",&pt_muon_3,"pt_muon_3/F");
    
    baselineTree->Branch("nMuons",&nMuons, "nMuons/I");
    baselineTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/F");
    baselineTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/F");
    baselineTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/F");
    baselineTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/F");
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
    
    // jets
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
    myTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
    myTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
    myTree->Branch("nJets_nonCSVL",&nJets_nonCSVL,"nJets_nonCSVL/I");
    myTree->Branch("nJets_nonCSVM",&nJets_nonCSVM,"nJets_nonCSVM/I");
    myTree->Branch("nJets_nonCSVT",&nJets_nonCSVT,"nJets_nonCSVT/I");
    myTree->Branch("nJets_nonCharmLCSVL",&nJets_nonCharmLCSVL,"nJets_nonCharmLCSVL/I");
    myTree->Branch("nJets_nonCharmLCSVM",&nJets_nonCharmLCSVM,"nJets_nonCharmLCSVM/I");
    myTree->Branch("nJets_nonCharmLCSVT",&nJets_nonCharmLCSVT,"nJets_nonCharmLCSVT/I");
    myTree->Branch("nJets_nonCharmMCSVL",&nJets_nonCharmMCSVL,"nJets_nonCharmMCSVL/I");
    myTree->Branch("nJets_nonCharmMCSVM",&nJets_nonCharmMCSVM,"nJets_nonCharmMCSVM/I");
    myTree->Branch("nJets_nonCharmMCSVT",&nJets_nonCharmMCSVT,"nJets_nonCharmMCSVT/I");
    myTree->Branch("nJets_nonCharmTCSVL",&nJets_nonCharmTCSVL,"nJets_nonCharmTCSVL/I");
    myTree->Branch("nJets_nonCharmTCSVM",&nJets_nonCharmTCSVM,"nJets_nonCharmTCSVM/I");
    myTree->Branch("nJets_nonCharmTCSVT",&nJets_nonCharmTCSVT,"nJets_nonCharmTCSVT/I");
    
    myTree->Branch("nJets_nonCharmLnonCSVL",&nJets_nonCharmLnonCSVL,"nJets_nonCharmLnonCSVL/I");
    myTree->Branch("nJets_nonCharmLnonCSVM",&nJets_nonCharmLnonCSVM,"nJets_nonCharmLnonCSVM/I");
    myTree->Branch("nJets_nonCharmLnonCSVT",&nJets_nonCharmLnonCSVT,"nJets_nonCharmLnonCSVT/I");
    myTree->Branch("nJets_nonCharmMnonCSVL",&nJets_nonCharmMnonCSVL,"nJets_nonCharmMnonCSVL/I");
    myTree->Branch("nJets_nonCharmMnonCSVM",&nJets_nonCharmMnonCSVM,"nJets_nonCharmMnonCSVM/I");
    myTree->Branch("nJets_nonCharmMnonCSVT",&nJets_nonCharmMnonCSVT,"nJets_nonCharmMnonCSVT/I");
    myTree->Branch("nJets_nonCharmTnonCSVL",&nJets_nonCharmTnonCSVL,"nJets_nonCharmTnonCSVL/I");
    myTree->Branch("nJets_nonCharmTnonCSVM",&nJets_nonCharmTnonCSVM,"nJets_nonCharmTnonCSVM/I");
    myTree->Branch("nJets_nonCharmTnonCSVT",&nJets_nonCharmTnonCSVT,"nJets_nonCharmTnonCSVT/I");
    
    myTree->Branch("nJets_CharmLnonCSVL",&nJets_CharmLnonCSVL,"nJets_CharmLnonCSVL/I");
    myTree->Branch("nJets_CharmLnonCSVM",&nJets_CharmLnonCSVM,"nJets_CharmLnonCSVM/I");
    myTree->Branch("nJets_CharmLnonCSVT",&nJets_CharmLnonCSVT,"nJets_CharmLnonCSVT/I");
    myTree->Branch("nJets_CharmMnonCSVL",&nJets_CharmMnonCSVL,"nJets_CharmMnonCSVL/I");
    myTree->Branch("nJets_CharmMnonCSVM",&nJets_CharmMnonCSVM,"nJets_CharmMnonCSVM/I");
    myTree->Branch("nJets_CharmMnonCSVT",&nJets_CharmMnonCSVT,"nJets_CharmMnonCSVT/I");
    myTree->Branch("nJets_CharmTnonCSVL",&nJets_CharmTnonCSVL,"nJets_CharmTnonCSVL/I");
    myTree->Branch("nJets_CharmTnonCSVM",&nJets_CharmTnonCSVM,"nJets_CharmTnonCSVM/I");
    myTree->Branch("nJets_CharmTnonCSVT",&nJets_CharmTnonCSVT,"nJets_CharmTnonCSVT/I");
    
    myTree->Branch("nJets_CharmLCSVL",&nJets_CharmLCSVL,"nJets_CharmLCSVL/I");
    myTree->Branch("nJets_CharmLCSVM",&nJets_CharmLCSVM,"nJets_CharmLCSVM/I");
    myTree->Branch("nJets_CharmLCSVT",&nJets_CharmLCSVT,"nJets_CharmLCSVT/I");
    myTree->Branch("nJets_CharmMCSVL",&nJets_CharmMCSVL,"nJets_CharmMCSVL/I");
    myTree->Branch("nJets_CharmMCSVM",&nJets_CharmMCSVM,"nJets_CharmMCSVM/I");
    myTree->Branch("nJets_CharmMCSVT",&nJets_CharmMCSVT,"nJets_CharmMCSVT/I");
    myTree->Branch("nJets_CharmTCSVL",&nJets_CharmTCSVL,"nJets_CharmTCSVL/I");
    myTree->Branch("nJets_CharmTCSVM",&nJets_CharmTCSVM,"nJets_CharmTCSVM/I");
    myTree->Branch("nJets_CharmTCSVT",&nJets_CharmTCSVT,"nJets_CharmTCSVT/I");
    
    
    myTree->Branch("nJets_CharmL",&nJets_CharmL,"nJets_CharmL/I");
    myTree->Branch("nJets_CharmM",&nJets_CharmM,"nJets_CharmM/I");
    myTree->Branch("nJets_CharmT",&nJets_CharmT,"nJets_CharmT/I");
    myTree->Branch("nJets_nonCharmL",&nJets_nonCharmL,"nJets_nonCharmL/I");
    myTree->Branch("nJets_nonCharmM",&nJets_nonCharmM,"nJets_nonCharmM/I");
    myTree->Branch("nJets_nonCharmT",&nJets_nonCharmT,"nJets_nonCharmT/I");
    
    
    myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/F");
    myTree->Branch("px_jet",px_jet,"px_jet[nJets]/F");
    myTree->Branch("py_jet",py_jet,"py_jet[nJets]/F");
    myTree->Branch("pz_jet",pz_jet,"pz_jet[nJets]/F");
    myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/F");
    myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/F");
    myTree->Branch("E_jet",E_jet,"E_jet[nJets]/F");
    myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/F");
    myTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/F");
    myTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/F");
    myTree->Branch("cdiscCvsL_jet_1",&cdiscCvsL_jet_1,"cdiscCvsL_jet_1/F");
    myTree->Branch("cdiscCvsB_jet_1",&cdiscCvsB_jet_1,"cdiscCvsB_jet_1/F");
    myTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/F");
    myTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/F");
    myTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/F");
    
    baselineTree->Branch("nJets",&nJets,"nJets/I");
    baselineTree->Branch("cdiscCvsL_jet_1",&cdiscCvsL_jet_1,"cdiscCvsL_jet_1/F");
    baselineTree->Branch("cdiscCvsB_jet_1",&cdiscCvsB_jet_1,"cdiscCvsB_jet_1/F");
    baselineTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
    baselineTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
    baselineTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
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
    
    
    baselineTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/F");
    baselineTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/F");
    baselineTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/F");
    baselineTree->Branch("E_jet",E_jet,"E_jet[nJets]/F");
    baselineTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    baselineTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/F");
    baselineTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/F");
    baselineTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/F");
    baselineTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/F");
    baselineTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/F");
    baselineTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/F");
    
    // Zboson
    myTree->Branch("Zboson_M",&Zboson_M,"Zboson_M/F");
    baselineTree->Branch("Zboson_M",&Zboson_M,"Zboson_M/F");
    myTree->Branch("Zboson2_M",&Zboson2_M,"Zboson2_M/F");
    baselineTree->Branch("Zboson2_M",&Zboson2_M,"Zboson2_M/F");
    myTree->Branch("mWt",&mWt,"mWt/F");
    baselineTree->Branch("mWt",&mWt,"mWt/F");
    myTree->Branch("FCNCtop_M",&FCNCtop_M,"FCNCtop_M/F");
    myTree->Branch("FCNCtop_tagger",&FCNCtop_M_tagger,"FCNCtop_M_tagger/F");
    baselineTree->Branch("FCNCtop_tagger",&FCNCtop_M_tagger,"FCNCtop_M_tagger/F");
    myTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/F");
    baselineTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/F");
    myTree->Branch("Zboson_Px",&Zboson_Px,"Zboson_Px/F");
    myTree->Branch("Zboson_Py",&Zboson_Py,"Zboson_Py/F");
    myTree->Branch("Zboson_Pz",&Zboson_Pz,"Zboson_Pz/F");
    myTree->Branch("Zboson_Energy",&Zboson_Energy,"Zboson_Energy/F");
    myTree->Branch("cjet_Pt",&cjet_Pt,"cjet_Pt/F");
    baselineTree->Branch("cjet_Pt",&cjet_Pt,"cjet_Pt/F");
    myTree->Branch("cjet_Pt_tagger",&cjet_Pt_tagger,"cjet_Pt_tagger/F");
    baselineTree->Branch("cjet_Pt_tagger",&cjet_Pt_tagger,"cjet_Pt_tagger/F");
    myTree->Branch("mlb",&mlb,"mlb/F");
    baselineTree->Branch("mlb",&mlb,"mlb/F");
    myTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/F");
    myTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/F");
    baselineTree->Branch("dRSMFCNCtop_tagger",&dRSMFCNCtop_tagger,"dRSMFCNCtop_tagger/F");
    myTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/F");
    myTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/F");
    myTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/F");
    baselineTree->Branch("dRWlepc_tagger",&dRWlepc_tagger,"dRWlepc_tagger/F");
    myTree->Branch("dRZb",&dRZb,"dRZb/F");
    myTree->Branch("dRZc",&dRZc,"dRZc/F");
    myTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/F");
    baselineTree->Branch("dRZc_tagger",&dRZc_tagger,"dRZc_tagger/F");
    myTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/F");
    myTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/F");
    baselineTree->Branch("dPhiSMFCNCtop_tagger",&dPhiSMFCNCtop_tagger,"dPhiSMFCNCtop_tagger/F");
    myTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/F");
    myTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/F");
    baselineTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/F");
    myTree->Branch("dPhiWlepc_tagger",&dPhiWlepc_tagger,"dPhiWlepc_tagger/F");
    myTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/F");
    myTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/F");
    myTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/F");
    baselineTree->Branch("dPhiZc_tagger",&dPhiZc_tagger,"dPhiZc_tagger/F");
    baselineTree->Branch("dRSMFCNCtop",&dRSMFCNCtop,"dRSMFCNCtop/F");
    baselineTree->Branch("dRWlepb",&dRWlepb,"dRWlepb/F");
    baselineTree->Branch("dRWlepc",&dRWlepc,"dRWlepc/F");
    baselineTree->Branch("dRZb",&dRZb,"dRZb/F");
    baselineTree->Branch("dRZc",&dRZc,"dRZc/F");
    baselineTree->Branch("dPhiSMFCNCtop",&dPhiSMFCNCtop,"dPhiSMFCNCtop/F");
    baselineTree->Branch("dPhiWlepb",&dPhiWlepb,"dPhiWlepb/F");
    baselineTree->Branch("dPhiWlepc",&dPhiWlepc,"dPhiWlepc/F");
    baselineTree->Branch("dPhiZb",&dPhiZb,"dPhiZb/F");
    baselineTree->Branch("dPhiZc",&dPhiZc,"dPhiZc/F");
    
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
    vector<TRootElectron*> selectedElectrons;
    vector<TRootElectron*> selectedVetoElectrons;
    vector<TRootPFJet*>    selectedJets;
    vector<TRootPFJet*>    PreselectedJets;
    vector<TRootMuon*>     selectedMuons;
    vector<TRootMuon*>     selectedLooseMuons;
    vector<TRootPFJet*>      selectedCSVLBJets;
    vector<TRootPFJet*>      selectedCSVMBJets;
    vector<TRootPFJet*>      selectedCSVTBJets;
    
    
    vector<TRootMCParticle*> mcParticles;
    vector <TRootPFJet*>     selectednonCSVLJets;
    vector<TRootPFJet*>      selectedCharmLJets;
    vector<TRootPFJet*>     selectedCharmMJets;
    vector<TRootPFJet*>     selectedCharmTJets;
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
    cutstep_string.clear();
    cutstep.clear();
    cutstep_eee.clear();
    cutstep_eeu.clear();
    cutstep_uue.clear();
    cutstep_uuu.clear();
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
      lep2 = false;
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
      mc_Zmass = -5.;
      reco_Zmass = -5.;
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
      PassedGoodPV = false;
      HBHEnoise = false;
      HBHEIso = false;
      CSCTight = false;
      badchan = false;
      badmu = false;
      EcalDead = false;
      //eeBad = false;
      eventweight = 1;
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
      METon = 0;
      orig_met_px = mets[0]->Px();
      orig_met_py = mets[0]->Py();
      orig_met_pt = sqrt(orig_met_px*orig_met_px + orig_met_py*orig_met_py);
      if((applyJES ||  applyJER) && isData)
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
      
      if(matching && !istZq) {Matcher(mcParticles)}
      
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
       
       if(unclusteredUp){
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
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first){   selectedCharmLJets.push_back(selectedJets[iJ]);   }
        else{   selectednonCharmLJets.push_back(selectedJets[iJ]);}
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first){   selectedCharmMJets.push_back(selectedJets[iJ]);   }
        else{   selectednonCharmMJets.push_back(selectedJets[iJ]);    }
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first){   selectedCharmTJets.push_back(selectedJets[iJ]);   }
        else{   selectednonCharmTJets.push_back(selectedJets[iJ]);    }
        
        
        /// c and b loose combi
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose ) {   selectedCLBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose){  selectednonCLBLJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first ){
          selectedCLnonBLJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCLnonBLJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose ) {   selectedCMBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose){  selectednonCMBLJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first ){
          selectedCMnonBLJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCMnonBLJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose ) {   selectedCTBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose){  selectednonCTBLJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first ){
          selectedCTnonBLJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCTnonBLJets.push_back(selectedJets[iJ]);}
        
        /// c and b medium combi
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium ) {   selectedCLBMJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium){  selectednonCLBMJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first ){
          selectedCLnonBMJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCLnonBMJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium ) {   selectedCMBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium){  selectednonCMBMJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first ){
          selectedCMnonBMJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCMnonBMJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium ) {   selectedCTBMJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium){  selectednonCTBMJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first ){
          selectedCTnonBMJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCTnonBMJets.push_back(selectedJets[iJ]);}
        
        /// c and b tight combi
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight ) {   selectedCLBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight){  selectednonCLBTJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Loose.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Loose.first ){
          selectedCLnonBTJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCLnonBTJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight ) {   selectedCMBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight){  selectednonCMBTJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Medium.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Medium.first ){
          selectedCMnonBTJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCMnonBTJets.push_back(selectedJets[iJ]);}
        
        if( selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first && selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight ) {   selectedCTBLJets.push_back(selectedJets[iJ]);  }
        else if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight){  selectednonCTBTJets.push_back(selectedJets[iJ]);   }
        else if(selectedJets[iJ]->ctag_pfCombinedCvsBJetTags() > c_workingpointvalue_Tight.second && selectedJets[iJ]->ctag_pfCombinedCvsLJetTags() > c_workingpointvalue_Tight.first ){
          selectedCTnonBTJets.push_back(selectedJets[iJ]);
        }
        else { selectednonCTnonBTJets.push_back(selectedJets[iJ]);}
        
        
        
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
        //PUweight = LumiWeights.ITweight((int)event->nTruePU());
        //not ready yet !
        PUweight = 1;
        
      }
      
      
      ///////////////////////////////
      //// Matching
      //////////////////////////////
      
      //cout << "matching " << matching << endl;
      
     /* if(matching){
        //cout << "in matching" << endl;
        int pdgID_charm = 4;
        int pdgID_bottom = 5;
        int pdgId_Z = 23;
        int pdgID_Z = 23;
        int pdgId_tau = 15;
        int pdgID_tau = 15;
        int pdgId_top = 6;
        int pdgId_W = 24;
        int pdgID_electron = 11;
        int pdgID_muon = 13;
        int pdgId_muon = 13;
        vector<TRootMCParticle*> mcParticlesMatching_;
        
        
        
        mcParticlesTLV_charm.clear();
        selectedJetsTLV.clear();
        mcParticlesTLV_bottom.clear();
        selectedElectronsTLV.clear();
        selectedMuonsTLV.clear();
        mcParticlesTLV_electrons.clear();
        mcParticlesTLV_Zboson.clear();
        mcParticlesTLV_Wboson.clear();
        mcParticlesTLV_Wtaus.clear();
        mcParticlesTLV_Wtaumuons.clear();
        mcParticlesTLV_Wtauelectrons.clear();
        mcParticlesTLV_Ztauelectrons.clear();
        mcParticlesTLV_Ztaumuons.clear();
        mcParticlesTLV_muons.clear();
        mcParticlesTLV_Welectrons.clear();
        mcParticlesTLV_Wmuons.clear();
        //cout << " ------ evt_num " << evt_num << " -------"<< endl;
        if(istZq){
          for (unsigned int i = 0; i < mcParticles.size(); i++)
          {
            
            //if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
            //  if(true)  cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
            
            
            if (mcParticles[i]->status()== 1 || ( mcParticles[i]->status() > 20 && mcParticles[i]->status() < 30 ) ){
              // cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
              
              //FCNC DECAY
              if ( fabs(mcParticles[i]->type()) == pdgID_charm && fabs(mcParticles[i]->motherType()) == pdgId_top){
                mcParticlesTLV_charm.push_back(*mcParticles[i]);
                //mcParticlesMatching_.push_back(mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              
              
              else if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_tau && fabs(mcParticles[i]->grannyType()) == pdgId_Z ){
                mcParticlesTLV_Ztaumuons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_tau && fabs(mcParticles[i]->grannyType()) == pdgId_Z ){
                mcParticlesTLV_Ztauelectrons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              
              else if ( fabs(mcParticles[i]->type()) == pdgID_bottom && fabs(mcParticles[i]->motherType()) == pdgId_top){
                mcParticlesTLV_bottom.push_back(*mcParticles[i]);
                //mcParticlesMatching_.push_back(mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_tau && fabs(mcParticles[i]->motherType()) == pdgId_Z && fabs(mcParticles[i]->grannyType()) == pdgId_top ){
                mcParticlesTLV_Ztaus.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_Z && fabs(mcParticles[i]->grannyType()) == pdgId_top){
                mcParticlesTLV_electrons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_Z && fabs(mcParticles[i]->grannyType()) == pdgId_top){
                mcParticlesTLV_muons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_Z && fabs(mcParticles[i]->motherType()) == pdgId_top  && fabs(mcParticles[i]->dauOneId()) == pdgId_muon  && fabs(mcParticles[i]->dauTwoId())== pdgId_muon ){
                mcParticlesTLV_Zboson.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_Z && fabs(mcParticles[i]->motherType()) == pdgId_top &&   fabs(mcParticles[i]->dauOneId()) == pdgID_electron  &&    fabs(mcParticles[i]->dauTwoId()) == pdgID_electron ){
                mcParticlesTLV_Zboson.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              
              
              // SM decay
              
              else if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_tau && fabs(mcParticles[i]->grannyType()) == pdgId_W ){
                mcParticlesTLV_Wtaumuons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_tau && fabs(mcParticles[i]->grannyType()) == pdgId_W ){
                mcParticlesTLV_Wtauelectrons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              
              else if(fabs(mcParticles[i]->type()) == pdgID_tau && fabs(mcParticles[i]->motherType()) == pdgId_W && fabs(mcParticles[i]->grannyType()) == pdgId_top ){
                mcParticlesTLV_Wtaus.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
             
              else if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_W && fabs(mcParticles[i]->grannyType()) == pdgId_top ){
                mcParticlesTLV_Welectrons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_W && fabs(mcParticles[i]->grannyType()) == pdgId_top){
                mcParticlesTLV_Wmuons.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              
              else if(fabs(mcParticles[i]->type()) == pdgId_W && fabs(mcParticles[i]->motherType()) == pdgId_top&&  (fabs(mcParticles[i]->dauOneId()) == pdgID_muon||  fabs(mcParticles[i]->dauTwoId()) == pdgID_muon) ){
                mcParticlesTLV_Wboson.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
              else if(fabs(mcParticles[i]->type()) == pdgId_W && fabs(mcParticles[i]->motherType()) == pdgId_top &&  ( fabs(mcParticles[i]->dauOneId()) == pdgID_electron ||  fabs(mcParticles[i]->dauTwoId()) == pdgID_electron) ){
                mcParticlesTLV_Wboson.push_back(*mcParticles[i]);
                histo1D["pdgID"]->Fill(mcParticles[i]->type());
                
              }
            }
          }
        }
        
    /*   if(istZq){
          if(fabs(mcParticles[i]->type()) == pdgID_Z && mcParticles[i]->nDau() == 2 &&  ( fabs(mcParticles[i]->dauOneId()) == pdgId_muon)  &&  ( fabs(mcParticles[i]->dauTwoId()) == pdgId_muon ) ){
            mcParticlesTLV_Zboson.push_back(*mcParticles[i]);
          }
          else if(fabs(mcParticles[i]->type()) == pdgID_Z && mcParticles[i]->nDau() == 2 &&  ( fabs(mcParticles[i]->dauOneId()) == pdgID_electron)  &&  (  fabs(mcParticles[i]->dauTwoId()) == pdgID_electron) ){
            mcParticlesTLV_Zboson.push_back(*mcParticles[i]);
          }
          if(fabs(mcParticles[i]->type()) == pdgId_W && mcParticles[i]->nDau() == 2 &&  ( fabs(mcParticles[i]->dauOneId()) == pdgID_muon||  fabs(mcParticles[i]->dauTwoId()) == pdgID_muon) ){
            mcParticlesTLV_Wboson.push_back(*mcParticles[i]);
          }
          else if(fabs(mcParticles[i]->type()) == pdgId_W && mcParticles[i]->nDau() == 2 &&  ( fabs(mcParticles[i]->dauOneId()) == pdgID_electron ||  fabs(mcParticles[i]->dauTwoId()) == pdgID_electron) ){
            mcParticlesTLV_Wboson.push_back(*mcParticles[i]);
          }

          if(fabs(mcParticles[i]->type()) == pdgID_electron && fabs(mcParticles[i]->motherType()) == pdgId_Z ){
            mcParticlesTLV_electrons.push_back(*mcParticles[i]);
            histo1D["pdgID"]->Fill(mcParticles[i]->type());
          }
          else if(fabs(mcParticles[i]->type()) == pdgID_muon && fabs(mcParticles[i]->motherType()) == pdgId_Z ){
            mcParticlesTLV_muons.push_back(*mcParticles[i]);
            histo1D["pdgID"]->Fill(mcParticles[i]->type());
          }
        }
*/
        
        //fill control histo's
        
        histo1D["mc_nZ"]->Fill(mcParticlesTLV_Zboson.size());
        histo1D["mc_nZLep"]->Fill(mcParticlesTLV_muons.size()+mcParticlesTLV_electrons.size());
        histo1D["mc_nZEl"]->Fill(mcParticlesTLV_electrons.size());
        histo1D["mc_nZMu"]->Fill(mcParticlesTLV_muons.size());
        histo1D["mc_nTZelectrons"]->Fill(mcParticlesTLV_Ztauelectrons.size());
        histo1D["mc_nTZmuons"]->Fill(mcParticlesTLV_Ztaumuons.size());
        histo1D["mc_nTWelectrons"]->Fill(mcParticlesTLV_Wtauelectrons.size());
        histo1D["mc_nTWmuons"]->Fill(mcParticlesTLV_Wtaumuons.size());
        histo1D["mc_nTZ"]->Fill(mcParticlesTLV_Ztaus.size());
        histo1D["mc_nTW"]->Fill(mcParticlesTLV_Wtaus.size());
        histo1D["mc_nWEl"]->Fill(mcParticlesTLV_Welectrons.size());
        histo1D["mc_nWMu"]->Fill(mcParticlesTLV_Wmuons.size());
        histo1D["mc_nWLep"]->Fill(mcParticlesTLV_Wmuons.size()+mcParticlesTLV_Welectrons.size());
        histo1D["mc_nW"]->Fill(mcParticlesTLV_Wboson.size());
        
        
        histo2D["nWboson"]->Fill(mcParticlesTLV_Wboson.size(), mcParticlesTLV_Wmuons.size()+mcParticlesTLV_Welectrons.size());
        histo2D["nZboson"]->Fill(mcParticlesTLV_Zboson.size(), mcParticlesTLV_muons.size()+mcParticlesTLV_electrons.size());
        histo2D["nZbosonnWboson"]->Fill(mcParticlesTLV_Zboson.size(), mcParticlesTLV_Wboson.size());
        
        if(mcParticlesTLV_Zboson.size()==1 && mcParticlesTLV_Wboson.size() == 1){
          if(mcParticlesTLV_Welectrons.size() == 1 || mcParticlesTLV_Wmuons.size() == 1){
            if(mcParticlesTLV_electrons.size() == 2 || mcParticlesTLV_muons.size()==2){
              cout << "FOUND AN EVENT" << endl;
            }
          }
        }
        
      //  if(mcParticlesTLV_Zboson.size()==1 && mcParticlesTLV_Wboson.size() == 1){
        if(mcParticlesTLV_Zboson.size()==1 ){
        //  if(mcParticlesTLV_Welectrons.size() == 1 || mcParticlesTLV_Wmuons.size() == 1){
            if(mcParticlesTLV_electrons.size()==2 ){
              
              /* for (unsigned int i = 0; i < mcParticles.size(); i++)
               {
               
               if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
               if(true)  cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: " << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
               
               }*/
              
              histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_electrons[0],mcParticlesTLV_electrons[1]));
              histo1D["dR_elec"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_electrons[0],mcParticlesTLV_electrons[1]));
              
              if(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_electrons[0],mcParticlesTLV_electrons[1]) < 0.1){
                for (unsigned int i = 0; i < mcParticles.size(); i++)
                {
                cout << setw(3) << right << "event " << evt_num << " " <<  i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: "  << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
                }
              }
              
              histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticlesTLV_electrons[0],mcParticlesTLV_electrons[1]));
              histo1D["dPhi_elec"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticlesTLV_electrons[0],mcParticlesTLV_electrons[1]));
              
              
              histo1D["mass_lep1"]->Fill(mcParticlesTLV_electrons[0].M());
              histo1D["mass_lep2"]->Fill(mcParticlesTLV_electrons[1].M());
              histo1D["Zmass_Zlep"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).M());
              histo1D["Zmass_Zbos"]->Fill(mcParticlesTLV_Zboson[0].M());
              histo2D["mass_lep"]->Fill(mcParticlesTLV_electrons[0].M(), mcParticlesTLV_electrons[1].M());
              histo2D["Zmass_Zbos_Zlep"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).M(), mcParticlesTLV_Zboson[0].M() );
              
              histo1D["pt_lep1"]->Fill(mcParticlesTLV_electrons[0].Pt());
              histo1D["pt_lep2"]->Fill(mcParticlesTLV_electrons[1].Pt());
              histo1D["pt_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Pt());
              histo1D["pt_Zlep"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Pt());
              histo2D["pt_lep"]->Fill(mcParticlesTLV_electrons[0].Pt(), mcParticlesTLV_electrons[1].Pt());
              histo2D["pt_Z"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Pt(), mcParticlesTLV_Zboson[0].Pt());
              
              histo1D["phi_lep1"]->Fill(mcParticlesTLV_electrons[0].Phi());
              histo1D["phi_lep2"]->Fill(mcParticlesTLV_electrons[1].Phi());
              histo1D["phi_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Phi());
              histo1D["phi_Zlep"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Phi());
              histo2D["phi_lep"]->Fill(mcParticlesTLV_electrons[0].Phi(), mcParticlesTLV_electrons[1].Phi());
              histo2D["phi_Z"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Phi(), mcParticlesTLV_Zboson[0].Phi());
              
              
              histo1D["eta_lep1"]->Fill(mcParticlesTLV_electrons[0].Eta());
              histo1D["eta_lep2"]->Fill(mcParticlesTLV_electrons[1].Eta());
              histo1D["eta_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Eta());
              histo1D["eta_Zlep"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Eta());
              histo2D["eta_lep"]->Fill(mcParticlesTLV_electrons[0].Eta(), mcParticlesTLV_electrons[1].Eta());
              histo2D["eta_Z"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Eta(), mcParticlesTLV_Zboson[0].Eta());
              
              histo1D["mass_elec1"]->Fill(mcParticlesTLV_electrons[0].M());
              histo1D["mass_elec2"]->Fill(mcParticlesTLV_electrons[1].M());
              histo1D["Zmass_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).M());
              histo1D["Zmass_Zbos_elec"]->Fill(mcParticlesTLV_Zboson[0].M());
              histo2D["mass_elec"]->Fill(mcParticlesTLV_electrons[0].M(), mcParticlesTLV_electrons[1].M());
              histo2D["Zmass_Zbos_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).M(), mcParticlesTLV_Zboson[0].M() );
              
              histo1D["pt_elec1"]->Fill(mcParticlesTLV_electrons[0].Pt());
              histo1D["pt_elec2"]->Fill(mcParticlesTLV_electrons[1].Pt());
              histo1D["pt_Zbos_elec"]->Fill(mcParticlesTLV_Zboson[0].Pt());
              histo1D["pt_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Pt());
              histo2D["pt_elec"]->Fill(mcParticlesTLV_electrons[0].Pt(), mcParticlesTLV_electrons[1].Pt());
              histo2D["pt_Zbos_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Pt(), mcParticlesTLV_Zboson[0].Pt());
              
              histo1D["phi_elec1"]->Fill(mcParticlesTLV_electrons[0].Phi());
              histo1D["phi_elec2"]->Fill(mcParticlesTLV_electrons[1].Phi());
              histo1D["phi_Zbos_elec"]->Fill(mcParticlesTLV_Zboson[0].Phi());
              histo1D["phi_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Phi());
              histo2D["phi_elec"]->Fill(mcParticlesTLV_electrons[0].Phi(), mcParticlesTLV_electrons[1].Phi());
              histo2D["phi_Zbos_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Phi(), mcParticlesTLV_Zboson[0].Phi());
              
              
              histo1D["eta_elec1"]->Fill(mcParticlesTLV_electrons[0].Eta());
              histo1D["eta_elec2"]->Fill(mcParticlesTLV_electrons[1].Eta());
              histo1D["eta_Zbos_elec"]->Fill(mcParticlesTLV_Zboson[0].Eta());
              histo1D["eta_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Eta());
              histo2D["eta_elec"]->Fill(mcParticlesTLV_electrons[0].Eta(), mcParticlesTLV_electrons[1].Eta());
              histo2D["eta_Zbos_Zelec"]->Fill((mcParticlesTLV_electrons[0] + mcParticlesTLV_electrons[1]).Eta(), mcParticlesTLV_Zboson[0].Eta());
              
            }
            else if(mcParticlesTLV_muons.size()==2){
              histo1D["dR_lep"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_muons[0],mcParticlesTLV_muons[1]));
              histo1D["dR_mu"]->Fill(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_muons[0],mcParticlesTLV_muons[1]));
              if(ROOT::Math::VectorUtil::DeltaR(mcParticlesTLV_muons[0],mcParticlesTLV_muons[1]) < 0.1){
                for (unsigned int i = 0; i < mcParticles.size(); i++)
                {
                cout << setw(3) << right << "event " << evt_num << " " <<  i << "  Status: " << setw(2) << mcParticles[i]->status() << "  pdgId: " << setw(3) << mcParticles[i]->type() << "  Mother: "  << setw(4) << mcParticles[i]->motherType() << "  Granny: " << setw(4) << mcParticles[i]->grannyType() << "  Pt: " << setw(7) << left << mcParticles[i]->Pt() << "  Eta: " << mcParticles[i]->Eta() << endl;
                }
                
              }
              
              
              histo1D["dPhi_lep"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticlesTLV_muons[0],mcParticlesTLV_muons[1]));
              histo1D["dPhi_mu"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(mcParticlesTLV_muons[0],mcParticlesTLV_muons[1]));
              
              histo1D["mass_lep1"]->Fill(mcParticlesTLV_muons[0].M());
              histo1D["mass_lep2"]->Fill(mcParticlesTLV_muons[1].M());
              histo1D["Zmass_Zlep"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).M());
              histo1D["Zmass_Zbos"]->Fill(mcParticlesTLV_Zboson[0].M());
              histo2D["mass_lep"]->Fill(mcParticlesTLV_muons[0].M(), mcParticlesTLV_muons[1].M());
              histo2D["Zmass_Zbos_Zlep"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).M(), mcParticlesTLV_Zboson[0].M() );
              
              histo1D["pt_lep1"]->Fill(mcParticlesTLV_muons[0].Pt());
              histo1D["pt_lep2"]->Fill(mcParticlesTLV_muons[1].Pt());
              histo1D["pt_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Pt());
              histo1D["pt_Zlep"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Pt());
              histo2D["pt_lep"]->Fill(mcParticlesTLV_muons[0].Pt(), mcParticlesTLV_muons[1].Pt());
              histo2D["pt_Z"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Pt(), mcParticlesTLV_Zboson[0].Pt());
              
              histo1D["phi_lep1"]->Fill(mcParticlesTLV_muons[0].Phi());
              histo1D["phi_lep2"]->Fill(mcParticlesTLV_muons[1].Phi());
              histo1D["phi_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Phi());
              histo1D["phi_Zlep"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Phi());
              histo2D["phi_lep"]->Fill(mcParticlesTLV_muons[0].Phi(), mcParticlesTLV_muons[1].Phi());
              histo2D["phi_Z"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Phi(), mcParticlesTLV_Zboson[0].Phi());
              
              
              histo1D["eta_lep1"]->Fill(mcParticlesTLV_muons[0].Eta());
              histo1D["eta_lep2"]->Fill(mcParticlesTLV_muons[1].Eta());
              histo1D["eta_Zbos"]->Fill(mcParticlesTLV_Zboson[0].Eta());
              histo1D["eta_Zlep"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Eta());
              histo2D["eta_lep"]->Fill(mcParticlesTLV_muons[0].Eta(), mcParticlesTLV_muons[1].Eta());
              histo2D["eta_Z"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Eta(), mcParticlesTLV_Zboson[0].Eta());
              
              histo1D["mass_mu1"]->Fill(mcParticlesTLV_muons[0].M());
              histo1D["mass_mu2"]->Fill(mcParticlesTLV_muons[1].M());
              histo1D["Zmass_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).M());
              histo1D["Zmass_Zbos_mu"]->Fill(mcParticlesTLV_Zboson[0].M());
              histo2D["mass_mu"]->Fill(mcParticlesTLV_muons[0].M(), mcParticlesTLV_muons[1].M());
              histo2D["Zmass_Zbos_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).M(), mcParticlesTLV_Zboson[0].M() );
              
              histo1D["pt_mu1"]->Fill(mcParticlesTLV_muons[0].Pt());
              histo1D["pt_mu2"]->Fill(mcParticlesTLV_muons[1].Pt());
              histo1D["pt_Zbos_mu"]->Fill(mcParticlesTLV_Zboson[0].Pt());
              histo1D["pt_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Pt());
              histo2D["pt_mu"]->Fill(mcParticlesTLV_muons[0].Pt(), mcParticlesTLV_muons[1].Pt());
              histo2D["pt_Zbos_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Pt(), mcParticlesTLV_Zboson[0].Pt());
              
              histo1D["phi_mu1"]->Fill(mcParticlesTLV_muons[0].Phi());
              histo1D["phi_mu2"]->Fill(mcParticlesTLV_muons[1].Phi());
              histo1D["phi_Zbos_mu"]->Fill(mcParticlesTLV_Zboson[0].Phi());
              histo1D["phi_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Phi());
              histo2D["phi_mu"]->Fill(mcParticlesTLV_muons[0].Phi(), mcParticlesTLV_muons[1].Phi());
              histo2D["phi_Zbos_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Phi(), mcParticlesTLV_Zboson[0].Phi());
              
              
              histo1D["eta_mu1"]->Fill(mcParticlesTLV_muons[0].Eta());
              histo1D["eta_mu2"]->Fill(mcParticlesTLV_muons[1].Eta());
              histo1D["eta_Zbos_mu"]->Fill(mcParticlesTLV_Zboson[0].Eta());
              histo1D["eta_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Eta());
              histo2D["eta_mu"]->Fill(mcParticlesTLV_muons[0].Eta(), mcParticlesTLV_muons[1].Eta());
              histo2D["eta_Zbos_Zmu"]->Fill((mcParticlesTLV_muons[0] + mcParticlesTLV_muons[1]).Eta(), mcParticlesTLV_Zboson[0].Eta());
            }
         // }
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
      
      // to be ok with triggers
      /*  if(trigged && (selectedElectrons.size()!=0 || selectedMuons.size() !=0)){
       cout << "----------------------------" << endl;
       cout << "check to be ok with trig" << endl;
       cout <<"Init Number of Muons, Electrons, Jets  ===>  " << init_muons.size() <<" "  << init_electrons.size()<<" "<< init_jets.size()   << endl;
       for(int iMu = 0 ; iMu < init_muons.size() ; iMu++){
       cout << "mu " << init_muons[iMu]->Pt() << " " << init_muons[iMu]->Eta()  << endl;
       }
       for(int iEl = 0 ; iEl < init_electrons.size() ; iEl++){
       cout << "el " << init_electrons[iEl]->Pt() << " " << init_electrons[iEl]->Eta() << endl;
       }
       cout <<"Number of Muons, Electrons, Jets  ===>  " << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< selectedJets.size()   << endl;
       for(int iMu = 0 ; iMu < selectedMuons.size() ; iMu++){
       cout << "mu " <<  selectedMuons[iMu]->Pt() << " " << selectedMuons[iMu]->Eta() << endl;
       }
       for(int iEl = 0 ; iEl < selectedElectrons.size() ; iEl++){
       cout << "el " << selectedElectrons[iEl]->Pt() << " " << selectedElectrons[iEl]->Eta() << endl;
       }
       cout << "Name, Trigged, Flow" << dName << " " << trigged << " " << continueFlow << endl;
       }*/
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
          baseSelected = true;
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
        vector <TLorentzVector> leptons = LeptonAssignerv2(selectedElectrons, selectedMuons);
        if(Assigned){
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
          Zboson2_Energy = ( lep0 + lep1).Energy();}
        
      }
      
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
          eventSelected = true;
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_3++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_3++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_3++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_3++; }
          histo1D["cutFlow"]->Fill(3., eventweight);
          // baseSelected = true;
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
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_4++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_4++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_4++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_4++; }
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
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_5++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_5++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_5++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_5++; }
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
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_6++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_6++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_6++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_6++; }
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
          //cout << "ncuts " << nCuts << endl;
          
          if(selectedMuons.size() == 3) {nbEvents_uuu_7++;}
          else if(selectedElectrons.size() == 3) {nbEvents_eee_7++;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_7++; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_7++; }
        }
      }
      
      //	   if(continueFlow)  eventSelected = true;
      //	   else eventSelected = false;
      if(passedMET && continueFlow){
        histo1D["cutFlow"]->Fill(8., eventweight);
        nCuts++;
        nbEvents_8++;
        
        if(selectedMuons.size() == 3) {nbEvents_uuu_8++;}
        else if(selectedElectrons.size() == 3) {nbEvents_eee_8++;}
        else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {nbEvents_eeu_8++; }
        else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){nbEvents_uue_8++; }
        
        
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
        // bjets
        nJets_CSVT =  selectedCSVTBJets.size();
        nJets_CSVM =  selectedCSVMBJets.size();
        nJets_CSVL =  selectedCSVLBJets.size();
        nJets_nonCSVL = selectednonCSVLJets.size();
        nJets_nonCSVM = selectednonCSVMJets.size();
        nJets_nonCSVT = selectednonCSVTJets.size();
        
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
          
          
          cout << " nMatched_charm " << nMatched_charm << " nNonMatched_charm" << nNonMatched_charm << endl;
          
          
        }
        if(matching && JetPartonPair_charm.size()>0){
          
          if(JetPartonPair_charm[0].first == cjetindex_tagger)  nMatched_charm_tag++;
          else nNonMatched_charm_tag++;
          
          cout << " nMatched_charm_tag " << nMatched_charm_tag << " nNonMatched_charm_tag " << nNonMatched_charm_tag << endl;
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
          
          if( ((elecIndices[0] == JetPartonPair_electron[0].first) && (elecIndices[1] == JetPartonPair_electron[1].first)) || ((elecIndices[1] == JetPartonPair_electron[0].first) && (elecIndices[0] == JetPartonPair_electron[1].first))){
            mc_Zmass= (mcParticlesTLV_electrons[JetPartonPair_electron[0].second] + mcParticlesTLV_electrons[JetPartonPair_electron[1].second]).M();
            reco_Zmass= (*selectedElectrons[elecIndices[0]]+*selectedElectrons[elecIndices[1]]).M();
          }
        }
        if(matching && muIndices.size() == 2){
          cout << " muIndices " << muIndices[0] << " " << muIndices[1] << endl; 
          if( (muIndices[0] == JetPartonPair_muon[0].first) && (muIndices[1] == JetPartonPair_muon[1].first)) nMatched_Zmu++;
          else if((muIndices[1] == JetPartonPair_muon[0].first) && (muIndices[0] == JetPartonPair_muon[1].first)) nMatched_Zmu++;
          else nNonMatched_Zmu++;
          if(((muIndices[0] == JetPartonPair_muon[0].first) && (muIndices[1] == JetPartonPair_muon[1].first) )|| ((muIndices[1] == JetPartonPair_muon[0].first) && (muIndices[0] == JetPartonPair_muon[1].first))){
            mc_Zmass= (mcParticlesTLV_muons[JetPartonPair_muon[0].second]+ mcParticlesTLV_muons[JetPartonPair_muon[1].second]).M();
            reco_Zmass= (*selectedMuons[muIndices[0]]+*selectedMuons[muIndices[1]]).M();
          }
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

vector <TLorentzVector> LeptonAssignerv2(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons)
{
  // cout << " in assigner " << endl;
  vector<TLorentzVector> ReturnColl;
  Assigned = false;
  
  
  
  elecbool = false;
  mubool = false;
  elecIndices.clear();
  muIndices.clear();
  
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
      elecbool = true;
      elecIndices.push_back(0);
      elecIndices.push_back(1);
      //cout << "assigned " <<endl;
    }
  }
  else if(muons.size() == 2){
    //    cout << "2 muons" << endl;
    if(muons[0]->charge() != muons[1]->charge()){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      
      Assigned = true;
      mubool = true;
      muIndices.push_back(0);
      muIndices.push_back(1);
      
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
      Assigned = true;
      elecIndices.push_back(0); elecIndices.push_back(1);
    }
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(electrons[0]->Px(), electrons[0]->Py(),electrons[0]->Pz(),electrons[0]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Assigned = true;
      elecIndices.push_back(0); elecIndices.push_back(2);
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(electrons[1]->Px(), electrons[1]->Py(),electrons[1]->Pz(),electrons[1]->Energy());
      Zlepcan1.SetPxPyPzE(electrons[2]->Px(), electrons[2]->Py(),electrons[2]->Pz(),electrons[2]->Energy());
      Assigned = true;
      elecIndices.push_back(1); elecIndices.push_back(2);
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
      Assigned = true;
      muIndices.push_back(0); muIndices.push_back(1);
    }
    else if(mass02 <= mass12 && mass02 < mass01){
      Zlepcan0.SetPxPyPzE(muons[0]->Px(), muons[0]->Py(),muons[0]->Pz(),muons[0]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Assigned = true;
      muIndices.push_back(0); muIndices.push_back(2);
    }
    else if(mass12 < mass01 && mass12 < mass02){
      Zlepcan0.SetPxPyPzE(muons[1]->Px(), muons[1]->Py(),muons[1]->Pz(),muons[1]->Energy());
      Zlepcan1.SetPxPyPzE(muons[2]->Px(), muons[2]->Py(),muons[2]->Pz(),muons[2]->Energy());
      Assigned = true;
      muIndices.push_back(1); muIndices.push_back(2);
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



int Matcher(vector<TRootMCParticle*> mcParticles_){
  int nMCP = mcParticles_.size();
  TLorentzVector mcpart;
  vector <TLorentzVector> mcpartTLV;

  
  int topQ = -999;
  int antitopQ = -999;
  int SMmu = -999;
  int SMel = -999;
  bool SMmuTop = false;
  bool SMmuATop = false;
  bool FCNCZATop = false;
  bool FCNCZTop = false;

  

  
  for (int iMC = 0; iMC < nMCP; iMC++)
  {
    mcpart.Clear();
    mcpart.SetPtEtaPhiE(mcParticles_[iMC]->Pt(), mcParticles_[iMC]->Eta(), mcParticles_[iMC]->Phi(), mcParticles_[iMC]->E());
    mcpartTLV.push_back(mcpart);
  }
  if(mcpartTLV.size() != mcParticles_.size()){cout << "mcP not filled correctly" << endl; return 0; }
  
  for (unsigned int i = 0; i < mcpartTLV.size(); i++)
  {
    if (verbose > 4)
      cout << setw(3) << right << i << "  Status: " << setw(2) << mcParticles_[iMC]->status() << "  pdgId: " << setw(3) << mcParticles_[iMC]->type() << "  Mother: " << setw(4) << mcParticles_[iMC]->motherType() << "  Granny: " << setw(4) << mcParticles_[iMC]->grannyType() << "  Pt: " << setw(7) << left << mcParticles_[iMC]->Pt() << "  Eta: " << mcParticles_[iMC]->Eta() << endl;
    
    
    if ( (mcParticles_[iMC]->status() > 1 && mcParticles_[iMC]->status() <= 20) || mcParticles_[iMC]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process
    
    if( mcParticles_[iMC]->type() == 6 ){ topQ = iMC; }
    else if( mcParticles_[iMC]->type() == -6 ){ antitopQ = iMC; }
    
    else if( mcParticles_[iMC]->type() ==  13 && mc_mother[i] == -24 && mc_granny[i] == -6 ){ SMmu = iMC;  SMmuATop = true; } // mu - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mc_mother[i] == 24 && mc_granny[i] == 6 ){ SMmu = iMC;  SMmuTop = true; } // mu+ from W+ from top
    
    else if( mcParticles_[iMC]->type() ==  11 && mc_mother[i] == -24 && mc_granny[i] == -6 ){ SMel = iMC;  SMelATop = true; } // el - from W - from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mc_mother[i] == 24 && mc_granny[i] == 6 ){ SMel = iMC;  SMelTop = true; } // el+ from W+ from top
    
    else if( mcParticles_[iMC]->type() ==  13 && mc_mother[i] == 23 && mc_granny[i] == -6 ){ FCNCmuMin = iMC;  FCNCZATop = true; } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mc_mother[i] == 23 && mc_granny[i] == -6 ){ FCNCmuPlus = iMC;  FCNCZATop = true; } // mu - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -13 && mc_mother[i] == 23 && mc_granny[i] == 6 ){ FCNCmuPlus = iMC;  FCNCZTop = true; } // mu+ from Z from top
    else if( mcParticles_[iMC]->type() ==  13 && mc_mother[i] == 23 && mc_granny[i] == 6 ){ FCNCmuMin = iMC;  FCNCZTop = true; } // mu+ from Z from top
    
    else if( mcParticles_[iMC]->type() ==  11 && mc_mother[i] == 23 && mc_granny[i] == -6 ){ FCNCmuMin = iMC;  FCNCZATop = true; } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mc_mother[i] == 23 && mc_granny[i] == -6 ){ FCNCmuPlus = iMC;  FCNCZATop = true; } // el - from Z  from tbar
    else if( mcParticles_[iMC]->type() ==  -11 && mc_mother[i] == 23 && mc_granny[i] == 6 ){ FCNCmuPlus = iMC;  FCNCZTop = true; } // el+ from Z from top
    else if( mcParticles_[iMC]->type() ==  11 && mc_mother[i] == 23 && mc_granny[i] == 6 ){ FCNCmuMin = iMC;  FCNCZTop = true; } // el+ from Z from top
    
    
  }
  if(SMmuTop && FCNCZTop) { cout << "bad event 1" << endl; }
  if(SMmuATop && FCNCZATop) { cout << "bad event 2" << endl; }
  if(SMelTop && FCNCZTop) { cout << "bad event 3" << endl; }
  if(SMelATop && FCNCZATop) { cout << "bad event 4" << endl; }
  if(SMmuTop && SMelTop ) { cout << "bad event 5 " << endl; }
  if(SMmuATop && SMelATop ) { cout << "bad event 6" << endl; }
  if(SMmuTop && FCNCZATop) { cout << "SM: t->W->mu   FCNCZ" << endl; }
  if(SMmuATop && FCNCZTop) { cout << "SM: tbar->W->mu   FCNCZ" << endl; }
  if(SMelTop && FCNCZATop) { cout << "SM: t->W->el   FCNCZ" << endl; }
  if(SMelATop && FCNCZTop) { cout << "SM: tbar->W->el   FCNCZ" << endl; }
  
  return 0;
}





