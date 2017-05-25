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

#include "Math/PtEtaPhiE4D.h"
#include "rochester/RoccoR.cc"
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


double lum_RunsBCDEF = 19.680129;//15.658183109;// /fb
double lum_RunsGH = 16.146177; // 15.199167277;// /fb



// home made functions
std::pair <Double_t,Double_t> CosThetaCalculation(TLorentzVector lepton, TLorentzVector Neutrino, TLorentzVector leptonicBJet, bool geninfo);
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
  cout << "nb of arguments " << argc << endl;
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
  Int_t nbSelectedEvents = 0;
  Int_t nbSelectedEvents_3L = 0;
  Int_t nbSelectedEvents_3Lveto = 0;
   Int_t nbSelectedEvents_3Lveto1J = 0;
  Int_t wnbSelectedEvents = 0;
  Int_t wnbSelectedEvents_3L = 0;
  Int_t wnbSelectedEvents_3Lveto = 0;
  Int_t wnbSelectedEvents_3Lveto1J = 0;
  int nbEvents = 0;
  double dataLumi = 0; //pb
  
  // to put  on here
  bool runHLT = true;
  bool applyJetLeptonCleaning = true;
  
  bool printTrigger = false;
  bool printLeptonSF = false;
  bool applyPU = true;
  
  int maxMCParticles = -1;
  
  // to put on with agrs
  bool applyJER = false;
  bool applyJES = false;
  
  bool doFakeLepton  =false;
  
  
  //////////////////////////////////////////////
  /// Set up everything for local submission ////
  ///////////////////////////////////////////////
  // check the arguments passed
  if(verbose>0)
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
  const int JES                 =  strtol(argv[argc-8], NULL,10);
  const int JER                 =  strtol(argv[argc-7], NULL,10);
  const int doFakes	 =  strtol(argv[argc-6], NULL,10);
  const int doJESJERshiftarg = strtol(argv[argc-5], NULL,10);
  const int doFakeshiftarg = strtol(argv[argc-4], NULL,10);
  const int JobNum		  = strtol(argv[argc-3], NULL, 10);
  const int startEvent  	  = strtol(argv[argc-2], NULL, 10);
  const int endEvent		  = strtol(argv[argc-1], NULL, 10);
  
  //cout << "argv[argc-1] " << argv[argc-1] << endl;
  
  applyJES = JES;
  applyJER = JER;
  doFakeLepton= doFakes;
  int doJESJERshift = doJESJERshiftarg;
  int doFakeshift = doFakeshiftarg;
  const int Usettbar = 1;
  
  // all the files are stored from arg 11 to argc-8
  vector<string> vecfileNames;
  
  for(int args = 11; args < argc-8; args++)
  {
   // cout << "pushing back " << argv[args] << endl;
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
  if(!isData)  anaEnv.METCollection = "PFMET_slimmedMETs";
  else if(isData) anaEnv.METCollection = "PFMET_slimmedMETs";
  //anaEnv.METCollection = "PFMET_slimmedMETs";
  anaEnv.MuonCollection = "Muons_slimmedMuons";
  anaEnv.ElectronCollection = "Electrons_selectedElectrons";
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
  Trigger* trigger_mumu  = new Trigger(1, 0, 0, 1,0,0,0); // mu , el, single, double, tri, met
  Trigger* trigger_ee  = new Trigger(0, 1, 0, 1,0,0,0);
  Trigger* trigger_emu  = new Trigger(1, 1, 0, 1,0,0,0) ;
  Trigger* trigger_mumumu  = new Trigger(1, 0, 0, 0,1,0,0);
  Trigger* trigger_eee  = new Trigger(0, 1, 0, 0,1,0,0);
  Trigger* trigger_emumu_mumue  = new Trigger(1, 1, 0, 0,1,0,0) ;
  Trigger* trigger_mu  = new Trigger(1, 0, 1, 0,0,0,0);
  Trigger* trigger_e  = new Trigger(0, 1, 1, 0,0,0,0);
  Trigger* trigger_met = new Trigger(0, 0, 0, 0,0,1,0);
  Trigger* trigger_jet = new Trigger(0, 0, 0, 0,0,0,1);
  
  ////////////////////////
  // intialize  Calibrations      //
  ///////////////////////

  BTagCalibrationReader * reader_csvv2;
  BTagCalibrationReader * reader_csvv2_JESdown;
  BTagCalibrationReader * reader_csvv2_JESup;
  BTagCalibrationReader * reader_csvv2_LFUp;
  BTagCalibrationReader * reader_csvv2_LFDown;
  BTagCalibrationReader * reader_csvv2_HFUp;
  BTagCalibrationReader * reader_csvv2_HFDown;
  BTagCalibrationReader * reader_csvv2_HFStats1Up;
  BTagCalibrationReader * reader_csvv2_HFStats1Down;
  BTagCalibrationReader * reader_csvv2_HFStats2Up;
  BTagCalibrationReader * reader_csvv2_HFStats2Down;
  BTagCalibrationReader * reader_csvv2_LFStats1Up;
  BTagCalibrationReader * reader_csvv2_LFStats1Down;
  BTagCalibrationReader * reader_csvv2_LFStats2Up;
  BTagCalibrationReader * reader_csvv2_LFStats2Down;
  BTagCalibrationReader * reader_csvv2_CFErr1Up;
  BTagCalibrationReader * reader_csvv2_CFErr1Down;
  BTagCalibrationReader * reader_csvv2_CFErr2Up;
  BTagCalibrationReader * reader_csvv2_CFErr2Down;
  
  BTagCalibration *btagcalib_c;
  BTagCalibrationReader *btagreader_c;
  BTagWeightTools *btwt_c = 0;
  BTagCalibrationReader * reader_charm;
  //TFile *histoFileHandle = 0;
  // for pu
  LumiReWeighting LumiWeights, LumiWeights_up, LumiWeights_down;
  
  // JER / JEC
  vector<JetCorrectorParameters> vCorrParam;
  
  
  ///////////////////////////////
  //  Set up Output ROOT file  ///
  //////////////////////////////
  stringstream ss;
  ss << JobNum;
  string strJobNum = ss.str();
 /* string histo_dir = "NtupleMakerOutput/TriLepton_histos/";
  string histo_dirdecay = histo_dir ;
  string histo_dir_date = histo_dirdecay+"/TriLepton_histos_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
  mkdir(histo_dirdecay.c_str(),0777);
  mkdir(histo_dir_date.c_str(),0777);*/
  
  /*string rootFileName (histo_dir_date+"/FCNC_3L_"+dName+".root");
  if (strJobNum != "0")
  {
    if(verbose > 0) cout << "strJobNum is " << strJobNum << endl;
    rootFileName = histo_dir_date+"/FCNC_3L_"+dName + "_"+strJobNum+".root";
  }
  cout << "Histofile: " << rootFileName << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  */
  ///////////////////////////
  /// Global variables ////
  //////////////////////////
  TRootEvent* event = 0;
  // TRootRun *runInfos = new TRootRun();
  
  /////////////////////////////
  /// Object ID              ///
  /////////////////////////////
  // electron
  float el_pt_cut =30.; // 42
  float el_eta_cut = 2.5;
  float el_iso_cone  = 0.3;
  // reliso cut fabs(eta supercluster) <= 1.479 --> 0.107587 // (fabs(eta supercluster) > 1.479 && fabs(eta supercluster) < 2.5) --> 0.113254
  // muon
  float mu_pt_cut = 30.; // 40
  float mu_eta_cut = 2.4;
  float mu_iso_cut = 0.15;
  float mu_iso_cut_loose = 0.25;
  float mu_iso_cut_fakeloose = 0.25;
  float mu_iso_cut_faketight= 0.15;
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
  LumiWeights_up = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");
  LumiWeights_down = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");
  
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
    double normfactor = datasets[d]->NormFactor();
    cout <<"found sample " << daName.c_str() << " with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    if(daName.find("Data")!=string::npos || daName.find("data")!=string::npos || daName.find("DATA")!=string::npos){
      isData = true;
    }
    if(daName.find("amc")!=string::npos) nlo = true;
    /////////////////////////////////////////
    ///    Calibrations                  ///
    ////////////////////////////////////////
    string histfile = "BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_comb_central.root";
    string histreadfile = "BTagHistosPtEta/BTagFile/Btag.root";
    if(!isData )
    {
      // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
        BTagCalibration calib_csvv2("csvv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
      reader_csvv2 = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                               BTagEntry::OP_RESHAPING, // operating point
                                               "iterativefit", // measurement type
                                               "central"); // systematics type  --> depending on JES up/Down andother reader is needed
      reader_csvv2_JESup = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                     BTagEntry::OP_RESHAPING, // operating point
                                                     "iterativefit", // measurement type
                                                     "up_jes"); // systematics type  --> depending on JES up/Down andother reader is needed
      reader_csvv2_JESdown = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                       BTagEntry::OP_RESHAPING, // operating point
                                                       "iterativefit", // measurement type
                                                       "down_jes"); // systematics type  --> depending on JES up/Down andother reader is needed
      
      
      // LFUp
      reader_csvv2_LFUp = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                    BTagEntry::OP_RESHAPING, // operating point
                                                    "iterativefit", // measurement type
                                                    "up_lf"); // systematics type
      // LFDown
      reader_csvv2_LFDown = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                      BTagEntry::OP_RESHAPING, // operating point
                                                      "iterativefit", // measurement type
                                                      "down_lf"); // systematics type
      
      // HFUp
      reader_csvv2_HFUp = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                    BTagEntry::OP_RESHAPING, // operating point
                                                    "iterativefit", // measurement type
                                                    "up_hf"); // systematics type
      // HFDown
      reader_csvv2_HFDown = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                      BTagEntry::OP_RESHAPING, // operating point
                                                      "iterativefit", // measurement type
                                                      "down_hf"); // systematics type
      
      // HFStats1Up
      reader_csvv2_HFStats1Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "up_hfstats1"); // systematics type
      // HFStats1Down
      reader_csvv2_HFStats1Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                            BTagEntry::OP_RESHAPING, // operating point
                                                            "iterativefit", // measurement type
                                                            "down_hfstats1"); // systematics type
      
      // HFStats2Up
      reader_csvv2_HFStats2Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "up_hfstats2"); // systematics type
      // HFStats2Down
      reader_csvv2_HFStats2Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                            BTagEntry::OP_RESHAPING, // operating point
                                                            "iterativefit", // measurement type
                                                            "down_hfstats2"); // systematics type
      
      // LFStats1Up
      reader_csvv2_LFStats1Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "up_lfstats1"); // systematics type
      // LFStats1Down
      reader_csvv2_LFStats1Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                            BTagEntry::OP_RESHAPING, // operating point
                                                            "iterativefit", // measurement type
                                                            "down_lfstats1"); // systematics type
      
      // LFStats2Up
      reader_csvv2_LFStats2Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "up_lfstats2"); // systematics type
      // LFStats2Down
      reader_csvv2_LFStats2Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                            BTagEntry::OP_RESHAPING, // operating point
                                                            "iterativefit", // measurement type
                                                            "down_lfstats2"); // systematics type
      
      // CFErr1Up
      reader_csvv2_CFErr1Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                        BTagEntry::OP_RESHAPING, // operating point
                                                        "iterativefit", // measurement type
                                                        "up_cferr1"); // systematics type
      // CFErr1Down
      reader_csvv2_CFErr1Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "down_cferr1"); // systematics type
      
      // CFErr2Up
      reader_csvv2_CFErr2Up = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                        BTagEntry::OP_RESHAPING, // operating point
                                                        "iterativefit", // measurement type
                                                        "up_cferr2"); // systematics type
      // CFErr2Down
      reader_csvv2_CFErr2Down = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                          BTagEntry::OP_RESHAPING, // operating point
                                                          "iterativefit", // measurement type
                                                          "down_cferr2"); // systematics type
      
      
    }
    
    if(verbose>1) cout << "btag done" << endl;
    
    TFile *muontrackfile = new TFile("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/Tracking_EfficienciesAndSF_BCDEFGH.root","read");
    TGraph* h_muonSFWeightTrack = (TGraph*) muontrackfile->Get("ratio_eff_eta3_dr030e030_corr")->Clone();//Tracking efficiency as function of eta

    
    
    MuonSFWeight* muonSFWeightID_BCDEF;
    MuonSFWeight* muonSFWeightID_GH;
    MuonSFWeight* muonSFWeightIso_BCDEF;
    MuonSFWeight* muonSFWeightIso_GH;
   //MuonSFWeight* muonSFWeightTrig_BCDEF;
   // MuonSFWeight* muonSFWeightTrig_GH;
    
    
    ElectronSFWeight* electronSFWeightID;
    ElectronSFWeight* electronSFWeightReco;
    
    muonSFWeightID_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/IDEfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
    muonSFWeightID_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF//20170413/IDEfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
    muonSFWeightIso_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/IsoEfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
    muonSFWeightIso_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/IsoEfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
    // muonSFWeightTrig_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsBCDEF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
    //  muonSFWeightTrig_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsGH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
    RoccoR rc("rochester/rcdata.2016.v3");
    
    
    
    if(verbose>1) cout << "muon SF loaded" << endl;
    
    electronSFWeightID = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/20170413/egammaEffi.txt_EGM2D_IDcutbTight_20170413.root","EGamma_SF2D",true,false,false, false); // extend range, debug,  print warnings
    electronSFWeightReco = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/20170413/egammaEffi.txt_EGM2D_reco_20170413.root","EGamma_SF2D",true,false,false, true);
    
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
    if(!isData) doFakeLepton = false;
    
    string channel_dir = "NtupleMakerOutput/Ntuples/" ;
    if(doFakeLepton) channel_dir = "NtupleMakerOutput/Ntuples_fakes/" ;
    string date_dir = channel_dir+ "/Ntuples_" + dateString +"/";
    mkdir("NtupleMakerOutput",0777);
    mkdir(channel_dir.c_str(),0777);
    mkdir(date_dir.c_str(),0777);
    
    string postfix = "";
    if(doJESJERshift == 1) postfix = "_JESdown" ;
    if(doJESJERshift == 2) postfix = "_JESup" ;
    if(doJESJERshift == 3) postfix = "_JERdown" ;
    if(doJESJERshift == 4) postfix = "_JERup" ;
    if(doFakeshift == 1) postfix ="_FakeShift";
    string Ntupname = date_dir +"FCNC_3L_" + dName + "_"+  strJobNum + postfix + ".root";
    cout << "Ntuple " << Ntupname << " created " << endl;
    
    TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
    tupfile->cd();
    cout << "outputfile " << Ntupname.c_str()<< " made and accessed " <<endl ;
    
    TTree* myTree = new TTree(("tree"+postfix).c_str(),("tree"+postfix).c_str());
    TTree* globalTree = new TTree(("globaltree"+postfix).c_str(),("globaltree"+postfix).c_str());
    ///////////////////////////
    /// output tree
    ///////////////////////////
    
    // event related variables
    Int_t run_num;
    double i_channel;
    Double_t x1;
    Double_t x2;
    Int_t id1;
    Int_t id2;
    Double_t q;
    Long64_t evt_num;
    Int_t lumi_num;
    Int_t nvtx;
    Int_t npu;
    Int_t PassedMETFilter;
    
    Int_t PassedTrigger;
    Int_t PassedTriggerMET;
    Int_t PassedTriggerJET;
    Int_t PassedTriggerNoLogic;
    Int_t PassedTriggerNoLogic2;
    Int_t PassedGoodPV;
    Int_t withfakes;
    Double_t hdamp_up;
    Double_t hdamp_down;
    Double_t puSF;
    Double_t puSF_up;
    Double_t puSF_down;

    Double_t btagSFshape = 1.;
    Double_t btagSFshape_down_cferr1 = 1.;
    Double_t btagSFshape_down_cferr2 = 1.;
    Double_t btagSFshape_down_hf= 1.;
    Double_t btagSFshape_down_hfstats1 = 1.;
    Double_t btagSFshape_down_hfstats2 = 1.;
    Double_t btagSFshape_down_lf = 1.;
    Double_t btagSFshape_down_lfstats1 = 1.;
    Double_t btagSFshape_down_lfstats2 = 1.;
    
    Double_t btagSFshape_up_cferr1 = 1.;
    Double_t btagSFshape_up_cferr2 = 1.;
    Double_t btagSFshape_up_hf= 1.;
    Double_t btagSFshape_up_hfstats1 = 1.;
    Double_t btagSFshape_up_hfstats2 = 1.;
    Double_t btagSFshape_up_lf = 1.;
    Double_t btagSFshape_up_lfstats1 = 1.;
    Double_t btagSFshape_up_lfstats2 = 1.;
    
    
    Double_t MuonIDSF[10];
    Double_t MuonIsoSF[10];
    Double_t MuonIDSF_BCDEF[10];
    Double_t MuonIsoSF_BCDEF[10];
    Double_t MuonIDSF_GH[10];
    Double_t MuonIsoSF_GH[10];
    Double_t MuonIDSF_BCDEF_up[10];
    Double_t MuonIsoSF_BCDEF_up[10];
    Double_t MuonIDSF_GH_up[10];
    Double_t MuonIsoSF_GH_up[10];
    Double_t MuonIDSF_BCDEF_down[10];
    Double_t MuonIsoSF_BCDEF_down[10];
    Double_t MuonIDSF_GH_down[10];
    Double_t MuonIsoSF_GH_down[10];
    
    
    Double_t MuonTrackSF[10];
    Double_t MuonTrackSF_up[10];
    Double_t MuonTrackSF_down[10];
    Double_t MuonIDSF_up[10];
    Double_t MuonIsoSF_up[10];
    Double_t MuonIDSF_down[10];
    Double_t MuonIsoSF_down[10];
    Double_t MuonTrigSFv2[10];
    Double_t MuonTrigSFv3[10];
    Double_t ElectronSF[10];
    Double_t ElectronSF_up[10];
    Double_t ElectronSF_down[10];
    Double_t ElectronSFID[10];
    Double_t ElectronSFID_up[10];
    Double_t ElectronSFID_down[10];
    Double_t ElectronSFReco[10];
    Double_t ElectronSFReco_up[10];
    Double_t ElectronSFReco_down[10];
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
    
    // variables for electrons
    Int_t nbOfLooseElectrons;
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
    Int_t nbOfLooseMuons;
    Int_t nMuons;
    Int_t rejecteventBadPFmuon;
   // Int_t PFmuon;
    Int_t badmueventmu[10];
    Int_t badmueventclonemu[10];
    Double_t pt_muon[10];
    Double_t pt_muon_corrected[10];
    Double_t ptSF_muon[10];
    Double_t phi_muon[10];
    Double_t eta_muon[10];
    Double_t E_muon[10];
    Double_t TrackLayers_muon[10];
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
    Double_t cdiscCvsB_jet[20];
    Double_t jet_Pt_before_JER[20];
    Double_t jet_Pt_after_JER[20];
    Double_t jet_Pt_before_JES[20];
    Double_t jet_Pt_after_JES[20];
    
    
    
    // met
    Double_t met_Pt;
    Double_t met_Px;
    Double_t met_Py;
    // Double_t met_Pz;
    Double_t met_Phi;
    Double_t met_Eta;
    Double_t met_before_JES;
    Double_t met_after_JES;
    
    
    /*int ZmuIndiceF_0 = -999;
     int ZmuIndiceF_1 = -999;
     int ZelecIndiceF_0 = -999;
     int ZelecIndiceF_1= -999;
     int WmuIndiceF = -999;
     int WelecIndiceF = -999;
     */
    
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
    Bool_t mc_isLastCopy[200];
    Bool_t mc_isPromptFinalState[200];
    Bool_t mc_isHardProcess[200];
    Bool_t mc_fromHardProcessFinalState[200];
    
    
    int nIniRecoLeptons;
    int nIniRecoElectrons;
    int nIniRecoMuons;
    
    globalTree->Branch("el_pt_cut", &el_pt_cut, "el_pt_cut/D");
    globalTree->Branch("el_eta_cut", &el_eta_cut, "el_eta_cut/D");
    globalTree->Branch("el_iso_cone", &el_iso_cone, "el_iso_cone/D");
    globalTree->Branch("mu_pt_cut", &mu_pt_cut , "mu_pt_cut/D");
    globalTree->Branch("mu_eta_cut", &mu_eta_cut, "mu_eta_cut/D");
    globalTree->Branch("mu_iso_cut", &mu_iso_cut, "mu_iso_cut/D");
    globalTree->Branch("jet_eta_cut", &jet_eta_cut, "jet_eta_cut/D");
    globalTree->Branch("jet_pt_cut", &jet_pt_cut,"jet_pt_cut/D");
    
    
    globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
    globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
    globalTree->Branch("nEv" , &nEv, "nEv/I");
    globalTree->Branch("nbSelectedEvents" , &nbSelectedEvents, "nbSelectedEvents/I");
    globalTree->Branch("nbSelectedEvents_3L" , &nbSelectedEvents_3L, "nbSelectedEvents_3L/I");
    globalTree->Branch("nbSelectedEvents_3Lveto" , &nbSelectedEvents_3Lveto, "nbSelectedEvents_3Lveto/I");
    globalTree->Branch("nbSelectedEvents_3Lveto1J" , &nbSelectedEvents_3Lveto1J, "nbSelectedEvents_3Lveto1J/I");
    globalTree->Branch("wnbSelectedEvents" , &wnbSelectedEvents, "wnbSelectedEvents/I");
    globalTree->Branch("wnbSelectedEvents_3L" , &wnbSelectedEvents_3L, "wnbSelectedEvents_3L/I");
    globalTree->Branch("wnbSelectedEvents_3Lveto" , &wnbSelectedEvents_3Lveto, "wnbSelectedEvents_3Lveto/I");
    globalTree->Branch("wnbSelectedEvents_3Lveto1J" , &wnbSelectedEvents_3Lveto1J, "wnbSelectedEvents_3Lveto1J/I");
    globalTree->Branch("sumW", &sumW, "sumW/I");
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
    
    // event related variables
    
    myTree->Branch("x1",&x1,"x1/D");
    myTree->Branch("x2",&x2,"x2/D");
    myTree->Branch("id1",&id1,"id1/I");
    myTree->Branch("id2",&id2,"id2/I");
    myTree->Branch("q",&q,"q/D");
    myTree->Branch("hdamp_up",&hdamp_up,"hdamp_up/D");
    myTree->Branch("hdamp_down",&hdamp_down,"hdamp_down/D");
    myTree->Branch("rejecteventBadPFmuon", &rejecteventBadPFmuon, "rejecteventBadPFmuon/I");
    
    myTree->Branch("channelInt", &channelInt, "channelInt/I");
    myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
    myTree->Branch("run_num",&run_num,"run_num/I");
    myTree->Branch("evt_num",&evt_num,"evt_num/L");
    myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
    myTree->Branch("nvtx",&nvtx,"nvtx/I");
    myTree->Branch("npu",&npu,"npu/I");
    myTree->Branch("puSF",&puSF,"puSF/D");
    myTree->Branch("puSF_up",&puSF_up,"puSF_up/D");
    myTree->Branch("puSF_down",&puSF_down,"puSF_down/D");
    myTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
    myTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
    myTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
    myTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
    myTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
    myTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
    myTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
    myTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
    myTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
    
    myTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
    myTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
    myTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
    myTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
    myTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
    myTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
    myTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
    myTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");
    
    myTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I");
    myTree->Branch("PassedTriggerMET", &PassedTriggerMET, "PassedTriggerMET/I");
    myTree->Branch("PassedTrigger", &PassedTrigger, "PassedTrigger/I");
    myTree->Branch("PassedTriggerJET", &PassedTriggerJET, "PassedTriggerJET/I");
    myTree->Branch("PassedTriggerNoLogic", &PassedTriggerNoLogic, "PassedTriggerNoLogic/I");
    myTree->Branch("PassedTriggerNoLogic2", &PassedTriggerNoLogic2, "PassedTriggerNoLogic2/I");
    myTree->Branch("PassedGoodPV", &PassedGoodPV,"PassedGoodPV/I");
    myTree->Branch("withfakes",&withfakes,"withfakes/I");
    /* myTree->Branch("WmuIndiceF", &WmuIndiceF, "WmuIndiceF/I");
     myTree->Branch("WelecIndiceF", &WelecIndiceF, "WelecIndiceF/I");
     myTree->Branch("ZelecIndiceF_0", &ZelecIndiceF_0, "ZelecIndiceF_0/I");
     myTree->Branch("ZelecIndiceF_1", &ZelecIndiceF_1, "ZelecIndiceF_1/I");
     myTree->Branch("ZmuIndiceF_0", &ZmuIndiceF_0, "ZmuIndiceF_0/I");
     myTree->Branch("ZmuIndiceF_1", &ZmuIndiceF_1, "ZmuIndiceF_1/I");
     
     */
    // electrons
    myTree->Branch("nbOfLooseElectrons",&nbOfLooseElectrons, "nbOfLooseElectrons/I");
    myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//
    myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
    myTree->Branch("ElectronSF_up",&ElectronSF_up,"ElectronSF_up[nElectrons]/D");
    myTree->Branch("ElectronSF_down",&ElectronSF_down,"ElectronSF_down[nElectrons]/D");
    myTree->Branch("ElectronSFID",&ElectronSFID,"ElectronSFID[nElectrons]/D");
    myTree->Branch("ElectronSFID_up",&ElectronSFID_up,"ElectronSFID_up[nElectrons]/D");
    myTree->Branch("ElectronSFID_down",&ElectronSFID_down,"ElectronSFID_down[nElectrons]/D");
    myTree->Branch("ElectronSFReco",&ElectronSFReco,"ElectronSFReco[nElectrons]/D");
    myTree->Branch("ElectronSFReco_up",&ElectronSFReco_up,"ElectronSFReco_up[nElectrons]/D");
    myTree->Branch("ElectronSFReco_down",&ElectronSFReco_down,"ElectronSFReco_down[nElectrons]/D");
    
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
    
    // muons
    myTree->Branch("nbOfLooseMuons", &nbOfLooseMuons, "nbOfLooseMuons/I");
    myTree->Branch("nMuons",&nMuons, "nMuons/I");
    myTree->Branch("badmueventclonemu", &badmueventclonemu, "badmueventclonemu[nMuons]/I");
    myTree->Branch("badmueventmu", &badmueventmu, "badmueventmu[nMuons]/I");
  //  myTree->Branch("PFmuon", &PFmuon, "PFMuon/I");
    
    myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
    myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
    
    myTree->Branch("MuonIDSF_BCDEF",&MuonIDSF_BCDEF,"MuonIDSF_BCDEF[nMuons]/D");
    myTree->Branch("MuonIsoSF_BCDEF",&MuonIsoSF_BCDEF, "MuonIsoSF_BCDEF[nMuons]/D");
    myTree->Branch("MuonIDSF_GH",&MuonIDSF_GH,"MuonIDSF_GH[nMuons]/D");
    myTree->Branch("MuonIsoSF_GH",&MuonIsoSF_GH, "MuonIsoSF_GH[nMuons]/D");
    myTree->Branch("MuonIDSF_BCDEF_up",&MuonIDSF_BCDEF_up,"MuonIDSF_BCDEF_up[nMuons]/D");
    myTree->Branch("MuonIsoSF_BCDEF_up",&MuonIsoSF_BCDEF_up, "MuonIsoSF_BCDEF_up[nMuons]/D");
    myTree->Branch("MuonIDSF_GH_up",&MuonIDSF_GH_up,"MuonIDSF_GH_up[nMuons]/D");
    myTree->Branch("MuonIsoSF_GH_up",&MuonIsoSF_GH_up, "MuonIsoSF_GH_up[nMuons]/D");
    myTree->Branch("MuonIDSF_BCDEF_down",&MuonIDSF_BCDEF_down,"MuonIDSF_BCDEF_down[nMuons]/D");
    myTree->Branch("MuonIsoSF_BCDEF_down",&MuonIsoSF_BCDEF_down, "MuonIsoSF_BCDEF_down[nMuons]/D");
    myTree->Branch("MuonIDSF_GH_down",&MuonIDSF_GH_down,"MuonIDSF_GH_down[nMuons]/D");
    myTree->Branch("MuonIsoSF_GH_down",&MuonIsoSF_GH_down, "MuonIsoSF_GH_down[nMuons]/D");
    
    
     myTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
    myTree->Branch("MuonTrackSF_up",&MuonTrackSF_up, "MuonTrackSF_up[nMuons]/D");
    myTree->Branch("MuonTrackSF_down",&MuonTrackSF_down, "MuonTrackSF_down[nMuons]/D");
    myTree->Branch("MuonIDSF_up",&MuonIDSF_up,"MuonIDSF_up[nMuons]/D");
    myTree->Branch("MuonIsoSF_up",&MuonIsoSF_up, "MuonIsoSF_up[nMuons]/D");
    myTree->Branch("MuonIDSF_down",&MuonIDSF_down,"MuonIDSF_down[nMuons]/D");
    myTree->Branch("MuonIsoSF_down",&MuonIsoSF_down, "MuonIsoSF_down[nMuons]/D");
    myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
    myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
    myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
    myTree->Branch("ptSF_muon",ptSF_muon,"ptSF_muon[nMuons]/D");
    myTree->Branch("pt_muon_corrected",pt_muon_corrected,"pt_muon_corrected[nMuons]/D");
    myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
    myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
    myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
    myTree->Branch("TrackLayers_muon", TrackLayers_muon,"TrackLayers_muon[nMuons]/D");
    myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
    myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
    myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
    myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
    myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
    myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
    myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
    myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
    myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
    
    
    
    // jets
    myTree->Branch("nJets",&nJets,"nJets/I");
    myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
    myTree->Branch("px_jet",px_jet,"px_jet[nJets]/D");
    myTree->Branch("py_jet",py_jet,"py_jet[nJets]/D");
    myTree->Branch("pz_jet",pz_jet,"pz_jet[nJets]/D");
    myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
    myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
    myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
    myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
    myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
    myTree->Branch("jet_Pt_before_JER",jet_Pt_before_JER,"jet_Pt_before_JER[nJets]/D");
    myTree->Branch("jet_Pt_before_JES",jet_Pt_before_JES,"jet_Pt_before_JES[nJets]/D");
    myTree->Branch("jet_Pt_after_JER",jet_Pt_after_JER,"jet_Pt_after_JER[nJets]/D");
    myTree->Branch("jet_Pt_after_JES",jet_Pt_after_JES,"jet_Pt_after_JES[nJets]/D");
    myTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
    myTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
    
    if(!isData){
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
      myTree->Branch("mc_isLastCopy", &mc_isLastCopy, "mc_isLastCopy[nMCParticles]/O");
      myTree->Branch("mc_isPromptFinalState", &mc_isPromptFinalState, "mc_isPromptFinalState[nMCParticles]/O");
      myTree->Branch("mc_isHardProcess", &mc_isHardProcess, "mc_isHardProcess[nMCParticles]/O");
      myTree->Branch("mc_fromHardProcessFinalState", &mc_fromHardProcessFinalState, "mc_fromHardProcessFinalState[nMCParticles]/O");
    }
    
    
    // met
    myTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
    myTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
    myTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
    myTree->Branch("met_Px", &met_Px, "met_Px/D");
    myTree->Branch("met_Py", &met_Py, "met_Py/D");
    //myTree->Branch("met_Pz", &met_Pz, "met_Pz/D");
    myTree->Branch("met_before_JES", &met_before_JES, "met_before_JES/D");
    myTree->Branch("met_after_JES", &met_after_JES, "met_after_JES/D");
    
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
    trigger_met->bookTriggers(isData, dName);
    trigger_jet->bookTriggers(isData, dName);
    
    
    
    //// define trigger branches
    /*
    std::array<Int_t,200> triggers_container;
    // cout << "trigger_mu->triggerList.size() " << trigger_mu->triggerList.size() << endl;
    // for(int iter_trig=0; iter_trig< (isData?trigger->triggerListData.size():trigger->triggerListMC.size()) && iter_trig<200; iter_trig++){
    for(int iter_trig=0; iter_trig< trigger_mu->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_mu->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_mumu->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_mumu->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      // std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_mumumu->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_mumumu->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
      //myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_e->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_e->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_ee->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_ee->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
      //myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_eee->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_eee->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_emu->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_emu->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_emumu_mumue->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_emumu_mumue->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
      //myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_met->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_met->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
     // myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
    for(int iter_trig=0; iter_trig< trigger_jet->triggerList.size() && iter_trig<200; iter_trig++){
      TString trigname;
      trigname = trigger_jet->triggerList[iter_trig];
      //cout << trigname << endl;
      trigname.ReplaceAll("_v*","");
      trigname.ReplaceAll("_v1","");
      trigname.ReplaceAll("_v2","");
      trigname.ReplaceAll("_v3","");
      trigname.ReplaceAll("_v4","");
      trigname.ReplaceAll("_v5","");
      trigname.ReplaceAll("_v6","");
      trigname.ReplaceAll("_v7","");
      trigname.ReplaceAll("_v8","");
      trigname.ReplaceAll("_v9","");
      trigname.ReplaceAll("_v10","");
      trigname.ReplaceAll("_v11","");
      trigname.ReplaceAll("_v12","");
      TString branchname = trigname+"/I";
      //std::cout << "adding trigger to trees " << trigname << " mapped to element " << iter_trig << " " << branchname << std::endl;
      //myTree->Branch(trigname,&(triggers_container[iter_trig]),branchname);
    }
*/
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
    vector <TRootMuon*>       init_muons_corrected;
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
    vector<TRootMuon*>        selectedFakeMuons;
    vector<TRootMuon*>        selectedFakeTightMuons;
    vector<TRootMuon*>        selectedFakeLooseMuons;
    vector<TRootElectron*>    selectedFakeElectrons;
    vector<TRootElectron*>    selectedFakeTightElectrons;
    vector<TRootElectron*>    selectedFakeLooseElectrons;
    
    
    vector<TRootMCParticle*> mcParticles;
    

    //////////////////////////////////////
    // Begin Event Loop
    //////////////////////////////////////
    nbEvents = 0;
    nofPosWeights = 0;
    nofNegWeights = 0;
    double eventweight = 1.;
    bool continueFlow ;
    nbSelectedEvents = 0;
    nbSelectedEvents_3L = 0;
    nbSelectedEvents_3Lveto = 0;
    nbSelectedEvents_3Lveto1J = 0;
    wnbSelectedEvents = 0;
    wnbSelectedEvents_3L = 0;
    wnbSelectedEvents_3Lveto = 0;
    wnbSelectedEvents_3Lveto1J = 0;
    
    bool debug = false;
    vector <int> selections;
    bool   passedMET = false;
    PassedGoodPV = false;
    bool   HBHEnoise = false;
    bool   HBHEIso = false;
    bool   CSCTight = false;
    bool badchan = false;
    bool badmu = false;
    
    bool   EcalDead = false;
    bool    eeBad = false;
    bool   lep3 = false;
    bool lep3veto = false;
    bool lep3veto1J = false;
    bool lep2 = false;
    TLorentzVector metTLV;
    TLorentzVector metTLVbf;
    string TriggBits;
    double met;
    bool leptonsAssigned ;
    
    
    
    for (unsigned int ievt = event_start; ievt < end_d; ievt++)
    {
      rejecteventBadPFmuon = 0;
      eventSelected = false;
      continueFlow = true;
      lep3 = false;
      lep2 = false;
      lep3veto = false;
      lep3veto1J = false;
      met = 0.;
      
      TriggBits = "";
      
      metTLV.Clear();
      metTLVbf.Clear();
      metTLV.SetPxPyPzE(0,0,0,0);
      selections.clear();
      bool lepsel = false;
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
      
      passedMET = false;
      PassedGoodPV = false;
      HBHEnoise = false;
      HBHEIso = false;
      CSCTight = false;
      badchan = false;
      badmu = false;
      EcalDead = false;
      eeBad = false;

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
      
      double scaleFactor = 1.;  // scale factor for the event
      
      
      event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
      if(debug) cout << "event loaded" << endl;
      genjets.clear();
      if(!isData){
        genjets = treeLoader.LoadGenJet(ievt,false);  //needed for JER
      }
      init_jets_corrected = init_jets;
      
      if(false)
      {
        cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
        cout <<"Number of Muons Loaded: " << init_muons.size() <<endl;
        cout << "Number of Jets  Loaded: " << init_jets.size() << endl;
	cout << "met size " << mets.size() << endl; 
        cout << "Met px / py loaded: "<< mets[0]->Px() << " / " << mets[0]->Py() << endl;
      }
      
      
      nIniRecoLeptons=0;
      nIniRecoElectrons = 0;
      nIniRecoMuons = 0;
      double dataMuPTSF = 1.;
      double mcMuPTSF = 1.;
      init_muons_corrected = init_muons;
      for(int iMu = 0 ; iMu < init_muons.size(); iMu++){
        if(init_muons[iMu]->Pt() > 10.0) nIniRecoMuons++;
        if(init_muons[iMu]->isBad80X()) rejecteventBadPFmuon = 1;
        if(init_muons[iMu]->isClone80X()) rejecteventBadPFmuon = 1;
       // if(!init_muons[iMu]->isPFMuon()) { cout << "rejected muon " << endl; rejecteventBadPFmuon = 1;}
       // PFmuon =   init_muons[0]->isPFMuon();
        
        
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
      eeBad = event->getEEBadScFilter();
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
      if(!isData )
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
      
      //save variables for pdf uncertainties and factorisation scale calculation
      hdamp_down = 1.;
      hdamp_up = 1.;
      Double_t weight0, weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8;
      weight0 = weight1 = weight2 = weight3 = weight4 = weight5 = weight6 = weight7 = weight8 = 1.;
      if(!isData){
        if(event->getWeight(1)!= -9999)
        {
          weight0 = (event->getWeight(1))/(abs(event->originalXWGTUP()));
          weight1 = (event->getWeight(2))/(abs(event->originalXWGTUP()));
          weight2 = (event->getWeight(3))/(abs(event->originalXWGTUP()));
          weight3 = (event->getWeight(4))/(abs(event->originalXWGTUP()));
          weight4 = (event->getWeight(5))/(abs(event->originalXWGTUP()));
          weight5 = (event->getWeight(6))/(abs(event->originalXWGTUP()));
          weight6 = (event->getWeight(7))/(abs(event->originalXWGTUP()));
          weight7 = (event->getWeight(8))/(abs(event->originalXWGTUP()));
          weight8 = (event->getWeight(9))/(abs(event->originalXWGTUP()));
        }
        else if (event->getWeight(1001)!= -9999)
        {
          weight0 = (event->getWeight(1001))/(abs(event->originalXWGTUP()));
          weight1 = (event->getWeight(1002))/(abs(event->originalXWGTUP()));
          weight2 = (event->getWeight(1003))/(abs(event->originalXWGTUP()));
          weight3 = (event->getWeight(1004))/(abs(event->originalXWGTUP()));
          weight4 = (event->getWeight(1005))/(abs(event->originalXWGTUP()));
          weight5 = (event->getWeight(1006))/(abs(event->originalXWGTUP()));
          weight6 = (event->getWeight(1007))/(abs(event->originalXWGTUP()));
          weight7 = (event->getWeight(1008))/(abs(event->originalXWGTUP()));
          weight8 = (event->getWeight(1009))/(abs(event->originalXWGTUP()));
        }
        
        // hdamp variation
        if (event->getWeight(1001)!= -9999)
        {
          hdamp_up = event->getWeight(5019)/fabs(event->originalXWGTUP());
          hdamp_down = event->getWeight(5010)/fabs(event->originalXWGTUP());
        }
        
        x1 = event->xParton1();
        x2 = event->xParton2();
        id1 = event->idParton1();
        id2 = event->idParton2();
        q = event->factorizationScale();
        
      }
    
    
    
      ///////////////////////////////////////////
      //  Trigger
      ///////////////////////////////////////////
      PassedTriggerNoLogic = false;
      PassedTriggerNoLogic2 = false;
      PassedTrigger = false;
      PassedTriggerMET = false;
      PassedTriggerJET = false;
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
      bool trigged_met = false;
      bool trigged_jet = false;
      if(runHLT )
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
        trigger_met->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_met =  trigger_met->checkIfFired();
        trigger_jet->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
        trigged_jet =  trigger_jet->checkIfFired();
        
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
        bool EM2 = false;
        bool MM2 = false;
        bool EE2 = false;
        
        EM2 = trigged_emu;
        MM2 = trigged_mumu;
        EE2 = trigged_ee;
        
        EM = (trigged_emumu_mumue|| trigged_emu);
        MM = (trigged_mumu || trigged_mumumu ) ;
        EE = (trigged_ee || trigged_eee );
        M  = ( trigged_mu );
        E  = (trigged_e);
        
        
        // testing TO FIX
        if ( EM2 &&                 (emdataset) ) PassedTriggerNoLogic2 = 1;
        if ( MM2 &&                 (mmdataset) ) PassedTriggerNoLogic2 = 1;
        if ( EE2 &&                 (eedataset) ) PassedTriggerNoLogic2 = 1;
        
        if ( EM &&                  (emdataset) ) PassedTriggerNoLogic = 1;
        if ( MM &&                  (mmdataset) ) PassedTriggerNoLogic = 1;
        if ( EE &&                  (eedataset) ) PassedTriggerNoLogic = 1;
        if ( M  &&                  (mdataset ) ) PassedTriggerNoLogic = 1;
        if ( E  &&                  (edataset ) ) PassedTriggerNoLogic = 1;
        
        //for data
        if ( EM &&                               (emdataset) ) result_trigger = 1;
        if ( MM && !EM &&                        (mmdataset) ) result_trigger = 1;
        if ( EE && !EM && !MM &&                 (eedataset) ) result_trigger = 1;
        if ( M  && !EM && !MM && !EE &&          (mdataset ) ) result_trigger = 1;
        if ( E  && !EM && !MM && !EE && !M &&    (edataset ) ) result_trigger = 1;
        // for MC
        if ( ( EM || MM || EE || E || M )    &&   !isData ) result_trigger = 1;
        if ( ( EM || MM || EE || E || M )    &&   (dName.find("data_MET")!=string::npos || dName.find("data_Jet")!=string::npos) ) result_trigger = 1;
        
        trigged = result_trigger;
        // if(dName.find("NP")!=string::npos) trigged = true; // needs to be fixed with the new MC
        
        
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
        trigged_jet = true;
        trigged_met = true;
      }
      
      if(verbose > 1) cout << "Apply trigger? " << runHLT << " trigged? " << trigged << endl;
      
      PassedTrigger = trigged;
      PassedTriggerMET = trigged_met;
      PassedTriggerJET = trigged_jet;
      ////////////////////////////
      ///// JES - JER smearing     ////
      //////////////////////////
      JERon = 0;
      
      for(int iJ = 0; iJ < init_jets_corrected.size(); iJ++){
        jet_Pt_before_JER[iJ] = init_jets_corrected[iJ]->Pt();
      }
      if(applyJER && !isData)
      {
        // cout << "applying JER" << endl;
        if(doJESJERshift == 3) jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "minus", false);
        else if(doJESJERshift == 4) jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "plus", false);
        else jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal", false);
        JERon = 1;
      }
      for(int iJ = 0; iJ < init_jets_corrected.size(); iJ++){
        jet_Pt_after_JER[iJ] = init_jets_corrected[iJ]->Pt();
        jet_Pt_before_JES[iJ] = init_jets_corrected[iJ]->Pt();
      }
      met_before_JES = mets[0]->Pt();
      JESon = 0;
      METon = 0;
      if(applyJES)
      {
        // cout << "applying JES" << endl;
        if(doJESJERshift == 1)  jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "minus");
        else if(doJESJERshift == 2)  jetTools->correctJetJESUnc(init_jets_corrected, mets[0], "plus");
        
        
         jetTools->correctJetsMet(init_jets_corrected,mets[0],event->fixedGridRhoFastjetAll() ,isData);
        
        
        METon = 1;
        JESon = 1;
      }
      met_after_JES = mets[0]->Pt();
      for(int iJ = 0; iJ < init_jets_corrected.size(); iJ++){
        jet_Pt_after_JES[iJ] = init_jets_corrected[iJ]->Pt();
        
      }
           ///////////////////////////////////////////////////////////
      // Event selection
      ///////////////////////////////////////////////////////////
      
      // Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons_corrected, init_electrons, mets,event->fixedGridRhoFastjetAll());
      PreselectedJets.clear();
      PreselectedJets  = selection.GetSelectedJets(jet_pt_cut,jet_eta_cut, true, "Loose");
      selectedMuons.clear();
      selectedLooseMuons.clear();
      selectedFakeMuons.clear();
      selectedFakeLooseMuons.clear();
      selectedFakeTightMuons.clear();
      
      
      
      
      selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Tight", "Summer16");   // spring 15 still counts for 2016
      selectedLooseMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut_loose, "Loose", "Summer16"); // spring 15 still counts for 2016
      selectedFakeTightMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut_faketight, "Fake", "Summer16");
      selectedFakeLooseMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut_fakeloose, "Fake", "Summer16");
      
      if(!doFakeshift) selectedFakeMuons = selectedFakeTightMuons;
      else if(doFakeshift) selectedFakeMuons = selectedFakeLooseMuons;
      // pt, eta, iso // run normally
      selectedElectrons.clear();
      selectedVetoElectrons.clear();
      selectedFakeElectrons.clear();
      selectedFakeLooseElectrons.clear();
      selectedFakeTightElectrons.clear();
      selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Tight","Spring16_80X",true,true);// pt, eta, WP point, campaign, cutbased, VID EA
      selectedVetoElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Veto","Spring16_80X",true,true);// pt, eta
      selectedFakeTightElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "FakeTight","Spring16_80X",true,false);
      selectedFakeLooseElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "FakeLoose","Spring16_80X",true,false);
      
      if(!doFakeshift) selectedFakeElectrons = selectedFakeTightElectrons;
      else if(doFakeshift) selectedFakeElectrons = selectedFakeLooseElectrons;
      
      /// For MC Information
      mcParticles.clear();
      if(!isData) treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
      if(!isData) sort(mcParticles.begin(),mcParticles.end(),HighestPt());
      
      
      
      
      
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
      
      
      
      
      ////////////////////////////////////////////////
      // Pre cut operations
      ////////////////////////////////////////////////
      // Apply primary vertex selection
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
      // Met filters
      if(isData && HBHEnoise && HBHEIso && CSCTight && EcalDead && eeBad && isGoodPV && badchan && badmu) passedMET = true;
      else if(!isData && HBHEnoise && HBHEIso && CSCTight && EcalDead  && isGoodPV && badchan && badmu ) passedMET = true;
      PassedMETFilter = passedMET;
      PassedGoodPV = isGoodPV;
      
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
     
      double btagWeightShape = 1.;
      Double_t btagWeight_shape_up_lf= 1.;
      Double_t btagWeight_shape_down_lf= 1.;
      Double_t btagWeight_shape_up_hf= 1.;
      Double_t btagWeight_shape_down_hf= 1.;
      Double_t btagWeight_shape_up_hfstats1= 1.;
      Double_t btagWeight_shape_down_hfstats1= 1.;
      Double_t btagWeight_shape_up_hfstats2= 1.;
      Double_t btagWeight_shape_down_hfstats2= 1.;
      Double_t btagWeight_shape_up_lfstats1= 1.;
      Double_t btagWeight_shape_down_lfstats1= 1.;
      Double_t btagWeight_shape_up_lfstats2= 1.;
      Double_t btagWeight_shape_down_lfstats2= 1.;
      Double_t btagWeight_shape_up_cferr1= 1.;
      Double_t btagWeight_shape_down_cferr1= 1.;
      Double_t btagWeight_shape_up_cferr2= 1.;
      Double_t btagWeight_shape_down_cferr2= 1.;
      double bTagEff = 1. , bTagEff_LFUp = 1., bTagEff_LFDown = 1., bTagEff_HFUp = 1., bTagEff_HFDown = 1., bTagEff_HFStats1Up = 1.,
      bTagEff_HFStats1Down = 1., bTagEff_HFStats2Up = 1., bTagEff_HFStats2Down = 1., bTagEff_LFStats1Up = 1., bTagEff_LFStats1Down = 1.,
      bTagEff_LFStats2Up = 1., bTagEff_LFStats2Down = 1., bTagEff_CFErr1Up = 1., bTagEff_CFErr1Down = 1., bTagEff_CFErr2Up = 1., bTagEff_CFErr2Down = 1.;
      
      if( !isData )
      {
        bool isBFlav = false;
        bool isLFlav = false;
        bool isCFlav = false;
        double jetpt ;
        double jeteta;
        double jetdisc ;
        int jetpartonflav ;
        for(int intJet = 0; intJet < selectedJets.size(); intJet++)
        {
          jetpt = selectedJets[intJet]->Pt();
          if(jetpt > 1000.) jetpt = 999.;
          jeteta = selectedJets[intJet]->Eta();
          jetdisc = selectedJets[intJet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
          if(jetdisc<0.0) jetdisc = -0.05;
          if(jetdisc>1.0) jetdisc = 1.0;
          
          isBFlav = false;
          isCFlav = false;
          isLFlav = false;
          BTagEntry::JetFlavor jflav;
          jetpartonflav = std::abs(selectedJets[intJet]->partonFlavour());
          if(debug) cout<<"parton flavour: "<<jetpartonflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
          
          if(jetpartonflav == 5){
            jflav = BTagEntry::FLAV_B;
            isBFlav = true;
          }
          else if(jetpartonflav == 4){
            jflav = BTagEntry::FLAV_C;
            isCFlav = true;
          }
          else{
            jflav = BTagEntry::FLAV_UDSG;
            isLFlav = true;
          }
          
          
          if(doJESJERshift == 1 && !isCFlav) bTagEff = reader_csvv2_JESdown->eval(jflav, jeteta, jetpt, jetdisc);
          else if(doJESJERshift == 2 && !isCFlav) bTagEff = reader_csvv2_JESup->eval(jflav, jeteta, jetpt, jetdisc);
          else bTagEff = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          
          
          
          // other systematics
          if( isBFlav ) bTagEff_LFUp = reader_csvv2_LFUp->eval(jflav, jeteta, jetpt, jetdisc);
          if( isBFlav ) bTagEff_LFDown = reader_csvv2_LFDown->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_HFUp = reader_csvv2_HFUp->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_HFDown = reader_csvv2_HFDown->eval(jflav, jeteta, jetpt, jetdisc);
          if( isBFlav ) bTagEff_HFStats1Up = reader_csvv2_HFStats1Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isBFlav ) bTagEff_HFStats1Down = reader_csvv2_HFStats1Down->eval(jflav, jeteta, jetpt, jetdisc);
          if( isBFlav ) bTagEff_HFStats2Up = reader_csvv2_HFStats2Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isBFlav ) bTagEff_HFStats2Down = reader_csvv2_HFStats2Down->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_LFStats1Up = reader_csvv2_LFStats1Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_LFStats1Down = reader_csvv2_LFStats1Down->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_LFStats2Up = reader_csvv2_LFStats2Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isLFlav ) bTagEff_LFStats2Down = reader_csvv2_LFStats2Down->eval(jflav, jeteta, jetpt, jetdisc);
          if( isCFlav ) bTagEff_CFErr1Up = reader_csvv2_CFErr1Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isCFlav ) bTagEff_CFErr1Down = reader_csvv2_CFErr1Down->eval(jflav, jeteta, jetpt, jetdisc);
          if( isCFlav ) bTagEff_CFErr2Up = reader_csvv2_CFErr2Up->eval(jflav, jeteta, jetpt, jetdisc);
          if( isCFlav ) bTagEff_CFErr2Down = reader_csvv2_CFErr2Down->eval(jflav, jeteta, jetpt, jetdisc);
          
          //If jet is not the appropriate flavor for that systematic, use the nominal reader so that all weights will be on the same
          //jet multiplicity footing. see code Kevin
          if( !isBFlav ) bTagEff_LFUp = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isBFlav ) bTagEff_LFDown = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_HFUp = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_HFDown = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isBFlav ) bTagEff_HFStats1Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isBFlav ) bTagEff_HFStats1Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isBFlav ) bTagEff_HFStats2Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isBFlav ) bTagEff_HFStats2Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_LFStats1Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_LFStats1Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_LFStats2Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isLFlav ) bTagEff_LFStats2Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isCFlav ) bTagEff_CFErr1Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isCFlav ) bTagEff_CFErr1Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isCFlav ) bTagEff_CFErr2Up = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          if( !isCFlav ) bTagEff_CFErr2Down = reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);
          
          
          
          // fill btagweights
          btagWeightShape *= bTagEff;
          
          btagWeight_shape_up_lf *= bTagEff_LFUp;
          btagWeight_shape_down_lf *= bTagEff_LFDown;
          btagWeight_shape_up_hf *= bTagEff_HFUp;
          btagWeight_shape_down_hf *= bTagEff_HFDown;
          btagWeight_shape_up_hfstats1 *= bTagEff_HFStats1Up;
          btagWeight_shape_down_hfstats1 *= bTagEff_HFStats1Down;
          btagWeight_shape_up_hfstats2 *= bTagEff_HFStats2Up;
          btagWeight_shape_down_hfstats2 *= bTagEff_HFStats2Down;
          btagWeight_shape_up_lfstats1 *= bTagEff_LFStats1Up;
          btagWeight_shape_down_lfstats1 *= bTagEff_LFStats1Down;
          btagWeight_shape_up_lfstats2 *= bTagEff_LFStats2Up;
          btagWeight_shape_down_lfstats2 *= bTagEff_LFStats2Down;
          btagWeight_shape_up_cferr1 *= bTagEff_CFErr1Up;
          btagWeight_shape_down_cferr1 *= bTagEff_CFErr1Down;
          btagWeight_shape_up_cferr2 *= bTagEff_CFErr2Up;
          btagWeight_shape_down_cferr2 *= bTagEff_CFErr2Down;
        }
        
      }
      
      
      
      
      ////////////////////////////////////
      //   Determine eventweight        ///
      /////////////////////////////////
      
         //////////////////////////////////////////////////////
      // Applying baseline selection
      //////////////////////////////////////////////////////
      continueFlow = true;
      nbEvents++;

      
      if(!doFakeLepton){
        if(((selectedMuons.size() + selectedElectrons.size()) < 2)){
          selections.push_back(0);
          continueFlow = false;
        }
        else if((selectedMuons.size() + selectedElectrons.size()) >1){
          selections.push_back(1);
          
          if(selectedMuons.size() + selectedElectrons.size() ==3 ){ lep3 = true; nbSelectedEvents_3L++;}
         
          if(selectedMuons.size() == 3) {channelInt = 0; i_channel = 0;}
          else if(selectedElectrons.size() == 3) {channelInt = 3; i_channel = 3;}
          else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1) {channelInt = 2; i_channel = 2; }
          else if(selectedMuons.size() == 2 && selectedElectrons.size() == 1){channelInt = 1; i_channel = 1; }
          else if(selectedMuons.size() == 2) {channelInt = 4; i_channel = 4;}
          else if(selectedElectrons.size() == 2) {channelInt = 5; i_channel = 5;}
          //else {cout << "ERROR no channel selected" << endl; break; }
        }
        
        
        
        
      }
      else  if(doFakeLepton){
        if((selectedMuons.size() + selectedElectrons.size()) == 2 && (selectedFakeMuons.size() + selectedFakeElectrons.size()) == 1 ){
          selections.push_back(1);
          nbSelectedEvents_3L++;
          lep3 = true;
          if(selectedMuons.size() == 2 && selectedFakeMuons.size() == 1) {channelInt = 0; i_channel = 0;}
          else if(selectedElectrons.size() == 2 && selectedFakeElectrons.size() == 1) {channelInt = 3; i_channel = 3;}
          else if(selectedElectrons.size() == 1 && selectedMuons.size() == 1 && selectedFakeElectrons.size() == 1 ) {channelInt = 2; i_channel = 2; }
          else if(selectedElectrons.size() == 2 && selectedFakeMuons.size() == 1  ) {channelInt = 2; i_channel = 2; }
          else if(selectedMuons.size() == 2 && selectedFakeElectrons.size() == 1){channelInt = 1; i_channel = 1; }
          else if(selectedMuons.size() == 1 && selectedElectrons.size() == 1  && selectedFakeMuons.size() == 1){channelInt = 1; i_channel = 1; }
          //else {cout << "ERROR no channel selected" << endl; break; }
        }
        else{
          selections.push_back(0);
          continueFlow = false;
        }
        
      }
      
      nbOfLooseMuons = selectedLooseMuons.size();
      nbOfLooseElectrons = selectedVetoElectrons.size(); 
      
     /* if((selectedMuons.size() == 3 && selectedLooseMuons.size() != 3) || (selectedVetoElectrons.size() != 3 &&  selectedElectrons.size()== 3)){
        selections.push_back(0);
        continueFlow = false;
      }
      else {
        selections.push_back(1);
        
        if(continueFlow) {
          lepsel = true;
          lep3veto = true;
          nbSelectedEvents_3Lveto++;
        }
      }*/
      
      double met_px = mets[0]->Px();
      double met_py = mets[0]->Py();
      met_Pt = sqrt(met_px*met_px + met_py*met_py);
      met = met_Pt;
      met_Phi = mets[0]->Phi();
      met_Eta = mets[0]->Eta();
      
      if(!isData)puSF =  LumiWeights.ITweight((int)event->nTruePU());
      if(!isData)puSF_up =  LumiWeights_up.ITweight((int)event->nTruePU());
      if(!isData)puSF_down =  LumiWeights_down.ITweight((int)event->nTruePU());
      
      
      
      btagSFshape = 1.;
      btagSFshape_down_cferr1 = 1.;
      btagSFshape_down_cferr2 = 1.;
      btagSFshape_down_hf= 1.;
      btagSFshape_down_hfstats1 = 1.;
      btagSFshape_down_hfstats2 = 1.;
      btagSFshape_down_lf = 1.;
      btagSFshape_down_lfstats1 = 1.;
      btagSFshape_down_lfstats2 = 1.;
      
      btagSFshape_up_cferr1 = 1.;
      btagSFshape_up_cferr2 = 1.;
      btagSFshape_up_hf= 1.;
      btagSFshape_up_hfstats1 = 1.;
      btagSFshape_up_hfstats2 = 1.;
      btagSFshape_up_lf = 1.;
      btagSFshape_up_lfstats1 = 1.;
      btagSFshape_up_lfstats2 = 1.;
      
      
      
      
      if(!isData){
        
        btagSFshape = btagWeightShape;
        btagSFshape_up_lf = btagWeight_shape_up_lf;
        btagSFshape_down_lf = btagWeight_shape_down_lf;
        btagSFshape_up_hf = btagWeight_shape_up_hf;
        btagSFshape_down_hf = btagWeight_shape_down_hf;
        btagSFshape_up_hfstats1 = btagWeight_shape_up_hfstats1;
        btagSFshape_down_hfstats1 = btagWeight_shape_down_hfstats1;
        btagSFshape_up_hfstats2 = btagWeight_shape_up_hfstats2;
        btagSFshape_down_hfstats2 = btagWeight_shape_down_hfstats2;
        btagSFshape_up_lfstats1 = btagWeight_shape_up_lfstats1;
        btagSFshape_down_lfstats1 = btagWeight_shape_down_lfstats1;
        btagSFshape_up_lfstats2 = btagWeight_shape_up_lfstats2;
        btagSFshape_down_lfstats2 = btagWeight_shape_down_lfstats2;
        btagSFshape_up_cferr1 = btagWeight_shape_up_cferr1;
        btagSFshape_down_cferr1 = btagWeight_shape_down_cferr1;
        btagSFshape_up_cferr2 = btagWeight_shape_up_cferr2;
        btagSFshape_down_cferr2 = btagWeight_shape_down_cferr2;
      }
      else if(isData) {btagSFshape = 1.;   puSF_down = 1.; puSF_up = 1.; puSF = 1.;}
      
      
      if( selectedJets.size() >0){
        if(continueFlow){ eventSelected = true; nbSelectedEvents_3Lveto1J++; lep3veto1J = true; }
      }
      
      withfakes = doFakeLepton;
      
      
      //////////////////////////////////////
      //  DO STUFF WITH SELECTED EVENTS ////
      /////////////////////////////////////
      // fill the tree
      
      if(eventSelected){
        eventweight = 1.;
        eventweight *= puSF;
        eventweight *= btagWeightShape;
        /*
        triggers_container.fill(0);
        for(int iter_trig=0; iter_trig< trigger_mu->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_mu->triggermap.find(trigger_mu->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_mumu->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_mumu->triggermap.find(trigger_mumu->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_mumumu->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_mumumu->triggermap.find(trigger_mumumu->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_e->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_e->triggermap.find(trigger_e->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_ee->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_ee->triggermap.find(trigger_eee->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_eee->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_eee->triggermap.find(trigger_eee->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_emu->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_emu->triggermap.find(trigger_emu->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_emumu_mumue->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_emumu_mumue->triggermap.find(trigger_emumu_mumue->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_met->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_met->triggermap.find(trigger_met->triggerList[iter_trig])->second.second;
        }
        for(int iter_trig=0; iter_trig< trigger_jet->triggerList.size() && iter_trig<200; iter_trig++){
          triggers_container[iter_trig] = trigger_jet->triggermap.find(trigger_jet->triggerList[iter_trig])->second.second;
        }
        */
        
        
        nMuons = 0;
        double muonSFtemp = 1.;
        
        for (Int_t selmu =0; selmu < selectedMuons.size() ; selmu++ )
        {
          pt_muon[nMuons]=selectedMuons[selmu]->Pt();
          phi_muon[nMuons]=selectedMuons[selmu]->Phi();
          eta_muon[nMuons]=selectedMuons[selmu]->Eta();
          E_muon[nMuons]=selectedMuons[selmu]->E();
          TrackLayers_muon[nMuons] = selectedMuons[selmu]->nofTrackerLayersWithMeasurement();
          if(isData) ptSF_muon[nMuons] = rc.kScaleDT(selectedMuons[selmu]->charge(), selectedMuons[selmu]->Pt(), selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Phi(), 0, 0);
          else ptSF_muon[nMuons] = rc.kScaleAndSmearMC(selectedMuons[selmu]->charge(), selectedMuons[selmu]->Pt(), selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Phi(), selectedMuons[selmu]->nofTrackerLayersWithMeasurement(),gRandom->Rndm(),gRandom->Rndm(),0, 0);
          
          pt_muon_corrected[nMuons]=selectedMuons[selmu]->Pt()*ptSF_muon[nMuons];
          
         
          
          badmueventmu[nMuons] = selectedMuons[selmu]->isBad80X();
          badmueventclonemu[nMuons] = selectedMuons[selmu]->isClone80X();
          
          
          pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
          if(!isData)
          {
            
            MuonIDSF_BCDEF[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
            MuonIsoSF_BCDEF[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
            MuonIDSF_GH[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
            MuonIsoSF_GH[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
            
            MuonIDSF_BCDEF_down[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
            MuonIsoSF_BCDEF_down[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
            MuonIDSF_GH_down[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
            MuonIsoSF_GH_down[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
            
            MuonIDSF_BCDEF_up[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
            MuonIsoSF_BCDEF_up[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
            MuonIDSF_GH_up[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
            MuonIsoSF_GH_up[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
            
            
            
            MuonIDSF[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            MuonIsoSF[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            
            MuonIDSF_up[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            MuonIsoSF_up[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            
            MuonIDSF_down[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            MuonIsoSF_down[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
            MuonTrackSF[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta());
            MuonTrackSF_up[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*1.01;
            MuonTrackSF_down[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*0.99;
            
            
            eventweight *= MuonIDSF[nMuons] * MuonIsoSF[nMuons] * MuonTrackSF[nMuons];
          }
          else
          {
            MuonIDSF[nMuons] = 1.;
            MuonIsoSF[nMuons] = 1.;
            MuonIDSF_up[nMuons] = 1.;
            MuonIsoSF_up[nMuons] = 1.;
            MuonIDSF_down[nMuons] = 1.;
            MuonIsoSF_down[nMuons] = 1.;
            MuonTrackSF[nMuons] = 1.;
            MuonTrackSF_up[nMuons] = 1.;
            MuonTrackSF_down[nMuons] = 1.;
          }
          if(MuonIDSF[nMuons]*MuonIsoSF[nMuons] == 0 ) cout << "  MuonIDSF[nMuons] " <<  MuonIDSF[nMuons] << " MuonIsoSF[nMuons] " << MuonIsoSF[nMuons] << "  MuonIDSF[nMuons]*MuonIsoSF[nMuons] " <<    MuonIDSF[nMuons]*MuonIsoSF[nMuons]     << endl;
          if(muonSFtemp == 0) cout << " muon SF " << muonSFtemp * MuonIDSF[nMuons]*MuonIsoSF[nMuons] << endl;
          charge_muon[nMuons]=selectedMuons[selmu]->charge();
          nMuons++;
        }
        if(doFakeLepton){
          for (Int_t selmu =0; selmu < selectedFakeMuons.size() ; selmu++ )
          {
            pt_muon[nMuons]=selectedFakeMuons[selmu]->Pt();
            phi_muon[nMuons]=selectedFakeMuons[selmu]->Phi();
            eta_muon[nMuons]=selectedFakeMuons[selmu]->Eta();
            E_muon[nMuons]=selectedFakeMuons[selmu]->E();
            badmueventmu[nMuons] = selectedFakeMuons[selmu]->isBad80X();
            badmueventclonemu[nMuons] = selectedFakeMuons[selmu]->isClone80X();
            TrackLayers_muon[nMuons] = selectedFakeMuons[selmu]->nofTrackerLayersWithMeasurement();
            if(isData) ptSF_muon[nMuons] = rc.kScaleDT(selectedFakeMuons[selmu]->charge(), selectedFakeMuons[selmu]->Pt(), selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Phi(), 0, 0);
            else ptSF_muon[nMuons] = rc.kScaleAndSmearMC(selectedFakeMuons[selmu]->charge(), selectedFakeMuons[selmu]->Pt(), selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Phi(), selectedFakeMuons[selmu]->nofTrackerLayersWithMeasurement(),gRandom->Rndm(),gRandom->Rndm(),0, 0);
            
            pt_muon_corrected[nMuons]=selectedFakeMuons[selmu]->Pt()*ptSF_muon[nMuons];
            
            
            pfIso_muon[nMuons]=selectedFakeMuons[selmu]->relPfIso(4,0);
            if(!isData)
            {
              MuonIDSF_BCDEF[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
              MuonIsoSF_BCDEF[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
              MuonIDSF_GH[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
              MuonIsoSF_GH[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
              
              MuonIDSF_BCDEF_down[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
              MuonIsoSF_BCDEF_down[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
              MuonIDSF_GH_down[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
              MuonIsoSF_GH_down[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1);
              
              MuonIDSF_BCDEF_up[nMuons] = muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
              MuonIsoSF_BCDEF_up[nMuons] = muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
              MuonIDSF_GH_up[nMuons] = muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
              MuonIsoSF_GH_up[nMuons] = muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1);
              
              
              MuonIDSF[nMuons]  = (muonSFWeightID_BCDEF->at(selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              MuonIsoSF[nMuons] = (muonSFWeightIso_BCDEF->at(selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedFakeMuons[selmu]->Eta(), selectedFakeMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              
              MuonIDSF_up[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              MuonIsoSF_up[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              
              MuonIDSF_down[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              MuonIsoSF_down[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
              MuonTrackSF[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta());
              MuonTrackSF_up[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*1.01;
              MuonTrackSF_down[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*0.99;
              
              eventweight *= MuonIDSF[nMuons] * MuonIsoSF[nMuons] *MuonTrackSF[nMuons];
            }
            else
            {
              MuonIDSF[nMuons] = 1.;
              MuonIsoSF[nMuons] = 1.;
              MuonIDSF_up[nMuons] = 1.;
              MuonIsoSF_up[nMuons] = 1.;
              MuonIDSF_down[nMuons] = 1.;
              MuonIsoSF_down[nMuons] = 1.;
              MuonTrackSF[nMuons] = 1.;
              MuonTrackSF_up[nMuons] = 1.;
              MuonTrackSF_down[nMuons] = 1.;

            }
            if(MuonIDSF[nMuons]*MuonIsoSF[nMuons] == 0 ) cout << "  MuonIDSF[nMuons] " <<  MuonIDSF[nMuons] << " MuonIsoSF[nMuons] " << MuonIsoSF[nMuons] << "  MuonIDSF[nMuons]*MuonIsoSF[nMuons] " <<    MuonIDSF[nMuons]*MuonIsoSF[nMuons]     << endl;
            if(muonSFtemp == 0) cout << " muon SF " << muonSFtemp * MuonIDSF[nMuons]*MuonIsoSF[nMuons] << endl;
            charge_muon[nMuons]=selectedFakeMuons[selmu]->charge();
            nMuons++;
          }
        }
        
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
          if(!isData){
            
            ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
            ElectronSF[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
            ElectronSF_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
            ElectronSF_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
            
            ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
            ElectronSFID[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
            ElectronSFID_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
            ElectronSFID_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
            
            ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
            ElectronSFReco[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
            ElectronSFReco_up[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
            ElectronSFReco_down[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
           
          }
          else {
            ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
            ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
            ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
            
          
          }
          
          eventweight *= ElectronSF[nElectrons];
          nElectrons++;
        }
        if(doFakeLepton){
          for (Int_t selel =0; selel < selectedFakeElectrons.size() ; selel++ )
          {
            pt_electron[nElectrons]=selectedFakeElectrons[selel]->Pt();
            phi_electron[nElectrons]=selectedFakeElectrons[selel]->Phi();
            eta_electron[nElectrons]=selectedFakeElectrons[selel]->Eta();
            eta_superCluster_electron[nElectrons]=selectedFakeElectrons[selel]->superClusterEta();
            E_electron[nElectrons]=selectedFakeElectrons[selel]->E();
            pfIso_electron[nElectrons]=selectedFakeElectrons[selel]->relPfIso(3,0);
            charge_electron[nElectrons]=selectedFakeElectrons[selel]->charge();
            if(!isData){
              ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
              
              ElectronSF[nElectrons] = electronSFWeightID->at(selectedFakeElectrons[selel]->Eta(),selectedFakeElectrons[selel]->Pt(),0)*electronSFWeightReco->at(selectedFakeElectrons[selel]->Eta(),selectedFakeElectrons[selel]->Pt(),0);
              ElectronSF_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
              ElectronSF_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1)*electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
              
              ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
              ElectronSFID[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
              ElectronSFID_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
              ElectronSFID_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
              
              ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
              ElectronSFReco[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
              ElectronSFReco_up[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),1);
              ElectronSFReco_down[nElectrons] = electronSFWeightReco->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),-1);
              
              
            }
            else {
              ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
              ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
              ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
              
            }
            
            eventweight *= ElectronSF[nElectrons];
            nElectrons++;
          }
          
        }
        if(lep3veto && !isData) wnbSelectedEvents_3Lveto += eventweight;
        if(lep3veto1J && !isData) wnbSelectedEvents_3Lveto1J *= eventweight;
          
          
          if (! isData)
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
              mc_isLastCopy[iMC] = mcParticles[iMC]->isLastCopy();
              mc_isPromptFinalState[iMC] = mcParticles[iMC]->isPromptFinalState();
              mc_isHardProcess[iMC] = mcParticles[iMC]->isHardProcess();
              mc_fromHardProcessFinalState[iMC] = mcParticles[iMC]->fromHardProcessFinalState();
            }
          }
      }
      if(eventSelected){
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
        
      }
      
      
      
      if(eventSelected){
        nbSelectedEvents++;
        if(!isData) wnbSelectedEvents += eventweight;
        myTree->Fill();
      }
      
      
      
      
    } // end eventloop
    
    
    
    
    
    
    //	for(int j_eeu = 0; j < 9; j++){       cout << cutstep[j] << endl; }
    sumW = (int) sumWeights;
    nEv = (int) nEvents;
    int  nEvPassed  = (int) nbSelectedEvents;
    globalTree->Fill();
    if(verbose > 0) cout << "end eventloop" << endl;
    
    cout << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl;
    //cout << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
    //cout << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
    
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
    
    tupfile->Close();
    
    delete tupfile;
    
    treeLoader.UnLoadDataset();
  } //End Loop on Datasets
  
  
  
  /////////////
  // Writing //
  /////////////
  
  cout << " - Writing outputs to the files ..." << endl;
  
  
  /*
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
  
  */
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









