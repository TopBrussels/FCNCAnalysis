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



struct HighestCSVBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

//Initializing CSVv2 b-tag WP
float workingpointvalue_Loose = 0.460;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Medium = 0.800;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Tight = 0.935;//working points updated to 2015 BTV-POG recommendations.


TLorentzVector FCNCjetCalculator(std::vector<TRootJet*> nonBJets,std::vector<TRootJet*> BJets, TLorentzVector recoZ ,int verb); 

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
  int nbTrig = 0; 
  int nbBaseline = 0; 
  int nbGPV = 0; 
  int nbSelectedEvents = 0; 
  int nbEvents = 0; 
  double dataLumi = 0; //pb
  bool eee = false; 
  bool eemu = false; 
  bool mumue = false; 
  bool mumumu = true;  
  bool runHLT = true; 
  bool hasMu = false; 
  bool hasEl = false; 
  bool dilep =false; 
  bool singlelep = false;
  bool applyJetCleaning = true; 
  bool fillBtagHisto = false; 
  bool printTrigger = false;
  bool printLeptonSF = false; 
  bool applyJER = false; 
  bool applyJES = false; 
  bool applyNegWeightCorrection = false; 
  bool applyPU = true; 
  bool applyLeptonSF = false;  
  bool btagShape = true; 
  string Channel = ""; 
  string xmlFileName = ""; 
  if(mumumu)
  {
      cout << " --> Using the TriMuon channel <-- " << endl; 
      Channel = "MuMuMu"; 
      xmlFileName = "config/Run2TriLepton_MuMuMu.xml" ; 
      dataLumi = 2700; //pb
      hasMu = true; 
      dilep = true; 
  }
  else
  {
      cerr << " ERROR: no channel specified " << endl; 
      exit(1); 
  }
  
  
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
  const int JES                 =  strtol(argv[argc-7], NULL,10);
  const int JER                 =  strtol(argv[argc-6], NULL,10);
  const int FillBtagHisto	 =  strtol(argv[argc-5], NULL,10);
  string chanName		  = argv[argc-4];
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
  
  // Print information to a textfile
  ofstream infoFile;
  string info_dir = "Information/"+Channel +"/";
  string info_date_dir = info_dir +  dateString +"/";
  mkdir(info_dir.c_str(),0777);
  mkdir(info_date_dir.c_str(),0777); 
  string infoName = info_date_dir + "information"; 
  infoName += "_"+ Channel;
  infoName += "_" + dName;
  infoName += "_" + JobNum; 
  infoName += ".txt"; 
  infoFile.open(infoName.c_str());
  
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
     
  
  infoFile << "---Dataset accepted from command line---" << endl;
  infoFile << "Dataset Name: " << dName << " data? " << isData << endl;
  infoFile << "Dataset Title: " << dTitle << endl;
  infoFile << "Dataset color: " << color << endl;
  infoFile << "Dataset ls: " << ls << endl;
  infoFile << "Dataset lw: " << lw << endl;
  infoFile << "Dataset normf: " << normf << endl;
  infoFile << "Dataset EqLumi: " << EqLumi << endl;
  infoFile << "Dataset xSect: " << xSect << endl;
  infoFile << "Dataset File Name: " << vecfileNames[0] << endl;
  infoFile << "Beginning Event: " << startEvent << endl;
  infoFile << "Ending Event: " << endEvent << endl;
  infoFile << "JobNum: " << JobNum << endl;
  infoFile << "Trigger: " << runHLT << " mu/e/single/di " << hasMu << "/"<< hasEl << "/"<< singlelep << "/" << dilep << endl; 
  infoFile << "Channel: mumumu/mumue/eee/eemu " << mumumu << "/" << mumue << "/" << eee << "/" <<
  eemu << endl; 
  infoFile << "xmlfile: " << xmlFileName.c_str()  << endl; 
  infoFile << "Jetcleaning on? " <<  applyJetCleaning << endl; 
  infoFile << "BtagReweighting  FillHisto? " << fillBtagHisto << endl; 
  infoFile << "JES? " << applyJES << " JER? " << applyJER << endl; 
  infoFile << "Neg Weight correction? " << applyNegWeightCorrection << endl;
  infoFile  << "Lepton SF? " << applyLeptonSF	 << endl; 

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
  Trigger* trigger = new Trigger(hasMu, hasEl, singlelep, dilep);

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
  string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";


  ///////////////////////////////
  //  Set up Output ROOT file  ///
  //////////////////////////////
  stringstream ss;
  ss << JobNum;
  string strJobNum = ss.str();
  string histo_dir = "NtupleMakerOutput/TriLepton_histos_"+ Channel;
  string histo_dir_date = histo_dir+"/TriLepton_histos_" + dateString +"/";
  mkdir(histo_dir.c_str(),0777);
  mkdir(histo_dir_date.c_str(),0777);
  
  string rootFileName (histo_dir_date+"/FCNC_3L_"+Channel+"_"+dName+".root");
  if (strJobNum != "0")
  {
    if(verbose == 0) cout << "strJobNum is " << strJobNum << endl;
    rootFileName = histo_dir_date+"/FCNC_3L_"+Channel+"_"+dName + "_"+strJobNum+".root";
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
  bool TightEl = false; 
  bool MediumEl = true; 
  bool LooseEl = false;  
  // muon
  float mu_pt_cut = 20.; // 40
  float mu_eta_cut = 2.4;
  float mu_iso_cut = 0.15;
  bool TightMu = false; 
  bool MediumMu = false; 
  bool LooseMu = true;  
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
    
  infoFile << "El: pt = "  << el_pt_cut_str << " - eta = " << el_eta_cut_str << " tight/medium/loose " << TightEl << "/" << MediumEl << "/" << LooseEl << endl; 
  infoFile << "Mu: pt = "  << mu_pt_cut_str << " - eta = " << mu_eta_cut_str << " - iso " << mu_iso_cut_str << " tight/medium/loose " << TightMu << "/" << MediumMu<< "/" << LooseMu << endl; 
  infoFile << "Jet: pt = "  << jet_pt_cut_str << " - eta = " << jet_eta_cut_str <<  endl; 
  
  
  

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
/*
    //Muons
    histo1D["MuonPt"]                                        = new TH1F( "MuonPt", "PT_{#mu}", 30, 0, 300);
    histo1D["LeptonPt"]                                        = new TH1F( "LeptonPt", "PT_{lep}", 30, 0, 300);
    histo1D["MuonRelIsolation"]                              = new TH1F( "MuonRelIsolation", "RelIso", 10, 0, .25);
    //Electrons
    histo1D["ElectronRelIsolation"]                          = new TH1F( "ElectronRelIsolation", "RelIso", 10, 0, .25);
    histo1D["ElectronPt"]                                    = new TH1F( "ElectronPt", "PT_{e}", 30, 0, 300);
    //Init Electron Plots

    histo1D["InitElectronPt"]                                = new TH1F( "InitElectronPt", "PT_{e}", 30, 0, 300);
    histo1D["InitElectronEta"]                               = new TH1F( "InitElectronEta", "#eta", 40, -4, 4);
    histo1D["NbOfElectronsInit"]                             = new TH1F( "NbOfElectronsInit", "Nb. of electrons", 10, 0, 10);
    histo1D["InitElectronRelIsolation"]                      = new TH1F( "InitElectronRelIsolation", "RelIso", 10, 0, .25);
    histo1D["InitElectronSuperClusterEta"]                   = new TH1F( "InitElectronSuperClusterEta", "#eta", 10, 0, 2.5);
    histo1D["InitElectrondEtaI"]                             = new TH1F( "InitElectrondEtaI", "#eta", 20, 0, .05);
    histo1D["InitElectrondPhiI"]                             = new TH1F( "InitElectrondPhiI", "#phi", 20, 0, .2);
    histo1D["InitElectronHoverE"]                            = new TH1F( "InitElectronHoverE", "H/E", 10, 0, .15);
    histo1D["InitElectrond0"]                                = new TH1F( "InitElectrond0", "d0", 20, 0, .1);
    histo1D["InitElectrondZ"]                                = new TH1F( "InitElectrondZ", "dZ", 10, 0, .25);
    histo1D["InitElectronEminusP"]                           = new TH1F( "InitElectronEminusP", "1/GeV", 10, 0, .25);
    histo1D["InitElectronConversion"]                        = new TH1F( "InitElectronConversion", "Conversion Pass", 2, 0, 2);
    histo1D["InitElectronMissingHits"]                       = new TH1F( "InitElectronMissingHits", "MissingHits", 10, 0, 10);
    histo1D["InitElectronCutFlow"]                           = new TH1F( "InitElectronCutFlow", "CutNumber", 12, 0, 12);

    //B-tagging discriminators
    histo1D["Bdisc_CSV_jet1"]                             = new TH1F( "Bdisc_CSV_jet1", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2"]                             = new TH1F( "Bdisc_CSV_jet2", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3"]                             = new TH1F( "Bdisc_CSV_jet3", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet1"]                             = new TH1F( "Bdisc_CSV_Bjet1", "CSV b-disc._{bjet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet2"]                             = new TH1F( "Bdisc_CSV_Bjet2", "CSV b-disc._{bjet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet3"]                             = new TH1F( "Bdisc_CSV_Bjet3", "CSV b-disc._{bjet3}", 30, 0, 1);
    //Jets
    histo1D["JetEta"]                                        = new TH1F( "JetEta", "Jet #eta", 40,-4, 4);
    histo1D["NbJets"]                                        = new TH1F( "NbJets", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets"]                                        = new TH1F( "NbCSVLJets", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets"]                                        = new TH1F( "NbCSVMJets", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets"]                                        = new TH1F( "NbCSVTJets", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["1stJetPt"]                                      = new TH1F( "1stJetPt", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt"]                                      = new TH1F( "2ndJetPt", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt"]                                      = new TH1F( "3rdJetPt", "PT_{jet3}", 30, 0, 300);
    histo1D["1stBJetPt"]                                      = new TH1F( "1stBJetPt", "PT_{bjet1}", 30, 0, 300);
    histo1D["2ndBJetPt"]                                      = new TH1F( "2ndBJetPt", "PT_{bjet2}", 30, 0, 300);
    histo1D["3rdBJetPt"]                                      = new TH1F( "3rdBJetPt", "PT_{bjet3}", 30, 0, 300);
    histo1D["HT_SelectedJets"]                               = new TH1F( "HT_SelectedJets", "HT", 30, 0, 1500);
    //MET
    histo1D["MET_preCut"]                                           = new TH1F( "MET_preCut", "MET", 70, 0, 700);
    histo1D["MT_LepMET_preCut"]                                           = new TH1F( "MET_LepMET_preCut", "MT(lep,MET)", 70, 0, 700);
    histo1D["MET"]                                           = new TH1F( "MET", "MET", 70, 0, 700);
    histo1D["MT_LepMET"]                                           = new TH1F( "MT_LepMET", "MT(lep,MET)", 70, 0, 700);

    ///////////////////
    // 2D histograms //
    ///////////////////
    histo2D["NJet_vs_Nbjet"] = new TH2F("NJet_vs_Nbjet","NJet:Nbjet",12,-0.5,11.5, 61, -0.5,11.5);
    histo2D["JetID_vs_pdgID"] = new TH2F("JetID_vs_pdgID","parton pdgID:jet number",12,-0.5,11.5, 61, -30.5,30.5);
*/


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
	infoFile <<"found sample " << daName.c_str() << " with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
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
//	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv"); 
           btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_76X_combToMujets.csv"); 
	   btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets","central");         	
	   if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
	   {
		btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_mujets_central.root",30,999,2.4);
	   }
	   else
	   {
//                btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+"_mujets_central.root",false,30,999,2.4);
                 btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
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
//        LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert.root", "pileup", "pileup");      

        LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data76X_25ns-Run246908-260627Cert.root", "pileup", "pileup");	
//          LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data76X_25ns-Run246908-260627Cert.root", "pileup", "pileup");  

       //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)

        MuonSFWeight* muonSFWeightID_T = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);
        MuonSFWeight* muonSFWeightID_M = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true,  printLeptonSF, printLeptonSF);
        MuonSFWeight* muonSFWeightID_L = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, printLeptonSF, printLeptonSF);
        MuonSFWeight* muonSFWeightIso_TT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);  // Tight RelIso, Tight ID
        MuonSFWeight* muonSFWeightIso_TM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Tight RelIso, Medium ID
        MuonSFWeight* muonSFWeightIso_LT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Tight ID
        MuonSFWeight* muonSFWeightIso_LM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Medium ID
//        double weightMuonHLTv2, weightMuonHLTv3 ; // for run C should also something like this be done
//        MuonSFWeight *muonSFWeightTrigHLTv4p2 = new MuonSFWeight(CaliPath+"LeptonSF/"+"SingleMuonTrigger_Z_RunCD_Reco76X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio", true, false, false);
//        MuonSFWeight *muonSFWeightTrigHLTv4p3 = new MuonSFWeight(CaliPath+"LeptonSF/"+"SingleMuonTrigger_Z_RunCD_Reco76X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio", true, false, false);
  


            
        string electronFile= "Elec_SF_TopEA.root";
        ElectronSFWeight* electronSFWeight = new ElectronSFWeight (CaliPath+"LeptonSF/"+electronFile,"GlobalSF", true,printLeptonSF, printLeptonSF); // (... , ... , debug, print warning)  

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

        string channel_dir = "NtupleMakerOutput/Ntuples_"+Channel;
        string date_dir = channel_dir+"/Ntuples_" + dateString +"/";
        mkdir(channel_dir.c_str(),0777);
        mkdir(date_dir.c_str(),0777);

        
        string Ntupname = date_dir +"FCNC_3L_" +Channel + "_" + dName + "_"+  strJobNum + ".root";

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
       Int_t evt_num;
       Int_t lumi_num;
       Int_t nvtx;
       Int_t npu;
       Int_t PassedMETFilter; 
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
	Double_t sf_electron[10];

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
	Double_t sf_muon[10];
        Int_t charge_muon[10];
  
        //variable for jets 
        Int_t nJets;
	Int_t nJets_CSVL; 
	Int_t nJets_CSVM; 
	Int_t nJets_CSVT;
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Int_t charge_jet[20];
        Double_t bdisc_jet[20];
        Double_t cdiscCvsL_jet[20]; 
	Double_t cdiscCvsB_jet[20]; 


        // variables for Zboson
        Double_t Zboson_M; 
/*	Double_t Zboson_Px; 
        Double_t Zboson_Py;
        Double_t Zboson_Pz;
	Double_t Zboson_Energy;
*/
        // met 
        Double_t met_Pt; 
	Double_t met_Phi; 
	Double_t met_Eta; 

	Double_t mWt; 
        Double_t FCNCtop_M; 
        Double_t SMtop_M;
        // global data set variables
	Int_t nofEventsHLTv2; 
	Int_t nofEventsHLTv3; 
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

       // event related variables
       myTree->Branch("nloWeight",&nloWeight,"nloWeight/D"); 
       myTree->Branch("run_num",&run_num,"run_num/I");
       myTree->Branch("evt_num",&evt_num,"evt_num/I");
       myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
       myTree->Branch("nvtx",&nvtx,"nvtx/I");
       myTree->Branch("npu",&npu,"npu/I");
       myTree->Branch("puSF",&puSF,"puSF/D");  
       myTree->Branch("btagSF",&btagSF,"btagSF/D");         
       myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
       myTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I"); 

       baselineTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I");
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
       myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
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
       baselineTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
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
       myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
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
       baselineTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
       baselineTree->Branch("pt_muon_1",&pt_muon_1,"pt_muon_1/D");
       baselineTree->Branch("pt_muon_2",&pt_muon_2,"pt_muon_2/D");
       baselineTree->Branch("pt_muon_3",&pt_muon_3,"pt_muon_3/D");

       // jets
       myTree->Branch("nJets",&nJets,"nJets/I");
       myTree->Branch("nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
       myTree->Branch("nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
       myTree->Branch("nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
       myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
       myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
       myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
       myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
       myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");	    
       myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
       myTree->Branch("cdiscCvsL_jet",cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
       myTree->Branch("cdiscCvsB_jet",cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
       myTree->Branch("pt_jet_1",&pt_jet_1,"pt_jet_1/D");
       myTree->Branch("pt_jet_2",&pt_jet_2,"pt_jet_2/D");
       myTree->Branch("pt_jet_3",&pt_jet_3,"pt_jet_3/D");
       
       baselineTree->Branch("nJets",&nJets,"nJets/I");
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
       myTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/D"); 
       baselineTree->Branch("SMtop_M",&SMtop_M, "SMtop_M/D");
 /*      myTree->Branch("Zboson_Px",&Zboson_Px,"Zboson_Px/D"); 
       myTree->Branch("Zboson_Py",&Zboson_Py,"Zboson_Py/D");
       myTree->Branch("Zboson_Pz",&Zboson_Pz,"Zboson_Pz/D");
       myTree->Branch("Zboson_Energy",&Zboson_Energy,"Zboson_Energy/D");
*/

        // met 
       myTree->Branch("met_Pt", &met_Pt, "met_Pt/D"); 
       myTree->Branch("met_Eta", &met_Eta,"met_Eta/D"); 
       myTree->Branch("met_Phi", &met_Phi, "met_Phi/D"); 
     
       baselineTree->Branch("met_Pt", &met_Pt, "met_Pt/D"); 
       baselineTree->Branch("met_Eta", &met_Eta,"met_Eta/D"); 
       baselineTree->Branch("met_Phi", &met_Phi, "met_Phi/D"); 

       

       /////////////////////////
       //// Corrections/trigger ///
       ///////////////////////////

       /// book triggers
       if (runHLT) { trigger->bookTriggers(isData);}

       



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
	vector < TRootJet* >      init_fatjets;
        vector < TRootJet* >      init_jets_corrected;
        vector < TRootGenJet* >   genjets;
        vector < TRootMET* >      mets;
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootJet*>      selectedCSVLBJets;
        vector<TRootJet*>      selectedCSVMBJets;
        vector<TRootJet*>      selectedCSVTBJets;
	vector<TRootJet*>      selectedCSVLLJets;
        vector<TRootJet*>      selectedCSVMLJets;
        vector<TRootJet*>      selectedCSVTLJets;
        vector<TRootMCParticle*> mcParticles;
        vector <TRootJet*>     selectednonCSVLJets; 

       TLorentzVector Zboson;
       TLorentzVector Zlep0;
       TLorentzVector Zlep1;
       TLorentzVector Wlep;
       TLorentzVector SMbjet;
       TLorentzVector cjet;
        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
	nbEvents = 0; 
        nofEventsHLTv2 = 0; 
        nofEventsHLTv3 = 0;
        nofPosWeights = 0; 
        nofNegWeights = 0; 
        float eventweight = 1;
        int nbEvents_0 = 0; 
        int nbEvents_1 = 0;
        int nbEvents_2 = 0;
        int nbEvents_3 = 0;
        int nbEvents_4 = 0;
        int nbEvents_5 = 0;
        int nbEvents_6 = 0;
	int nbEvents_7 = 0; 
	int nbEvents_8 = 0; 
	int nbEvents_9 = 0; 
        bool debug = false; 

        bool   passedMET = false;
        bool   HBHEnoise = false;
        bool   HBHEIso = false;
        bool   CSCTight = false;
        bool   EcalDead = false;
        bool    eeBad = false; 
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
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
	    genjets.clear();
	    if(!isData) genjets = treeLoader.LoadGenJet(ievt,false);  //needed for JER
	    
	    
	    if(verbose == 0)
	    {
	    	cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
	        cout <<"Number of Muons Loaded: " << init_muons.size() <<endl; 
	        cout << "Number of Jets Loaded: " << init_jets_corrected.size() << endl;  
	    }

  
             //  take the event          
            datasets[d]->eventTree()->LoadTree(ievt);
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            int currentRun = event->runId();
	    run_num = event->runId(); 
	    evt_num = event->eventId();
/*          // to be applied from 76X v1 in our ttp
            HBHEnoise = event->getHBHENoiseFilter();
            HBHEIso = event->getHBHENoiseIsoFilter();
            CSCTight = event->getCSCTightHalo2015Filter(); 
            EcalDead = event->getEcalDeadCellTriggerPrimitiveFilter();
            eeBad = event->getEEBadScFilter();  
*/
            
//            cout << "eeBadSc " << eeBadSc << endl;
	    lumi_num=event->lumiBlockId(); 
	    nvtx = vertex.size();
	    npu = (int) event->nTruePU(); 

/*           if(isData) // run C should be added as third counter
           {
                 if(currentRun >= 256630 && currentRun <= 257819 )  // run nbrs need to be checked
        	{
        	  nofEventsHLTv2++;
        	}
           	else
        	{
        	  nofEventsHLTv3++;
        	}

          }

*/
           /////////////////////////////////////
           //  fix negative weights for amc@nlo/// 
           /////////////////////////////////////
	   double hasNegWeight = false; 
           double mc_baseweight = 1; 
	   if(!isData && (event->getWeight(1001) != -9999.))
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
	   if( !isData && (event->getWeight(1) != -9999. )) 
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
	   if(!isData)
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
            bool filechanged = false; 
            bool runchanged = false; 
            
            if(runHLT)
            {
               trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
               trigged = trigger->checkIfFired();

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
            if(dName.find("NP")!=string::npos) trigged = true; 
 
            if(verbose==0) cout << "Apply trigger? " << runHLT << " trigged? " << trigged << endl; 

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
            Run2Selection selection(init_jets,init_fatjets, init_muons, init_electrons, mets,event->fixedGridRhoFastjetAll());
	    selectedJets.clear(); 
	    selectedJets  = selection.GetSelectedJets(jet_pt_cut,jet_eta_cut, true, "Tight"); 
	    selectedMuons.clear();
            if(TightMu)  selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Tight", "Spring15"); 
            if(MediumMu)  selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Medium", "Spring15");
            if(LooseMu)  selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Loose", "Spring15");
	    // pt, eta, iso // run normally
	    selectedElectrons.clear();
	    if(TightEl) selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Tight","Spring15_25ns",true);// pt, eta
            if(MediumEl) selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Medium","Spring15_25ns",true);// pt, eta
            if(LooseEl) selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Loose","Spring15_25ns",true);// pt, eta
            /// For MC Information
            mcParticles.clear();
            treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
            sort(mcParticles.begin(),mcParticles.end(),HighestPt());
	    // void TTreeLoader::LoadMCEvent(int, TopTree::TRootNPGenEvent*, std::vector<TopTree::TRootMCParticle*>&, bool) 
	    if (verbose == 0) cout <<"Number of Muons, Electrons, Jets  ===>  " << endl << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< selectedJets.size()   << endl;

            
            ////////////////////////////////////////////////
            // Pre cut operations
            ////////////////////////////////////////////////
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
	    // Met filters if(HBHEnoise && HBHEIso && CSCTight && EcalDead && eeBad && isGoodPV) passedMET = true;
	     passedMET = true; 
	     PassedMETFilter = passedMET; 
	    
	    if(applyJetCleaning){
             if(verbose == 0) cout << "Applying jet cleaning " << endl; 
             int OrigSize = selectedJets.size();
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
              if(verbose == 0){
              if( OrigSize != selectedJets.size()) cout << "--> original = " << OrigSize  << " after cleaning = " << selectedJets.size() << endl;
              else cout << "--> no change" << endl; 
              }
           }


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
 		btagWeight =  btwt->getMCEventWeight(selectedJets);

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
            if(hasNegWeight && applyNegWeightCorrection && !isData) eventweight *= -1.; 
	    histo1D["init_nPVs_before"]->Fill(vertex.size(), eventweight); 
            if(applyPU && !isData)  eventweight *= PUweight;
	    histo1D["init_nPVs_after"]->Fill(vertex.size(), eventweight);

            //////////////////////////////////////////////////////
            // Applying baseline selection
            //////////////////////////////////////////////////////
	    nbEvents++;
	    eventweight = 1.;  
	    if(!isGoodPV) continue;
            nbGPV++; 
            if(verbose == 0) cout << "good pv" << endl; 
            if(!passedMET) continue;  
	    if(!trigged) continue; 
            nbTrig++; 
            if(verbose == 0 ) cout << "trigger" << endl; 
            histo1D["cutFlow"]->Fill(0., eventweight);
	    nCuts++;
            nbEvents_0++; 
//            cout << " after " << nCuts << " " << nbEvents_0 << endl;
	    if(mumumu &&  selectedMuons.size() < 2) continue; 
	    if(mumue &&  selectedMuons.size() < 2) continue; 
	    if(eemu  && selectedElectrons.size() < 2) continue; 
	    if(eee  &&  selectedElectrons.size() < 2) continue;
	    if(verbose == 0 ) cout << "baseline" << endl;
	    histo1D["cutFlow"]->Fill(1., eventweight);   
            nCuts++;
            nbEvents_1++; 
	    
            nElectrons=0;
            for (Int_t selel =0; selel < selectedElectrons.size() ; selel++ )
	    {
	      
              pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
	      phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
	      eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
	      eta_superCluster_electron[nElectrons]=selectedElectrons[selel]->superClusterEta();
	      E_electron[nElectrons]=selectedElectrons[selel]->E();
	      d0_electron[nElectrons]=selectedElectrons[selel]->d0();
	      d0BeamSpot_electron[nElectrons]=selectedElectrons[selel]->d0BeamSpot();
	      chargedHadronIso_electron[nElectrons]=selectedElectrons[selel]->chargedHadronIso(3);
	      neutralHadronIso_electron[nElectrons]=selectedElectrons[selel]->neutralHadronIso(3);
	      photonIso_electron[nElectrons]=selectedElectrons[selel]->photonIso(3);
	      pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
	      charge_electron[nElectrons]=selectedElectrons[selel]->charge();
	      sigmaIEtaIEta_electron[nElectrons]=selectedElectrons[selel]->sigmaIEtaIEta();
	      deltaEtaIn_electron[nElectrons]=selectedElectrons[selel]->deltaEtaIn();
	      deltaPhiIn_electron[nElectrons]=selectedElectrons[selel]->deltaPhiIn();
	      hadronicOverEm_electron[nElectrons]=selectedElectrons[selel]->hadronicOverEm();
	      missingHits_electron[nElectrons]=selectedElectrons[selel]->missingHits();
	      passConversion_electron[nElectrons]=selectedElectrons[selel]->passConversion();
	      isEBEEGap[nElectrons]=selectedElectrons[selel]->isEBEEGap();
	      if(!isData) sf_electron[nElectrons]=electronSFWeight->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0); 
	      else sf_electron[nElectrons] = 1.; 
	      if(!isData) ElectronSF[nElectrons] = electronSFWeight->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
              else ElectronSF[nElectrons] = 1.; 
	     
              nElectrons++;
            }
	    if(selectedElectrons.size()>0) pt_electron_1 = selectedElectrons[0]->Pt(); 
	    if(selectedElectrons.size()>1) pt_electron_2 = selectedElectrons[1]->Pt(); 
	    if(selectedElectrons.size()>2) pt_electron_3 = selectedElectrons[2]->Pt(); 

            //////////////////////
            // Muon Based Plots //
            //////////////////////
            nMuons = 0; 
            for (Int_t selmu =0; selmu < selectedMuons.size() ; selmu++ )
            {
              
              pt_muon[nMuons]=selectedMuons[selmu]->Pt();
	      phi_muon[nMuons]=selectedMuons[selmu]->Phi();
	      eta_muon[nMuons]=selectedMuons[selmu]->Eta();
	      E_muon[nMuons]=selectedMuons[selmu]->E();
	      d0_muon[nMuons]=selectedMuons[selmu]->d0();
	      d0BeamSpot_muon[nMuons]=selectedMuons[selmu]->d0BeamSpot();
	      chargedHadronIso_muon[nMuons]=selectedMuons[selmu]->chargedHadronIso(4);
	      neutralHadronIso_muon[nMuons]=selectedMuons[selmu]->neutralHadronIso(4);
	      photonIso_muon[nMuons]=selectedMuons[selmu]->photonIso(4);
              pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
	      charge_muon[nMuons]=selectedMuons[selmu]->charge();
	      if(!isData) sf_muon[nMuons]= muonSFWeightIso_TT->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID_T->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
	      else sf_muon[nMuons] = 1.;
	      if(!isData)
	      {
		MuonIDSF[nMuons] = muonSFWeightID_T->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
		MuonIsoSF[nMuons] =  muonSFWeightIso_TT->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0); 
//		MuonTrigSFv2[nMuons] = muonSFWeightTrigHLTv4p2->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0); 
//		MuonTrigSFv3[nMuons] = muonSFWeightTrigHLTv4p3->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0); 
	      }
	      else
	      {
		MuonIDSF[nMuons] = 1.; 
		MuonIsoSF[nMuons] = 1.; 
//		MuonTrigSFv2[nMuons] = 1.;
//		MuonTrigSFv3[nMuons] = 1.; 
              }
              nMuons++;
            }
            if(selectedMuons.size()>0) pt_muon_1 = selectedMuons[0]->Pt(); 
	    if(selectedMuons.size()>1) pt_muon_2 = selectedMuons[1]->Pt(); 
	    if(selectedMuons.size()>2) pt_muon_3 = selectedMuons[2]->Pt();
            nLeptons = nMuons + nElectrons; 
	    ///////////////////////
	    //   Jet based plots //
	    //////////////////////
	    nJets = 0; 
            for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
            {
                 
                pt_jet[nJets]=selectedJets[seljet]->Pt(); 
                phi_jet[nJets]=selectedJets[seljet]->Phi();
                eta_jet[nJets]=selectedJets[seljet]->Eta();
                E_jet[nJets]=selectedJets[seljet]->E();
                charge_jet[nJets]=selectedJets[seljet]->charge();
                bdisc_jet[nJets]=selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() ;
                cdiscCvsB_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsBJetTags() ;
                cdiscCvsL_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsLJetTags() ;
                nJets++;

            }
	    if(selectedJets.size()>0) pt_jet_1 = selectedJets[0]->Pt(); 
	    if(selectedJets.size()>1) pt_jet_2 = selectedJets[1]->Pt(); 
	    if(selectedJets.size()>2) pt_jet_3 = selectedJets[2]->Pt();
            nJets_CSVT =  selectedCSVTBJets.size(); 
	    nJets_CSVM =  selectedCSVMBJets.size();
            nJets_CSVL =  selectedCSVLBJets.size();
	    double met_px = mets[0]->Px();
	    double met_py = mets[0]->Py();
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
	    met_Phi = mets[0]->Phi(); 
	    met_Eta = mets[0]->Eta();
	    puSF = PUweight;
	    btagSF = btagWeight;  
	    
            histo1D["cutFlow"]->Fill(2., eventweight); 
            nCuts++;
            nbEvents_2++;  
//            cout << " after " << nCuts << " " << nbEvents_2 << endl;
            if(selectedJets.size() < 2) continue; 
	    histo1D["cutFlow"]->Fill(3., eventweight);
            nCuts++;
            nbEvents_3++;  
//            cout << " after " << nCuts << " " << nbEvents_3 << endl;
//            if(selectedCSVLBJets.size() < 1) continue; 
	    histo1D["cutFlow"]->Fill(4., eventweight);
            nCuts++;
            nbEvents_4++;  
//            cout << " after " << nCuts << " " << nbEvents_4 << endl;
            if(selectedMuons.size() + selectedElectrons.size() <3)    baselineTree->Fill();
             if(selectedMuons.size() +selectedElectrons.size() <3)   nbBaseline++;
	    //check flavour
            histo1D["nbMuons"]->Fill(selectedMuons.size(), eventweight);
            histo1D["nbElectrons"]->Fill(selectedElectrons.size(), eventweight);
            histo1D["nbJets"]->Fill(selectedJets.size(), eventweight);

            if(selectedElectrons.size() + selectedMuons.size() <3) continue;
            histo1D["cutFlow"]->Fill(5., eventweight);
            nCuts++;
            nbEvents_5++;
 
            if(mumumu && selectedMuons.size() <3) continue;
            histo1D["cutFlow"]->Fill(6., eventweight);
            nCuts++;
            nbEvents_6++;  
//            cout << " after " << nCuts << " " << nbEvents_5 << endl;            
            if(selectedCSVLBJets.size() > 0) continue;
	    Zlep0.Clear(); 
	    Zlep1.Clear();
            Wlep.Clear();  

	    // check sign
	    bool OS = false; 
	    if(eemu && (selectedElectrons[0]->charge() == selectedElectrons[1]->charge())) continue;
            if(mumue && (selectedMuons[0]->charge() == selectedMuons[1]->charge())) continue;
	    if(mumumu)
	    {
	       if(selectedMuons[0]->charge() != selectedMuons[1]->charge()){
 			 OS = true;
 			 Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
  			 Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
			 if(selectedMuons.size() > 2) Wlep.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
  			 else Wlep.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
               } 
	       if(selectedMuons.size() > 2) {
		if(selectedMuons[2]->charge() != selectedMuons[1]->charge()){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
              		 Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
		}
                 else if(selectedMuons[0]->charge() != selectedMuons[2]->charge()){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
			 Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
		}
	       }
            }
            if(eee)
            {
               if(selectedElectrons[0]->charge() != selectedElectrons[1]->charge()){
 			 OS = true;
 			 Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
  			 Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
			 if(selectedElectrons.size() > 2) Wlep.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
  			 else Wlep.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
               } 
	       if(selectedElectrons.size() > 2) {
		if(selectedElectrons[2]->charge() != selectedElectrons[1]->charge()){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
              		 Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
		}
                 else if(selectedElectrons[0]->charge() != selectedElectrons[2]->charge()){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
			 Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
		}
	       }   
            }
	

             eventSelected = true; 
            if(mumue) Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
            if(mumue) Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
            if(mumue) Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
           
            if(eemu)  Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
            if(eemu) Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
            if(eemu) Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
            if(mumumu && !OS) continue;
            if(eee && !OS) continue;
            histo1D["cutFlow"]->Fill(7., eventweight); 
            nCuts++;
            nbEvents_7++;  
            eventSelected = true;


            //Make event variables
	    Zboson.Clear(); 
            Zboson.SetPxPyPzE(( Zlep0 + Zlep1).Px() ,( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Energy()) ;
            Zboson_M = (Zlep0+Zlep1).M();

 	    SMbjet.Clear(); 
	    SMbjet.SetPxPyPzE(selectedCSVLBJets[0]->Px(),selectedCSVLBJets[0]->Py(),selectedCSVLBJets[0]->Pz(),selectedCSVLBJets[0]->Energy());

	    cjet.Clear(); 
	    cjet = FCNCjetCalculator(selectedCSVLLJets,selectedCSVLBJets, Zboson ,3);

            FCNCtop_M = (Zboson+cjet).M();
	    SMtop_M = (Wlep+SMbjet).M();  

            mWt = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_px)*(Wlep.Px() + met_px) - (Wlep.Py() + met_py)* (Wlep.Py() + met_py)); 

	    if(fabs((Zlep0+Zlep1).M() - 90.0 ) > 15) continue; 
	    histo1D["cutFlow"]->Fill(8., eventweight);
            nCuts++;
	    nbEvents_8++; 
            if(fabs((Wlep+SMbjet).M() - 173.0) > 35 ) continue; 
            histo1D["cutFlow"]->Fill(9., eventweight);
            nCuts++;
            nbEvents_9++;	    
             //////////////////////////////////////
	    //  DO STUFF WITH SELECTED EVENTS ////
	    ////////////////////////////////////// 
	    if(!eventSelected) continue;
	    nbSelectedEvents++; 
	    myTree->Fill(); 
	       
	    
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
	if(debug)	for(int j = 0; j < 7; j++){       cout << cutstep[j] << endl; }
        sumW = (int) sumWeights; 
        nEv = (int) nEvents; 
	globalTree->Fill(); 
        if(verbose == 0) cout << "end eventloop" << endl; 
	infoFile << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl;
        infoFile << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
        infoFile << nbBaseline << " baseline events out of trigged " << nbTrig <<  " selected " << endl;
        infoFile << setprecision(2) << ((double)nbGPV/(double)nbEvents)*100 << " % of the initial events stay after Good PV" << endl;
        infoFile << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
        infoFile << setprecision(2) << ((double)nbTrig/(double)nbGPV)*100 << " % of the GPV  events stay after Trigger" << endl;
        cout << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl; 
        cout << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
        cout << nbBaseline << " baseline events out of trigged " << nbTrig <<  " selected " << endl;
        cout << setprecision(2) << ((double)nbGPV/(double)nbEvents)*100 << " % of the initial events stay after Good PV" << endl; 
	cout << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
        cout << setprecision(2) << ((double)nbTrig/(double)nbGPV)*100 << " % of the GPV  events stay after Trigger" << endl;
	if (! isData  ) 
        {
      		cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
      		cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
      		cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents) << endl;
      
                // Determine scale factor due to negative weights
            	nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
                cout << "This corresponds to an event scale factor of " << nloSF  << endl; 
        }
	infoFile.close();
	tupfile->Write();   
    	tupfile->Close();
        delete tupfile;
        if(!isData && !btagShape) delete btwt; 
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




TLorentzVector FCNCjetCalculator(std::vector<TRootJet*> nonBJets,std::vector<TRootJet*> BJets, TLorentzVector recoZ ,int verb)
{
    TLorentzVector FCNCjet; 
    FCNCjet.Clear(); 


     double TempMinMass = 100000.00;
     double TopMass = 172.9;
     TLorentzVector Jetcandidate;
     int NbInColl = -1;
     if(nonBJets.size() != 0){

       for(unsigned int iJ = 0; iJ < nonBJets.size(); iJ++)
       {
     	  TLorentzVector Jet;
     	  Jet.SetPxPyPzE(nonBJets[iJ]->Px(),nonBJets[iJ]->Py(),nonBJets[iJ]->Pz(),nonBJets[iJ]->Energy());

     	  if(fabs((recoZ+Jet).M() - TopMass) < TempMinMass)
     	  {
     	    TempMinMass = fabs((recoZ+Jet).M() - TopMass);
     	    Jetcandidate.SetPxPyPzE(Jet.Px(), Jet.Py(), Jet.Pz(), Jet.E());
     	    NbInColl = iJ;

     	  }


       }
       FCNCjet.SetPxPyPzE(nonBJets[NbInColl]->Px(),nonBJets[NbInColl]->Py(),nonBJets[NbInColl]->Pz(),nonBJets[NbInColl]->Energy());
     }
     else {
       for(unsigned int iJ = 1; iJ < BJets.size(); iJ++)
       {
     	  TLorentzVector Jet;
     	  Jet.SetPxPyPzE(BJets[iJ]->Px(),BJets[iJ]->Py(),BJets[iJ]->Pz(),BJets[iJ]->Energy());

     	  if(fabs((recoZ+Jet).M() - TopMass) < TempMinMass)
     	  {
     	    TempMinMass = fabs((recoZ+Jet).M() - TopMass);
     	    Jetcandidate.SetPxPyPzE(Jet.Px(), Jet.Py(), Jet.Pz(), Jet.E());
     	    NbInColl = iJ;

     	  }

       }

       FCNCjet.SetPxPyPzE(BJets[NbInColl]->Px(),BJets[NbInColl]->Py(),BJets[NbInColl]->Pz(),BJets[NbInColl]->Energy());
     }


    return FCNCjet;
}
