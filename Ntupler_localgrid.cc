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
float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.



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
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
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
  string Channel = ""; 
  string xmlFileName = ""; 
  if(mumumu)
  {
      cout << " --> Using the TriMuon channel <-- " << endl; 
      Channel = "MuMuMu"; 
      xmlFileName = "config/Run2TriLepton_MuMuMu.xml" ; 
      dataLumi = 1200; //pb
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
  if(argc < 15)
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
  const int FillBtagHisto	 =  strtol(argv[argc-5], NULL,10);
  string chanName		  = argv[argc-4];
  const int JobNum		  = strtol(argv[argc-3], NULL, 10);
  const int startEvent  	  = strtol(argv[argc-2], NULL, 10);
  const int endEvent		  = strtol(argv[argc-1], NULL, 10);

  fillBtagHisto = FillBtagHisto; 
  // all the files are stored from arg 11 to argc-2
  vector<string> vecfileNames;
  for(int args = 11; args < argc-5; args++)
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
  string infoName =  "Information/information"; 
  infoName += "_"+ Channel;
  infoName += "_" + dName;
  infoName += "_" + JobNum; 
  infoName += "_" + dateString;
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

  // for pu
  LumiReWeighting LumiWeights;
 

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
  
  string rootFileName (histo_dir_date+"/FCNC_3L_"+Channel+".root");
  if (strJobNum != "0")
  {
    if(verbose == 0) cout << "strJobNum is " << strJobNum << endl;
    rootFileName = histo_dir_date+"/FCNC_3L_"+Channel+"_"+strJobNum+".root";
  }
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
  float el_eta_cut = 2.4;
  bool TightEl = true; 
  bool MediumEl = false; 
  bool LooseEl = false; 
  // muon
  float mu_pt_cut = 20.; // 40
  float mu_eta_cut = 2.4;
  float mu_iso_cut = 0.15;
  bool TightMu = true; 
  bool MediumMu = false; 
  bool LooseMu = false;  
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
/*
    histo1D["NbOfVertices"]                                  = new TH1F("NbOfVertices", "Nb. of vertices", 60, 0, 60);
    histo1D["cutFlow"]                                  = new TH1F( "cutFlow", "cutFlow", 15, -0.5, 14.5);
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


    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    
	treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset	

        string daName = datasets[d]->Name();
        float normfactor = datasets[d]->NormFactor();
	cout <<"found sample " << daName.c_str() << " with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
	infoFile <<"found sample " << daName.c_str() << " with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
        if(daName.find("Data")!=string::npos || daName.find("data")!=string::npos || daName.find("DATA")!=string::npos){
   	   isData = true;
         }	

        /////////////////////////////////////////
        ///    Calibrations                  ///
        ////////////////////////////////////////
        string CaliPath = "../TopTreeAnalysisBase/Calibrations/"; 
        string BCaliPath = CaliPath + "BTagging/CSVv2_13TeV_25ns_combToMujets.csv";
        if(!isData)
	{
           // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv"); 
	   btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_MEDIUM, "mujets","central");         	
	   if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
	   {
		btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+daName+ "_" + strJobNum +"_mujets_central.root",30,999,2.4);
	   }
	   else
	   {
                btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);

	   }


        }

        LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root", "pileup60", "pileup");      	


       //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)

        MuonSFWeight* muonSFWeightID_T = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);
        MuonSFWeight* muonSFWeightID_M = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio",true,  printLeptonSF, printLeptonSF);
        MuonSFWeight* muonSFWeightID_L = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, printLeptonSF, printLeptonSF);
        MuonSFWeight* muonSFWeightIso_TT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio",true, printLeptonSF,printLeptonSF);  // Tight RelIso, Tight ID
        MuonSFWeight* muonSFWeightIso_TM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Tight RelIso, Medium ID
        MuonSFWeight* muonSFWeightIso_LT = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Tight ID
        MuonSFWeight* muonSFWeightIso_LM = new MuonSFWeight(CaliPath+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true,printLeptonSF, printLeptonSF);  // Loose RelIso, Medium ID
            
        string electronFile= "Elec_SF_TopEA.root";
        ElectronSFWeight* electronSFWeight = new ElectronSFWeight (CaliPath+"LeptonSF/"+electronFile,"GlobalSF", true,printLeptonSF, printLeptonSF); // (... , ... , debug, print warning)  
      

        ////////////////////////////////////////////////////////////
        // Setup Date string and nTuple for output  
        ///////////////////////////////////////////////////////////

        string channel_dir = "NtupleMakerOutput/Ntuples_"+Channel;
        string date_dir = channel_dir+"/Ntuples_" + dateString +"/";
        mkdir(channel_dir.c_str(),0777);
        mkdir(date_dir.c_str(),0777);

        
        string Ntupname = date_dir +"FCNC_3L_" +Channel + "_" + strJobNum + ".root";

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
	tupfile->cd();
	TTree* myTree = new TTree("tree","tree");
                    
	///////////////////////////
       /// output tree
       ///////////////////////////
       // event related variables
       Int_t run_num;
       Int_t evt_num;
       Int_t lumi_num;
       Int_t nvtx;
       Int_t npu;
       Double_t puSF;
       // event related variables
       myTree->Branch("run_num",&run_num,"run_num/I");
       myTree->Branch("evt_num",&evt_num,"evt_num/I");
       myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
       myTree->Branch("nvtx",&nvtx,"nvtx/I");
       myTree->Branch("npu",&npu,"npu/I");
       myTree->Branch("puSF",&puSF,"puSF/D");  

       //////////////////////////
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



        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
	nbEvents = 0; 
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
           if(verbose == 0 ) cout << "new event " << ievt << endl; 
           double ievt_d = ievt;
	    
	    bool debug = false; 
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
	    if(!isData) genjets = treeLoader.LoadGenJet(ievt,false);
	    
	    
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
	    lumi_num=event->lumiBlockId(); 
	    nvtx = vertex.size();
	    npu = (int) event->nTruePU(); 

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


            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Declare selection instance
            Run2Selection selection(init_jets,init_fatjets, init_muons, init_electrons, mets);
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
	    for(unsigned int iJ = 0; iJ < selectedJets.size(); iJ++)
	    {
		if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose) selectedCSVLBJets.push_back(selectedJets[iJ]); 
	        else selectedCSVLLJets.push_back(selectedJets[iJ]);
                if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium) selectedCSVMBJets.push_back(selectedJets[iJ]);
                else selectedCSVMLJets.push_back(selectedJets[iJ]);
                if(selectedJets[iJ]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight) selectedCSVTBJets.push_back(selectedJets[iJ]);
                else selectedCSVTLJets.push_back(selectedJets[iJ]);

	    }

	   ////////////////////////////////////
	   //   Event Weights               ///
	   ///////////////////////////////////
	   float btagWeight  =  1;
           if( fillBtagHisto && !isData)
           {
		btwt->FillMCEfficiencyHistos(selectedJets);

	   } 
           else if( !fillBtagHisto && !isData)
	   {
 		btagWeight =  btwt->getMCEventWeight(selectedJets);

            }


           float PUweight = 1; 
	   if(!isData)
           {
                PUweight = LumiWeights.ITweight((int)event->nTruePU()); 


           }

            //////////////////////////////////////////////////////
            // Applying baseline selection
            //////////////////////////////////////////////////////
	    nbEvents++; 
	    if(!isGoodPV) continue;
            nbGPV++; 
            if(verbose == 0) cout << "good pv" << endl;  
	    if(!trigged) continue; 
            nbTrig++; 
            if(verbose == 0 ) cout << "trigger" << endl; 
	    if(mumumu && !hasMu &&  selectedMuons.size() < 2) continue; 
	    if(mumue && hasMu && !hasEl && selectedMuons.size() < 2) continue; 
	    if(eemu && hasEl && !hasMu &&selectedElectrons.size() < 2) continue; 
	    if(mumue && hasMu && hasEl && (selectedMuons.size() < 1 || selectedElectrons.size() < 1) ) continue; 
	    if(eemu && hasEl && !hasMu && (selectedElectrons.size() < 1 || selectedMuons.size() <1)) continue; 
	    if(eee && hasEl &&  selectedElectrons.size() < 2) continue;
            nbBaseline++;  
	    if(verbose == 0 ) cout << "baseline" << endl;
	    eventSelected = true; 
	    
	    
	    
	   
	    if(eventSelected) 
	    {

           	float MUweight = 1;
           	if(!isData)
           	{
                	for(unsigned int iMu =0 ; iMu < selectedMuons.size(); iMu++)
                	{
             		       	if(TightMu) MUweight *= muonSFWeightIso_TT->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0)* muonSFWeightID_T->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0);
                    		if(MediumMu) MUweight *= muonSFWeightIso_TM->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0)* muonSFWeightID_M->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0); // needs to be checked                     
                    		if(LooseMu) MUweight *= muonSFWeightIso_LM->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0)* muonSFWeightID_L->at(selectedMuons[iMu]->Eta(), selectedMuons[iMu]->Pt(), 0); // needs to be checked
                	}
           	}
           	float ELweight = 1;
           	if(!isData)
           	{
                	for(unsigned int iEl = 0; iEl < selectedElectrons.size(); iEl++)
                	{
                    	    ELweight *= electronSFWeight->at(selectedElectrons[iEl]->Eta(),selectedElectrons[iEl]->Pt(),0);

                	}	
           	}
	       nbSelectedEvents++; 
	       myTree->Fill(); 
	       
	    }
	} // end eventloop
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
	infoFile.close(); 
	 tupfile->Write();   
    	tupfile->Close();
        delete tupfile;
        delete btwt; 
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

  

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;


/*

  for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
    cout << "1D Plot: " << it->first << endl;
   // TCanvas ctemp = 
    
    TH1F *temp = it->second;
    temp->Draw();  
  }
  for (map<string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {
     cout << "2D Plot: " << it->first << endl;
   
     TH2F *temp = it->second;
     temp->Draw();
  }


*/


    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
