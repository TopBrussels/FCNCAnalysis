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
double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy);
float EffectiveAreaRho(TRootElectron *el, float _rho) ;
float EffectiveArea(TRootElectron *el) ; 
float relPfIsoEl(TRootElectron *el, float _rho);
float IsoDBeta(TRootMuon *mu); 

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
  bool mumumu =false;  
  bool all = false;
  bool runHLT = true; 
  bool hasMu = false; 
  bool hasEl = false; 
  bool dilep =false; 
  bool singlelep = false;
  bool applyJetCleaning =false;
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
  string Channel = ""; 
  string xmlFileName = ""; 
  
  
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
  /// define channels
  //
  if(chanName.find("mumumu")!=string::npos) mumumu = true; 
  if(chanName.find("eemu")!=string::npos) eemu = true;
  if(chanName.find("mumue")!=string::npos) mumue = true;
  if(chanName.find("eee")!=string::npos) eee = true;
  if(chanName.find("all")!=string::npos) all = true; 
   if(mumumu)
  {
      cout << " --> Using the TriMuon channel <-- " << endl;
      Channel = "MuMuMu";
      xmlFileName = "config/Run2TriLepton.xml" ;
      dataLumi = 2700; //pb
      hasMu = true;
      hasEl = false;
      dilep = true;
      singlelep = false;
  }
 if(eee)
  {
      cout << " --> Using the TriElectron channel <-- " << endl;
      Channel = "ElElEl";
      xmlFileName = "config/Run2TriLepton.xml" ;
      dataLumi = 2700; //pb
      hasMu = false;
      hasEl = true;
      dilep = true;
      singlelep =false;
  }
  if(mumue)
  {
      cout << " --> Using the MuMuEl channel <-- " << endl;
      Channel = "MuMuEl";
      xmlFileName = "config/Run2TriLepton.xml" ;
      dataLumi = 2700; //pb
      hasMu = true;
      hasEl = true;
      dilep = true;
      singlelep =false;
  }
  if(eemu)
  {
      cout << " --> Using the ElElMu channel <-- " << endl;
      Channel = "ElElMu";
      xmlFileName = "config/Run2TriLepton.xml" ;
      dataLumi = 2700; //pb
      hasMu = true;
      hasEl = true;
      dilep = true;
      singlelep =false;
  }
  if(all)
  {
      cout << " --> Using the all  channel <-- " << endl;
      Channel = "All";
      xmlFileName = "config/Run2TriLepton.xml" ;
      dataLumi = 2700; //pb
      hasMu = true;
      hasEl = true;
      dilep = true;
      singlelep = true;
  }
/*  else
  {
      cerr << " ERROR: no channel specified " << endl;
      exit(1);
  }*/
  
  // Print information to a textfile
  ofstream infoFile;
  ofstream isoFile;
  ofstream jetFile;
  ofstream jetJECFile;
  ofstream jetSelFile;
  ofstream topFile; 
  ofstream mWtFile; 
  ofstream muSelFile; 
  ofstream muIniFile;
  string info_dir = "Information/"+Channel +"/";
  string iso_dir = "Isolation/"+Channel +"/";
  
  string info_date_dir = info_dir +  dateString +"/";
  string iso_date_dir = iso_dir +  dateString +"/";
  cout << "info dir " << info_dir.c_str() << endl; 
  mkdir(info_dir.c_str(),0777);
  mkdir(info_date_dir.c_str(),0777); 
  mkdir(iso_dir.c_str(),0777);
  mkdir(iso_date_dir.c_str(),0777);
  string infoName = info_date_dir + "information"; 
  infoName += "_"+ Channel;
  infoName += "_" + dName;
  infoName += "_" + JobNum; 
  infoName += ".txt"; 
  infoFile.open(infoName.c_str());
  infoFile.precision(3);  
  string isoName = iso_date_dir + "isolation";
  isoName += "_"+ Channel;
  isoName += "_" + dName;
  isoName += "_" + JobNum;
  isoName += ".txt";
  isoFile.open(isoName.c_str());
//  isoFile.precision(3);
  string jetName = info_date_dir + "jetinfo"; 
  jetName += "_"+ Channel;
  jetName += "_" + dName;
  jetName += "_" + JobNum;
  jetName += ".txt";
  jetFile.open(jetName.c_str());
  string jetJECName = info_date_dir + "jetinfoJEC"; 
  jetJECName += "_"+ Channel;
  jetJECName += "_" + dName;
  jetJECName += "_" + JobNum;
  jetJECName += ".txt";
  jetJECFile.open(jetJECName.c_str());
  string jetSelName = info_date_dir + "jetinfoSel";
  jetSelName += "_"+ Channel;
  jetSelName += "_" + dName;
  jetSelName += "_" + JobNum;
  jetSelName += ".txt";
  jetSelFile.open(jetSelName.c_str());
  string topName = info_date_dir + "topinfo";
  topName += "_"+ Channel;
  topName += "_" + dName;
  topName += "_" + JobNum;
  topName += ".txt";
  topFile.open(topName.c_str());
  string mWtName = info_date_dir + "mWtinfo";
  mWtName += "_"+ Channel;
  mWtName += "_" + dName;
  mWtName += "_" + JobNum;
  mWtName += ".txt";
  mWtFile.open(mWtName.c_str());
  string muSelName = info_date_dir + "muSelinfo";
  muSelName += "_"+ Channel;
  muSelName += "_" + dName;
  muSelName += "_" + JobNum;
  muSelName += ".txt";
  muSelFile.open(muSelName.c_str());
  string muIniName = info_date_dir + "muIniinfo";
  muIniName += "_"+ Channel;
  muIniName += "_" + dName;
  muIniName += "_" + JobNum;
  muIniName += ".txt";
  muIniFile.open(muIniName.c_str());
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
  Trigger* trigger_mumu  = new Trigger(1, 0, 0, 1);
  Trigger* trigger_ee  = new Trigger(0, 1, 0, 1);
  Trigger* trigger_emu  = new Trigger(1, 1, 0, 1) ;

  ///////////////////////
  // MET calculator /// 
  /////////////////////
  MEzCalculator* MEzCalculator; 
  

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

        LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data76X_25ns-Run246908-260627Cert.root", "pileup", "pileup");	

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
           trigger_mumu->bookTriggers(isData);
           trigger_ee->bookTriggers(isData);
           trigger_emu->bookTriggers(isData);
       



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
        bool   Wmu = false; 
        bool   Wel = false;
        bool   lep3 = false; 
        TLorentzVector metTLV;  
        string TriggBits; 
        string channel;
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
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
	   continueFlow = true; 
           lep3 = false; 
	   leading_jetPt = 0.; 
	   met = 0.;
	   leading_jet_btagDiscr = 0.; 
           TriggBits = ""; 
           channel = ""; 
	   pt_lept1 = pt_lept2 = pt_lept3 = 0. ; 
	   metTLV.Clear();
           metTLV.SetPxPyPzE(0,0,0,0); 
           selections.clear();
           
           selectionsnb.clear(); 
           selectionsnb.str(std::string());
           nCuts = 0;
           Wmu = false; 
           Wel = false;
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
	    init_jets_corrected = init_jets;  
	    
	    if(verbose>3)
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
          // to be applied from 76X v1 in our ttp
            HBHEnoise = event->getHBHENoiseFilter();
            HBHEIso = event->getHBHENoiseIsoFilter();
            CSCTight = event->getCSCTightHalo2015Filter(); 
            EcalDead = event->getEcalDeadCellTriggerPrimitiveFilter();
            eeBad = event->getEEBadScFilter(); 
 
           for(int iEl = 0 ; iEl < init_electrons.size() ; iEl ++){
              isoFile << evt_num << " sumChargedHadronPt=" << init_electrons[iEl]->chargedHadronIso(3) << ", sumNeutralHadronEt=" << init_electrons[iEl]->neutralHadronIso(3) << ", sumPhotonEt=" << init_electrons[iEl]->photonIso(3)<< ", effArea=" << EffectiveArea(init_electrons[iEl]) << endl;  
	   }
	   for(int iJet = 0; iJet < init_jets.size(); iJet++){
              TRootPFJet* tempJet = (TRootPFJet*) init_jets[iJet];
                double ptTemp = sqrt(tempJet->Px()*tempJet->Px()+tempJet->Py()*tempJet->Py());
//               jetFile << "EvtNb="<< evt_num << " jet_pt=" << tempJet->Pt() << " " << ptTemp << endl; 
                  jetFile << "EvtNb="<< evt_num << " jet_pt=" << tempJet->Pt() <<" jet_eta=" << tempJet->Eta() << " jet_phi=" << tempJet->Phi() << " NEMfraction="  << tempJet->neutralEmEnergyFraction() << " CEMfraction=" << tempJet->chargedEmEnergyFraction() << " NHfraction=" << tempJet->neutralHadronEnergyFraction() << " CHfraction=" << tempJet->chargedHadronEnergyFraction() << " Cmult=" << tempJet->chargedMultiplicity() << " nConst=" << tempJet->nConstituents() << endl;

//              jetFile << "EvtNb="<< evt_num << " jet_pt=" << tempJet->Pt() << " jet_eta=" << tempJet->Eta() << " jet_phi=" << tempJet->Phi()  << endl;
//              jetFile << "EvtNb="<< evt_num << " jet_pt=" << tempJet->Pt() << " jet_eta=" << init_jets[iJet]->Eta() << " jet_phi=" << init_jets[iJet]->Phi()  << endl;
           }
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
            bool trigged_mumu = false; 
	    bool trigged_ee = false; 
	    bool trigged_emu = false; 
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
                 if(all && isData)
		 {
		    if(dName.find("DoubleEG")!=string::npos && trigged_ee) trigged = true;
	            if(dName.find("DoubleMu")!=string::npos && trigged_mumu) trigged = true; 
		    if(dName.find("MuonEG")!=string::npos && !trigged_ee && !trigged_mumu && trigged_emu) trigged = true;
                 }
		 else if((trigged_emu || trigged_ee || trigged_mumu) && all) trigged = true;
                 if( trigged_ee && eee) trigged = true;
                 if( trigged_mumu && mumumu ) {trigged = true; nbTrig++;}
                 if( trigged_emu && !trigged_ee && !trigged_mumu  && (eemu || mumue)) trigged = true;
		 if(trigged_emu &&  trigged_ee && trigged_mumu) TriggBits = "111"; 
                 else if(!trigged_emu && !trigged_ee && !trigged_mumu) TriggBits = "000";
                 else if(!trigged_emu && !trigged_ee && trigged_mumu) TriggBits = "100";
                 else if(!trigged_emu && trigged_ee && !trigged_mumu) TriggBits = "010";
                 else if(trigged_emu && !trigged_ee && !trigged_mumu) TriggBits = "001";
                 else if(!trigged_emu && trigged_ee && trigged_mumu) TriggBits = "110";
                 else if(trigged_emu && trigged_ee && !trigged_mumu) TriggBits = "011";
                 else if(trigged_emu && !trigged_ee && trigged_mumu) TriggBits = "101";
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
 	   for(int iJet = 0; iJet < init_jets_corrected.size(); iJet++){
              TRootPFJet* tempJ = (TRootPFJet*) init_jets_corrected[iJet];
//              jetJECFile << "EvtNb="<< evt_num << " jet_pt=" << tempJ->Pt() <<" jet_eta=" << tempJ->Eta() << " jet_phi=" << tempJ->Phi() << " NEMfraction="  << tempJ->neutralEmEnergyFraction() << " CEMfraction=" << tempJ->chargedEmEnergyFraction() << " NHfraction=" << tempJ->neutralHadronEnergyFraction() << " CHfraction=" << tempJ->chargedHadronEnergyFraction() << " Cmult=" << tempJ->chargedMultiplicity() << " nConst=" << tempJ->nConstituents() << endl;
               jetJECFile << "EvtNb="<< evt_num << " jet_pt=" << tempJ->Pt() <<" jet_eta=" << tempJ->Eta() << " jet_phi=" << tempJ->Phi() << " jet_bDis=" << tempJ->btag_combinedInclusiveSecondaryVertexV2BJetTags() << endl;
   
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
            treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
            sort(mcParticles.begin(),mcParticles.end(),HighestPt());
	    // void TTreeLoader::LoadMCEvent(int, TopTree::TRootNPGenEvent*, std::vector<TopTree::TRootMCParticle*>&, bool) 
	    if (verbose>4) cout <<"Number of Muons, Electrons, Jets  ===>  " << endl << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< PreselectedJets.size()   << endl;
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
//	    cout << evt_num << " init " << init_jets_corrected.size() << " sel "  << selectedJets.size() << " bf cleaning " << PreselectedJets.size() << endl; 
            
            ////////////////////////////////////////////////
            // Pre cut operations
            ////////////////////////////////////////////////
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
	    // Met filters 
	    if(HBHEnoise && HBHEIso && CSCTight && EcalDead && eeBad && isGoodPV) passedMET = true;
	    PassedMETFilter = passedMET; 
	    
            for(int iJet = 0; iJet < selectedJets.size(); iJet++){
              TRootPFJet* tempJ = (TRootPFJet*) selectedJets[iJet];
              jetSelFile << "EvtNb="<< evt_num << " jet_pt=" << tempJ->Pt() <<" jet_eta=" << tempJ->Eta() << " jet_phi=" << tempJ->Phi() << " jet_bDis=" << tempJ->btag_combinedInclusiveSecondaryVertexV2BJetTags() << endl; 
     	    }
           for(int iMu = 0; iMu < selectedMuons.size(); iMu++){
	      muSelFile << "EvtNb="<< evt_num << " mu_pt=" << selectedMuons[iMu]->Pt() <<" mu_eta=" << selectedMuons[iMu]->Eta() << " mu_phi=" << selectedMuons[iMu]->Phi() << " mu_iso=" << IsoDBeta(selectedMuons[iMu]) << endl;
           }
           for(int iMu = 0; iMu < init_muons.size(); iMu++){
              muIniFile << "EvtNb="<< evt_num << " mu_pt=" << init_muons[iMu]->Pt() <<" mu_eta=" << init_muons[iMu]->Eta() << " mu_phi=" << init_muons[iMu]->Phi() << " mu_iso=" << IsoDBeta(init_muons[iMu]) << endl;
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
            // determine channels for synch
            //////////////////////////////////////////////////////
            if(selectedJets.size() > 0){
		leading_jetPt = selectedJets[0]->Pt(); 
                leading_jet_btagDiscr = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();		
	    }
            
            if(selectedMuons.size() > 2 && selectedElectrons.size() <= 2) { 
		channel = "mmm";
		pt_lept1 = selectedMuons[0]->Pt(); 
		pt_lept2 = selectedMuons[1]->Pt();
  		pt_lept3 = selectedMuons[2]->Pt();
	        iso_lept1 = IsoDBeta( selectedMuons[0]);
                iso_lept2 =  IsoDBeta(selectedMuons[1]);
                iso_lept3 = IsoDBeta(selectedMuons[2]);
	    }
	    else if(selectedElectrons.size() > 2 && selectedMuons.size() <= 2){
		channel = "eee"; 
		pt_lept1 = selectedElectrons[0]->Pt();
                pt_lept2 = selectedElectrons[1]->Pt();
                pt_lept3 = selectedElectrons[2]->Pt();
                iso_lept1 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                iso_lept2 = relPfIsoEl(selectedElectrons[1],event->fixedGridRhoFastjetAll());
                iso_lept3 = relPfIsoEl(selectedElectrons[2],event->fixedGridRhoFastjetAll());

	    }
            else if(selectedElectrons.size() > 2 && selectedMuons.size() > 2){
                cout << "SOMETHING IS WRONG " << endl; 
           }
	    else if(selectedElectrons.size() == 2 && selectedMuons.size() == 1){ 
		channel = "eem"; 
		if(selectedMuons[0]->Pt() > selectedElectrons[0]->Pt()){
	    		pt_lept1 = selectedMuons[0]->Pt();
        	        pt_lept2 = selectedElectrons[0]->Pt();
                	pt_lept3 = selectedElectrons[1]->Pt();
                        iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept2 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept3 = relPfIsoEl(selectedElectrons[1],event->fixedGridRhoFastjetAll());
		}
		else  if(selectedMuons[0]->Pt() < selectedElectrons[1]->Pt()){
                        pt_lept3 = selectedMuons[0]->Pt();
                        pt_lept1 = selectedElectrons[0]->Pt();
                        pt_lept2 = selectedElectrons[1]->Pt();
                        iso_lept3 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept1 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept2 = relPfIsoEl(selectedElectrons[1],event->fixedGridRhoFastjetAll());
                }
		else {
                        pt_lept2 = selectedMuons[0]->Pt();
                        pt_lept1 = selectedElectrons[0]->Pt();
                        pt_lept3 = selectedElectrons[1]->Pt();
                        iso_lept2 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept1 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept3 = relPfIsoEl(selectedElectrons[1],event->fixedGridRhoFastjetAll());
                }

	    }
	    else if(selectedElectrons.size() == 1 && selectedMuons.size() == 2){
                channel = "mme"; 
                if(selectedElectrons[0]->Pt() > selectedMuons[0]->Pt()){
                        pt_lept1 = selectedElectrons[0]->Pt();
                        pt_lept2 = selectedMuons[0]->Pt();
                        pt_lept3 = selectedMuons[1]->Pt();
                        iso_lept1 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept2 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept3 = (selectedMuons[1]->chargedHadronIso(4) + std::max(0.0, selectedMuons[1]->neutralHadronIso(4) + selectedMuons[1]->photonIso(4) - 0.5*selectedMuons[1]->puChargedHadronIso(4)))/selectedMuons[1]->Pt();


                }
                else  if(selectedElectrons[0]->Pt() < selectedMuons[1]->Pt()){
                        pt_lept3 = selectedElectrons[0]->Pt();
                        pt_lept1 = selectedMuons[0]->Pt();
                        pt_lept2 = selectedMuons[1]->Pt();
                        iso_lept3 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept2 = (selectedMuons[1]->chargedHadronIso(4) + std::max(0.0, selectedMuons[1]->neutralHadronIso(4) + selectedMuons[1]->photonIso(4) - 0.5*selectedMuons[1]->puChargedHadronIso(4)))/selectedMuons[1]->Pt();
                }
                else {
                        pt_lept2 = selectedElectrons[0]->Pt();
                        pt_lept1 = selectedMuons[0]->Pt();
                        pt_lept3 = selectedMuons[1]->Pt();
                        iso_lept2 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
                        iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept3 = (selectedMuons[1]->chargedHadronIso(4) + std::max(0.0, selectedMuons[1]->neutralHadronIso(4) + selectedMuons[1]->photonIso(4) - 0.5*selectedMuons[1]->puChargedHadronIso(4)))/selectedMuons[1]->Pt();
                }
            }
	    else{
              if(selectedMuons.size() == 2)
	      {
		id_lept3 = 0; 
                pt_lept1 = selectedMuons[0]->Pt();
                pt_lept2 = selectedMuons[1]->Pt();
                iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt(); // TO BE CHECKED
                iso_lept2 = (selectedMuons[1]->chargedHadronIso(4) + std::max(0.0, selectedMuons[1]->neutralHadronIso(4) + selectedMuons[1]->photonIso(4) - 0.5*selectedMuons[1]->puChargedHadronIso(4)))/selectedMuons[1]->Pt();
              }
	     else if(selectedElectrons.size() == 2){
               id_lept3 = 0;
                pt_lept1 = selectedElectrons[0]->Pt();
                pt_lept2 = selectedElectrons[1]->Pt();
                iso_lept1 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll()); // TO BE CHECKED
                iso_lept2 = relPfIsoEl(selectedElectrons[1],event->fixedGridRhoFastjetAll());
            }
            else if(selectedMuons.size() == 1 && selectedElectrons.size() == 1){
                id_lept3 = 0;
	     	if(selectedMuons[0]->Pt() > selectedElectrons[0]->Pt()){
			pt_lept1 = selectedMuons[0]->Pt();
			pt_lept2 = selectedElectrons[0]->Pt();
			iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
			iso_lept2 = relPfIsoEl(selectedElectrons[0],event->fixedGridRhoFastjetAll());
		}
		else{
                        pt_lept2 = selectedMuons[0]->Pt();
                        pt_lept1 = selectedElectrons[0]->Pt();
                        iso_lept2 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
                        iso_lept1 = relPfIsoEl(selectedElectrons[0], event->fixedGridRhoFastjetAll());
                }
	     }
	     else if(selectedMuons.size() == 1){
		 id_lept3 = 0;
		 id_lept2 = 0; 
		 pt_lept1 = selectedMuons[0]->Pt();
		 iso_lept1 = (selectedMuons[0]->chargedHadronIso(4) + std::max(0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.5*selectedMuons[0]->puChargedHadronIso(4)))/selectedMuons[0]->Pt();
             }
	     else if(selectedElectrons.size() == 1 && selectedMuons.size()==0){
		 id_lept3 = 0;
                 id_lept2 = 0;
		 pt_lept1 = selectedElectrons[0]->Pt();
		 iso_lept1 = relPfIsoEl(selectedElectrons[0], event->fixedGridRhoFastjetAll());

	     }
              channel = "nan";
            }
	    
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
               // continue;
            }
             if(dName.find("DoubleEG")!=string::npos && selectedElectrons.size() < 2) { continueFlow = false; }
             else if(dName.find("DoubleEG")!=string::npos) { nbEvents_test++ ;}
	     if(dName.find("DoubleMu")!=string::npos && selectedMuons.size() < 2) { continueFlow = false; }
             else if(dName.find("DoubleMu")!=string::npos) { nbEvents_test++ ;}
             if(dName.find("MuonEG")!=string::npos && (selectedElectrons.size() < 1 || selectedMuons.size() < 1)) { continueFlow = false; }
             else if(dName.find("MuonEG")!=string::npos){ nbEvents_test++ ;}
	   
            if(all && ((selectedMuons.size() + selectedElectrons.size()) != 3)){
               selections.push_back(0);
	       continueFlow = false; 
            }
	    else if(all && ((selectedMuons.size() + selectedElectrons.size()) == 3)){
               selections.push_back(1); 
	       if(continueFlow){
		 histo1D["cutFlow"]->Fill(1., eventweight);   
                 nCuts++;
                 nbEvents_1++;
               }
               lep3 = true;  
            }
            if(mumumu && (selectedMuons.size() != 3)){
               selections.push_back(0);
               continueFlow = false;
            }
            else if(mumumu && (selectedMuons.size() == 3)){
               selections.push_back(1);
               if(continueFlow){
                 histo1D["cutFlow"]->Fill(1., eventweight);
                 nCuts++;
                 nbEvents_1++;
               }
               lep3 = true;
            }
            if(eee && (selectedElectrons.size() != 3)){
               selections.push_back(0);
               continueFlow = false;
            }
            else if(eee && (selectedElectrons.size() == 3)){
               selections.push_back(1);
               if(continueFlow){
                 histo1D["cutFlow"]->Fill(1., eventweight);
                 nCuts++;
                 nbEvents_1++;
               }
               lep3 = true;
            }
            if(eemu && (selectedMuons.size() != 1|| selectedElectrons.size() != 2)){
               selections.push_back(0);
               continueFlow = false;
            }
            else if(eemu && (selectedMuons.size() == 1 && selectedElectrons.size() == 2)){
               selections.push_back(1);
               if(continueFlow){
                 histo1D["cutFlow"]->Fill(1., eventweight);
                 nCuts++;
                 nbEvents_1++;
               }
               lep3 = true;
            }
            if(mumue && (selectedMuons.size() != 2 || selectedElectrons.size() != 1)){
               selections.push_back(0);
               continueFlow = false;
            }
            else if(mumue && (selectedMuons.size() == 2 && selectedElectrons.size() == 1)){
               selections.push_back(1);
               if(continueFlow){
                 histo1D["cutFlow"]->Fill(1., eventweight);
                 nCuts++;
                 nbEvents_1++;
               }
               lep3 = true;
            }
            
            

            if(selectedMuons.size() == selectedLooseMuons.size() && continueFlow) nbEvents_1m++;
	    else continueFlow = false; 
            if(selectedLooseElectrons.size() == selectedElectrons.size() && continueFlow)  nbEvents_2m++;
	    else continueFlow = false; 
	    if((selectedMuons.size() != selectedLooseMuons.size()) || (selectedLooseElectrons.size() != selectedElectrons.size())){
              selections.push_back(0); 
              continueFlow = false; 
            }
	    else{
	     selections.push_back(1); 
            }

	    double met_px = mets[0]->Px();
	    double met_py = mets[0]->Py();
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
	    met = met_Pt; 
	    met_Phi = mets[0]->Phi(); 
	    met_Eta = mets[0]->Eta();
            
	    puSF = PUweight;
	    btagSF = btagWeight;  


	    Zlep0.Clear(); 
	    Zlep1.Clear();
            Wlep.Clear(); 
            Wlep.SetPxPyPzE(0,0,0,0); 

	    // check sign
	    bool OS = false; 
	    if(selectedElectrons.size() == 2){ 
	       if(selectedElectrons[0]->charge() == selectedElectrons[1]->charge()){ OS = false;  }
	       else {
	          OS = true;
	          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
  	          Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
		  if(lep3) {
		    Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
		    Wmu = true ; 	
		    
                  }
	       }   
	    }
            else if(selectedMuons.size() == 2){ 
	    	if(selectedMuons[0]->charge() == selectedMuons[1]->charge()){ OS = false; }
		else{
		  OS = true;
		  Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
  		  Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
		  if(lep3) {
		     Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
		     Wel = true; 
                  }
		}
	    }
	    else if(selectedMuons.size()==3)
	    {
              Wmu = true; 
              bool first = false; 
              bool second = false; 
              bool third = false; 
	      if(selectedMuons[0]->charge() != selectedMuons[1]->charge()) first = true; 
  	      if(selectedMuons[2]->charge() != selectedMuons[1]->charge()) second = true; 
              if(selectedMuons[0]->charge() != selectedMuons[2]->charge()) third = true; 
	      if(first || second || third) OS = true; 
	      else OS = false; ; 
	      if(first && !second && !third){
 			 Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
  			 Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
			 Wlep.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
               } 
               else if(third && !second && !first){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
			 Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
	       }
               else if(second && !first && !third){
                         Zlep0.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
               }
	       else if(OS)
               {
                     TLorentzVector tempMu0; 
		     tempMu0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy()); 
		     TLorentzVector tempMu1; 
	             tempMu1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy()); 
	             TLorentzVector tempMu2; 
		     tempMu2.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());  
		     double mass01 = (tempMu0 + tempMu1).M(); 
                     double mass02 = (tempMu0 + tempMu2).M();
                     double mass12 = (tempMu2 + tempMu1).M();
		    if(first && second && !third){
			if(fabs(mass01-90.0) < fabs(mass12-90.0) ){
			  Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
			}
			else{
    		         Zlep0.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
			}
                    }
		    else if(first && third && !second)
		    {
                        if(fabs(mass01-90.0) < fabs(mass02-90.0) ){
                          Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                        }
                        else{
                          Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
			}
		    }
		    else if(third && second && !first)
		    {
			if(fabs(mass02-90.0) < fabs(mass12-90.0) ){
                          Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                        }
                        else{
                         Zlep0.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                        }
                    }
		    else if (first && second && third){
			if(fabs(mass01-90.0) < fabs(mass12-90.0) && fabs(mass01-90.0) < fabs(mass02-90.0) ){
                          Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                        }
                        else if( fabs(mass12-90.0) < fabs(mass01-90.0) && fabs(mass12-90.0) < fabs(mass02-90.0) ){
                         Zlep0.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                        }			 
			else{
                          Zlep0.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(), selectedMuons[0]->Pz(), selectedMuons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(), selectedMuons[2]->Pz(), selectedMuons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(), selectedMuons[1]->Pz(), selectedMuons[1]->Energy());
                        }
                    }
                }
            }
	    else if(selectedElectrons.size()==3)
	    {
              Wel = true; 
              bool first = false; 
              bool second = false; 
              bool third = false; 
	      if(selectedElectrons[0]->charge() != selectedElectrons[1]->charge()) first = true; 
  	      if(selectedElectrons[2]->charge() != selectedElectrons[1]->charge()) second = true; 
              if(selectedElectrons[0]->charge() != selectedElectrons[2]->charge()) third = true; 
	      if(first || second || third) OS = true; 
	      else continue;
	      if(first && !second && !third){
 			 Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
  			 Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
			 Wlep.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
               } 
                 else if(third && !second && !first){
			 OS = true; 
			 Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
			 Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
			 Wlep.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
		}
                else if(second && !first && !third){
                         Zlep0.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                }
	        else if(OS)
                {
                     TLorentzVector tempMu0; 
		     tempMu0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy()); 
		     TLorentzVector tempMu1; 
	             tempMu1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy()); 
	             TLorentzVector tempMu2; 
		     tempMu2.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());  
		     double mass01 = (tempMu0 + tempMu1).M(); 
                     double mass02 = (tempMu0 + tempMu2).M();
                     double mass12 = (tempMu2 + tempMu1).M();
		    if(first && second && !third){
			if(fabs(mass01-90.0) < fabs(mass12-90.0) ){
			  Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
			}
			else{
    		         Zlep0.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
			}
                    }
		    else if(first && third && !second)
		    {
                        if(fabs(mass01-90.0) < fabs(mass02-90.0) ){
                          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                        }
                        else{
                          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
			}
		    }
		    else if(third && second && !first)
		    {
			if(fabs(mass02-90.0) < fabs(mass12-90.0) ){
                          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                        }
                        else{
                         Zlep0.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                        }
                    }
		    else if (first && second && third){
			if(fabs(mass01-90.0) < fabs(mass12-90.0) && fabs(mass01-90.0) < fabs(mass02-90.0) ){
                          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                        }
                        else if( fabs(mass12-90.0) < fabs(mass01-90.0) && fabs(mass12-90.0) < fabs(mass02-90.0) ){
                         Zlep0.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                         Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                         Wlep.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                        }			 
			else{
                          Zlep0.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(), selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy());
                          Zlep1.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(), selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy());
                          Wlep.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(), selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy());
                        }
                    }
                }
            }

            if(!OS){
                continueFlow = false; 
 		//coninue;
            }
            else if(continueFlow){ 
              histo1D["cutFlow"]->Fill(2., eventweight); 
              nCuts++;
              nbEvents_2++;  
	    }
	    Zboson.Clear(); 
            if(OS) Zboson.SetPxPyPzE(( Zlep0 + Zlep1).Px() ,( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Py(),( Zlep0 + Zlep1).Energy()) ;
            if(OS) Zboson_M = (Zlep0+Zlep1).M();
            else if(!OS) Zboson_M = 0;
//            cout << " Zmass" << Zboson_M << endl;  
	    if(Zboson_M < 76 || Zboson_M > 106)
            {
                selections.push_back(0); 
		continueFlow = false; 
                eventSelected = false;
                // continue; 
            }  
            else{
               selections.push_back(1);
	       if(continueFlow){
		  nCuts++; 
	          nbEvents_3++; 
	          histo1D["cutFlow"]->Fill(3., eventweight);
                 eventSelected = true;
		}
            }
            if(selectedJets.size() == 0){
		selections.push_back(0);
		continueFlow = false; 
               // continue; 
            }
            else{
	      selections.push_back(1);
	      if(continueFlow){ 
	       histo1D["cutFlow"]->Fill(4., eventweight);
               nCuts++;
               nbEvents_4++;  
	      }
            }
//            cout << " after " << nCuts << " " << nbEvents_3 << endl;
            if(selectedCSVLBJets.size() != 1){
               selections.push_back(0); 
	       continueFlow = false; 
	       //continue; 
            }
	    else{
		selections.push_back(1); 
		if(continueFlow){
	          histo1D["cutFlow"]->Fill(5., eventweight);
                  nCuts++;
                  nbEvents_5++;  
		}
	    }
            //double mWtsecond = 0.;
             
	    if(Wel|| Wmu){
                 //double phis = Wlep.Phi() - mets[0]->Phi();
                 //double cosphis = TMath::Cos(phis);
		 mWt = TMath::Sqrt((Wlep.Pt() + met_Pt)*(Wlep.Pt() +met_Pt)-(Wlep.Px() + met_px)*(Wlep.Px() + met_px) - (Wlep.Py() + met_py)* (Wlep.Py() + met_py));
                 //mWtsecond = TMath::Sqrt(2*Wlep.Pt() * met_Pt*(1-cosphis));
	    }
       
	    else mWt = 0.;
            mWtFile << "EvtNb="<< evt_num << " mWt=" << mWt << " met_Pt=" << met_Pt << " WlepPt=" << Wlep.Pt() << "CosPhi=" << TMath::Cos(Wlep.Phi() - met_Phi) << endl; //" second=" <<  mWtsecond << endl;   
            if(mWt < 20){
               selections.push_back(0); 
	       continueFlow = false; 
//               continue; 
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
            double met_pz = 0.; // has to be adapted !!! 
            metTLV.SetPxPyPzE(met_px,met_py,met_pz,TMath::Sqrt(met_px*met_px+met_py*met_py+met_pz*met_pz));
 	    SMbjet.Clear(); 
            if(selectedCSVLBJets.size() > 0){
                 SMbjet.SetPxPyPzE(selectedCSVLBJets[0]->Px(),selectedCSVLBJets[0]->Py(),selectedCSVLBJets[0]->Pz(),selectedCSVLBJets[0]->Energy());            
	         if(Wel|| Wmu)  SMtop_M = (Wlep+SMbjet+metTLV).M();
		  else SMtop_M = 0.; 
            }
	    else  SMtop_M = 0. ;
            if(continueFlow) topFile << "EvtNb="<< evt_num << " Bjet_pt=" << SMbjet.Pt() <<" Bjet_px=" << SMbjet.Px() << " Bjet_py=" << SMbjet.Py() << " Bjet_pz()=" << SMbjet.Pz() << " Bjet_Energy=" << SMbjet.Energy() << " Wlep_pt=" << Wlep.Pt() <<" Wlep_px=" << Wlep.Px() << " Wlep_py=" << Wlep.Py() << " Wlep_pz()=" << Wlep.Pz() << " Wlep_Energy=" << Wlep.Energy() << " met_Pt=" << metTLV.Pt() <<" met_px=" << metTLV.Px() << " met_py=" << metTLV.Py() << " met_pz()=" << metTLV.Pz() << " met_Energy=" << metTLV.Energy() << " topmass= " << SMtop_M << endl;  

//	    cjet.Clear(); 
//	    cjet = FCNCjetCalculator(selectedCSVLLJets,selectedCSVLBJets, Zboson ,3);

 //           FCNCtop_M = (Zboson+cjet).M();

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
	    if(isGoodPV && passedMET && continueFlow){ 
	      histo1D["cutFlow"]->Fill(8., eventweight);
              nCuts++;
               nbEvents_8++;
	    }
             //////////////////////////////////////
	    //  DO STUFF WITH SELECTED EVENTS ////
	    ////////////////////////////////////// 
	    if(eventSelected){ 
	      nbSelectedEvents++; 
	      myTree->Fill(); 
 	    }
	    if(selections.size() != 8) cout << "ERROR SOMETHING WENT WRONG WITH THE SELECTIONS " << endl; 
            for(int inb = 0; inb <selections.size(); inb++)
            {
                 selectionsnb << selections[inb];
            }
//            infoFile << "|" << evt_num << "|"  << TriggBits << "|"  <<channel << "|"  << pt_lept1 << "|" << pt_lept2 << "|" << pt_lept3 << "|" << iso_lept1 << "|" << iso_lept2 << "|" << iso_lept3 << "|" << id_lept1 << "|" << id_lept2 << "|" << id_lept3 << "|" << leading_jetPt << "|" << leading_jet_btagDiscr << "|" << met << "|" << selectionsnb.str() << endl;    
	    infoFile << "|" << evt_num << "|"  << TriggBits << "|"  <<channel << "|"  << pt_lept1 << "|" << pt_lept2 << "|" << pt_lept3 << "|" << iso_lept1 << "|" << iso_lept2 << "|" << iso_lept3 << "|" << leading_jetPt << "|" << leading_jet_btagDiscr << "|" << met << "|" << selectionsnb.str() << endl;
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
	cout << "2 lep " << nbEvents_test << endl; 
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
 //	for(int j = 0; j < 9; j++){       cout << cutstep[j] << endl; }
        sumW = (int) sumWeights; 
        nEv = (int) nEvents; 
	globalTree->Fill(); 
        if(verbose == 0) cout << "end eventloop" << endl; 

/*	infoFile << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl;
        infoFile << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
        infoFile << nbBaseline << " baseline events out of trigged " << nbTrig <<  " selected " << endl;
        infoFile << setprecision(2) << ((double)nbGPV/(double)nbEvents)*100 << " % of the initial events stay after Good PV" << endl;
        nfoFile << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
        infoFile << setprecision(2) << ((double)nbTrig/(double)nbGPV)*100 << " % of the GPV  events stay after Trigger" << endl;
*/        cout << nbSelectedEvents << " events out of initial " << nbEvents <<  " selected " << endl; 
        cout << nbSelectedEvents << " events out of trigged  " << nbTrig <<  " selected " << endl;
       // cout << nbBaseline << " baseline events out of trigged " << nbTrig <<  " selected " << endl;
       // cout << setprecision(2) << ((double)nbGPV/(double)nbEvents)*100 << " % of the initial events stay after Good PV" << endl; 
	cout << setprecision(2) << ((double)nbTrig/(double)nbEvents)*100 << " % of the initial events stay after Trigger" << endl;
       // cout << setprecision(2) << ((double)nbTrig/(double)nbGPV)*100 << " % of the GPV  events stay after Trigger" << endl;
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
        isoFile.close();
	topFile.close();
        jetFile.close();
	jetJECFile.close();
	jetSelFile.close();
	muSelFile.close();
        mWtFile.close();  
        muIniFile.close(); 
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


double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy)
{
  double M_W  = 80.4;
  double M_mu =  0.10566; // 105.66 MeV/c^2
  double M_el = 0.000510999; // 0.510998910 Mev/c^2
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


