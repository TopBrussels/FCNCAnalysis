#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>

#include <cmath>
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
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;


Bool_t debug = true;



/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;




int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO

  if (debug)
    {
      cout << "list of arguments are ..." << endl;
      for (int n_arg=1; n_arg<argc; n_arg++)
	{
	  std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl;
	}
    }
    


  if(argc < 14 )
    {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
	for (int n_arg=1; n_arg<argc; n_arg++)
	  {
	    std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl; 
	  }
        return 1;
    }

    //Placing arguments in properly typed variables for Dataset creation


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
    // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list) 
    const int JobNum                = strtol(argv[argc-3], NULL, 10);
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);


    // all the files are stored from arg 11 to argc-2
    vector<string> vecfileNames;
    for(int args = 11; args < argc-3; args++)
      {
	vecfileNames.push_back(argv[args]);
      }
    
    if (debug){
      cout << "The list of file to run over will be printed..." << endl;
      for ( int nfiles = 0; nfiles < vecfileNames.size(); nfiles++)
	{
	  cout << "file number " << nfiles << " is " << vecfileNames[nfiles] << endl;
	}
    }




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
    }
     
    cout << "----------------------------------------" << endl;
   

   // ofstream eventlist;
  //  eventlist.open ("interesting_events_mu.txt");


    int passed = 0;
    int eventCount = 0;



    clock_t start = clock();

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the Displaced Top search ! "           << endl;
    cout << "*************************************************************" << endl;


    string postfix = "_Run2_TopTree_Study_" + dName; // to relabel the names of the output file

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels 
    bool ee = false; 
    bool emu = true; 
    bool mumu = false; 
    bool runHLT = true;
    bool printTrigger = false;  
    bool applyJetCleaning = true; 

    if(emu)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/Run2TriLepton_samples.xml";
    }
    else if(mumu)
    {
        cout << " --> Using the Muon-Muon channel..." << endl;
        channelpostfix = "_MuMu";
        xmlFileName = "config/Run2TriLepton_samples.xml";
    }
    else if(ee)
    {
        cout << " --> Using the Electron-Electron channel..." << endl;
        channelpostfix = "_ElEl";
        xmlFileName = "config/Run2TriLepton_samples.xml";
    }
    else
    {
        cerr<<"Correct Di-lepton Channel not selected."<<endl;
        exit(1);
    }



    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
//    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = false;
    anaEnv.loadGenJetCollection = true;
    anaEnv.loadGenEventCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.loadTrackMETCollection = false;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;
    int verbose = 2;//anaEnv.Verbose;



    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);


    string dataSetName;

    ////////////////////////////////
    //  Event Scale Factor
    ////////////////////////////////

    string pathToCaliDir="../TopTreeAnalysisBase/Calibrations/";


    /// Leptons

    // Muon SF
    double muonSFID, muonSFIso;
    MuonSFWeight *muonSFWeight_ = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"Muon_SF_TopEA.root","SF_totErr", false, false); // (... , ... , debug, print warning)
    MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
    //    MuonSFWeight *muonSFWeightID_M = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
    //  MuonSFWeight *muonSFWeightID_L = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
  
    MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Tight RelIso, Tight ID
    //   MuonSFWeight *muonSFWeightIso_TM = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Tight RelIso, Medium ID
    //   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Tight ID
    //   MuonSFWeight *muonSFWeightIso_LM = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_MediumID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Medium ID
    //   MuonSFWeight *muonSFWeightIso_LT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_LooseRelIso_DEN_LooseID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Loose RelIso, Loose ID


    // Electron SF
    string electronFile= "Elec_SF_TopEA.root";
    ElectronSFWeight *electronSFWeight_ = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electronFile,"GlobalSF", false, false); // (... , ... , debug, print warning)

    // initialize trigger
    //                 
    Trigger* trigger = new Trigger(1, 1, 0, 1);

  


    // PU SF
   LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert.root","pileup60","pileup");
//   LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root","pileup","pileup");
    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    int ndatasets = datasets.size() - 1 ;

    

    //Output ROOT file
    string outputDirectory("/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/FCNCAnalysis/MACRO_Output_"+channelpostfix);
    mkdir(outputDirectory.c_str(),0777);

    // add jobs number at the end of file if 
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();

    string rootFileName (outputDirectory+"/Trial"+postfix+channelpostfix+".root");
    if (strJobNum != "0")
      {
	cout << "strJobNum is " << strJobNum << endl;
        rootFileName = outputDirectory+"/Trial"+postfix+channelpostfix+"_"+strJobNum+".root";
      }
    
    
    

    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");


    //Global variable
    TRootEvent* event = 0;

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////
    /*
    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    */

    ///////////////////
    // 1D histograms //
    ///////////////////

    map <string,TH1F*> histo1D;
    
    std::string  titlePlot = ""; 
    titlePlot = "cutFlow"+channelpostfix; 
    histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 13,-0.5,12.5);



    ///////////////////
    // 2D histograms //
    ///////////////////
    //    histo2D["HTLepSep"] = new TH2F("HTLepSep","dR_{ll}:HT",50,0,1000, 20, 0,4);

    //Plots
    string pathPNG = "/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/FCNCAnalysis/MSPlots_"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);


    // define cuts here

    // electron
    float el_pt_cut =20.; // 42
    float el_eta_cut = 2.5;


    // muon
    float mu_pt_cut = 20.; // 40
    float mu_eta_cut = 2.1;
    float mu_iso_cut = 0.15;
 
    //jets
    float jet_pt_cut = 30.;
    float jet_eta_cut = 2.4;

    // convert into string

/*    std::ostringstream el_pt_cut_strs, el_eta_cut_strs, mu_pt_cut_strs, mu_eta_cut_strs, mu_iso_cut_strs, jet_pt_cut_strs, jet_eta_cut_strs;
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
*/
    // check
    cout <<"Making directory :"<< pathPNG  <<endl;
    

    
    



    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {

        cout << "Load Dataset" << endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset	
        string previousFilename = "";
        int iFile = -1;
        dataSetName = datasets[d]->Name();


	cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
	

        //////////////////////////////////////////////
        // Setup Date string and nTuple for output  //
        //////////////////////////////////////////////

        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );

        int year = now->tm_year + 1900;
        int month =  now->tm_mon + 1;
        int day = now->tm_mday;
        int hour = now->tm_hour;
        int min = now->tm_min;
        int sec = now->tm_sec;

        string year_str;
        string month_str;
        string day_str;
        string hour_str;
        string min_str;
        string sec_str;

        ostringstream convert;   // stream used for the conversion
        convert << year;      // insert the textual representation of 'Number' in the characters in the stream
        year_str = convert.str();
        convert.str("");
        convert.clear();
        convert << month;      // insert the textual representation of 'Number' in the characters in the stream
        month_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        day_str = convert.str();
        convert.str("");
        convert.clear();
        convert << hour;      // insert the textual representation of 'Number' in the characters in the stream
        hour_str = convert.str();
        convert.str("");
        convert.clear();
        convert << min;      // insert the textual representation of 'Number' in the characters in the stream
        min_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        sec_str = convert.str();
        convert.str("");
        convert.clear();


        string date_str = day_str + "_" + month_str + "_" + year_str;

        cout <<"DATE STRING   "<<date_str << endl;

        string channel_dir = "/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/FCNCAnalysis/Output"+channelpostfix;
        string date_dir = channel_dir+"/Output" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);

       if(debug) cout << "Make Ntuple name" << endl; 


 //       string Ntupname = "/user/ivanpari/CMSSW_7_4_15/src/TopBrussels/FCNCAnalysis/Output"+channelpostfix+"/Output"+ date_str  +"/Output_" + dataSetName +postfix + ".root";
//        string Ntuptitle = "Output" + channelpostfix;
//        if(debug) cout << Ntupname.c_str() << " " << Ntuptitle.c_str() << endl; 
//        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
//        if(debug) cout << "made ntuple root file " << endl; 

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
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Int_t charge_jet[20];
        Double_t bdisc_jet[20];

	// event related variables
	Int_t run_num;
	Int_t evt_num;
	Int_t lumi_num;
	Int_t nvtx;
	Int_t npu;
	Double_t puSF;


	// define the output tree                                                                                       
	if(debug) cout << "accessing fout" << endl;   
	fout->cd();
        if(debug) cout << "accessed fout " << endl; 
        TTree* myTree = new TTree("tree","tree");
        if(debug) cout << "made  myTree" << endl; 

      	// electrons
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");//                                                            
        myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
      //  myTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
      //  myTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
      //  myTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
      //  myTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
      //  myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
      //  myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
	//myTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
	//myTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
	//myTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
	//myTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
	//myTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
	//myTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
	//myTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
	//myTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
	//myTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
	//myTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
	myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
	

	// muons
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        //myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        //myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        //myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        //myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        //myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        //myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        //myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
	//myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
	myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");

        // jets
        myTree->Branch("nJets",&nJets,"nJets/I");
        myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");       
        myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");

	// event related variables
	myTree->Branch("run_num",&run_num,"run_num/I");
	myTree->Branch("evt_num",&evt_num,"evt_num/I");
	myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
	myTree->Branch("nvtx",&nvtx,"nvtx/I");
	myTree->Branch("npu",&npu,"npu/I");
	myTree->Branch("puSF",&puSF,"puSF/D");	



        


       // book trigger

       if(runHLT){  trigger->bookTriggers(isData); };	






        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////


        int start = 0;
        if(debug) cout << " ending " << endl;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();
        if(debug) cout << "ending" << endl; 
        cout <<"Number of events in total dataset = "<<  ending  <<endl;


        int event_start = startEvent;
        if (verbose > 1) cout << " - Loop over events " << endl;


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
            return 1;
        }


        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootMuon*>     selectedMuons;
        
        // initial variables
        vector < TRootVertex* >   vertex;
        vector < TRootMuon* >     init_muons;
        vector < TRootElectron* > init_electrons;
        vector < TRootJet* >      init_jets;
        vector < TRootMET* >      mets;

        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
	
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
	  nElectrons = 0;
	  nMuons = 0; 
	  nJets = 0;
          selectedElectrons.clear();
          selectedMuons.clear();
	  selectedJets.clear();  
          vertex.clear(); 
	  init_muons.clear(); 
          init_jets.clear(); 
          mets.clear(); 
	  
	  double ievt_d = ievt;
	  currentfrac = ievt_d/end_d;
	  if (debug)cout << endl << endl << "Starting a new event loop!"<<endl;
	  
	  if(ievt%1000 == 0)
            {
	      std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }
	  
	  event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets,  mets, debug);  //load event
	  
	  
	  if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
	  if(debug) cout <<"Number of Muons Loaded: " << init_muons.size() <<endl; 
	  if(debug) cout << "Number of Jets Loaded: " << init_jets.size() << endl;  

	  
	  //////////////////
	  //Loading Gen jets
	  //////////////////
	  
	  
	  datasets[d]->eventTree()->LoadTree(ievt); 
	  string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	  int currentRun = event->runId();
	  run_num=event->runId();
	  evt_num=event->eventId();
	  lumi_num=event->lumiBlockId();
	  nvtx=vertex.size();
	  npu=(int)event->nTruePU();
	  
	  

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


	    ///////////////////////////////////////////
	    //  Pile up Scale Factor
	    ///////////////////////////////////////////

	    double puWeight = LumiWeights.ITweight( npu ); // simplest reweighting, just use reconstructed number of PV
	    puSF=puWeight;
	    if (isData) puSF =1;


            ///////////////////////
            // JER smearing
            //////////////////////

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////





            // Declare selection instance
            Run2Selection selection(init_jets, init_muons, init_electrons, mets);

            // Define object selection cuts

	    if (debug)cout<<"Getting Jets"<<endl; 
	    selectedJets                                        = selection.GetSelectedJets(jet_pt_cut,jet_eta_cut, true, "Tight"); // Relying solely on cuts defined in setPFJetCuts()

	    // make a new collections of muons
	    if (debug)cout<<"Getting Muons"<<endl;
	    selectedMuons = selection.GetSelectedMuons(mu_pt_cut, mu_eta_cut, mu_iso_cut, "Tight", "Spring15"); // pt, eta, iso // run normally

	    // make a new collections of electrons
	    if (debug)cout<<"Getting Electrons"<<endl;
	    selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut, "Tight","Spring15_25ns",true);// pt, eta






            //////////////////////////////////
            // Preselection Lepton Operations //
            //////////////////////////////////
           if(applyJetCleaning){
            if(debug) cout << "Applying jet cleaning " << endl; 
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
             if(debug){
             if( OrigSize != selectedJets.size()) cout << "--> original = " << OrigSize  << " after cleaning = " << selectedJets.size() << endl;
             else cout << "--> no change" << endl; 
             }
           }	    

          



            //////////////////////
            // Sync'ing cutflow //
            //////////////////////


            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;

            
             // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);



            eventCount++;


            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection
           //  if (!isGoodPV) continue; // Check that there is a good Primary Vertex

	    if (debug) cout <<"Number of Muons, Electrons, Jets  ===>  " << endl << selectedMuons.size() <<" "  << selectedElectrons.size()<<" "<< selectedJets.size()   << endl;


            if (debug)	cout <<"applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections




            //////////////////////////
            // Electron Based Plots //
            //////////////////////////
	    if (debug) cout << "before electrons loop" << endl;

	    
	    nElectrons=0;
            for (Int_t selel =0; selel < selectedElectrons.size() ; selel++ )
	    {
	      
              pt_electron[nElectrons]=selectedElectrons[selel]->Pt();
	      phi_electron[nElectrons]=selectedElectrons[selel]->Phi();
	      eta_electron[nElectrons]=selectedElectrons[selel]->Eta();
	      eta_superCluster_electron[nElectrons]=selectedElectrons[selel]->superClusterEta();
	      E_electron[nElectrons]=selectedElectrons[selel]->E();
//	      d0_electron[nElectrons]=selectedElectrons[selel]->d0();
//	      d0BeamSpot_electron[nElectrons]=selectedElectrons[selel]->d0BeamSpot();
//	      chargedHadronIso_electron[nElectrons]=selectedElectrons[selel]->chargedHadronIso(3);
//	      neutralHadronIso_electron[nElectrons]=selectedElectrons[selel]->neutralHadronIso(3);
//	      photonIso_electron[nElectrons]=selectedElectrons[selel]->photonIso(3);
//	      pfIso_electron[nElectrons]=selectedElectrons[selel]->relPfIso(3,0);
	      charge_electron[nElectrons]=selectedElectrons[selel]->charge();
//	      sigmaIEtaIEta_electron[nElectrons]=selectedElectrons[selel]->sigmaIEtaIEta();
//	      deltaEtaIn_electron[nElectrons]=selectedElectrons[selel]->deltaEtaIn();
//	      deltaPhiIn_electron[nElectrons]=selectedElectrons[selel]->deltaPhiIn();
//	      hadronicOverEm_electron[nElectrons]=selectedElectrons[selel]->hadronicOverEm();
//	      missingHits_electron[nElectrons]=selectedElectrons[selel]->missingHits();
//	      passConversion_electron[nElectrons]=selectedElectrons[selel]->passConversion();
//	      isEBEEGap[nElectrons]=selectedElectrons[selel]->isEBEEGap();
	      sf_electron[nElectrons]=electronSFWeight_->at(selectedElectrons[selel]->Eta(),selectedElectrons[selel]->Pt(),0);
              nElectrons++;
            }


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
//	      d0_muon[nMuons]=selectedMuons[selmu]->d0();
//	      d0BeamSpot_muon[nMuons]=selectedMuons[selmu]->d0BeamSpot();
//	      chargedHadronIso_muon[nMuons]=selectedMuons[selmu]->chargedHadronIso(4);
//	      neutralHadronIso_muon[nMuons]=selectedMuons[selmu]->neutralHadronIso(4);
//	      photonIso_muon[nMuons]=selectedMuons[selmu]->photonIso(4);
//               pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
	      charge_muon[nMuons]=selectedMuons[selmu]->charge();
	      sf_muon[nMuons]=muonSFWeightIso_TT->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID_T->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
              nMuons++;
            }
            

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
                nJets++;

            } 




 

            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();



            //////////////////
            //Filling Tree//
            //////////////////


	    // filling the cutflow and the histo 1 D

	    // preselection cut
	    histo1D["h_cutFlow"]->Fill(0., 1);
//            if(runHLT && !trigged) continue;
            if(runHLT && trigged){
              histo1D["h_cutFlow"]->Fill(1., 1); 
//              if (!isGoodPV) continue; 
              if(isGoodPV){
                 histo1D["h_cutFlow"]->Fill(2., 1);
//                 if(selectedElectrons.size() + selectedMuons.size()==0) continue; 
                 if(selectedElectrons.size() + selectedMuons.size()>0){ 
                    histo1D["h_cutFlow"]->Fill(3., 1);
//                    if(selectedElectrons.size() + selectedMuons.size()<2) continue;
                    if(selectedElectrons.size() + selectedMuons.size()>1){
                        histo1D["h_cutFlow"]->Fill(4., 1);
                        myTree->Fill(); 
                       passed++;			    
//                       if(selectedElectrons.size() + selectedMuons.size()<3) continue;
                       if(selectedElectrons.size() + selectedMuons.size()>2){
                          histo1D["h_cutFlow"]->Fill(5., 1);

			   if (debug) {
				cout << "npu vtx " << nvtx << endl;
				cout << "npu is " << npu << endl;
				cout << "puSF is" << puSF << endl;
	      		   }
		       }
		     }
		  }
               }
            }

	    
            

	    

        } //End Loop on Events

	fout->cd();
	///Write histogram
	for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
	  {
	    cout << "1D Plot: " << it->first << endl;
	    // TCanvas ctemp = 
	    
	    TH1F *temp = it->second;
	    string name = it->first;
	    cout << name << endl; 
	    temp->Draw();  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
	    //c1->SaveAs(pathPNG+"1DPlot/"+name modeString[mode] + "_" + cutLabel + ".png");
	    //	    temp->Write(fout, name, true, pathPNG+"1DPlot/", "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
	  }
	
	if (debug) cout << "Done writing the Tree" << endl;
	fout->Write();   
	fout->Close();
//        tupfile->Close();
        cout <<"n events after all the cuts is  =  "<< passed <<endl;
        cout << "Event Count: " << eventCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    //eventlist.close();

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;





    delete fout;
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}


