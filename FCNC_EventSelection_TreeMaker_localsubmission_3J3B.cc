//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.            ////
//////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:                              //ttbar @ NNLO 8 TeV:
//all-had ->679 * .46 = 312.34                      //all-had -> 245.8 * .46 = 113.068
//semi-lep ->679 *.45 = 305.55                      //semi-lep-> 245.8 * .45 = 110.61
//di-lep-> 679* .09 = 61.113                        //di-lep ->  245.8 * .09 = 22.122

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
//#include "TopTreeAnalysisBase/Selection/interface/FCNC_1L3BSelectionTable.h"
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
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

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


int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 14)
    {
        std::cerr << "INVALID INPUT FROM XMLFILE.  CHECK XML IMPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
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

    vector<string> vecfileNames;
    for(int args = 11; args < argc-2; args++)
    {
        vecfileNames.push_back(argv[args]);
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
    cout << "----------------------------------------" << endl;

    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();

    //Initializing event counts etc.
    int passed = 0;
    int preTrig = 0;
    int postTrig = 0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;

    //Initializing CSVv2 b-tag WP
	  float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
	  float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
	  float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.

    clock_t start = clock();
    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FCNC_1L3B search ! "           << endl;
    cout << "*************************************************************" << endl;



    ///////////////////////////////////////////////////////////////
    // Initialize scale&reweight-handlings
    //////////////////////////////////////////////////////////////
    bool bTagReweight = true;
    bool bLeptonSF = true;
    
    int doJESShift = 0; // 0: off 1: minus 2: plus
    cout << "doJESShift: " << doJESShift << endl;

    int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
    cout << "doJERShift: " << doJERShift << endl;

    int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
    cout << "dobTagEffShift: " << dobTagEffShift << endl;

    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;

    string postfix = "_Run2_TopTree_Study_" + dName; // to relabel the names of the output file

    if (doJESShift == 1)
        postfix= postfix+"_JESMinus";
    if (doJESShift == 2)
        postfix= postfix+"_JESPlus";
    if (doJERShift == 1)
        postfix= postfix+"_JERMinus";
    if (doJERShift == 2)
        postfix= postfix+"_JERPlus";
    if (dobTagEffShift == -1)
        postfix= postfix+"_bTagMinus";
    if (dobTagEffShift == 1)
        postfix= postfix+"_bTagPlus";
    if(domisTagEffShift == -1)
        postfix= postfix+"_misTagMinus";
    if(domisTagEffShift == 1)
        postfix= postfix+"_misTagPlus";


    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    bool debug = false;
    bool Muon = true;
    bool Electron = false;
    bool bTagReweight_PreReweighting = false; //Needs to be set only once to true in order to produce the EtaPtHistos
    string btagger = "CSVM";
    bool printTriggers = false;
    bool applyTriggers = true;
	  float Luminosity = 2094.087; //pb^-1 Muon  = 2196.422335, Electron = 2094.087
    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels
   	FILE* eventlist;

    if(Muon && !Electron)
    {
        cout << " --> Using the Muon channel..." << endl;
        channelpostfix = "_Mu";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
    	eventlist = fopen("EventInfo_mu.txt","w");
    }
    else if(!Muon && Electron)
    {
        cout << " --> Using the Electron channel..." << endl;
        channelpostfix = "_El";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
    	eventlist = fopen("EventInfo_El.txt","w");
    }
    else
    {
        cerr<<"Correct lepton Channel not selected."<<endl;
        exit(1);
    }

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

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
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = true;
    anaEnv.loadGenJetCollection = false;
    anaEnv.loadGenEventCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.loadTrackMETCollection = false;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;

    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);


    //////////////////////////////////////////////
    // Btag and lepton scale factors
    //////////////////////////////////////////////
    BTagCalibration * bTagCalib;   
//    BTagCalibrationReader * bTagReader_comb_central;
//    BTagCalibrationReader * bTagReader_comb_up;
//    BTagCalibrationReader * bTagReader_comb_down;
    BTagCalibrationReader * bTagReader_mujets_central;
//    BTagCalibrationReader * bTagReader_mujets_up;
//    BTagCalibrationReader * bTagReader_mujets_down;
//    BTagCalibrationReader * bTagReader_ttbar_central;
//    BTagCalibrationReader * bTagReader_ttbar_up;
//    BTagCalibrationReader * bTagReader_ttbar_down;
//    BTagWeightTools *btwt_comb_central;
//    BTagWeightTools *btwt_comb_up;
//    BTagWeightTools *btwt_comb_down;
    BTagWeightTools *btwt_mujets_central = 0;
//    BTagWeightTools *btwt_mujets_up;
//    BTagWeightTools *btwt_mujets_down;
//    BTagWeightTools *btwt_ttbar_central;
//    BTagWeightTools *btwt_ttbar_up;
//    BTagWeightTools *btwt_ttbar_down;

    if(bTagReweight)
    {
        if(dName.find("Data")==string::npos)        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
        {
            bTagCalib = new BTagCalibration("CSVv2","../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv");
            bTagReader_mujets_central = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","central"); //mujets
//            bTagReader_mujets_up = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","up"); //mujets
//            bTagReader_mujets_down = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","down"); //mujets
//            bTagReader_comb_central = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"comb","central"); //mujets
//            bTagReader_comb_up = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"comb","up"); //mujets
//            bTagReader_comb_down = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"comb","down"); //mujets
//            bTagReader_ttbar_central = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"ttbar","central"); //mujets
//            bTagReader_ttbar_up = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"ttbar","up"); //mujets
//            bTagReader_ttbar_down = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"ttbar","down"); //mujets
            if(bTagReweight_PreReweighting)// Need to differentiate BTagWeightTools according to filling the histos and just reading, because of overwriting possibilities in grid submission
            {
//                btwt_comb_central = new BTagWeightTools(bTagReader_comb_central,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum + "_comb_central.root",false,30,999,2.4);
//                btwt_comb_up = new BTagWeightTools(bTagReader_comb_up,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_comb_up.root",false,30,999,2.4);
//                btwt_comb_down = new BTagWeightTools(bTagReader_comb_down,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_comb_down.root",false,30,999,2.4);
                btwt_mujets_central = new BTagWeightTools(bTagReader_mujets_central,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_central.root",false,30,999,2.4);
//                btwt_mujets_up = new BTagWeightTools(bTagReader_mujets_up,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_up.root",false,30,999,2.4);
//                btwt_mujets_down = new BTagWeightTools(bTagReader_mujets_down,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_down.root",false,30,999,2.4);
//                btwt_ttbar_central = new BTagWeightTools(bTagReader_ttbar_central,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_ttbar_central.root",false,30,999,2.4);
//                btwt_ttbar_up = new BTagWeightTools(bTagReader_ttbar_up,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_ttbar_up.root",false,30,999,2.4);
//                btwt_ttbar_down = new BTagWeightTools(bTagReader_ttbar_down,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_ttbar_down.root",false,30,999,2.4);
            }
            else
            {
//                btwt_comb_central = new BTagWeightTools(bTagReader_comb_central,"BTagHistosPtEta/HistosPtEta_"+dName + "_comb_central.root",false,30,999,2.4);
//                btwt_comb_up = new BTagWeightTools(bTagReader_comb_up,"BTagHistosPtEta/HistosPtEta_"+dName +"_comb_up.root",false,30,999,2.4);
//                btwt_comb_down = new BTagWeightTools(bTagReader_comb_down,"BTagHistosPtEta/HistosPtEta_"+dName+"_comb_down.root",false,30,999,2.4);
					cout << "CAVEAT!!! Using the BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root as standard PtEta histo for b-tag reweighing" << endl;
                btwt_mujets_central = new BTagWeightTools(bTagReader_mujets_central,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
//                btwt_mujets_up = new BTagWeightTools(bTagReader_mujets_up,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_up.root",false,30,999,2.4);
//                btwt_mujets_down = new BTagWeightTools(bTagReader_mujets_down,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_down.root",false,30,999,2.4);
//                btwt_ttbar_central = new BTagWeightTools(bTagReader_ttbar_central,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_central.root",false,30,999,2.4);
//                btwt_ttbar_up = new BTagWeightTools(bTagReader_ttbar_up,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_up.root",false,30,999,2.4);
//                btwt_ttbar_down = new BTagWeightTools(bTagReader_ttbar_down,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_down.root",false,30,999,2.4);
            }
        }       
    }

    MuonSFWeight* muonSFWeight;   
    ElectronSFWeight* electronSFWeight; 
    if(bLeptonSF){
        if(Muon){
            muonSFWeight = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/Muon_SF_TopEA.root","SF_totErr",false,false);
        }
        else if(Electron){
            electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/Elec_SF_TopEA.root","GlobalSF",false,false);    
        }
    }

    LumiReWeighting LumiWeights;
//    LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert_norm.root", "pileup50", "pileup");    
    LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root", "pileup60", "pileup");    
//    LumiWeights = LumiReWeighting("/user/lbeck/CMSSW_7_6_0/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "/user/lbeck/CMSSW_7_6_0/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert_Silver.root", "pileup50", "pileup");    

    ////////////////////////////
    ///  Initialise trigger  ///
    ////////////////////////////

    if (debug) cout << "Initializing trigger" << endl;    
    //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);
    Trigger* trigger = 0;
    if(applyTriggers) trigger = new Trigger(Muon, Electron, true, false);

    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////
/*    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    if(dName.find("Data")==string::npos)
    {
        Luminosity = theDataset->EquivalentLumi();
        cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
*/    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    // add jobs number at the end of file if 
    //Output ROOT file
    SourceDate *strdate = new SourceDate();
    string date_str = strdate->ReturnDateStr();
    string histo_dir = "TreeMakerOutput/MACRO_histos"+channelpostfix;
    string histo_dir_date = histo_dir+"/MACRO_histos_" + date_str +"/";
    int mkdirstatus_histos = mkdir(histo_dir.c_str(),0777);
    mkdirstatus_histos = mkdir(histo_dir_date.c_str(),0777);
    
    string rootFileName (histo_dir_date+"/FCNC_1L3B_"+postfix+channelpostfix+".root");
    if (strJobNum != "0")
      {
	      cout << "strJobNum is " << strJobNum << endl;
        rootFileName = histo_dir_date+"/FCNC_1L3B_"+postfix+channelpostfix+"_"+strJobNum+".root";
      }
    
    
    

    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_fatjets;
    vector < TRootMET* >      mets;

    //Global variable
    TRootEvent* event = 0;
    TRootRun *runInfos = new TRootRun();

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////

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
/*
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
*/
    //B-tagging discriminators
    histo1D["Bdisc_CSV_jet1"]                             = new TH1F( "Bdisc_CSV_jet1", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2"]                             = new TH1F( "Bdisc_CSV_jet2", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3"]                             = new TH1F( "Bdisc_CSV_jet3", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet4"]                             = new TH1F( "Bdisc_CSV_jet4", "CSV b-disc._{jet4}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet5"]                             = new TH1F( "Bdisc_CSV_jet5", "CSV b-disc._{jet5}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet6"]                             = new TH1F( "Bdisc_CSV_jet6", "CSV b-disc._{jet6}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet1"]                             = new TH1F( "Bdisc_CSV_Bjet1", "CSV b-disc._{bjet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet2"]                             = new TH1F( "Bdisc_CSV_Bjet2", "CSV b-disc._{bjet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet3"]                             = new TH1F( "Bdisc_CSV_Bjet3", "CSV b-disc._{bjet3}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet4"]                             = new TH1F( "Bdisc_CSV_Bjet4", "CSV b-disc._{bjet4}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet5"]                             = new TH1F( "Bdisc_CSV_Bjet5", "CSV b-disc._{bjet5}", 30, 0, 1);
    histo1D["Bdisc_CSV_Bjet6"]                             = new TH1F( "Bdisc_CSV_Bjet6", "CSV b-disc._{bjet6}", 30, 0, 1);
    //Jets
    histo1D["JetEta"]                                        = new TH1F( "JetEta", "Jet #eta", 40,-4, 4);
    //0J0B
    histo1D["NbJets_0J0B"]                                        = new TH1F( "NbJets_0J0B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_0J0B"]                                        = new TH1F( "NbCSVLJets_0J0B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_0J0B"]                                        = new TH1F( "NbCSVMJets_0J0B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_0J0B"]                                        = new TH1F( "NbCSVTJets_0J0B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_0J0B"]                                        = new TH1F( "NbLightJets_0J0B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_0J0B"]                             = new TH1F( "Bdisc_CSV_jet1_0J0B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_0J0B"]                             = new TH1F( "Bdisc_CSV_jet2_0J0B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_0J0B"]                             = new TH1F( "Bdisc_CSV_jet3_0J0B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_0J0B"]                                        = new TH1F( "MuonPt_0J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_0J0B"]                                        = new TH1F( "ElectronPt_0J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_0J0B"]                                      = new TH1F( "1stJetPt_0J0B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_0J0B"]                                      = new TH1F( "2ndJetPt_0J0B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_0J0B"]                                      = new TH1F( "3rdJetPt_0J0B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_0J0B"]                                      = new TH1F( "4thJetPt_0J0B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_0J0B"]                               = new TH1F( "HT_SelectedJets_0J0B", "HT", 30, 0, 1500);
    histo1D["MET_0J0B"]                                           = new TH1F( "MET_0J0B", "MET", 70, 0, 700);
    //0J1B
    histo1D["NbJets_0J1B"]                                        = new TH1F( "NbJets_0J1B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_0J1B"]                                        = new TH1F( "NbCSVLJets_0J1B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_0J1B"]                                        = new TH1F( "NbCSVMJets_0J1B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_0J1B"]                                        = new TH1F( "NbCSVTJets_0J1B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_0J1B"]                                        = new TH1F( "NbLightJets_0J1B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_0J1B"]                             = new TH1F( "Bdisc_CSV_jet1_0J1B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_0J1B"]                             = new TH1F( "Bdisc_CSV_jet2_0J1B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_0J1B"]                             = new TH1F( "Bdisc_CSV_jet3_0J1B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_0J1B"]                                        = new TH1F( "MuonPt_0J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_0J1B"]                                        = new TH1F( "ElectronPt_0J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_0J1B"]                                      = new TH1F( "1stJetPt_0J1B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_0J1B"]                                      = new TH1F( "2ndJetPt_0J1B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_0J1B"]                                      = new TH1F( "3rdJetPt_0J1B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_0J1B"]                                      = new TH1F( "4thJetPt_0J1B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_0J1B"]                               = new TH1F( "HT_SelectedJets_0J1B", "HT", 30, 0, 1500);
    histo1D["MET_0J1B"]                                           = new TH1F( "MET_0J1B", "MET", 70, 0, 700);
    //0J2B
    histo1D["NbJets_0J2B"]                                        = new TH1F( "NbJets_0J2B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_0J2B"]                                        = new TH1F( "NbCSVLJets_0J2B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_0J2B"]                                        = new TH1F( "NbCSVMJets_0J2B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_0J2B"]                                        = new TH1F( "NbCSVTJets_0J2B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_0J2B"]                                        = new TH1F( "NbLightJets_0J2B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_0J2B"]                             = new TH1F( "Bdisc_CSV_jet1_0J2B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_0J2B"]                             = new TH1F( "Bdisc_CSV_jet2_0J2B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_0J2B"]                             = new TH1F( "Bdisc_CSV_jet3_0J2B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_0J2B"]                                        = new TH1F( "MuonPt_0J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_0J2B"]                                        = new TH1F( "ElectronPt_0J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_0J2B"]                                      = new TH1F( "1stJetPt_0J2B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_0J2B"]                                      = new TH1F( "2ndJetPt_0J2B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_0J2B"]                                      = new TH1F( "3rdJetPt_0J2B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_0J2B"]                                      = new TH1F( "4thJetPt_0J2B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_0J2B"]                               = new TH1F( "HT_SelectedJets_0J2B", "HT", 30, 0, 1500);
    histo1D["MET_0J2B"]                                           = new TH1F( "MET_0J2B", "MET", 70, 0, 700);
    //0J3B
    histo1D["NbJets_0J3B"]                                        = new TH1F( "NbJets_0J3B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_0J3B"]                                        = new TH1F( "NbCSVLJets_0J3B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_0J3B"]                                        = new TH1F( "NbCSVMJets_0J3B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_0J3B"]                                        = new TH1F( "NbCSVTJets_0J3B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_0J3B"]                                        = new TH1F( "NbLightJets_0J3B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_0J3B"]                             = new TH1F( "Bdisc_CSV_jet1_0J3B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_0J3B"]                             = new TH1F( "Bdisc_CSV_jet2_0J3B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_0J3B"]                             = new TH1F( "Bdisc_CSV_jet3_0J3B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_0J3B"]                                        = new TH1F( "MuonPt_0J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_0J3B"]                                        = new TH1F( "ElectronPt_0J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_0J3B"]                                      = new TH1F( "1stJetPt_0J3B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_0J3B"]                                      = new TH1F( "2ndJetPt_0J3B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_0J3B"]                                      = new TH1F( "3rdJetPt_0J3B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_0J3B"]                                      = new TH1F( "4thJetPt_0J3B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_0J3B"]                               = new TH1F( "HT_SelectedJets_0J3B", "HT", 30, 0, 1500);
    histo1D["MET_0J3B"]                                           = new TH1F( "MET_0J3B", "MET", 70, 0, 700);
    //1J0B
    histo1D["NbJets_1J0B"]                                        = new TH1F( "NbJets_1J0B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_1J0B"]                                        = new TH1F( "NbCSVLJets_1J0B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_1J0B"]                                        = new TH1F( "NbCSVMJets_1J0B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_1J0B"]                                        = new TH1F( "NbCSVTJets_1J0B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_1J0B"]                                        = new TH1F( "NbLightJets_1J0B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_1J0B"]                             = new TH1F( "Bdisc_CSV_jet1_1J0B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_1J0B"]                             = new TH1F( "Bdisc_CSV_jet2_1J0B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_1J0B"]                             = new TH1F( "Bdisc_CSV_jet3_1J0B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_1J0B"]                                        = new TH1F( "MuonPt_1J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_1J0B"]                                        = new TH1F( "ElectronPt_1J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_1J0B"]                                      = new TH1F( "1stJetPt_1J0B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_1J0B"]                                      = new TH1F( "2ndJetPt_1J0B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_1J0B"]                                      = new TH1F( "3rdJetPt_1J0B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_1J0B"]                                      = new TH1F( "4thJetPt_1J0B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_1J0B"]                               = new TH1F( "HT_SelectedJets_1J0B", "HT", 30, 0, 1500);
    histo1D["MET_1J0B"]                                           = new TH1F( "MET_1J0B", "MET", 70, 0, 700);
    //1J1B
    histo1D["NbJets_1J1B"]                                        = new TH1F( "NbJets_1J1B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_1J1B"]                                        = new TH1F( "NbCSVLJets_1J1B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_1J1B"]                                        = new TH1F( "NbCSVMJets_1J1B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_1J1B"]                                        = new TH1F( "NbCSVTJets_1J1B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_1J1B"]                                        = new TH1F( "NbLightJets_1J1B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_1J1B"]                             = new TH1F( "Bdisc_CSV_jet1_1J1B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_1J1B"]                             = new TH1F( "Bdisc_CSV_jet2_1J1B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_1J1B"]                             = new TH1F( "Bdisc_CSV_jet3_1J1B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_1J1B"]                                        = new TH1F( "MuonPt_1J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_1J1B"]                                        = new TH1F( "ElectronPt_1J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_1J1B"]                                      = new TH1F( "1stJetPt_1J1B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_1J1B"]                                      = new TH1F( "2ndJetPt_1J1B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_1J1B"]                                      = new TH1F( "3rdJetPt_1J1B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_1J1B"]                                      = new TH1F( "4thJetPt_1J1B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_1J1B"]                               = new TH1F( "HT_SelectedJets_1J1B", "HT", 30, 0, 1500);
    histo1D["MET_1J1B"]                                           = new TH1F( "MET_1J1B", "MET", 70, 0, 700);
    //1J2B
    histo1D["NbJets_1J2B"]                                        = new TH1F( "NbJets_1J2B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_1J2B"]                                        = new TH1F( "NbCSVLJets_1J2B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_1J2B"]                                        = new TH1F( "NbCSVMJets_1J2B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_1J2B"]                                        = new TH1F( "NbCSVTJets_1J2B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_1J2B"]                                        = new TH1F( "NbLightJets_1J2B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_1J2B"]                             = new TH1F( "Bdisc_CSV_jet1_1J2B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_1J2B"]                             = new TH1F( "Bdisc_CSV_jet2_1J2B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_1J2B"]                             = new TH1F( "Bdisc_CSV_jet3_1J2B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_1J2B"]                                        = new TH1F( "MuonPt_1J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_1J2B"]                                        = new TH1F( "ElectronPt_1J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_1J2B"]                                      = new TH1F( "1stJetPt_1J2B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_1J2B"]                                      = new TH1F( "2ndJetPt_1J2B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_1J2B"]                                      = new TH1F( "3rdJetPt_1J2B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_1J2B"]                                      = new TH1F( "4thJetPt_1J2B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_1J2B"]                               = new TH1F( "HT_SelectedJets_1J2B", "HT", 30, 0, 1500);
    histo1D["MET_1J2B"]                                           = new TH1F( "MET_1J2B", "MET", 70, 0, 700);
    //1J3B
    histo1D["NbJets_1J3B"]                                        = new TH1F( "NbJets_1J3B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_1J3B"]                                        = new TH1F( "NbCSVLJets_1J3B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_1J3B"]                                        = new TH1F( "NbCSVMJets_1J3B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_1J3B"]                                        = new TH1F( "NbCSVTJets_1J3B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_1J3B"]                                        = new TH1F( "NbLightJets_1J3B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_1J3B"]                             = new TH1F( "Bdisc_CSV_jet1_1J3B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_1J3B"]                             = new TH1F( "Bdisc_CSV_jet2_1J3B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_1J3B"]                             = new TH1F( "Bdisc_CSV_jet3_1J3B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_1J3B"]                                        = new TH1F( "MuonPt_1J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_1J3B"]                                        = new TH1F( "ElectronPt_1J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_1J3B"]                                      = new TH1F( "1stJetPt_1J3B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_1J3B"]                                      = new TH1F( "2ndJetPt_1J3B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_1J3B"]                                      = new TH1F( "3rdJetPt_1J3B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_1J3B"]                                      = new TH1F( "4thJetPt_1J3B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_1J3B"]                               = new TH1F( "HT_SelectedJets_1J3B", "HT", 30, 0, 1500);
    histo1D["MET_1J3B"]                                           = new TH1F( "MET_1J3B", "MET", 70, 0, 700);
    //2J0B
    histo1D["NbJets_2J0B"]                                        = new TH1F( "NbJets_2J0B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_2J0B"]                                        = new TH1F( "NbCSVLJets_2J0B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_2J0B"]                                        = new TH1F( "NbCSVMJets_2J0B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_2J0B"]                                        = new TH1F( "NbCSVTJets_2J0B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_2J0B"]                                        = new TH1F( "NbLightJets_2J0B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_2J0B"]                             = new TH1F( "Bdisc_CSV_jet1_2J0B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_2J0B"]                             = new TH1F( "Bdisc_CSV_jet2_2J0B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_2J0B"]                             = new TH1F( "Bdisc_CSV_jet3_2J0B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_2J0B"]                                        = new TH1F( "MuonPt_2J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_2J0B"]                                        = new TH1F( "ElectronPt_2J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_2J0B"]                                      = new TH1F( "1stJetPt_2J0B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_2J0B"]                                      = new TH1F( "2ndJetPt_2J0B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_2J0B"]                                      = new TH1F( "3rdJetPt_2J0B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_2J0B"]                                      = new TH1F( "4thJetPt_2J0B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_2J0B"]                               = new TH1F( "HT_SelectedJets_2J0B", "HT", 30, 0, 1500);
    histo1D["MET_2J0B"]                                           = new TH1F( "MET_2J0B", "MET", 70, 0, 700);
    //2J1B
    histo1D["NbJets_2J1B"]                                        = new TH1F( "NbJets_2J1B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_2J1B"]                                        = new TH1F( "NbCSVLJets_2J1B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_2J1B"]                                        = new TH1F( "NbCSVMJets_2J1B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_2J1B"]                                        = new TH1F( "NbCSVTJets_2J1B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_2J1B"]                                        = new TH1F( "NbLightJets_2J1B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_2J1B"]                             = new TH1F( "Bdisc_CSV_jet1_2J1B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_2J1B"]                             = new TH1F( "Bdisc_CSV_jet2_2J1B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_2J1B"]                             = new TH1F( "Bdisc_CSV_jet3_2J1B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_2J1B"]                                        = new TH1F( "MuonPt_2J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_2J1B"]                                        = new TH1F( "ElectronPt_2J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_2J1B"]                                      = new TH1F( "1stJetPt_2J1B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_2J1B"]                                      = new TH1F( "2ndJetPt_2J1B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_2J1B"]                                      = new TH1F( "3rdJetPt_2J1B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_2J1B"]                                      = new TH1F( "4thJetPt_2J1B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_2J1B"]                               = new TH1F( "HT_SelectedJets_2J1B", "HT", 30, 0, 1500);
    histo1D["MET_2J1B"]                                           = new TH1F( "MET_2J1B", "MET", 70, 0, 700);
    //2J2B
    histo1D["NbJets_2J2B"]                                        = new TH1F( "NbJets_2J2B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_2J2B"]                                        = new TH1F( "NbCSVLJets_2J2B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_2J2B"]                                        = new TH1F( "NbCSVMJets_2J2B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_2J2B"]                                        = new TH1F( "NbCSVTJets_2J2B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_2J2B"]                                        = new TH1F( "NbLightJets_2J2B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_2J2B"]                             = new TH1F( "Bdisc_CSV_jet1_2J2B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_2J2B"]                             = new TH1F( "Bdisc_CSV_jet2_2J2B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_2J2B"]                             = new TH1F( "Bdisc_CSV_jet3_2J2B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_2J2B"]                                        = new TH1F( "MuonPt_2J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_2J2B"]                                        = new TH1F( "ElectronPt_2J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_2J2B"]                                      = new TH1F( "1stJetPt_2J2B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_2J2B"]                                      = new TH1F( "2ndJetPt_2J2B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_2J2B"]                                      = new TH1F( "3rdJetPt_2J2B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_2J2B"]                                      = new TH1F( "4thJetPt_2J2B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_2J2B"]                               = new TH1F( "HT_SelectedJets_2J2B", "HT", 30, 0, 1500);
    histo1D["MET_2J2B"]                                           = new TH1F( "MET_2J2B", "MET", 70, 0, 700);
    //2J3B
    histo1D["NbJets_2J3B"]                                        = new TH1F( "NbJets_2J3B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_2J3B"]                                        = new TH1F( "NbCSVLJets_2J3B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_2J3B"]                                        = new TH1F( "NbCSVMJets_2J3B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_2J3B"]                                        = new TH1F( "NbCSVTJets_2J3B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_2J3B"]                                        = new TH1F( "NbLightJets_2J3B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_2J3B"]                             = new TH1F( "Bdisc_CSV_jet1_2J3B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_2J3B"]                             = new TH1F( "Bdisc_CSV_jet2_2J3B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_2J3B"]                             = new TH1F( "Bdisc_CSV_jet3_2J3B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_2J3B"]                                        = new TH1F( "MuonPt_2J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_2J3B"]                                        = new TH1F( "ElectronPt_2J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_2J3B"]                                      = new TH1F( "1stJetPt_2J3B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_2J3B"]                                      = new TH1F( "2ndJetPt_2J3B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_2J3B"]                                      = new TH1F( "3rdJetPt_2J3B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_2J3B"]                                      = new TH1F( "4thJetPt_2J3B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_2J3B"]                               = new TH1F( "HT_SelectedJets_2J3B", "HT", 30, 0, 1500);
    histo1D["MET_2J3B"]                                           = new TH1F( "MET_2J3B", "MET", 70, 0, 700);
    //3J0B
    histo1D["NbJets_3J0B"]                                        = new TH1F( "NbJets_3J0B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_3J0B"]                                        = new TH1F( "NbCSVLJets_3J0B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_3J0B"]                                        = new TH1F( "NbCSVMJets_3J0B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_3J0B"]                                        = new TH1F( "NbCSVTJets_3J0B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_3J0B"]                                        = new TH1F( "NbLightJets_3J0B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_3J0B"]                             = new TH1F( "Bdisc_CSV_jet1_3J0B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_3J0B"]                             = new TH1F( "Bdisc_CSV_jet2_3J0B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_3J0B"]                             = new TH1F( "Bdisc_CSV_jet3_3J0B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_3J0B"]                                        = new TH1F( "MuonPt_3J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_3J0B"]                                        = new TH1F( "ElectronPt_3J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_3J0B"]                                      = new TH1F( "1stJetPt_3J0B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_3J0B"]                                      = new TH1F( "2ndJetPt_3J0B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_3J0B"]                                      = new TH1F( "3rdJetPt_3J0B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_3J0B"]                                      = new TH1F( "4thJetPt_3J0B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_3J0B"]                               = new TH1F( "HT_SelectedJets_3J0B", "HT", 30, 0, 1500);
    histo1D["MET_3J0B"]                                           = new TH1F( "MET_3J0B", "MET", 70, 0, 700);
    //3J1B
    histo1D["NbJets_3J1B"]                                        = new TH1F( "NbJets_3J1B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_3J1B"]                                        = new TH1F( "NbCSVLJets_3J1B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_3J1B"]                                        = new TH1F( "NbCSVMJets_3J1B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_3J1B"]                                        = new TH1F( "NbCSVTJets_3J1B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_3J1B"]                                        = new TH1F( "NbLightJets_3J1B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_3J1B"]                             = new TH1F( "Bdisc_CSV_jet1_3J1B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_3J1B"]                             = new TH1F( "Bdisc_CSV_jet2_3J1B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_3J1B"]                             = new TH1F( "Bdisc_CSV_jet3_3J1B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_3J1B"]                                        = new TH1F( "MuonPt_3J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_3J1B"]                                        = new TH1F( "ElectronPt_3J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_3J1B"]                                      = new TH1F( "1stJetPt_3J1B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_3J1B"]                                      = new TH1F( "2ndJetPt_3J1B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_3J1B"]                                      = new TH1F( "3rdJetPt_3J1B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_3J1B"]                                      = new TH1F( "4thJetPt_3J1B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_3J1B"]                               = new TH1F( "HT_SelectedJets_3J1B", "HT", 30, 0, 1500);
    histo1D["MET_3J1B"]                                           = new TH1F( "MET_3J1B", "MET", 70, 0, 700);
    //3J2B
    histo1D["NbJets_3J2B"]                                        = new TH1F( "NbJets_3J2B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_3J2B"]                                        = new TH1F( "NbCSVLJets_3J2B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_3J2B"]                                        = new TH1F( "NbCSVMJets_3J2B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_3J2B"]                                        = new TH1F( "NbCSVTJets_3J2B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_3J2B"]                                        = new TH1F( "NbLightJets_3J2B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_3J2B"]                             = new TH1F( "Bdisc_CSV_jet1_3J2B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_3J2B"]                             = new TH1F( "Bdisc_CSV_jet2_3J2B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_3J2B"]                             = new TH1F( "Bdisc_CSV_jet3_3J2B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_3J2B"]                                        = new TH1F( "MuonPt_3J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_3J2B"]                                        = new TH1F( "ElectronPt_3J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_3J2B"]                                      = new TH1F( "1stJetPt_3J2B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_3J2B"]                                      = new TH1F( "2ndJetPt_3J2B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_3J2B"]                                      = new TH1F( "3rdJetPt_3J2B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_3J2B"]                                      = new TH1F( "4thJetPt_3J2B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_3J2B"]                               = new TH1F( "HT_SelectedJets_3J2B", "HT", 30, 0, 1500);
    histo1D["MET_3J2B"]                                           = new TH1F( "MET_3J2B", "MET", 70, 0, 700);
    //3J3B
    histo1D["NbJets_3J3B"]                                        = new TH1F( "NbJets_3J3B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_3J3B"]                                        = new TH1F( "NbCSVLJets_3J3B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_3J3B"]                                        = new TH1F( "NbCSVMJets_3J3B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_3J3B"]                                        = new TH1F( "NbCSVTJets_3J3B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_3J3B"]                                        = new TH1F( "NbLightJets_3J3B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_3J3B"]                             = new TH1F( "Bdisc_CSV_jet1_3J3B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_3J3B"]                             = new TH1F( "Bdisc_CSV_jet2_3J3B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_3J3B"]                             = new TH1F( "Bdisc_CSV_jet3_3J3B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_3J3B"]                                        = new TH1F( "MuonPt_3J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_3J3B"]                                        = new TH1F( "ElectronPt_3J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_3J3B"]                                      = new TH1F( "1stJetPt_3J3B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_3J3B"]                                      = new TH1F( "2ndJetPt_3J3B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_3J3B"]                                      = new TH1F( "3rdJetPt_3J3B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_3J3B"]                                      = new TH1F( "4thJetPt_3J3B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_3J3B"]                               = new TH1F( "HT_SelectedJets_3J3B", "HT", 30, 0, 1500);
    histo1D["MET_3J3B"]                                           = new TH1F( "MET_3J3B", "MET", 70, 0, 700);
    //4J0B
    histo1D["NbJets_4J0B"]                                        = new TH1F( "NbJets_4J0B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_4J0B"]                                        = new TH1F( "NbCSVLJets_4J0B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_4J0B"]                                        = new TH1F( "NbCSVMJets_4J0B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_4J0B"]                                        = new TH1F( "NbCSVTJets_4J0B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_4J0B"]                                        = new TH1F( "NbLightJets_4J0B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_4J0B"]                             = new TH1F( "Bdisc_CSV_jet1_4J0B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_4J0B"]                             = new TH1F( "Bdisc_CSV_jet2_4J0B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_4J0B"]                             = new TH1F( "Bdisc_CSV_jet3_4J0B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_4J0B"]                                        = new TH1F( "MuonPt_4J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_4J0B"]                                        = new TH1F( "ElectronPt_4J0B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_4J0B"]                                      = new TH1F( "1stJetPt_4J0B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_4J0B"]                                      = new TH1F( "2ndJetPt_4J0B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_4J0B"]                                      = new TH1F( "3rdJetPt_4J0B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_4J0B"]                                      = new TH1F( "4thJetPt_4J0B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_4J0B"]                               = new TH1F( "HT_SelectedJets_4J0B", "HT", 30, 0, 1500);
    histo1D["MET_4J0B"]                                           = new TH1F( "MET_4J0B", "MET", 70, 0, 700);
    //4J1B
    histo1D["NbJets_4J1B"]                                        = new TH1F( "NbJets_4J1B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_4J1B"]                                        = new TH1F( "NbCSVLJets_4J1B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_4J1B"]                                        = new TH1F( "NbCSVMJets_4J1B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_4J1B"]                                        = new TH1F( "NbCSVTJets_4J1B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_4J1B"]                                        = new TH1F( "NbLightJets_4J1B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_4J1B"]                             = new TH1F( "Bdisc_CSV_jet1_4J1B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_4J1B"]                             = new TH1F( "Bdisc_CSV_jet2_4J1B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_4J1B"]                             = new TH1F( "Bdisc_CSV_jet3_4J1B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_4J1B"]                                        = new TH1F( "MuonPt_4J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_4J1B"]                                        = new TH1F( "ElectronPt_4J1B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_4J1B"]                                      = new TH1F( "1stJetPt_4J1B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_4J1B"]                                      = new TH1F( "2ndJetPt_4J1B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_4J1B"]                                      = new TH1F( "3rdJetPt_4J1B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_4J1B"]                                      = new TH1F( "4thJetPt_4J1B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_4J1B"]                               = new TH1F( "HT_SelectedJets_4J1B", "HT", 30, 0, 1500);
    histo1D["MET_4J1B"]                                           = new TH1F( "MET_4J1B", "MET", 70, 0, 700);
    //4J2B
    histo1D["NbJets_4J2B"]                                        = new TH1F( "NbJets_4J2B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_4J2B"]                                        = new TH1F( "NbCSVLJets_4J2B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_4J2B"]                                        = new TH1F( "NbCSVMJets_4J2B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_4J2B"]                                        = new TH1F( "NbCSVTJets_4J2B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_4J2B"]                                        = new TH1F( "NbLightJets_4J2B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_4J2B"]                             = new TH1F( "Bdisc_CSV_jet1_4J2B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_4J2B"]                             = new TH1F( "Bdisc_CSV_jet2_4J2B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_4J2B"]                             = new TH1F( "Bdisc_CSV_jet3_4J2B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_4J2B"]                                        = new TH1F( "MuonPt_4J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_4J2B"]                                        = new TH1F( "ElectronPt_4J2B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_4J2B"]                                      = new TH1F( "1stJetPt_4J2B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_4J2B"]                                      = new TH1F( "2ndJetPt_4J2B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_4J2B"]                                      = new TH1F( "3rdJetPt_4J2B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_4J2B"]                                      = new TH1F( "4thJetPt_4J2B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_4J2B"]                               = new TH1F( "HT_SelectedJets_4J2B", "HT", 30, 0, 1500);
    histo1D["MET_4J2B"]                                           = new TH1F( "MET_4J2B", "MET", 70, 0, 700);
    //4J3B
    histo1D["NbJets_4J3B"]                                        = new TH1F( "NbJets_4J3B", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets_4J3B"]                                        = new TH1F( "NbCSVLJets_4J3B", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets_4J3B"]                                        = new TH1F( "NbCSVMJets_4J3B", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets_4J3B"]                                        = new TH1F( "NbCSVTJets_4J3B", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["NbLightJets_4J3B"]                                        = new TH1F( "NbLightJets_4J3B", "nb. Light tags", 15,-0.5, 14.5);
    histo1D["Bdisc_CSV_jet1_4J3B"]                             = new TH1F( "Bdisc_CSV_jet1_4J3B", "CSV b-disc._{jet1}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet2_4J3B"]                             = new TH1F( "Bdisc_CSV_jet2_4J3B", "CSV b-disc._{jet2}", 30, 0, 1);
    histo1D["Bdisc_CSV_jet3_4J3B"]                             = new TH1F( "Bdisc_CSV_jet3_4J3B", "CSV b-disc._{jet3}", 30, 0, 1);
    histo1D["MuonPt_4J3B"]                                        = new TH1F( "MuonPt_4J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["ElectronPt_4J3B"]                                        = new TH1F( "ElectronPt_4J3B", "PT_{#mu}", 30, 0, 300);
    histo1D["1stJetPt_4J3B"]                                      = new TH1F( "1stJetPt_4J3B", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt_4J3B"]                                      = new TH1F( "2ndJetPt_4J3B", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt_4J3B"]                                      = new TH1F( "3rdJetPt_4J3B", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt_4J3B"]                                      = new TH1F( "4thJetPt_4J3B", "PT_{jet4}", 30, 0, 300);
    histo1D["HT_SelectedJets_4J3B"]                               = new TH1F( "HT_SelectedJets_4J3B", "HT", 30, 0, 1500);
    histo1D["MET_4J3B"]                                           = new TH1F( "MET_4J3B", "MET", 70, 0, 700);
    histo1D["NbJets"]                                        = new TH1F( "NbJets", "nb. jets", 15,-0.5, 14.5);
    histo1D["NbCSVLJets"]                                        = new TH1F( "NbCSVLJets", "nb. CSVL tags", 15,-0.5, 14.5);
    histo1D["NbCSVMJets"]                                        = new TH1F( "NbCSVMJets", "nb. CSVM tags", 15,-0.5, 14.5);
    histo1D["NbCSVTJets"]                                        = new TH1F( "NbCSVTJets", "nb. CSVT tags", 15,-0.5, 14.5);
    histo1D["1stJetPt"]                                      = new TH1F( "1stJetPt", "PT_{jet1}", 30, 0, 300);
    histo1D["2ndJetPt"]                                      = new TH1F( "2ndJetPt", "PT_{jet2}", 30, 0, 300);
    histo1D["3rdJetPt"]                                      = new TH1F( "3rdJetPt", "PT_{jet3}", 30, 0, 300);
    histo1D["4thJetPt"]                                      = new TH1F( "4thJetPt", "PT_{jet4}", 30, 0, 300);
    histo1D["5thJetPt"]                                      = new TH1F( "5thJetPt", "PT_{jet5}", 30, 0, 300);
    histo1D["6thJetPt"]                                      = new TH1F( "6thJetPt", "PT_{jet6}", 30, 0, 300);
    histo1D["1stBJetPt"]                                      = new TH1F( "1stBJetPt", "PT_{bjet1}", 30, 0, 300);
    histo1D["2ndBJetPt"]                                      = new TH1F( "2ndBJetPt", "PT_{bjet2}", 30, 0, 300);
    histo1D["3rdBJetPt"]                                      = new TH1F( "3rdBJetPt", "PT_{bjet3}", 30, 0, 300);
    histo1D["4thBJetPt"]                                      = new TH1F( "4thBJetPt", "PT_{bjet4}", 30, 0, 300);
    histo1D["5thBJetPt"]                                      = new TH1F( "5thBJetPt", "PT_{bjet5}", 30, 0, 300);
    histo1D["6thBJetPt"]                                      = new TH1F( "6thBJetPt", "PT_{bjet6}", 30, 0, 300);
    histo1D["HT_SelectedJets"]                               = new TH1F( "HT_SelectedJets", "HT", 30, 0, 1500);
    //MET
    histo1D["MET_preCut"]                                           = new TH1F( "MET_preCut", "MET", 70, 0, 700);
    histo1D["MT_LepMET_preCut"]                                           = new TH1F( "MET_LepMET_preCut", "MT(lep,MET)", 70, 0, 700);
    histo1D["MET"]                                           = new TH1F( "MET", "MET", 70, 0, 700);
    histo1D["MT_LepMET"]                                           = new TH1F( "MT_LepMET", "MT(lep,MET)", 70, 0, 700);
	//MC-info plots
    histo1D["JetID_Hmother"]                                           = new TH1F( "JetID_Hmother", "jet number", 12,-0.5,11.5);

    ///////////////////
    // 2D histograms //
    ///////////////////
    histo2D["NJet_vs_Nbjet"] = new TH2F("NJet_vs_Nbjet","NJet:Nbjet",12,-0.5,11.5, 61, -0.5,11.5);
    histo2D["JetID_vs_pdgID"] = new TH2F("JetID_vs_pdgID","parton pdgID:jet number",12,-0.5,11.5, 61, -30.5,30.5);



    /////////////////////////////////
    //       Loop on datasets      //
    /////////////////////////////////
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    string btagHisto_str = "BTagHistosPtEta/HistosPtEta"+dName+"_mujets_central.root";
    TFile * btaghistos = new TFile(btagHisto_str.c_str());

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset

        bool nlo = false;
        bool isData = false;


        if(dName.find("NLO") != std::string::npos || dName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;

        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;
        float normfactor = datasets[d]->NormFactor();

        ////////////////////////////////////////////////////////////
        // Setup Date string and nTuple for output  
        ///////////////////////////////////////////////////////////

//        SourceDate *strdate = new SourceDate();
//        string date_str = strdate->ReturnDateStr();
        if(debug)cout<<"date print"<<endl;

        string channel_dir = "TreeMakerOutput/Ntuples"+channelpostfix;
        string date_dir = channel_dir+"/Ntuples_" + date_str +"/";
        int mkdirstatus__ = mkdir(channel_dir.c_str(),0777);
        mkdirstatus__ = mkdir(date_dir.c_str(),0777);

        string jobNumString = static_cast<ostringstream*>( &(ostringstream() << JobNum) )->str();
        string Ntupname = date_dir +"FCNC_1L3B_" +postfix + channelpostfix + "_" + jobNumString + ".root";
        string Ntuptitle_AdvancedVars = "AdvancedVarsTree";
        string Ntuptitle_ObjectVars = "ObjectVarsTree";
        string Ntuptitle_EventInfo = "EventInfoTree";
        string Ntuptitle_Weights = "Weights";

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
        TNtuple * tup      = new TNtuple(Ntuptitle_AdvancedVars.c_str(), Ntuptitle_AdvancedVars.c_str(), "MTlepmet:MLepTop_GenMatch:MHadTop_GenMatch:EtaLepTop_GenMatch:EtaHadTop_GenMatch:MassW_GenMatch:EtaW_GenMatch:dR_lepJet_min:MLepTop:MHadTop:EtaLepTop:EtaHadTop:MassW:EtaW");
        TNtuple * tup_ObjectVars      = new TNtuple(Ntuptitle_ObjectVars.c_str(), Ntuptitle_ObjectVars.c_str(), "qlepton:leptonpt:leptoneta:leptonX:leptonY:leptonZ:leptonE:bdisc1:bdisc2:bdisc3:bdisc4:bdisc5:jet1_Pt:jet2_Pt:jet3_Pt:jet4_Pt:jet5_Pt:jet1_Eta:jet2_Eta:jet3_Eta:jet4_Eta:jet5_Eta:jet1_x:jet2_x:jet3_x:jet4_x:jet5_x:jet1_y:jet2_y:jet3_y:jet4_y:jet5_y:jet1_z:jet2_z:jet3_z:jet4_z:jet5_z:jet1_E:jet2_E:jet3_E:jet4_E:jet5_E:MissingEt");
        TNtuple * tup_EventInfo      = new TNtuple(Ntuptitle_EventInfo.c_str(), Ntuptitle_EventInfo.c_str(), "nbVertices:nb_jets:nb_bjets");
        TNtuple * tup_Weights      = new TNtuple(Ntuptitle_Weights.c_str(), Ntuptitle_Weights.c_str(), "lumiWeight:fleptonSF:btagWeight_comb_central:btagWeight_comb_up:btagWeight_comb_down:btagWeight_mujets_central:btagWeight_mujets_up:btagWeight_mujets_down:btagWeight_ttbar_central:btagWeight_ttbar_up:btagWeight_ttbar_down");
                    
        
        
        if(debug)cout<<"created craneens"<<endl;
        ///////////////////////////////////////////////////////////////
        // JEC
        ///////////////////////////////////////////////////////////////
        vector<JetCorrectorParameters> vCorrParam;

        if(dName.find("Data")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
        }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer15_25nsV6_DATA_Uncertainty_AK4PFchs.txt");

        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);

        /// book triggers
        if (applyTriggers) { trigger->bookTriggers(isData);}


        //////////////////////////////////////////////////
        // Pre-event loop definitions
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1, start = 0;
        int currentRun;
        int iFile = -1;
        int isdata_int = 0;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;
        string previousFilename = "";

        int event_start = startEvent;

        if (debug) cout << " - Loop over events " << endl;

        float qlepton,leptonpt,bdisc1,bdisc2,bdisc3,bdisc4,bdisc5,nb_jets,nb_bjets,jet1_Pt,jet2_Pt,jet3_Pt,MissingEt,leptoneta,MTlepmet,MLepTop_GenMatch,MHadTop_GenMatch,EtaLepTop_GenMatch,EtaHadTop_GenMatch,MassW_GenMatch,EtaW_GenMatch,dR_lepJet_min,MLepTop,MHadTop,EtaLepTop,EtaHadTop,MassW,EtaW,nbVertices;
        float leptonX,leptonY,leptonZ,leptonE,jet4_Pt,jet5_Pt,jet1_Eta,jet2_Eta,jet3_Eta,jet4_Eta,jet5_Eta,jet1_x,jet2_x,jet3_x,jet4_x,jet5_x,jet1_y,jet2_y,jet3_y,jet4_y,jet5_y,jet1_z,jet2_z,jet3_z,jet4_z,jet5_z,jet1_E,jet2_E,jet3_E,jet4_E,jet5_E;
        float lumiWeight, fleptonSF;
        float btagWeight_comb_central,btagWeight_comb_up,btagWeight_comb_down,btagWeight_mujets_central,btagWeight_mujets_up,btagWeight_mujets_down,btagWeight_ttbar_central,btagWeight_ttbar_up,btagWeight_ttbar_down;


        double end_d = ending;
        if(endEvent > ending) end_d = ending;
        else end_d = endEvent;

        double EqLumi = 1.;
        if(dName.find("Data")==string::npos) EqLumi = end_d/xSect;
        cout <<"Equivalent Lumi: "<<  EqLumi  << endl;
        cout <<"Will run over "<<  end_d<< " events..."<<endl;    cout <<"Starting event = = = = "<< event_start  << endl;


        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
	        qlepton = 0; leptonpt = -1.; bdisc1 = -1.; bdisc2 = -1.; bdisc3 = -1.; bdisc4 = -1.; bdisc5 = -1.; nb_jets = -1; nb_bjets = -1.; jet1_Pt = -1; jet2_Pt = -1.;jet3_Pt = -1.; MissingEt = -1.;
            leptoneta = -99999.; MTlepmet = -1.; MLepTop_GenMatch = -1.; MHadTop_GenMatch = -1.; EtaLepTop_GenMatch = -99999.;
            EtaHadTop_GenMatch = -99999.; MassW_GenMatch = -1.; EtaW_GenMatch = -99999.; dR_lepJet_min = 99999.;
            MLepTop = -1.; MHadTop = -1.; EtaLepTop = -99999.; EtaHadTop = -99999.; MassW = -1.; EtaW = -99999.; nbVertices = -1;
            jet4_Pt = -1.; jet5_Pt = -1.; jet1_Eta = -99999.; jet2_Eta = -99999.; jet3_Eta = -99999.; jet4_Eta = -99999.; jet5_Eta = -99999.;
            jet1_x = -99999.; jet2_x = -99999.; jet3_x = -99999.; jet4_x = -99999.; jet5_x = -99999.;
            jet1_y = -99999.; jet2_y = -99999.; jet3_y = -99999.; jet4_y = -99999.; jet5_y = -99999.;
            jet1_z = -99999.; jet2_z = -99999.; jet3_z = -99999.; jet4_z = -99999.; jet5_z = -99999.; jet1_E = -1.; jet2_E = -1.; jet3_E = -1.; jet4_E = -1.; jet5_E = -1.;
            leptonX = -99999.; leptonY = -99999.; leptonZ = -99999.; leptonE = -1.;
            lumiWeight = 1.; fleptonSF = 1.;
            btagWeight_comb_central = 1.;btagWeight_comb_up = 1.;btagWeight_comb_down = 1.;btagWeight_mujets_central = 1.;
            btagWeight_mujets_up = 1.;btagWeight_mujets_down = 1.;btagWeight_ttbar_central = 1.;btagWeight_ttbar_up = 1.;btagWeight_ttbar_down = 1.;

            double ievt_d = ievt;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%10000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            if(debug)cout<<"before tree load"<<endl;
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            if(debug)cout<<"after tree load"<<endl;

            ////////////////////////////
            ///  Include trigger set up here when using data
            ////////////////////////////
            datasets[d]->eventTree()->LoadTree(ievt);
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            int currentRun = event->runId();

            bool fileChanged = false;
            bool runChanged = false;
            if( event->runId() > 10000) isdata_int=1;


            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;
            //////////////////////////////////
            //Loading Gen jets //
            //////////////////////////////////

            vector<TRootGenJet*> genjets;
            if(dName.find("Data")==string::npos)//genjets only available in non-data samples
            {
                if(debug) cout << "loading genjets" << endl;
                // loading GenJets as I need them for JER
                genjets = treeLoader.LoadGenJet(ievt);
            }

            ///////////////////////
            // JER smearing
            //////////////////////

            if(dName.find("Data")==string::npos)//JER smearing only feasible for non-data samples
            {
                //JER
                doJERShift == 0;
                if(doJERShift == 1)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
                else if(doJERShift == 2)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
                else
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");

                //     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");

                // JES sysematic!
                if (doJESShift == 1)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
                else if (doJESShift == 2)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "plus");

                //            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");

            }
            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            //define object containers
            vector<TRootElectron*> selectedElectrons;
            vector<TRootPFJet*>    selectedJets;
            vector<TRootMuon*>     selectedMuons;
            vector<TRootElectron*> selectedExtraElectrons;
            vector<TRootMuon*>     selectedExtraMuons;
            selectedElectrons.reserve(10);
            selectedMuons.reserve(10);
            vector<TRootJet*>      selectedLBJets;
            vector<TRootJet*>      selectedMBJets;
            vector<TRootJet*>      selectedTBJets;
            vector<TRootJet*>      selectedLightJets_LWP;
            vector<TRootJet*>      selectedLightJets_MWP;
            vector<TRootJet*>      selectedLightJets_TWP;
            vector<TRootJet*>      selectedLightJets;
		        vector<TLorentzVector> selectedMuonsTLV, selectedElectronsTLV, metsTLV, selectedJetsTLV, selectedBJetsTLV, selectedLeptonsTLV;
            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC, selectedLooseIsoMuonsTLV;
            vector<TLorentzVector> mcParticlesTLV, mcMuonsTLV, mcPartonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_,mcParticles;
            vector<int> mcMuonIndex, mcPartonIndex;
            JetPartonMatching muonMatching, jetMatching;



            // Declare selection instance
            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets);

            // Define object selection cuts

            int nMu, nEl, nLooseMu, nLooseEl; //number of (loose) muons/electrons
            nMu = 0, nEl = 0, nLooseMu=0, nLooseEl=0;

            /////////////////////////////////////////////
            // Define object selection cuts
            /////////////////////////////////////////////
            if (Muon)
            {
				if (debug)cout<<"Getting Jets"<<endl;
				selectedJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Tight"); // ApplyJetId
				if (debug)cout<<"Getting Tight Muons"<<endl;
				selectedMuons                                       = r2selection.GetSelectedMuons(30,2.4,0.15, "Tight", "Spring15"); //Selected
				if (debug)cout<<"Getting Loose Electrons"<<endl;
				selectedElectrons                                   = r2selection.GetSelectedElectrons(20,2.5,"Loose", "Spring15_25ns", true); //Vetoed  
				if (debug)cout<<"Getting Loose Muons"<<endl;
				selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20,"Loose","Spring15"); //Vetoed         
            }
            if (Electron)
            {
				if (debug)cout<<"Getting Jets"<<endl;
				selectedJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Tight"); // ApplyJetId
				if (debug)cout<<"Getting Loose Muons"<<endl;
				selectedMuons                                       = r2selection.GetSelectedMuons(20, 2.4, 0.20,"Loose","Spring15"); //Vetoed
				if (debug)cout<<"Getting Medium Electrons"<<endl;
				selectedElectrons                                   = r2selection.GetSelectedElectrons(30,2.5,"Medium", "Spring15_25ns", true); //Selected                       
				if (debug)cout<<"Getting Loose Electrons"<<endl;
				selectedExtraElectrons                              = r2selection.GetSelectedElectrons(20,2.5,"Loose", "Spring15_25ns", true); //Vetoed
            }
            if(Muon)
            {
              		nMu = selectedMuons.size(); //Number of Muons in Event (Tight only)
              		nEl = selectedElectrons.size(); //Number of Electrons in Event   (tight and loose)
              		nLooseMu = selectedExtraMuons.size();   //Number of loose muons      (loose only)
            }
            else if(Electron)
            {
              		nMu = selectedMuons.size(); //Number of Muons in Event (tight and loose)
              		nEl = selectedElectrons.size(); //Number of Electrons in Event (Tight only)
              		nLooseEl = selectedExtraElectrons.size(); //Number of loose electrons  (loose only)
            }



            ////////////////////////////////////////////////
            // PV selection and trigger
            ////////////////////////////////////////////////
            
            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV <<endl;
//            if (debug) cin.get();

            histo1D["cutFlow"]->Fill(0., 1. );

 //           weightCount += scaleFactor;
            eventCount++;

            bool trigged = false;
            if ( ! applyTriggers && previousFilename != currentFilename )
            {
              fileChanged = true;
              previousFilename = currentFilename;
              iFile++;
              cout << "File changed!!! => iFile = " << iFile << endl;
              trigged = true;
            }
            
            if (applyTriggers)
            {
              trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTriggers);
              trigged = trigger->checkIfFired();
              
            }
             if(dName.find("NP")!=string::npos)//JER smearing only feasible for non-data samples
            {
                trigged = true;
            }
            //////////////////////////////////////////////////////
            // Applying baseline lepton selection
            //////////////////////////////////////////////////////


          if (!isGoodPV) continue; // Check that there is a good Primary Vertex
          histo1D["cutFlow"]->Fill(1., 1. );//goodPV

          if(!trigged) continue;
          histo1D["cutFlow"]->Fill(2., 1. );//trigger

           if (debug)
          {
          	cout <<" applying baseline event selection..."<<endl;
          	cout <<"number of muons: " << nMu <<endl;
          	cout <<"number of electrons: " << nEl <<endl;
          }
          //Apply the lepton, btag and HT selections
          if (Muon && !Electron)
          {
            if  (  !( nMu ==1 && nEl == 0)) continue; // Muon Channel Selection
            if (selectedMuons[0]->Pt() < 30) continue;
            if (debug)	cout <<"Muon selection passed..."<<endl;
          }
          else if (!Muon && Electron)
          {
              if  (  !( nMu == 0 && nEl == 1)) continue; // Electron Channel Selection
              if (selectedElectrons[0]->Pt() < 30) continue;
           if (debug)	cout <<"Electron selection passed..."<<endl;
          }
          else
          {
              cerr<<"Correct Channel not selected."<<endl;
              exit(1);
          }
          histo1D["cutFlow"]->Fill(3., 1. );//1lepton

			if(Muon && !Electron)
			{
				if(nLooseMu != 1) continue;
	            if (debug)	cout <<"Vetoed extra muons..."<<endl;
			}
			if(!Muon && Electron)
			{
				if(nLooseEl != 1) continue;
	            if (debug)	cout <<"Vetoed extra electrons..."<<endl;
			}
			histo1D["cutFlow"]->Fill(4., 1. ); //LooseLepton removal
			
			////////////////////////////////////////////////////////////////////////////////
			// Clean jet collection from jets overlapping with lepton
      ////////////////////////////////////////////////////////////////////////////////
			float clean_dR = 0.4;
			if(Muon)
			{
			    for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			    {
			        float tmp_dR = 5.;
			        tmp_dR = sqrt(pow(selectedMuons[0]->Phi()-selectedJets[iJet]->Phi(),2)+pow(selectedMuons[0]->Eta()-selectedJets[iJet]->Eta(),2));
			        if(tmp_dR < clean_dR)
			        {
			            selectedJets.erase(selectedJets.begin() + iJet);
			            if(debug) cout << "Cleaned jet collection from overlapping lepton" << endl;
			        }
			    }
			}
			if(Electron)
			{
			    for(unsigned int iJet = 0; iJet < selectedJets.size(); iJet++)
			    {
			        float tmp_dR = 5.;
			        tmp_dR = sqrt(pow(selectedElectrons[0]->Phi()-selectedJets[iJet]->Phi(),2)+pow(selectedElectrons[0]->Eta()-selectedJets[iJet]->Eta(),2));
			        if(tmp_dR < clean_dR)
			        {
			            selectedJets.erase(selectedJets.begin() + iJet);
			            if(debug) cout << "Cleaned jet collection from overlapping lepton" << endl;
			        }
			    }
			}
			histo1D["cutFlow"]->Fill(5., 1. ); // JetCleaning
			
			/////////////////////////////////////////////
			// Make TLorentzVectors //
			////////////////////////////////////////////
			for(unsigned int iMuon=0; iMuon<selectedMuons.size(); iMuon++)
       		{
				selectedMuonsTLV.push_back( *selectedMuons[iMuon]);			
			}
			for(unsigned int iElectron=0; iElectron<selectedElectrons.size(); iElectron++)
       		{
				selectedElectronsTLV.push_back( *selectedElectrons[iElectron]);		
			}
			if(mets.size() > 0) metsTLV.push_back( *mets[0] );
			else cout<<"No MET??"<<endl;
			
			////////////////////////////////////////////////////////////////
			//Calculations for MT(lep,MET) selection cut
			////////////////////////////////////////////////////////////////
			float MT = -999;
			if(Muon)
			{
			    MT = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
			    selectedLeptonsTLV.push_back(*selectedMuons[0]);
			}
			else if(Electron)
			{
			    MT = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
			    selectedLeptonsTLV.push_back(*selectedElectrons[0]);
			}
			else cout << "Wrong channel (1)" << endl;

			histo1D["MT_LepMET_preCut"] ->Fill(MT, Luminosity/EqLumi*scaleFactor );
			histo1D["MET_preCut"] ->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor );

//			if(MT <= 50) continue;
//			histo1D["cutFlow"]->Fill(???., 1. );
//	        if (debug)	cout <<"Cut on MT..."<<endl;

            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output
			////////////////////////////////////
		    //Fill b-jet collections
            ////////////////////////////////////
		    for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
		    {
		     	if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose   )
		        {
		        	selectedLBJets.push_back(selectedJets[seljet]);
		            if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium)
		            {
		            	selectedMBJets.push_back(selectedJets[seljet]);
		                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight)
		                {
		                	selectedTBJets.push_back(selectedJets[seljet]);
		                }
		   				else selectedLightJets_TWP.push_back(selectedJets[seljet]);
		        	}
		   			else
		   			{
		   				selectedLightJets_MWP.push_back(selectedJets[seljet]);
		   				selectedLightJets_TWP.push_back(selectedJets[seljet]);
		   			}
		   		}
		   		else
		   		{
		   			selectedLightJets_LWP.push_back(selectedJets[seljet]);
		   			selectedLightJets_MWP.push_back(selectedJets[seljet]);
		   			selectedLightJets_TWP.push_back(selectedJets[seljet]);
		   		}
		  	}
		  	
      	histo2D["NJet_vs_Nbjet"]->Fill(selectedJets.size(),selectedMBJets.size());

            //////////////////////////////////////////////////
            // Apply scale factors
            //////////////////////////////////////////////////
            if(dName.find("Data")!=string::npos) //If sample is data, no PU reweighting
            {
                lumiWeight=1;
            }
            else
            {
                lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
            }
            if(debug) cout << "lumiweight: " << lumiWeight << endl;
            //cout << "lumiweight: " << lumiWeight << ", nTruePU: " << (int)event->nTruePU() << endl;
            if(bLeptonSF){
                if(Muon){
                    fleptonSF = muonSFWeight->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                }
                else if(Electron){
                    fleptonSF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                }
            }
            if(debug) cout<<"lepton SF:  "<<fleptonSF<<endl;
            if(bTagReweight && !bTagReweight_PreReweighting)
            {
                if(dName.find("Data")==string::npos) //If sample is data, no PU reweighting
                {
                    if(debug) cout << "Applying b-tag weights " << endl;
//                    btagWeight_comb_central =  btwt_comb_central->getMCEventWeight(selectedJets, false);
//                    btagWeight_comb_up =  btwt_comb_up->getMCEventWeight(selectedJets, false);
//                    btagWeight_comb_down =  btwt_comb_down->getMCEventWeight(selectedJets, false);
                    btagWeight_mujets_central =  btwt_mujets_central->getMCEventWeight(selectedJets, false);
//                    btagWeight_mujets_up =  btwt_mujets_up->getMCEventWeight(selectedJets, false);
//                    btagWeight_mujets_down =  btwt_mujets_down->getMCEventWeight(selectedJets, false);
//                    btagWeight_ttbar_central =  btwt_ttbar_central->getMCEventWeight(selectedJets, false);
//                    btagWeight_ttbar_up =  btwt_ttbar_up->getMCEventWeight(selectedJets, false);
//                    btagWeight_ttbar_down =  btwt_ttbar_down->getMCEventWeight(selectedJets, false);
                }
            }
            if(debug) cout<<"btag SF:  "<< endl;

            scaleFactor = scaleFactor * lumiWeight * fleptonSF;




			//////////////////////////////////////
			// Pre-jet histograms //
			/////////////////////////////////////
			float bdisc_jet1=-1., bdisc_jet2=-1., bdisc_jet3=-1.;
			float jet1Pt=-1., jet2Pt=-1., jet3Pt=-1., jet4Pt=-1.;
			if(selectedJets.size()>=1){ bdisc_jet1 = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();jet1Pt = selectedJets[0]->Pt();}
			if(selectedJets.size()>=2){ bdisc_jet2 = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();jet2Pt = selectedJets[1]->Pt();}
			if(selectedJets.size()>=3){ bdisc_jet3 = selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags();jet3Pt = selectedJets[2]->Pt();}
			if(selectedJets.size()>=4){ jet4Pt = selectedJets[3]->Pt();}
   		float HT = 0, H = 0;


   		for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
   		{
        		//Event-level variables
        		HT = HT + selectedJets[seljet1]->Pt();
        		H = H +  selectedJets[seljet1]->P();
   		}



			if(selectedJets.size() >= 0 && selectedMBJets.size() >=0)
			{
				histo1D["NbJets_0J0B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_0J0B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_0J0B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_0J0B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_0J0B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_0J0B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_0J0B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_0J0B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_0J0B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_0J0B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_0J0B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_0J0B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_0J0B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_0J0B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_0J0B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_0J0B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 0 && selectedMBJets.size() >=1)
			{
				histo1D["NbJets_0J1B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_0J1B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_0J1B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_0J1B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_0J1B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_0J1B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_0J1B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_0J1B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_0J1B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_0J1B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_0J1B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_0J1B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_0J1B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_0J1B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_0J1B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_0J1B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 0 && selectedMBJets.size() >=2)
			{
				histo1D["NbJets_0J2B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_0J2B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_0J2B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_0J2B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_0J2B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_0J2B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_0J2B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_0J2B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_0J2B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_0J2B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_0J2B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_0J2B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_0J2B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_0J2B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_0J2B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_0J2B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 0 && selectedMBJets.size() >=3)
			{
				histo1D["NbJets_0J3B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_0J3B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_0J3B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_0J3B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_0J3B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_0J3B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_0J3B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_0J3B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_0J3B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_0J3B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_0J3B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_0J3B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_0J3B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_0J3B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_0J3B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_0J3B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 1 && selectedMBJets.size() >=0)
			{
				histo1D["NbJets_1J0B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_1J0B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_1J0B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_1J0B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_1J0B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_1J0B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_1J0B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_1J0B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_1J0B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_1J0B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_1J0B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_1J0B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_1J0B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_1J0B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_1J0B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_1J0B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 1 && selectedMBJets.size() >=1)
			{
				histo1D["NbJets_1J1B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_1J1B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_1J1B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_1J1B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_1J1B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_1J1B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_1J1B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_1J1B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_1J1B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_1J1B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_1J1B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_1J1B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_1J1B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_1J1B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_1J1B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_1J1B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 1 && selectedMBJets.size() >=2)
			{
				histo1D["NbJets_1J2B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_1J2B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_1J2B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_1J2B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_1J2B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_1J2B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_1J2B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_1J2B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_1J2B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_1J2B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_1J2B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_1J2B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_1J2B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_1J2B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_1J2B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_1J2B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 1 && selectedMBJets.size() >=3)
			{
				histo1D["NbJets_1J3B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_1J3B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_1J3B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_1J3B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_1J3B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_1J3B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_1J3B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_1J3B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_1J3B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_1J3B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_1J3B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_1J3B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_1J3B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_1J3B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_1J3B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_1J3B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 2 && selectedMBJets.size() >=0)
			{
				histo1D["NbJets_2J0B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_2J0B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_2J0B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_2J0B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_2J0B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_2J0B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_2J0B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_2J0B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_2J0B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_2J0B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_2J0B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_2J0B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_2J0B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_2J0B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_2J0B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_2J0B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 2 && selectedMBJets.size() >=1)
			{
				histo1D["NbJets_2J1B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_2J1B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_2J1B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_2J1B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_2J1B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_2J1B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_2J1B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_2J1B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_2J1B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_2J1B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_2J1B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_2J1B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_2J1B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_2J1B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_2J1B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_2J1B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 2 && selectedMBJets.size() >=2)
			{
				histo1D["NbJets_2J2B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_2J2B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_2J2B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_2J2B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_2J2B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_2J2B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_2J2B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_2J2B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_2J2B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_2J2B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_2J2B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_2J2B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_2J2B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_2J2B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_2J2B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_2J2B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 2 && selectedMBJets.size() >=3)
			{
				histo1D["NbJets_2J3B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_2J3B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_2J3B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_2J3B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_2J3B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_2J3B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_2J3B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_2J3B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_2J3B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_2J3B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_2J3B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_2J3B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_2J3B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_2J3B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_2J3B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_2J3B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 3 && selectedMBJets.size() >=0)
			{
				histo1D["NbJets_3J0B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_3J0B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_3J0B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_3J0B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_3J0B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_3J0B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_3J0B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_3J0B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_3J0B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_3J0B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_3J0B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_3J0B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_3J0B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_3J0B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_3J0B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_3J0B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 3 && selectedMBJets.size() >=1)
			{
				histo1D["NbJets_3J1B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_3J1B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_3J1B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_3J1B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_3J1B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_3J1B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_3J1B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_3J1B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_3J1B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_3J1B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_3J1B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_3J1B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_3J1B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_3J1B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_3J1B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_3J1B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 3 && selectedMBJets.size() >=2)
			{
				histo1D["NbJets_3J2B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_3J2B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_3J2B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_3J2B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_3J2B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_3J2B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_3J2B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_3J2B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_3J2B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_3J2B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_3J2B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_3J2B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_3J2B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_3J2B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_3J2B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_3J2B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 3 && selectedMBJets.size() >=3)
			{
				histo1D["NbJets_3J3B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_3J3B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_3J3B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_3J3B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_3J3B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_3J3B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_3J3B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_3J3B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_3J3B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_3J3B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_3J3B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_3J3B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_3J3B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_3J3B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_3J3B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_3J3B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 4 && selectedMBJets.size() >=0)
			{
				histo1D["NbJets_4J0B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_4J0B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_4J0B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_4J0B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_4J0B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_4J0B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_4J0B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_4J0B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_4J0B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_4J0B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_4J0B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_4J0B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_4J0B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_4J0B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_4J0B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_4J0B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 4 && selectedMBJets.size() >=1)
			{
				histo1D["NbJets_4J1B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_4J1B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_4J1B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_4J1B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_4J1B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_4J1B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_4J1B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_4J1B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_4J1B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_4J1B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_4J1B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_4J1B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_4J1B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_4J1B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_4J1B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_4J1B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 4 && selectedMBJets.size() >=2)
			{
				histo1D["NbJets_4J2B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_4J2B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_4J2B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_4J2B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_4J2B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_4J2B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_4J2B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_4J2B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_4J2B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_4J2B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_4J2B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_4J2B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_4J2B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_4J2B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_4J2B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_4J2B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 4 && selectedMBJets.size() >=3)
			{
				histo1D["NbJets_4J3B"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVLJets_4J3B"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVMJets_4J3B"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbCSVTJets_4J3B"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);
				histo1D["NbLightJets_4J3B"] ->Fill(selectedLightJets_MWP.size(), Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet1_4J3B"]->Fill(bdisc_jet1, Luminosity/EqLumi*scaleFactor);
    			histo1D["Bdisc_CSV_jet2_4J3B"]->Fill(bdisc_jet2, Luminosity/EqLumi*scaleFactor);
		    	histo1D["Bdisc_CSV_jet3_4J3B"]->Fill(bdisc_jet3, Luminosity/EqLumi*scaleFactor);
				if(Electron) histo1D["ElectronPt_4J3B"]->Fill(selectedElectrons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
				if(Muon) histo1D["MuonPt_4J3B"]->Fill(selectedMuons[0]->Pt(),  Luminosity/EqLumi*scaleFactor);
    			histo1D["1stJetPt_4J3B"]->Fill(jet1Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["2ndJetPt_4J3B"]->Fill(jet2Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["3rdJetPt_4J3B"]->Fill(jet3Pt, Luminosity/EqLumi*scaleFactor);
    			histo1D["4thJetPt_4J3B"]->Fill(jet4Pt, Luminosity/EqLumi*scaleFactor);
				histo1D["MET_4J3B"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["HT_SelectedJets_4J3B"]->Fill(HT, Luminosity/EqLumi*scaleFactor);
			}

            //////////////////////////////////////////////////////////////////////
            // Cut on nb of jets and b-jets
            //////////////////////////////////////////////////////////////////////
			if(selectedJets.size() < 3)  continue;
			histo1D["cutFlow"]->Fill(6., 1. );//n Jets
	        if (debug)	cout <<"Cut on nb jets..."<<endl;

            ///////////////////////////////////////////////////
            // Fill b-tag histos for scale factors
            // info: http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v2.pdf
            ////////////////////////////////////////////////////
            if(bTagReweight_PreReweighting)
            {
                if(dName.find("Data")==string::npos)        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
                {
//                    btwt_comb_central->FillMCEfficiencyHistos(selectedJets);
//                    btwt_comb_up->FillMCEfficiencyHistos(selectedJets);
//                    btwt_comb_down->FillMCEfficiencyHistos(selectedJets);
                    btwt_mujets_central->FillMCEfficiencyHistos(selectedJets);
//                    btwt_mujets_up->FillMCEfficiencyHistos(selectedJets);
//                    btwt_mujets_down->FillMCEfficiencyHistos(selectedJets);
//                    btwt_ttbar_central->FillMCEfficiencyHistos(selectedJets);
//                    btwt_ttbar_up->FillMCEfficiencyHistos(selectedJets);
//                    btwt_ttbar_down->FillMCEfficiencyHistos(selectedJets);
                }
            }       





		  	if(selectedMBJets.size() < 3) continue;
	        if (debug)	cout <<"Cut on nb b-jets..."<<endl;
			histo1D["cutFlow"]->Fill(7., 1. ); //n BJets

			histo1D["cutFlow"]->Fill(8., 1. );

            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
            for (unsigned int i = 0; i < selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
            for (unsigned int i = 0; i < selectedMBJets.size(); i++) selectedBJetsTLV.push_back(*selectedMBJets[i]);

            //////////////////////////////////////
            // Peeking at the MC info 
            /////////////////////////////////////
            pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);// First one is jet number, second one is mcParticle number
            pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
                  
            int pdgID_top = 6; //top quark
                  
            bool Posleptonmatched = false;
            bool Negleptonmatched = false;
            
            if(dName != "data" && dName != "Data" && dName != "Data" && dName != "D_ata")
            {
                treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles, false);
                sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
                  
                mcParticlesMatching_.clear();
                    
                    
                for (unsigned int i = 0; i < mcParticles.size(); i++)
                {
                    if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process                   
                      
                    if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == 13 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )		// mu-, W-, tbar
                    {
                        Negleptonmatched = true;
                    }
                    else if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == -13 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )		// mu-, W-, tbar
                    {
                        Posleptonmatched = true;
                    }
                    else if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == 11 && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -pdgID_top )		// mu-, W-, tbar
                    {
                        Negleptonmatched = true;
                    }
                    else if ( mcParticles[i]->status() == 1 && mcParticles[i]->type() == -11 && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == pdgID_top )		// mu-, W-, tbar
                    {
                        Posleptonmatched = true;
                    }
                      
                    if ( abs(mcParticles[i]->type()) < 6 || abs(mcParticles[i]->type()) == 21 )  //light/b quarks, 6 should stay hardcoded, OR gluon
                    {
                        mcParticlesTLV.push_back(*mcParticles[i]);
                        mcParticlesMatching_.push_back(mcParticles[i]);
                    }
                      
                }
                    
                    
                JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
                if (matching.getNumberOfAvailableCombinations() != 1) cerr << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<" .  This should be equal to 1 !!!"<<endl;
                    
                    
                vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
                    
                for (unsigned int i = 0; i < mcParticlesTLV.size(); i++)
                {
                    int matchedJetNumber = matching.getMatchForParton(i, 0);
                    if (matchedJetNumber > -1) JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
                }
                    
                  
                for (unsigned int i = 0; i < JetPartonPair.size(); i++)
                {
                    unsigned int j = JetPartonPair[i].second;
                      
                    if ( fabs(mcParticlesMatching_[j]->type()) < 6 )  //light/b quarks, 6 should stay hardcoded
                    {
                        if ( ( Posleptonmatched && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -pdgID_top )
                               || ( Negleptonmatched && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == pdgID_top ) )  // if mu+, check if mother of particle is W- and granny tbar --> then it is a quark from W- decay
                        {
                            if (hadronicWJet1_.first == 9999)
                            {
                                hadronicWJet1_ = JetPartonPair[i];
                            }
                            else if (hadronicWJet2_.first == 9999)
                            {
                                hadronicWJet2_ = JetPartonPair[i];
                            }
                            else
                            {
                                cerr << "Found a third jet coming from a W boson which comes from a top quark..." << endl;
                                cerr << " -- pdgId: " << mcParticlesMatching_[j]->type() << " mother: " << mcParticlesMatching_[j]->motherType() << " granny: " << mcParticlesMatching_[j]->grannyType() << " Pt: " << mcParticlesMatching_[j]->Pt() << endl;
                                cerr << " -- ievt: " << ievt << endl;
                                exit(1);
                            }
                        }
                    }
                    if ( fabs(mcParticlesMatching_[j]->type()) == 5 )
                    {
                        if ( ( Posleptonmatched && mcParticlesMatching_[j]->motherType() == -pdgID_top )
                            || ( Negleptonmatched && mcParticlesMatching_[j]->motherType() == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
                        {
                            hadronicBJet_ = JetPartonPair[i];
                        }
                        else if ( ( Posleptonmatched && mcParticlesMatching_[j]->motherType() == pdgID_top )
                                        || ( Negleptonmatched && mcParticlesMatching_[j]->motherType() == -pdgID_top ) )
                        {
                            leptonicBJet_ = JetPartonPair[i];
                        }
                    }
                }  /// End loop over Jet Parton Pairs
            
            }// End MC matching (end of data-if-statement)




            ///////////////////////////////////////////////////
            // Filling histograms / plotting //
            //////////////////////////////////////////////////

            histo1D["NbOfVertices"]->Fill(vertex.size(), Luminosity/EqLumi*scaleFactor);




            /////////////////////////////////////
            // Muon Based Plots //
            /////////////////////////////////////
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                histo1D["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), Luminosity/EqLumi*scaleFactor);
                histo1D["LeptonPt"]->Fill(selectedMuons[selmu]->Pt(), Luminosity/EqLumi*scaleFactor);
                float reliso = selectedMuons[selmu]->relPfIso(4, 0.5);
                histo1D["MuonRelIsolation"]->Fill(reliso, Luminosity/EqLumi*scaleFactor);
            }

            //////////////////////////////////////////
            // Electron Based Plots //
            /////////////////////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                histo1D["ElectronRelIsolation"]->Fill(reliso, Luminosity/EqLumi*scaleFactor);
                histo1D["ElectronPt"]->Fill(selectedElectrons[selel]->Pt(), Luminosity/EqLumi*scaleFactor);
                histo1D["LeptonPt"]->Fill(selectedElectrons[selel]->Pt(), Luminosity/EqLumi*scaleFactor);
            }

            //////////////////////////////////
            // Jets Based Plots //
            //////////////////////////////////
			histo1D["NbJets"]->Fill(selectedJets.size(), Luminosity/EqLumi*scaleFactor);
			histo1D["NbCSVLJets"]->Fill(selectedLBJets.size(), Luminosity/EqLumi*scaleFactor);
			histo1D["NbCSVMJets"]->Fill(selectedMBJets.size(), Luminosity/EqLumi*scaleFactor);
			histo1D["NbCSVTJets"]->Fill(selectedTBJets.size(), Luminosity/EqLumi*scaleFactor);

			if(selectedJets.size() >= 1)
			{
				histo1D["Bdisc_CSV_jet1"]->Fill(selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["1stJetPt"]->Fill(selectedJets[0]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 2)
			{
				histo1D["Bdisc_CSV_jet2"]->Fill(selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["2ndJetPt"]->Fill(selectedJets[1]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 3)
			{
				histo1D["Bdisc_CSV_jet3"]->Fill(selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["3rdJetPt"]->Fill(selectedJets[2]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 4)
			{
				histo1D["Bdisc_CSV_jet4"]->Fill(selectedJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["4thJetPt"]->Fill(selectedJets[3]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 5)
			{
				histo1D["Bdisc_CSV_jet5"]->Fill(selectedJets[4]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["5thJetPt"]->Fill(selectedJets[4]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedJets.size() >= 6)
			{
				histo1D["Bdisc_CSV_jet6"]->Fill(selectedJets[5]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["6thJetPt"]->Fill(selectedJets[5]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			//B-jets
			if(debug) cout << "selectedMBJets.size() = "<< selectedMBJets.size() << endl;
			if(selectedMBJets.size() >= 1)
			{
				histo1D["Bdisc_CSV_Bjet1"]->Fill(selectedMBJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["1stBJetPt"]->Fill(selectedMBJets[0]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedMBJets.size() >= 2)
			{
				histo1D["Bdisc_CSV_Bjet2"]->Fill(selectedMBJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["2ndBJetPt"]->Fill(selectedMBJets[1]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedMBJets.size() >= 3)
			{
				histo1D["Bdisc_CSV_Bjet3"]->Fill(selectedMBJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["3rdBJetPt"]->Fill(selectedMBJets[2]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedMBJets.size() >= 4)
			{
				histo1D["Bdisc_CSV_Bjet4"]->Fill(selectedMBJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["4thBJetPt"]->Fill(selectedMBJets[3]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedMBJets.size() >= 5)
			{
				histo1D["Bdisc_CSV_Bjet5"]->Fill(selectedMBJets[4]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["5thBJetPt"]->Fill(selectedMBJets[4]->Pt(), Luminosity/EqLumi*scaleFactor);
			}
			if(selectedMBJets.size() >= 6)
			{
				histo1D["Bdisc_CSV_Bjet6"]->Fill(selectedMBJets[5]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), Luminosity/EqLumi*scaleFactor);
				histo1D["6thBJetPt"]->Fill(selectedMBJets[5]->Pt(), Luminosity/EqLumi*scaleFactor);
			}

			for(unsigned int i = 0; i < selectedJets.size(); i++)
			{
				histo1D["JetEta"]->Fill(selectedJets[i]->Eta(), Luminosity/EqLumi*scaleFactor);
			}


            histo1D["HT_SelectedJets"]->Fill(HT, Luminosity/EqLumi*scaleFactor);


            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();

            /////////////////////////////////
            //MET Based Plots//
            /////////////////////////////////

            histo1D["MET"]->Fill(mets[0]->Et(), Luminosity/EqLumi*scaleFactor);
            histo1D["MT_LepMET"]->Fill(MT, Luminosity/EqLumi*scaleFactor);

            /////////////////////////////////////////////////////////
            //Topology Reconstructions (truth level)
            ////////////////////////////////////////////////////////
            if( int(leptonicBJet_.first) != 9999) //Leptonic decays ///////////////////// CHANGE TO LEPTONIC TOP
            {
                int ind_jet = leptonicBJet_.first;///////////////////// CHANGE TO LEPTONIC TOP
                MLepTop_GenMatch = (selectedJetsTLV[ind_jet] + selectedLeptonsTLV[0] + metsTLV[0]).M();
                EtaLepTop_GenMatch = (selectedJetsTLV[ind_jet] + selectedLeptonsTLV[0] + metsTLV[0]).Eta();
            }
            if(int(hadronicWJet1_.first) != 9999 && int(hadronicWJet2_.first) != 9999) //Hadronic decays
            {
                int ind_jet1 = hadronicWJet1_.first;
                int ind_jet2 = hadronicWJet2_.first;
                MassW_GenMatch = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2]).M();
                EtaW_GenMatch = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2]).Eta();

                if(int(hadronicBJet_.first) != 9999)/////////////////////// CHANGE TO HADRONIC TOP
                {
                    int ind_jet_top = hadronicBJet_.first;///////////////////// CHANGE TO LEPTONIC TOP
                    MHadTop_GenMatch = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2] + selectedJetsTLV[ind_jet_top]).M();
                    EtaHadTop_GenMatch = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2] + selectedJetsTLV[ind_jet_top]).Eta();
                }
                
            }
            
            ///////////////////////////
            //Defining ntuple variables
            ///////////////////////////
            float chi_LepTop = 99999.;float chi_HadTop = 99999.;float chi_W = 99999.;float chi_comb = 99999.;
            float TopMass = 172.5;
            float WMass = 80.4;
            float Topsigma = 40;//ref https://indico.cern.ch/event/366968/session/9/contribution/10/attachments/729442/1000907/kskovpenFCNC20150530_tH.pdf
            float Wsigma = 15;//ref https://indico.cern.ch/event/366968/session/9/contribution/10/attachments/729442/1000907/kskovpenFCNC20150530_tH.pdf
            float Hsigma = 30;//ref https://indico.cern.ch/event/366968/session/9/contribution/10/attachments/729442/1000907/kskovpenFCNC20150530_tH.pdf
            //Reconstruction of top objects
            for(unsigned int  iBJet1 = 0; iBJet1 < selectedBJetsTLV.size(); iBJet1++)
            {
                float tmp_LepTopMass = -1.;float tmp_HadTopMass = -1.;float tmp_WMass = -1.;float tmp_LepTopEta = -99999.;float tmp_HadTopEta = -99999.;float tmp_WEta = -99999.;
                float chi_lep = 9999.;
                float chi_W = 9999.;
                float chi_hadTop = 9999.;
                for(unsigned int  iJet1 = 0; iJet1 < selectedJetsTLV.size(); iJet1++)
                {
                    if(selectedJetsTLV[iJet1].Pt() == selectedBJetsTLV[iBJet1].Pt() && selectedJetsTLV[iJet1].Eta() == selectedBJetsTLV[iBJet1].Eta() )
                    {
                        if(debug) cout << "Skipped same jet as in b-jet coll." << endl;
                        continue;
                    }
                    for(unsigned int  iJet2 = 0; iJet2 < selectedJetsTLV.size(); iJet2++)
                    {
                        if(selectedJetsTLV[iJet2].Pt() == selectedBJetsTLV[iBJet1].Pt() && selectedJetsTLV[iJet2].Eta() == selectedBJetsTLV[iBJet1].Eta() ) continue;
                        if(iJet1 == iJet2) continue;
                        
                        tmp_WMass = (selectedJetsTLV[iJet1]+selectedJetsTLV[iJet2]).M();
                        tmp_HadTopMass = (selectedJetsTLV[iJet1]+selectedJetsTLV[iJet2]+selectedBJetsTLV[iBJet1]).M();
                        tmp_WEta = (selectedJetsTLV[iJet1]+selectedJetsTLV[iJet2]).Eta();
                        tmp_HadTopEta = (selectedJetsTLV[iJet1]+selectedJetsTLV[iJet2]+selectedBJetsTLV[iBJet1]).Eta();
                        
                        chi_hadTop = pow(tmp_HadTopMass-TopMass,2)/pow(Topsigma,2);
                        chi_W = pow(tmp_WMass-WMass,2)/pow(Wsigma,2);
                        if(chi_hadTop + chi_W < chi_comb)
                        {
                            chi_comb = chi_hadTop + chi_W;
                            MHadTop = tmp_HadTopMass;
                            EtaHadTop = tmp_HadTopEta;
                            MassW = tmp_WMass;
                            EtaW = tmp_WEta;
                        }
                        
                    }//Loop iJet2
                }//Loop iJet1
            }//Loop iBJet1
            for(unsigned int  iJet = 0; iJet < selectedJetsTLV.size(); iJet++)
            {
                float tmp_dR_lepJet = 99999.;
                if(Muon) tmp_dR_lepJet = selectedMuonsTLV[0].DeltaR(selectedJetsTLV[iJet]);
                if(Electron) tmp_dR_lepJet= selectedElectronsTLV[0].DeltaR(selectedJetsTLV[iJet]);
                
                if(tmp_dR_lepJet < dR_lepJet_min) dR_lepJet_min = tmp_dR_lepJet;
            }
			if(Muon && !Electron)
			{
			    leptonpt = selectedMuons[0]->Pt();
			    leptoneta = selectedMuons[0]->Eta();
			    leptonX = selectedMuons[0]->X();
			    leptonY = selectedMuons[0]->Y();
			    leptonZ = selectedMuons[0]->Z();
			    leptonE = selectedMuons[0]->E();
			    qlepton = selectedMuons[0]->charge();
			}
			else if(!Muon && Electron)
			{
			    leptonpt = selectedElectrons[0]->Pt();
			    leptoneta = selectedElectrons[0]->Eta();
			    leptonX = selectedElectrons[0]->X();
			    leptonY = selectedElectrons[0]->Y();
			    leptonZ = selectedElectrons[0]->Z();
			    leptonE = selectedElectrons[0]->E();
			    qlepton = selectedElectrons[0]->charge();
			}
			bdisc1 = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			bdisc2 = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			if(selectedJets.size() >= 3)bdisc3 = selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			nb_jets = selectedJets.size();
			nb_bjets = selectedMBJets.size();
			jet1_Pt = selectedJets[0]->Pt();
			jet2_Pt = selectedJets[1]->Pt();
			if(selectedJets.size() >= 3)			jet3_Pt = selectedJets[2]->Pt();
			jet1_Eta = selectedJets[0]->Eta();
			jet2_Eta = selectedJets[1]->Eta();
			if(selectedJets.size() >= 3)			jet3_Eta = selectedJets[2]->Eta();
			jet1_x = selectedJets[0]->X();
			jet2_x = selectedJets[1]->X();
			if(selectedJets.size() >= 3)			jet3_x = selectedJets[2]->X();
			jet1_y = selectedJets[0]->Y();
			jet2_y = selectedJets[1]->Y();
			if(selectedJets.size() >= 3)			jet3_y = selectedJets[2]->Y();
			jet1_z = selectedJets[0]->Z();
			jet2_z = selectedJets[1]->Z();
			if(selectedJets.size() >= 3)			jet3_z = selectedJets[2]->Z();
			jet1_E = selectedJets[0]->E();
			jet2_E = selectedJets[1]->E();
			if(selectedJets.size() >= 3)			jet3_E = selectedJets[2]->E();
			if(selectedJets.size() >= 4)
			{
			    jet4_Pt = selectedJets[3]->Pt();
			    jet4_Eta = selectedJets[3]->Eta();
			    jet4_x = selectedJets[3]->X();
			    jet4_y = selectedJets[3]->Y();
			    jet4_z = selectedJets[3]->Z();
			    jet4_E = selectedJets[3]->E();
			    bdisc4 = selectedJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			}
			if(selectedJets.size() >= 5)
			{
			    jet5_Pt = selectedJets[4]->Pt();
			    jet5_Eta = selectedJets[4]->Eta();
			    jet5_x = selectedJets[4]->X();
			    jet5_y = selectedJets[4]->Y();
			    jet5_z = selectedJets[4]->Z();
			    jet5_E = selectedJets[4]->E();
			    bdisc5 = selectedJets[4]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			}
			MissingEt = mets[0]->Et();
			MTlepmet = MT;
			nbVertices = vertex.size();



            float vals[14] = {MTlepmet,MLepTop_GenMatch,MHadTop_GenMatch,EtaLepTop_GenMatch,EtaHadTop_GenMatch,MassW_GenMatch,EtaW_GenMatch,dR_lepJet_min,MLepTop,MHadTop,EtaLepTop,EtaHadTop,MassW,EtaW};
            float vals_ObjectVars[43] = {qlepton,leptonpt,leptoneta,leptonX,leptonY,leptonZ,leptonE,bdisc1,bdisc2,bdisc3,bdisc4,bdisc5,jet1_Pt,jet2_Pt,jet3_Pt,jet4_Pt,jet5_Pt,jet1_Eta,jet2_Eta,jet3_Eta,jet4_Eta,jet5_Eta,jet1_x,jet2_x,jet3_x,jet4_x,jet5_x,jet1_y,jet2_y,jet3_y,jet4_y,jet5_y,jet1_z,jet2_z,jet3_z,jet4_z,jet5_z,jet1_E,jet2_E,jet3_E,jet4_E,jet5_E,MissingEt};
            float vals_EventInfo[3] = {nbVertices,nb_jets,nb_bjets};
            float vals_Weights[11] = {lumiWeight,fleptonSF,btagWeight_comb_central,btagWeight_comb_up,btagWeight_comb_down,btagWeight_mujets_central,btagWeight_mujets_up,btagWeight_mujets_down,btagWeight_ttbar_central,btagWeight_ttbar_up,btagWeight_ttbar_down};
            tup->Fill(vals);
            tup_ObjectVars->Fill(vals_ObjectVars);
            tup_EventInfo->Fill(vals_EventInfo);
            tup_Weights->Fill(vals_Weights);


			///////////////////////////
			// Event info //////
			///////////////////////////
			if(Muon && !Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedMuons[0]->type(), selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), selectedMuons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());
			else if(!Muon && Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->type(), selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), selectedElectrons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());


        } //End Loop on Events

		tup->Write();
		tup_ObjectVars->Write();
		tup_EventInfo->Write();
       	tup_Weights->Write();
       	tupfile->Close();
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
//        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    fclose (eventlist);

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;

  fout->cd();

/* 	TDirectory* _0J0Bdir = fout->mkdir("0J0B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _0J1Bdir = fout->mkdir("0J1B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _0J2Bdir = fout->mkdir("0J2B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _0J3Bdir = fout->mkdir("0J3B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _1J0Bdir = fout->mkdir("1J0B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _1J1Bdir = fout->mkdir("1J1B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _1J2Bdir = fout->mkdir("1J2B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _1J3Bdir = fout->mkdir("1J3B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _2J0Bdir = fout->mkdir("2J0B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _2J1Bdir = fout->mkdir("2J1B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _2J2Bdir = fout->mkdir("2J2B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _2J3Bdir = fout->mkdir("2J3B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _3J0Bdir = fout->mkdir("3J0B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _3J1Bdir = fout->mkdir("3J1B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _3J2Bdir = fout->mkdir("3J2B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _3J3Bdir = fout->mkdir("3J3B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _4J0Bdir = fout->mkdir("4J0B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _4J1Bdir = fout->mkdir("4J1B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _4J2Bdir = fout->mkdir("4J2B"); //Makes different subdirectories in outputfile for every variable
 	TDirectory* _4J3Bdir = fout->mkdir("4J3B"); //Makes different subdirectories in outputfile for every variable
*/

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
//    delete btwt_comb_central;
//    delete btwt_comb_up;
//    delete btwt_comb_down;
    delete btwt_mujets_central;
//    delete btwt_mujets_up;
//    delete btwt_mujets_down;
//    delete btwt_ttbar_central;
//    delete btwt_ttbar_up;
//    delete btwt_ttbar_down;

/*    for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {

        TH2F *temp = it->second;
        temp->Write();
    }
*/
    fout->Write();   
    fout->Close();
    delete fout;




    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}


