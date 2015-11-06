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
#include "TopTreeAnalysisBase/Tools/interface/SourceDate.h"

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


// General flags
bool debug = false;
bool Muon = false;
bool Electron = true;
string btagger = "CSVM";

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;
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
    bool bTagReweight = false;
    bool bLeptonSF = false;
    
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
    BTagCalibrationReader * bTagReader;
    BTagWeightTools *btwt;
    // BTagCalibrationReader * bTagReadercomb;

    if(bTagReweight && dName.find("Data")==string::npos){
        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
        bTagCalib = new BTagCalibration("CSVv2","CSVv2_13TeV_combTomujets.csv");
        bTagReader = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","central"); //mujets
        btwt = new BTagWeightTools("HistosPtEta_"+dName+".root",true);        
        // bTagReadercomb = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"comb","central"); //comb
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


    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////
	float Luminosity = 1.; //pb^-1
    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    if(dName.find("Data")<=0 || dName.find("data")<=0 || dName.find("DATA")<=0)
    {
        Luminosity = theDataset->EquivalentLumi();
        cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    // add jobs number at the end of file if 
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();
    //Output ROOT file
    string outputDirectory("MACRO_Output"+channelpostfix);
    mkdir(outputDirectory.c_str(),0777);
    
    string rootFileName (outputDirectory+"/DisplacedTop"+postfix+channelpostfix+".root");
    if (strJobNum != "0")
      {
	cout << "strJobNum is " << strJobNum << endl;
        rootFileName = outputDirectory+"/DisplacedTop"+postfix+channelpostfix+"_"+strJobNum+".root";
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

    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    MSPlot["cutFlow"]                                  = new MultiSamplePlot(datasets, "cutFlow", 15, -0.5, 14.5, "cutFlow");
    //Muons
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
    MSPlot["LeptonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{lep}");
    MSPlot["MuonRelIsolation"]                              = new MultiSamplePlot(datasets, "MuonRelIsolation", 10, 0, .25, "RelIso");
    //Electrons
    MSPlot["ElectronRelIsolation"]                          = new MultiSamplePlot(datasets, "ElectronRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["ElectronPt"]                                    = new MultiSamplePlot(datasets, "ElectronPt", 30, 0, 300, "PT_{e}");
    MSPlot["NbOfElectronsPreJetSel"]                           = new MultiSamplePlot(datasets, "NbOfElectronsPreJetSel", 10, 0, 10, "Nb. of electrons");
    //Init Electron Plots
    MSPlot["InitElectronPt"]                                = new MultiSamplePlot(datasets, "InitElectronPt", 30, 0, 300, "PT_{e}");
    MSPlot["InitElectronEta"]                               = new MultiSamplePlot(datasets, "InitElectronEta", 40, -4, 4, "#eta");
    MSPlot["NbOfElectronsInit"]                             = new MultiSamplePlot(datasets, "NbOfElectronsInit", 10, 0, 10, "Nb. of electrons");
    MSPlot["InitElectronRelIsolation"]                      = new MultiSamplePlot(datasets, "InitElectronRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["InitElectronSuperClusterEta"]                   = new MultiSamplePlot(datasets, "InitElectronSuperClusterEta", 10, 0, 2.5, "#eta");
    MSPlot["InitElectrondEtaI"]                             = new MultiSamplePlot(datasets, "InitElectrondEtaI", 20, 0, .05, "#eta");
    MSPlot["InitElectrondPhiI"]                             = new MultiSamplePlot(datasets, "InitElectrondPhiI", 20, 0, .2, "#phi");
    MSPlot["InitElectronHoverE"]                            = new MultiSamplePlot(datasets, "InitElectronHoverE", 10, 0, .15, "H/E");
    MSPlot["InitElectrond0"]                                = new MultiSamplePlot(datasets, "InitElectrond0", 20, 0, .1, "d0");
    MSPlot["InitElectrondZ"]                                = new MultiSamplePlot(datasets, "InitElectrondZ", 10, 0, .25, "dZ");
    MSPlot["InitElectronEminusP"]                           = new MultiSamplePlot(datasets, "InitElectronEminusP", 10, 0, .25, "1/GeV");
    MSPlot["InitElectronConversion"]                        = new MultiSamplePlot(datasets, "InitElectronConversion", 2, 0, 2, "Conversion Pass");
    MSPlot["InitElectronMissingHits"]                       = new MultiSamplePlot(datasets, "InitElectronMissingHits", 10, 0, 10, "MissingHits");
    MSPlot["InitElectronCutFlow"]                           = new MultiSamplePlot(datasets, "InitElectronCutFlow", 12, 0, 12, "CutNumber");
    //B-tagging discriminators
    MSPlot["Bdisc_CSV_jet1"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet1", 30, 0, 1, "CSV b-disc._{jet1}");
    MSPlot["Bdisc_CSV_jet2"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet2", 30, 0, 1, "CSV b-disc._{jet2}");
    MSPlot["Bdisc_CSV_jet3"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet3", 30, 0, 1, "CSV b-disc._{jet3}");
    MSPlot["Bdisc_CSV_jet4"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet4", 30, 0, 1, "CSV b-disc._{jet4}");
    MSPlot["Bdisc_CSV_jet5"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet5", 30, 0, 1, "CSV b-disc._{jet5}");
    MSPlot["Bdisc_CSV_jet6"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet6", 30, 0, 1, "CSV b-disc._{jet6}");
    MSPlot["Bdisc_CSV_Bjet1"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet1", 30, 0, 1, "CSV b-disc._{bjet1}");
    MSPlot["Bdisc_CSV_Bjet2"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet2", 30, 0, 1, "CSV b-disc._{bjet2}");
    MSPlot["Bdisc_CSV_Bjet3"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet3", 30, 0, 1, "CSV b-disc._{bjet3}");
    MSPlot["Bdisc_CSV_Bjet4"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet4", 30, 0, 1, "CSV b-disc._{bjet4}");
    MSPlot["Bdisc_CSV_Bjet5"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet5", 30, 0, 1, "CSV b-disc._{bjet5}");
    MSPlot["Bdisc_CSV_Bjet6"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_Bjet6", 30, 0, 1, "CSV b-disc._{bjet6}");
    //Jets
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["NbJetsPreJetSel"]                                        = new MultiSamplePlot(datasets, "NbJetsPreJetSel", 15,-0.5, 14.5, "nb. jets");
    MSPlot["NbCSVLJetsPreJetSel_0"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_0", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_0"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_0", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_0"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_0", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_0"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_0", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbCSVLJetsPreJetSel_1"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_1", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_1"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_1", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_1"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_1", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_1"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_1", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbCSVLJetsPreJetSel_2"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_2", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_2"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_2", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_2"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_2", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_2"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_2", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbCSVLJetsPreJetSel_3"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_3", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_3"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_3", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_3"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_3", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_3"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_3", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbCSVLJetsPreJetSel_4"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_4", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_4"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_4", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_4"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_4", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_4"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_4", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbCSVLJetsPreJetSel_5"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel_5", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel_5"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel_5", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel_5"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel_5", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["NbLightJetsPreJetSel_5"]                                        = new MultiSamplePlot(datasets, "NbLightJetsPreJetSel_5", 15,-0.5, 14.5, "nb. Light tags");
    MSPlot["NbJets"]                                        = new MultiSamplePlot(datasets, "NbJets", 15,-0.5, 14.5, "nb. jets");
    MSPlot["NbCSVLJets"]                                        = new MultiSamplePlot(datasets, "NbCSVLJets", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJets"]                                        = new MultiSamplePlot(datasets, "NbCSVMJets", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJets"]                                        = new MultiSamplePlot(datasets, "NbCSVTJets", 15,-0.5, 14.5, "nb. CSVT tags");
    MSPlot["1stJetPt"]                                      = new MultiSamplePlot(datasets, "1stJetPt", 30, 0, 300, "PT_{jet1}");
    MSPlot["2ndJetPt"]                                      = new MultiSamplePlot(datasets, "2ndJetPt", 30, 0, 300, "PT_{jet2}");
    MSPlot["3rdJetPt"]                                      = new MultiSamplePlot(datasets, "3rdJetPt", 30, 0, 300, "PT_{jet3}");
    MSPlot["4thJetPt"]                                      = new MultiSamplePlot(datasets, "4thJetPt", 30, 0, 300, "PT_{jet4}");
    MSPlot["5thJetPt"]                                      = new MultiSamplePlot(datasets, "5thJetPt", 30, 0, 300, "PT_{jet5}");
    MSPlot["6thJetPt"]                                      = new MultiSamplePlot(datasets, "6thJetPt", 30, 0, 300, "PT_{jet6}");
    MSPlot["1stBJetPt"]                                      = new MultiSamplePlot(datasets, "1stBJetPt", 30, 0, 300, "PT_{bjet1}");
    MSPlot["2ndBJetPt"]                                      = new MultiSamplePlot(datasets, "2ndBJetPt", 30, 0, 300, "PT_{bjet2}");
    MSPlot["3rdBJetPt"]                                      = new MultiSamplePlot(datasets, "3rdBJetPt", 30, 0, 300, "PT_{bjet3}");
    MSPlot["4thBJetPt"]                                      = new MultiSamplePlot(datasets, "4thBJetPt", 30, 0, 300, "PT_{bjet4}");
    MSPlot["5thBJetPt"]                                      = new MultiSamplePlot(datasets, "5thBJetPt", 30, 0, 300, "PT_{bjet5}");
    MSPlot["6thBJetPt"]                                      = new MultiSamplePlot(datasets, "6thBJetPt", 30, 0, 300, "PT_{bjet6}");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    //MET
    MSPlot["MET_preCut"]                                           = new MultiSamplePlot(datasets, "MET_preCut", 70, 0, 700, "MET");
    MSPlot["MT_LepMET_preCut"]                                           = new MultiSamplePlot(datasets, "MET_LepMET_preCut", 70, 0, 700, "MT(lep,MET)");
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");
    MSPlot["MT_LepMET"]                                           = new MultiSamplePlot(datasets, "MT_LepMET", 70, 0, 700, "MT(lep,MET)");
	//MC-info plots
    MSPlot["JetID_Hmother"]                                           = new MultiSamplePlot(datasets, "JetID_Hmother", 12,-0.5,11.5, "jet number");

    ///////////////////
    // 1D histograms //
    ///////////////////

    ///////////////////
    // 2D histograms //
    ///////////////////
    histo2D["JetID_vs_pdgID"] = new TH2F("JetID_vs_pdgID","parton pdgID:jet number",12,-0.5,11.5, 61, -30.5,30.5);


    LumiReWeighting LumiWeights;
    LumiWeights = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root", "pileup", "pileup");    

    /////////////////////////////////
    //       Loop on datasets      //
    /////////////////////////////////
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    TFile * btaghistos = new TFile("histos/HistosPtEta.root");

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename2 = "";
        string currentfilename2 = "";

        bool nlo = false;

        if(dName.find("NLO") != std::string::npos || dName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;

        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;
        float normfactor = datasets[d]->NormFactor();

        ////////////////////////////////////////////////////////////
        // Setup Date string and nTuple for output  
        ///////////////////////////////////////////////////////////

        SourceDate *strdate = new SourceDate();
        string date_str = strdate->ReturnDateStr();
        if(debug)cout<<"date print"<<endl;

        string channel_dir = "Trees_SelectionOutput"+channelpostfix;
        string date_dir = channel_dir+"/Trees_SelectionOutput_" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);

        string Ntupname = "Trees_SelectionOutput"+channelpostfix+"/Trees_SelectionOutput_"+ date_str  +"/FCNC_1L3B_" +postfix + channelpostfix + ".root";
        string Ntuptitle = "tree";

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
        TNtuple * tup      = new TNtuple(Ntuptitle.c_str(), Ntuptitle.c_str(), "leptonpt:bdisc1:bdisc2:bdisc3:nb_jets:nb_bjets:jet1_Pt:jet2_Pt:jet3_Pt:MissingEt:leptoneta:MTlepmet:MLepTop:MHadTop:EtaLepTop:EtaHadTop:MassW:EtaW");
        
        
        if(debug)cout<<"created craneens"<<endl;
        //////////////////////////////////////////////////
        // Pre-event loop definitions
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1, start = 0;
        int currentRun;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;

        int event_start = startEvent;

        if (debug) cout << " - Loop over events " << endl;

        float leptonpt,bdisc1,bdisc2,bdisc3,nb_jets,nb_bjets,jet1_Pt,jet2_Pt,jet3_Pt,MissingEt,leptoneta,MTlepmet,MLepTop,MHadTop,EtaLepTop,EtaHadTop,MassW,EtaW;

        double end_d = ending;
        if(endEvent > ending) end_d = ending;
        else end_d = endEvent;

        cout <<"Will run over "<<  end_d<< " events..."<<endl;    cout <<"Starting event = = = = "<< event_start  << endl;


        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
	        leptonpt = -1., bdisc1 = -1., bdisc2 = -1., bdisc3 = -1., nb_jets = -1, nb_bjets = -1., jet1_Pt = -1, jet2_Pt = -1.,jet3_Pt = -1., MissingEt = -1.;
            leptoneta = -99999.; MTlepmet = -1.; MLepTop = -1.; MHadTop = -1.; EtaLepTop = -99999.; EtaHadTop = -99999.; MassW = -1.; EtaW = -99999.;
           

            double ievt_d = ievt;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            if(debug)cout<<"before tree load"<<endl;
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            if(debug)cout<<"after tree load"<<endl;

            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;

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
		    vector<TLorentzVector> selectedMuonsTLV, selectedElectronsTLV, metsTLV, selectedJetsTLV, selectedLeptonsTLV;
            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC, selectedLooseIsoMuonsTLV;
            vector<TLorentzVector> mcParticlesTLV, mcMuonsTLV, mcPartonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_,mcParticles;
            vector<int> mcMuonIndex, mcPartonIndex;
            JetPartonMatching muonMatching, jetMatching;


            currentRun = event->runId();

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
				selectedJets                                        = r2selection.GetSelectedJets(30,2.5,true,"Tight"); // ApplyJetId
				if (debug)cout<<"Getting Tight Muons"<<endl;
				selectedMuons                                       = r2selection.GetSelectedMuons(30,2.1,0.15, "Tight", "Spring15"); //Selected
				if (debug)cout<<"Getting Loose Electrons"<<endl;
				selectedElectrons                                   = r2selection.GetSelectedElectrons(20,2.5,"Loose", "Spring15_25ns", true); //Vetoed  
				if (debug)cout<<"Getting Loose Muons"<<endl;
				selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20,"Loose","Spring15"); //Vetoed         
            }
            if (Electron)
            {
				if (debug)cout<<"Getting Jets"<<endl;
				selectedJets                                        = r2selection.GetSelectedJets(30,2.5,true,"Tight"); // ApplyJetId
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
            bool trigged = true;
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();
            MSPlot["cutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor );
            weightCount += scaleFactor;
            eventCount++;




            ////////////////////////////////////////////
            // Applying baseline selection //
            ///////////////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection

            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            MSPlot["cutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor );


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
            MSPlot["cutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor );

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
			MSPlot["cutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor );

			
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
			
			
			//Calculations for MT(lep,MET) selection cut
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

			MSPlot["MT_LepMET_preCut"] ->Fill(MT, datasets[d], true, Luminosity*scaleFactor );
			MSPlot["MET_preCut"] ->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor );

			if(MT <= 50) continue;
			MSPlot["cutFlow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor );
	        if (debug)	cout <<"Cut on MT..."<<endl;

            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output
			
		    //Fill b-jet collections
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
		  	

			//////////////////////////////////////
			// Pre-jet histograms //
			/////////////////////////////////////
			MSPlot["NbJetsPreJetSel"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVLJetsPreJetSel_0"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVMJetsPreJetSel_0"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVTJetsPreJetSel_0"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbLightJetsPreJetSel_0"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			if(selectedJets.size() >= 1)
			{
				MSPlot["NbCSVLJetsPreJetSel_1"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVMJetsPreJetSel_1"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVTJetsPreJetSel_1"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbLightJetsPreJetSel_1"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 2)
			{
				MSPlot["NbCSVLJetsPreJetSel_2"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVMJetsPreJetSel_2"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVTJetsPreJetSel_2"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbLightJetsPreJetSel_2"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 3)
			{
				MSPlot["NbCSVLJetsPreJetSel_3"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVMJetsPreJetSel_3"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVTJetsPreJetSel_3"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbLightJetsPreJetSel_3"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 4)
			{
				MSPlot["NbCSVLJetsPreJetSel_4"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVMJetsPreJetSel_4"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVTJetsPreJetSel_4"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbLightJetsPreJetSel_4"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 5)
			{
				MSPlot["NbCSVLJetsPreJetSel_5"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVMJetsPreJetSel_5"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbCSVTJetsPreJetSel_5"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbLightJetsPreJetSel_5"] ->Fill(selectedLightJets_MWP.size(), datasets[d], true, Luminosity*scaleFactor);
			}


			if(selectedJets.size() < 3)  continue;
			MSPlot["cutFlow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor );
	        if (debug)	cout <<"Cut on nb jets..."<<endl;


		  	if(selectedMBJets.size() < 3) continue;
	        if (debug)	cout <<"Cut on nb b-jets..."<<endl;
			MSPlot["cutFlow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor );

			MSPlot["cutFlow"]->Fill(7, datasets[d], true, Luminosity*scaleFactor );

            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;


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
                    
                // take all the selectedJets_ to study the radiation stuff, selectedJets_ are already ordened in decreasing Pt()
                for (unsigned int i = 0; i < selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
                    
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

            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);




            /////////////////////////////////////
            // Muon Based Plots //
            /////////////////////////////////////
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["LeptonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                float reliso = selectedMuons[selmu]->relPfIso(4, 0.5);
                MSPlot["MuonRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////////////////////////
            // Electron Based Plots //
            /////////////////////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronPt"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["LeptonPt"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////////////////
            // Jets Based Plots //
            //////////////////////////////////
			MSPlot["NbJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVLJets"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVMJets"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbCSVTJets"]->Fill(selectedTBJets.size(), datasets[d], true, Luminosity*scaleFactor);

			if(selectedJets.size() >= 1)
			{
				MSPlot["Bdisc_CSV_jet1"]->Fill(selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["1stJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 2)
			{
				MSPlot["Bdisc_CSV_jet2"]->Fill(selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["2ndJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 3)
			{
				MSPlot["Bdisc_CSV_jet3"]->Fill(selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["3rdJetPt"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 4)
			{
				MSPlot["Bdisc_CSV_jet4"]->Fill(selectedJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 5)
			{
				MSPlot["Bdisc_CSV_jet5"]->Fill(selectedJets[4]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedJets.size() >= 6)
			{
				MSPlot["Bdisc_CSV_jet6"]->Fill(selectedJets[5]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			//B-jets
			if(debug) cout << "selectedMBJets.size() = "<< selectedMBJets.size() << endl;
			if(selectedMBJets.size() >= 1)
			{
				MSPlot["Bdisc_CSV_Bjet1"]->Fill(selectedMBJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["1stBJetPt"]->Fill(selectedMBJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedMBJets.size() >= 2)
			{
				MSPlot["Bdisc_CSV_Bjet2"]->Fill(selectedMBJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["2ndBJetPt"]->Fill(selectedMBJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedMBJets.size() >= 3)
			{
				MSPlot["Bdisc_CSV_Bjet3"]->Fill(selectedMBJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["3rdBJetPt"]->Fill(selectedMBJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedMBJets.size() >= 4)
			{
				MSPlot["Bdisc_CSV_Bjet4"]->Fill(selectedMBJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["4thBJetPt"]->Fill(selectedMBJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedMBJets.size() >= 5)
			{
				MSPlot["Bdisc_CSV_Bjet5"]->Fill(selectedMBJets[4]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["5thBJetPt"]->Fill(selectedMBJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}
			if(selectedMBJets.size() >= 6)
			{
				MSPlot["Bdisc_CSV_Bjet6"]->Fill(selectedMBJets[5]->btag_combinedInclusiveSecondaryVertexV2BJetTags(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["6thBJetPt"]->Fill(selectedMBJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
			}

			for(unsigned int i = 0; i < selectedJets.size(); i++)
			{
				MSPlot["JetEta"]->Fill(selectedJets[i]->Eta(), datasets[d], true, Luminosity*scaleFactor);
			}

            float HT = 0, H = 0;


            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                //Event-level variables
                HT = HT + selectedJets[seljet1]->Pt();
                H = H +  selectedJets[seljet1]->P();
            }

            MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);


            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();

            /////////////////////////////////
            //MET Based Plots//
            /////////////////////////////////

            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["MT_LepMET"]->Fill(MT, datasets[d], true, Luminosity*scaleFactor);

            /////////////////////////////////////////////////////////
            //Topology Reconstructions (truth level)
            ////////////////////////////////////////////////////////
            if( int(leptonicBJet_.first) != 9999) //Leptonic decays ///////////////////// CHANGE TO LEPTONIC TOP
            {
                int ind_jet = leptonicBJet_.first;///////////////////// CHANGE TO LEPTONIC TOP
                MLepTop = (selectedJetsTLV[ind_jet] + selectedLeptonsTLV[0] + metsTLV[0]).M();
                EtaLepTop = (selectedJetsTLV[ind_jet] + selectedLeptonsTLV[0] + metsTLV[0]).Eta();
            }
            if(int(hadronicWJet1_.first) != 9999 && int(hadronicWJet2_.first) != 9999) //Hadronic decays
            {
                int ind_jet1 = hadronicWJet1_.first;
                int ind_jet2 = hadronicWJet2_.first;
                MassW = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2]).M();
                EtaW = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2]).Eta();

                if(int(hadronicBJet_.first) != 9999)/////////////////////// CHANGE TO HADRONIC TOP
                {
                    int ind_jet_top = hadronicBJet_.first;///////////////////// CHANGE TO LEPTONIC TOP
                    MHadTop = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2] + selectedJetsTLV[ind_jet_top]).M();
                    EtaHadTop = (selectedJetsTLV[ind_jet1] + selectedJetsTLV[ind_jet2] + selectedJetsTLV[ind_jet_top]).Eta();
                }
                
            }
            
            ///////////////////////////
            //Filling nTuple//
            ///////////////////////////
			if(Muon && !Electron)
			{
			    leptonpt = selectedMuons[0]->Pt();
			    leptoneta = selectedMuons[0]->Eta();
			}
			else if(!Muon && Electron)
			{
			    leptonpt = selectedElectrons[0]->Pt();
			    leptoneta = selectedElectrons[0]->Eta();
			}
			bdisc1 = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			bdisc2 = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			bdisc3 = selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
			nb_jets = selectedJets.size();
			nb_bjets = selectedMBJets.size();
			jet1_Pt = selectedJets[0]->Pt();
			jet2_Pt = selectedJets[1]->Pt();
			jet3_Pt = selectedJets[2]->Pt();
			MissingEt = mets[0]->Et();
			MTlepmet = MT;

            float vals[18] = {leptonpt,bdisc1,bdisc2,bdisc3,nb_jets,nb_bjets,jet1_Pt,jet2_Pt,jet3_Pt,MissingEt,leptoneta,MTlepmet,MLepTop,MHadTop,EtaLepTop,EtaHadTop,MassW,EtaW};
            tup->Fill(vals);


			///////////////////////////
			// Event info //////
			///////////////////////////
			if(Muon && !Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedMuons[0]->type(), selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), selectedMuons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());
			else if(!Muon && Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->type(), selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), selectedElectrons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());


        } //End Loop on Events

		tup->Write();
       	tupfile->Close();
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    fclose (eventlist);

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;

    fout->cd();

	//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Write(fout, name, false, outputDirectory, "png");
    }

    TDirectory* th2dir = fout->mkdir("Histos2D");
    th2dir->cd();


    for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {

        TH2F *temp = it->second;
        temp->Write();
    }
    delete fout;



/*    if(Muon){

        cout<<"TRIGGGG"<<endl;

        cout<<"preTrig: "<<preTrig<<"   postTrig: "<<postTrig<<endl;
        cout<<"********"<<endl;
    }
*/



    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}





/*

            float HTb = 0.;
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244   )
                {
                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.679)
                    {
                        HTb += selectedJets[seljet]->Pt();
                        selectedMBJets.push_back(selectedJets[seljet]);
                        if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.898)
                        {
                            selectedTBJets.push_back(selectedJets[seljet]);
                        }
                    }
                }
                else
                {
                    selectedLightJets.push_back(selectedJets[seljet]);
                }
            }
            float nJets = selectedJets.size(); //Number of Jets in Event
            float nLtags = selectedLBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM & CSVT)
            float nMtags = selectedMBJets.size(); //Number of CSVM tags in Event (includes jets that pass CSVT)
            float nTtags = selectedTBJets.size(); //Number of CSVT tags in Event 
            float nLights=selectedLightJets.size();

            selectedLBJets.clear();
            selectedMBJets.clear();
            selectedTBJets.clear();
            selectedLightJets.clear();

            float temp_HT = 0.;
            float HT_leading = 0.;
            float HT_lagging = 0.;
            float HTRat = 0;
            for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
                temp_HT += selectedJets[seljet0]->Pt();
                if (seljet0 < 4){
                    HT_leading += selectedJets[seljet0]->Pt();
                }else{
                    HT_lagging += selectedJets[seljet0]->Pt();
                }
            }

            //HTRat = HT_leading/HT_lagging;



            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////


            float weight_0 = 1; //nominal
            float weight_1 = 1;
            float weight_2 = 1;
            float weight_3 = 1;
            float weight_4 = 1;
            float weight_5 = 1;
            float weight_6 = 1;
            float weight_7 = 1;
            float weight_8 = 1;



            // Apply trigger selection
            //trigged = treeLoader.EventTrigged (itrigger);
            bool trigged = false;  // Disabling the HLT requirement

            trigged=true;
            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();

            preTrig++;
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(dName.find("Data") != string::npos || dName.find("data") != string::npos || dName.find("DATA") != string::npos){
            if (!trigged)          continue;  //If an HLT condition is not present, skip this event in the loop.       
            }
            postTrig++;

            if(dName.find("TTJets") != string::npos){
                weight_0 = (runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1")/(abs(event->originalXWGTUP()))); //nominal                   
            }
            weightCount += scaleFactor;
            eventCount++;


            if(nlo)
            {
                if(weight_0 < 0.0)
                {
                    scaleFactor = -1.0;  //Taking into account negative weights in NLO Monte Carlo
                    negWeights++;
                }
            }


            //Filling Histogram of the number of vertices before Event Selection
            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

            if (!isGoodPV) continue; // Check that there is a good Primary Vertex

            if (debug) cout <<"Number of Muons = "<< nMu <<"    Electrons =  "  <<nEl<<"     Jets = "<< selectedJets.size()   <<" loose BJets = "<<  nLtags   <<
                "  MuonChannel = "<<Muon<<" Electron Channel"<<Electron<<endl;

            //Apply the lepton, jet, btag and HT & MET selections
            if (Muon)
            {
                if  (  !( nMu == 1 && nEl == 0 && nLooseMu == 1 && nJets>=6 && nMtags >=2)) continue; // Muon Channel Selection

            }
            else if(Electron){
                if  (  !( nMu == 0 && nEl == 1 && nLooseEl == 1 && nJets>=6 && nMtags >=2)) continue; // Electron Channel Selection
            }
            else{
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            if(dName.find("TTJets")!=string::npos){
                weight_1 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 2")/(abs(event->originalXWGTUP())));                
                weight_2 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 3")/(abs(event->originalXWGTUP())));                
                weight_3 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 4")/(abs(event->originalXWGTUP())));                
                weight_4 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 5")/(abs(event->originalXWGTUP())));                
                weight_5 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 6")/(abs(event->originalXWGTUP())));                
                weight_6 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 7")/(abs(event->originalXWGTUP())));                
                weight_7 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 8")/(abs(event->originalXWGTUP())));                
                weight_8 = (runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 9")/(abs(event->originalXWGTUP())));                  
            }
            //cout<<"weight0: "<<weight_0<<" weight1: "<<weight_1<<" weight2: "<<weight_2<<endl;

            if(dName.find("Data") !=string::npos || dName.find("data") != string::npos || dName.find("DATA") != string::npos)
            {
                lumiWeight=1;
            }
            else{
                lumiWeight = LumiWeights.ITweight( vertex.size() );
            }
            scaleFactor = scaleFactor * lumiWeight;
            float bTagEff(-1);
            // if(bTagReweight && dName.find("Data")==string::npos){
            // //get btag weight info
            //     for(int jetbtag = 0; jetbtag<selectedJets.size(); jetbtag++){
            //         float jetpt = selectedJets[jetbtag]->Pt();
            //         float jeteta = selectedJets[jetbtag]->Eta();
            //         float jetdisc = selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
            //         BTagEntry::JetFlavor jflav;
            //         int jetpartonflav = std::abs(selectedJets[jetbtag]->partonFlavour());
            //         if(debug) cout<<"parton flavour: "<<jetpartonflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
            //         if(jetpartonflav == 5){
            //             jflav = BTagEntry::FLAV_B;
            //         }
            //         else if(jetpartonflav == 4){
            //             jflav = BTagEntry::FLAV_C;
            //         }
            //         else{
            //             jflav = BTagEntry::FLAV_UDSG;
            //         }
            //         bTagEff = bTagReader->eval(jflav, jeteta, jetpt, jetdisc);                     
            //         if(debug)cout<<"btag efficiency = "<<bTagEff<<endl;       
            //     }      
            //     btwt->FillMCEfficiencyHistos(selectedJets); 
            // }

            float btagWeight = 1;
            if(bTagReweight && dName.find("Data")==string::npos){
                btagWeight =  btwt->getMCEventWeight(selectedJets, false);
                // cout<<"btag weight "<<btagWeight<<endl;
            }

            float fleptonSF = 0;
            if(bLeptonSF){
                if(Muon){
                    fleptonSF = muonSFWeight->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                }
                else if(Electron){
                    fleptonSF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                }
            }
            if(debug) cout<<"lepton SF:  "<<fleptonSF<<endl;

            scaleFactor *= fleptonSF;
            scaleFactor *= btagWeight;

            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            TRootGenEvent* genEvt = 0;

            sort(selectedJets.begin(),selectedJets.end(),HighestCSVBtag());


            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////

            //////////////////////
            // Muon Based Plots //
            //////////////////////
            float leptonIso=0;
            float chargedHIso = 0;
            float neutralHIso = 0;
            float photonIso = 0;
            float PUIso = 0;
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                float relisomu = (selectedMuons[0]->chargedHadronIso(4) + max( 0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.) ) / selectedMuons[0]->Pt();
                chargedHIso = selectedMuons[0]->chargedHadronIso(4);
                neutralHIso = selectedMuons[0]->neutralHadronIso(4);
                photonIso = selectedMuons[0]->photonIso(4);
                PUIso = selectedMuons[0]->puChargedHadronIso(4);

                
                MSPlot["leptonIso"]->Fill(relisomu, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                leptonIso = relisomu;
            }

            //////////////////////////
            // Electron Based Plots //
            //////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                leptonIso = reliso;
            }

            //////////////////////
            // Jets Based Plots //
            //////////////////////

            HT = 0;
            float HT1M2L=0, H1M2L=0, HTbjets=0, HT2M=0, H2M=0, HT2L2J=0;
            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                //Event-level variables
                jetpt = selectedJets[seljet1]->Pt();
                HT = HT + jetpt;
                H = H +  selectedJets[seljet1]->P();
                if (seljet1 > 4  )  HTHi +=  selectedJets[seljet1]->Pt();
            }

            HTH = HT/H;
            HTRat = HTHi/HT;

            float selectedLeptonPt = 0 ;
            if(Muon){
                selectedLeptonPt = selectedMuons[0]->Pt();
            }
            else if(Electron){
                selectedLeptonPt = selectedElectrons[0]->Pt();
            }

       

            MSPlot["HTb_SelectedJets"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["NbOfSelectedBJets"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
            float met = mets[0]->Et();
            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);

            if(nJets>5){
                MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }

            MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);


            if(Muon){
                muonpt  = selectedMuons[0]->Pt();
                muoneta = selectedMuons[0]->Eta();
                leptonphi = selectedMuons[0]->Phi();
            }
            else if(Electron){
                muonpt  = selectedElectrons[0]->Pt();
                muoneta = selectedElectrons[0]->Eta();
                leptonphi = selectedElectrons[0]->Phi();

            }

            bjetpt = selectedMBJets[0]->Pt();


            float firstjetpt = selectedJets[0]->Pt();
            float secondjetpt = selectedJets[1]->Pt();
            
            float nvertices = vertex.size();
            float angletoplep = 0;
            float angletop1top2 = 0;
            float diTopness = 0;

            float vals[39] = {BDTScore,nJets,nLtags,nMtags,nTtags,HT,muonpt,muoneta,bjetpt,HT2M,HTb,HTH,HTRat,diTopness,scaleFactor,nvertices,normfactor,Luminosity,weight_0,weight_1,weight_2,weight_3,weight_4,weight_5,weight_6,weight_7,weight_8,met,angletop1top2,angletoplep, firstjetpt,secondjetpt, leptonIso, leptonphi, chargedHIso,neutralHIso,photonIso,PUIso,lumiWeight};
            tupfile->cd();
            tup->Fill(vals);

			///////////////////////////
			// Event info //////
			///////////////////////////
			if(Muon && !Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedMuons[0]->type(), selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), selectedMuons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());
			else if(!Muon && Electron)
			fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->type(), selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), selectedElectrons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());


        } //End Loop on Events

        tupfile->cd();
        tup->Write();
        tupfile->Close();
        
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    fclose (eventlist);
*/
