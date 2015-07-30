//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:
//all-had ->679 * .46 = 312.34
//semi-lep ->679 *.45 = 305.55
//di-lep-> 679* .09 = 61.11

//ttbar @ NNLO 8 TeV:
//all-had -> 245.8 * .46 = 113.068
//semi-lep-> 245.8 * .45 = 110.61
//di-lep ->  245.8 * .09 = 22.122

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
//#include "TopTreeAnalysisBase/Selection/interface/FCNC_1L3BSelectionTable.h"

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

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = false;
bool split_leptonflavours = false;
bool debug = false;
float topness;
string btagger = "CSVM";

int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TProfile*> histoProfile;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

struct HighestCSVBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
int Factorial(int N);
float Sphericity(vector<TLorentzVector> parts );
float Centrality(vector<TLorentzVector> parts);

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 14)
    {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
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
    const float PreJetSelEff           = strtod(argv[10], NULL);
    string fileName                 = argv[11];
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
    cout << "----------------------------------------" << endl;
//    cin.get();



    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");

    int passed = 0;
    int ndefs =0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;

    float scalefactorbtageff, mistagfactor;
	float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
	float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
	float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.

    clock_t start = clock();

    BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_NOttbar_payload_EPS13.txt", btagger) ;

    int doJESShift = 0; // 0: off 1: minus 2: plus
    cout << "doJESShift: " << doJESShift << endl;

    int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
    cout << "doJERShift: " << doJERShift << endl;

    int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
    cout << "dobTagEffShift: " << dobTagEffShift << endl;

    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FCNC_1L3B search ! "           << endl;
    cout << "*************************************************************" << endl;


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
    bool Muon = true;
    bool Electron = false;

    if(Muon && !Electron)
    {
        cout << " --> Using the Muon channel..." << endl;
        channelpostfix = "_Mu";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
    }
    else if(!Muon && Electron)
    {
        cout << " --> Using the Electron channel..." << endl;
        channelpostfix = "_El";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
    }
    else
    {
        cerr<<"Correct lepton Channel not selected."<<endl;
        exit(1);
    }

	if(!split_leptonflavours) channelpostfix = "";

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

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
    int verbose = 2;//anaEnv.Verbose;



    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi*normf);
    datasets.push_back(theDataset);
    float Luminosity = 15000.0; //pb^-1??


    string dataSetName;

    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    dataSetName = theDataset->Name();
    if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
    {
        Luminosity = theDataset->EquivalentLumi();
        cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    //Output ROOT file
    string outputDirectory("MACRO_Output"+channelpostfix);
    mkdir(outputDirectory.c_str(),0777);
    string rootFileName (outputDirectory+"/FCNC_1L3B"+postfix+channelpostfix+".root");
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

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////

    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    //Muons
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
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
    MSPlot["Bdisc_CSV_jet1"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet1", 20, 0, 1, "CSV b-disc._{jet1}");
    MSPlot["Bdisc_CSV_jet2"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet2", 20, 0, 1, "CSV b-disc._{jet2}");
    MSPlot["Bdisc_CSV_jet3"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet3", 20, 0, 1, "CSV b-disc._{jet3}");
    MSPlot["Bdisc_CSV_jet4"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet4", 20, 0, 1, "CSV b-disc._{jet4}");
    MSPlot["Bdisc_CSV_jet5"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet5", 20, 0, 1, "CSV b-disc._{jet5}");
    MSPlot["Bdisc_CSV_jet6"]                             = new MultiSamplePlot(datasets, "Bdisc_CSV_jet6", 20, 0, 1, "CSV b-disc._{jet6}");
    //Jets
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["NbJetsPreJetSel"]                                        = new MultiSamplePlot(datasets, "NbJetsPreJetSel", 15,-0.5, 14.5, "nb. jets");
    MSPlot["NbCSVLJetsPreJetSel"]                                        = new MultiSamplePlot(datasets, "NbCSVLJetsPreJetSel", 15,-0.5, 14.5, "nb. CSVL tags");
    MSPlot["NbCSVMJetsPreJetSel"]                                        = new MultiSamplePlot(datasets, "NbCSVMJetsPreJetSel", 15,-0.5, 14.5, "nb. CSVM tags");
    MSPlot["NbCSVTJetsPreJetSel"]                                        = new MultiSamplePlot(datasets, "NbCSVTJetsPreJetSel", 15,-0.5, 14.5, "nb. CSVT tags");
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
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    //MET
    MSPlot["MET_preCut"]                                           = new MultiSamplePlot(datasets, "MET_preCut", 70, 0, 700, "MET");
    MSPlot["MT_LepMET_preCut"]                                           = new MultiSamplePlot(datasets, "MET_LepMET_preCut", 70, 0, 700, "MT(lep,MET)");
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");
    MSPlot["MT_LepMET"]                                           = new MultiSamplePlot(datasets, "MT_LepMET", 70, 0, 700, "MT(lep,MET)");

    ///////////////////
    // 1D histograms //
    ///////////////////

    ///////////////////
    // 2D histograms //
    ///////////////////

    //Plots
    string pathPNG = "MSPlots_FCNC_1L3B"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);

    cout <<"Making directory :"<< pathPNG  <<endl;
    vector<string> CutsselecTable;
    /////////////////////////////////
    // Selection table: Lepton + jets
    /////////////////////////////////
    if(Muon && !Electron)
    {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("Event cleaning and Trigger"));
        CutsselecTable.push_back(string("Exactly 1 Tight Isolated Muon")); //Including pt-cut of 30 GeV
        CutsselecTable.push_back(string("Loose leptons veto"));
        CutsselecTable.push_back(string("Mt(lep,MET) > 50 GeV"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("At least 3 CSVM Jets"));
        CutsselecTable.push_back(string("Pt cuts on jets (BL)")); //End of baseline selection
    }
    if(!Muon && Electron)
    {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("Event cleaning and Trigger"));
        CutsselecTable.push_back(string("Exactly 1 Tight Isolated Electron")); //Including pt-cut of 30 GeV
        CutsselecTable.push_back(string("Loose leptons veto"));
        CutsselecTable.push_back(string("Mt(lep,MET) > 50 GeV"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("At least 3 CSVM Jets"));
        CutsselecTable.push_back(string("Pt cuts on jets (BL)")); //End of baseline selection
    }


    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;
        bool nlo = false;
//    	bool bx25 = false;
        dataSetName = datasets[d]->Name();
//        if(dataSetName.find("bx50") != std::string::npos) bx25 = false;
//        else bx25 = true;

        if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;

//        if(bx25) cout << "Dataset with 25ns Bunch Spacing!" <<endl;
//        else cout << "Dataset with 50ns Bunch Spacing!" <<endl;
        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;


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

        //string dataSetName = datasets[d]->Name();
        string channel_dir = "SelectionOutput"+channelpostfix;
        string date_dir = channel_dir+"/SelectionOutput_" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);

        string Ntupname = "SelectionOutput"+channelpostfix+"/SelectionOutput_"+ date_str  +"/FCNC_1L3B_" + dataSetName +postfix + ".root";
        string Ntuptitle = "FCNC_1L3B_" + channelpostfix;

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");

//        TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"Faco:BDT:nJets:nFatJets:nWTags:nTopTags:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2L:HTb:HTH:HTRat:topness:EventSph:EventCen:DiLepSph:DiLepCen:TopDiLepSph:TopDiLepCen:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");


        //////////////////////////////////////////////////
        /// Initialize JEC factors ///////////
        //////////////////////////////////////////////////

        vector<JetCorrectorParameters> vCorrParam;

        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
        }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
//    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalMC")));
//    JetCorrectionUncertainty *jecUncTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));

        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);

        //////////////////////////////////////////////////
        // Loop on events ///////////////////
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1;

        int start = 0;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events in total dataset = "<<  ending  <<endl;

        int event_start = startEvent;
        if (verbose > 1) cout << " - Loop over events " << endl;

//        float faco, BDTScore, MHT, MHTSig, STJet,muoneta, muonpt,electronpt,bjetpt, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX ,HTHi,HTRat, HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;

        double currentfrac =0.;
        double end_d;
        if(endEvent > ending)
            end_d = ending;
        else
            end_d = endEvent;

        //end_d = 10000; //artifical ending for debug
        int nEvents = end_d - event_start;
        cout <<"Will run over "<<  (end_d - event_start) << " events..."<<endl;
        cout <<"Starting event = = = = "<< event_start  << endl;

        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootSubstructureJet*>    selectedFatJets;
        vector<TRootPFJet*>    MVASelJets1;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);

        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////

	for (unsigned int ievt = event_start; ievt < end_d; ievt++)
	  //        for (unsigned int ievt = event_start; ievt < 1000; ievt++)
        {
//	  faco=99.0, BDTScore= -99999.0, MHT = 0.,MHTSig = 0.,muoneta = 0., muonpt =0., electronpt=0., bjetpt =0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0., HTHi =0., HTRat = 0;
//            H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, init_fatjets,  mets, debug);  //load event
            if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
            MSPlot["NbOfElectronsInit"]->Fill(init_electrons.size(), datasets[d], true, Luminosity*scaleFactor );
            for (Int_t initel =0; initel < init_electrons.size(); initel++ )
            {
                float initreliso = init_electrons[initel]->relPfIso(3, 0.5);
                MSPlot["InitElectronPt"]->Fill(init_electrons[initel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronEta"]->Fill(init_electrons[initel]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronRelIsolation"]->Fill(initreliso, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronSuperClusterEta"]->Fill(fabs(init_electrons[initel]->superClusterEta()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectrondEtaI"]->Fill(fabs(init_electrons[initel]->deltaEtaIn()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectrondPhiI"]->Fill(fabs(init_electrons[initel]->deltaPhiIn()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronHoverE"]->Fill(init_electrons[initel]->hadronicOverEm(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectrond0"]->Fill(fabs(init_electrons[initel]->d0()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectrondZ"]->Fill(fabs(init_electrons[initel]->dz()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronEminusP"]->Fill(fabs(1/init_electrons[initel]->E() - 1/init_electrons[initel]->P()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronConversion"]->Fill(init_electrons[initel]->passConversion(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronMissingHits"]->Fill(init_electrons[initel]->missingHits(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["InitElectronCutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor);
                if(init_electrons[initel]->Pt() > 30)
                {
                    MSPlot["InitElectronCutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
                    if(fabs(init_electrons[initel]->Eta()) < 2.5)
                    {
                        MSPlot["InitElectronCutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
                        if(fabs(init_electrons[initel]->deltaEtaIn()) < 0.009277)
                        {
                            MSPlot["InitElectronCutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
                            if(fabs(init_electrons[initel]->deltaPhiIn()) < 0.094739)
                            {
                                MSPlot["InitElectronCutFlow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
                                if(fabs(init_electrons[initel]->hadronicOverEm()) < 0.093068)
                                {
                                    MSPlot["InitElectronCutFlow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
                                    if(fabs(init_electrons[initel]->d0()) < 0.035904)
                                    {
                                        MSPlot["InitElectronCutFlow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
                                        if(fabs(init_electrons[initel]->dz()) < 0.075496)
                                        {
                                            MSPlot["InitElectronCutFlow"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
                                            if(fabs((1/init_electrons[initel]->E()) - (1/init_electrons[initel]->P())) < 0.189968)
                                            {
                                                MSPlot["InitElectronCutFlow"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
                                                if(fabs(init_electrons[initel]->relPfIso(3, 0.5)) < 0.130136)
                                                {
                                                    MSPlot["InitElectronCutFlow"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                                                    if(init_electrons[initel]->passConversion())
                                                    {
                                                        MSPlot["InitElectronCutFlow"]->Fill(10, datasets[d], true, Luminosity*scaleFactor);
                                                        if(fabs(init_electrons[initel]->missingHits()) <= 1)
                                                        {
                                                            MSPlot["InitElectronCutFlow"]->Fill(11, datasets[d], true, Luminosity*scaleFactor);

                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }
            float weight_0 = event->weight0();
            if (debug)cout <<"Weight0: " << weight_0 <<endl;
            if(nlo)
            {
                if(weight_0 < 0.0)
                {
                    scaleFactor = -1.0;  //Taking into account negative weights in NLO Monte Carlo
                    negWeights++;
                }
            }

            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;
            string graphName;

            //////////////////////////////////
            //Loading Gen jets //
            //////////////////////////////////

            vector<TRootGenJet*> genjets;
            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
            {
                // loading GenJets as I need them for JER
                genjets = treeLoader.LoadGenJet(ievt);
            }

            ///////////////////////
            // JER smearing
            //////////////////////

//            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
//            {
//                //JER
//                doJERShift == 0;
//                if(doJERShift == 1)
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
//                else if(doJERShift == 2)
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
//                else
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
//
//                //     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
//
//                // JES sysematic!
//                if (doJESShift == 1)
//                    jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
//                else if (doJESShift == 2)
//                    jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
//
//                //            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
//
//            }

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
//            trigged = treeLoader.EventTrigged (itrigger);
            bool trigged = true;  // Disabling the HLT requirement
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(!trigged)		   continue;  //If an HLT condition is not present, skip this event in the loop.
            // Declare selection instance
            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets);

            // Define object selection cuts
            if (Muon && !Electron)
            {
				if (debug)cout<<"Getting Jets"<<endl;
				selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
				if (debug)cout<<"Getting Tight Muons"<<endl;
				selectedMuons                                       = r2selection.GetSelectedMuons();
				if (debug)cout<<"Getting Tight Electrons"<<endl;
				selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14",true); // VBTF ID    
				if (debug)cout<<"Getting Loose Muons"<<endl;
				selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20);                                   
            }
            if (!Muon && Electron)
            {
				if (debug)cout<<"Getting Jets"<<endl;
				selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
				if (debug)cout<<"Getting Tight Muons"<<endl;
				selectedMuons                                       = r2selection.GetSelectedMuons();
				if (debug)cout<<"Getting Tight Electrons"<<endl;
				selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14", true); // VBTF ID                       
				if (debug)cout<<"Getting Loose Electrons"<<endl;
				selectedExtraElectrons                              = r2selection.GetSelectedElectrons("Loose", "PHYS14", true);
            }


            vector<TRootJet*>      selectedLBJets;
            vector<TRootJet*>      selectedMBJets;
            vector<TRootJet*>      selectedTBJets;
            vector<TRootJet*>      selectedLightJets;

        	int JetCut =0;
            int nMu, nEl, nLooseMu, nLooseEl; //number of (loose) muons/electrons
            if(Muon && !Electron)
            {
              		nMu = selectedMuons.size(); //Number of Muons in Event (Tight only)
              		nEl = selectedElectrons.size(); //Number of Electrons in Event   (tight and loose)
              		nLooseMu = selectedExtraMuons.size();   //Number of loose muons      (loose only)
            }
            else if(Electron && !Muon)
            {
              		nMu = selectedMuons.size(); //Number of Muons in Event (tight and loose)
              		nEl = selectedElectrons.size(); //Number of Electrons in Event (Tight only)
              		nLooseEl = selectedExtraElectrons.size(); //Number of loose electrons  (loose only)
            }
            else cout << "Wrong channel (1)" << endl;

            MSPlot["NbOfElectronsPreJetSel"]->Fill(nEl, datasets[d], true, Luminosity*scaleFactor );


            bool isTagged =false;

			/////////////////////////////////////////////
			// Make TLorentzVectors //
			////////////////////////////////////////////
			vector<TLorentzVector> selectedMuonsTLV, selectedElectronsTLV, metsTLV, selectedJetsTLV;
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
			

//            vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV, mcPartonsTLV;
//            vector<TRootMCParticle*> mcParticlesMatching_;
//            vector<int> mcMuonIndex, mcPartonIndex;
//            JetPartonMatching muonMatching, jetMatching;

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // PreJetSelection looping over Jet Collection                                      /////////////////////////
            // Summing HT and calculating leading, lagging, and ratio for Selected and BJets //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();
            selecTable.Fill(d,0,scaleFactor);
            weightCount += scaleFactor;
            eventCount++;




            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection

//            if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            selecTable.Fill(d,1,scaleFactor);


            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections
            if (Muon && !Electron)
            {
                if  (  !( nMu ==1 && nEl == 0)) continue; // Muon Channel Selection
                if (selectedMuons[0]->Pt() < 30) continue;
            }
            else if (!Muon && Electron)
            {
                if  (  !( nMu == 0 && nEl == 1)) continue; // Electron Channel Selection
                if (selectedElectrons[0]->Pt() < 30) continue;
            }
            else
            {
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            selecTable.Fill(d,2,scaleFactor);

			if(Muon && !Electron)
			{
				if(nLooseMu != 0) continue;
			}
			if(!Muon && Electron)
			{
				if(nLooseEl != 0) continue;
			}
			selecTable.Fill(d,3,scaleFactor);

			//Calculations for MT(lep,MET) selection cut
			float MT = -999;
			if(Muon && !Electron) MT = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
			else if(Electron && !Muon) MT = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
			else cout << "Wrong channel (1)" << endl;

			MSPlot["MT_LepMET_preCut"] ->Fill(MT, datasets[d], true, Luminosity*scaleFactor );
			MSPlot["MET_preCut"] ->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor );

			if(MT <= 50) continue;
			selecTable.Fill(d,4,scaleFactor);

            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output
			if(selectedJets.size() < 3)  continue;
			selecTable.Fill(d,5,scaleFactor);
			
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
		        	}
		   		}
		   		else
		        {
		        	selectedLightJets.push_back(selectedJets[seljet]);
		        }
		  	}
		  	
		  	if(selectedMBJets.size() < 3) continue;
			selecTable.Fill(d,6,scaleFactor);

			selecTable.Fill(d,7,scaleFactor);

            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            ///////////////////////
            // Getting Gen Event //
            ///////////////////////

/*            TRootGenEvent* genEvt = 0;

            if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data")
            {
                vector<TRootMCParticle*> mcParticles;
                vector<TRootMCParticle*> mcTops;
                mcParticlesMatching_.clear();
                mcParticlesTLV.clear();
                //selectedJetsTLV.clear();
                mcParticles.clear();
                mcTops.clear();

                int leptonPDG, muonPDG = 13, electronPDG = 11;
                leptonPDG = muonPDG;

                genEvt = treeLoader.LoadGenEvent(ievt,false);
                treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticlesMatching_,false);
                if (debug) cout <<"size   "<< mcParticlesMatching_.size()<<endl;
            }
*/

            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////

            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);




            //////////////////////
            // Muon Based Plots //
            //////////////////////
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                float reliso = selectedMuons[selmu]->relPfIso(4, 0.5);
                MSPlot["MuonRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////////
            // Electron Based Plots //
            //////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronPt"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////
            // Jets Based Plots //
            //////////////////////

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

            ///////////////////
            //MET Based Plots//
            ///////////////////

            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["MT_LepMET"]->Fill(MT, datasets[d], true, Luminosity*scaleFactor);

            //////////////////
            //Topology Plots//
            /////////////////



            //////////////////
            //Filling nTuple//
            //////////////////



        } //End Loop on Events

       	tupfile->Close();
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    eventlist.close();

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;

    //////////////////////
    // Selection tables //
    //////////////////////

    //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    selecTable.TableCalculator(  true, true, true, true, true);

    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTable.Write(  outputDirectory+"/FCNC_1L3B"+postfix+"_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);

    fout->cd();

//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin();
            it != MSPlot.end();
            it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Write(fout, name, true, pathPNG, "pdf");
    }


    TDirectory* th2dir = fout->mkdir("Histos2D");
    th2dir->cd();


    for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {

        TH2F *temp = it->second;
        temp->Write();
    }
    delete fout;
    cout << "It took " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}

int Factorial(int N = 1)
{
    int fact = 1;
    for( int i=1; i<=N; i++ )
        fact = fact * i;  // OR fact *= i;
    return fact;
}

float Sphericity(vector<TLorentzVector> parts )
{
    if(parts.size()>0)
    {
        double spTensor[3*3] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
        int counter = 0;
        float tensorNorm = 0, y1 = 0, y2 = 0, y3 = 0;

        for(int tenx = 0; tenx < 3; tenx++)
        {
            for(int teny = 0; teny < 3; teny++)
            {
                for(int selpart = 0; selpart < parts.size(); selpart++)
                {

                    spTensor[counter] += ((parts[selpart][tenx])*(parts[selpart][teny]));
//                    if((tenx == 0 && teny == 2) || (tenx == 2 && teny == 1))
//                    {
//                    cout << "nan debug term " << counter+1 << ": " << (parts[selpart][tenx])*(parts[selpart][teny]) << endl;
//                    cout << "Tensor Building Term " << counter+1 << ": " << spTensor[counter] << endl;
//                    }
                    if(tenx ==0 && teny == 0)
                    {
                        tensorNorm += parts[selpart].Vect().Mag2();
                    }
                }
                if((tenx == 0 && teny == 2) || (tenx == 2 && teny == 1))
                {
//                    cout << "Tensor term pre-norm " << counter+1 << ": " << spTensor[counter] << endl;
                }
                spTensor[counter] /= tensorNorm;
//                cout << "Tensor Term " << counter+1 << ": " << spTensor[counter] << endl;
                counter++;
            }
        }
        TMatrixDSym m(3, spTensor);
        //m.Print();
        TMatrixDSymEigen me(m);
        TVectorD eigenval = me.GetEigenValues();
        vector<float> eigenVals;
        eigenVals.push_back(eigenval[0]);
        eigenVals.push_back(eigenval[1]);
        eigenVals.push_back(eigenval[2]);
        sort(eigenVals.begin(), eigenVals.end());
        //cout << "EigenVals: "<< eigenVals[0] << ", " << eigenVals[1] << ", " << eigenVals[2] << ", " << endl;
        float sp = 3.0*(eigenVals[0] + eigenVals[1])/2.0;
        //cout << "Sphericity: " << sp << endl;
        return sp;
    }
    else
    {
        return 0;
    }
}
float Centrality(vector<TLorentzVector> parts)
{
    float E = 0, ET = 0;
    for(int selpart = 0; selpart < parts.size(); selpart++)
    {
        E += parts[selpart].E();
        ET += parts[selpart].Et();
    }
    return ET/E;
}
