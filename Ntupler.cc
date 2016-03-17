//////////////////////////////////////////////////////////////////////////////////////////////
////         Ntupler for FCNC search in top decays at 13 TEV            ///
////              --- Kevin Deroover                                                        ///
//////////////////////////////////////////////////////////////////////////////////////////////

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
    const float Luminosity        =  strtod(argv[11], NULL);
    string fileName                 = argv[12];
    // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list) 
    const string channel            = argv[argc-4];
    const int JobNum                = strtol(argv[argc-3], NULL, 10);
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);

    vector<string> vecfileNames;
    for(int args = 12; args < argc-4; args++)
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
    cout << "Channel: " << channel << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "JobNum: " << JobNum << endl;
    cout << "----------------------------------------" << endl;

    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();

    //Initializing event counts etc.
    int passed = 0;
    int eventCount = 0;

    //Initializing CSVv2 b-tag WP
	  float workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations.
	  float workingpointvalue_Medium = 0.800;//working points updated to 2016 BTV-POG recommendations.
	  float workingpointvalue_Tight = 0.935;//working points updated to 2016 BTV-POG recommendations.

    clock_t start = clock();
    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FCNC_1L3B search ! "           << endl;
    cout << "*************************************************************" << endl;



    ///////////////////////////////////////////////////////////////
    // Initialize scale&reweight-handlings
    //////////////////////////////////////////////////////////////
    bool bTagReweight = true;
    bool bLeptonSF = true;
    bool bTagReweight_PreReweighting = true; //Needs to be set only once to true in order to produce the EtaPtHistos
    bool applyJES = true;
    bool applyJER = true;
    
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
    bool Muon = false;
    bool Electron = false;
    string electronID = "";
    string btagger = "CSVM";
    bool printTriggers = false;
    bool applyTriggers = true;
    string channelpostfix = "";

    //Setting Lepton Channels
    if(channel == "Mu") Muon = true;
    else if (channel == "El") Electron = true;
    else 
    {
        cerr<<"Correct lepton Channel not selected."<<endl;
        exit(1);
    }
   	FILE* eventlist;

    if(Muon && !Electron)
    {
        cout << " --> Using the Muon channel..." << endl;
        channelpostfix = "_Mu";
    	  eventlist = fopen("EventInfo_mu.txt","w");
    	  electronID = "Loose";
    }
    else if(!Muon && Electron)
    {
        cout << " --> Using the Electron channel..." << endl;
        channelpostfix = "_El";
    	  eventlist = fopen("EventInfo_El.txt","w");
    	  electronID = "Medium";
    }
    else 
    {
        cerr<<"Correct lepton Channel not selected."<<endl;
        exit(1);
    }

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
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
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
            bTagCalib = new BTagCalibration("CSVv2","../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_76X_combToMujets.csv");
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
                btwt_mujets_central = new BTagWeightTools(bTagReader_mujets_central,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",false,30,999,2.4);
//                btwt_mujets_up = new BTagWeightTools(bTagReader_mujets_up,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_up.root",false,30,999,2.4);
//                btwt_mujets_down = new BTagWeightTools(bTagReader_mujets_down,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_down.root",false,30,999,2.4);
//                btwt_ttbar_central = new BTagWeightTools(bTagReader_ttbar_central,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_central.root",false,30,999,2.4);
//                btwt_ttbar_up = new BTagWeightTools(bTagReader_ttbar_up,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_up.root",false,30,999,2.4);
//                btwt_ttbar_down = new BTagWeightTools(bTagReader_ttbar_down,"BTagHistosPtEta/HistosPtEta_"+dName+"_ttbar_down.root",false,30,999,2.4);
            }
        }       
    }

    /////////////////////////////////////////////////
    //                   Lepton SF                 //
    /////////////////////////////////////////////////
    MuonSFWeight* muonSFWeightID_TT;   
    MuonSFWeight* muonSFWeightIso_TT;
    MuonSFWeight* muonSFWeightTrigC_TT;
    MuonSFWeight* muonSFWeightTrigD1_TT;
    MuonSFWeight* muonSFWeightTrigD2_TT;


    ElectronSFWeight* electronSFWeight; 
    if(bLeptonSF)
    {
        if(Muon)
        { 
            muonSFWeightID_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
            muonSFWeightIso_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
            muonSFWeightTrigC_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runC_IsoMu20_OR_IsoTkMu20_PtEtaBins/abseta_pt_ratio", true, false, false);
            muonSFWeightTrigD1_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio", true, false, false);
            muonSFWeightTrigD2_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio", true, false, false);
        }
        else if(Electron)
        {
                if(electronID == "Loose") electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/CutBasedID_LooseWP_76X_18Feb.txt_SF2D.root","EGamma_SF2D",true,false);
                if(electronID == "Medium") electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/CutBasedID_MediumWP_76X_18Feb.txt_SF2D.root","EGamma_SF2D",true,false);
                if(electronID == "Tight") electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root","EGamma_SF2D",true,false);
                if(electronID == "Veto") electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/CutBasedID_VetoWP_76X_18Feb.txt_SF2D.root","EGamma_SF2D",true,false);
        }
    }

    LumiReWeighting puSFs;
    puSFs = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data76X_25ns-Run246908-260627Cert.root", "pileup", "pileup");    

    ////////////////////////////
    /// Initialise trigger /
    ////////////////////////////

    if (debug) cout << "Initializing trigger" << endl;    
    //Trigger* trigger = new Trigger(hasMuon, hasElectron, trigSingleLep, trigDoubleLep);
    Trigger* trigger = 0;
    if(applyTriggers) trigger = new Trigger(Muon, Electron, true, false);

    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////
    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    SourceDate *strdate = new SourceDate();
    string date_str = strdate->ReturnDateStr();

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


    /////////////////////////////////
    //       Loop on datasets      //
    /////////////////////////////////
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset

        bool nlo = false;
        bool isData = false;
        double nloSF = 1; 
        double sumWeights = 0;

        if(dName.find("NLO") != std::string::npos || dName.find("nlo") !=std::string::npos || dName.find("amc") !=std::string::npos) nlo = true;
        else nlo = false;

        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;
        float normfactor = datasets[d]->NormFactor();

        ////////////////////////////////////////////////////////////
        // Setup Date string and nTuple for output  
        ///////////////////////////////////////////////////////////
        if(debug)cout<<"date print"<<endl;

        string channel_dir = "NtuplerOutput/Ntuples"+channelpostfix;
        string date_dir = channel_dir+"/Ntuples_" + date_str +"/";
        int mkdirstatus__ = mkdir("NtuplerOutput",0777);
        mkdirstatus__ = mkdir(channel_dir.c_str(),0777);
        mkdirstatus__ = mkdir(date_dir.c_str(),0777);

        string jobNumString = static_cast<ostringstream*>( &(ostringstream() << JobNum) )->str();
        string Ntupname = date_dir +"FCNC_1L3B_" +postfix + channelpostfix + "_" + jobNumString + ".root";
        string Ntuptitle_AdvancedVars = "AdvancedVarsTree";
        string Ntuptitle_ObjectVars = "ObjectVarsTree";
        string Ntuptitle_NtupleInfo = "NtupleInfoTree";
        string Ntuptitle_Weights = "Weights";

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
	      tupfile->cd();
	      TTree* tup = new TTree(Ntuptitle_AdvancedVars.c_str(), Ntuptitle_AdvancedVars.c_str());
        TTree* tup_ObjectVars      = new TTree(Ntuptitle_ObjectVars.c_str(), Ntuptitle_ObjectVars.c_str());
	      TTree* tup_ntupleinfo      = new TTree(Ntuptitle_NtupleInfo.c_str(), Ntuptitle_NtupleInfo.c_str());
        TTree* tup_Weights      = new TTree(Ntuptitle_Weights.c_str(), Ntuptitle_Weights.c_str());

        //////////////////////////////////////////////////
        // Pre-event loop definitions
        /////////////////////////////////////////////////
        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;
        string previousFilename = "";
        int event_start = startEvent;

        if (debug) cout << " - Loop over events " << endl;

        double end_d = ending;
        if(endEvent > ending) end_d = ending;
        else end_d = endEvent;
        double EqLumi = 1.;
        int nEvents = end_d - event_start;
        if(dName.find("Data")==string::npos) EqLumi = nEvents/xSect;
        cout <<"Equivalent Lumi: "<<  EqLumi  << endl;

        ///////////////////////////////////////
        // Variables for ntuples
        ///////////////////////////////////////
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        Int_t cutstep[10]; //0: no cut, 1: PV cleaning, 2: trigger, 3: lepton selection, 4: loose lepton veto, 5: nb jets, 6: nb b-jets
        Int_t nCuts = 7; //REDEFINE if ncuts change
        Int_t nofPosWeights = 0;
        Int_t nofNegWeights = 0;
        Int_t sumW = 0; 
        Int_t JERon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Int_t JESon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Double_t Luminosity_ = Luminosity;

        //Weights
        Double_t puSF;
        Double_t fleptonSF;
        Double_t btagWeight_mujets_central;
        Double_t btagWeight_mujets_up;
        Double_t btagWeight_mujets_down;
        Double_t nloWeight;// for amc@nlo samples
        Double_t weight1;
        Double_t weight2;
        Double_t weight3;
        Double_t weight4;
        Double_t weight5;
        Double_t weight6;
        Double_t weight7;
        Double_t weight8; 
        Double_t MuonIDSF; 
        Double_t MuonIsoSF; 
        Double_t MuonTrigSF;

 
      
	      // variables for electrons
        Int_t nElectrons;
        Double_t pt_electron;
        Double_t phi_electron;
        Double_t eta_electron;
        Double_t eta_superCluster_electron;
        Double_t E_electron;
        Double_t d0_electron;
        Double_t d0BeamSpot_electron;
        Double_t chargedHadronIso_electron;
        Double_t neutralHadronIso_electron;
        Double_t photonIso_electron;
        Double_t pfIso_electron;
        Double_t charge_electron;
        Double_t sigmaIEtaIEta_electron;
	      Double_t deltaEtaIn_electron;
	      Double_t deltaPhiIn_electron;
	      Double_t hadronicOverEm_electron;
	      Int_t missingHits_electron;
	      Bool_t passConversion_electron;
	      Bool_t isId_electron;
	      Bool_t isIso_electron;      
        Bool_t isEBEEGap; 
	      Double_t sf_electron;
      
        //variable for muons
        Int_t nMuons;
        Double_t pt_muon;
        Double_t phi_muon;
        Double_t eta_muon;
        Double_t E_muon;
        Double_t d0_muon;
        Double_t d0BeamSpot_muon;
        Double_t chargedHadronIso_muon;
        Double_t neutralHadronIso_muon;
        Double_t photonIso_muon;
        Double_t relIso_muon;
	      Bool_t isId_muon;
	      Bool_t isIso_muon;
        Double_t pfIso_muon;
        Double_t charge_muon;
  
        //variable for jets 
        Int_t nJets;
	      Int_t nJets_CSVL; 
	      Int_t nJets_CSVM; 
	      Int_t nJets_CSVT;
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Double_t charge_jet[20];
        Double_t bdisc_jet[20];
        Double_t cdiscCvsL_jet[20]; 
	      Double_t cdiscCvsB_jet[20];
	      Int_t jet_matchedMC_pdgID[20];
	      Int_t jet_matchedMC_motherpdgID[20];
	      Int_t jet_matchedMC_grannypdgID[20];
      
        // met 
        Double_t met_Pt; 
	      Double_t met_Phi; 
	      Double_t met_Eta;
	      
	      //Reconstructed variables
	      Double_t Mbb;
	      Double_t MTlepmet;
	      Double_t MLepTop_GenMatch;
	      Double_t MHadTop_GenMatch;
	      Double_t EtaLepTop_GenMatch;
	      Double_t EtaHadTop_GenMatch;
	      Double_t MassW_GenMatch;
	      Double_t EtaW_GenMatch;
	      Double_t dR_lepJet_min;
	      Double_t MLepTop;
	      Double_t MHadTop;
	      Double_t EtaLepTop;
	      Double_t EtaHadTop;
	      Double_t MassW;
	      Double_t EtaW;

        // global data set variables
        tup_ntupleinfo->Branch("Luminosity_",&Luminosity_,"Luminosity_/D");  
        tup_ntupleinfo->Branch("I_nofPosWeights",&nofPosWeights,"nofPosWeights/I");  
	      tup_ntupleinfo->Branch("I_nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        tup_ntupleinfo->Branch("I_nEvents" , &nEvents, "nEvents/I"); 
        tup_ntupleinfo->Branch("I_sumW", &sumW, "sumW/I");
        tup_ntupleinfo->Branch("I_nCuts",&nCuts, "nCuts/I"); 
        tup_ntupleinfo->Branch("I_cutstep",&cutstep,"cutstep[nCuts]/I");
        tup_ntupleinfo->Branch("I_JERon",&JERon,"JERon/I"); 
        tup_ntupleinfo->Branch("I_JESon", &JESon, "JESon/I");
        tup_ntupleinfo->Branch("workingpointvalue_Loose", &workingpointvalue_Loose, "workingpointvalue_Loose/D"); 
        tup_ntupleinfo->Branch("workingpointvalue_Medium", &workingpointvalue_Medium, "workingpointvalue_Medium/D");
        tup_ntupleinfo->Branch("workingpointvalue_Tight", &workingpointvalue_Tight, "workingpointvalue_Tight/D");
        tup_ntupleinfo->Branch("I_run_num",&run_num,"run_num/I");
        tup_ntupleinfo->Branch("I_evt_num",&evt_num,"evt_num/I");
        tup_ntupleinfo->Branch("I_lumi_num",&lumi_num,"lumi_num/I");
        tup_ObjectVars->Branch("I_nvtx",&nvtx,"nvtx/I");
        tup_ObjectVars->Branch("I_npu",&npu,"npu/I");

        // Weights
        tup_Weights->Branch("fleptonSF",&fleptonSF,"fleptonSF/D"); //Contains, if muon, the  isoSF, idSF & trigSF
        tup_Weights->Branch("puSF",&puSF,"puSF/D");  
        tup_Weights->Branch("btagWeight_mujets_central",&btagWeight_mujets_central,"btagWeight_mujets_central/D"); 
        tup_Weights->Branch("btagWeight_mujets_up",&btagWeight_mujets_up,"btagWeight_mujets_up/D");  
        tup_Weights->Branch("btagWeight_mujets_down",&btagWeight_mujets_down,"btagWeight_mujets_down/D"); 
        tup_Weights->Branch("nloWeight",&nloWeight,"nloWeight/D"); 
        tup_Weights->Branch("weight1",&weight1,"weight1/D");  
        tup_Weights->Branch("weight2",&weight2,"weight2/D");  
        tup_Weights->Branch("weight3",&weight3,"weight3/D"); 
        tup_Weights->Branch("weight4",&weight4,"weight4/D");  
        tup_Weights->Branch("weight5",&weight5,"weight5/D"); 
        tup_Weights->Branch("weight6",&weight6,"weight6/D");  
        tup_Weights->Branch("weight7",&weight7,"weight7/D"); 
        tup_Weights->Branch("weight8",&weight8,"weight8/D");  
        tup_Weights->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF/D");  
        tup_Weights->Branch("MuonIsoSF",&MuonIsoSF,"MuonIsoSF/D");  
        tup_Weights->Branch("MuonTrigSF",&MuonTrigSF,"MuonTrigSF/D");  

        // electrons
        tup_ObjectVars->Branch("pt_electron",&pt_electron,"pt_electron/D");
        tup_ObjectVars->Branch("phi_electron",&phi_electron,"phi_electron/D");
        tup_ObjectVars->Branch("eta_electron",&eta_electron,"eta_electron/D");
        tup_ObjectVars->Branch("eta_superCluster_electron",&eta_superCluster_electron,"eta_superCluster_electron/D");
        tup_ObjectVars->Branch("E_electron",&E_electron,"E_electron/D");
        tup_ObjectVars->Branch("chargedHadronIso_electron",&chargedHadronIso_electron,"chargedHadronIso_electron/D");
        tup_ObjectVars->Branch("neutralHadronIso_electron",&neutralHadronIso_electron,"neutralHadronIso_electron/D");
        tup_ObjectVars->Branch("photonIso_electron",&photonIso_electron,"photonIso_electron/D");
        tup_ObjectVars->Branch("pfIso_electron",&pfIso_electron,"pfIso_electron/D");
        tup_ObjectVars->Branch("charge_electron",&charge_electron,"charge_electron/D");
        tup_ObjectVars->Branch("d0_electron",&d0_electron,"d0_electron/D");
        tup_ObjectVars->Branch("d0BeamSpot_electron",&d0BeamSpot_electron,"d0BeamSpot_electron/D");
        tup_ObjectVars->Branch("sigmaIEtaIEta_electron",&sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron/D");
        tup_ObjectVars->Branch("deltaEtaIn_electron",&deltaEtaIn_electron,"deltaEtaIn_electron/D");
        tup_ObjectVars->Branch("deltaPhiIn_electron",&deltaPhiIn_electron,"deltaPhiIn_electron/D");
        tup_ObjectVars->Branch("hadronicOverEm_electron",&hadronicOverEm_electron,"hadronicOverEm_electron/D");
        tup_ObjectVars->Branch("I_missingHits_electron",&missingHits_electron,"missingHits_electron/I");
        tup_ObjectVars->Branch("I_passConversion_electron",&passConversion_electron,"passConversion_electron/O)");
        tup_ObjectVars->Branch("I_isId_electron",&isId_electron,"isId_electron/O)");
        tup_ObjectVars->Branch("I_isIso_electron",&isIso_electron,"isIso_electron/O)");
        tup_ObjectVars->Branch("I_isEBEEGap",&isEBEEGap,"isEBEEGap/O)");
      

        // muons
        tup_ObjectVars->Branch("pt_muon",&pt_muon,"pt_muon/D");
        tup_ObjectVars->Branch("phi_muon",&phi_muon,"phi_muon/D");
        tup_ObjectVars->Branch("eta_muon",&eta_muon,"eta_muon/D");
        tup_ObjectVars->Branch("E_muon",&E_muon,"E_muon/D");
        tup_ObjectVars->Branch("chargedHadronIso_muon",&chargedHadronIso_muon,"chargedHadronIso_muon/D");
        tup_ObjectVars->Branch("neutralHadronIso_muon",&neutralHadronIso_muon,"neutralHadronIso_muon/D");
        tup_ObjectVars->Branch("photonIso_muon",&photonIso_muon,"photonIso_muon/D");
        tup_ObjectVars->Branch("I_isId_muon",&isId_muon,"isId_muon/O");
        tup_ObjectVars->Branch("I_isIso_muon",&isIso_muon,"isIso_muon/O");
        tup_ObjectVars->Branch("pfIso_muon",&pfIso_muon,"pfIso_muon/D");
        tup_ObjectVars->Branch("charge_muon",&charge_muon,"charge_muon/D");
        tup_ObjectVars->Branch("d0_muon",&d0_muon,"d0_muon/D");
        tup_ObjectVars->Branch("d0BeamSpot_muon",&d0BeamSpot_muon,"d0BeamSpot_muon/D");

        // jets
        tup_ObjectVars->Branch("I_nJets",&nJets,"nJets/I");
        tup_ObjectVars->Branch("I_nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
        tup_ObjectVars->Branch("I_nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
        tup_ObjectVars->Branch("I_nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
        tup_ObjectVars->Branch("pt_jet",&pt_jet,"pt_jet[nJets]/D");
        tup_ObjectVars->Branch("phi_jet",&phi_jet,"phi_jet[nJets]/D");
        tup_ObjectVars->Branch("eta_jet",&eta_jet,"eta_jet[nJets]/D");
        tup_ObjectVars->Branch("E_jet",&E_jet,"E_jet[nJets]/D");
        tup_ObjectVars->Branch("charge_jet",&charge_jet,"charge_jet[nJets]/D");	    
        tup_ObjectVars->Branch("bdisc_jet",&bdisc_jet,"bdisc_jet[nJets]/D");
        tup_ObjectVars->Branch("cdiscCvsL_jet",&cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
        tup_ObjectVars->Branch("cdiscCvsB_jet",&cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
        tup_ObjectVars->Branch("jet_matchedMC_pdgID",&jet_matchedMC_pdgID,"jet_matchedMC_pdgID[nJets]/D");
        tup_ObjectVars->Branch("jet_matchedMC_motherpdgID",&jet_matchedMC_motherpdgID,"jet_matchedMC_motherpdgID/D");
        tup_ObjectVars->Branch("jet_matchedMC_grannypdgID",&jet_matchedMC_grannypdgID,"jet_matchedMC_grannypdgID[nJets]/D");
       
        // met 
        tup_ObjectVars->Branch("met_Pt", &met_Pt, "met_Pt/D"); 
        tup_ObjectVars->Branch("met_Eta", &met_Eta,"met_Eta/D"); 
        tup_ObjectVars->Branch("met_Phi", &met_Phi, "met_Phi/D"); 

        // Advanced variables
        tup->Branch("Mbb",&Mbb,"Mbb/D");
        tup->Branch("MTlepmet",&MTlepmet,"MTlepmet/D");
        tup->Branch("MLepTop_GenMatch",&MLepTop_GenMatch,"MLepTop_GenMatch/D");
        tup->Branch("MHadTop_GenMatch",&MHadTop_GenMatch,"MHadTop_GenMatch/D");
        tup->Branch("EtaLepTop_GenMatch",&EtaLepTop_GenMatch,"EtaLepTop_GenMatch/D");
        tup->Branch("EtaHadTop_GenMatch",&EtaHadTop_GenMatch,"EtaHadTop_GenMatch/D");
        tup->Branch("MassW_GenMatch",&MassW_GenMatch,"MassW_GenMatch/D");
        tup->Branch("EtaW_GenMatch",&EtaW_GenMatch,"EtaW_GenMatch/D");
        tup->Branch("dR_lepJet_min",&dR_lepJet_min,"dR_lepJet_min/D");
        tup->Branch("MLepTop",&MLepTop,"MLepTop/D");
        tup->Branch("MHadTop",&MHadTop,"MHadTop/D");
        tup->Branch("EtaLepTop",&EtaLepTop,"EtaLepTop/D");
        tup->Branch("EtaHadTop",&EtaHadTop,"EtaHadTop/D");
        tup->Branch("MassW",&MassW,"MassW/D");
        tup->Branch("EtaW",&EtaW,"EtaW/D");

        if(debug)cout<<"created ntuples"<<endl;
        ///////////////////////////////////////////////////////////////
        // JEC
        ///////////////////////////////////////////////////////////////
        vector<JetCorrectorParameters> vCorrParam;

        if(dName.find("Data")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
        }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt");

        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);

        /// book triggers
        if (applyTriggers) { trigger->bookTriggers(isData);}




        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
        cout <<"Will run over "<<  nEvents << " events..."<<endl;    cout <<"Starting event = = = = "<< event_start  << endl;
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            if(debug)
            {
                cout << " " << endl;
                cout << "------------NEW EVENT: " << ievt << " --------------" << endl;
            }

            double ievt_d = ievt;

            if(ievt%10000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            puSF = 1;
            fleptonSF = 1;
            btagWeight_mujets_central = 1;
            btagWeight_mujets_up = 1;
            btagWeight_mujets_down = 1;
            nloWeight = 1;// for amc@nlo samples
            weight1 = 1;
            weight2 = 1;
            weight3 = 1;
            weight4 = 1;
            weight5 = 1;
            weight6 = 1;
            weight7 = 1;
            weight8 = 1; 
            MuonIDSF = 1; 
            MuonIsoSF = 1; 
            MuonTrigSF = 1;

            if(debug)cout<<"before tree load"<<endl;
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            if(debug)cout<<"after tree load"<<endl;

            //////////////////////////////////////////////////////////////////
            ///  Include trigger set up here when using data
            //////////////////////////////////////////////////////////////////
            datasets[d]->eventTree()->LoadTree(ievt);
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            int currentRun = event->runId();
	          run_num = event->runId(); 
	          evt_num = event->eventId();
	          lumi_num=event->lumiBlockId(); 
	          nvtx = vertex.size();
	          npu = (int) event->nTruePU(); 

            bool fileChanged = false;
            bool runChanged = false;


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
            if(applyJER)
            {
                if(dName.find("Data")==string::npos)//JER smearing only feasible for non-data samples
                {
                    //JER
                    if(applyJER)
                    {
                        JERon = doJERShift;
                        if(doJERShift == 1)
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
                        else if(doJERShift == 2)
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
                        else
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
                    }
                    else JERon = -1;

                    // JES
                    if(applyJES)
                    {
                        
                        JESon = doJESShift;
                        if (doJESShift == 1)
                            jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
                        else if (doJESShift == 2)
                            jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
                        
                            jetTools->correctJets(init_jets,event->fixedGridRhoFastjetAll() ,false);
                    }
                    else JESon = -1;

                }
            }
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
		            }
		            else 
		            {
		                if(nlo) hasNegWeight = true;
		                nofNegWeights++; 
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
                }
                 else
                 {
                    if(nlo) hasNegWeight = true;
                    nofNegWeights++;
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
          	    sumWeights += mc_baseweight;

            }
            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            //define object containers
            vector<TRootElectron*> selectedElectrons;
            vector<TRootPFJet*>    selectedJets;
            vector<TRootPFJet*>    selectedOrigJets;
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
		        vector<TLorentzVector> selectedMuonsTLV, selectedElectronsTLV, metsTLV, selectedJetsTLV, selectedBJetsTLV, selectedLightJetsTLV, selectedLeptonsTLV;
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
				        selectedOrigJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Tight"); // ApplyJetId
				        if (debug)cout<<"Getting Tight Muons"<<endl;
				        selectedMuons                                       = r2selection.GetSelectedMuons(30,2.1,0.15, "Tight", "Spring15"); //Selected
				        if (debug)cout<<"Getting Loose Electrons"<<endl;
				        selectedElectrons                                   = r2selection.GetSelectedElectrons(20,2.5,"Loose", "Spring15_25ns", true); //Vetoed  
				        if (debug)cout<<"Getting Loose Muons"<<endl;
				        selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20,"Loose","Spring15"); //Vetoed         
            }
            else if (Electron)
            {
				        if (debug)cout<<"Getting Jets"<<endl;
				        selectedOrigJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Tight"); // ApplyJetId
				        if (debug)cout<<"Getting Loose Muons"<<endl;
				        selectedMuons                                       = r2selection.GetSelectedMuons(20, 2.4, 0.20,"Loose","Spring15"); //Vetoed
				        if (debug)cout<<"Getting Electrons"<<endl;
				        selectedElectrons                                   = r2selection.GetSelectedElectrons(30,2.1,electronID, "Spring15_25ns", true); //Selected                       
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
            
            if(debug)
            {
                cout << "Number of jets: " << selectedOrigJets.size() << endl;
                cout << "Number of selected muons: " << nMu << endl;
                cout << "Number of selected electrons: " << nEl << endl;
                cout << "Number of loose muons: " << nLooseMu << endl;
                cout << "Number of loose electrons: " << nLooseEl << endl;
            }

            //////////////////////////////////////////////////
            // Apply scale factors
            //////////////////////////////////////////////////
            if(dName.find("Data")!=string::npos) //If sample is data, no PU reweighting
            {
                puSF=1;
            }
            else
            {
                puSF = puSFs.ITweight( (int)event->nTruePU() );
            }
            if(debug) cout << "puSF: " << puSF << endl;
            /////////////////////////////////////////////////
            //                   Lepton SF                 //
            /////////////////////////////////////////////////
            float fleptonSF = 1;
            if(bLeptonSF)
            {
                if(Muon && nMu>0){
                    fleptonSF = muonSFWeightID_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0) * muonSFWeightIso_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    MuonIDSF = muonSFWeightID_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    MuonIsoSF = muonSFWeightIso_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                }
                else if(Electron && nEl>0){
                    fleptonSF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                }
            }

            float trigSFC = 1;
            float trigSFD1 = 1;
            float trigSFD2 = 1;
            if(bLeptonSF)
            {
                if(dName.find("Data")==string::npos && Muon && nMu>0)
                {
                    trigSFC = muonSFWeightTrigC_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    trigSFD1 = muonSFWeightTrigD1_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    trigSFD2 = muonSFWeightTrigD2_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);       
                    MuonTrigSF =( (trigSFC*17.2) + (trigSFD1*923.88) + (trigSFD2*1639.20) )/Luminosity;  
                }
                fleptonSF*=MuonTrigSF;
            }

            if(debug) cout<<"lepton SF:  "<<fleptonSF<<endl;

            /////////////////////////////////////////////////
            //                   Btag SF                    //
            /////////////////////////////////////////////////
            if(bTagReweight && !bTagReweight_PreReweighting)
            {
                if(dName.find("Data")==string::npos) //If sample is data, no PU reweighting
                {
                    if(debug) cout << "Applying b-tag weights " << endl;
                    //                    btagWeight_comb_central =  btwt_comb_central->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_comb_up =  btwt_comb_up->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_comb_down =  btwt_comb_down->getMCEventWeight(selectedOrigJets, false);
                    btagWeight_mujets_central =  btwt_mujets_central->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_mujets_up =  btwt_mujets_up->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_mujets_down =  btwt_mujets_down->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_ttbar_central =  btwt_ttbar_central->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_ttbar_up =  btwt_ttbar_up->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_ttbar_down =  btwt_ttbar_down->getMCEventWeight(selectedOrigJets, false);
                }
            }
            if(debug) cout<<"btag SF:  "<< btagWeight_mujets_central << endl;

            scaleFactor = scaleFactor * puSF * fleptonSF * btagWeight_mujets_central;
            ////////////////////////////////////////////////
            // Pre-baseline initializations
            ////////////////////////////////////////////////           
            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV <<endl;
            //            if (debug) cin.get();

            if(debug) cout << "Past cut 0: NO CUTS" << endl;
            cutstep[0]++; //Order of appearance of cutstep & nCuts is important here

            eventCount++;

            bool trigged = false;
            if ( ! applyTriggers && previousFilename != currentFilename )
            {
                fileChanged = true;
                previousFilename = currentFilename;
                trigged = true;
            }
            
            if (applyTriggers)
            {
              trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTriggers);
              trigged = trigger->checkIfFired();
              
            }
            else trigged = true;
             if(dName.find("NP")!=string::npos)//JER smearing only feasible for non-data samples
            {
                trigged = true;
            }
            //////////////////////////////////////////////////////
            // Applying baseline lepton selection
            //////////////////////////////////////////////////////


            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            if(debug) cout << "Past cut 1: good PV selection" << endl;
            cutstep[1]++; //Order of appearance of cutstep & nCuts is important here

            if(!trigged) continue;
            if(debug) cout << "Past cut 2: Passed trigger" << endl;
            cutstep[2]++; //Order of appearance of cutstep & nCuts is important here

            if (debug)
            {
              	cout <<" applying baseline event selection..."<<endl;
              	cout <<"number of muons: " << nMu <<endl;
              	cout <<"number of electrons: " << nEl <<endl;
            }
            //Apply the lepton, btag and HT selections
            if (Muon && !Electron)
            {
                if  (  !( nMu ==1 && nEl == 0)) continue; // Muon Channel Selection + veto on electrons
                if (selectedMuons[0]->Pt() < 30) continue;
                if (debug)	cout <<"Muon selection passed..."<<endl;
            }
            else if (!Muon && Electron)
            {
                  if  (  !( nMu == 0 && nEl == 1)) continue; // Electron Channel Selection + veto on muons
                  if (selectedElectrons[0]->Pt() < 30) continue;
                  if (debug)	cout <<"Electron selection passed..."<<endl;
            }
            else
            {
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            if(debug) cout << "Past cut 3: Single lepton selected" << endl;
            cutstep[3]++; //Order of appearance of cutstep & nCuts is important here

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
            if(debug) cout << "Past cut 4: Vetoed extra loose leptons" << endl;
            cutstep[4]++; //Order of appearance of cutstep & nCuts is important here
            
            if(Electron)
            {
              pt_electron=selectedElectrons[0]->Pt();
	            phi_electron=selectedElectrons[0]->Phi();
	            eta_electron=selectedElectrons[0]->Eta();
	            eta_superCluster_electron=selectedElectrons[0]->superClusterEta();
	            E_electron=selectedElectrons[0]->E();
	            d0_electron=selectedElectrons[0]->d0();
	            d0BeamSpot_electron=selectedElectrons[0]->d0BeamSpot();
	            chargedHadronIso_electron=selectedElectrons[0]->chargedHadronIso(3);
	            neutralHadronIso_electron=selectedElectrons[0]->neutralHadronIso(3);
	            photonIso_electron=selectedElectrons[0]->photonIso(3);
	            pfIso_electron=selectedElectrons[0]->relPfIso(3,0);
	            charge_electron=selectedElectrons[0]->charge();
	            sigmaIEtaIEta_electron=selectedElectrons[0]->sigmaIEtaIEta();
	            deltaEtaIn_electron=selectedElectrons[0]->deltaEtaIn();
	            deltaPhiIn_electron=selectedElectrons[0]->deltaPhiIn();
	            hadronicOverEm_electron=selectedElectrons[0]->hadronicOverEm();
	            missingHits_electron=selectedElectrons[0]->missingHits();
	            passConversion_electron=selectedElectrons[0]->passConversion();
	            isEBEEGap=selectedElectrons[0]->isEBEEGap();
	          }
	          else if (Muon)
	          {
	              pt_muon=selectedMuons[0]->Pt();
	              phi_muon=selectedMuons[0]->Phi();
	              eta_muon=selectedMuons[0]->Eta();
	              E_muon=selectedMuons[0]->E();
	              d0_muon=selectedMuons[0]->d0();
	              d0BeamSpot_muon=selectedMuons[0]->d0BeamSpot();
	              chargedHadronIso_muon=selectedMuons[0]->chargedHadronIso(4);
	              neutralHadronIso_muon=selectedMuons[0]->neutralHadronIso(4);
	              photonIso_muon=selectedMuons[0]->photonIso(4);
                pfIso_muon=selectedMuons[0]->relPfIso(4,0);
	              charge_muon=selectedMuons[0]->charge();

	          }

            /////////////////////////////////////////////////
            //            Jet lepton cleaning          //
            /////////////////////////////////////////////////
            selectedJets.clear();
		       
		        float clean_dR = 0.4;
            for (int origJets=0; origJets<selectedOrigJets.size(); origJets++)
            {
                if(Electron)
                {
                    if(selectedOrigJets[origJets]->DeltaR(*selectedExtraElectrons[0])>clean_dR)
                    {
                        selectedJets.push_back(selectedOrigJets[origJets]);
                    }
                }                    
                else if(Muon)
                {
                    if(selectedOrigJets[origJets]->DeltaR(*selectedExtraMuons[0])>clean_dR)
                    {
                        selectedJets.push_back(selectedOrigJets[origJets]);
                    }
                }                    
			      }
            if(debug) cout << "Cleaned jets from isolated leptons" << endl;
		
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

			      //////////////////////////////////////
			      // Jet variables //
			      /////////////////////////////////////
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
            nJets_CSVT =  selectedTBJets.size(); 
	          nJets_CSVM =  selectedMBJets.size();
            nJets_CSVL =  selectedLBJets.size();
	          double met_px = mets[0]->Px();
	          double met_py = mets[0]->Py();
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
	          met_Phi = mets[0]->Phi(); 
	          met_Eta = mets[0]->Eta();

         		float HT = 0, H = 0;
         		for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
         		{
              		//Event-level variables
              		HT = HT + selectedJets[seljet1]->Pt();
              		H = H +  selectedJets[seljet1]->P();
         		}


            //////////////////////////////////////////////////////////////////////
            // Cut on nb of jets and b-jets
            //////////////////////////////////////////////////////////////////////
			      if(selectedJets.size() < 2)  continue;
            if(debug) cout << "Past cut 5: Passed number of jets cut" << endl;
            cutstep[5]++; //Order of appearance of cutstep & nCuts is important here

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

		  	    if(selectedLBJets.size() < 2) continue;
	          if (debug)	cout <<"Cut on nb b-jets..."<<endl;
            if(debug) cout << "Past cut 6: Passed cut on number of b-jets" << endl;
            cutstep[6]++; //Order of appearance of cutstep & nCuts is important here


            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            for (unsigned int i = 0; i < selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
            for (unsigned int i = 0; i < selectedMBJets.size(); i++) selectedBJetsTLV.push_back(*selectedMBJets[i]);
            for (unsigned int i = 0; i < selectedLightJets_MWP.size(); i++) selectedLightJetsTLV.push_back(*selectedLightJets_MWP[i]);

            //////////////////////////////////////
            // Peeking at the MC info 
            /////////////////////////////////////
            for( unsigned int nj = 0; nj < selectedJets.size(); nj++)
            {
            	      jet_matchedMC_pdgID[nj] = -999;
	                  jet_matchedMC_motherpdgID[nj] = -999;
	                  jet_matchedMC_grannypdgID[nj] = -999;
            }
            
            pair<unsigned int, unsigned int> leptonicBJet_ = pair<unsigned int,unsigned int>(9999,9999);// First one is jet number, second one is mcParticle number
            pair<unsigned int, unsigned int> hadronicBJet_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
                  
            int pdgID_top = 6; //top quark
                  
            bool Posleptonmatched = false;
            bool Negleptonmatched = false;
            
            if(dName != "data" && dName != "Data" && dName != "Data" && dName != "D_ata")
            {
                treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
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
                    
                    //Storing matched MC pdgid's in ntuples
            	      jet_matchedMC_pdgID[JetPartonPair[i].first] = mcParticlesMatching_[j]->type();
	                  jet_matchedMC_motherpdgID[JetPartonPair[i].first] = mcParticlesMatching_[j]->motherType();
	                  jet_matchedMC_grannypdgID[JetPartonPair[i].first] = mcParticlesMatching_[j]->grannyType();
                      
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
            
            /////////////////////////////////////////////////////
            //Defining complexer ntuple variables
            /////////////////////////////////////////////////////
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
			      MTlepmet = MT;


            if(selectedBJetsTLV.size()>=2)   Mbb = (selectedBJetsTLV[0]+selectedBJetsTLV[1]).M();
            if(selectedLightJetsTLV.size()>=2) MassW = (selectedLightJetsTLV[0]+selectedLightJetsTLV[1]).M();


	          tup->Fill();
            tup_ObjectVars->Fill();
            tup_Weights->Fill();


			      ///////////////////////////
			      // Event info //////
			      ///////////////////////////
			      if(Muon && !Electron)
			      fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedMuons[0]->type(), selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), selectedMuons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());
			      else if(!Muon && Electron)
			      fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->type(), selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), selectedElectrons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());



        } //End Loop on Events
        
        sumW = (int) sumWeights;
        
	      if (! isData  ) 
        {
            cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
            cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
            cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents) << endl;
            
            // Determine scale factor due to negative weights
            nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
            cout << "This corresponds to an event scale factor of " << nloSF  << endl; 
        }

        tup_ntupleinfo->Fill();
        tup_ntupleinfo->Print("all");	          
        tup->Print("all");
        tup_ObjectVars->Print("all");
        tup_Weights->Print("all");

	      tupfile->Write();   
       	tupfile->Close();
        delete tupfile;
        cout <<"n events passed  =  "<<passed <<endl;
        cout << "Event Count: " << eventCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();

    } //End Loop on Datasets

    fclose (eventlist);

//    delete btwt_comb_central;
//    delete btwt_comb_up;
//    delete btwt_comb_down;
    delete btwt_mujets_central;
//    delete btwt_mujets_up;
//    delete btwt_mujets_down;
//    delete btwt_ttbar_central;
//    delete btwt_ttbar_up;
//    delete btwt_ttbar_down;



    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}


