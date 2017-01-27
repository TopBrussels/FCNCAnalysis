//////////////////////////////////////////////////////////////////////////////////////////////
////         Ntupler for FCNC search in top decays at 13 TEV            ///
////              --- Kevin Deroover                                                        ///
//////////////////////////////////////////////////////////////////////////////////////////////

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

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"

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

    //Checking passed Arguments to ensure proper execution of MACRO
    if(argc < 14)
    {
        std::cerr << "INVALID INPUT FROM XMLFILE.  CHECK XML IMPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN passed." << std::endl;
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
    int passed_FinalSelection = 0;
    int passed_Step2 = 0;
    int passed_Step3 = 0;
    int passed_Step4 = 0;
    int passed_Step5 = 0;
    int passed_Step6 = 0;
    int passed_Step7 = 0;
    int passed_Step8 = 0;
    int eventCount = 0;

    //Initializing b-tag WP
	  float CSVv2_workingpointvalue_Loose = 0.5426;//working points updated to 2016ReReco BTV-POG recommendations.
	  float CSVv2_workingpointvalue_Medium = 0.8484;//working points updated to 2016ReReco BTV-POG recommendations.
	  float CSVv2_workingpointvalue_Tight = 0.9535;//working points updated to 2016ReReco BTV-POG recommendations.
	  float cMVA_workingpointvalue_Loose = -0.5884;//working points updated to 2016ReReco BTV-POG recommendations.
	  float cMVA_workingpointvalue_Medium = 0.4432;//working points updated to 2016ReReco BTV-POG recommendations.
	  float cMVA_workingpointvalue_Tight = 0.9432;//working points updated to 2016ReReco BTV-POG recommendations.

    clock_t start = clock();
    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FCNC_1L3B search ! "           << endl;
    cout << "*************************************************************" << endl;



    ///////////////////////////////////////////////////////////////
    // Initialize scale&reweight-handlings
    //////////////////////////////////////////////////////////////
    bool bLeptonSF = true;
    bool bTagReweight_FillMChistos = false; //Needs to be set only once to true in order to produce the BTagEtaPtHistos
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
    {
        postfix += "JESMinus";
    }
    if (doJESShift == 2)
    {
        postfix += "JESPlus";
    }
    if (doJERShift == 1)
    {
        postfix += "JERMinus";
    }
    if (doJERShift == 2)
    {
        postfix += "JERPlus";
    }
/*    if (dobTagEffShift == -1)
        postfix= postfix+"_bTagMinus";
    if (dobTagEffShift == 1)
        postfix= postfix+"_bTagPlus";
    if(domisTagEffShift == -1)
        postfix= postfix+"_misTagMinus";
    if(domisTagEffShift == 1)
        postfix= postfix+"_misTagPlus";
*/

    ///////////////////////////////////////
    // Configuration - ALWAYS CHECK THESE!!!!
    ///////////////////////////////////////
    bool debug = false;
    bool Muon = false;
    bool Electron = false;
    string btagger = "CSVv2M"; //Define which b-tagger + WP is used in the SF for the cutflow-table// valable: CSVv2M, cMVAM
    bool printTriggers = false;
    bool applyTriggers = true;
    bool CSVv2nonshape = false;//Put to false if you don't want to use the regular CSVv2 SF's
    bool applyMVAJetComb = true;
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
    }
    else if(!Muon && Electron)
    {
        cout << " --> Using the Electron channel..." << endl;
        channelpostfix = "_El";
    	  eventlist = fopen("EventInfo_El.txt","w");
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
    anaEnv.ElectronCollection = "Electrons_calibratedPatElectrons";
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
    int mkdirstatus_btag = mkdir("BTagHistosPtEta",0777);

    BTagCalibration * bTagCalib_CSVv2;   
//    BTagCalibration * bTagCalib_cMVA;

    BTagCalibrationReader * bTagReader_CSVv2M_mujets_central;
    BTagCalibrationReader * bTagReader_shape;
    BTagCalibrationReader * reader_JESUp;
    BTagCalibrationReader * reader_JESDown;
    BTagCalibrationReader * reader_LFUp;
    BTagCalibrationReader * reader_LFDown;
    BTagCalibrationReader * reader_HFUp;
    BTagCalibrationReader * reader_HFDown;
    BTagCalibrationReader * reader_HFStats1Up;
    BTagCalibrationReader * reader_HFStats1Down;
    BTagCalibrationReader * reader_HFStats2Up;
    BTagCalibrationReader * reader_HFStats2Down;
    BTagCalibrationReader * reader_LFStats1Up;
    BTagCalibrationReader * reader_LFStats1Down;
    BTagCalibrationReader * reader_LFStats2Up;
    BTagCalibrationReader * reader_LFStats2Down;
    BTagCalibrationReader * reader_CFErr1Up;
    BTagCalibrationReader * reader_CFErr1Down;
    BTagCalibrationReader * reader_CFErr2Up;
    BTagCalibrationReader * reader_CFErr2Down;

    BTagWeightTools *btwt_CSVv2M_mujets_central = 0;


        if(dName.find("Data")==string::npos)
        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
        {
            bTagCalib_CSVv2 = new BTagCalibration("CSVv2","../TopTreeAnalysisBase/Calibrations/BTagging/ttH_BTV_CSVv2_13TeV_2016All_36p5_2017_1_10.csv");

            if(CSVv2nonshape) bTagReader_CSVv2M_mujets_central = new BTagCalibrationReader(bTagCalib_CSVv2,BTagEntry::OP_MEDIUM,"mujets","central"); //mujets
            bTagReader_shape = new BTagCalibrationReader(bTagCalib_CSVv2,BTagEntry::OP_RESHAPING,"iterativefit","central"); //reshaping


            // JESUp
            reader_JESUp = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_jes"); // systematics type
            // JESDown
            reader_JESDown = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_jes"); // systematics type

            // LFUp
            reader_LFUp = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_lf"); // systematics type
            // LFDown
            reader_LFDown = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_lf"); // systematics type

            // HFUp
            reader_HFUp = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_hf"); // systematics type
            // HFDown
            reader_HFDown = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_hf"); // systematics type

            // HFStats1Up
            reader_HFStats1Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_hfstats1"); // systematics type
            // HFStats1Down
            reader_HFStats1Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_hfstats1"); // systematics type

            // HFStats2Up
            reader_HFStats2Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_hfstats2"); // systematics type
            // HFStats2Down
            reader_HFStats2Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_hfstats2"); // systematics type

            // LFStats1Up
            reader_LFStats1Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_lfstats1"); // systematics type
            // LFStats1Down
            reader_LFStats1Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_lfstats1"); // systematics type

            // LFStats2Up
            reader_LFStats2Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_lfstats2"); // systematics type
            // LFStats2Down
            reader_LFStats2Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_lfstats2"); // systematics type

            // CFErr1Up
            reader_CFErr1Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_cferr1"); // systematics type
            // CFErr1Down
            reader_CFErr1Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_cferr1"); // systematics type

            // CFErr2Up
            reader_CFErr2Up = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "up_cferr2"); // systematics type
            // CFErr2Down
            reader_CFErr2Down = new BTagCalibrationReader(bTagCalib_CSVv2, // calibration instance
                           BTagEntry::OP_RESHAPING, // operating point
                           "iterativefit", // measurement type
                           "down_cferr2"); // systematics type        


            if(bTagReweight_FillMChistos && CSVv2nonshape)// Need to differentiate BTagWeightTools according to filling the histos and just reading, because of overwriting possibilities in grid submission
            {
                btwt_CSVv2M_mujets_central = new BTagWeightTools(bTagReader_CSVv2M_mujets_central,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_central.root",false,30,670,2.4);
            }
            else if(CSVv2nonshape)
            {
                btwt_CSVv2M_mujets_central = new BTagWeightTools(bTagReader_CSVv2M_mujets_central,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",true,30,670,2.4);
           }
        }       

    /////////////////////////////////////////////////
    //                   Lepton SF                 //
    /////////////////////////////////////////////////
    MuonSFWeight* muonSFWeightID_BCDEF;   
    MuonSFWeight* muonSFWeightID_GH;   
    MuonSFWeight* muonSFWeightIso_BCDEF;
    MuonSFWeight* muonSFWeightIso_GH;
    MuonSFWeight* muonSFWeightTrig_BCDEF;
    MuonSFWeight* muonSFWeightTrig_GH;


    ElectronSFWeight* electronSFWeightID; 
    ElectronSFWeight* electronSFWeightReco; 
    if(bLeptonSF)
    {
        if(Muon)
        { 
            muonSFWeightID_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
            muonSFWeightID_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonID_EfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false);
            muonSFWeightIso_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonIso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
            muonSFWeightIso_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/MuonIso_EfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
            muonSFWeightTrig_BCDEF = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsBCDEF.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
            muonSFWeightTrig_GH = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/SingleMuonTrigger_EfficienciesAndSF_RunsGH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio", true, false, false);
        }
        else if(Electron)
        {
                electronSFWeightID = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/egammaEffi.txt_SF2D_CutBasedMediumID.root","EGamma_SF2D",true,false,false);
                electronSFWeightReco = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/ElectronSF/egammaEffi.txt_SF2D_GsfTrackingEff.root","EGamma_SF2D",true,false,false);
        }
    }

    LumiReWeighting W_puSFs, W_puSFs_Minus, W_puSFs_Plus;
    W_puSFs = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");    
    W_puSFs_Minus = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");    
    W_puSFs_Plus = LumiReWeighting("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");    
    
    
    ////////////////////////////
    /// Initialise trigger /
    ////////////////////////////
    bool isData = false;
    if(dName.find("Data")!=string::npos)
    {
        isData = true;
        applyTriggers = true; //Always apply trigger on data
    }
    else if(dName.find("noTrigger")!=string::npos) applyTriggers = false;

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
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset

        bool nlo = false;
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
        string Ntuptitle_ObjectVars = "ObjectVarsTree";
        string Ntuptitle_NtupleInfo = "NtupleInfoTree";

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");
	      tupfile->cd();
        TTree* tup_ObjectVars      = new TTree(Ntuptitle_ObjectVars.c_str(), Ntuptitle_ObjectVars.c_str());
	      TTree* tup_ntupleinfo      = new TTree(Ntuptitle_NtupleInfo.c_str(), Ntuptitle_NtupleInfo.c_str());

        //////////////////////////////////////////////////
        // Pre-event loop definitions
        /////////////////////////////////////////////////
        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;
        int event_start = startEvent;

        if (debug) cout << " - Loop over events " << endl;

        double end_d = ending;
        if(endEvent > ending) end_d = ending;
        else end_d = endEvent;
        int nEvents = end_d - event_start;

        ///////////////////////////////////////
        // Variables for ntuples
        ///////////////////////////////////////
        Int_t run_num;
        Int_t evt_num;
        Int_t genTTX;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        Double_t cutstep[10]; //0: no cut, 1: PV cleaning, 2:event cleaning, 3: trigger, 4: lepton selection, 5: loose other-flavoured lepton veto, 6: veto extra loose leptons, 7: nb jets, 8: nb b-jets
        Int_t nCuts = 9; //REDEFINE if ncuts change
        Int_t nofPosWeights = 0;
        Int_t nofNegWeights = 0;
        Int_t sumW = 0; 
        Int_t JERon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Int_t JESon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Double_t Luminosity_ = Luminosity;

        //Weights
        Double_t W_puSF;
        Double_t W_puSF_Minus;
        Double_t W_puSF_Plus;
        Double_t W_fleptonSF;
        Double_t W_fleptonSF_Plus;
        Double_t W_fleptonSF_Minus;
        Double_t W_btagWeight_CSVv2M_mujets_central;
        Double_t W_btagWeight_CSVv2M_mujets_up;
        Double_t W_btagWeight_CSVv2M_mujets_down;
        Double_t W_btagWeight_shape;
        Double_t W_btagWeight_shape_up_lf; 
        Double_t W_btagWeight_shape_down_lf; 
        Double_t W_btagWeight_shape_up_hf; 
        Double_t W_btagWeight_shape_down_hf; 
        Double_t W_btagWeight_shape_up_hfstats1; 
        Double_t W_btagWeight_shape_down_hfstats1; 
        Double_t W_btagWeight_shape_up_hfstats2; 
        Double_t W_btagWeight_shape_down_hfstats2; 
        Double_t W_btagWeight_shape_up_lfstats1; 
        Double_t W_btagWeight_shape_down_lfstats1; 
        Double_t W_btagWeight_shape_up_lfstats2; 
        Double_t W_btagWeight_shape_down_lfstats2; 
        Double_t W_btagWeight_shape_up_cferr1; 
        Double_t W_btagWeight_shape_down_cferr1; 
        Double_t W_btagWeight_shape_up_cferr2; 
        Double_t W_btagWeight_shape_down_cferr2; 
        Double_t W_nloWeight;// for amc@nlo samples
        Double_t W_weight1;
        Double_t W_weight2;
        Double_t W_weight3;
        Double_t W_weight4;
        Double_t W_weight5;
        Double_t W_weight6;
        Double_t W_weight7;
        Double_t W_weight8; 
        Double_t W_MuonIDSF; //One of the 3 components for the total muon SF
        Double_t W_MuonIsoSF; //One of the 3 components for the total muon SF
        Double_t W_MuonTrigSF;//One of the 3 components for the total muon SF
        Double_t W_ElectronIDSF; //One of the 2 components for the total electron SF
        Double_t W_ElectronRecoSF; //One of the 2 components for the total electron SF
        Double_t W_MuonIDSF_Plus; //One of the 3 components for the total muon SF
        Double_t W_MuonIsoSF_Plus; //One of the 3 components for the total muon SF
        Double_t W_MuonTrigSF_Plus;//One of the 3 components for the total muon SF
        Double_t W_ElectronIDSF_Plus; //One of the 2 components for the total electron SF
        Double_t W_ElectronRecoSF_Plus; //One of the 2 components for the total electron SF
        Double_t W_MuonIDSF_Minus; //One of the 3 components for the total muon SF
        Double_t W_MuonIsoSF_Minus; //One of the 3 components for the total muon SF
        Double_t W_MuonTrigSF_Minus;//One of the 3 components for the total muon SF
        Double_t W_ElectronIDSF_Minus; //One of the 2 components for the total electron SF
        Double_t W_ElectronRecoSF_Minus; //One of the 2 components for the total electron SF
        Double_t W_TopPtReweighing;
 
      
	      // variables for electrons
        Double_t eta_superCluster_electron;
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
        Bool_t isEBEEGap; 
      
        //variable for muons
        Double_t d0_muon;
        Double_t d0BeamSpot_muon;
        Double_t chargedHadronIso_muon;
        Double_t neutralHadronIso_muon;
        Double_t photonIso_muon;
        Double_t pfIso_muon;
        Double_t charge_muon;
        
        //variable for  leptons
        Double_t pt_lepton;
        Double_t eta_lepton;
        Double_t phi_lepton;
        Double_t E_lepton;
        Int_t LepCharge;
        
        //MC particle variables (affected by TopPtReweighing)
        Double_t MC_TopPt;
        Double_t MC_AntiTopPt;


        //variable for jets 
        Int_t nJets;
	      Int_t nJets_CSVL; 
	      Int_t nJets_CSVM; 
	      Int_t nJets_CSVT;
	      Int_t nJets_cMVAL; 
	      Int_t nJets_cMVAM; 
	      Int_t nJets_cMVAT;
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Double_t charge_jet[20];
        Double_t incl_charge_jet[20];
        Double_t CSVv2[20];
        Double_t cMVA[20];
        Double_t cdiscCvsL_jet[20]; 
	      Double_t cdiscCvsB_jet[20];
	      Double_t jet_matchedMC_pdgID[20];
	      Double_t jet_matchedMC_motherpdgID[20];
	      Double_t jet_matchedMC_grannypdgID[20];
      
        // met 
        Double_t met_Px; 
        Double_t met_Py; 
        Double_t met_Pt; 
	      Double_t met_Phi; 
	      Double_t met_Eta;	      
	      
	      //JetIndices_correctJetComb
	      Int_t TOPTOPLEPHAD_JetIdx_LepTop;
	      Int_t TOPTOPLEPHAD_JetIdx_HadTop;
	      Int_t TOPTOPLEPHAD_JetIdx_W1;
	      Int_t TOPTOPLEPHAD_JetIdx_W2;
	      Int_t TOPTOPLEPHBB_JetIdx_LepTop;
	      Int_t TOPTOPLEPHBB_JetIdx_HadTop;
	      Int_t TOPTOPLEPHBB_JetIdx_H1;
	      Int_t TOPTOPLEPHBB_JetIdx_H2;
	      Int_t TOPHLEPBB_JetIdx_LepTop_hut;
	      Int_t TOPHLEPBB_JetIdx_H1_hut;
	      Int_t TOPHLEPBB_JetIdx_H2_hut;
	      Int_t TOPHLEPBB_JetIdx_LepTop_hct;
	      Int_t TOPHLEPBB_JetIdx_H1_hct;
	      Int_t TOPHLEPBB_JetIdx_H2_hct;
        Double_t MVA_TOPTOPLEPHAD;
        Double_t MVA_TOPTOPLEPHBB;
        Double_t MVA_TOPHLEPBB_hut;
        Double_t MVA_TOPHLEPBB_hct;

        //Variables for signal/background training
	      Double_t HiggsMass_TOPHLEPBB_hut;
	      Double_t HiggsMass_TOPHLEPBB_hct;
	      Double_t HiggsEta_TOPHLEPBB_hut;
	      Double_t HiggsEta_TOPHLEPBB_hct;
	      Double_t TopLepMass_TOPHLEPBB_hut;
	      Double_t TopLepMass_TOPHLEPBB_hct;
        Double_t TopLepPt_TOPHLEPBB_hut;
        Double_t TopLepPt_TOPHLEPBB_hct;
        Double_t TopLepEta_TOPHLEPBB_hut;
        Double_t TopLepEta_TOPHLEPBB_hct;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
        Double_t TopLepHiggsDr_TOPHLEPBB_hut;
        Double_t TopLepHiggsDr_TOPHLEPBB_hct;
        Double_t HiggsBJet1CSVv2_TOPHLEPBB_hut;
        Double_t HiggsBJet1CSVv2_TOPHLEPBB_hct;
        Double_t HiggsBJet2CSVv2_TOPHLEPBB_hut;
        Double_t HiggsBJet2CSVv2_TOPHLEPBB_hct;
        Double_t TopLepBJetCSVv2_TOPHLEPBB_hut;
        Double_t TopLepBJetCSVv2_TOPHLEPBB_hct;
        Double_t TopHadMass_TOPTOPLEPHAD;
        Double_t TopLepMass_TOPTOPLEPHAD;
        Double_t TopLepTopHadDr_TOPTOPLEPHAD;
        Double_t TopLepBJetCSVv2_TOPTOPLEPHAD;
        Double_t TopHadBJetCSVv2_TOPTOPLEPHAD;
        Double_t TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
        Double_t TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
        Double_t HiggsMass_TOPTOPLEPHBB;
        Double_t TopLepMass_TOPTOPLEPHBB;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
        Double_t TopLepHiggsDr_TOPTOPLEPHBB;
        Double_t HiggsBJet1CSVv2_TOPTOPLEPHBB;
        Double_t HiggsBJet2CSVv2_TOPTOPLEPHBB;
        Double_t TopLepBJetCSVv2_TOPTOPLEPHBB;
        Double_t TopHadNonBJetCSVv2_TOPTOPLEPHBB;

        // global data set variables
        tup_ntupleinfo->Branch("I_nofPosWeights",&nofPosWeights,"nofPosWeights/I");  
	      tup_ntupleinfo->Branch("I_nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        tup_ntupleinfo->Branch("I_nEvents" , &nEvents, "nEvents/I"); 
        tup_ntupleinfo->Branch("I_sumW", &sumW, "sumW/I");
        tup_ntupleinfo->Branch("I_nCuts",&nCuts, "nCuts/I"); 
        tup_ntupleinfo->Branch("cutstep",&cutstep,"cutstep[nCuts]/D");
/*        tup_ntupleinfo->Branch("I_JERon",&JERon,"JERon/I"); 
        tup_ntupleinfo->Branch("I_JESon", &JESon, "JESon/I");
        tup_ntupleinfo->Branch("CSVv2_workingpointvalue_Loose", &CSVv2_workingpointvalue_Loose, "CSVv2_workingpointvalue_Loose/D"); 
        tup_ntupleinfo->Branch("CSVv2_workingpointvalue_Medium", &CSVv2_workingpointvalue_Medium, "CSVv2_workingpointvalue_Medium/D");
        tup_ntupleinfo->Branch("CSVv2_workingpointvalue_Tight", &CSVv2_workingpointvalue_Tight, "CSVv2_workingpointvalue_Tight/D");
        tup_ntupleinfo->Branch("cMVA_workingpointvalue_Loose", &cMVA_workingpointvalue_Loose, "cMVA_workingpointvalue_Loose/D"); 
        tup_ntupleinfo->Branch("cMVA_workingpointvalue_Medium", &cMVA_workingpointvalue_Medium, "cMVA_workingpointvalue_Medium/D");
        tup_ntupleinfo->Branch("cMVA_workingpointvalue_Tight", &cMVA_workingpointvalue_Tight, "cMVA_workingpointvalue_Tight/D");
*/
        // Weights
        tup_ObjectVars->Branch("W_fleptonSF",&W_fleptonSF,"W_fleptonSF/D"); //Contains, if muon, the  isoSF, idSF & trigSF
        tup_ObjectVars->Branch("W_fleptonSF_Plus",&W_fleptonSF_Plus,"W_fleptonSF_Plus/D"); //Contains, if muon, the  isoSF, idSF & trigSF
        tup_ObjectVars->Branch("W_fleptonSF_Minus",&W_fleptonSF_Minus,"W_fleptonSF_Minus/D"); //Contains, if muon, the  isoSF, idSF & trigSF
        tup_ObjectVars->Branch("W_puSF",&W_puSF,"W_puSF/D");
        tup_ObjectVars->Branch("W_puSF_Minus",&W_puSF_Minus,"W_puSF_Minus/D");
        tup_ObjectVars->Branch("W_puSF_Plus",&W_puSF_Plus,"W_puSF_Plus/D");
        tup_ObjectVars->Branch("W_btagWeight_CSVv2M_mujets_central",&W_btagWeight_CSVv2M_mujets_central,"W_btagWeight_CSVv2M_mujets_central/D"); 
        tup_ObjectVars->Branch("W_btagWeight_CSVv2M_mujets_up",&W_btagWeight_CSVv2M_mujets_up,"W_btagWeight_CSVv2M_mujets_up/D");  
        tup_ObjectVars->Branch("W_btagWeight_CSVv2M_mujets_down",&W_btagWeight_CSVv2M_mujets_down,"W_btagWeight_CSVv2M_mujets_down/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape",&W_btagWeight_shape,"W_btagWeight_shape/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_lf",&W_btagWeight_shape_up_lf,"W_btagWeight_shape_up_lf/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_lf",&W_btagWeight_shape_down_lf,"W_btagWeight_shape_down_lf/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_hf",&W_btagWeight_shape_up_hf,"W_btagWeight_shape_up_hf/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_hf",&W_btagWeight_shape_down_hf,"W_btagWeight_shape_down_hf/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_hfstats1",&W_btagWeight_shape_up_hfstats1,"W_btagWeight_shape_up_hfstats1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_hfstats1",&W_btagWeight_shape_down_hfstats1,"W_btagWeight_shape_down_hfstats1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_hfstats2",&W_btagWeight_shape_up_hfstats2,"W_btagWeight_shape_up_hfstats2/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_hfstats2",&W_btagWeight_shape_down_hfstats2,"W_btagWeight_shape_down_hfstats2/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_lfstats1",&W_btagWeight_shape_up_lfstats1,"W_btagWeight_shape_up_lfstats1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_lfstats1",&W_btagWeight_shape_down_lfstats1,"W_btagWeight_shape_down_lfstats1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_lfstats2",&W_btagWeight_shape_up_lfstats2,"W_btagWeight_shape_up_lfstats2/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_lfstats2",&W_btagWeight_shape_down_lfstats2,"W_btagWeight_shape_down_lfstats2/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_cferr1",&W_btagWeight_shape_up_cferr1,"W_btagWeight_shape_up_cferr1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_cferr1",&W_btagWeight_shape_down_cferr1,"W_btagWeight_shape_down_cferr1/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_up_cferr2",&W_btagWeight_shape_up_cferr2,"W_btagWeight_shape_up_cferr2/D"); 
        tup_ObjectVars->Branch("W_btagWeight_shape_down_cferr2",&W_btagWeight_shape_down_cferr2,"W_btagWeight_shape_down_cferr2/D"); 
        tup_ObjectVars->Branch("W_nloWeight",&W_nloWeight,"W_nloWeight/D"); 
        tup_ObjectVars->Branch("W_weight1",&W_weight1,"W_weight1/D");  
        tup_ObjectVars->Branch("W_weight2",&W_weight2,"W_weight2/D");  
        tup_ObjectVars->Branch("W_weight3",&W_weight3,"W_weight3/D"); 
        tup_ObjectVars->Branch("W_weight4",&W_weight4,"W_weight4/D");  
        tup_ObjectVars->Branch("W_weight5",&W_weight5,"W_weight5/D"); 
        tup_ObjectVars->Branch("W_weight6",&W_weight6,"W_weight6/D");  
        tup_ObjectVars->Branch("W_weight7",&W_weight7,"W_weight7/D"); 
        tup_ObjectVars->Branch("W_weight8",&W_weight8,"W_weight8/D");  
        tup_ObjectVars->Branch("W_TopPtReweighing",&W_TopPtReweighing,"W_TopPtReweighing/D");  

        tup_ObjectVars->Branch("I_run_num",&run_num,"run_num/I");
        tup_ObjectVars->Branch("I_genTTX",&genTTX,"genTTX/I");
        tup_ObjectVars->Branch("I_evt_num",&evt_num,"evt_num/I");
        tup_ObjectVars->Branch("I_lumi_num",&lumi_num,"lumi_num/I");
        tup_ObjectVars->Branch("I_nvtx",&nvtx,"nvtx/I");
        tup_ObjectVars->Branch("I_npu",&npu,"npu/I");


        // electrons
        tup_ObjectVars->Branch("eta_superCluster_electron",&eta_superCluster_electron,"eta_superCluster_electron/D");
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
        tup_ObjectVars->Branch("I_isEBEEGap",&isEBEEGap,"isEBEEGap/O)");
      

        // muons
        tup_ObjectVars->Branch("chargedHadronIso_muon",&chargedHadronIso_muon,"chargedHadronIso_muon/D");
        tup_ObjectVars->Branch("neutralHadronIso_muon",&neutralHadronIso_muon,"neutralHadronIso_muon/D");
        tup_ObjectVars->Branch("photonIso_muon",&photonIso_muon,"photonIso_muon/D");
        tup_ObjectVars->Branch("pfIso_muon",&pfIso_muon,"pfIso_muon/D");
        tup_ObjectVars->Branch("charge_muon",&charge_muon,"charge_muon/D");
        tup_ObjectVars->Branch("d0_muon",&d0_muon,"d0_muon/D");
        tup_ObjectVars->Branch("d0BeamSpot_muon",&d0BeamSpot_muon,"d0BeamSpot_muon/D");

        //SelectedLepton
        tup_ObjectVars->Branch("pt_lepton",&pt_lepton,"pt_lepton/D");
        tup_ObjectVars->Branch("phi_lepton",&phi_lepton,"phi_lepton/D");
        tup_ObjectVars->Branch("eta_lepton",&eta_lepton,"eta_lepton/D");
        tup_ObjectVars->Branch("E_lepton",&E_lepton,"E_lepton/D");
        tup_ObjectVars->Branch("I_LepCharge",&LepCharge,"LepCharge/I");
        

        // jets
        tup_ObjectVars->Branch("I_nJets",&nJets,"nJets/I");
        tup_ObjectVars->Branch("I_nJets_CSVL",&nJets_CSVL,"nJets_CSVL/I");
        tup_ObjectVars->Branch("I_nJets_CSVM",&nJets_CSVM,"nJets_CSVM/I");
        tup_ObjectVars->Branch("I_nJets_CSVT",&nJets_CSVT,"nJets_CSVT/I");
        tup_ObjectVars->Branch("I_nJets_cMVAL",&nJets_cMVAL,"nJets_cMVAL/I");
        tup_ObjectVars->Branch("I_nJets_cMVAM",&nJets_cMVAM,"nJets_cMVAM/I");
        tup_ObjectVars->Branch("I_nJets_cMVAT",&nJets_cMVAT,"nJets_cMVAT/I");
        tup_ObjectVars->Branch("pt_jet",&pt_jet,"pt_jet[nJets]/D");
        tup_ObjectVars->Branch("phi_jet",&phi_jet,"phi_jet[nJets]/D");
        tup_ObjectVars->Branch("eta_jet",&eta_jet,"eta_jet[nJets]/D");
        tup_ObjectVars->Branch("E_jet",&E_jet,"E_jet[nJets]/D");
        tup_ObjectVars->Branch("charge_jet",&charge_jet,"charge_jet[nJets]/D");	    
        tup_ObjectVars->Branch("incl_charge_jet",&incl_charge_jet,"incl_charge_jet[nJets]/D");	    
        tup_ObjectVars->Branch("CSVv2",&CSVv2,"CSVv2[nJets]/D");
        tup_ObjectVars->Branch("cMVA",&cMVA,"cMVA[nJets]/D");
        tup_ObjectVars->Branch("cdiscCvsL_jet",&cdiscCvsL_jet,"cdiscCvsL_jet[nJets]/D");
        tup_ObjectVars->Branch("cdiscCvsB_jet",&cdiscCvsB_jet,"cdiscCvsB_jet[nJets]/D");
        tup_ObjectVars->Branch("jet_matchedMC_pdgID",&jet_matchedMC_pdgID,"jet_matchedMC_pdgID[nJets]/D");
        tup_ObjectVars->Branch("jet_matchedMC_motherpdgID",&jet_matchedMC_motherpdgID,"jet_matchedMC_motherpdgID[nJets]/D");
        tup_ObjectVars->Branch("jet_matchedMC_grannypdgID",&jet_matchedMC_grannypdgID,"jet_matchedMC_grannypdgID[nJets]/D");

        //MC variables (affected by TopPtReweighing
        tup_ObjectVars->Branch("MC_TopPt",&MC_TopPt,"MC_TopPt/D");
        tup_ObjectVars->Branch("MC_AntiTopPt",&MC_AntiTopPt,"MC_AntiTopPt/D");

       
        // met 
        tup_ObjectVars->Branch("met_Px", &met_Px, "met_Px/D"); 
        tup_ObjectVars->Branch("met_Py", &met_Py, "met_Py/D"); 
        tup_ObjectVars->Branch("met_Pt", &met_Pt, "met_Pt/D"); 
        tup_ObjectVars->Branch("met_Eta", &met_Eta,"met_Eta/D"); 
        tup_ObjectVars->Branch("met_Phi", &met_Phi, "met_Phi/D"); 

        // Advanced variables
       tup_ObjectVars->Branch("I_TOPTOPLEPHAD_JetIdx_LepTop",&TOPTOPLEPHAD_JetIdx_LepTop,"TOPTOPLEPHAD_JetIdx_LepTop/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHAD_JetIdx_HadTop",&TOPTOPLEPHAD_JetIdx_HadTop,"TOPTOPLEPHAD_JetIdx_HadTop/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHAD_JetIdx_W1",&TOPTOPLEPHAD_JetIdx_W1,"TOPTOPLEPHAD_JetIdx_W1/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHAD_JetIdx_W2",&TOPTOPLEPHAD_JetIdx_W2,"TOPTOPLEPHAD_JetIdx_W2/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHBB_JetIdx_LepTop",&TOPTOPLEPHBB_JetIdx_LepTop,"TOPTOPLEPHBB_JetIdx_LepTop/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHBB_JetIdx_HadTop",&TOPTOPLEPHBB_JetIdx_HadTop,"TOPTOPLEPHBB_JetIdx_HadTop/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHBB_JetIdx_H1",&TOPTOPLEPHBB_JetIdx_H1,"TOPTOPLEPHBB_JetIdx_H1/I");
       tup_ObjectVars->Branch("I_TOPTOPLEPHBB_JetIdx_H2",&TOPTOPLEPHBB_JetIdx_H2,"TOPTOPLEPHBB_JetIdx_H2/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_LepTop_hut",&TOPHLEPBB_JetIdx_LepTop_hut,"TOPHLEPBB_JetIdx_LepTop_hut/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_H1_hut",&TOPHLEPBB_JetIdx_H1_hut,"TOPHLEPBB_JetIdx_H1_hut/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_H2_hut",&TOPHLEPBB_JetIdx_H2_hut,"TOPHLEPBB_JetIdx_H2_hut/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_LepTop_hct",&TOPHLEPBB_JetIdx_LepTop_hct,"TOPHLEPBB_JetIdx_LepTop_hct/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_H1_hct",&TOPHLEPBB_JetIdx_H1_hct,"TOPHLEPBB_JetIdx_H1_hct/I");
       tup_ObjectVars->Branch("I_TOPHLEPBB_JetIdx_H2_hct",&TOPHLEPBB_JetIdx_H2_hct,"TOPHLEPBB_JetIdx_H2_hct/I");
       tup_ObjectVars->Branch("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD,"MVA_TOPTOPLEPHAD/D");
       tup_ObjectVars->Branch("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB,"MVA_TOPTOPLEPHBB/D");
       tup_ObjectVars->Branch("MVA_TOPHLEPBB_hut",&MVA_TOPHLEPBB_hut,"MVA_TOPHLEPBB_hut/D");
       tup_ObjectVars->Branch("MVA_TOPHLEPBB_hct",&MVA_TOPHLEPBB_hct,"MVA_TOPHLEPBB_hct/D");

        //Variables for signal/background training
	    tup_ObjectVars->Branch("HiggsMass_TOPHLEPBB_hut",&HiggsMass_TOPHLEPBB_hut,"HiggsMass_TOPHLEPBB_hut/D");
	    tup_ObjectVars->Branch("HiggsMass_TOPHLEPBB_hct",&HiggsMass_TOPHLEPBB_hct,"HiggsMass_TOPHLEPBB_hct/D");
	    tup_ObjectVars->Branch("HiggsEta_TOPHLEPBB_hut",&HiggsEta_TOPHLEPBB_hut,"HiggsEta_TOPHLEPBB_hut/D");
	    tup_ObjectVars->Branch("HiggsEta_TOPHLEPBB_hct",&HiggsEta_TOPHLEPBB_hct,"HiggsEta_TOPHLEPBB_hct/D");
	    tup_ObjectVars->Branch("TopLepMass_TOPHLEPBB_hut",&TopLepMass_TOPHLEPBB_hut,"TopLepMass_TOPHLEPBB_hut/D");
	    tup_ObjectVars->Branch("TopLepMass_TOPHLEPBB_hct",&TopLepMass_TOPHLEPBB_hut,"TopLepMass_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("TopLepPt_TOPHLEPBB_hut",&TopLepPt_TOPHLEPBB_hut,"TopLepPt_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("TopLepPt_TOPHLEPBB_hct",&TopLepPt_TOPHLEPBB_hct,"TopLepPt_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("TopLepEta_TOPHLEPBB_hut",&TopLepEta_TOPHLEPBB_hut,"TopLepEta_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("TopLepEta_TOPHLEPBB_hct",&TopLepEta_TOPHLEPBB_hct,"TopLepEta_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut,"HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct,"HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("TopLepHiggsDr_TOPHLEPBB_hut",&TopLepHiggsDr_TOPHLEPBB_hut,"TopLepHiggsDr_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("TopLepHiggsDr_TOPHLEPBB_hct",&TopLepHiggsDr_TOPHLEPBB_hct,"TopLepHiggsDr_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("HiggsBJet1CSVv2_TOPHLEPBB_hut",&HiggsBJet1CSVv2_TOPHLEPBB_hut,"HiggsBJet1CSVv2_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("HiggsBJet1CSVv2_TOPHLEPBB_hct",&HiggsBJet1CSVv2_TOPHLEPBB_hct,"HiggsBJet1CSVv2_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("HiggsBJet2CSVv2_TOPHLEPBB_hut",&HiggsBJet2CSVv2_TOPHLEPBB_hut,"HiggsBJet2CSVv2_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("HiggsBJet2CSVv2_TOPHLEPBB_hct",&HiggsBJet2CSVv2_TOPHLEPBB_hct,"HiggsBJet2CSVv2_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("TopLepBJetCSVv2_TOPHLEPBB_hut",&TopLepBJetCSVv2_TOPHLEPBB_hut,"TopLepBJetCSVv2_TOPHLEPBB_hut/D");
      tup_ObjectVars->Branch("TopLepBJetCSVv2_TOPHLEPBB_hct",&TopLepBJetCSVv2_TOPHLEPBB_hct,"TopLepBJetCSVv2_TOPHLEPBB_hct/D");
      tup_ObjectVars->Branch("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD,"TopHadMass_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD,"TopLepMass_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD,"TopLepTopHadDr_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD,"TopLepBJetCSVv2_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD,"TopHadBJetCSVv2_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD,"TopHadWNonBJet1CSVv2_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD,"TopHadWNonBJet2CSVv2_TOPTOPLEPHAD/D");
      tup_ObjectVars->Branch("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB,"HiggsMass_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB,"TopLepMass_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB,"HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB,"TopLepHiggsDr_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB,"HiggsBJet1CSVv2_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB,"HiggsBJet2CSVv2_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB,"TopLepBJetCSVv2_TOPTOPLEPHBB/D");
      tup_ObjectVars->Branch("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB,"TopHadNonBJetCSVv2_TOPTOPLEPHBB/D");

        if(debug)cout<<"created ntuples"<<endl;



        //////////////////////////////////////////////////////////////////////////////////////////////
        // Initializing TopKinFit + MVA for correct jet-comb selection
        //////////////////////////////////////////////////////////////////////////////////////////////
        TRandom3 *rnd = new TRandom3(666);

        int nToys = 50;
        std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
        std::string pdfFileName_TTHypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
        std::string pdfFileName_STHypo 
        = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

        KINFIT::kfit *kfit_SMttHypo = new KINFIT::kfit();
        KINFIT::kfit *kfit_TTHypo = new KINFIT::kfit();
        KINFIT::kfit *kfit_STHypo = new KINFIT::kfit();

        kfit_SMttHypo->Init(TOPTOPLEPHAD);
        kfit_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
        kfit_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
        kfit_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
        kfit_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
        kfit_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
        kfit_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
        kfit_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
        kfit_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
        kfit_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
        kfit_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
        kfit_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
        kfit_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
        kfit_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
        kfit_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
        kfit_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
        kfit_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
        kfit_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
        kfit_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
        kfit_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
        kfit_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
        kfit_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
        kfit_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
        kfit_SMttHypo->SetNToy(nToys);

        kfit_TTHypo->Init(TOPTOPLEPHBB);
        kfit_TTHypo->SetPDF("TopWMass",pdfFileName_TTHypo.c_str(),"TopLepWM_Fit");
        kfit_TTHypo->SetPDF("TopMass",pdfFileName_TTHypo.c_str(),"TopLepRecM_Fit");
        kfit_TTHypo->SetPDF("HiggsMass",pdfFileName_TTHypo.c_str(),"HiggsRecM_Fit");
        kfit_TTHypo->SetPDF("TopHadMass",pdfFileName_TTHypo.c_str(),"TopHadRecM_Fit");
        kfit_TTHypo->SetPDF("MetPx",pdfFileName_TTHypo.c_str(),"dMetPx_Gaus");
        kfit_TTHypo->SetPDF("MetPy",pdfFileName_TTHypo.c_str(),"dMetPy_Gaus");
        kfit_TTHypo->SetPDF("BJetPx",pdfFileName_TTHypo.c_str(),"dBJetPx_Fit");
        kfit_TTHypo->SetPDF("BJetPy",pdfFileName_TTHypo.c_str(),"dBJetPy_Fit");
        kfit_TTHypo->SetPDF("BJetPz",pdfFileName_TTHypo.c_str(),"dBJetPz_Fit");
        kfit_TTHypo->SetPDF("BJetE",pdfFileName_TTHypo.c_str(),"dBJetE_Fit");
        kfit_TTHypo->SetPDF("ElecPx",pdfFileName_TTHypo.c_str(),"dElecPx_Fit");
        kfit_TTHypo->SetPDF("ElecPy",pdfFileName_TTHypo.c_str(),"dElecPy_Fit");
        kfit_TTHypo->SetPDF("ElecPz",pdfFileName_TTHypo.c_str(),"dElecPz_Fit");
        kfit_TTHypo->SetPDF("ElecE",pdfFileName_TTHypo.c_str(),"dElecE_Fit");
        kfit_TTHypo->SetPDF("MuonPx",pdfFileName_TTHypo.c_str(),"dMuonPx_Fit");
        kfit_TTHypo->SetPDF("MuonPy",pdfFileName_TTHypo.c_str(),"dMuonPy_Fit");
        kfit_TTHypo->SetPDF("MuonPz",pdfFileName_TTHypo.c_str(),"dMuonPz_Fit");
        kfit_TTHypo->SetPDF("MuonE",pdfFileName_TTHypo.c_str(),"dMuonE_Fit");
        kfit_TTHypo->SetPDF("NonBJetPx",pdfFileName_TTHypo.c_str(),"dNonBJetPx_Fit");
        kfit_TTHypo->SetPDF("NonBJetPy",pdfFileName_TTHypo.c_str(),"dNonBJetPy_Fit");
        kfit_TTHypo->SetPDF("NonBJetPz",pdfFileName_TTHypo.c_str(),"dNonBJetPz_Fit");
        kfit_TTHypo->SetPDF("NonBJetE",pdfFileName_TTHypo.c_str(),"dNonBJetE_Fit");
        kfit_TTHypo->SetNToy(nToys);

        kfit_STHypo->Init(TOPHLEPBB);
        kfit_STHypo->SetPDF("TopWMass",pdfFileName_STHypo.c_str(),"TopLepWM_Fit");
        kfit_STHypo->SetPDF("TopMass",pdfFileName_STHypo.c_str(),"TopLepRecM_Fit");
        kfit_STHypo->SetPDF("HiggsMass",pdfFileName_STHypo.c_str(),"HiggsRecM_Fit");
        kfit_STHypo->SetPDF("MetPx",pdfFileName_STHypo.c_str(),"dMetPx_Gaus");
        kfit_STHypo->SetPDF("MetPy",pdfFileName_STHypo.c_str(),"dMetPy_Gaus");
        kfit_STHypo->SetPDF("BJetPx",pdfFileName_STHypo.c_str(),"dBJetPx_Fit");
        kfit_STHypo->SetPDF("BJetPy",pdfFileName_STHypo.c_str(),"dBJetPy_Fit");
        kfit_STHypo->SetPDF("BJetPz",pdfFileName_STHypo.c_str(),"dBJetPz_Fit");
        kfit_STHypo->SetPDF("BJetE",pdfFileName_STHypo.c_str(),"dBJetE_Fit");
        kfit_STHypo->SetPDF("ElecPx",pdfFileName_STHypo.c_str(),"dElecPx_Fit");
        kfit_STHypo->SetPDF("ElecPy",pdfFileName_STHypo.c_str(),"dElecPy_Fit");
        kfit_STHypo->SetPDF("ElecPz",pdfFileName_STHypo.c_str(),"dElecPz_Fit");
        kfit_STHypo->SetPDF("ElecE",pdfFileName_STHypo.c_str(),"dElecE_Fit");
        kfit_STHypo->SetPDF("MuonPx",pdfFileName_STHypo.c_str(),"dMuonPx_Fit");
        kfit_STHypo->SetPDF("MuonPy",pdfFileName_STHypo.c_str(),"dMuonPy_Fit");
        kfit_STHypo->SetPDF("MuonPz",pdfFileName_STHypo.c_str(),"dMuonPz_Fit");
        kfit_STHypo->SetPDF("MuonE",pdfFileName_STHypo.c_str(),"dMuonE_Fit");
        kfit_STHypo->SetNToy(nToys);
        

        //Initialize variables used in the MVAreader
        string bMethod = "CSVv2M";
        float MVAFullReco_TopHadRecM_;
        float MVAFullReco_HiggsRecM_;
        float MVAFullReco_TopLepRecM_;
        float MVAFullReco_HiggsTopLepRecDr_;
        float MVAFullReco_TopLepTopHadRecDr_;
        float MVAFullReco_TopLepRecPt_;
        float MVAPartReco_TopHadRecM_;
        float MVAPartReco_HiggsRecM_;
				float MVAPartReco_TopLepRecMT_;
        float MVAPartReco_TopLepTopHadRecDphiT_;
				float MVAPartReco_HiggsTopLepRecDphiT_;
				float MVAPartReco_TopLepRecPtT_;
        
        //Initialize readers (FullReco + PartReco) for SMttHypo
        TMVA::Reader* reader_FullReco_SMttHypo = new TMVA::Reader("!Color:!Silent");
        TMVA::Reader* reader_PartReco_SMttHypo = new TMVA::Reader("!Color:!Silent");

	      reader_FullReco_SMttHypo->AddVariable("TopHadRecM"+bMethod,&MVAFullReco_TopHadRecM_);
	      reader_FullReco_SMttHypo->AddVariable("TopLepRecM"+bMethod,&MVAFullReco_TopLepRecM_);
	      reader_FullReco_SMttHypo->AddVariable("TopLepTopHadRecDr"+bMethod,&MVAFullReco_TopLepTopHadRecDr_);
	      reader_FullReco_SMttHypo->AddVariable("TopLepRecPt"+bMethod,&MVAFullReco_TopLepRecPt_);

	      reader_PartReco_SMttHypo->AddVariable("TopHadRecM"+bMethod,&MVAPartReco_TopHadRecM_);
	      reader_PartReco_SMttHypo->AddVariable("TopLepRecMT"+bMethod,&MVAPartReco_TopLepRecMT_);
	      reader_PartReco_SMttHypo->AddVariable("TopLepTopHadRecDphiT"+bMethod,&MVAPartReco_TopLepTopHadRecDphiT_);
	      reader_PartReco_SMttHypo->AddVariable("TopLepRecPtT"+bMethod,&MVAPartReco_TopLepRecPtT_);

	      std::string weightsFile_FullReco_SMttHypo= "TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullReco"+bMethod+"_BDT.weights.xml";
	      std::string weightsFile_PartReco_SMttHypo= "TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartReco"+bMethod+"_BDT.weights.xml";
	      reader_FullReco_SMttHypo->BookMVA("BDTG method",weightsFile_FullReco_SMttHypo.c_str());
	      reader_PartReco_SMttHypo->BookMVA("BDTG method",weightsFile_PartReco_SMttHypo.c_str());
        //Initialize readers (FullReco + PartReco) for TTHypo
        TMVA::Reader* reader_FullReco_TTHypo = new TMVA::Reader("!Color:!Silent");
        TMVA::Reader* reader_PartReco_TTHypo = new TMVA::Reader("!Color:!Silent");

	      reader_FullReco_TTHypo->AddVariable("HiggsRecM"+bMethod,&MVAFullReco_HiggsRecM_);
	      reader_FullReco_TTHypo->AddVariable("TopLepRecM"+bMethod,&MVAFullReco_TopLepRecM_);
	      reader_FullReco_TTHypo->AddVariable("HiggsTopLepRecDr"+bMethod,&MVAFullReco_HiggsTopLepRecDr_);
	      reader_FullReco_TTHypo->AddVariable("TopLepRecPt"+bMethod,&MVAFullReco_TopLepRecPt_);

	      reader_PartReco_TTHypo->AddVariable("HiggsRecM"+bMethod,&MVAPartReco_HiggsRecM_);
	      reader_PartReco_TTHypo->AddVariable("TopLepRecMT"+bMethod,&MVAPartReco_TopLepRecMT_);
	      reader_PartReco_TTHypo->AddVariable("HiggsTopLepRecDphiT"+bMethod,&MVAPartReco_HiggsTopLepRecDphiT_);
	      reader_PartReco_TTHypo->AddVariable("TopLepRecPtT"+bMethod,&MVAPartReco_TopLepRecPtT_);

	      std::string weightsFile_FullReco_TTHypo= "TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullReco"+bMethod+"_BDT.weights.xml";
	      std::string weightsFile_PartReco_TTHypo= "TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartReco"+bMethod+"_BDT.weights.xml";
	      reader_FullReco_TTHypo->BookMVA("BDTG method",weightsFile_FullReco_TTHypo.c_str());
	      reader_PartReco_TTHypo->BookMVA("BDTG method",weightsFile_PartReco_TTHypo.c_str());
        //Initialize readers (FullReco + PartReco) for STHypo Hut
        TMVA::Reader* reader_FullReco_STHypo_hut = new TMVA::Reader("!Color:!Silent");
        TMVA::Reader* reader_PartReco_STHypo_hut = new TMVA::Reader("!Color:!Silent");

	      reader_FullReco_STHypo_hut->AddVariable("HiggsRecM"+bMethod,&MVAFullReco_HiggsRecM_);
	      reader_FullReco_STHypo_hut->AddVariable("TopLepRecM"+bMethod,&MVAFullReco_TopLepRecM_);
	      reader_FullReco_STHypo_hut->AddVariable("HiggsTopLepRecDr"+bMethod,&MVAFullReco_HiggsTopLepRecDr_);
	      reader_FullReco_STHypo_hut->AddVariable("TopLepRecPt"+bMethod,&MVAFullReco_TopLepRecPt_);

	      reader_PartReco_STHypo_hut->AddVariable("HiggsRecM"+bMethod,&MVAPartReco_HiggsRecM_);
	      reader_PartReco_STHypo_hut->AddVariable("TopLepRecMT"+bMethod,&MVAPartReco_TopLepRecMT_);
	      reader_PartReco_STHypo_hut->AddVariable("HiggsTopLepRecDphiT"+bMethod,&MVAPartReco_HiggsTopLepRecDphiT_);
	      reader_PartReco_STHypo_hut->AddVariable("TopLepRecPtT"+bMethod,&MVAPartReco_TopLepRecPtT_);

	      std::string weightsFile_FullReco_STHypo_hut= "TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAFullRecoHut"+bMethod+"_BDT.weights.xml";
	      std::string weightsFile_PartReco_STHypo_hut= "TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAPartRecoHut"+bMethod+"_BDT.weights.xml";
	      reader_FullReco_STHypo_hut->BookMVA("BDTG method",weightsFile_FullReco_STHypo_hut.c_str());
	      reader_PartReco_STHypo_hut->BookMVA("BDTG method",weightsFile_PartReco_STHypo_hut.c_str());
        //Initialize readers (FullReco + PartReco) for STHypo Hct
        TMVA::Reader* reader_FullReco_STHypo_hct = new TMVA::Reader("!Color:!Silent");
        TMVA::Reader* reader_PartReco_STHypo_hct = new TMVA::Reader("!Color:!Silent");

	      reader_FullReco_STHypo_hct->AddVariable("HiggsRecM"+bMethod,&MVAFullReco_HiggsRecM_);
	      reader_FullReco_STHypo_hct->AddVariable("TopLepRecM"+bMethod,&MVAFullReco_TopLepRecM_);
	      reader_FullReco_STHypo_hct->AddVariable("HiggsTopLepRecDr"+bMethod,&MVAFullReco_HiggsTopLepRecDr_);
	      reader_FullReco_STHypo_hct->AddVariable("TopLepRecPt"+bMethod,&MVAFullReco_TopLepRecPt_);

	      reader_PartReco_STHypo_hct->AddVariable("HiggsRecM"+bMethod,&MVAPartReco_HiggsRecM_);
	      reader_PartReco_STHypo_hct->AddVariable("TopLepRecMT"+bMethod,&MVAPartReco_TopLepRecMT_);
	      reader_PartReco_STHypo_hct->AddVariable("HiggsTopLepRecDphiT"+bMethod,&MVAPartReco_HiggsTopLepRecDphiT_);
	      reader_PartReco_STHypo_hct->AddVariable("TopLepRecPtT"+bMethod,&MVAPartReco_TopLepRecPtT_);

	      std::string weightsFile_FullReco_STHypo_hct= "TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAFullRecoHct"+bMethod+"_BDT.weights.xml";
	      std::string weightsFile_PartReco_STHypo_hct= "TopKinFit/test/Validation/TopHLepbb/MVA/weights/TMVAPartRecoHct"+bMethod+"_BDT.weights.xml";
	      reader_FullReco_STHypo_hct->BookMVA("BDTG method",weightsFile_FullReco_STHypo_hct.c_str());
	      reader_PartReco_STHypo_hct->BookMVA("BDTG method",weightsFile_PartReco_STHypo_hct.c_str());


        ///////////////////////////////////////////////////////////////
        // JEC
        ///////////////////////////////////////////////////////////////
        vector<JetCorrectorParameters> vCorrParam;
        JetCorrectionUncertainty *jecUnc;

        if(dName.find("Data_Run2016B")!=string::npos || dName.find("Data_Run2016C")!=string::npos || dName.find("Data_Run2016D")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016BCDV2_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016BCDV2_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016BCDV2_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016BCDV2_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016BCDV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("Data_Run2016E")!=string::npos || dName.find("Data_Run2016F")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016EFV2_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016EFV2_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016EFV2_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016EFV2_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016EFV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("Data_Run2016G")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016GV2_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016GV2_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016GV2_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016GV2_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016GV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("Data_Run2016H")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016HV2_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016HV2_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016HV2_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016HV2_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_23Sep2016V2/Spring16_23Sep2016HV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Spring16_25nsV10/Spring16_25nsV10_MC_Uncertainty_AK4PFchs.txt");
        }

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
//                cin.get();
                cout << " " << endl;
                cout << "------------NEW EVENT: " << ievt << " --------------" << endl;
            }

            double ievt_d = ievt;

            if(ievt%10000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-event_start)/(end_d-event_start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            W_puSF =1;
            W_fleptonSF =1;
            W_fleptonSF_Plus =1;
            W_fleptonSF_Minus =1;
            W_btagWeight_CSVv2M_mujets_central =1;
            W_btagWeight_CSVv2M_mujets_up =1;
            W_btagWeight_CSVv2M_mujets_down =1;
            W_btagWeight_shape =1;
            W_btagWeight_shape_up_lf =1; 
            W_btagWeight_shape_down_lf =1; 
            W_btagWeight_shape_up_hf =1; 
            W_btagWeight_shape_down_hf =1; 
            W_btagWeight_shape_up_hfstats1 =1; 
            W_btagWeight_shape_down_hfstats1 =1; 
            W_btagWeight_shape_up_hfstats2 =1; 
            W_btagWeight_shape_down_hfstats2 =1; 
            W_btagWeight_shape_up_lfstats1 =1; 
            W_btagWeight_shape_down_lfstats1 =1; 
            W_btagWeight_shape_up_lfstats2 =1; 
            W_btagWeight_shape_down_lfstats2 =1; 
            W_btagWeight_shape_up_cferr1 =1; 
            W_btagWeight_shape_down_cferr1 =1; 
            W_btagWeight_shape_up_cferr2 =1; 
            W_btagWeight_shape_down_cferr2 =1; 
            W_nloWeight =1;// for amc@nlo samples
            W_weight1 =1;
            W_weight2 =1;
            W_weight3 =1;
            W_weight4 =1;
            W_weight5 =1;
            W_weight6 =1;
            W_weight7 =1;
            W_weight8 =1; 
            W_MuonIDSF =1; //One of the 3 components for the total muon SF
            W_MuonIsoSF =1; //One of the 3 components for the total muon SF
            W_MuonTrigSF =1;//One of the 3 components for the total muon SF
            W_MuonTrigSF = 1;//Used in calculation for W_MuonTrigSF
            W_ElectronIDSF =1; //One of the 2 components for the total electron SF
            W_ElectronRecoSF =1;     //One of the 2 components for the total electron SF
            W_MuonIDSF_Plus =1; //One of the 3 components for the total muon SF
            W_MuonIsoSF_Plus =1; //One of the 3 components for the total muon SF
            W_MuonTrigSF_Plus =1;//One of the 3 components for the total muon SF
            W_MuonTrigSF_Plus = 1;//Used in calculation for W_MuonTrigSF
            W_ElectronIDSF_Plus =1; //One of the 2 components for the total electron SF
            W_ElectronRecoSF_Plus =1;     //One of the 2 components for the total electron SF
            W_MuonIDSF_Minus =1; //One of the 3 components for the total muon SF
            W_MuonIsoSF_Minus =1; //One of the 3 components for the total muon SF
            W_MuonTrigSF_Minus =1;//One of the 3 components for the total muon SF
            W_MuonTrigSF_Minus = 1;//Used in calculation for W_MuonTrigSF
            W_ElectronIDSF_Minus =1; //One of the 2 components for the total electron SF
            W_ElectronRecoSF_Minus =1;     //One of the 2 components for the total electron SF
            W_TopPtReweighing =1;
            genTTX = -666;

            MC_TopPt = -99.;
            MC_AntiTopPt = -99.;
            
	          //JetIndices_correctJetComb
	          TOPTOPLEPHAD_JetIdx_LepTop = -99;
	          TOPTOPLEPHAD_JetIdx_HadTop = -99;
	          TOPTOPLEPHAD_JetIdx_W1 = -99;
	          TOPTOPLEPHAD_JetIdx_W2 = -99;
	          TOPTOPLEPHBB_JetIdx_LepTop = -99;
	          TOPTOPLEPHBB_JetIdx_HadTop = -99;
	          TOPTOPLEPHBB_JetIdx_H1 = -99;
	          TOPTOPLEPHBB_JetIdx_H2 = -99;
	          TOPHLEPBB_JetIdx_LepTop_hut = -99;
	          TOPHLEPBB_JetIdx_H1_hut = -99;
	          TOPHLEPBB_JetIdx_H2_hut = -99;
	          TOPHLEPBB_JetIdx_LepTop_hct = -99;
	          TOPHLEPBB_JetIdx_H1_hct = -99;
	          TOPHLEPBB_JetIdx_H2_hct = -99;
            MVA_TOPTOPLEPHAD = -999.;
            MVA_TOPTOPLEPHBB = -999.;
            MVA_TOPHLEPBB_hut = -999.;
            MVA_TOPHLEPBB_hct = -999.;

            HiggsMass_TOPHLEPBB_hut = -999.;
            HiggsMass_TOPHLEPBB_hct = -999.;
            HiggsEta_TOPHLEPBB_hut = -999.;
            HiggsEta_TOPHLEPBB_hct = -999.;
            TopLepMass_TOPHLEPBB_hut = -999.;
            TopLepMass_TOPHLEPBB_hct = -999.;
            TopLepPt_TOPHLEPBB_hut = -999.;
            TopLepPt_TOPHLEPBB_hct = -999.;
            TopLepEta_TOPHLEPBB_hut = -999.;
            TopLepEta_TOPHLEPBB_hct = -999.;
            HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut = -999.;
            HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct = -999.;
            TopLepHiggsDr_TOPHLEPBB_hut = -999.;
            TopLepHiggsDr_TOPHLEPBB_hct = -999.;
            HiggsBJet1CSVv2_TOPHLEPBB_hut = -999.;
            HiggsBJet1CSVv2_TOPHLEPBB_hct = -999.;
            HiggsBJet2CSVv2_TOPHLEPBB_hut = -999.;
            HiggsBJet2CSVv2_TOPHLEPBB_hct = -999.;
            TopLepBJetCSVv2_TOPHLEPBB_hut = -999.;
            TopLepBJetCSVv2_TOPHLEPBB_hct = -999.;
            TopHadMass_TOPTOPLEPHAD = -999.;
            TopLepMass_TOPTOPLEPHAD = -999.;
            TopLepTopHadDr_TOPTOPLEPHAD = -999.;
            TopLepBJetCSVv2_TOPTOPLEPHAD = -999.;
            TopHadBJetCSVv2_TOPTOPLEPHAD = -999.;
            TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = -999.;
            TopHadWNonBJet2CSVv2_TOPTOPLEPHAD = -999.;
            HiggsMass_TOPTOPLEPHBB = -999.;
            TopLepMass_TOPTOPLEPHBB = -999.;
            HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = -999.;
            TopLepHiggsDr_TOPTOPLEPHBB = -999.;
            HiggsBJet1CSVv2_TOPTOPLEPHBB = -999.;
            HiggsBJet2CSVv2_TOPTOPLEPHBB = -999.;
            TopLepBJetCSVv2_TOPTOPLEPHBB = -999.;
            TopHadNonBJetCSVv2_TOPTOPLEPHBB = -999.;

            if(debug)cout<<"before tree load"<<endl;
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            if(debug)cout<<"after tree load"<<endl;

            //////////////////////////////////////////////////////////////////
            ///  Include trigger set up here when using data
            //////////////////////////////////////////////////////////////////
            datasets[d]->eventTree()->LoadTree(ievt);
            if(datasets[d]->eventTree()->GetBranch("genTTX_id_")) genTTX = event->getgenTTX_id(); //Safety check to make sure the code doesn't crash when using trees that do not have the genTTX_id variable

		        if(dName.find("TTJets")!=string::npos)//Split ttbar samples into ttbb, ttcc and ttlf
		        {
                bool isttbb = (genTTX == 051 || genTTX == 151 || genTTX == 251 ||
		              genTTX == 052 || genTTX == 152 || genTTX == 252 ||
		              genTTX == 053 || genTTX == 153 || genTTX == 253 ||
		              genTTX == 054 || genTTX == 154 || genTTX == 254 ||
		              genTTX == 055 || genTTX == 155 || genTTX == 255);
               
                bool isttcc = (genTTX == 041 || genTTX == 141 || genTTX == 241 ||
		              genTTX == 042 || genTTX == 142 || genTTX == 242 ||
		              genTTX == 043 || genTTX == 143 || genTTX == 243 ||
		              genTTX == 044 || genTTX == 144 || genTTX == 244 ||
		              genTTX == 045 || genTTX == 145 || genTTX == 245);
               
                bool isttlf = (!isttbb && !isttcc);
            }

	          run_num = event->runId(); 
	          evt_num = event->eventId();
	          lumi_num=event->lumiBlockId(); 
	          nvtx = vertex.size();
	          npu = (int) event->nTruePU(); 
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
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "minus", false);
                        else if(doJERShift == 2)
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "plus", false);
                        else
                            jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
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
                    
                    jetTools->correctMETTypeOne(init_jets, mets[0], isData);


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
          
          	    W_nloWeight = mc_baseweight;
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
            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets, rho);

            // Define object selection cuts

            int nMu, nEl, nLooseMu, nLooseEl; //number of (loose) muons/electrons
            nMu = 0, nEl = 0, nLooseMu=0, nLooseEl=0;

            /////////////////////////////////////////////
            // Define object selection cuts: https://twiki.cern.ch/twiki/bin/view/CMS/TTbarXSecSynchronization#Lepton_jets
            /////////////////////////////////////////////
            if (Muon)
            {
				        if (debug)cout<<"Getting Jets"<<endl;
				        selectedOrigJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Loose"); // ApplyJetId
				        if (debug)cout<<"Getting Tight Muons"<<endl;
				        selectedMuons                                       = r2selection.GetSelectedMuons(27,2.1,0.15, "Tight", "Spring15"); //Selected - Trigger: HLT_Iso(Tk)Mu24_v*
				        if (debug)cout<<"Getting Loose Electrons"<<endl;
				        selectedElectrons                                   = r2selection.GetSelectedElectrons(10,2.5,"Loose", "Spring16_80X", true, true); //Vetoed  
				        if (debug)cout<<"Getting Loose Muons"<<endl;
				        selectedExtraMuons                                  = r2selection.GetSelectedMuons(10, 2.4, 0.25,"Loose","Spring15"); //Vetoed         
            }
            else if (Electron)
            {
				        if (debug)cout<<"Getting Jets"<<endl;
				        selectedOrigJets                                        = r2selection.GetSelectedJets(30,2.4,true,"Loose"); // ApplyJetId
				        if (debug)cout<<"Getting Loose Muons"<<endl;
				        selectedMuons                                       = r2selection.GetSelectedMuons(10, 2.4, 0.25,"Loose","Spring15"); //Vetoed
				        if (debug)cout<<"Getting Electrons"<<endl;
				        selectedElectrons                                   = r2selection.GetSelectedElectrons(35,2.1,"Medium", "Spring16_80X", true, true); //Selected - Trigger:  HLT_Ele32_eta2p1_WPTight_Gsf_v*                
				        if (debug)cout<<"Getting Loose Electrons"<<endl;
				        selectedExtraElectrons                              = r2selection.GetSelectedElectrons(10,2.5,"Loose", "Spring16_80X", true, true); //Vetoed
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
                W_puSF=1;
                W_puSF_Minus=1;
                W_puSF_Plus=1;
            }
            else
            {
                W_puSF = W_puSFs.ITweight( (int)event->nTruePU() );
                W_puSF_Minus = W_puSFs_Minus.ITweight( (int)event->nTruePU() );
                W_puSF_Plus = W_puSFs_Plus.ITweight( (int)event->nTruePU() );
            }
            if(debug) cout << "W_puSF: " << W_puSF << endl;
            /////////////////////////////////////////////////
            //                   Lepton SF                 //
            /////////////////////////////////////////////////
            float lum_RunsBCDEF = 15.658183109;// /fb
            float lum_RunsGH = 15.199167277;// /fb
            if(bLeptonSF && !isData)
            {
                if(Muon && nMu>0){
                    W_MuonIDSF = (muonSFWeightID_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonIsoSF = (muonSFWeightIso_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonTrigSF = (muonSFWeightTrig_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightTrig_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);

                    W_MuonIDSF_Plus = (muonSFWeightID_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonIsoSF_Plus = (muonSFWeightIso_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonTrigSF_Plus = (muonSFWeightTrig_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightTrig_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);

                    W_MuonIDSF_Minus = (muonSFWeightID_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonIsoSF_Minus = (muonSFWeightIso_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    W_MuonTrigSF_Minus = (muonSFWeightTrig_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightTrig_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                    
                    W_fleptonSF_Plus = W_MuonIDSF_Plus * W_MuonIsoSF_Plus * W_MuonTrigSF_Plus;
                    W_fleptonSF_Minus = W_MuonIDSF_Minus * W_MuonIsoSF_Minus * W_MuonTrigSF_Minus;
                    W_fleptonSF = W_MuonIDSF * W_MuonIsoSF * W_MuonTrigSF;
                }
                else if(Electron && nEl>0){
                    W_ElectronIDSF = electronSFWeightID->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                    W_ElectronRecoSF = electronSFWeightReco->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                    W_fleptonSF = W_ElectronIDSF * W_ElectronRecoSF;
                    W_ElectronIDSF_Plus = electronSFWeightID->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),1);
                    W_ElectronRecoSF_Plus = electronSFWeightReco->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),1);
                    W_fleptonSF_Plus = W_ElectronIDSF_Plus * W_ElectronRecoSF_Plus;
                    W_ElectronIDSF_Minus = electronSFWeightID->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),-1);
                    W_ElectronRecoSF_Minus = electronSFWeightReco->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),-1);
                    W_fleptonSF_Minus = W_ElectronIDSF_Minus * W_ElectronRecoSF_Minus;
                }
            }


            if(debug) cout<<"lepton SF:  "<<W_fleptonSF  << endl;

            /////////////////////////////////////////////////
            //                   Btag SF                    //
            /////////////////////////////////////////////////
            if(!bTagReweight_FillMChistos && CSVv2nonshape)
            {
                if(dName.find("Data")==string::npos) //If sample is data, no b-tag-reweighting
                {
                    if(debug) cout << "Applying b-tag weights " << endl;
                    W_btagWeight_CSVv2M_mujets_central =  btwt_CSVv2M_mujets_central->getMCEventWeight(selectedOrigJets, false);
                    //                    W_btagWeight_CSVv2M_mujets_up =  btwt_CSVv2M_mujets_up->getMCEventWeight(selectedOrigJets, false);
                    //                    W_btagWeight_CSVv2M_mujets_down =  btwt_CSVv2M_mujets_down->getMCEventWeight(selectedOrigJets, false);
                    //btagWeight_cMVAM_mujets_central =  btwt_cMVAM_mujets_central->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_cMVAM_mujets_up =  btwt_CSVv2M_mujets_up->getMCEventWeight(selectedOrigJets, false);
                    //                    btagWeight_cMVAM_mujets_down =  btwt_CSVv2M_mujets_down->getMCEventWeight(selectedOrigJets, false);
                }
            }
            ////////////////////////////////////////////////
            // Pre-baseline initializations
            ////////////////////////////////////////////////           
            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV <<endl;

            if(debug) cout << "Past cut 0: NO CUTS" << endl;
            cutstep[0]=cutstep[0]+scaleFactor; //Order of appearance of cutstep & nCuts is important here

            eventCount++;
            //////////////////////////////////////////////////////
            // Applying baseline lepton selection
            //////////////////////////////////////////////////////


            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            if(debug) cout << "Past cut 1: good PV selection" << endl;
            cutstep[1]=cutstep[1]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
            passed_Step2++;

            bool trigged = false;
            if (applyTriggers)
            {
              trigger->checkAvail(run_num, datasets, d, &treeLoader, event, printTriggers);
              trigged = trigger->checkIfFired();
              
            }
            else trigged = true;
//             if(dName.find("NP")!=string::npos)
//            {
//                trigged = true;
//            }

            //Event cleaning filters
            if(isData)
            {
                if(!event->getHBHENoiseFilter() || !event->getHBHENoiseIsoFilter() || !event->getEEBadScFilter() || !event->getglobalTightHalo2016Filter()
                 || !event->getEcalDeadCellTriggerPrimitiveFilter() || !event->getPVFilter() || !event->getBadChCandFilter() || !event->getBadPFMuonFilter()) continue;
            }
            else if(!isData)
            {
                 if(!event->getHBHENoiseFilter() || !event->getHBHENoiseIsoFilter() || !event->getglobalTightHalo2016Filter()
                 || !event->getEcalDeadCellTriggerPrimitiveFilter() || !event->getPVFilter() || !event->getBadChCandFilter() || !event->getBadPFMuonFilter()) continue;
            }
            cutstep[2]=cutstep[2]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
                
            passed_Step3++;


            if(!trigged) continue;
            if(debug) cout << "Past cut 2: passed_FinalSelection trigger" << endl;
            cutstep[3]=cutstep[3]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
            passed_Step4++;

            if (debug)
            {
              	cout <<" applying baseline event selection..."<<endl;
              	cout <<"number of muons: " << nMu <<endl;
              	cout <<"number of electrons: " << nEl <<endl;
            }
            //Apply the lepton, btag and HT selections
            if (Muon && !Electron)
            {
                if  (  !( nMu ==1)) continue; // Muon Channel Selection
                if(debug) cout << "Past cut 3: Single Muon selected" << endl;
            }
            else if (!Muon && Electron)
            {
                if  (  !( nEl ==1)) continue; // Electron Channel Selection
                if(debug) cout << "Past cut 3: Single Electron selected" << endl;
            }
            else
            {
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            cutstep[4]=cutstep[4]+scaleFactor; //Order of appearance of cutstep & nCuts is important here

            passed_Step5++;

			      if(Muon && !Electron)
			      {
                  if( !(nEl == 0)) continue;
                  cutstep[5]=cutstep[5]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
                  passed_Step6++;
				          if(nLooseMu != 1) continue;
                  cutstep[6]=cutstep[6]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
                  passed_Step7++;
	                if (debug)	cout <<"Vetoed extra muons..."<<endl;
			      }
			      if(!Muon && Electron)
			      {
                  if( !(nMu == 0)) continue;
                  cutstep[5]=cutstep[5]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
                  passed_Step6++;
				          if(nLooseEl != 1) continue;
                  cutstep[6]=cutstep[6]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
                  passed_Step7++;
	                if (debug)	cout <<"Vetoed extra electrons..."<<endl;
			      }
            if(debug) cout << "Past cut 5: Vetoed extra loose leptons" << endl;
            cutstep[4]=cutstep[4]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
            
            if(Electron)
            {
              pt_lepton=selectedElectrons[0]->Pt();
	            phi_lepton=selectedElectrons[0]->Phi();
	            eta_lepton=selectedElectrons[0]->Eta();
	            E_lepton=selectedElectrons[0]->E();
	            LepCharge=selectedElectrons[0]->charge();
	            eta_superCluster_electron=selectedElectrons[0]->superClusterEta();
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
                pt_lepton=selectedMuons[0]->Pt();
	              phi_lepton=selectedMuons[0]->Phi();
	              eta_lepton=selectedMuons[0]->Eta();
	              E_lepton=selectedMuons[0]->E();
                LepCharge=selectedMuons[0]->charge();
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

            for(int jetbtag = 0; jetbtag<selectedJets.size(); jetbtag++)
  			    {

                if(!isData)
                {
                    float bTagEff = 1, bTagEff_LFUp = 1, bTagEff_LFDown = 1, bTagEff_HFUp = 1, bTagEff_HFDown = 1, bTagEff_HFStats1Up = 1,
                    bTagEff_HFStats1Down = 1, bTagEff_HFStats2Up = 1, bTagEff_HFStats2Down = 1, bTagEff_LFStats1Up = 1, bTagEff_LFStats1Down = 1,
                    bTagEff_LFStats2Up = 1, bTagEff_LFStats2Down = 1, bTagEff_CFErr1Up = 1, bTagEff_CFErr1Down = 1, bTagEff_CFErr2Up = 1, bTagEff_CFErr2Down = 1;

                    float jetpt = selectedJets[jetbtag]->Pt();
                    float jeteta = selectedJets[jetbtag]->Eta();
                    float jetdisc = selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                    bool isBFlav = false;
                    bool isLFlav = false;
                    bool isCFlav = false;
                    if(jetdisc<0.0) jetdisc = -0.05;
                    if(jetdisc>1.0) jetdisc = 1.0;
                    BTagEntry::JetFlavor jflav;
                    int jethadronflav = std::abs(selectedJets[jetbtag]->hadronFlavour());
                    if(debug) cout<<"hadron flavour: "<<jethadronflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
                    if(jethadronflav == 5){
                        jflav = BTagEntry::FLAV_B;
                        isBFlav =true;
                    }
                    else if(jethadronflav == 4){
                        jflav = BTagEntry::FLAV_C;
                        isCFlav=true;
                    }
                    else{
                        jflav = BTagEntry::FLAV_UDSG;
                        isLFlav = true;
                    }
                    if( doJESShift == 2)        bTagEff = reader_JESUp->eval(jflav, jeteta, jetpt, jetdisc);
                    else if( doJESShift == 1) bTagEff = reader_JESDown->eval(jflav, jeteta, jetpt, jetdisc);
                    else bTagEff = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);

                    if( isBFlav ) bTagEff_LFUp = reader_LFUp->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isBFlav ) bTagEff_LFDown = reader_LFDown->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_HFUp = reader_HFUp->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_HFDown = reader_HFDown->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isBFlav ) bTagEff_HFStats1Up = reader_HFStats1Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isBFlav ) bTagEff_HFStats1Down = reader_HFStats1Down->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isBFlav ) bTagEff_HFStats2Up = reader_HFStats2Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isBFlav ) bTagEff_HFStats2Down = reader_HFStats2Down->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_LFStats1Up = reader_LFStats1Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_LFStats1Down = reader_LFStats1Down->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_LFStats2Up = reader_LFStats2Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isLFlav ) bTagEff_LFStats2Down = reader_LFStats2Down->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isCFlav ) bTagEff_CFErr1Up = reader_CFErr1Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isCFlav ) bTagEff_CFErr1Down = reader_CFErr1Down->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isCFlav ) bTagEff_CFErr2Up = reader_CFErr2Up->eval(jflav, jeteta, jetpt, jetdisc);
                    if( isCFlav ) bTagEff_CFErr2Down = reader_CFErr2Down->eval(jflav, jeteta, jetpt, jetdisc);
                    
                    //If jet is not the appropriate flavor for that systematic, use the nominal reader so that all weights will be on the same
                    //jet multiplicity footing.
                    if( !isBFlav ) bTagEff_LFUp = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_LFDown = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_HFUp = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_HFDown = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats1Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats1Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats2Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats2Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats1Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats1Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats2Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats2Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr1Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr1Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr2Up = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr2Down = bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc);

                    W_btagWeight_shape_up_lf *= bTagEff_LFUp;
                    W_btagWeight_shape_down_lf *= bTagEff_LFDown;
                    W_btagWeight_shape_up_hf *= bTagEff_HFUp;
                    W_btagWeight_shape_down_hf *= bTagEff_HFDown;
                    W_btagWeight_shape_up_hfstats1 *= bTagEff_HFStats1Up;
                    W_btagWeight_shape_down_hfstats1 *= bTagEff_HFStats1Down;
                    W_btagWeight_shape_up_hfstats2 *= bTagEff_HFStats2Up;
                    W_btagWeight_shape_down_hfstats2 *= bTagEff_HFStats2Down;            
                    W_btagWeight_shape_up_lfstats1 *= bTagEff_LFStats1Up;
                    W_btagWeight_shape_down_lfstats1 *= bTagEff_LFStats1Down;
                    W_btagWeight_shape_up_lfstats2 *= bTagEff_LFStats2Up;
                    W_btagWeight_shape_down_lfstats2 *= bTagEff_LFStats2Down; 
                    W_btagWeight_shape_up_cferr1 *= bTagEff_CFErr1Up;
                    W_btagWeight_shape_down_cferr1 *= bTagEff_CFErr1Down;
                    W_btagWeight_shape_up_cferr2 *= bTagEff_CFErr2Up;
                    W_btagWeight_shape_down_cferr2 *= bTagEff_CFErr2Down; 

                    W_btagWeight_shape*=bTagEff;
                    if (debug){
                        cout<<"hadron flavour: "<<jethadronflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
                        cout << " isBFlav " << isBFlav;
                        cout << " isLFlav" << isLFlav;
                        cout << " isCFlav" << isCFlav << endl;
                        cout << " reader_csvv2->eval(jflav, jeteta, jetpt, jetdisc);  " << bTagReader_shape->eval(jflav, jeteta, jetpt, jetdisc) << endl;
                        cout << " bTagEff " << bTagEff;
                        cout << " bTagEff_LFUp " << bTagEff_LFUp;
                        cout << " bTagEff_LFDown " << bTagEff_LFDown << endl;
                        cout << " bTagEff_HFUp " << bTagEff_HFUp;
                        cout << " bTagEff_HFDown " << bTagEff_HFDown << endl;
                        cout << " bTagEff_HFStats1Up " << bTagEff_HFStats1Up;
                        cout << " bTagEff_HFStats1Down " << bTagEff_HFStats1Down << endl;
                        cout << " bTagEff_HFStats2Up " << bTagEff_HFStats2Up;
                        cout << " bTagEff_HFStats2Down " << bTagEff_HFStats2Down << endl;
                        cout << " bTagEff_LFStats1Up " << bTagEff_LFStats1Up;
                        cout << " bTagEff_LFStats1Down " << bTagEff_LFStats1Down << endl;
                        cout << " bTagEff_LFStats2Up " << bTagEff_LFStats2Up;
                        cout << " bTagEff_LFStats2Down " << bTagEff_LFStats2Down << endl;
                        cout << " bTagEff_CFErr1Up " << bTagEff_CFErr1Up;
                        cout << " bTagEff_CFErr1Down " << bTagEff_CFErr1Down << endl;
                        cout << " bTagEff_CFErr2Up " << bTagEff_CFErr2Up;
                        cout << " bTagEff_CFErr2Down " << bTagEff_CFErr2Down << endl;
                    }
                    if(debug)cout<<"btag efficiency = "<<bTagEff<<endl;       

                }
            }
            
            float btagWeight = 1;
            if(btagger == "CSVv2M") btagWeight = W_btagWeight_CSVv2M_mujets_central;
//            else if(btagger == "cMVAM") btagWeight = btagWeight_cMVAM_mujets_central;
            if(!CSVv2nonshape)  btagWeight = W_btagWeight_shape;


            if(debug) cout<<"btag SF:  "<< btagWeight << endl;
            scaleFactor = scaleFactor * W_puSF * W_fleptonSF * btagWeight;
            if(isData) scaleFactor = 1;

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
            for(int jetbtag = 0; jetbtag<selectedJets.size(); jetbtag++)
  			    {
            
		           	if (selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Loose   )
		            {
		              	  selectedLBJets.push_back(selectedJets[jetbtag]);
		                  if (selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Medium)
		                  {
		                  	  selectedMBJets.push_back(selectedJets[jetbtag]);
		                      if (selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Tight)
		                      {
		                      	selectedTBJets.push_back(selectedJets[jetbtag]);
		                      }
		         				      else selectedLightJets_TWP.push_back(selectedJets[jetbtag]);
		              	  }
		         			    else
		         			    {
		         				      selectedLightJets_MWP.push_back(selectedJets[jetbtag]);
		         				      selectedLightJets_TWP.push_back(selectedJets[jetbtag]);
		         			    }
		         		}
		         		else
		         		{
		         			selectedLightJets_LWP.push_back(selectedJets[jetbtag]);
		         			selectedLightJets_MWP.push_back(selectedJets[jetbtag]);
		         			selectedLightJets_TWP.push_back(selectedJets[jetbtag]);
		         		}
		        }

			      //////////////////////////////////////
			      // Jet variables //
			      /////////////////////////////////////
			      nJets = 0; 
            nJets_CSVT =  0; 
	          nJets_CSVM =  0;
            nJets_CSVL =  0;
            for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
            {
                      
                pt_jet[nJets]=selectedJets[seljet]->Pt(); 
                phi_jet[nJets]=selectedJets[seljet]->Phi();
                eta_jet[nJets]=selectedJets[seljet]->Eta();
                E_jet[nJets]=selectedJets[seljet]->E();
                charge_jet[nJets]=selectedJets[seljet]->charge();
                incl_charge_jet[nJets]=selectedJets[seljet]->inclusiveJetCharge();
                CSVv2[nJets]=selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() ;
                cMVA[nJets]=selectedJets[seljet]->btag_PFCombinedMVAV2BJetTags() ;
                cdiscCvsB_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsBJetTags() ;
                cdiscCvsL_jet[nJets]=selectedJets[seljet]->ctag_pfCombinedCvsLJetTags() ;
                if(cMVA[nJets] > cMVA_workingpointvalue_Loose) nJets_cMVAL++;
                if(cMVA[nJets] > cMVA_workingpointvalue_Medium) nJets_cMVAM++;
                if(cMVA[nJets] > cMVA_workingpointvalue_Tight) nJets_cMVAT++;
                if(CSVv2[nJets] > CSVv2_workingpointvalue_Loose) nJets_CSVL++;
                if(CSVv2[nJets] > CSVv2_workingpointvalue_Medium) nJets_CSVM++;
                if(CSVv2[nJets] > CSVv2_workingpointvalue_Tight) nJets_CSVT++;
                nJets++;
            }
	          double met_px = mets[0]->Px();
	          double met_py = mets[0]->Py();
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
	          met_Px = mets[0]->Px(); 
	          met_Py = mets[0]->Py(); 
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
			      if(selectedJets.size() < 3)  continue;
            if(debug) cout << "Past cut 5: passed_FinalSelection number of jets cut" << endl;
            cutstep[7]=cutstep[7]+scaleFactor; //Order of appearance of cutstep & nCuts is important here
            passed_Step8++;

            ///////////////////////////////////////////////////
            // Fill b-tag histos for scale factors
            // info: http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees_v2.pdf
            ////////////////////////////////////////////////////
            if(bTagReweight_FillMChistos && CSVv2nonshape)
            {
                if(dName.find("Data")==string::npos)        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
                {
                    btwt_CSVv2M_mujets_central->FillMCEfficiencyHistos(selectedJets);
//                    btwt_CSVv2M_mujets_up->FillMCEfficiencyHistos(selectedJets);
//                    btwt_CSVv2M_mujets_down->FillMCEfficiencyHistos(selectedJets);
//                    btwt_cMVAM_mujets_central->FillMCEfficiencyHistos(selectedJets);
//                    btwt_cMVAM_mujets_up->FillMCEfficiencyHistos(selectedJets);
//                    btwt_cMVAM_mujets_down->FillMCEfficiencyHistos(selectedJets);
                }
                continue;
            }

		  	    if(selectedMBJets.size() < 1) continue;
	          if (debug)	cout <<"Cut on nb b-jets..."<<endl;
            if(debug) cout << "Past cut 7: passed_FinalSelection cut on number of b-jets" << endl;
            cutstep[8]=cutstep[8]+scaleFactor; //Order of appearance of cutstep & nCuts is important here


            if(debug)
            {
                cout<<"Selection passed_FinalSelection."<<endl;
            }
            passed_FinalSelection++;

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
            pair<unsigned int, unsigned int> hadronicTopJet = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> HiggsBJet1_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> HiggsBJet2_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet1_ = pair<unsigned int,unsigned int>(9999,9999);
            pair<unsigned int, unsigned int> hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
            
            int HiggsBJetCounter = 0;                  
            int pdgID_top = 6; //top quark
                  
            bool Posleptonmatched = false;
            bool Negleptonmatched = false;
            
            double TopPtReweigh = 1;
            double AntiTopPtReweigh = 1;
            
            if(dName != "data" && dName != "Data" && dName != "Data" && dName != "D_ata")
            {
                treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
                sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
                  
                mcParticlesMatching_.clear();
                    
                    
                for (unsigned int i = 0; i < mcParticles.size(); i++)
                {
                    if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process                   

                    //Calculating event weight according to the TopPtReweighing: https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
                    if(mcParticles[i]->type() == 6)
                    {
                        if(mcParticles[i]->Pt() < 400) TopPtReweigh  = TMath::Exp(0.0615-0.0005*mcParticles[i]->Pt());
                        else TopPtReweigh  = TMath::Exp(0.0615-0.0005*400);
                        MC_TopPt = mcParticles[i]->Pt();
                    }
                    else if(mcParticles[i]->type() == -6)
                    {
                        if(mcParticles[i]->Pt() < 400) AntiTopPtReweigh  = TMath::Exp(0.159-0.00141*mcParticles[i]->Pt());
                        else AntiTopPtReweigh  = TMath::Exp(0.159-0.00141*400);
                        MC_AntiTopPt = mcParticles[i]->Pt();
                    }
                      
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
            	      if(debug) cout << "matched MCparticle of jet " << JetPartonPair[i].first << " has pdgID: " << jet_matchedMC_pdgID[JetPartonPair[i].first] << endl;
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
                    if ( fabs(mcParticlesMatching_[j]->type()) == 5 || fabs(mcParticlesMatching_[j]->type()) == 4 || fabs(mcParticlesMatching_[j]->type()) == 2 )
                    {
                        if ( ( Posleptonmatched && mcParticlesMatching_[j]->motherType() == -pdgID_top )
                            || ( Negleptonmatched && mcParticlesMatching_[j]->motherType() == pdgID_top ) )  // if mu+ (top decay leptonic) and mother is antitop ---> hadronic b
                        {
                            hadronicTopJet = JetPartonPair[i];
                        }
                        else if ( ( Posleptonmatched && mcParticlesMatching_[j]->motherType() == pdgID_top )
                                        || ( Negleptonmatched && mcParticlesMatching_[j]->motherType() == -pdgID_top ) )
                        {
                            leptonicBJet_ = JetPartonPair[i];
                        }
                        
                        if(mcParticlesMatching_[j]->motherType() == 25 && HiggsBJetCounter < 1)
                        {
                            HiggsBJet1_ = JetPartonPair[i];
                            HiggsBJetCounter++;
                        }
                        else if(mcParticlesMatching_[j]->motherType() == 25 && HiggsBJetCounter == 1)
                        {
                            HiggsBJet2_ = JetPartonPair[i];
                            HiggsBJetCounter++;
                        }

                    }
                }  /// End loop over Jet Parton Pairs
            
            }// End MC matching (end of data-if-statement)


            W_TopPtReweighing = TMath::Sqrt(TopPtReweigh*AntiTopPtReweigh);

            ////////////////////////////////////////////////////////////////////////////////
            // Jet-combination selection according to MVA-TopKinFit
            ////////////////////////////////////////////////////////////////////////////////
            if(applyMVAJetComb)
            {
                sort(selectedJets.begin(), selectedJets.end(), HighestCSVBtag());//Sort according to the highest b-tag value.
                
                //First define all the variables that go into the kinfit procedure
	              std::vector<float> BJetPt_SMttHypo;
	              std::vector<float> BJetEta_SMttHypo;
	              std::vector<float> BJetPhi_SMttHypo;
	              std::vector<float> BJetE_SMttHypo;
	              std::vector<float> BJetPt_TTHypo;
	              std::vector<float> BJetEta_TTHypo;
	              std::vector<float> BJetPhi_TTHypo;
	              std::vector<float> BJetE_TTHypo;
	              std::vector<float> BJetPt_STHypo;
	              std::vector<float> BJetEta_STHypo;
	              std::vector<float> BJetPhi_STHypo;
	              std::vector<float> BJetE_STHypo;

	              std::vector<float> NonBJetPt_SMttHypo;
	              std::vector<float> NonBJetEta_SMttHypo;
	              std::vector<float> NonBJetPhi_SMttHypo;
	              std::vector<float> NonBJetE_SMttHypo;
	              std::vector<float> NonBJetPt_TTHypo;
	              std::vector<float> NonBJetEta_TTHypo;
	              std::vector<float> NonBJetPhi_TTHypo;
	              std::vector<float> NonBJetE_TTHypo;
	              std::vector<float> NonBJetPt_STHypo;
	              std::vector<float> NonBJetEta_STHypo;
	              std::vector<float> NonBJetPhi_STHypo;
	              std::vector<float> NonBJetE_STHypo;

	              std::vector<float> ElectronPt;
	              std::vector<float> ElectronEta;
	              std::vector<float> ElectronPhi;
	              std::vector<float> ElectronE;
	              std::vector<float> MuonPt;
	              std::vector<float> MuonEta;
	              std::vector<float> MuonPhi;
	              std::vector<float> MuonE;
	              
                
/*
                for(int i_Jet = 0; i_Jet < selectedJets.size(); i_Jet++)
                {
                    if(i_Jet <= 4)
                }
*/
                for(int i_Jet = 0; i_Jet < selectedJets.size(); i_Jet++)
                {
                      if(i_Jet < 3 || selectedJets[i_Jet]->btag_combinedInclusiveSecondaryVertexV2BJetTags()>CSVv2_workingpointvalue_Medium)//The 3 jets with the highest CSVv2 value are used as b-jets.
                      {
                          if(i_Jet != selectedJets.size()-1)//If all jets are b-tagged, assign the last jet as non-b tagged jet
                          {
                              BJetPt_TTHypo.push_back(selectedJets[i_Jet]->Pt());
                              BJetEta_TTHypo.push_back(selectedJets[i_Jet]->Eta());
                              BJetPhi_TTHypo.push_back(selectedJets[i_Jet]->Phi());
                              BJetE_TTHypo.push_back(selectedJets[i_Jet]->E());
                          }
                          else
                          {
                              NonBJetPt_TTHypo.push_back(selectedJets[i_Jet]->Pt());
                              NonBJetEta_TTHypo.push_back(selectedJets[i_Jet]->Eta());
                              NonBJetPhi_TTHypo.push_back(selectedJets[i_Jet]->Phi());
                              NonBJetE_TTHypo.push_back(selectedJets[i_Jet]->E());
                          }
                          BJetPt_STHypo.push_back(selectedJets[i_Jet]->Pt());
                          BJetEta_STHypo.push_back(selectedJets[i_Jet]->Eta());
                          BJetPhi_STHypo.push_back(selectedJets[i_Jet]->Phi());
                          BJetE_STHypo.push_back(selectedJets[i_Jet]->E());
                      }
                      else
                      {
                          NonBJetPt_TTHypo.push_back(selectedJets[i_Jet]->Pt());
                          NonBJetEta_TTHypo.push_back(selectedJets[i_Jet]->Eta());
                          NonBJetPhi_TTHypo.push_back(selectedJets[i_Jet]->Phi());
                          NonBJetE_TTHypo.push_back(selectedJets[i_Jet]->E());
                          NonBJetPt_STHypo.push_back(selectedJets[i_Jet]->Pt());
                          NonBJetEta_STHypo.push_back(selectedJets[i_Jet]->Eta());
                          NonBJetPhi_STHypo.push_back(selectedJets[i_Jet]->Phi());
                          NonBJetE_STHypo.push_back(selectedJets[i_Jet]->E());
                      }

                      if(i_Jet < 2 || selectedJets[i_Jet]->btag_combinedInclusiveSecondaryVertexV2BJetTags()>CSVv2_workingpointvalue_Medium)//The 2 jets with the highest CSVv2 value are used as b-jets, in case the number of b-jets is samller than 2.
                      {
                          if(i_Jet != selectedJets.size()-2 && i_Jet != selectedJets.size()-1)//If all jets are b-tagged, assign the last 2 jets as non-b tagged jet
                          {
                              BJetPt_SMttHypo.push_back(selectedJets[i_Jet]->Pt());
                              BJetEta_SMttHypo.push_back(selectedJets[i_Jet]->Eta());
                              BJetPhi_SMttHypo.push_back(selectedJets[i_Jet]->Phi());
                              BJetE_SMttHypo.push_back(selectedJets[i_Jet]->E());
                          }
                          else
                          {
                              NonBJetPt_SMttHypo.push_back(selectedJets[i_Jet]->Pt());
                              NonBJetEta_SMttHypo.push_back(selectedJets[i_Jet]->Eta());
                              NonBJetPhi_SMttHypo.push_back(selectedJets[i_Jet]->Phi());
                              NonBJetE_SMttHypo.push_back(selectedJets[i_Jet]->E());
                          }
                      }
                      else
                      {
                          NonBJetPt_SMttHypo.push_back(selectedJets[i_Jet]->Pt());
                          NonBJetEta_SMttHypo.push_back(selectedJets[i_Jet]->Eta());
                          NonBJetPhi_SMttHypo.push_back(selectedJets[i_Jet]->Phi());
                          NonBJetE_SMttHypo.push_back(selectedJets[i_Jet]->E());
                      }
                }

                if(Electron)
                {
                    ElectronPt.push_back(selectedElectrons[0]->Pt());
                    ElectronEta.push_back(selectedElectrons[0]->Eta());
                    ElectronPhi.push_back(selectedElectrons[0]->Phi());
                    ElectronE.push_back(selectedElectrons[0]->E());
                    MuonPt.push_back(0.);
                    MuonEta.push_back(0.);
                    MuonPhi.push_back(0.);
                    MuonE.push_back(0.);
                }
                else if(Muon)
                {
                    MuonPt.push_back(selectedMuons[0]->Pt());
                    MuonEta.push_back(selectedMuons[0]->Eta());
                    MuonPhi.push_back(selectedMuons[0]->Phi());
                    MuonE.push_back(selectedMuons[0]->E());
                    ElectronPt.push_back(0.);
                    ElectronEta.push_back(0.);
                    ElectronPhi.push_back(0.);
                    ElectronE.push_back(0.);
                }
                
                
                kfit_STHypo->SetBJet(BJetPt_STHypo,BJetEta_STHypo,BJetPhi_STHypo,BJetE_STHypo);
                kfit_STHypo->SetNonBJet(NonBJetPt_STHypo,NonBJetEta_STHypo,NonBJetPhi_STHypo,NonBJetE_STHypo);
                kfit_STHypo->SetMet(met_px,met_py);
                kfit_STHypo->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
                kfit_STHypo->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
                
                kfit_STHypo->Run();
                
		            int NPerm_STHypo = kfit_STHypo->GetNPerm();

		            float TopLepWLepFitPt;
		            float TopLepWLepFitEta;
		            float TopLepWLepFitPhi;
		            float TopLepWLepFitE;
                if(Electron)
                {
                    TopLepWLepFitPt = selectedElectrons[0]->Pt();
                    TopLepWLepFitEta = selectedElectrons[0]->Eta();
                    TopLepWLepFitPhi = selectedElectrons[0]->Phi();
                    TopLepWLepFitE = selectedElectrons[0]->E();
                }
                else if(Muon)
                {
                    TopLepWLepFitPt = selectedMuons[0]->Pt();
                    TopLepWLepFitEta = selectedMuons[0]->Eta();
                    TopLepWLepFitPhi = selectedMuons[0]->Phi();
                    TopLepWLepFitE = selectedMuons[0]->E();
                }

		            
                if(selectedJets.size() >=  4)
                {
                    //Initialize the hypotheses
                    kfit_SMttHypo->SetBJet(BJetPt_SMttHypo,BJetEta_SMttHypo,BJetPhi_SMttHypo,BJetE_SMttHypo);
                    kfit_SMttHypo->SetNonBJet(NonBJetPt_SMttHypo,NonBJetEta_SMttHypo,NonBJetPhi_SMttHypo,NonBJetE_SMttHypo);
                    kfit_SMttHypo->SetMet(met_px,met_py);
                    kfit_SMttHypo->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
                    kfit_SMttHypo->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
                    kfit_TTHypo->SetBJet(BJetPt_TTHypo,BJetEta_TTHypo,BJetPhi_TTHypo,BJetE_TTHypo);
                    kfit_TTHypo->SetNonBJet(NonBJetPt_TTHypo,NonBJetEta_TTHypo,NonBJetPhi_TTHypo,NonBJetE_TTHypo);
                    kfit_TTHypo->SetMet(met_px,met_py);
                    kfit_TTHypo->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
                    kfit_TTHypo->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);

                    kfit_SMttHypo->Run();
                    kfit_TTHypo->Run();
		                int NPerm_SMttHypo = kfit_SMttHypo->GetNPerm();
		                int NPerm_TTHypo = kfit_TTHypo->GetNPerm();
                    //Run over SMttHypo
		                for(int ip=0;ip<NPerm_SMttHypo;ip++)
		                {
		                     float disc = kfit_SMttHypo->GetDisc(ip);

		                     int idxTopLepWElecFit = kfit_SMttHypo->GetIndex(ELECTRON_TOPTOPLEPHAD,ip);
		                     int idxTopLepWMuonFit = kfit_SMttHypo->GetIndex(MUON_TOPTOPLEPHAD,ip);
		                     int idxTopLepBJetFit = kfit_SMttHypo->GetIndex(BJETLEP_TOPTOPLEPHAD,ip);
		                     int idxTopHadBJetFit = kfit_SMttHypo->GetIndex(BJETHAD_TOPTOPLEPHAD,ip);
		                     int idxTopHadWNonBJet1Fit = kfit_SMttHypo->GetIndex(NONBJET1_TOPTOPLEPHAD,ip);
		                     int idxTopHadWNonBJet2Fit = kfit_SMttHypo->GetIndex(NONBJET2_TOPTOPLEPHAD,ip);

		                     float NuPx = kfit_SMttHypo->GetNuPx(ip,0);
		                     float NuPy = kfit_SMttHypo->GetNuPy(ip,0);
		                     float NuPz = kfit_SMttHypo->GetNuPz(ip,0);
		                     float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

		                     TLorentzVector *TopLepWNuFitP4 = new TLorentzVector();
		                     TopLepWNuFitP4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
		                     
		                     TLorentzVector *TopLepWLepFitP4 = new TLorentzVector();
		                     TopLepWLepFitP4->SetPtEtaPhiE(TopLepWLepFitPt,TopLepWLepFitEta,TopLepWLepFitPhi,TopLepWLepFitE);
		                     
		                     float TopLepBJetFitPt = BJetPt_SMttHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitEta = BJetEta_SMttHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitPhi = BJetPhi_SMttHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitE = BJetE_SMttHypo[idxTopLepBJetFit];
		                     
		                     TLorentzVector *TopLepBJetFitP4 = new TLorentzVector();
		                     TopLepBJetFitP4->SetPtEtaPhiE(TopLepBJetFitPt,TopLepBJetFitEta,TopLepBJetFitPhi,TopLepBJetFitE);
		                     
		                     float TopHadBJetFitPt = BJetPt_SMttHypo[idxTopHadBJetFit];
		                     float TopHadBJetFitEta = BJetEta_SMttHypo[idxTopHadBJetFit];
		                     float TopHadBJetFitPhi = BJetPhi_SMttHypo[idxTopHadBJetFit];
		                     float TopHadBJetFitE = BJetE_SMttHypo[idxTopHadBJetFit];
		                     
		                     TLorentzVector *TopHadBJetFitP4 = new TLorentzVector();
		                     TopHadBJetFitP4->SetPtEtaPhiE(TopHadBJetFitPt,TopHadBJetFitEta,TopHadBJetFitPhi,TopHadBJetFitE);
		                     
		                     float TopHadWNonBJet1FitPt = NonBJetPt_SMttHypo[idxTopHadWNonBJet1Fit];
		                     float TopHadWNonBJet1FitEta = NonBJetEta_SMttHypo[idxTopHadWNonBJet1Fit];
		                     float TopHadWNonBJet1FitPhi = NonBJetPhi_SMttHypo[idxTopHadWNonBJet1Fit];
		                     float TopHadWNonBJet1FitE = NonBJetE_SMttHypo[idxTopHadWNonBJet1Fit];
		                     
		                     TLorentzVector *TopHadWNonBJet1FitP4 = new TLorentzVector();
		                     TopHadWNonBJet1FitP4->SetPtEtaPhiE(TopHadWNonBJet1FitPt,TopHadWNonBJet1FitEta,TopHadWNonBJet1FitPhi,TopHadWNonBJet1FitE);
		                     
		                     float TopHadWNonBJet2FitPt = NonBJetPt_SMttHypo[idxTopHadWNonBJet2Fit];
		                     float TopHadWNonBJet2FitEta = NonBJetEta_SMttHypo[idxTopHadWNonBJet2Fit];
		                     float TopHadWNonBJet2FitPhi = NonBJetPhi_SMttHypo[idxTopHadWNonBJet2Fit];
		                     float TopHadWNonBJet2FitE = NonBJetE_SMttHypo[idxTopHadWNonBJet2Fit];

		                     TLorentzVector *TopHadWNonBJet2FitP4 = new TLorentzVector();
		                     TopHadWNonBJet2FitP4->SetPtEtaPhiE(TopHadWNonBJet2FitPt,TopHadWNonBJet2FitEta,TopHadWNonBJet2FitPhi,TopHadWNonBJet2FitE);

		                     TLorentzVector TopHadW = *TopHadWNonBJet1FitP4+*TopHadWNonBJet2FitP4;
		                     TLorentzVector TopLep = *TopLepWLepFitP4+*TopLepWNuFitP4+*TopLepBJetFitP4;
		                     TLorentzVector TopHad = TopHadW+*TopHadBJetFitP4;

		                     float VarTopHadWRecM = TopHadW.M();
		                     float VarTopLepRecM = TopLep.M();
		                     float VarTopHadRecM = TopHad.M();
		                     float VarTopLepTopHadRecDr = TopLep.DeltaR(TopHad);
		                     float VarTopLepRecPt = TopLep.Pt();
		                     float VarTopHadRecPt = TopHad.Pt();

		                     TLorentzVector *TopHadFitT = new TLorentzVector();
		                     TopHadFitT->SetPxPyPzE(TopHad.Px(),TopHad.Py(),0.,TopHad.Et());
		                     
		                     TLorentzVector *TopLepWLepFitT = new TLorentzVector();
		                     TopLepWLepFitT->SetPxPyPzE(TopLepWLepFitP4->Px(),TopLepWLepFitP4->Py(),0.,TopLepWLepFitP4->Et());

		                     TLorentzVector *TopLepWNuFitT = new TLorentzVector();
		                     TopLepWNuFitT->SetPxPyPzE(TopLepWNuFitP4->Px(),TopLepWNuFitP4->Py(),0.,TopLepWNuFitP4->Et());

		                     TLorentzVector *TopLepBJetFitT = new TLorentzVector();
		                     TopLepBJetFitT->SetPxPyPzE(TopLepBJetFitP4->Px(),TopLepBJetFitP4->Py(),0.,TopLepBJetFitP4->Et());
		                     
		                     TLorentzVector TopLepT = *TopLepWLepFitT+*TopLepWNuFitT+*TopLepBJetFitT;

		                     float VarTopLepRecMT = sqrt(2*(*TopLepWNuFitT+*TopLepBJetFitT).Pt() * TopLepWNuFitT->Pt() * (1-cos( (*TopLepWNuFitT+*TopLepBJetFitT).DeltaPhi( *TopLepWNuFitT )) ) );
		                     float VarTopLepTopHadRecDphiT = TopLepT.DeltaPhi(*TopHadFitT);
		                     float VarTopLepRecPtT = TopLepT.Pt();
		                     
		                     delete TopHadFitT;
		                     delete TopLepWLepFitT;
		                     delete TopLepWNuFitT;
		                     delete TopLepBJetFitT;
		                     
		                     delete TopLepWLepFitP4;
		                     delete TopLepWNuFitP4;
		                     delete TopLepBJetFitP4;
		                     delete TopHadBJetFitP4;
		                     delete TopHadWNonBJet1FitP4;
		                     delete TopHadWNonBJet2FitP4;

                         float MVA_tmp;
				                 if( disc < 10E+8 )
				                   {				 
				                      MVAFullReco_TopHadRecM_ = VarTopHadRecM;
				                      MVAFullReco_TopLepRecM_ = VarTopLepRecM;
				                      MVAFullReco_TopLepTopHadRecDr_ = VarTopLepTopHadRecDr;
				                      MVAFullReco_TopLepRecPt_ = VarTopLepRecPt;
				                      
				                      MVA_tmp = reader_FullReco_SMttHypo->EvaluateMVA("BDTG method");
				                   }
				                 else
				                   {
				                      MVAPartReco_TopHadRecM_ = VarTopHadRecM;
				                      MVAPartReco_TopLepRecMT_ = VarTopLepRecMT;
				                      MVAPartReco_TopLepTopHadRecDphiT_ = VarTopLepTopHadRecDphiT;
				                      MVAPartReco_TopLepRecPtT_ = VarTopLepRecPtT;
				                      
				                      MVA_tmp = reader_PartReco_SMttHypo->EvaluateMVA("BDTG method");
				                   }

				                   if(MVA_tmp > MVA_TOPTOPLEPHAD)
				                   {
				                        MVA_TOPTOPLEPHAD = MVA_tmp;

				                        for(int i_IndexMatch = 0; i_IndexMatch < selectedJets.size(); i_IndexMatch++)
				                        {
				                            if(float(pt_jet[i_IndexMatch]) == BJetPt_SMttHypo[idxTopLepBJetFit]) TOPTOPLEPHAD_JetIdx_LepTop = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) == BJetPt_SMttHypo[idxTopHadBJetFit]) TOPTOPLEPHAD_JetIdx_HadTop = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) == NonBJetPt_SMttHypo[idxTopHadWNonBJet1Fit]) TOPTOPLEPHAD_JetIdx_W1 = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) == NonBJetPt_SMttHypo[idxTopHadWNonBJet2Fit]) TOPTOPLEPHAD_JetIdx_W2 = i_IndexMatch;
        /*				                    else
				                            {
				                                cout << "An error occurred in the jet-index matching of the sorted jet collection to the original jet collection" << endl;
				                                return 1;
				                            }
        */				                
                                }
                                  
                                TopHadMass_TOPTOPLEPHAD = TopHad.M();
                                TopLepMass_TOPTOPLEPHAD = TopLep.M();
                                TopLepTopHadDr_TOPTOPLEPHAD = TopLep.DeltaR(TopHad);
                                TopLepBJetCSVv2_TOPTOPLEPHAD = CSVv2[TOPTOPLEPHAD_JetIdx_LepTop];
                                TopHadBJetCSVv2_TOPTOPLEPHAD = CSVv2[TOPTOPLEPHAD_JetIdx_HadTop];
                                TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = CSVv2[TOPTOPLEPHAD_JetIdx_W1];
                                TopHadWNonBJet2CSVv2_TOPTOPLEPHAD = CSVv2[TOPTOPLEPHAD_JetIdx_W2];

/*
    cout << " JetIndices: " << TOPTOPLEPHAD_JetIdx_W1 << " " << TOPTOPLEPHAD_JetIdx_W2 << " " << TOPTOPLEPHAD_JetIdx_LepTop << " " << TOPTOPLEPHAD_JetIdx_HadTop << endl;
    cout << " - pt_jet[TOPTOPLEPHAD_JetIdx_W1]: " << pt_jet[TOPTOPLEPHAD_JetIdx_W1] << endl;
    cout << " - pt_jet[TOPTOPLEPHAD_JetIdx_W2]: " << pt_jet[TOPTOPLEPHAD_JetIdx_W2] << endl;
    cout << " - pt_jet[TOPTOPLEPHAD_JetIdx_LepTop]: " << pt_jet[TOPTOPLEPHAD_JetIdx_LepTop] << endl;
    cout << " - pt_jet[TOPTOPLEPHAD_JetIdx_HadTop]: " << pt_jet[TOPTOPLEPHAD_JetIdx_HadTop] << endl;
    cout << " - eta_jet[TOPTOPLEPHAD_JetIdx_W1]: " << eta_jet[TOPTOPLEPHAD_JetIdx_W1] << endl;
    cout << " - eta_jet[TOPTOPLEPHAD_JetIdx_W2]: " << eta_jet[TOPTOPLEPHAD_JetIdx_W2] << endl;
    cout << " - eta_jet[TOPTOPLEPHAD_JetIdx_LepTop]: " << eta_jet[TOPTOPLEPHAD_JetIdx_LepTop] << endl;
    cout << " - eta_jet[TOPTOPLEPHAD_JetIdx_HadTop]: " << eta_jet[TOPTOPLEPHAD_JetIdx_HadTop] << endl;
    cout << " - phi_jet[TOPTOPLEPHAD_JetIdx_W1]: " << phi_jet[TOPTOPLEPHAD_JetIdx_W1] << endl;
    cout << " - phi_jet[TOPTOPLEPHAD_JetIdx_W2]: " << phi_jet[TOPTOPLEPHAD_JetIdx_W2] << endl;
    cout << " - phi_jet[TOPTOPLEPHAD_JetIdx_LepTop]: " << phi_jet[TOPTOPLEPHAD_JetIdx_LepTop] << endl;
    cout << " - phi_jet[TOPTOPLEPHAD_JetIdx_HadTop]: " << phi_jet[TOPTOPLEPHAD_JetIdx_HadTop] << endl;
    cout << " - E_jet[TOPTOPLEPHAD_JetIdx_W1]: " << E_jet[TOPTOPLEPHAD_JetIdx_W1] << endl;
    cout << " - E_jet[TOPTOPLEPHAD_JetIdx_W2]: " << E_jet[TOPTOPLEPHAD_JetIdx_W2] << endl;
    cout << " - E_jet[TOPTOPLEPHAD_JetIdx_LepTop]: " << E_jet[TOPTOPLEPHAD_JetIdx_LepTop] << endl;
    cout << " - E_jet[TOPTOPLEPHAD_JetIdx_HadTop]: " << E_jet[TOPTOPLEPHAD_JetIdx_HadTop] << endl;
    cout << "HadTopMass (TOPTOPLEPHAD): " << VarTopHadRecM << endl;
    cout << "HiggsMass (TOPTOPLEPHAD): " << VarTopHadWRecM << endl;
*/				                   }

                    }
                    //Run over TTHypo
		                for(int ip=0;ip<NPerm_TTHypo;ip++)
		                {
		                     float disc = kfit_TTHypo->GetDisc(ip);

		                     int idxTopLepWElecFit = kfit_TTHypo->GetIndex(ELECTRON_TOPTOPLEPHBB,ip);
		                     int idxTopLepWMuonFit = kfit_TTHypo->GetIndex(MUON_TOPTOPLEPHBB,ip);
		                     int idxTopLepBJetFit = kfit_TTHypo->GetIndex(BJETLEP_TOPTOPLEPHBB,ip);
		                     int idxTopHadNonBJetFit = kfit_TTHypo->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip);
		                     int idxHiggsBJet1Fit = kfit_TTHypo->GetIndex(BJET1_TOPTOPLEPHBB,ip);
		                     int idxHiggsBJet2Fit = kfit_TTHypo->GetIndex(BJET2_TOPTOPLEPHBB,ip);

		                     float NuPx = kfit_TTHypo->GetNuPx(ip,0);
		                     float NuPy = kfit_TTHypo->GetNuPy(ip,0);
		                     float NuPz = kfit_TTHypo->GetNuPz(ip,0);
		                     float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

		                     TLorentzVector *TopLepWNuFitP4 = new TLorentzVector();
		                     TopLepWNuFitP4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
		                     
		                     TLorentzVector *TopLepWLepFitP4 = new TLorentzVector();
		                     TopLepWLepFitP4->SetPtEtaPhiE(TopLepWLepFitPt,TopLepWLepFitEta,TopLepWLepFitPhi,TopLepWLepFitE);
		                     
		                     float TopLepBJetFitPt = BJetPt_TTHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitEta = BJetEta_TTHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitPhi = BJetPhi_TTHypo[idxTopLepBJetFit];
		                     float TopLepBJetFitE = BJetE_TTHypo[idxTopLepBJetFit];
		                     
		                     TLorentzVector *TopLepBJetFitP4 = new TLorentzVector();
		                     TopLepBJetFitP4->SetPtEtaPhiE(TopLepBJetFitPt,TopLepBJetFitEta,TopLepBJetFitPhi,TopLepBJetFitE);
		                     
		                     float HiggsBJet1FitPt = BJetPt_TTHypo[idxHiggsBJet1Fit];
		                     float HiggsBJet1FitEta = BJetEta_TTHypo[idxHiggsBJet1Fit];
		                     float HiggsBJet1FitPhi = BJetPhi_TTHypo[idxHiggsBJet1Fit];
		                     float HiggsBJet1FitE = BJetE_TTHypo[idxHiggsBJet1Fit];
		                     
		                     TLorentzVector *HiggsBJet1FitP4 = new TLorentzVector();
		                     HiggsBJet1FitP4->SetPtEtaPhiE(HiggsBJet1FitPt,HiggsBJet1FitEta,HiggsBJet1FitPhi,HiggsBJet1FitE);
		                     
		                     float HiggsBJet2FitPt = BJetPt_TTHypo[idxHiggsBJet2Fit];
		                     float HiggsBJet2FitEta = BJetEta_TTHypo[idxHiggsBJet2Fit];
		                     float HiggsBJet2FitPhi = BJetPhi_TTHypo[idxHiggsBJet2Fit];
		                     float HiggsBJet2FitE = BJetE_TTHypo[idxHiggsBJet2Fit];

		                     TLorentzVector *HiggsBJet2FitP4 = new TLorentzVector();
		                     HiggsBJet2FitP4->SetPtEtaPhiE(HiggsBJet2FitPt,HiggsBJet2FitEta,HiggsBJet2FitPhi,HiggsBJet2FitE);

		                     
		                     TLorentzVector Higgs = *HiggsBJet1FitP4+*HiggsBJet2FitP4;
		                     TLorentzVector TopLep = *TopLepWLepFitP4+*TopLepWNuFitP4+*TopLepBJetFitP4;

		                     float VarHiggsRecM = Higgs.M();
		                     float VarTopLepRecM = TopLep.M();
		                     float VarHiggsTopLepRecDr = Higgs.DeltaR(TopLep);
		                     float VarTopLepRecPt = TopLep.Pt();
		                     float VarHiggsRecPt = Higgs.Pt();

		                     TLorentzVector *HiggsFitT = new TLorentzVector();
		                     HiggsFitT->SetPxPyPzE(Higgs.Px(),Higgs.Py(),0.,Higgs.Et());
		                     
		                     TLorentzVector *TopLepWLepFitT = new TLorentzVector();
		                     TopLepWLepFitT->SetPxPyPzE(TopLepWLepFitP4->Px(),TopLepWLepFitP4->Py(),0.,TopLepWLepFitP4->Et());

		                     TLorentzVector *TopLepWNuFitT = new TLorentzVector();
		                     TopLepWNuFitT->SetPxPyPzE(TopLepWNuFitP4->Px(),TopLepWNuFitP4->Py(),0.,TopLepWNuFitP4->Et());

		                     TLorentzVector *TopLepBJetFitT = new TLorentzVector();
		                     TopLepBJetFitT->SetPxPyPzE(TopLepBJetFitP4->Px(),TopLepBJetFitP4->Py(),0.,TopLepBJetFitP4->Et());
		                     
		                     TLorentzVector TopLepT = *TopLepWLepFitT+*TopLepWNuFitT+*TopLepBJetFitT;

		                     float VarTopLepRecMT = sqrt(2*(*TopLepWNuFitT+*TopLepBJetFitT).Pt() * TopLepWNuFitT->Pt() * (1-cos( (*TopLepWNuFitT+*TopLepBJetFitT).DeltaPhi( *TopLepWNuFitT )) ) );
		                     float VarHiggsTopLepRecDphiT = HiggsFitT->DeltaPhi(TopLepT);
		                     float VarTopLepRecPtT = TopLepT.Pt();
		                     
		                     delete HiggsFitT;
		                     delete TopLepWLepFitT;
		                     delete TopLepWNuFitT;
		                     delete TopLepBJetFitT;
		                     
		                     delete TopLepWLepFitP4;
		                     delete TopLepWNuFitP4;
		                     delete TopLepBJetFitP4;

                         float MVA_tmp;
				                 if( disc < 10E+8 )
				                   {				 
				                      MVAFullReco_HiggsRecM_ = VarHiggsRecM;
				                      MVAFullReco_TopLepRecM_ = VarTopLepRecM;
				                      MVAFullReco_HiggsTopLepRecDr_ = VarHiggsTopLepRecDr;
				                      MVAFullReco_TopLepRecPt_ = VarTopLepRecPt;
				                      
				                      MVA_tmp = reader_FullReco_TTHypo->EvaluateMVA("BDTG method");
				                   }
				                 else
				                   {
				                      MVAPartReco_HiggsRecM_ = VarHiggsRecM;
				                      MVAPartReco_TopLepRecMT_ = VarTopLepRecMT;
				                      MVAPartReco_HiggsTopLepRecDphiT_ = VarHiggsTopLepRecDphiT;
				                      MVAPartReco_TopLepRecPtT_ = VarTopLepRecPtT;
				                      
				                      MVA_tmp = reader_PartReco_TTHypo->EvaluateMVA("BDTG method");
				                   }

				                   if(MVA_tmp > MVA_TOPTOPLEPHBB)
				                   {
				                        MVA_TOPTOPLEPHBB = MVA_tmp;
				                        for(int i_IndexMatch = 0; i_IndexMatch < selectedJets.size(); i_IndexMatch++)
				                        {
				                            if(float(pt_jet[i_IndexMatch]) == BJetPt_TTHypo[idxTopLepBJetFit]) TOPTOPLEPHBB_JetIdx_LepTop = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) == NonBJetPt_TTHypo[idxTopHadNonBJetFit]) TOPTOPLEPHBB_JetIdx_HadTop = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) == BJetPt_TTHypo[idxHiggsBJet1Fit]) TOPTOPLEPHBB_JetIdx_H1 = i_IndexMatch;
				                            else if(float(pt_jet[i_IndexMatch]) ==BJetPt_TTHypo[idxHiggsBJet2Fit]) TOPTOPLEPHBB_JetIdx_H2 = i_IndexMatch;
        /*				                    else
				                            {
				                                cout << "An error occurred in the jet-index matching of the sorted jet collection to the original jet collection" << endl;
				                                return 1;
				                            }
        */				                }
/*    cout << " JetIndices: " << TOPTOPLEPHBB_JetIdx_W1 << " " << TOPTOPLEPHBB_JetIdx_W2 << " " << TOPTOPLEPHBB_JetIdx_LepTop << " " << TOPTOPLEPHBB_JetIdx_HadTop << endl;
    cout << " - pt_jet[TOPTOPLEPHBB_JetIdx_W1]: " << pt_jet[TOPTOPLEPHBB_JetIdx_W1] << endl;
    cout << " - pt_jet[TOPTOPLEPHBB_JetIdx_W2]: " << pt_jet[TOPTOPLEPHBB_JetIdx_W2] << endl;
    cout << " - pt_jet[TOPTOPLEPHBB_JetIdx_LepTop]: " << pt_jet[TOPTOPLEPHBB_JetIdx_LepTop] << endl;
    cout << " - pt_jet[TOPTOPLEPHBB_JetIdx_HadTop]: " << pt_jet[TOPTOPLEPHBB_JetIdx_HadTop] << endl;
    cout << " - eta_jet[TOPTOPLEPHBB_JetIdx_W1]: " << eta_jet[TOPTOPLEPHBB_JetIdx_W1] << endl;
    cout << " - eta_jet[TOPTOPLEPHBB_JetIdx_W2]: " << eta_jet[TOPTOPLEPHBB_JetIdx_W2] << endl;
    cout << " - eta_jet[TOPTOPLEPHBB_JetIdx_LepTop]: " << eta_jet[TOPTOPLEPHBB_JetIdx_LepTop] << endl;
    cout << " - eta_jet[TOPTOPLEPHBB_JetIdx_HadTop]: " << eta_jet[TOPTOPLEPHBB_JetIdx_HadTop] << endl;
    cout << " - phi_jet[TOPTOPLEPHBB_JetIdx_W1]: " << phi_jet[TOPTOPLEPHBB_JetIdx_W1] << endl;
    cout << " - phi_jet[TOPTOPLEPHBB_JetIdx_W2]: " << phi_jet[TOPTOPLEPHBB_JetIdx_W2] << endl;
    cout << " - phi_jet[TOPTOPLEPHBB_JetIdx_LepTop]: " << phi_jet[TOPTOPLEPHBB_JetIdx_LepTop] << endl;
    cout << " - phi_jet[TOPTOPLEPHBB_JetIdx_HadTop]: " << phi_jet[TOPTOPLEPHBB_JetIdx_HadTop] << endl;
    cout << " - E_jet[TOPTOPLEPHBB_JetIdx_W1]: " << E_jet[TOPTOPLEPHBB_JetIdx_W1] << endl;
    cout << " - E_jet[TOPTOPLEPHBB_JetIdx_W2]: " << E_jet[TOPTOPLEPHBB_JetIdx_W2] << endl;
    cout << " - E_jet[TOPTOPLEPHBB_JetIdx_LepTop]: " << E_jet[TOPTOPLEPHBB_JetIdx_LepTop] << endl;
    cout << " - E_jet[TOPTOPLEPHBB_JetIdx_HadTop]: " << E_jet[TOPTOPLEPHBB_JetIdx_HadTop] << endl;
    cout << "HadTopMass (TOPTOPLEPHBB): " << VarTopHadRecM << endl;
    cout << "HiggsMass (TOPTOPLEPHBB): " << VarTopHadHRecM << endl;
*/
                              HiggsMass_TOPTOPLEPHBB = Higgs.M();
                              TopLepMass_TOPTOPLEPHBB = TopLep.M();
                              HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB = HiggsBJet1FitP4->DeltaR(*HiggsBJet2FitP4);
                              TopLepHiggsDr_TOPTOPLEPHBB = Higgs.DeltaR(TopLep);
                              HiggsBJet1CSVv2_TOPTOPLEPHBB = CSVv2[TOPTOPLEPHBB_JetIdx_H1];
                              HiggsBJet2CSVv2_TOPTOPLEPHBB = CSVv2[TOPTOPLEPHBB_JetIdx_H2];
                              TopLepBJetCSVv2_TOPTOPLEPHBB = CSVv2[TOPTOPLEPHBB_JetIdx_LepTop];
                              TopHadNonBJetCSVv2_TOPTOPLEPHBB = CSVv2[TOPTOPLEPHBB_JetIdx_HadTop];
				                   }

                    }
                }//End if-not-at-least-4-jets
                //Run over STHypo
		            for(int ip=0;ip<NPerm_STHypo;ip++)
		            {
		                 float disc = kfit_STHypo->GetDisc(ip);

		                 int idxTopLepWElecFit = kfit_STHypo->GetIndex(ELECTRON_TOPHLEPBB,ip);
		                 int idxTopLepWMuonFit = kfit_STHypo->GetIndex(MUON_TOPHLEPBB,ip);
		                 int idxTopLepBJetFit = kfit_STHypo->GetIndex(BJETLEP_TOPHLEPBB,ip);
		                 int idxTopHadNonBJetFit = 0;//The Had jet will be assigned as the non-b-jet in this specific hypothesis, which has index 0 in it's respective non-b-jet collection
		                 int idxHiggsBJet1Fit = kfit_STHypo->GetIndex(BJET1_TOPHLEPBB,ip);
		                 int idxHiggsBJet2Fit = kfit_STHypo->GetIndex(BJET2_TOPHLEPBB,ip);

		                 float NuPx = kfit_STHypo->GetNuPx(ip,0);
		                 float NuPy = kfit_STHypo->GetNuPy(ip,0);
		                 float NuPz = kfit_STHypo->GetNuPz(ip,0);
		                 float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

		                 TLorentzVector *TopLepWNuFitP4 = new TLorentzVector();
		                 TopLepWNuFitP4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
		                 
		                 TLorentzVector *TopLepWLepFitP4 = new TLorentzVector();
		                 TopLepWLepFitP4->SetPtEtaPhiE(TopLepWLepFitPt,TopLepWLepFitEta,TopLepWLepFitPhi,TopLepWLepFitE);
		                 
		                 float TopLepBJetFitPt = BJetPt_STHypo[idxTopLepBJetFit];
		                 float TopLepBJetFitEta = BJetEta_STHypo[idxTopLepBJetFit];
		                 float TopLepBJetFitPhi = BJetPhi_STHypo[idxTopLepBJetFit];
		                 float TopLepBJetFitE = BJetE_STHypo[idxTopLepBJetFit];
		                 
		                 TLorentzVector *TopLepBJetFitP4 = new TLorentzVector();
		                 TopLepBJetFitP4->SetPtEtaPhiE(TopLepBJetFitPt,TopLepBJetFitEta,TopLepBJetFitPhi,TopLepBJetFitE);
		                 
		                 float HiggsBJet1FitPt = BJetPt_STHypo[idxHiggsBJet1Fit];
		                 float HiggsBJet1FitEta = BJetEta_STHypo[idxHiggsBJet1Fit];
		                 float HiggsBJet1FitPhi = BJetPhi_STHypo[idxHiggsBJet1Fit];
		                 float HiggsBJet1FitE = BJetE_STHypo[idxHiggsBJet1Fit];
		                 
		                 TLorentzVector *HiggsBJet1FitP4 = new TLorentzVector();
		                 HiggsBJet1FitP4->SetPtEtaPhiE(HiggsBJet1FitPt,HiggsBJet1FitEta,HiggsBJet1FitPhi,HiggsBJet1FitE);
		                 
		                 float HiggsBJet2FitPt = BJetPt_STHypo[idxHiggsBJet2Fit];
		                 float HiggsBJet2FitEta = BJetEta_STHypo[idxHiggsBJet2Fit];
		                 float HiggsBJet2FitPhi = BJetPhi_STHypo[idxHiggsBJet2Fit];
		                 float HiggsBJet2FitE = BJetE_STHypo[idxHiggsBJet2Fit];

		                 TLorentzVector *HiggsBJet2FitP4 = new TLorentzVector();
		                 HiggsBJet2FitP4->SetPtEtaPhiE(HiggsBJet2FitPt,HiggsBJet2FitEta,HiggsBJet2FitPhi,HiggsBJet2FitE);

		                 
		                 TLorentzVector Higgs = *HiggsBJet1FitP4+*HiggsBJet2FitP4;
		                 TLorentzVector TopLep = *TopLepWLepFitP4+*TopLepWNuFitP4+*TopLepBJetFitP4;

		                 float VarHiggsRecM = Higgs.M();
		                 float VarTopLepRecM = TopLep.M();
		                 float VarHiggsTopLepRecDr = Higgs.DeltaR(TopLep);
		                 float VarTopLepRecPt = TopLep.Pt();
		                 float VarHiggsRecPt = Higgs.Pt();

		                 TLorentzVector *HiggsFitT = new TLorentzVector();
		                 HiggsFitT->SetPxPyPzE(Higgs.Px(),Higgs.Py(),0.,Higgs.Et());
		                 
		                 TLorentzVector *TopLepWLepFitT = new TLorentzVector();
		                 TopLepWLepFitT->SetPxPyPzE(TopLepWLepFitP4->Px(),TopLepWLepFitP4->Py(),0.,TopLepWLepFitP4->Et());

		                 TLorentzVector *TopLepWNuFitT = new TLorentzVector();
		                 TopLepWNuFitT->SetPxPyPzE(TopLepWNuFitP4->Px(),TopLepWNuFitP4->Py(),0.,TopLepWNuFitP4->Et());

		                 TLorentzVector *TopLepBJetFitT = new TLorentzVector();
		                 TopLepBJetFitT->SetPxPyPzE(TopLepBJetFitP4->Px(),TopLepBJetFitP4->Py(),0.,TopLepBJetFitP4->Et());
		                 
		                 TLorentzVector TopLepT = *TopLepWLepFitT+*TopLepWNuFitT+*TopLepBJetFitT;

		                 float VarTopLepRecMT = sqrt(2*(*TopLepWNuFitT+*TopLepBJetFitT).Pt() * TopLepWNuFitT->Pt() * (1-cos( (*TopLepWNuFitT+*TopLepBJetFitT).DeltaPhi( *TopLepWNuFitT )) ) );
		                 float VarHiggsTopLepRecDphiT = HiggsFitT->DeltaPhi(TopLepT);
		                 float VarTopLepRecPtT = TopLepT.Pt();
		                 
		                 delete HiggsFitT;
		                 delete TopLepWLepFitT;
		                 delete TopLepWNuFitT;
		                 delete TopLepBJetFitT;
		                 
		                 delete TopLepWLepFitP4;
		                 delete TopLepWNuFitP4;
		                 delete TopLepBJetFitP4;

                     float MVA_tmp_hut;
                     float MVA_tmp_hct;
				             if( disc < 10E+8 )
				               {				 
				                  MVAFullReco_HiggsRecM_ = VarHiggsRecM;
				                  MVAFullReco_TopLepRecM_ = VarTopLepRecM;
				                  MVAFullReco_HiggsTopLepRecDr_ = VarHiggsTopLepRecDr;
				                  MVAFullReco_TopLepRecPt_ = VarTopLepRecPt;
				                  
				                  MVA_tmp_hut = reader_FullReco_STHypo_hut->EvaluateMVA("BDTG method");
				                  MVA_tmp_hct = reader_FullReco_STHypo_hct->EvaluateMVA("BDTG method");
				               }
				             else
				               {
				                  MVAPartReco_HiggsRecM_ = VarHiggsRecM;
				                  MVAPartReco_TopLepRecMT_ = VarTopLepRecMT;
				                  MVAPartReco_HiggsTopLepRecDphiT_ = VarHiggsTopLepRecDphiT;
				                  MVAPartReco_TopLepRecPtT_ = VarTopLepRecPtT;
				                  
				                  MVA_tmp_hut = reader_PartReco_STHypo_hut->EvaluateMVA("BDTG method");
				                  MVA_tmp_hct = reader_PartReco_STHypo_hct->EvaluateMVA("BDTG method");
				               }

				               if(MVA_tmp_hut > MVA_TOPHLEPBB_hut)
				               {
				                    MVA_TOPHLEPBB_hut = MVA_tmp_hut;
				                    for(int i_IndexMatch = 0; i_IndexMatch < selectedJets.size(); i_IndexMatch++)
				                    {
				                        if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxTopLepBJetFit]) TOPHLEPBB_JetIdx_LepTop_hut = i_IndexMatch;
				                        else if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxHiggsBJet1Fit]) TOPHLEPBB_JetIdx_H1_hut = i_IndexMatch;
				                        else if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxHiggsBJet2Fit]) TOPHLEPBB_JetIdx_H2_hut = i_IndexMatch;
    /*				                    else
				                        {
				                            cout << "An error occurred in the jet-index matching of the sorted jet collection to the original jet collection" << endl;
				                            return 1;
				                        }
    */				                }
	                          HiggsMass_TOPHLEPBB_hut = Higgs.M();
	                          HiggsEta_TOPHLEPBB_hut = Higgs.Eta();
	                          TopLepMass_TOPHLEPBB_hut = TopLep.M();
                            TopLepPt_TOPHLEPBB_hut = TopLep.Pt();
                            TopLepEta_TOPHLEPBB_hut = TopLep.Eta();
                            HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut = HiggsBJet1FitP4->DeltaR(*HiggsBJet2FitP4);
                            TopLepHiggsDr_TOPHLEPBB_hut = Higgs.DeltaR(TopLep);
                            HiggsBJet1CSVv2_TOPHLEPBB_hut = CSVv2[TOPHLEPBB_JetIdx_H1_hct];
                            HiggsBJet2CSVv2_TOPHLEPBB_hut = CSVv2[TOPHLEPBB_JetIdx_H2_hct];
                            TopLepBJetCSVv2_TOPHLEPBB_hut = CSVv2[TOPHLEPBB_JetIdx_LepTop_hct];
				               }
				               if(MVA_tmp_hct > MVA_TOPHLEPBB_hct)
				               {
				                    MVA_TOPHLEPBB_hct = MVA_tmp_hct;
				                    for(int i_IndexMatch = 0; i_IndexMatch < selectedJets.size(); i_IndexMatch++)
				                    {
				                        if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxTopLepBJetFit]) TOPHLEPBB_JetIdx_LepTop_hct = i_IndexMatch;
				                        else if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxHiggsBJet1Fit]) TOPHLEPBB_JetIdx_H1_hct = i_IndexMatch;
				                        else if(float(pt_jet[i_IndexMatch]) == BJetPt_STHypo[idxHiggsBJet2Fit]) TOPHLEPBB_JetIdx_H2_hct = i_IndexMatch;
    /*				                    else
				                        {
				                            cout << "An error occurred in the jet-index matching of the sorted jet collection to the original jet collection" << endl;
				                            return 1;
				                        }
    */				                }

	                          HiggsMass_TOPHLEPBB_hct = Higgs.M();
	                          HiggsEta_TOPHLEPBB_hct = Higgs.Eta();
	                          TopLepMass_TOPHLEPBB_hct = TopLep.M();
                            TopLepPt_TOPHLEPBB_hct = TopLep.Pt();
                            TopLepEta_TOPHLEPBB_hct = TopLep.Eta();
                            HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct = HiggsBJet1FitP4->DeltaR(*HiggsBJet2FitP4);
                            TopLepHiggsDr_TOPHLEPBB_hct = Higgs.DeltaR(TopLep);
                            HiggsBJet1CSVv2_TOPHLEPBB_hct = CSVv2[TOPHLEPBB_JetIdx_H1_hct];
                            HiggsBJet2CSVv2_TOPHLEPBB_hct = CSVv2[TOPHLEPBB_JetIdx_H2_hct];
                            TopLepBJetCSVv2_TOPHLEPBB_hct = CSVv2[TOPHLEPBB_JetIdx_LepTop_hct];
				               }

                }
            }



            // Filling ntuples
            tup_ObjectVars->Fill();


			      ///////////////////////////
			      // Event info //////
			      ///////////////////////////
			      if(Muon && !Electron)
			      fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedMuons[0]->type(), selectedMuons[0]->Pt(), selectedMuons[0]->Eta(), selectedMuons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size());
			      else if(!Muon && Electron)
			      {
			        fprintf(eventlist,"%6d %6d %10d  %+2d  %6.2f %+4.2f %+4.2f   %6.1f  %+4.2f    %d %d \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->type(), selectedElectrons[0]->Pt(), selectedElectrons[0]->Eta(), selectedElectrons[0]->Phi(),mets[0]->Et(), mets[0]->Phi(),selectedJets.size(), selectedMBJets.size(), selectedElectrons[0]->superClusterEta(), selectedElectrons[0]->deltaEtaIn(), selectedElectrons[0]->deltaPhiIn(), selectedElectrons[0]->sigmaIEtaIEta_full5x5(), selectedElectrons[0]->hadronicOverEm(), selectedElectrons[0]->ioEmIoP(), r2selection.pfElectronIso(selectedElectrons[0]), selectedElectrons[0]->missingHits());
//               fprintf(eventlist,"%6d %6d %10d  %6.3f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6d %6.5f %6.5f %6.5f %6.5f %6.5f \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->superClusterEta(), selectedElectrons[0]->deltaEtaIn(), selectedElectrons[0]->deltaPhiIn(), selectedElectrons[0]->sigmaIEtaIEta_full5x5(), selectedElectrons[0]->hadronicOverEm(), selectedElectrons[0]->ioEmIoP(), r2selection.pfElectronIso(selectedElectrons[0]), selectedElectrons[0]->missingHits(), selectedElectrons[0]->chargedHadronIso(3), selectedElectrons[0]->neutralHadronIso(3), selectedElectrons[0]->photonIso(3), rho, r2selection.GetElectronIsoCorrType(selectedElectrons[0])/rho);
//                fprintf(eventlist,"%6d %6d %10d  %6.5f %6.5f %6.5f %6.5f %6.5f \n", event->runId(), event->lumiBlockId(), event->eventId(), selectedElectrons[0]->neutralHadronIso(3), selectedElectrons[0]->photonIso(3), selectedElectrons[0]->chargedHadronIso(3), rho, r2selection.GetElectronIsoCorrType(selectedElectrons[0],true)/rho);
            }

            if(debug)
            {
                cout << "End selected event" << endl;
            }

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
        //tup_ntupleinfo->Print("all");	          
        //tup_ObjectVars->Print("all");

        tupfile->cd();
        tup_ntupleinfo->Write();
        tup_ObjectVars->Write();

       	tupfile->Close();
//        delete tupfile;
//        delete tup_ntupleinfo;
//        delete tup_ObjectVars;
        cout <<"n events passed_Step2  =  "<<passed_Step2 <<endl;
        cout <<"n events passed_Step3  =  "<<passed_Step3 <<endl;
        cout <<"n events passed_Step4  =  "<<passed_Step4 <<endl;
        cout <<"n events passed_Step5  =  "<<passed_Step5 <<endl;
        cout <<"n events passed_Step6  =  "<<passed_Step6 <<endl;
        cout <<"n events passed_Step7  =  "<<passed_Step7 <<endl;
        cout <<"n events passed_Step8  =  "<<passed_Step8 <<endl;
        cout <<"n events passed_FinalSelection  =  "<<passed_FinalSelection <<endl;
        cout << "Event Count: " << eventCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();

      delete btwt_CSVv2M_mujets_central;

    } //End Loop on Datasets

    fclose (eventlist);


  

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}


