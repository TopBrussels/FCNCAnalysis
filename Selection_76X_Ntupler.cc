///////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////


#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <TLorentzVector.h>
#include <map>


#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>


#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <iostream>
#include <cstdlib>

//used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
//the header file for the triggers are taken from Lieselotte
#include "../TopTreeAnalysisBase/Tools/interface/Trigger.h"
//This header file is taken directly from the BTV wiki. It contains
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"



/////
using namespace std;
using namespace reweight;
using namespace TopTree;

map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;



int main (int argc, char *argv[])
{
    clock_t start = clock();
    
    /// Some variables from POG/PAG
    const float workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations for 76X.
    const float workingpointvalue_Medium = 0.800;//working points updated to 2015 BTV-POG recommendations.
    const float workingpointvalue_Tight = 0.935;//working points updated to 2015 BTV-POG recommendations.
    //// *** Working Conditions ////
    bool Elec_Elec, Mu_Mu, Elec_Mu, Apply_HLT_Triggers, eventSelected, Fake_Electrons, Charge_misID, ApplyElec_SF , ApplyMu_SF , ApplyPU_SF, Apply_btag_SF, Apply_JetCleaning, trigged,debug, printTrigger;
    Elec_Elec = false;
    Mu_Mu = true;
    Elec_Mu = false;
    Apply_HLT_Triggers = true;
    printTrigger = false;
    eventSelected= false;
    Fake_Electrons = false;
    Charge_misID = false;
    ApplyElec_SF = true;
    ApplyMu_SF = true;
    ApplyPU_SF = true;
    //Apply_btag_SF = false;
    Apply_JetCleaning = true;
    trigged = false;
    debug = false;
    bool applyJER = false;
    bool applyJES = false;
   // bool bTagReweight_fillBtagHisto = false;
    bool applyNegWeightCorrection = true;
    bool btagShape = true;
    bool fillBtagHisto = false;
    
    std::string channelpostfix = "";
    /////////////////////
    ///  Configuration
    /////////////////////
    
    /// xml file
    //string xmlFileName ="config/Run2SameSignDiLepton_76XSamples.xml";
    string xmlFileName ="";
    string Channel = "";
    if (argc > 1) xmlFileName = (string)argv[1];
    const char *xmlfile = xmlFileName.c_str();
    
    //Setting Lepton Channels
    if(Elec_Elec)
    {
        cout << " --> Using the Electron-Electron channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_76X_ElEl_Samples.xml";
        channelpostfix = "_ElEl_";
        Channel = "Dilepton_ElecElec";
    }
    else if(Mu_Mu)
    {
        cout << " --> Using the Muon-Muon channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_76X_MuMu_Samples.xml";
        channelpostfix = "_MuMu_";
        Channel = "Dilepton_MuMu";
    }
    else if(Elec_Mu)
    {
        cout << " --> Using the Electron-Muon channel..." << endl;
        channelpostfix = "_ElMu_";
        Channel = "Dilepton_ElecMu";
    }
    else
    {
        cerr<<"ERROR: Correct Di-lepton Channel not selected."<<endl;
        exit(1);
    }
    cout << " - Using config file " << xmlfile << endl;

    
    // placing argument in variables
    
    if(argc < 14)
    {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
        for (int n_arg=1; n_arg<argc; n_arg++)
        {
            std:: cerr << "arg number " << n_arg << " is " << argv[n_arg] << std::endl;
        }
        return 1;
    }
    
    
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
    const int JobNum                = strtol(argv[argc-3], NULL, 10);
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);
    
    vector<string> vecfileNames;
    for(int args = 11; args < argc-3; args++)
    {
        vecfileNames.push_back(argv[args]);
    }
    
    //info
    cout << "===Dataset accepted from command line===" << endl;
    cout << "Dataset Name (dName)  : " << dName << endl;
    cout << "Dataset Title (dTitle): " << dTitle << endl;
    cout << "Dataset color (color) : " << color << endl;
    cout << "Dataset Line Style (ls): " << ls << endl;
    cout << "Dataset Line Width (lw): " << lw << endl;
    cout << "Dataset Normalization Factor (normf): " << normf << endl;
    cout << "Dataset Equivalent Luminosity (EqLumi): " << EqLumi << endl;
    cout << "Dataset Xsection (xSect): " << xSect << endl;
    cout << "Dataset File Name: " << vecfileNames[0] << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "JobNum: " << JobNum << endl;
    cout << " =============================" <<endl;
    
    ////////////////////////////////////
    ///  AnalysisEnvironment
    ////////////////////////////////////
    
    //TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    
    AnalysisEnvironment anaEnv;
    cout << " - Loading environment ..." << endl;  // AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile); doesn't work on localgird
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = false;
    anaEnv.loadGenJetCollection = true;// changed on 31okt
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.JetType = 2; //0: TRootJet - 1: CaloJet - 2: PFJet - 3: JPTJet
    anaEnv.METType = 2; //0: TRootMET - 1: CaloMET - 2: PFMET - 3: TCMET

    cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
    
   // new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv); // to be understanded
    int verbose = 2; //anaEnv.Verbose;
    int nbEvents = 0;
    
    
    
    ////////////////////////////////////
    ///  DETERMINE EVENT SCALEFACTOR  ///
    /////////////////////////////////////
    
    
    
    ///// --- Define Scaling Factors ----- /////
  
    string pathToCaliDir = "../TopTreeAnalysisBase/Calibrations/";
    
    //////PU SF
    
   LumiReWeighting LumiWeights("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIIFall15DR76-Asympt25ns.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert.root", "pileup", "pileup");
    
    
    ///// lepton scaling factors
    
    //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)
    
    MuonSFWeight *muonSFWeightID = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false); // Tight ID
    
    MuonSFWeight *muonSFWeightIso = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunCD_Reco76X_Feb15.root", "MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso
    ////Triggers SF for muons to be added
    
    string electronFile= "Elec_SF_TopEA.root";
    ElectronSFWeight* electronSFWeight = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electronFile,"GlobalSF", true,false, false); // (... , ... , debug, print warning)

    
    ///// b-tagging scaling factor
    
    BTagCalibration * btagcalib;
    BTagCalibrationReader * btagreader;
    BTagCalibrationReader *reader_csvv2;
    BTagWeightTools *btwt = 0;

    
    
    
    ///-- dfeine Triggering ---//
    // /// HLT Triggers will used in the analysis is according to Top Trigger (Run2)
    // //// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
    
    //Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep);
    Trigger* DilepTrigger = new Trigger(1,0,0,1);
    
    // JER / JEC
    vector<JetCorrectorParameters> vCorrParam;
    string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";

    
    /////////////////////
    ///  Load Datasets
    /////////////////////
    
    TTreeLoader treeLoader;
    vector <Dataset*> datasets;
    cout << " - Loading datasets ..." << endl;
    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName,dTitle,true,color,ls,lw,normf,xSect,vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi*normf);
    datasets.push_back(theDataset);
    cout << "Number of datasets: " << datasets.size() << endl;
    
    
    
    
    /// //////////
    /// determine lumi
    ////////////////////////
    
    float oldLuminosity = 2612.180735004; ///2629.40667073 ; //RunD+C   //2612.180735004 Run D ;   //anaEnv.Luminosity;  // in 1/pb
    cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
    
    float Luminosity=oldLuminosity;
    
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
        string dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")==0 || dataSetName.find("data")==0 || dataSetName.find("DATA")==0)
        {
            Luminosity = datasets[d]->EquivalentLumi();
            cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        }
    }

    if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
    
    
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();
    
    //Global variable
    TRootEvent* event = 0;
    
    //nof selected events
    int nofSelectedEvents = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    
    
    
    ///////////////////////\\\\\\\\\\\\\\\\\\\\\\\
    ///// Create root file contains histograms \\\\
    /////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    string runDate = "Test_DiMu_18Apr";
    string runChannel = channelpostfix;
    //string histoPath = "OutPutHistos";
    string histoPath = "/user/sabuzeid/FCNC_Study/CMSSW_7_6_0/src/TopBrussels/FCNCAnalysis/OutputHistos";
    histoPath += runChannel+runDate+"/";
    mkdir(histoPath.c_str(),0777); // create the directory histoPath if it is not exist already
    string histoPathSampleName = channelpostfix+dName;
    mkdir((histoPath+histoPathSampleName+"/").c_str(),0777);
    string histoFileName = histoPath+histoPathSampleName+"/"+channelpostfix+dName+".root";
    if (strJobNum != "0")
    {
        if(verbose == 0) cout << "strJobNum is " << strJobNum << endl;
        histoFileName = histoPath+histoPathSampleName+"/"+channelpostfix+dName+"_"+strJobNum+".root";
    }
    TFile *fout = new TFile(histoFileName.c_str(),"RECREATE");
    
    // map <string,TH1F*> histo1D;
    std::string titlePlot = "";
    
    titlePlot = "cutFlow"+channelpostfix;
    histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 16,-0.5,15.5);
    
    
    
    //*** histos for Jets *** //
    titlePlot = "initial_Nb_Jets"+channelpostfix;
    histo1D["h_initial_Nb_Jets"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets",  16, - 0.5, 15.5 );
    
    titlePlot = "trigged_Nb_Jets"+channelpostfix;
    histo1D["h_trigged_Nb_Jets"] = new TH1F(titlePlot.c_str(), "trigged nb. of jets",  16, - 0.5, 15.5 );
    titlePlot = "2L_Nb_Jets"+channelpostfix;
    histo1D["h_2L_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of jets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_2J_1BJ_Nb_Jets"+channelpostfix;
    histo1D["h_2L_2J_1BJ_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After 2L + at least 2 jets + at least 1 CSVL jet cut: nb. of jets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_1st_Jet_Pt"+channelpostfix;
    histo1D["h_2L_1st_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 1st jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2nd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2nd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 2nd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_3rd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_3rd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 3rd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_4th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_4th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 4th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_5th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_5th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 5th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_6th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_6th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L cut: 6th jet P_{T}",  100, 0, 500 );
    titlePlot = "2L_2J_Nb_Jets"+channelpostfix;
    histo1D["h_2L_2J_Nb_Jets"] = new TH1F(titlePlot.c_str(), "After 2L + at least 2 jets cut: nb. of jets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_2J_1st_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1st_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 1st jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_2nd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_2nd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 2nd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_3rd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_3rd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 3rd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_4th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_4th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 4th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_5th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_5th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 5th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_6th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_6th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + 2J cut: 6th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_1st_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_1st_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + cut: 1st jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_2nd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_2nd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bJet cut: 2nd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_3rd_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_3rd_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bJet cut: 3rd jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_4th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_4th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bJet cut: 4th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_5th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_5th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bJet cut: 5th jet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_6th_Jet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_6th_Jet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bJet cut: 6th jet P_{T}",  100, 0, 500);
    
    titlePlot = "2L_2J_Jet_Eta"+channelpostfix;
    histo1D["h_2L_2J_Jet_Eta"]= new TH1F(titlePlot.c_str(), "After 2L + >= 2J #eta",  100, -3.0, 3.0);
    titlePlot = "2L_2J_Jet_Phi"+channelpostfix;
    histo1D["h_2L_2J_Jet_Phi"]= new TH1F(titlePlot.c_str(), "After 2L + >= 2J #Phi",  100, -4.0, 4.0);
    
    
    //////***** CSVLBJet *****////////
    titlePlot = "initial_Nb_CSVLBJets"+channelpostfix;
    histo1D["h_initial_Nb_CSVLBJets"] = new TH1F(titlePlot.c_str(), "Initial nb. of CSVLBJets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_Nb_CSVLBJets"+channelpostfix;
    histo1D["h_2L_Nb_CSVLBJets"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of CSVLBJets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_2J_Nb_CSVLBJets"+channelpostfix;
    histo1D["h_2L_2J_Nb_CSVLBJets"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of CSVLBJets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_2J_1BJ_Nb_CSVLBJets"+channelpostfix;
    histo1D["h_2L_2J_1BJ_Nb_CSVLBJets"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of CSVLBJets",  16, - 0.5, 15.5 );
    
    titlePlot = "2L_2J_1bJ_1st_CSVLBJet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_1st_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + cut: 1st CSVLBJet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_2nd_CSVLBJet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_2nd_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 2nd CSVLBJet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_3rd_CSVLBJet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_3rd_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 3rd CSVLBJet P_{T}",  100, 0, 500 );
    
    titlePlot = "2L_2J_1bJ_4th_CSVLBJet_Pt"+channelpostfix;
    histo1D["h_2L_2J_1bJ_4th_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 4th CSVLBJet P_{T}",  100, 0, 500 );
    
    
    
    //*** histos for Muons *** //
    if (Mu_Mu)
    {
        titlePlot = "initial_Nb_Mu"+channelpostfix;
        histo1D["h_initial_Nb_Mu"] = new TH1F(titlePlot.c_str(), "Initial nb. of Muons",  16, - 0.5, 15.5 );
        
        titlePlot = "trigged_Nb_Mu"+channelpostfix;
        histo1D["h_trigged_Nb_Mu"] = new TH1F(titlePlot.c_str(), "Initial nb. of Muons",  16, - 0.5, 15.5 );
        
        titlePlot = "2L_Nb_Mu"+channelpostfix;
        histo1D["h_2L_Nb_Mu"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of Muons",  16, - 0.5, 15.5 );
        
        titlePlot = "2L_1st_Mu_Pt"+channelpostfix;
        histo1D["h_2L_1st_Mu_Pt"] = new TH1F(titlePlot.c_str(), "After 2L : 1st Muon P_{T}",  100, 0, 500 );
        
        titlePlot = "2L_2nd_Mu_Pt"+channelpostfix;
        histo1D["h_2L_2nd_Mu_Pt"] = new TH1F(titlePlot.c_str(), "After 2L  cut: 2nd Muon P_{T}",  100, 0, 500);
        
        titlePlot = "2L_DeltaR_2Mu"+channelpostfix;
        histo1D["h_2L_DeltaR_2Mu"] = new TH1F(titlePlot.c_str(), "After ALL  cuts: DeltaR",  100, 0, 5);
        
        titlePlot = "2L_DeltaPhi_2Mu"+channelpostfix;
        histo1D["h_2L_DeltaPhi_2Mu"] = new TH1F(titlePlot.c_str(), "After ALL  cuts: #Delta #Phi",  100, -4., 4.);
        
        titlePlot = "2L_Mll_2Mu"+channelpostfix;
        histo1D["h_2L_Mll_2Mu"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: inv_mass_2Mu",  100, 0., 500.);
        
        titlePlot = "2L_2J_1b_Mll_2Mu"+channelpostfix;
        histo1D["h_2J_1b_2L_Mll_2Mu"] = new TH1F(titlePlot.c_str(), "After 2L+>=2J +1>=b cuts: inv_mass_2Mu",  100, 0., 500.);
        
        titlePlot = "2SSL_Mll_2Mu"+channelpostfix;
        histo1D["h_2SSL_Mll_2Mu"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: inv_mass_2SSMu",  100, 0., 500.);
        
        titlePlot = "2OSL_Mll_2Mu"+channelpostfix;
        histo1D["h_2OSL_Mll_2Mu"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: inv_mass_2OSMu",  100, 0., 500.);
        
        titlePlot = "2SSL_DeltaR_2Mu"+channelpostfix;
        histo1D["h_2SSL_DeltaR_2Mu"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta R",  100, 0, 5);
        
        titlePlot = "2SSL_DeltaPhi_2Mu"+channelpostfix;
        histo1D["h_2SSL_DeltaPhi_2Mu"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta #Phi",  100, -4., 4.);
        
        titlePlot = "2OSL_DeltaR_2Mu"+channelpostfix;
        histo1D["h_2OSL_DeltaR_2Mu"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R",  100, 0, 5);
        
        titlePlot = "2OSL_DeltaPhi_2Mu"+channelpostfix;
        histo1D["h_2OSL_DeltaPhi_2Mu"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta #Phi",  100, -4., 4.);
        
        
        titlePlot = "2L_SumPT_2Mu"+channelpostfix;
        histo1D["h_2L_SumPT_2Mu"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2Mu",  100, 0., 500.);
        
        titlePlot = "2SSL_SumPT_2Mu"+channelpostfix;
        histo1D["h_2SSL_SumPT_2Mu"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2SSMu",  100, 0., 500.);
        titlePlot = "2OSL_SumPT_2Mu"+channelpostfix;
        histo1D["h_2OSL_SumPT_2Mu"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2OSMu",  100, 0., 500.);
        
        titlePlot = "2L_DeltaR_b0_Mu0"+channelpostfix;
        histo1D["h_2L_DeltaR_b0_Mu0"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: #Delta R b_{0} Mu_{0}",  100, 0, 5);
        
        titlePlot = "2L_mMu0b0"+channelpostfix;
        histo1D["h_2L_mMu0b0"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Mass b_{0} Mu_{0}",  100, 0., 500);
        
        titlePlot = "2SSL_DeltaR_b0_Mu0"+channelpostfix;
        histo1D["h_2SSL_DeltaR_b0_Mu0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta R b_{0} Mu_{0}",  100, 0, 5);
        
        titlePlot = "2SSL_mMu0b0"+channelpostfix;
        histo1D["h_2SSL_mMu0b0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: Mass b_{0} Mu_{0}",  100, 0., 500);
        
        titlePlot = "2OSL_DeltaR_b0_Mu0"+channelpostfix;
        histo1D["h_2OSL_DeltaR_b0_Mu0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R b_{0} Mu_{0}",  100, 0, 5);
        
        titlePlot = "2OSL_mMu0b0"+channelpostfix;
        histo1D["h_2OSL_mMu0b0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: Mass b_{0} Mu_{0}",  100, 0., 500);
        
    }
    
    if (Elec_Elec)
    {
        titlePlot = "initial_Nb_Elec"+channelpostfix;
        histo1D["h_initial_Nb_Elec"] = new TH1F(titlePlot.c_str(), "Initial nb. of Electrons",  16, - 0.5, 15.5 );
        
        titlePlot = "trigged_Nb_Elec"+channelpostfix;
        histo1D["h_trigged_Nb_Elec"] = new TH1F(titlePlot.c_str(), "Trigged nb. of Electrons",  16, - 0.5, 15.5 );
        
        titlePlot = "2L_Nb_Elec"+channelpostfix;
        histo1D["h_2L_Nb_Elec"] = new TH1F(titlePlot.c_str(), "After 2L cut: nb. of Electrons",  16, - 0.5, 15.5 );
        
        titlePlot = "2L_1st_Elec_Pt"+channelpostfix;
        histo1D["h_2L_1st_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2L : 1st Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2L_2nd_Elec_Pt"+channelpostfix;
        histo1D["h_2L_2nd_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2L  cut: 2nd Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2SSL_1st_Elec_Pt"+channelpostfix;
        histo1D["h_2SSL_1st_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2SSL : 1st Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2SSL_2nd_Elec_Pt"+channelpostfix;
        histo1D["h_2SSL_2nd_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2SSL  cut: 2nd Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2OSL_1st_Elec_Pt"+channelpostfix;
        histo1D["h_2OSL_1st_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2OSL : 1st Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2OSL_2nd_Elec_Pt"+channelpostfix;
        histo1D["h_2OSL_2nd_Elec_Pt"] = new TH1F(titlePlot.c_str(), "After 2OSL  cut: 2nd Electron P_{T}",  100, 0, 300);
        
        titlePlot = "2L_DeltaR_2Elec"+channelpostfix;
        histo1D["h_2L_DeltaR_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: #Delta R",  100, 0, 5);
        
        titlePlot = "2L_DeltaPhi_2Elec"+channelpostfix;
        histo1D["h_2L_DeltaPhi_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: #Delta #Phi",  100, -4., 4.);
        
        titlePlot = "2SSL_DeltaR_2Elec"+channelpostfix;
        histo1D["h_2SSL_DeltaR_2Elec"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta R",  100, 0, 5);
        
        titlePlot = "2SSL_DeltaPhi_2Elec"+channelpostfix;
        histo1D["h_2SSL_DeltaPhi_2Elec"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta #Phi",  100, -4., 4.);
        
        titlePlot = "2OSL_DeltaR_2Elec"+channelpostfix;
        histo1D["h_2OSL_DeltaR_2Elec"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R",  100, 0, 5);
        
        titlePlot = "2OSL_DeltaPhi_2Elec"+channelpostfix;
        histo1D["h_2OSL_DeltaPhi_2Elec"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta #Phi",  100, -4., 4.);
        
        titlePlot = "2L_Mll_2Elec"+channelpostfix;
        histo1D["h_2L_Mll_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: inv_mass_2Elec",  100, 0., 500.);
        
        titlePlot = "2L_2J_1b_Mll_2Elec"+channelpostfix;
        histo1D["h_2J_1b_2L_Mll_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L+>=2J +1>=b cuts: inv_mass_2Elec",  100, 0., 500.);
        
        titlePlot = "2SSL_Mll_2Elec"+channelpostfix;
        histo1D["h_2SSL_Mll_2Elec"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: inv_mass_2SSElec",  100, 0., 500.);
        
        titlePlot = "2OSL_Mll_2Elec"+channelpostfix;
        histo1D["h_2OSL_Mll_2Elec"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: inv_mass_2OSElec",  100, 0., 500.);
        
        titlePlot = "2L_SumPT_2Elec"+channelpostfix;
        histo1D["h_2L_SumPT_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2Elec",  100, 0., 500.);
        
        titlePlot = "2SSL_SumPT_2Elec"+channelpostfix;
        histo1D["h_2SSL_SumPT_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2SSElec",  100, 0., 500.);
        titlePlot = "2OSL_SumPT_2Elec"+channelpostfix;
        histo1D["h_2OSL_SumPT_2Elec"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Sum P_{T} 2OSElec",  100, 0., 500.);
        
        titlePlot = "2L_DeltaR_b0_Elec0"+channelpostfix;
        histo1D["h_2L_DeltaR_b0_Elec0"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: #Delta R b_{0} Elec_{0}",  100, 0, 5);
        
        titlePlot = "2L_mElec0b0"+channelpostfix;
        histo1D["h_2L_mElec0b0"] = new TH1F(titlePlot.c_str(), "After 2L  cuts: Mass b_{0} Elec_{0}",  100, 0., 500);
        
        titlePlot = "2SSL_DeltaR_b0_Elec0"+channelpostfix;
        histo1D["h_2SSL_DeltaR_b0_Elec0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: #Delta R b_{0} Elec_{0}",  100, 0, 5);
        
        titlePlot = "2SSL_mElec0b0"+channelpostfix;
        histo1D["h_2SSL_mElec0b0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: Mass b_{0} Elec_{0}",  100, 0., 500);
        
        titlePlot = "2OSL_DeltaR_b0_Elec0"+channelpostfix;
        histo1D["h_2OSL_DeltaR_b0_Elec0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R b_{0} Elec_{0}",  100, 0, 5);
        
        titlePlot = "2OSL_mElec0b0"+channelpostfix;
        histo1D["h_2OSL_mElec0b0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: Mass b_{0} Elec_{0}",  100, 0., 500);
        
    }

    
    ///*** histos for mets ***///
    //        titlePlot = "initial_met"+channelpostfix;
    //        histo1D["h_initial_met"] = new TH1F(titlePlot.c_str(), " initial missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
    //
    //        titlePlot = "2L_met"+channelpostfix;
    //        histo1D["h_2L_met"] = new TH1F(titlePlot.c_str(), "After 2L missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
    //        titlePlot = "2L_2J_met"+channelpostfix;
    //        histo1D["h_2L_2J_met"] = new TH1F(titlePlot.c_str(), "After 2L + >=2J missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
    //        titlePlot = "2L_2J_1bJ_met"+channelpostfix;
    //        histo1D["h_2L_2J_1bJ_met"] = new TH1F(titlePlot.c_str(), "After 2L + >=2J + 1bJet missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
    
    titlePlot = "Ht_AllJets"+channelpostfix;
    histo1D["h_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts:  H_{T}",  150, 0, 1500 );
    
    titlePlot = "St_AllJets+Lep"+channelpostfix;
    histo1D["h_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts:  S_{T}",  150, 0, 1500 );
    
    titlePlot = "2SSL_Ht_AllJets"+channelpostfix;
    histo1D["h_2SSL_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts 2SSL:  H_{T}",  150, 0, 1500 );
    
    titlePlot = "2SSL_St_AllJets+Lep"+channelpostfix;
    histo1D["h_2SSL_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts 2SSL:  S_{T}",  150, 0, 1500 );
    
    titlePlot = "2OSL_Ht_AllJets"+channelpostfix;
    histo1D["h_2OSL_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts 2OSL:  H_{T}",  150, 0, 1500 );
    
    titlePlot = "2OSL_St_AllJets+Lep"+channelpostfix;
    histo1D["h_2OSL_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts 2OSL:  S_{T}",  150, 0, 1500 );
    
    
    
    
    ////////////////////////////////////
    ///  Loop on datasets
    ////////////////////////////////////
    
    if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
    
    
    ///////----Start looping over datasets -----/////
    bool nlo = false;
    
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
        string previousFilename = "";
        int iFile = -1;
        bool isData;
        string dataSetName = datasets[d]->Name();
        if (verbose > 0)
        cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
        if (verbose > 0)
        cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
        double nloSF = 1;
        double sumWeights = 0;
        nlo = false;
        
        double lumiWeight = -99.;
        if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        {
            lumiWeight=1;
            isData = true;
        }
        else
        {
            lumiWeight = Luminosity/datasets[d]->EquivalentLumi();
            cout << "the weight to apply for each event of this data set is " << "Lumi / (EquivalentLumi) = ( "  << Luminosity << " / " << datasets[d]->EquivalentLumi() << ") = " << Luminosity/datasets[d]->EquivalentLumi()  <<  endl;
            
            ///// -- for mc@Nlo samples the correction for negative weight should be applied and this have to be done before any other SF and
            //////// also before applying any selection cuts
            if(dataSetName.find("amc")!=string::npos) nlo = true;
            
        }
        
        if (verbose >1)cout << "the lumiweight of the dataset : "<< datasets[d]->Name() << " is " << lumiWeight << endl;
        
        
        //// --- Calibration - Applying Jet correction (Jec) --- ////
        
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
        
        
        // string BCaliPath = CaliPath + "BTagging/CSVv2_13TeV_25ns_combToMujets.csv";
        if(!isData && !btagShape)
        {
            // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
            //	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv");
            btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_76X_combToMujets.csv");
            btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_LOOSE, "mujets","central");
            if(fillBtagHisto)  // before btag reweighting can be apply, you first have to make the histograms
            {
                btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_central.root",30,999,2.4);
                //btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",false,30,999,2.4);
            }
            else
            {
                btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",false,30,999,2.4);
                // btwt = new BTagWeightTools(btagreader,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
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
        

        
        
        //// ***************** /////
        /// output Ntuples /////
        /// ***************** ////
        string rootTreesDir = "/user/sabuzeid/FCNC_Study/CMSSW_7_6_0/src/TopBrussels/FCNCAnalysis/Output_Ntuples/";
        mkdir(rootTreesDir.c_str(),0777);
        string rootTreespath = rootTreesDir+Channel+"_";
        rootTreespath+=runDate+"/";
        mkdir(rootTreespath.c_str(),0777);
        string roottreename = rootTreespath+"FCNC_2SSL_";
        roottreename+= datasets[d]->Name();
        roottreename += channelpostfix;
        roottreename+= strJobNum;
        roottreename+="_tree.root";
        
        cout << "  - Recreate outputfile ... " << roottreename.c_str() << endl;
        
        // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
       // fileout->cd();
        TTree *bookkeeping = new TTree("startevents","startevents");
        TTree* InitialTree = new TTree("Initialtree","Initialtree");
        TTree* myTree = new TTree("tree","tree");
        TTree* SSLeptonTree = new TTree("SSLeptonTree","SSLeptonTree");
        
        TTree* OSLeptonTree = new TTree("OSLeptonTree","OSLeptonTree");
        TTree* globalTree = new TTree("globaltree","globaltree");

        
        //////////////////////////////
        // My tree - variables //
        //////////////////////////////
        
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        Int_t nEv;
        // various weights
        Double_t puSF;
        Double_t btagSF;
        Int_t nofPosWeights;
        Int_t nofNegWeights;
        Double_t nloWeight; // for amc@nlo samples
        Int_t sumW;
        Double_t MuonIDSF[10];
        Double_t MuonIsoSF[10];
        Double_t sf_muon[10];
        Double_t sf_electron[10];
        
        //// Variables for met
        Double_t met_Pt;
        Double_t met_Eta;
        Double_t met_Phi;
        
        //// Variables for leptons
        Int_t nLeptons;
        
        // Variables for Electrons
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
        Double_t pt_electron_1;
        Double_t pt_electron_2;
        
        //variable for muons
        Int_t nMuons;
        Double_t pt_muon[10];
        Double_t phi_muon[10];
        Double_t eta_muon[10];
        Double_t E_muon[10];
//        Double_t d0_muon[10];
//        Double_t d0BeamSpot_muon[10];
//        Double_t chargedHadronIso_muon[10];
//        Double_t neutralHadronIso_muon[10];
//        Double_t photonIso_muon[10];
//        Double_t relIso_muon[10];
//        Bool_t isId_muon[10];
//        Bool_t isIso_muon[10];
//        Double_t pfIso_muon[10];
        
        Int_t charge_muon[10];
        Double_t pt_1st_Muon;
        Double_t pt_2nd_Muon;
        Double_t eta_1st_Muon;
        Double_t eta_2nd_Muon;
        Double_t phi_1st_Muon;
        Double_t phi_2nd_Muon;
        Double_t E_1st_Muon;
        Double_t E_2nd_Muon;
        Int_t charge_1st_Muon;
        Int_t charge_2nd_Muon;
        
        
        //variable for jets
        Int_t nJets;
        Int_t nCSVLbJets;
        Int_t nCSVMbJets;
        Int_t nCSVTbJets;
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Int_t charge_jet[20];
        
        
        
        
        ///// --- Variables used for MVA --- //
        Double_t Ht;
        Double_t St;
        Double_t DeltaR_2L;
        Double_t DeltaPhi_2L;
        Double_t invMass_2L;
        Double_t DeltaR_Mu0b0;
        Double_t Mass_WJetPair;
        Double_t Mass_JetPair;
    
        
        ///////////////////////////////
        // My trees                  //
        ///////////////////////////////
        
        bookkeeping->Branch("run_num",&run_num,"run_num/I");
        bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
        bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
        bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
        bookkeeping->Branch("npu",&npu,"npu/I");
        
        globalTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        globalTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        globalTree->Branch("nEv" , &nEv, "nEv/I");
        globalTree->Branch("sumW", &sumW, "sumW/I");
        
        // define the output trees may I need to make many trees depend on selection cuts
        /////(Integer variables)
        
        myTree->Branch("isData",&isData,"isData/I");
        myTree->Branch("run_num",&run_num,"run_num/I");
        myTree->Branch("evt_num",&evt_num,"evt_num/I");
        myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        myTree->Branch("nvtx",&nvtx,"nvtx/I");
        myTree->Branch("npu",&npu,"npu/I");
        myTree->Branch("puSF",&puSF,"puSF/D");
        myTree->Branch("btagSF",&btagSF,"btagSF/D");
        myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        myTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        myTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        myTree->Branch("sumW", &sumW, "sumW/I");
        
        
        InitialTree->Branch("isData",&isData,"isData/I");
        InitialTree->Branch("run_num",&run_num,"run_num/I");
        InitialTree->Branch("evt_num",&evt_num,"evt_num/I");
        InitialTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        InitialTree->Branch("nvtx",&nvtx,"nvtx/I");
        InitialTree->Branch("npu",&npu,"npu/I");
        InitialTree->Branch("puSF",&puSF,"puSF/D");
        InitialTree->Branch("btagSF",&btagSF,"btagSF/D");
        InitialTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        
        
        SSLeptonTree->Branch("isData",&isData,"isData/I");
        SSLeptonTree->Branch("run_num",&run_num,"run_num/I");
        SSLeptonTree->Branch("evt_num",&evt_num,"evt_num/I");
        SSLeptonTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        SSLeptonTree->Branch("nvtx",&nvtx,"nvtx/I");
        SSLeptonTree->Branch("npu",&npu,"npu/I");
        SSLeptonTree->Branch("puSF",&puSF,"puSF/D");
        SSLeptonTree->Branch("btagSF",&btagSF,"btagSF/D");
        SSLeptonTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        SSLeptonTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        SSLeptonTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        SSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        
        OSLeptonTree->Branch("isData",&isData,"isData/I");
        OSLeptonTree->Branch("run_num",&run_num,"run_num/I");
        OSLeptonTree->Branch("evt_num",&evt_num,"evt_num/I");
        OSLeptonTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        OSLeptonTree->Branch("nvtx",&nvtx,"nvtx/I");
        OSLeptonTree->Branch("npu",&npu,"npu/I");
        OSLeptonTree->Branch("puSF",&puSF,"puSF/D");
        OSLeptonTree->Branch("btagSF",&btagSF,"btagSF/D");
        OSLeptonTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        OSLeptonTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        OSLeptonTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        OSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//

        
        // Set branches for different Trees
        ////--> Electrons <----////
        InitialTree ->Branch("nElectrons",&nElectrons, "nElectrons/I");
        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");

       // myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
        myTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        myTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        myTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        myTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
//        myTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
//        myTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
//        myTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
//        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
//        myTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
//        myTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
//        myTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
//        myTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
//        myTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
//        myTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
//        myTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
//        myTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
//        myTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
//        myTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
//        myTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        myTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        myTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/D");
        myTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/D");
       
        
        SSLeptonTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        SSLeptonTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        SSLeptonTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        SSLeptonTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        SSLeptonTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        SSLeptonTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        SSLeptonTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        //        SSLeptonTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        //        SSLeptonTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        //        SSLeptonTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        //        SSLeptonTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        //        SSLeptonTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        //        SSLeptonTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        SSLeptonTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        SSLeptonTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/D");
        SSLeptonTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/D");
        
        
        OSLeptonTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        OSLeptonTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        OSLeptonTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        OSLeptonTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        OSLeptonTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        OSLeptonTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        OSLeptonTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        //        OSLeptonTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        //        OSLeptonTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        //        OSLeptonTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        //        OSLeptonTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        //        OSLeptonTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        //        OSLeptonTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        OSLeptonTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        OSLeptonTree->Branch("pt_electron_1",&pt_electron_1,"pt_electron_1/D");
        OSLeptonTree->Branch("pt_electron_2",&pt_electron_2,"pt_electron_2/D");
        
        
        
        ////--> muons <----////
        InitialTree->Branch("nMuons",&nMuons, "nMuons/I");
        InitialTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        InitialTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //InitialTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // InitialTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        InitialTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        InitialTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        InitialTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        InitialTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        //        InitialTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        //        InitialTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        //        InitialTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        //        InitialTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        //        InitialTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        //        InitialTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        InitialTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        //        InitialTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        //        InitialTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        InitialTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");

        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
       // myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        myTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        myTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        myTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
//        myTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
//        myTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
//        myTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
//        myTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
//        myTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
//        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
//        myTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
//        myTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        myTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        myTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/D");
        myTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/D");
        myTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/D");
        myTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/D");
        myTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/I");
        myTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/D");
        myTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/D");
        myTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/D");
        myTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/D");
        myTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/I");
       
        
        SSLeptonTree->Branch("nMuons",&nMuons, "nMuons/I");
        SSLeptonTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        SSLeptonTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //SSLeptonTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // SSLeptonTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        SSLeptonTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        SSLeptonTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        SSLeptonTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        SSLeptonTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        //        SSLeptonTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        //        SSLeptonTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        //        SSLeptonTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        //        SSLeptonTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        //        SSLeptonTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        //        SSLeptonTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        SSLeptonTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        //        SSLeptonTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        //        SSLeptonTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        SSLeptonTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        SSLeptonTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/D");
        SSLeptonTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/D");
        SSLeptonTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/D");
        SSLeptonTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/D");
        SSLeptonTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/I");
        SSLeptonTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/D");
        SSLeptonTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/D");
        SSLeptonTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/D");
        SSLeptonTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/D");
        SSLeptonTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/I");

        OSLeptonTree->Branch("nMuons",&nMuons, "nMuons/I");
        OSLeptonTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        OSLeptonTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //OSLeptonTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // OSLeptonTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        OSLeptonTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        OSLeptonTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        OSLeptonTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        OSLeptonTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        //        OSLeptonTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        //        OSLeptonTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        //        OSLeptonTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        //        OSLeptonTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        //        OSLeptonTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        //        OSLeptonTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        OSLeptonTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        //        OSLeptonTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        //        OSLeptonTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        OSLeptonTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        OSLeptonTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/D");
        OSLeptonTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/D");
        OSLeptonTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/D");
        OSLeptonTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/D");
        OSLeptonTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/I");
        OSLeptonTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/D");
        OSLeptonTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/D");
        OSLeptonTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/D");
        OSLeptonTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/D");
        OSLeptonTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/I");

        
        

        ////--> jets <----////
        
        myTree->Branch("nJets",&nJets,"nJets/I");
        myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        
        InitialTree->Branch("nJets",&nJets,"nJets/I");
        InitialTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        InitialTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        InitialTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        InitialTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        InitialTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        
        SSLeptonTree->Branch("nJets",&nJets,"nJets/I");
        SSLeptonTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        SSLeptonTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        SSLeptonTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        SSLeptonTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        SSLeptonTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        
        OSLeptonTree->Branch("nJets",&nJets,"nJets/I");
        OSLeptonTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        OSLeptonTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        OSLeptonTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        OSLeptonTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        OSLeptonTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        
        ////--> met <----////
        myTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
        myTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
        myTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
        
        InitialTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
        InitialTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
        InitialTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
        
        SSLeptonTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
        SSLeptonTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
        SSLeptonTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
        
        OSLeptonTree->Branch("met_Pt", &met_Pt, "met_Pt/D");
        OSLeptonTree->Branch("met_Eta", &met_Eta,"met_Eta/D");
        OSLeptonTree->Branch("met_Phi", &met_Phi, "met_Phi/D");
        
        ////--> bJets <----////
        myTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        myTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        myTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        
        InitialTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        InitialTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        InitialTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        
        SSLeptonTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        SSLeptonTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        SSLeptonTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        
        OSLeptonTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        OSLeptonTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        OSLeptonTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        
        
        ///// --> MVA variables <--- ////
        myTree->Branch("Ht",&Ht,"Ht/D");
        myTree->Branch("St",&St,"St/D");
        myTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/D");
        myTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/D");
        myTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/D");
        myTree->Branch("DeltaR_Mu0b0",&DeltaR_Mu0b0,"DeltaR_Mu0b0/D");
        myTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/D");
        myTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/D");
        
        SSLeptonTree->Branch("Ht",&Ht,"Ht/D");
        SSLeptonTree->Branch("St",&St,"St/D");
        SSLeptonTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/D");
        SSLeptonTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/D");
        SSLeptonTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/D");
        SSLeptonTree->Branch("DeltaR_Mu0b0",&DeltaR_Mu0b0,"DeltaR_Mu0b0/D");
        SSLeptonTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/D");
        SSLeptonTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/D");
        
        OSLeptonTree->Branch("Ht",&Ht,"Ht/D");
        OSLeptonTree->Branch("St",&St,"St/D");
        OSLeptonTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/D");
        OSLeptonTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/D");
        OSLeptonTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/D");
        OSLeptonTree->Branch("DeltaR_Mu0b0",&DeltaR_Mu0b0,"DeltaR_Mu0b0/D");
        OSLeptonTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/D");
        OSLeptonTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/D");


        

        ////////////////////////////////////
        ///  Loop on events
        ////////////////////////////////////
        //some bookkeeping variables
        //nEvents[d] = 0;
        nbEvents = 0;
        int previousRun = -1;
        int currentRun;
        bool Btagged = false;
        nofPosWeights = 0;
        nofNegWeights = 0;
        
        int itrigger = -1;
        
        
        if (Apply_HLT_Triggers) {
            DilepTrigger->bookTriggers(isData);
        }
        
        
        //// Vectors for objects that will be used during the analysis
        vector < TRootVertex* > vertex;
        vector < TRootMuon* > init_muons;
        vector < TRootElectron* > init_electrons;
        vector < TRootJet* > init_jets_corrected;
        vector < TRootJet* > init_jets;
        vector < TRootMET* > mets;
        vector < TRootGenJet* > genjets;
        vector<TRootMCParticle*> mcParticles;
        vector<TRootPFJet*> selectedJets;
        vector < TRootMuon* > selectedMuons;
        vector < TRootElectron* > selectedElectrons;
        vector<TRootJet*> selectedBCSVLJets;
        vector<TRootJet*> selectedBCSVMJets;
        vector<TRootJet*> selectedBCSVTJets;
        vector<TRootJet*> JetsExcludingHighestCSVLb;
        float rho;
        
        
        /////// ************* ///////////
        /// start looping on events ///
        ////////************** /////////
        
        int event_start = startEvent;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();
        cout <<"Number of events = "<<  ending  <<endl;
        
        //double currentfrac =0.;
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

        //for (unsigned int ievt = 0; ievt < 1000; ievt++) // used for test code to run over few numbers of events
        if (verbose >0) cout<< "looping over events ... " << endl;
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
        {
           // nEvents[d]++;
            selectedElectrons.clear();
            selectedMuons.clear();
            selectedJets.clear();
            vertex.clear();
            init_muons.clear();
            init_jets.clear();
            init_jets_corrected.clear();
            genjets.clear();
            mets.clear();
            selectedBCSVLJets.clear();
            selectedBCSVMJets.clear();
            selectedBCSVTJets.clear();
            
            
            
            if (ievt%1000 == 0)std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
            ////////////////////
            ///  LOAD EVENT  ///
            ///////////////////
            if (verbose >2) cout<< "load event" << endl;

            //TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
            event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets, debug); //updated by adding 5th argument to check analysis at 3 Feb
            if(!isData) genjets = treeLoader.LoadGenJet(ievt,false);  //needed for JER
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if (previousFilename != currentFilename)
            {
                previousFilename = currentFilename;
                iFile++;
                cout << "File changed!!! => iFile = " << iFile << endl;
            }
            
            currentRun = event->runId();
            rho = event->fixedGridRhoFastjetAll(); //// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2#Rho_effective_area_corrections
            
            if (verbose >3) cout << "The current Run ID  =  "<< currentRun << endl;
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
            bookkeeping->Fill();

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
                   // histo1D["weightIndex"]->Fill(1.,1.);
                    
                }
                else
                {
                    if(nlo) hasNegWeight = true;
                   // cout << "hasNegWeight" << endl;
                    nofNegWeights++;
                   // histo1D["weightIndex"]->Fill(-1.,1.);
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
                  //  histo1D["weightIndex"]->Fill(2.,1.);
                    
                }
                else
                {
                    if(nlo) hasNegWeight = true;
                    nofNegWeights++;
                  //  histo1D["weightIndex"]->Fill(-2.,1.);
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
               // histo1D["nloweight"]->Fill(mc_baseweight, 1.);
                sumWeights += mc_baseweight;
              //  if (verbose > 1) cout << "mc_baseweight  =  " << mc_baseweight << endl;
              //  if (verbose > 1) cout << "sumWeights  =  " << sumWeights << endl;
                
            }
            
           //// Trigger.checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
            if (Apply_HLT_Triggers)
            {
                DilepTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
                itrigger = DilepTrigger->checkIfFired();
            }
            
            ////////////////////////////
            ///// JES - JER smearing     ////
            //////////////////////////
            if(applyJER && !isData)
            {
                jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
            }
            if(applyJES && !isData)
            {
                jetTools->correctJets(init_jets,event->fixedGridRhoFastjetAll() ,false);
            }
           
            /////////////////////////
            ///  EVENT SELECTION  ///
            /////////////////////////
            
            //Declare selection instance
            
            Run2Selection selection(init_jets, init_muons, init_electrons, mets, rho); ///rho = event->fixedGridRhoFastjetAll(); and it is used in electron isolation
            
            
            //// choose good primary vertex
            bool isGoodPV = selection.isPVSelected(vertex, 4 , 24. ,2.); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
            
            // Jets Selection
            selectedJets = selection.GetSelectedJets(30., 2.4, true , "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
            
            
            /// --- Muons Selection -- ///
            selectedMuons = selection.GetSelectedMuons(15. , 2.1 , .15 ,"Tight","Spring15");  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
            
            //// --- Electron Selection --- ///
            selectedElectrons = selection.GetSelectedElectrons(15., 2.4 , "Tight" , "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
            
           // cout << "the nb of selected electrons is " << selectedElectrons.size()<<endl;
            /// For MC Information
            mcParticles.clear();
            treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
            
            

            //// sorting objects in the the event according to Pt
            
            sort(selectedJets.begin(), selectedJets.end(),HighestPt());
            sort(selectedMuons.begin(), selectedMuons.end(), HighestPt());
            sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt());
            sort(mcParticles.begin(),mcParticles.end(),HighestPt());
            
            ///// --- Jet Cleaning -- ////
            
//            
//            if(Apply_JetCleaning)
//            {
//                if(debug) cout << "Applying jet cleaning " << endl;
//                int OrigSize = selectedJets.size();
//                for (int origJets=0; origJets<selectedJets.size(); origJets++)
//                {
//                    bool erased = false;
//                    if(selectedMuons.size()>0){
//                        if(selectedJets[origJets]->DeltaR(*selectedMuons[0])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                    if(selectedMuons.size()>1 && !erased){
//                        if(selectedJets[origJets]->DeltaR(*selectedMuons[1])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                    if(selectedMuons.size()>2 && !erased){
//                        if(selectedJets[origJets]->DeltaR(*selectedMuons[2])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                    if(selectedElectrons.size()>0 && !erased){
//                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[0])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                    if(selectedElectrons.size()>1 && !erased){
//                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[1])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                    if(selectedElectrons.size()>2 && !erased){
//                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[2])<0.4){ selectedJets.erase(selectedJets.begin()+origJets); erased = true;}
//                    }
//                }
//                if(debug)
//                {
//                    if( OrigSize != selectedJets.size()) cout << "--> original jet collection size = " << OrigSize  << "  And after cleaning =  " << selectedJets.size() << endl;
//                    else cout << "--> no change" << endl;
//                }
//            }
            
            if(Apply_JetCleaning)
            {
                if(debug) cout << "Applying jet cleaning " << endl;
                int OrigSize = selectedJets.size();
                for (int origJets=selectedJets.size()-1; origJets>=0; origJets--)
                {
                    bool erased = false;
                    if(selectedMuons.size()>0){
                        if(selectedJets[origJets]->DeltaR(*selectedMuons[0])<0.4){selectedJets[origJets]=selectedJets.back();selectedJets.pop_back(); erased = true;}
                    }
                    if(selectedMuons.size()>1 && !erased){
                        if(selectedJets[origJets]->DeltaR(*selectedMuons[1])<0.4){ selectedJets[origJets]=selectedJets.back();selectedJets.pop_back(); erased = true;}
                    }
                    if(selectedMuons.size()>2 && !erased){
                        if(selectedJets[origJets]->DeltaR(*selectedMuons[2])<0.4){ selectedJets[origJets]=selectedJets.back();selectedJets.pop_back(); erased = true;}
                    }
                    if(selectedElectrons.size()>0 && !erased){
                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[0])<0.4){ selectedJets[origJets]=selectedJets.back();selectedJets.pop_back(); erased = true;}
                    }
                    if(selectedElectrons.size()>1 && !erased){
                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[1])<0.4){ selectedJets[origJets]=selectedJets.back();selectedJets.pop_back();erased = true;}
                    }
                    if(selectedElectrons.size()>2 && !erased){
                        if(selectedJets[origJets]->DeltaR(*selectedElectrons[2])<0.4){ selectedJets[origJets]=selectedJets.back();selectedJets.pop_back(); erased = true;}
                    }
                }
                if(debug)
                {
                    if( OrigSize != selectedJets.size()) cout << "--> original jet collection size = " << OrigSize  << "  And after cleaning =  " << selectedJets.size() << endl;
                    else cout << "--> no change" << endl;
                }
            }


            /////---- bTagging ----\\\\\\
            
            for(unsigned int i = 0; i < selectedJets.size() ; i++)
            {
                
                TRootJet* tempJet = (TRootJet*) selectedJets[i];
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose)//loose WP
                {
                    
                    selectedBCSVLJets.push_back(tempJet);
                    Btagged = true;
                }
                
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium)//medium WP
                {
                    selectedBCSVMJets.push_back(tempJet);
                }
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight)//tight WP
                {
                    selectedBCSVTJets.push_back(tempJet);
                }
            }
            
            if(verbose > 2) cout << "btagging done" << endl;
            /// sorting bTag Jets
            sort(selectedBCSVLJets.begin(),selectedBCSVLJets.end(),HighestPt());
            sort(selectedBCSVMJets.begin(),selectedBCSVMJets.end(),HighestPt());
            sort(selectedBCSVTJets.begin(),selectedBCSVTJets.end(),HighestPt());
            
            
            //// create a vector containing jets after excluding highest CSVL b-tagged jet
            JetsExcludingHighestCSVLb.clear();
            
            for (unsigned int ijet = 0; ijet < selectedJets.size() ; ijet++)
            {
                if (selectedBCSVLJets.size() != 0)
                {
                    if (fabs(selectedJets[ijet]->Pt() - selectedBCSVLJets[0]->Pt()) > .000000001)
                    {
                        JetsExcludingHighestCSVLb.push_back(selectedJets[ijet]);
                    }
                    
                }
            }
            
            ////===========================///////
            //// Applying Scale Factors    //////
            ////==========================//////
            
            // scale factor for the event
            Double_t scaleFactor = 1.;
            Double_t Elec_scaleFactor = 1.;
            Double_t Muon_scaleFactor = 1.;
            Double_t MuonID_scaleFactor = 1.;
            Double_t MuonIso_scaleFactor = 1.;
            Double_t PU_scaleFactor =1.;
            //Double_t bTag_scaleFactor =1.;
            float muon1SF, muon2SF,muonID1SF, muonID2SF,muonIso1SF, muonIso2SF,electron1SF, electron2SF, muon3SF;
            muon1SF =  muon2SF = electron1SF = electron2SF = muon3SF = 0.;
            
            histo1D["h_cutFlow"]->Fill(0., scaleFactor*lumiWeight); //// fill histogram before applying any scaling factors or triggers
            
            
            
            ////initial histograms after applying scaling factors or triggers
            //histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
            histo1D["h_initial_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
            //histo1D["h_initial_Nb_PV"]->Fill(vertex.size(),scaleFactor*lumiWeight);
            //histo1D["h_initial_Nb_Elec"]->Fill(selectedElectrons.size(),scaleFactor*lumiWeight);
           // histo1D["h_initial_Nb_Mu"]->Fill(selectedMuons.size(),scaleFactor*lumiWeight);
            histo1D["h_initial_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
            // histo1D["h_initial_met"]->Fill(met_pt, scaleFactor*lumiWeight);
           
            
            
            /////applyNegWeightCorrection
            if(hasNegWeight && applyNegWeightCorrection && !isData) scaleFactor *= -1.;
            
            ////// PU SF
            
            float puWeight = 1;
            // if (verbose>2) cout << " isData =  "<< isData << endl;
            if (ApplyPU_SF)
            {
                
                if(!isData)
                {
                    puWeight = LumiWeights.ITweight((int)event->nTruePU()); // simplest reweighting, just use reconstructed number of PV. faco
                    PU_scaleFactor =puWeight;
                    // if (verbose>1) cout << "while puWeight  =  "<< puWeight << endl;
                }
                else if (isData) {PU_scaleFactor =1;}
                // if (verbose>1) cout << "PU_scaleFactor is " << PU_scaleFactor << endl;
                scaleFactor *= PU_scaleFactor;
                
            }
            histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
            
            //if (verbose>1) cout << "puSF " << puSF << endl;

            
            ///// btagging SF
            
            float btagWeight  =  1.;
            float bTagEff = 1.;
            if(fillBtagHisto && !isData && !btagShape)
            {
                btwt->FillMCEfficiencyHistos(selectedJets);
                
            }
            else if(!fillBtagHisto && !isData && !btagShape)
            {
                btagWeight =  btwt->getMCEventWeight(selectedJets);
                
            }
            else if(!isData && btagShape)
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
            
            scaleFactor *= btagWeight;
            
            
            /////////////////////////////////
            //// *** Apply selections *** ///
            ////////////////////////////////
            nbEvents++;
            bool diElectron = false;
            bool diMuon = false;
            bool diEMu = false;
            bool SSdiLepton = false;
            bool OSdiLepton = false;
           // bool Event_passed = false;
            double InvMass_ll = 0.;
            double InvMass_lb = 0.;
            const double Zmass = 91.1876; // GeV SM xsections at 13 TeV https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
            const double Wmass = 80.398; // GeV SM xsections at 13 TeV https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
            TLorentzVector tempLepton_0;
            TLorentzVector tempLepton_1;
            TLorentzVector Lepton0;
            TLorentzVector Lepton1;
            TLorentzVector SM_bJet;
            //float mjj=0. ;
//            float massDiff_from_Wmass = 0.;
//            float MinimassDiff_from_Wmass = 999;
//            float WPairMass,PairMass;
            int qLepton0 , qLepton1, qElec0, qElec1, qMu0, qMu1;
            qLepton0 = qLepton1 = qElec0 = qElec1 =  qMu0 = qMu1 = 0;
            double sum_jet_PT = 0.;
            double Ht_AllJets = 0.;
            double St_AllJets_Leps = 0.;
            double Sum_Leptons_Pt =0.;
            double Zmass_Window = 0.;
            DeltaR_2L = 0.;
            DeltaPhi_2L = 0.;
            invMass_2L = 0.;
            DeltaR_Mu0b0= 0.;
            Mass_WJetPair = 0.;
            Mass_JetPair = 0.;
            
            

            /// Trigger
            if(Apply_HLT_Triggers) trigged = itrigger; //trigged = treeLoader.EventTrigged(itrigger);
            //else trigged = true;
           // cout << "trigged = " << trigged << endl;
            if (!trigged) continue;
            
            // if (trigged)
            // {
            histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
            histo1D["h_trigged_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
            //histo1D["h_trigged_Nb_PV"]->Fill(vertex.size(),scaleFactor*lumiWeight);
            //histo1D["h_initial_Nb_Elec"]->Fill(selectedElectrons.size(),scaleFactor*lumiWeight);
            //histo1D["h_trigged_Nb_Mu"]->Fill(selectedMuons.size(),scaleFactor*lumiWeight);
            // if (verbose>2) cout << "the nb. of initial muons after pass triggers is  =  " << init_muons.size() << endl;
            // if (verbose>2) cout << "the nb. of muons after pass triggers is  =  " << selectedMuons.size() << endl;
            
            //histo1D["h_Nb_PV_AfterPU_Weight"]->Fill(vertex.size(),scaleFactor*lumiWeight);
            
            if (!isGoodPV) continue;
            
            //  if (isGoodPV)
            // {
            if(verbose>2) cout << "GoodPV" << endl;
            histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
            
            int numLep = selectedElectrons.size()+ selectedMuons.size();
            nMuons = 0;
            nLeptons = 0;
            nElectrons = 0;
            for (Int_t selmu =0; selmu < selectedMuons.size() ; selmu++ )
            {
                
                pt_muon[nMuons]=selectedMuons[selmu]->Pt();
                phi_muon[nMuons]=selectedMuons[selmu]->Phi();
                eta_muon[nMuons]=selectedMuons[selmu]->Eta();
                E_muon[nMuons]=selectedMuons[selmu]->E();
                //                        d0_muon[nMuons]=selectedMuons[selmu]->d0();
                //                        d0BeamSpot_muon[nMuons]=selectedMuons[selmu]->d0BeamSpot();
                //                        chargedHadronIso_muon[nMuons]=selectedMuons[selmu]->chargedHadronIso(4);
                //                        neutralHadronIso_muon[nMuons]=selectedMuons[selmu]->neutralHadronIso(4);
                //                        photonIso_muon[nMuons]=selectedMuons[selmu]->photonIso(4);
                //                        pfIso_muon[nMuons]=selectedMuons[selmu]->relPfIso(4,0);
                charge_muon[nMuons]=selectedMuons[selmu]->charge();
                if(selectedMuons.size()>0) pt_1st_Muon = selectedMuons[0]->Pt();
                if(selectedMuons.size()>1) pt_2nd_Muon = selectedMuons[1]->Pt();
                if(selectedMuons.size()>0) phi_1st_Muon = selectedMuons[0]->Phi();
                if(selectedMuons.size()>1) phi_2nd_Muon = selectedMuons[1]->Phi();
                if(selectedMuons.size()>0) eta_1st_Muon = selectedMuons[0]->Eta();
                if(selectedMuons.size()>1) eta_2nd_Muon = selectedMuons[1]->Eta();
                if(selectedMuons.size()>0) charge_1st_Muon = selectedMuons[0]->charge();
                if(selectedMuons.size()>1) charge_2nd_Muon = selectedMuons[1]->charge();
                
                ///// Applying muon scaling factors
                
                if(ApplyMu_SF && !isData)
                {
                    
                    sf_muon[nMuons]= muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    Muon_scaleFactor = muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    MuonIDSF[nMuons] = muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    MuonIsoSF[nMuons] =  muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    //		MuonTrigSFv2[nMuons] = muonSFWeightTrigHLTv4p2->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    //		MuonTrigSFv3[nMuons] = muonSFWeightTrigHLTv4p3->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                }
                else
                {
                    sf_muon[nMuons] = 1.;
                    Muon_scaleFactor = 1.;
                    MuonIDSF[nMuons] = 1.;
                    MuonIsoSF[nMuons] = 1.;
                    //		MuonTrigSFv2[nMuons] = 1.;
                    //		MuonTrigSFv3[nMuons] = 1.;
                    
                }
//                if(ApplyMu_SF && !isData)
//                {
//                    MuonIDSF[nMuons] = muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
//                    MuonIsoSF[nMuons] =  muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
//                    //		MuonTrigSFv2[nMuons] = muonSFWeightTrigHLTv4p2->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
//                    //		MuonTrigSFv3[nMuons] = muonSFWeightTrigHLTv4p3->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
//                }
//                else
//                {
//                    MuonIDSF[nMuons] = 1.;
//                    MuonIsoSF[nMuons] = 1.;
//                    //		MuonTrigSFv2[nMuons] = 1.;
//                    //		MuonTrigSFv3[nMuons] = 1.;
//                }
                
                nMuons++;
            }
            scaleFactor *= Muon_scaleFactor;
            
            for (unsigned int selelec = 0 ; selelec > selectedElectrons.size() ; selelec++)
            {
                    pt_electron[nElectrons]=selectedElectrons[selelec]->Pt();
                    phi_electron[nElectrons]=selectedElectrons[selelec]->Phi();
                    eta_electron[nElectrons]=selectedElectrons[selelec]->Eta();
                    eta_superCluster_electron[nElectrons]=selectedElectrons[selelec]->superClusterEta();
                    E_electron[nElectrons]=selectedElectrons[selelec]->E();
//                    d0_electron[nElectrons]=selectedElectrons[selelec]->d0();
//                    d0BeamSpot_electron[nElectrons]=selectedElectrons[selelec]->d0BeamSpot();
//                    chargedHadronIso_electron[nElectrons]=selectedElectrons[selelec]->chargedHadronIso(3);
//                    neutralHadronIso_electron[nElectrons]=selectedElectrons[selelec]->neutralHadronIso(3);
//                    photonIso_electron[nElectrons]=selectedElectrons[selelec]->photonIso(3);
//                    pfIso_electron[nElectrons]=selectedElectrons[selelec]->relPfIso(3,0);
                    charge_electron[nElectrons]=selectedElectrons[selelec]->charge();
//                    sigmaIEtaIEta_electron[nElectrons]=selectedElectrons[selelec]->sigmaIEtaIEta();
//                    deltaEtaIn_electron[nElectrons]=selectedElectrons[selelec]->deltaEtaIn();
//                    deltaPhiIn_electron[nElectrons]=selectedElectrons[selelec]->deltaPhiIn();
//                    hadronicOverEm_electron[nElectrons]=selectedElectrons[selelec]->hadronicOverEm();
//                    missingHits_electron[nElectrons]=selectedElectrons[selelec]->missingHits();
//                    passConversion_electron[nElectrons]=selectedElectrons[selelec]->passConversion();
//                    isEBEEGap[nElectrons]=selectedElectrons[selelec]->isEBEEGap();
                    if(selectedElectrons.size()>0) pt_electron_1 = selectedElectrons[0]->Pt();
                    if(selectedElectrons.size()>1) pt_electron_2 = selectedElectrons[1]->Pt();
                
                if (ApplyElec_SF && !isData)
                {
                    sf_electron[nElectrons]= electronSFWeight->at(selectedElectrons[selelec]->Eta(),selectedElectrons[selelec]->Pt(),0);
                    
                }
                else sf_electron[nElectrons]= 1.;
                
                    nElectrons++;
            }
            
            nJets = 0;
            for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
            {
                
                pt_jet[nJets]=selectedJets[seljet]->Pt();
                phi_jet[nJets]=selectedJets[seljet]->Phi();
                eta_jet[nJets]=selectedJets[seljet]->Eta();
                E_jet[nJets]=selectedJets[seljet]->E();
                charge_jet[nJets]=selectedJets[seljet]->charge();
                // bdisc_jet[nJets]=selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() ;
                nJets++;
                
            }
            nCSVTbJets = selectedBCSVTJets.size();
            nCSVMbJets = selectedBCSVMJets.size();
            nCSVLbJets = selectedBCSVLJets.size();
            double met_px = mets[0]->Px();
            double met_py = mets[0]->Py();
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
            met_Phi = mets[0]->Phi();
            met_Eta = mets[0]->Eta();
            
            puSF = PU_scaleFactor;
            btagSF = btagWeight;
            
            nLeptons = nMuons + nElectrons;
            InitialTree -> Fill();
            
            if (numLep==2)
            {
                histo1D["h_cutFlow"]->Fill(4., scaleFactor*lumiWeight);
                
                if (Elec_Elec && !Mu_Mu && !Elec_Mu) // Electron-Electron Channel
                {
                    if (selectedElectrons.size()==2 && selectedMuons.size()==0 && selectedElectrons[0]->Pt()>25. && selectedElectrons[1]->Pt()>= 15.)
                    {
                        tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                        //qLepton0 = selectedElectrons[0]->charge();
                        qElec0 = selectedElectrons[0]->charge();
                        tempLepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
                        //qLepton1 = selectedElectrons[1]->charge();
                        qElec1=selectedElectrons[1]->charge();
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        if (InvMass_ll >12.)
                        {
                            diElectron = true;
                            histo1D["h_2L_1st_Elec_Pt"]->Fill(tempLepton_0.Pt(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2nd_Elec_Pt"]->Fill(tempLepton_1.Pt(),scaleFactor*lumiWeight);
                            //histo1D["h_2L_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                            //histo1D["h_2L_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                            histo1D["h_2L_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                            invMass_2L = InvMass_ll;
                        }
                        
                    }
                    
                }
                
                else if (Mu_Mu && !Elec_Elec && !Elec_Mu) // Muon-Muon channel
                {
                    
                    if (selectedMuons.size()==2 && selectedElectrons.size()==0 && selectedMuons[0]->Pt() >= 20. && selectedMuons[1]->Pt()>= 20.)
                    {
                        
                        histo1D["h_cutFlow"]->Fill(5., scaleFactor*lumiWeight);
                        tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                        qMu0 = selectedMuons[0]->charge();
                        tempLepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
                        qMu1 = selectedMuons[1]->charge();
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        
                        // Zmass_Window = fabs(Zmass - InvMass_ll);
                        if (InvMass_ll >12.) //&& Zmass_Window >= 15.)
                        {
                            diMuon = true;
                            histo1D["h_2L_1st_Mu_Pt"]->Fill(tempLepton_0.Pt(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2nd_Mu_Pt"]->Fill(tempLepton_1.Pt(),scaleFactor*lumiWeight);
                            // histo1D["h_2L_DeltaR_2Mu"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                            // histo1D["h_2L_DeltaPhi_2Mu"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                            histo1D["h_2L_Mll_2Mu"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                            histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
                            invMass_2L = InvMass_ll;
                            
                        }
                        
                        
                    }
                }
                
                //                        else if (Elec_Mu) // Electron- Muon channel
                //                        {
                //                            if (selectedElectrons.size()==1 && selectedMuons.size()==1)
                //                            {
                //                                if (selectedMuons[0]->Pt() > selectedElectrons[0]->Pt()&& selectedMuons[0]->Pt()>= 20. && selectedElectrons[0]->Pt()>=15.)
                //                                {
                //                                    tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                //                                    qLepton0 = selectedMuons[0]->charge();
                //                                    tempLepton_1.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                //                                    qLepton1 = selectedElectrons[0]->charge();
                //
                //                                }
                //                                if (selectedElectrons[0]->Pt()>selectedMuons[0]->Pt() && selectedElectrons[0]->Pt()>=26. && selectedMuons[0]->Pt())
                //                                {
                //                                    tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                //                                    qLepton0 = selectedElectrons[0]->charge();
                //                                    tempLepton_1.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                //                                    qLepton1 = selectedMuons[0]->charge();
                //                                }
                //                                InvMass_ll = (tempLepton_0+tempLepton_1).M();
                //                                if (InvMass_ll >12.){diEMu = true;}
                //
                //                            }
                //
                //                        }
               // cout << "dielectron  =  " << diElectron << endl;
                if (diElectron || diMuon || diEMu)
                {
                    histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                    histo1D["h_2L_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
                    histo1D["h_2L_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 0)histo1D["h_2L_1st_Jet_Pt"]->Fill(selectedJets[0]->Pt(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 1)histo1D["h_2L_2nd_Jet_Pt"]->Fill(selectedJets[1]->Pt(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 2)histo1D["h_2L_3rd_Jet_Pt"]->Fill(selectedJets[2]->Pt(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 3)histo1D["h_2L_4th_Jet_Pt"]->Fill(selectedJets[3]->Pt(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 4)histo1D["h_2L_5th_Jet_Pt"]->Fill(selectedJets[4]->Pt(),scaleFactor*lumiWeight);
                    if(selectedJets.size() > 5)histo1D["h_2L_6th_Jet_Pt"]->Fill(selectedJets[5]->Pt(),scaleFactor*lumiWeight);
                    if(Elec_Elec)histo1D["h_2L_Nb_Elec"]->Fill(selectedElectrons.size(),scaleFactor*lumiWeight);
                    if(Mu_Mu)histo1D["h_2L_Nb_Mu"]->Fill(selectedMuons.size(),scaleFactor*lumiWeight);
                    // histo1D["h_2L_met"]->Fill(met_pt, scaleFactor*lumiWeight);
                    //histo2D["h_2L_Nb_bjets_vs_Nb_jets"]->Fill(selectedBCSVLJets.size(), selectedJets.size(),scaleFactor*lumiWeight);
                    //if(Elec_Elec) histo2D["h_2L_Lep0Pt_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                    
                    if (selectedJets.size()>= 2)
                    {
                        histo1D["h_cutFlow"]->Fill(8., scaleFactor*lumiWeight);
                        histo1D["h_2L_2J_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
                        histo1D["h_2L_2J_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                        histo1D["h_2L_2J_1st_Jet_Pt"]->Fill(selectedJets[0]->Pt(),scaleFactor*lumiWeight);
                        histo1D["h_2L_2J_2nd_Jet_Pt"]->Fill(selectedJets[1]->Pt(),scaleFactor*lumiWeight);
                        if(selectedJets.size() > 2)histo1D["h_2L_2J_3rd_Jet_Pt"]->Fill(selectedJets[2]->Pt(),scaleFactor*lumiWeight);
                        if(selectedJets.size() > 3)histo1D["h_2L_2J_4th_Jet_Pt"]->Fill(selectedJets[3]->Pt(),scaleFactor*lumiWeight);
                        if(selectedJets.size() > 4)histo1D["h_2L_2J_5th_Jet_Pt"]->Fill(selectedJets[4]->Pt(),scaleFactor*lumiWeight);
                        if(selectedJets.size() > 5)histo1D["h_2L_2J_6th_Jet_Pt"]->Fill(selectedJets[5]->Pt(),scaleFactor*lumiWeight);
                        // histo1D["h_2L_2J_met"]->Fill(met_pt, scaleFactor*lumiWeight);
                        // histo2D["h_2L_2J_Nb_bjets_vs_Nb_jets"]->Fill(selectedBCSVLJets.size(), selectedJets.size(),scaleFactor*lumiWeight);
                        
                        
                        //eventSelected = true;
                        
                        if (selectedBCSVLJets.size()>= 1)
                        {
                            //myTree->Fill();
                            histo1D["h_cutFlow"]->Fill(9., scaleFactor*lumiWeight);
                            SM_bJet.SetPxPyPzE(selectedBCSVLJets[0]->Px(),selectedBCSVLJets[0]->Py(),selectedBCSVLJets[0]->Pz(),selectedBCSVLJets[0]->Energy());
                            histo1D["h_2L_2J_1BJ_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2J_1BJ_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2J_1bJ_1st_Jet_Pt"]->Fill(selectedJets[0]->Pt(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2J_1bJ_2nd_Jet_Pt"]->Fill(selectedJets[1]->Pt(),scaleFactor*lumiWeight);
                            if(selectedJets.size() > 2)histo1D["h_2L_2J_1bJ_3rd_Jet_Pt"]->Fill(selectedJets[2]->Pt(),scaleFactor*lumiWeight);
                            if(selectedJets.size() > 3)histo1D["h_2L_2J_1bJ_4th_Jet_Pt"]->Fill(selectedJets[3]->Pt(),scaleFactor*lumiWeight);
                            if(selectedJets.size() > 4)histo1D["h_2L_2J_1bJ_5th_Jet_Pt"]->Fill(selectedJets[4]->Pt(),scaleFactor*lumiWeight);
                            if(selectedJets.size() > 5)histo1D["h_2L_2J_1bJ_6th_Jet_Pt"]->Fill(selectedJets[5]->Pt(),scaleFactor*lumiWeight);
                            histo1D["h_2L_2J_1bJ_1st_CSVLBJet_Pt"]->Fill(selectedBCSVLJets[0]->Pt(),scaleFactor*lumiWeight);
                            if(selectedBCSVLJets.size()>1) histo1D["h_2L_2J_1bJ_2nd_CSVLBJet_Pt"]->Fill(selectedBCSVLJets[1]->Pt(),scaleFactor*lumiWeight);
                            if(selectedBCSVLJets.size()>2) histo1D["h_2L_2J_1bJ_3rd_CSVLBJet_Pt"]->Fill(selectedBCSVLJets[2]->Pt(),scaleFactor*lumiWeight);
                            if(selectedBCSVLJets.size()>3) histo1D["h_2L_2J_1bJ_4th_CSVLBJet_Pt"]->Fill(selectedBCSVLJets[3]->Pt(),scaleFactor*lumiWeight);
                            //histo1D["h_2L_2J_1bJ_met"]->Fill(met_pt, scaleFactor*lumiWeight);
                            //histo2D["h_2L_2j_1b_Nb_bjets_vs_Nb_jets"]->Fill(selectedBCSVLJets.size(), selectedJets.size(),scaleFactor*lumiWeight);
                            
                           // TLorentzVector temp_W_1st_Jet;
                            //TLorentzVector temp_W_2nd_Jet;
                            float PairMass = 0.;
                            bool W_RecoMass = false;
                            
                            if (selectedJets.size()>=3)
                            {
                            
                                if (JetsExcludingHighestCSVLb.size() >= 2)
                            {
                               // cout <<" selectedJets.size() =  " << selectedJets.size() << endl;
                              //  cout << "JetsExcludingHighestCSVLb.size()  =  " << JetsExcludingHighestCSVLb.size() << endl;
                                //cout<< "the size of JetsExcludingHighestCSVLb  =   " << JetsExcludingHighestCSVLb.size() <<endl;
                                TLorentzVector temp_W_1st_Jet;
                                TLorentzVector temp_W_2nd_Jet;
                                double mjj = 0.;
                                double massDiff_from_Wmass = 0.;
                                double MinimassDiff_from_Wmass = 999;
                                //float WPairMass,PairMass;
                                for (unsigned int xjet =JetsExcludingHighestCSVLb.size()-1; xjet > 0 ; xjet--)
                                {
                                   // cout<< "xjet =  " << xjet <<endl;
                                  //  temp_W_1st_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[xjet]->Px(),JetsExcludingHighestCSVLb[xjet]->Py(),JetsExcludingHighestCSVLb[xjet]->Pz(),JetsExcludingHighestCSVLb[xjet]->Energy());
                                    for (unsigned int kjet = 0 ; kjet < xjet; kjet++)
                                    {
                                        //cout<< "xjet =  " << xjet <<endl;
                                       // cout << " kjet =  " << kjet << endl;
                                        temp_W_1st_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[xjet]->Px(),JetsExcludingHighestCSVLb[xjet]->Py(),JetsExcludingHighestCSVLb[xjet]->Pz(),JetsExcludingHighestCSVLb[xjet]->Energy());
                                        temp_W_2nd_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[kjet]->Px(),JetsExcludingHighestCSVLb[kjet]->Py(),JetsExcludingHighestCSVLb[kjet]->Pz(),JetsExcludingHighestCSVLb[kjet]->Energy());
                                        mjj = (temp_W_1st_Jet+ temp_W_2nd_Jet).M();
                                        Mass_JetPair = mjj;
                                      //  cout << "the value of mjj is =  " << mjj << endl;
                                       // cout << "the value of Mass_JetPair is =  " << Mass_JetPair << endl;
                                        massDiff_from_Wmass = fabs(Wmass - mjj);
                                        if (massDiff_from_Wmass < MinimassDiff_from_Wmass)
                                        {
                                            PairMass = mjj;
                                            MinimassDiff_from_Wmass = massDiff_from_Wmass;
                                            W_RecoMass = true;
                                           // cout << "mass of W jets pair (inside loop ) =  "<< mjj << endl;
                                        }
                                    }
                                }
                               // Mass_WJetPair = PairMass;
                               //cout << "mass of W jets pair (outside loop ) =  "<< PairMass << endl;
                            }
                            }
                            
                            //if (W_RecoMass) {
                                Mass_WJetPair = PairMass;
                           // }
                           // cout << "mass of W jets pair (outside loop ) =  "<< PairMass << endl;
                            
                            if (diElectron && !diMuon)
                            {
                                histo1D["h_2J_1b_2L_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2L_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                                InvMass_lb = (tempLepton_0+SM_bJet).M();
                                //cout << " InvMass_lb =  " << InvMass_lb;
                                histo1D["h_2L_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaR_b0_Elec0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                                DeltaR_2L = tempLepton_0.DeltaR(tempLepton_1);
                                DeltaPhi_2L = tempLepton_0.DeltaPhi(tempLepton_1);
                                DeltaR_Mu0b0 = tempLepton_0.DeltaR(SM_bJet);
                            }
                            if (diMuon && !diElectron)
                            {
                                histo1D["h_2J_1b_2L_Mll_2Mu"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaPhi_2Mu"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaR_2Mu"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2L_SumPT_2Mu"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                                InvMass_lb = (tempLepton_0+SM_bJet).M();
                                //cout << " InvMass_lb =  " << InvMass_lb;
                                histo1D["h_2L_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2L_DeltaR_b0_Mu0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                                DeltaR_2L = tempLepton_0.DeltaR(tempLepton_1);
                                DeltaPhi_2L = tempLepton_0.DeltaPhi(tempLepton_1);
                                DeltaR_Mu0b0 = tempLepton_0.DeltaR(SM_bJet);
                            }
                            
                            for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                            {
                                histo1D["h_2L_2J_Jet_Eta"]->Fill(selectedJets[ijet]->Eta(),scaleFactor*lumiWeight);
                                histo1D["h_2L_2J_Jet_Phi"]->Fill(selectedJets[ijet]->Phi(),scaleFactor*lumiWeight);
                                sum_jet_PT+= selectedJets[ijet]->Pt();
                            }
                            Ht_AllJets=sum_jet_PT;
                            St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                            histo1D["h_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                            histo1D["h_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                            
                            myTree->Fill();
                            
                            // eventSelected = true;
                            
                            //if (Elec_Elec && selectedElectrons[0]->charge()== selectedElectrons[1]->charge())// in case of Same Sign dilepton
                            if (diElectron && qElec0 == qElec1)
                            {
                                SSdiLepton = true;
                                histo1D["h_2SSL_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2SSL_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                                histo1D["h_2SSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaR_b0_Elec0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                            }
                            if (diElectron && qElec0 != qElec1) // in case of Opposite Sign dilepton
                            {
                                OSdiLepton = true;
                                histo1D["h_2OSL_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2OSL_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                                //InvMass_lb = (tempLepton_0+SM_bJet).M();
                                // cout << " InvMass_lb =  " << InvMass_lb;
                                histo1D["h_2OSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaR_b0_Elec0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                                
                            }
                            if (diMuon && !diElectron && qMu0 == qMu1)
                            {
                                SSdiLepton = true;
                                histo1D["h_2SSL_Mll_2Mu"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaPhi_2Mu"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaR_2Mu"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2SSL_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2SSL_DeltaR_b0_Mu0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                                histo1D["h_2SSL_SumPT_2Mu"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                            }
                            if (diMuon && !diElectron && qMu0 != qMu1) // in case of Opposite Sign dilepton
                            {
                                OSdiLepton = true;
                                histo1D["h_2OSL_Mll_2Mu"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaPhi_2Mu"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaR_2Mu"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                histo1D["h_2OSL_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                histo1D["h_2OSL_DeltaR_b0_Mu0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
                                histo1D["h_2OSL_SumPT_2Mu"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                            }
                            //                                    if (Elec_Mu && qLepton0 == qLepton1)
                            //                                    {
                            //                                        SSdiLepton = true;
                            //                                    }
                            //                                    if (Elec_Mu && qLepton0 != qLepton1) // in case of Opposite Sign dilepton
                            //                                    {
                            //                                        OSdiLepton = true;
                            //                                    }
                            
                            if (SSdiLepton && !OSdiLepton)
                            {
                                histo1D["h_cutFlow"]->Fill(10., scaleFactor*lumiWeight);
                                //histo1D["h_2SSL_1st_Elec_Pt"]->Fill(tempLepton_0.Pt(),scaleFactor*lumiWeight);
                                //histo1D["h_2SSL_2nd_Elec_Pt"]->Fill(tempLepton_1.Pt(),scaleFactor*lumiWeight);
                                //histo2D["h_2SSL_2j_1b_Nb_bjets_vs_Nb_jets"]->Fill(selectedBCSVLJets.size(), selectedJets.size(),scaleFactor*lumiWeight);
                                for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                {
                                    sum_jet_PT+= selectedJets[ijet]->Pt();
                                }
                                Ht_AllJets=sum_jet_PT;
                                St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                                histo1D["h_2SSL_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                                histo1D["h_2SSL_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                                SSLeptonTree->Fill();
                            }
                            if (OSdiLepton && !SSdiLepton)
                            {
                                histo1D["h_cutFlow"]->Fill(11., scaleFactor*lumiWeight);
                                //histo1D["h_2OSL_1st_Elec_Pt"]->Fill(tempLepton_0.Pt(),scaleFactor*lumiWeight);
                                //histo1D["h_2OSL_2nd_Elec_Pt"]->Fill(tempLepton_1.Pt(),scaleFactor*lumiWeight);
                                //histo2D["h_2OSL_2j_1b_Nb_bjets_vs_Nb_jets"]->Fill(selectedBCSVLJets.size(), selectedJets.size(),scaleFactor*lumiWeight);
                                for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                {
                                    sum_jet_PT+= selectedJets[ijet]->Pt();
                                }
                                Ht_AllJets=sum_jet_PT;
                                St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                                histo1D["h_2OSL_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                                histo1D["h_2OSL_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                                
                                OSLeptonTree->Fill();
                            }
                            eventSelected = true;
                            
                            
                        }//2L+>=2Jets+>=1CSVLB
                        
                        //                                if (selectedJets.size()>= 3 && selectedBCSVLJets.size()>= 1) {
                        //                                    histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
                        //                                }
                        //                                if (selectedJets.size()>= 4 && selectedBCSVLJets.size()>= 1) {
                        //                                    histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                        //                                }
                        
                    }//2L+>=2Jets
                    
                }//2L
                
            }//2L
            
            // }//GoodPV
            // }//trigged
            
            
            
            
            
//            if (! eventSelected )
//            {
//                continue;
//            }
            
            if(verbose>3)cout << "filling the tree" << endl;
            //myTree->Fill();
            nofSelectedEvents++;
            
            
            if (verbose > 2)
                cout << "  Event " << ievt << " is selected" << endl;
        } // end the loop over events
        sumW = (int) sumWeights;
        nEv = (int) nEvents;
        
        globalTree->Fill();
        
        //////////////////////
        ///  END OF EVENT  ///
        //////////////////////
        
        cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
        if (! isData  )
        {
            cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
            cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
            cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents) << endl;
            
            // Determine scale factor due to negative weights
            nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
            cout << "This corresponds to an event scale factor of " << nloSF  << endl;
        }
        myTree->Write();
        fileout->cd();
        fileout->Write();
        fileout->Close();
        delete fileout;  //   delete myTree;
        
        ///*****************///
        ///   Write plots   ///
        ///*****************///
        
        ///////////////////
        /// CLEANING
        /////////////////
        //delete fileout;  //   delete myTree;
        
        if(!isData && !btagShape) delete btwt;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }///end the loop over the datasets
    
//    string pathPNG = "OutPutHistos/";
//    mkdir(pathPNG.c_str(),0777);
//    mkdir((pathPNG+"1DPlot/").c_str(),0777); // 0777 if it doesn't exist already, make it
    fout->cd();
    for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
        if(verbose>3) cout << "1D Plot: " << it->first << endl;
        TCanvas *ctemp = new TCanvas();
        ctemp->cd();
        TH1F *temp = it->second;
        string name = it->first;
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
    
    //        for (map<string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    //        {
    //            if(verbose>1) cout << "2D Plot: " << it->first << endl;
    //            TH2F *temp = it->second;
    //            string name = it->first;
    //            temp->Draw();
    //        }
    //
    fout->Write();
    fout->Close();

    delete fout;
    
    
   // delete tcAnaEnv;
    
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;

   
}