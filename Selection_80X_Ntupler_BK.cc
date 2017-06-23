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

struct HighestCSVBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

////// user defined functions
void ElectronChargeMisId(vector<TRootElectron*> SelElectron);

int main (int argc, char *argv[])
{
    clock_t start = clock();
    
    
    //Initializing b-tag WP ref https://indico.cern.ch/event/535758/contributions/2177471/attachments/1279035/1899097/BTagPerf_160525_80XWPs.pdf
    

    //// CSVv2 tagger // updated according to latest recommendation https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    float CSVv2_workingpointvalue_Loose = 0.5426;//working points updated to 2016 BTV-POG recommendations.
    float CSVv2_workingpointvalue_Medium = 0.8484;//working points updated to 2016 BTV-POG recommendations.
    float CSVv2_workingpointvalue_Tight = 0.9535;//working points updated to 2016 BTV-POG recommendations.
    //// supercombined tagger
    float cMVA_workingpointvalue_Loose = -0.5884;//working points updated to 2016 BTV-POG recommendations.
    float cMVA_workingpointvalue_Medium = 0.4432;//working points updated to 2016 BTV-POG recommendations.
    float cMVA_workingpointvalue_Tight = 0.9432;//working points updated to 2016 BTV-POG recommendations.
    
    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FCNC_2SSL search ! "           << endl;
    cout << "*************************************************************" << endl;
    
    ///// ** Define the name of Directory that contains output *** ////
    
    string runDate = "Test_DDChargeMisIdCorrected_13June17";
    
    /////////////////////////////////
    //// *** Working Conditions //// have to be Checked everytime before running the code !!!
    ////////////////////////////////
    bool Elec_Elec, Mu_Mu, Elec_Mu, Mu_Elec, Apply_HLT_Triggers, eventSelected, Fake_Electrons, ApplyCharge_misID, ApplyElec_SF , ApplyMu_SF , ApplyPU_SF, Apply_btag_SF, Apply_JetCleaning, trigged,debug, printTrigger, All_lep, RemoveBadMuon;
    
    Elec_Elec = true;
    ApplyElec_SF = true;
    ApplyCharge_misID = true;
    Mu_Mu =false;
    ApplyMu_SF = false;
    RemoveBadMuon = true;
    Elec_Mu = false;
    Mu_Elec = false;
    All_lep = false;
    Apply_HLT_Triggers = true;
    printTrigger = false;
    eventSelected= false;
    Fake_Electrons = false;
    ApplyPU_SF = true;
    Apply_JetCleaning = true;
    trigged = false;
    debug = false;
    bool applyJER = true;
    bool applyJES = true;
    bool applyNegWeightCorrection = true;
    bool btagShape = true;
    bool fillBtagHisto = false;
   // bool CSVv2nonshape = false;//Put to false if you don't want to use the regular CSVv2 SF's
    std::string channelpostfix = "";
   
    //// ** use this for run JER/JES systematics ////
    int doJESShift = 0; // 0: off 1: minus 2: plus
    cout << "doJESShift: " << doJESShift << endl;
    
    int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
    cout << "doJERShift: " << doJERShift << endl;
    
    int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
    cout << "dobTagEffShift: " << dobTagEffShift << endl;
    
    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;
    
   // string postfix = "_Run2_TopTree_Study_" + dName; // to relabel the names of the output file
    
//    if (doJESShift == 1)
//    {
//        postfix += "JESMinus";
//    }
//    if (doJESShift == 2)
//    {
//        postfix += "JESPlus";
//    }
//    if (doJERShift == 1)
//    {
//        postfix += "JERMinus";
//    }
//    if (doJERShift == 2)
//    {
//        postfix += "JERPlus";
//    }
    
    /////////////////////
    ///  Configuration
    /////////////////////
    
    /// xml file
    string xmlFileName ="";
    string Channel = "";
    if (argc > 1) xmlFileName = (string)argv[1];
    const char *xmlfile = xmlFileName.c_str();
    double dataLumi = 0;
    float lum_RunsBCDEF = 19.68;// /fb
    float lum_RunsGH = 16.146;// /fb
    
    //Setting Lepton Channels
    if(Elec_Elec)
    {
        cout << " --> Using the Electron-Electron channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_ElEl_V10_Samples.xml";
        dataLumi = 35826.306833965; //Runs from B to H
        channelpostfix = "_ElEl_";
        Channel = "Dilepton_ElecElec";
    }
    else if(Mu_Mu)
    {
        cout << " --> Using the Muon-Muon channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_MuMu_V10_Samples.xml";
        dataLumi = 35826.306833965; //Run B+C+D+E+F+G+H
        channelpostfix = "_MuMu_";
        Channel = "Dilepton_MuMu";
    }
    else if(Elec_Mu)
    {
        cout << " --> Using the Electron-Muon channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_ElMu_V10_Samples.xml";
        dataLumi = 35826.306833965;
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
    // if there only two arguments after the fileName, the jobNum will be set to 0 by default as an integer is expected and it will get a string (lastfile of the list)

    const int JobNum                = strtol(argv[argc-3], NULL, 10);
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);
    
    vector<string> vecfileNames;
    for(int args = 11; args < argc-3; args++)
    {
        vecfileNames.push_back(argv[args]);
    }
    
    
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();
    
    ofstream infoFile;
    
    string info_dir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_26_patch1/src/TopBrussels/FCNCAnalysis/Information/"+Channel +"/";
    mkdir(info_dir.c_str(),0777);
    string info_date_dir = info_dir + runDate +"/";
    cout << "info dir " << info_dir.c_str() << endl;
    mkdir(info_date_dir.c_str(),0777);
    string infoName = info_date_dir + "information";
    infoName += "_"+ Channel;
    infoName += "_" + dName;
    infoName += "_" + strJobNum;
    infoName += ".txt";
    infoFile.open(infoName.c_str());
    infoFile.precision(3);

    
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
    
    bool isData = false;
    // Print information to a textfile
    if(dName.find("Data")!=string::npos || dName.find("data")!=string::npos || dName.find("DATA")!=string::npos){
        isData = true;
        cout << "running on data !!!!" << endl;
        cout << "luminosity is " << dataLumi << endl;
    }
    ////////////////////////////////////
    ///  AnalysisEnvironment
    ////////////////////////////////////
    
    //TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    
    AnalysisEnvironment anaEnv;
    cout << " - Loading environment ..." << endl;
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    if (!isData) anaEnv.METCollection = "PFMET_slimmedMETs";
    if(isData) anaEnv.METCollection = "PFMET_slimmedMETsMuEGClean";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    ////anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    //////anaEnv.ElectronCollection = "Electrons_calibratedPatElectrons"; ///     used from tag 80X_v2 to version v4
    anaEnv.ElectronCollection = "Electrons_selectedElectrons"; /// used from v7 where reminiaod samples used
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
    
    
    /////////////////////////////////////////////
    ///// --- Define Scaling Factors ----- /////
    ////////////////////////////////////////////
  
    string pathToCaliDir = "../TopTreeAnalysisBase/Calibrations/";
    
    ////***PU SF + Systematicss added ***/////
    
    LumiReWeighting LumiWeights("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
    
    LumiReWeighting LumiWeights_UP("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysPlus.root", "pileup", "pileup");
    
    LumiReWeighting LumiWeights_Down("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet_sysMinus.root", "pileup", "pileup");

    ///// *** lepton scaling factors *** ////
    
    //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)
    
    MuonSFWeight *muonSFWeightID;
    MuonSFWeight *muonSFWeightIso;
    MuonSFWeight *muonSFWeightID_BCDEF;
    MuonSFWeight *muonSFWeightIso_BCDEF;
    MuonSFWeight *muonSFWeightID_GH;
    MuonSFWeight *muonSFWeightIso_GH;
    
    TFile *muontrackfile = new TFile("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonSF/20170413/Tracking_EfficienciesAndSF_BCDEFGH.root","read");
    TGraph* h_muonSFWeightTrack = (TGraph*) muontrackfile->Get("ratio_eff_eta3_dr030e030_corr")->Clone();//Tracking efficiency as function of eta

    
    if (ApplyMu_SF && (Mu_Mu || Elec_Mu))
    {
        muonSFWeightID_BCDEF = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/20170413/"+"IDEfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
        muonSFWeightIso_BCDEF = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/20170413/"+"IsoEfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
        
        muonSFWeightID_GH = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/20170413/"+"IDEfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
        muonSFWeightIso_GH = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/20170413/"+"IsoEfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
        
        
    }


    ////Triggers SF for muons to be added
    

    string electronFile= "egammaEffi.txt_EGM2D_IDcutbTight_20170413.root";
    string electronRecoFile = "egammaEffi.txt_EGM2D_reco_20170413.root";
    string elecHistName = "EGamma_SF2D";
    
    //ElectronSFWeight::ElectronSFWeight(std::string const&, std::string const&, bool const&, bool const&, bool const&)
    
    ElectronSFWeight* electronSFWeightID = new ElectronSFWeight (pathToCaliDir+"LeptonSF/ElectronSF/20170413/"+electronFile,elecHistName, true,false, false); // (... , ... , debug, print warning)
    ElectronSFWeight* electronSFWeightReco = new ElectronSFWeight(pathToCaliDir+"LeptonSF/ElectronSF/20170413/"+electronRecoFile,elecHistName, true,false, false,true);

    
    ///-- dfeine Triggering ---//
    // /// HLT Triggers will used in the analysis is according to Top Trigger (Run2)
    // //// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
    
    if(verbose == 0) cout << "Initializing trigger" << endl;
    //Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep);
   
    Trigger* DiElecTrigger = new Trigger(0,1,0,1);
    Trigger* DiMuTrigger = new Trigger(1,0,0,1);
    Trigger* DiElMuTrigger = new Trigger(1,1,0,1);
    
    // JER / JEC
    vector<JetCorrectorParameters> vCorrParam;
    
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
    
    float oldLuminosity = dataLumi ;  //anaEnv.Luminosity;  // in 1/pb
    cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
    
    float Luminosity=oldLuminosity;
    
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
        string dataSetName = datasets[d]->Name();
        //if(dataSetName.find("Data")==0 || dataSetName.find("data")==0 || dataSetName.find("DATA")==0)
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            Luminosity = datasets[d]->EquivalentLumi();
            cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        }
    }

    if (Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
    
    
//    stringstream ss;
//    ss << JobNum;
//    string strJobNum = ss.str();
    
    //Global variable
    TRootEvent* event = 0;
    
    //nof selected events
    int nofSelectedEvents = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    
    /////////////////////////////
    /// Object ID              ///
    /////////////////////////////
    // electron
    float el_pt_cut =15.;
    float el_eta_cut = 2.5;
    float el_iso_cone  = 0.3;
    //// reliso cut fabs(eta supercluster) <= 1.479 --> 0.107587 // (fabs(eta supercluster) > 1.479 && fabs(eta supercluster) < 2.5) --> 0.113254
    // muon
    float mu_pt_cut = 12.;
    float mu_eta_cut = 2.1;
    float mu_iso_cut = 0.15;
    //jets
    float jet_pt_cut = 25.;
    float jet_eta_cut = 2.4;
    
    
    ///////////////////////\\\\\\\\\\\\\\\\\\\\\\\
    ///// Create root file contains histograms \\\\
    /////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\

    string histoDir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_26_patch1/src/TopBrussels/FCNCAnalysis/Output_Histos/";
    mkdir(histoDir.c_str(),0777);
    string histoPath = histoDir + Channel+"_";
    histoPath += runDate+"/";
    mkdir(histoPath.c_str(),0777); // create the directory histoPath if it is not exist already
    string histoPathSampleName = histoPath+dName;
    mkdir((histoPathSampleName+"/").c_str(),0777);
    string histoFileName = histoPathSampleName+"/"+channelpostfix+dName+".root";
    if (strJobNum != "0")
    {
        if(verbose == 0) cout << "strJobNum is " << strJobNum << endl;
        histoFileName = histoPathSampleName+"/"+"FCNC_2SSL_"+dName+"_"+strJobNum+".root";
    }
    TFile *fout = new TFile(histoFileName.c_str(),"RECREATE");
    
    // map <string,TH1F*> histo1D;
    std::string titlePlot = "";
    
    titlePlot = "cutFlow"+channelpostfix;
    histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 16,-0.5,15.5);
    
    titlePlot = "Nb_Events_Lumi"+channelpostfix;
    histo1D["h_Nb_Events_Lumi"] = new TH1F(titlePlot.c_str(), "Nb_Events_Lumi",  6, - 0.5, 5.5 );

    
    titlePlot = "cutFlow_LeptonCuts"+channelpostfix;
    histo1D["h_cutFlow_LeptonCuts"] = new TH1F(titlePlot.c_str(), "#events after Applied Cuts ", 16,-0.5,15.5);
    
    titlePlot = "cutFlow_ChargeMisId"+channelpostfix;
    histo1D["h_cutFlow_ChargeMisId"] = new TH1F(titlePlot.c_str(), "#events before & after Charge mis-id ", 16,-0.5,15.5);
    
    titlePlot = "cutFlow_DD_ChargeMisId"+channelpostfix;
    histo1D["h_cutFlow_DD_ChargeMisId"] = new TH1F(titlePlot.c_str(), "#events DD before & after Charge mis-id ", 16,-0.5,15.5);

    
    titlePlot = "NSSL_barrel_cutFlow"+channelpostfix;
    histo1D["h_NSSL_barrel_cutFlow"] = new TH1F(titlePlot.c_str(), "NSSL_barrel_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NOSL_barrel_cutFlow"+channelpostfix;
    histo1D["h_NOSL_barrel_cutFlow"] = new TH1F(titlePlot.c_str(), "NOSL_barrel_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NSSL_EndCap_cutFlow"+channelpostfix;
    histo1D["h_NSSL_EndCap_cutFlow"] = new TH1F(titlePlot.c_str(), "NSSL_EndCap_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NOSL_EndCap_cutFlow"+channelpostfix;
    histo1D["h_NOSL_EndCap_cutFlow"] = new TH1F(titlePlot.c_str(), "NOSL_EndCap_cutflow", 6,-0.5,5.5);
    
    titlePlot = "Nb_EventsWith3ChargeAgr"+channelpostfix;
    histo1D["h_EventsWith3ChargeAgr"] = new TH1F(titlePlot.c_str(), "#Events After Selective method Applied", 6,-0.5,5.5);
    
    titlePlot = "Nb_EventsBefore3ChargeAgr"+channelpostfix;
    histo1D["h_EventsBefore3ChargeAgr"] = new TH1F(titlePlot.c_str(), "#Events After Selective method Applied", 6,-0.5,5.5);
    
    titlePlot = "SelctiveElecCharge_CutFlow"+channelpostfix;
    histo1D["h_SelctiveElecCharge_CutFlow"] = new TH1F(titlePlot.c_str(), "#Events CutFlow for Selctive Electron Charge ", 6,-0.5,5.5);
    
    titlePlot = "ElecChargeFlips_2atBarrel_CutFlow"+channelpostfix;
    histo1D["h_ElecChargeFlips_2atBarrel_CutFlow"] = new TH1F(titlePlot.c_str(), "#Events CutFlow for Electron Charge flips barrel ", 6,-0.5,5.5);
    
    titlePlot = "ElecChargeFlips_2atEndCap_CutFlow"+channelpostfix;
    histo1D["h_ElecChargeFlips_2atEndCap_CutFlow"] = new TH1F(titlePlot.c_str(), "#Events CutFlow for Electron Charge flips EndCap ", 6,-0.5,5.5);
    
    titlePlot = "ElecChargeFlips_1atBarrelEndCap_CutFlow"+channelpostfix;
    histo1D["h_ElecChargeFlips_1atBarrelEndCap_CutFlow"] = new TH1F(titlePlot.c_str(), "#Events CutFlow for Electron Charge flips 1atBarrelEndCap ", 6,-0.5,5.5);
    
    
//    titlePlot = "IsSelectiveCharge"+channelpostfix;
//    histo1D["h_IsSelectiveCharge"] = new TH1F(titlePlot.c_str(), "NOSL_EndCap_cutflow", 3,-0.5,2.5);
//    
//    titlePlot = "IsChargeFilpped"+channelpostfix;
//    histo1D["h_IsChargeFilpped"] = new TH1F(titlePlot.c_str(), "NOSL_EndCap_cutflow", 3,-0.5,2.5);
    
    
    
    
    //*** histos for Jets *** //
    
    titlePlot = "initial_Nb_Jets"+channelpostfix;
    histo1D["h_initial_Nb_Jets"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets",  16, - 0.5, 15.5 );
    

    
    /////////////////////////////////
    //////*** 2D Histograms ***/////
    ///////////////////////////////
    titlePlot = "2L_Nb_jets_vs_CSVLbjets"+channelpostfix;
    histo2D["2L_Nb_jets_vs_CSVLbjets"] = new TH2F(titlePlot.c_str(),"After 2L #Jets Vs #CSVLbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_Nb_jets_vs_CSVMbjets"+channelpostfix;
    histo2D["2L_Nb_jets_vs_CSVMbjets"] = new TH2F(titlePlot.c_str(),"After 2L #Jets Vs #CSVMbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_Nb_jets_vs_CSVTbjets"+channelpostfix;
    histo2D["2L_Nb_jets_vs_CSVTbjets"] = new TH2F(titlePlot.c_str(),"After 2L #Jets Vs #CSVTbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_3Jets_Nb_jets_vs_CSVLbjets"+channelpostfix;
    histo2D["2L_3Jets_Nb_jets_vs_CSVLbjets"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets #Jets Vs #CSVLbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_3Jets_Nb_jets_vs_CSVMbjets"+channelpostfix;
    histo2D["2L_3Jets_Nb_jets_vs_CSVMbjets"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets #Jets Vs #CSVMbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_3Jets_Nb_jets_vs_CSVTbjets"+channelpostfix;
    histo2D["2L_3Jets_Nb_jets_vs_CSVTbjets"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets #Jets Vs #CSVTbJets",15,-0.5,14.5, 15, -0.5,14.5);
    
    titlePlot = "2L_Lep0_vs_Lep1_Pt"+channelpostfix;
    histo2D["h_2L_Lep0_vs_Lep1Pt"] = new TH2F(titlePlot.c_str(),"After 2L Lep0 P_{T} Vs Lep1 P_{T}",100,0,500, 100, 0,500);
    
    titlePlot = "2L_3Jets_Lep0_vs_Lep1_Pt"+channelpostfix;
    histo2D["h_2L_3Jets_Lep0_vs_Lep1Pt"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets Lep0 P_{T} Vs Lep1 P_{T}",100,0,500, 100, 0,500);
    
    titlePlot = "2L_3Jets1b_Lep0_vs_Lep1_Pt"+channelpostfix;
    histo2D["h_2L_3Jets1b_Lep0_vs_Lep1Pt"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets1b Lep0 P_{T} Vs Lep1 P_{T}",100,0,500, 100, 0,500);

    titlePlot = "2L_3Jets1b_Ht_vs_metPt"+channelpostfix;
    histo2D["h_2L_3Jets1b_Ht_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2L+3Jets1b H_{T} Vs met P_{T}",100,0,500, 100, 0,500);
    
    titlePlot = "2SSL_3Jets1b_Ht_vs_metPt"+channelpostfix;
    histo2D["h_2SSL_3Jets1b_Ht_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2SSL+3Jets1b H_{T} Vs met P_{T}",100,0,500, 100, 0,500);
    
    titlePlot = "2OSL_3Jets1b_Ht_vs_metPt"+channelpostfix;
    histo2D["h_2OSL_3Jets1b_Ht_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2OSL+3Jets1b H_{T} Vs met P_{T}",100,0,500, 100, 0,500);
    
    titlePlot = "2L_2lDeltaPhi_vs_metPt"+channelpostfix;
    histo2D["h_2L_2lDeltaPhi_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2L #Delta #Phi  Vs E_{T}^{mis} ",100,0,500, 100, 0,500);
    
    titlePlot = "2SSL_2lDeltaPhi_vs_metPt"+channelpostfix;
    histo2D["h_2SSL_2lDeltaPhi_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2SSL #Delta #Phi  Vs E_{T}^{mis} ",100,0,500, 100, 0,500);
    
    titlePlot = "2OSL_2lDeltaPhi_vs_metPt"+channelpostfix;
    histo2D["h_2OSL_2lDeltaPhi_vs_metPt"] = new TH2F(titlePlot.c_str(),"After 2SSL #Delta #Phi  Vs E_{T}^{mis} ",100,0,500, 100, 0,500);
    
    
    titlePlot = "NSSL_Pt_eta"+channelpostfix;
    histo2D["h_NSSL_Pt_eta"] = new TH2F(titlePlot.c_str()," NSSL P_{T} vs #eta  ",6,-0.5,5.5 , 6,-0.5,5.5);
    
    titlePlot = "NOSL_Pt_eta"+channelpostfix;
    histo2D["h_NOSL_Pt_eta"] = new TH2F(titlePlot.c_str()," NOSL P_{T} vs #eta  ",6,-0.5,5.5 , 6,-0.5,5.5);
    
    titlePlot = "DD_ChargemisId_Pt_eta"+channelpostfix;
    histo2D["h_DD_ChargemisId_Pt_eta"] = new TH2F(titlePlot.c_str()," DD charge mis-Id P_{T} vs #eta  ",6,-0.5,5.5 , 6,-0.5,5.5);

    
    ////////////////////////////////////
    ///  Loop on datasets
    ////////////////////////////////////
    
    if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
    bool nlo = false;
    bool isSignal = false;
    bool isDY = false;
    
    ///////----Start looping over datasets -----/////

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
        string previousFilename = "";
        int iFile = -1;
       
        string dataSetName = datasets[d]->Name();
        if (verbose > 0)cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
        if (verbose > 0)cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
        double nloSF = 1;
        double sumWeights = 0;
        nlo = false;
        double lumiWeight = -99.;
        //if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        if (dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
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
            if(dataSetName.find("FCNC")!=string::npos) isSignal = true;
            if(dataSetName.find("DYJets")!=string::npos) isDY = true;
            
        }
        
        if (verbose >1)cout << "the lumiweight of the dataset : "<< datasets[d]->Name() << " is " << lumiWeight << endl;
        //SampleLumiWeight = lumiWeight;
        
        //// --- Calibration - Applying Jet correction (Jec) --- ////
        
        vCorrParam.clear();
        JetCorrectionUncertainty *jecUnc;
        
        if(dName.find("DataRun2016B")!=string::npos || dName.find("DataRun2016C")!=string::npos || dName.find("DataRun2016D")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("DataRun2016E")!=string::npos || dName.find("DataRun2016F")!=string::npos )
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("DataRun2016G")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt");
        }
        else if(dName.find("DataRun2016H_v2")!=string::npos || dName.find("DataRun2016H_v3")!=string::npos)
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            isData = true;
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_DATA/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt");
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt");
        }
        
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
        
        ///////*** btag reweighing **** ///////
        ///// b-tagging scaling factor
        int mkdirstatus_btag = mkdir("BTagHistosPtEta",0777);
//        BTagCalibration * bTagCalib_CSVv2;
//        BTagCalibrationReader * b_reader_CSVv2;
//        BTagCalibrationReader *reader_csvv2;
//
        //    ///// b-tagging scaling factor
        //
        BTagCalibration * bTagCalib_CSVv2;
        BTagCalibrationReader * b_reader_CSVv2;
        BTagCalibrationReader *reader_csvv2;

        BTagWeightTools *btwt_CSVv2L_mujets_central = 0;
        BTagCalibrationReader * b_reader_CSVv2_shape;
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
        
        
        if(!isData)
        {
            // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
            //	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv");
            
            bTagCalib_CSVv2 = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
            
            if(!btagShape)b_reader_CSVv2 = new BTagCalibrationReader(bTagCalib_CSVv2,// calibration instance
                                                       BTagEntry::OP_LOOSE, // operating point
                                                       "mujets",// measurement type
                                                       "central"); // systematics type  --> depending on JES up/Down andother reader is needed
            
            b_reader_CSVv2_shape = new BTagCalibrationReader(bTagCalib_CSVv2,BTagEntry::OP_RESHAPING,"iterativefit","central"); //reshaping
            
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
            
            
            if(fillBtagHisto && !btagShape)  // before btag reweighting can be applied, you first have to make the histograms
            {
                btwt_CSVv2L_mujets_central = new BTagWeightTools(b_reader_CSVv2,"BTagHistosPtEta/HistosPtEta_"+dName+ "_" + strJobNum +"_mujets_central.root",false,30,999,2.4);
                //btwt_CSVv2L_mujets_central = new BTagWeightTools(b_reader_CSVv2,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",false,30,999,2.4);
            }
            else if(!btagShape)
            {
                btwt_CSVv2L_mujets_central = new BTagWeightTools(b_reader_CSVv2,"BTagHistosPtEta/HistosPtEta_"+dName+"_mujets_central.root",true,30,999,2.4);
                // btwt_CSVv2L_mujets_central = new BTagWeightTools(b_reader_CSVv2,"BTagHistosPtEta/HistosPtEta_TTJets_mujets_central.root",false,30,999,2.4);
            }
            
            
        }
        
        //// ***************** /////
        /// output Ntuples /////
        /// ***************** ////
        string rootTreesDir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_26_patch1/src/TopBrussels/FCNCAnalysis/Output_Ntuples/";
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
        /// Decleration of the Trees
        
        TTree *bookkeeping = new TTree("startevents","startevents");
        TTree* InitialTree = new TTree("Initialtree","Initialtree");
        TTree* NoZmassVetoTree = new TTree("NoZmassVetoTree","NoZmassVetoTree");
        TTree* SSLNoZmassVetoTree = new TTree("SSLNoZmassVetoTree","SSLNoZmassVetoTree");
        TTree* OSLNoZmassVetoTree = new TTree("OSLNoZmassVetoTree","OSLNoZmassVetoTree");
        TTree* myTree = new TTree("tree","tree");
        TTree* SSLeptonTree = new TTree("SSLeptonTree","SSLeptonTree");
        
        TTree* OSLeptonTree = new TTree("OSLeptonTree","OSLeptonTree");
        TTree* globalTree = new TTree("globaltree","globaltree");
        TTree* MVATree = new TTree("MVATree","MVATree");

        
        //////////////////////////////
        // My tree - variables //
        //////////////////////////////
        
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        Int_t nEv;
        Int_t nEv_lumi;
        Double_t Count_cut[12];
        Int_t nCuts = 11;
        // various weights
        
        Int_t nofPosWeights;
        Int_t nofNegWeights;
        Int_t sumW;
        Int_t JERon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Int_t JESon = 0; // -1: not on, 0: nominal, 1: minus, 2: plus
        Double_t nloWeight; // for amc@nlo samples
        
        Double_t sf_muon[10];
        
        Double_t puSF;
        Double_t puSF_UP;
        Double_t puSF_Down;
        Double_t btagSF;
        Float_t SampleLumiWeight;
        
        Double_t btagSFshape = 1.;
        Double_t btagSFshape_down_cferr1 = 1.;
        Double_t btagSFshape_down_cferr2 = 1.;
        Double_t btagSFshape_down_hf= 1.;
        Double_t btagSFshape_down_hfstats1 = 1.;
        Double_t btagSFshape_down_hfstats2 = 1.;
        Double_t btagSFshape_down_lf = 1.;
        Double_t btagSFshape_down_lfstats1 = 1.;
        Double_t btagSFshape_down_lfstats2 = 1.;
        
        Double_t btagSFshape_up_cferr1 = 1.;
        Double_t btagSFshape_up_cferr2 = 1.;
        Double_t btagSFshape_up_hf= 1.;
        Double_t btagSFshape_up_hfstats1 = 1.;
        Double_t btagSFshape_up_hfstats2 = 1.;
        Double_t btagSFshape_up_lf = 1.;
        Double_t btagSFshape_up_lfstats1 = 1.;
        Double_t btagSFshape_up_lfstats2 = 1.;
        
        
        //// Variables for met
        Int_t PassedMETFilter;
        Float_t met_Pt;
        Float_t met_Eta;
        Float_t met_Phi;
        Float_t corrected_met_Eta ;
        Float_t corrected_met_Phi ;
        Float_t UncorrJES_met_Pt ;
        Float_t corrJES_met_Pt ;
        Float_t MT_lep1;
        Float_t MT_lep2;
        Float_t MT_Elec1;
        Float_t MT_Elec2;
        Float_t MT_Mu1;
        Float_t MT_Mu2;
        Float_t DeltaPhi_met_lep1;
        Float_t DeltaPhi_met_lep2;
        
        
        
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
        Bool_t isSelectiveElec[10];
        Float_t pt_1st_Electron;
        Float_t pt_2nd_Electron;
        Float_t eta_1st_Electron;
        Float_t eta_2nd_Electron;
        Float_t phi_1st_Electron;
        Float_t phi_2nd_Electron;
        Float_t E_1st_Electron;
        Float_t E_2nd_Electron;
        Float_t charge_1st_Electron;
        Float_t charge_2nd_Electron;
        
        
        Double_t ElectronSF[10];
        Double_t ElectronSF_up[10];
        Double_t ElectronSF_down[10];
        Double_t ElectronSFID[10];
        Double_t ElectronSFID_up[10];
        Double_t ElectronSFID_down[10];
        Double_t ElectronSFReco[10];
        Double_t ElectronSFReco_up[10];
        Double_t ElectronSFReco_down[10];
        
        //Kinamatic variables for Muons
        Int_t nMuons;
        Double_t MuonIDSF[10];
        Double_t MuonIsoSF[10];
        Double_t MuonTrackSF[10];
        Double_t MuonTrackSF_up[10];
        Double_t MuonTrackSF_down[10];
        Double_t MuonIDSF_up[10];
        Double_t MuonIsoSF_up[10];
        Double_t MuonIDSF_down[10];
        Double_t MuonIsoSF_down[10];
        
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
        
        
        Int_t charge_muon[10];
        Float_t pt_1st_Muon;
        Float_t pt_2nd_Muon;
        Float_t eta_1st_Muon;
        Float_t eta_2nd_Muon;
        Float_t phi_1st_Muon;
        Float_t phi_2nd_Muon;
        Float_t E_1st_Muon;
        Float_t E_2nd_Muon;
        Float_t charge_1st_Muon;
        Float_t charge_2nd_Muon;
        
        
        //Kinamatic variables for jets
        Int_t nJets;
        Int_t nJetpair;
        Int_t nCSVLbJets;
        Int_t nCSVMbJets;
        Int_t nCSVTbJets;
        Int_t nNonCSVLbJets;
        Int_t nNonCSVMbJets;
        Int_t nNonCSVTbJets;
        Float_t bdiscCSVv2_2nd_bjet;
        Float_t bdiscCSVv2_1st_bjet;
        Float_t bdiscCSVv2_2nd_jet;
        Float_t bdiscCSVv2_1st_jet;
        Float_t bdiscCSVv2_3rd_jet;
        Float_t bdiscCSVv2_4th_jet;
        
        Double_t pt_jet[20];
        Double_t phi_jet[20];
        Double_t eta_jet[20];
        Double_t E_jet[20];
        Double_t charge_jet[20];
        Double_t bdisc_jet[20];
        Double_t Mass_JetPair_arr[20];
        
        Float_t pt_1st_jet;
        Float_t pt_2nd_jet;
        Float_t pt_3rd_jet;
        Float_t pt_4th_jet;
        Float_t eta_1st_jet;
        Float_t eta_2nd_jet;
        Float_t eta_3rd_jet;
        Float_t eta_4th_jet;
        Float_t phi_1st_jet;
        Float_t phi_2nd_jet;
        Float_t phi_3rd_jet;
        Float_t phi_4th_jet;
        Float_t E_1st_jet;
        Float_t E_2nd_jet;
        Float_t E_3rd_jet;
        Float_t E_4th_jet;
        
        
        ///// --- Variables used for MVA --- //
        Float_t Ht;
        Float_t St;
        Float_t DeltaR_2L;
        Float_t DeltaPhi_2L;
        Float_t invMass_2L;
        Float_t invMass_2SSL_Zmass;
        Float_t invMass_2OSL_Zmass;
        Float_t DeltaR_Mu0b0_DiMu;
        Float_t DeltaR_Mu1b0_DiMu;
        Float_t DeltaR_Elec0b0_DiElec;
        Float_t DeltaR_Elec1b0_DiElec;
        Float_t DeltaR_Mu0b0_DiMuEl;
        Float_t DeltaR_Mu1b0_DiElMu;
        Float_t DeltaR_Elec1b0_DiMuEl;
        Float_t DeltaR_Elec0b0_DiElMu;
        Float_t Mass_WJetPair;
        Float_t Mass_JetPair;
        Float_t MVA_nJets;
        Float_t MVA_nCSVLbtagJets;
        
        /*
        Float_t Match_MC_Ht;
        Float_t Match_MC_St;
        Float_t Match_MC_DeltaR_2L;
        Float_t Match_MC_DeltaPhi_2L;
        Float_t Match_MC_invMass_2L;
        Float_t Match_SM_MC_DeltaR_Mu0b0_DiMu;
        Float_t Match_MC_DeltaR_Mu1b0_DiMu;
        Float_t Match_SM_MC_DeltaR_Elec0b0_DiElec;
        Float_t Match_MC_DeltaR_Elec1b0_DiElec;
        Float_t Match_SM_MC_DeltaR_Mu0b0_DiMuEl;
        Float_t Match_MC_DeltaR_Mu1b0_DiElMu;
        Float_t Match_MC_DeltaR_Elec1b0_DiMuEl;
        Float_t Match_SM_MC_DeltaR_Elec0b0_DiElMu;
        Float_t Match_FCNC_MC_Mass_WJetPair;
        Float_t Match_FCNC_MC_Mass_JetPair;
        */
        ///// ---  variables used in MC ---- ////
        Double_t Jet_matchMC_pdgId[20];
        Double_t Jet_matchMC_motherpdgId[20];
        Double_t Jet_matchMC_grannypdgId[20];
        
    
        
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
        globalTree->Branch("nEv_lumi" , &nEv_lumi, "nEv_lumi/I");
        globalTree->Branch("sumW", &sumW, "sumW/I");
        globalTree->Branch("nCuts",&nCuts, "nCuts/I");
        globalTree->Branch("Count_cut",&Count_cut,"Count_cut[nCuts]/D");
        
        // define the output trees may I need to make many trees depend on selection cuts
        /////(Integer variables)
        ////** Tree before applying Z_Mass Veto ** ///
        
        NoZmassVetoTree->Branch("isData",&isData,"isData/I");
        NoZmassVetoTree->Branch("run_num",&run_num,"run_num/I");
        NoZmassVetoTree->Branch("evt_num",&evt_num,"evt_num/I");
        NoZmassVetoTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        NoZmassVetoTree->Branch("nvtx",&nvtx,"nvtx/I");
        NoZmassVetoTree->Branch("npu",&npu,"npu/I");
        NoZmassVetoTree->Branch("puSF",&puSF,"puSF/D");
        NoZmassVetoTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        NoZmassVetoTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        
        NoZmassVetoTree->Branch("btagSF",&btagSF,"btagSF/D");
        NoZmassVetoTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        NoZmassVetoTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        NoZmassVetoTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        NoZmassVetoTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        NoZmassVetoTree->Branch("sumW", &sumW, "sumW/I");
        NoZmassVetoTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
       // NoZmassVetoTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
        NoZmassVetoTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("ElectronSF",ElectronSF,"ElectronSF[nElectrons]/D");
        
        NoZmassVetoTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        NoZmassVetoTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        NoZmassVetoTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        NoZmassVetoTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        NoZmassVetoTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        NoZmassVetoTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        NoZmassVetoTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        NoZmassVetoTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        NoZmassVetoTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        NoZmassVetoTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        NoZmassVetoTree->Branch("nMuons",&nMuons, "nMuons/I");
        NoZmassVetoTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        NoZmassVetoTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        NoZmassVetoTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        
        NoZmassVetoTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        NoZmassVetoTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        NoZmassVetoTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        NoZmassVetoTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        NoZmassVetoTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        NoZmassVetoTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        NoZmassVetoTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        NoZmassVetoTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        NoZmassVetoTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        NoZmassVetoTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        NoZmassVetoTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        NoZmassVetoTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        NoZmassVetoTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        NoZmassVetoTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        NoZmassVetoTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");
        
        NoZmassVetoTree->Branch("nJets",&nJets,"nJets/I");
        NoZmassVetoTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        NoZmassVetoTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        NoZmassVetoTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        NoZmassVetoTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        NoZmassVetoTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        NoZmassVetoTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
        NoZmassVetoTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        NoZmassVetoTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        
        NoZmassVetoTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        NoZmassVetoTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        NoZmassVetoTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        
        NoZmassVetoTree->Branch("Ht",&Ht,"Ht/F");
        NoZmassVetoTree->Branch("St",&St,"St/F");
        NoZmassVetoTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        NoZmassVetoTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        NoZmassVetoTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        NoZmassVetoTree->Branch("invMass_2SSL_Zmass",&invMass_2SSL_Zmass,"invMass_2SSL_Zmass/F");
        NoZmassVetoTree->Branch("invMass_2OSL_Zmass",&invMass_2OSL_Zmass,"invMass_2OSL_Zmass/F");
        NoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        NoZmassVetoTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        NoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        NoZmassVetoTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        NoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        NoZmassVetoTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        NoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        NoZmassVetoTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        NoZmassVetoTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        NoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        
        NoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        
        NoZmassVetoTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        NoZmassVetoTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        NoZmassVetoTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        NoZmassVetoTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        NoZmassVetoTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        NoZmassVetoTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        
        NoZmassVetoTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        NoZmassVetoTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        NoZmassVetoTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        NoZmassVetoTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        NoZmassVetoTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        NoZmassVetoTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        NoZmassVetoTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        NoZmassVetoTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        NoZmassVetoTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        
        NoZmassVetoTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        NoZmassVetoTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        NoZmassVetoTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        NoZmassVetoTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        NoZmassVetoTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        NoZmassVetoTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        NoZmassVetoTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        NoZmassVetoTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");

    //    NoZmassVetoTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        
        SSLNoZmassVetoTree->Branch("isData",&isData,"isData/I");
        SSLNoZmassVetoTree->Branch("run_num",&run_num,"run_num/I");
        SSLNoZmassVetoTree->Branch("evt_num",&evt_num,"evt_num/I");
        SSLNoZmassVetoTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        SSLNoZmassVetoTree->Branch("nvtx",&nvtx,"nvtx/I");
        SSLNoZmassVetoTree->Branch("npu",&npu,"npu/I");
        SSLNoZmassVetoTree->Branch("puSF",&puSF,"puSF/D");
        SSLNoZmassVetoTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        SSLNoZmassVetoTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        SSLNoZmassVetoTree->Branch("btagSF",&btagSF,"btagSF/D");
        SSLNoZmassVetoTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        SSLNoZmassVetoTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        SSLNoZmassVetoTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        SSLNoZmassVetoTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        SSLNoZmassVetoTree->Branch("sumW", &sumW, "sumW/I");
        SSLNoZmassVetoTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        SSLNoZmassVetoTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        SSLNoZmassVetoTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        SSLNoZmassVetoTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        SSLNoZmassVetoTree->Branch("ElectronSF",ElectronSF,"ElectronSF[nElectrons]/D");
        
        SSLNoZmassVetoTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        SSLNoZmassVetoTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        SSLNoZmassVetoTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        SSLNoZmassVetoTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        SSLNoZmassVetoTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        SSLNoZmassVetoTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        SSLNoZmassVetoTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        SSLNoZmassVetoTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        SSLNoZmassVetoTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        SSLNoZmassVetoTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        SSLNoZmassVetoTree->Branch("nMuons",&nMuons, "nMuons/I");
        SSLNoZmassVetoTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        SSLNoZmassVetoTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        SSLNoZmassVetoTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        SSLNoZmassVetoTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        SSLNoZmassVetoTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        SSLNoZmassVetoTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        SSLNoZmassVetoTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        SSLNoZmassVetoTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        SSLNoZmassVetoTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        SSLNoZmassVetoTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        SSLNoZmassVetoTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        SSLNoZmassVetoTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        SSLNoZmassVetoTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        SSLNoZmassVetoTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        SSLNoZmassVetoTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        SSLNoZmassVetoTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        SSLNoZmassVetoTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        SSLNoZmassVetoTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");
        
        SSLNoZmassVetoTree->Branch("nJets",&nJets,"nJets/I");
        SSLNoZmassVetoTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        SSLNoZmassVetoTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        SSLNoZmassVetoTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        SSLNoZmassVetoTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        SSLNoZmassVetoTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        SSLNoZmassVetoTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
        SSLNoZmassVetoTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        SSLNoZmassVetoTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        
        SSLNoZmassVetoTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        SSLNoZmassVetoTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        SSLNoZmassVetoTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        
        SSLNoZmassVetoTree->Branch("Ht",&Ht,"Ht/F");
        SSLNoZmassVetoTree->Branch("St",&St,"St/F");
        SSLNoZmassVetoTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        SSLNoZmassVetoTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        SSLNoZmassVetoTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        SSLNoZmassVetoTree->Branch("invMass_2SSL_Zmass",&invMass_2SSL_Zmass,"invMass_2SSL_Zmass/F");
        SSLNoZmassVetoTree->Branch("invMass_2OSL_Zmass",&invMass_2OSL_Zmass,"invMass_2OSL_Zmass/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        SSLNoZmassVetoTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        SSLNoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        SSLNoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
       
        SSLNoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        
        SSLNoZmassVetoTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        SSLNoZmassVetoTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        SSLNoZmassVetoTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        SSLNoZmassVetoTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        SSLNoZmassVetoTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        SSLNoZmassVetoTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        
        SSLNoZmassVetoTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        
        SSLNoZmassVetoTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        SSLNoZmassVetoTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");

        ////////////
        OSLNoZmassVetoTree->Branch("isData",&isData,"isData/I");
        OSLNoZmassVetoTree->Branch("run_num",&run_num,"run_num/I");
        OSLNoZmassVetoTree->Branch("evt_num",&evt_num,"evt_num/I");
        OSLNoZmassVetoTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        OSLNoZmassVetoTree->Branch("nvtx",&nvtx,"nvtx/I");
        OSLNoZmassVetoTree->Branch("npu",&npu,"npu/I");
        OSLNoZmassVetoTree->Branch("puSF",&puSF,"puSF/D");
        OSLNoZmassVetoTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        OSLNoZmassVetoTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        OSLNoZmassVetoTree->Branch("btagSF",&btagSF,"btagSF/D");
        OSLNoZmassVetoTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        OSLNoZmassVetoTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        OSLNoZmassVetoTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        OSLNoZmassVetoTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        OSLNoZmassVetoTree->Branch("sumW", &sumW, "sumW/I");
        OSLNoZmassVetoTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        OSLNoZmassVetoTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        OSLNoZmassVetoTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        OSLNoZmassVetoTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        OSLNoZmassVetoTree->Branch("ElectronSF",ElectronSF,"ElectronSF[nElectrons]/D");
        
        OSLNoZmassVetoTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        OSLNoZmassVetoTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        OSLNoZmassVetoTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        OSLNoZmassVetoTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        OSLNoZmassVetoTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        OSLNoZmassVetoTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        OSLNoZmassVetoTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        OSLNoZmassVetoTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        OSLNoZmassVetoTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        OSLNoZmassVetoTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        OSLNoZmassVetoTree->Branch("nMuons",&nMuons, "nMuons/I");
        OSLNoZmassVetoTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        OSLNoZmassVetoTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        OSLNoZmassVetoTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        OSLNoZmassVetoTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        OSLNoZmassVetoTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        OSLNoZmassVetoTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        OSLNoZmassVetoTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        OSLNoZmassVetoTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        OSLNoZmassVetoTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        OSLNoZmassVetoTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        OSLNoZmassVetoTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        OSLNoZmassVetoTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        OSLNoZmassVetoTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        OSLNoZmassVetoTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        OSLNoZmassVetoTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        OSLNoZmassVetoTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        OSLNoZmassVetoTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        OSLNoZmassVetoTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");
        
        OSLNoZmassVetoTree->Branch("nJets",&nJets,"nJets/I");
        OSLNoZmassVetoTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        OSLNoZmassVetoTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        OSLNoZmassVetoTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        OSLNoZmassVetoTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        OSLNoZmassVetoTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        OSLNoZmassVetoTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
        OSLNoZmassVetoTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        OSLNoZmassVetoTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        
        OSLNoZmassVetoTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        OSLNoZmassVetoTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        OSLNoZmassVetoTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        
        OSLNoZmassVetoTree->Branch("Ht",&Ht,"Ht/F");
        OSLNoZmassVetoTree->Branch("St",&St,"St/F");
        OSLNoZmassVetoTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        OSLNoZmassVetoTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        OSLNoZmassVetoTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        OSLNoZmassVetoTree->Branch("invMass_2SSL_Zmass",&invMass_2SSL_Zmass,"invMass_2SSL_Zmass/F");
        OSLNoZmassVetoTree->Branch("invMass_2OSL_Zmass",&invMass_2OSL_Zmass,"invMass_2OSL_Zmass/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        OSLNoZmassVetoTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        OSLNoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        OSLNoZmassVetoTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        OSLNoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        
        OSLNoZmassVetoTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        OSLNoZmassVetoTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        OSLNoZmassVetoTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        OSLNoZmassVetoTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        OSLNoZmassVetoTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        OSLNoZmassVetoTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        
        OSLNoZmassVetoTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        
        OSLNoZmassVetoTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        OSLNoZmassVetoTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");


        
        /////////////////////////////////////////////////////
        
        myTree->Branch("isData",&isData,"isData/I");
        myTree->Branch("run_num",&run_num,"run_num/I");
        myTree->Branch("evt_num",&evt_num,"evt_num/I");
        myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        myTree->Branch("nvtx",&nvtx,"nvtx/I");
        myTree->Branch("npu",&npu,"npu/I");
        myTree->Branch("puSF",&puSF,"puSF/D");
        myTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        myTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        myTree->Branch("btagSF",&btagSF,"btagSF/D");
        myTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        myTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        myTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        myTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        myTree->Branch("sumW", &sumW, "sumW/I");
        myTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");
        
        myTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        myTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        myTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        myTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        myTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        myTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        myTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        myTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        myTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        
        myTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        myTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        myTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        myTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        myTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        myTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        myTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        myTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");

        SSLeptonTree->Branch("isData",&isData,"isData/I");
        SSLeptonTree->Branch("run_num",&run_num,"run_num/I");
        SSLeptonTree->Branch("evt_num",&evt_num,"evt_num/I");
        SSLeptonTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        SSLeptonTree->Branch("nvtx",&nvtx,"nvtx/I");
        SSLeptonTree->Branch("npu",&npu,"npu/I");
        SSLeptonTree->Branch("puSF",&puSF,"puSF/D");
        SSLeptonTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        SSLeptonTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        SSLeptonTree->Branch("btagSF",&btagSF,"btagSF/D");
        SSLeptonTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        SSLeptonTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        SSLeptonTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        SSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
        SSLeptonTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");
        SSLeptonTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        SSLeptonTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        SSLeptonTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        SSLeptonTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        SSLeptonTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        SSLeptonTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        SSLeptonTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        SSLeptonTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        SSLeptonTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        SSLeptonTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        SSLeptonTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        SSLeptonTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        SSLeptonTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        SSLeptonTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        SSLeptonTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        SSLeptonTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        SSLeptonTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");
        
        OSLeptonTree->Branch("isData",&isData,"isData/I");
        OSLeptonTree->Branch("run_num",&run_num,"run_num/I");
        OSLeptonTree->Branch("evt_num",&evt_num,"evt_num/I");
        OSLeptonTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        OSLeptonTree->Branch("nvtx",&nvtx,"nvtx/I");
        OSLeptonTree->Branch("npu",&npu,"npu/I");
        OSLeptonTree->Branch("puSF",&puSF,"puSF/D");
        OSLeptonTree->Branch("puSF_UP",&puSF_UP,"puSF_UP/D");
        OSLeptonTree->Branch("puSF_Down",&puSF_Down,"puSF_Down/D");
        OSLeptonTree->Branch("btagSF",&btagSF,"btagSF/D");
        OSLeptonTree->Branch("nofPosWeights",&nofPosWeights,"nofPosWeights/I");
        OSLeptonTree->Branch("nofNegWeights",&nofNegWeights,"nofNegWeights/I");
        OSLeptonTree->Branch("nloWeight",&nloWeight,"nloWeight/D");
        OSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
        OSLeptonTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");
        OSLeptonTree->Branch("btagSFshape",&btagSFshape,"btagSFshape/D");
        OSLeptonTree->Branch("btagSFshape_down_cferr1",&btagSFshape_down_cferr1,"btagSFshape_down_cferr1/D");
        OSLeptonTree->Branch("btagSFshape_down_cferr2",&btagSFshape_down_cferr2,"btagSFshape_down_cferr2/D");
        OSLeptonTree->Branch("btagSFshape_down_hf",&btagSFshape_down_hf,"btagSFshape_down_hf/D");
        OSLeptonTree->Branch("btagSFshape_down_hfstats1",&btagSFshape_down_hfstats1,"btagSFshape_down_hfstats1/D");
        OSLeptonTree->Branch("btagSFshape_down_hfstats2",&btagSFshape_down_hfstats2,"btagSFshape_down_hfstats2/D");
        OSLeptonTree->Branch("btagSFshape_down_lf",&btagSFshape_down_lf,"btagSFshape_down_lf/D");
        OSLeptonTree->Branch("btagSFshape_down_lfstats1",&btagSFshape_down_lfstats1,"btagSFshape_down_lfstats1/D");
        OSLeptonTree->Branch("btagSFshape_down_lfstats2",&btagSFshape_down_lfstats2,"btagSFshape_down_lfstats2/D");
        OSLeptonTree->Branch("btagSFshape_up_cferr1",&btagSFshape_up_cferr1,"btagSFshape_up_cferr1/D");
        OSLeptonTree->Branch("btagSFshape_up_cferr2",&btagSFshape_up_cferr2,"btagSFshape_up_cferr2/D");
        OSLeptonTree->Branch("btagSFshape_up_hf",&btagSFshape_up_hf,"btagSFshape_up_hf/D");
        OSLeptonTree->Branch("btagSFshape_up_hfstats1",&btagSFshape_up_hfstats1,"btagSFshape_up_hfstats1/D");
        OSLeptonTree->Branch("btagSFshape_up_hfstats2",&btagSFshape_up_hfstats2,"btagSFshape_up_hfstats2/D");
        OSLeptonTree->Branch("btagSFshape_up_lf",&btagSFshape_up_lf,"btagSFshape_up_lf/D");
        OSLeptonTree->Branch("btagSFshape_up_lfstats1",&btagSFshape_up_lfstats1,"btagSFshape_up_lfstats1/D");
        OSLeptonTree->Branch("btagSFshape_up_lfstats2",&btagSFshape_up_lfstats2,"btagSFshape_up_lfstats2/D");

        
        // Set branches for different Trees
        ////--> Electrons <----////
        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
        myTree->Branch("ElectronSF_up",&ElectronSF_up,"ElectronSF_up[nElectrons]/D");
        myTree->Branch("ElectronSF_down",&ElectronSF_down,"ElectronSF_down[nElectrons]/D");
        myTree->Branch("ElectronSFID",&ElectronSFID,"ElectronSFID[nElectrons]/D");
        myTree->Branch("ElectronSFID_up",&ElectronSFID_up,"ElectronSFID_up[nElectrons]/D");
        myTree->Branch("ElectronSFID_down",&ElectronSFID_down,"ElectronSFID_down[nElectrons]/D");
        myTree->Branch("ElectronSFReco",&ElectronSFReco,"ElectronSFReco[nElectrons]/D");
        myTree->Branch("ElectronSFReco_up",&ElectronSFReco_up,"ElectronSFReco_up[nElectrons]/D");
        myTree->Branch("ElectronSFReco_down",&ElectronSFReco_down,"ElectronSFReco_down[nElectrons]/D");
        
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
        myTree->Branch("isSelectiveElec",isSelectiveElec,"isSelectiveElec[nElectrons]/O)");
        
        
        myTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        myTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        myTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        myTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        myTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        myTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        myTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        myTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        myTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        myTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        myTree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
        myTree->Branch("MVA_nCSVLbtagJets", &MVA_nCSVLbtagJets, "MVA_nCSVLbtagJets/F");
       
        
        SSLeptonTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        SSLeptonTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSF_up",&ElectronSF_up,"ElectronSF_up[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSF_down",&ElectronSF_down,"ElectronSF_down[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFID",&ElectronSFID,"ElectronSFID[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFID_up",&ElectronSFID_up,"ElectronSFID_up[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFID_down",&ElectronSFID_down,"ElectronSFID_down[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFReco",&ElectronSFReco,"ElectronSFReco[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFReco_up",&ElectronSFReco_up,"ElectronSFReco_up[nElectrons]/D");
        SSLeptonTree->Branch("ElectronSFReco_down",&ElectronSFReco_down,"ElectronSFReco_down[nElectrons]/D");
        SSLeptonTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        SSLeptonTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        SSLeptonTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        SSLeptonTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        SSLeptonTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        SSLeptonTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        SSLeptonTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        SSLeptonTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        SSLeptonTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        SSLeptonTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        SSLeptonTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        SSLeptonTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        SSLeptonTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        SSLeptonTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        SSLeptonTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        SSLeptonTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        SSLeptonTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        SSLeptonTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        SSLeptonTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        SSLeptonTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        SSLeptonTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        SSLeptonTree->Branch("ElectronSF",ElectronSF,"ElectronSF[nElectrons]/D");
        SSLeptonTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        SSLeptonTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        SSLeptonTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        SSLeptonTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        SSLeptonTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        SSLeptonTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        SSLeptonTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        SSLeptonTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        SSLeptonTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        SSLeptonTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        SSLeptonTree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
        SSLeptonTree->Branch("MVA_nCSVLbtagJets", &MVA_nCSVLbtagJets, "MVA_nCSVLbtagJets/F");

        
        OSLeptonTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        OSLeptonTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSF_up",&ElectronSF_up,"ElectronSF_up[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSF_down",&ElectronSF_down,"ElectronSF_down[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFID",&ElectronSFID,"ElectronSFID[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFID_up",&ElectronSFID_up,"ElectronSFID_up[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFID_down",&ElectronSFID_down,"ElectronSFID_down[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFReco",&ElectronSFReco,"ElectronSFReco[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFReco_up",&ElectronSFReco_up,"ElectronSFReco_up[nElectrons]/D");
        OSLeptonTree->Branch("ElectronSFReco_down",&ElectronSFReco_down,"ElectronSFReco_down[nElectrons]/D");
        OSLeptonTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        OSLeptonTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        OSLeptonTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        OSLeptonTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        OSLeptonTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        OSLeptonTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        OSLeptonTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        OSLeptonTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        OSLeptonTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        OSLeptonTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        OSLeptonTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        OSLeptonTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        OSLeptonTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        OSLeptonTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        OSLeptonTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        OSLeptonTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        OSLeptonTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        OSLeptonTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        OSLeptonTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        OSLeptonTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        OSLeptonTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        OSLeptonTree->Branch("ElectronSF",ElectronSF,"ElectronSF[nElectrons]/D");
        OSLeptonTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        OSLeptonTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        OSLeptonTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        OSLeptonTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        OSLeptonTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        OSLeptonTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        OSLeptonTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        OSLeptonTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        OSLeptonTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        OSLeptonTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        OSLeptonTree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
        OSLeptonTree->Branch("MVA_nCSVLbtagJets", &MVA_nCSVLbtagJets, "MVA_nCSVLbtagJets/F");

        
        
        ////--> muons <----////
        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        myTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        myTree->Branch("MuonTrackSF_up",&MuonTrackSF_up, "MuonTrackSF_up[nMuons]/D");
        myTree->Branch("MuonTrackSF_down",&MuonTrackSF_down, "MuonTrackSF_down[nMuons]/D");
        myTree->Branch("MuonIDSF_up",&MuonIDSF_up,"MuonIDSF_up[nMuons]/D");
        myTree->Branch("MuonIsoSF_up",&MuonIsoSF_up, "MuonIsoSF_up[nMuons]/D");
        myTree->Branch("MuonIDSF_down",&MuonIDSF_down,"MuonIDSF_down[nMuons]/D");
        myTree->Branch("MuonIsoSF_down",&MuonIsoSF_down, "MuonIsoSF_down[nMuons]/D");
        

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
        
        myTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        myTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        myTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        myTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        myTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        myTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        myTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        myTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        myTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        myTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");
       
        
        SSLeptonTree->Branch("nMuons",&nMuons, "nMuons/I");
        SSLeptonTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        SSLeptonTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        SSLeptonTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        SSLeptonTree->Branch("MuonTrackSF_up",&MuonTrackSF_up, "MuonTrackSF_up[nMuons]/D");
        SSLeptonTree->Branch("MuonTrackSF_down",&MuonTrackSF_down, "MuonTrackSF_down[nMuons]/D");
        SSLeptonTree->Branch("MuonIDSF_up",&MuonIDSF_up,"MuonIDSF_up[nMuons]/D");
        SSLeptonTree->Branch("MuonIsoSF_up",&MuonIsoSF_up, "MuonIsoSF_up[nMuons]/D");
        SSLeptonTree->Branch("MuonIDSF_down",&MuonIDSF_down,"MuonIDSF_down[nMuons]/D");
        SSLeptonTree->Branch("MuonIsoSF_down",&MuonIsoSF_down, "MuonIsoSF_down[nMuons]/D");        //SSLeptonTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // SSLeptonTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        SSLeptonTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        SSLeptonTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        SSLeptonTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        SSLeptonTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        SSLeptonTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        SSLeptonTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        SSLeptonTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        SSLeptonTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        SSLeptonTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        SSLeptonTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        SSLeptonTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/F");
        //        SSLeptonTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        //        SSLeptonTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        SSLeptonTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        SSLeptonTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        SSLeptonTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        SSLeptonTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        SSLeptonTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        SSLeptonTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        SSLeptonTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        SSLeptonTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        SSLeptonTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        SSLeptonTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        SSLeptonTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");

        OSLeptonTree->Branch("nMuons",&nMuons, "nMuons/I");
        OSLeptonTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        OSLeptonTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        OSLeptonTree->Branch("MuonTrackSF",&MuonTrackSF, "MuonTrackSF[nMuons]/D");
        OSLeptonTree->Branch("MuonTrackSF_up",&MuonTrackSF_up, "MuonTrackSF_up[nMuons]/D");
        OSLeptonTree->Branch("MuonTrackSF_down",&MuonTrackSF_down, "MuonTrackSF_down[nMuons]/D");
        OSLeptonTree->Branch("MuonIDSF_up",&MuonIDSF_up,"MuonIDSF_up[nMuons]/D");
        OSLeptonTree->Branch("MuonIsoSF_up",&MuonIsoSF_up, "MuonIsoSF_up[nMuons]/D");
        OSLeptonTree->Branch("MuonIDSF_down",&MuonIDSF_down,"MuonIDSF_down[nMuons]/D");
        OSLeptonTree->Branch("MuonIsoSF_down",&MuonIsoSF_down, "MuonIsoSF_down[nMuons]/D");
        //OSLeptonTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // OSLeptonTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        OSLeptonTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        OSLeptonTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        OSLeptonTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        OSLeptonTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        OSLeptonTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        OSLeptonTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        OSLeptonTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        OSLeptonTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        OSLeptonTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        OSLeptonTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        OSLeptonTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/D");
        OSLeptonTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        OSLeptonTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        OSLeptonTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        
        OSLeptonTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        OSLeptonTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        OSLeptonTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        OSLeptonTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        OSLeptonTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        OSLeptonTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        OSLeptonTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        OSLeptonTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        OSLeptonTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        OSLeptonTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");
        

        ////--> jets <----////
        
        myTree->Branch("nJets",&nJets,"nJets/I");
        myTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        myTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        myTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        myTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        myTree->Branch("bdisc_jet",bdisc_jet,"bdisc_jet[nJets]/D");
        myTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        myTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        myTree->Branch("bdiscCSVv2_1st_jet", &bdiscCSVv2_1st_jet, "bdiscCSVv2_1st_jet/F");
        myTree->Branch("bdiscCSVv2_2nd_jet", &bdiscCSVv2_2nd_jet, "bdiscCSVv2_2nd_jet/F");
        myTree->Branch("bdiscCSVv2_3rd_jet", &bdiscCSVv2_3rd_jet, "bdiscCSVv2_3rd_jet/F");
        myTree->Branch("bdiscCSVv2_4th_jet", &bdiscCSVv2_4th_jet, "bdiscCSVv2_4th_jet/F");
        
        
        myTree->Branch("pt_1st_jet", &pt_1st_jet, "pt_1st_jet/F");
        myTree->Branch("pt_2nd_jet", &pt_2nd_jet, "pt_2nd_jet/F");
        myTree->Branch("pt_3rd_jet", &pt_3rd_jet, "pt_3rd_jet/F");
        myTree->Branch("pt_4th_jet", &pt_4th_jet, "pt_4th_jet/F");
     
        SSLeptonTree->Branch("nJets",&nJets,"nJets/I");
        SSLeptonTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        SSLeptonTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        SSLeptonTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        SSLeptonTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        SSLeptonTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        SSLeptonTree->Branch("pt_1st_jet", &pt_1st_jet, "pt_1st_jet/F");
        SSLeptonTree->Branch("pt_2nd_jet", &pt_2nd_jet, "pt_2nd_jet/F");
        SSLeptonTree->Branch("pt_3rd_jet", &pt_3rd_jet, "pt_3rd_jet/F");
        SSLeptonTree->Branch("pt_4th_jet", &pt_4th_jet, "pt_4th_jet/F");
        SSLeptonTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        SSLeptonTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        SSLeptonTree->Branch("bdiscCSVv2_1st_jet", &bdiscCSVv2_1st_jet, "bdiscCSVv2_1st_jet/F");
        SSLeptonTree->Branch("bdiscCSVv2_2nd_jet", &bdiscCSVv2_2nd_jet, "bdiscCSVv2_2nd_jet/F");
        SSLeptonTree->Branch("bdiscCSVv2_3rd_jet", &bdiscCSVv2_3rd_jet, "bdiscCSVv2_3rd_jet/F");
        SSLeptonTree->Branch("bdiscCSVv2_4th_jet", &bdiscCSVv2_4th_jet, "bdiscCSVv2_4th_jet/F");
        
        OSLeptonTree->Branch("nJets",&nJets,"nJets/I");
        OSLeptonTree->Branch("pt_jet",pt_jet,"pt_jet[nJets]/D");
        OSLeptonTree->Branch("phi_jet",phi_jet,"phi_jet[nJets]/D");
        OSLeptonTree->Branch("eta_jet",eta_jet,"eta_jet[nJets]/D");
        OSLeptonTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
        OSLeptonTree->Branch("charge_jet",charge_jet,"charge_jet[nJets]/I");
        OSLeptonTree->Branch("pt_1st_jet", &pt_1st_jet, "pt_1st_jet/F");
        OSLeptonTree->Branch("pt_2nd_jet", &pt_2nd_jet, "pt_2nd_jet/F");
        OSLeptonTree->Branch("pt_3rd_jet", &pt_3rd_jet, "pt_3rd_jet/F");
        OSLeptonTree->Branch("pt_4th_jet", &pt_4th_jet, "pt_4th_jet/F");
        OSLeptonTree->Branch("bdiscCSVv2_1st_bjet", &bdiscCSVv2_1st_bjet, "bdiscCSVv2_1st_bjet/F");
        OSLeptonTree->Branch("bdiscCSVv2_2nd_bjet", &bdiscCSVv2_2nd_bjet, "bdiscCSVv2_2nd_bjet/F");
        OSLeptonTree->Branch("bdiscCSVv2_1st_jet", &bdiscCSVv2_1st_jet, "bdiscCSVv2_1st_jet/F");
        OSLeptonTree->Branch("bdiscCSVv2_2nd_jet", &bdiscCSVv2_2nd_jet, "bdiscCSVv2_2nd_jet/F");
        OSLeptonTree->Branch("bdiscCSVv2_3rd_jet", &bdiscCSVv2_3rd_jet, "bdiscCSVv2_3rd_jet/F");
        OSLeptonTree->Branch("bdiscCSVv2_4th_jet", &bdiscCSVv2_4th_jet, "bdiscCSVv2_4th_jet/F");
        
        
        ////--> met <----////
        myTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        myTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        myTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        myTree->Branch("corrJES_met_Pt", &corrJES_met_Pt, "corrJES_met_Pt/F");
        myTree->Branch("UncorrJES_met_Pt", &UncorrJES_met_Pt, "UncorrJES_met_Pt/F");
        myTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        myTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        myTree->Branch("PassedMETFilter", &PassedMETFilter,"PassedMETFilter/I");
        
        myTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        myTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        myTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        myTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        myTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        myTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        myTree->Branch("DeltaPhi_met_lep1", &DeltaPhi_met_lep1, "DeltaPhi_met_lep1/F");
        myTree->Branch("DeltaPhi_met_lep2", &DeltaPhi_met_lep2, "DeltaPhi_met_lep2/F");

        SSLeptonTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        SSLeptonTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        SSLeptonTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        SSLeptonTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        SSLeptonTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        SSLeptonTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        SSLeptonTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        SSLeptonTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        SSLeptonTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        SSLeptonTree->Branch("DeltaPhi_met_lep1", &DeltaPhi_met_lep1, "DeltaPhi_met_lep1/F");
        SSLeptonTree->Branch("DeltaPhi_met_lep2", &DeltaPhi_met_lep2, "DeltaPhi_met_lep2/F");

        
        OSLeptonTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        OSLeptonTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        OSLeptonTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        OSLeptonTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        OSLeptonTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        OSLeptonTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        OSLeptonTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        OSLeptonTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        OSLeptonTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");
        OSLeptonTree->Branch("DeltaPhi_met_lep1", &DeltaPhi_met_lep1, "DeltaPhi_met_lep1/F");
        OSLeptonTree->Branch("DeltaPhi_met_lep2", &DeltaPhi_met_lep2, "DeltaPhi_met_lep2/F");
        
        InitialTree->Branch("UncorrJES_met_Pt", &UncorrJES_met_Pt, "UncorrJES_met_Pt/F");
        InitialTree->Branch("corrJES_met_Pt", &corrJES_met_Pt, "corrJES_met_Pt/F");
        InitialTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        InitialTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        SSLeptonTree->Branch("UncorrJES_met_Pt", &UncorrJES_met_Pt, "UncorrJES_met_Pt/F");
        SSLeptonTree->Branch("corrJES_met_Pt", &corrJES_met_Pt, "corrJES_met_Pt/F");
        SSLeptonTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        SSLeptonTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        OSLeptonTree->Branch("UncorrJES_met_Pt", &UncorrJES_met_Pt, "UncorrJES_met_Pt/F");
        OSLeptonTree->Branch("corrJES_met_Pt", &corrJES_met_Pt, "corrJES_met_Pt/F");
        OSLeptonTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        OSLeptonTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        ////--> bJets <----////
        myTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        myTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        myTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        myTree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        myTree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        myTree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        
        
        SSLeptonTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        SSLeptonTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        SSLeptonTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        SSLeptonTree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        SSLeptonTree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        SSLeptonTree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        
        OSLeptonTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        OSLeptonTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        OSLeptonTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        OSLeptonTree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        OSLeptonTree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        OSLeptonTree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        
        ////--> MC Particles <---/////
//        
        
        ///// --> MVA variables <--- ////
        myTree->Branch("Ht",&Ht,"Ht/F");
        myTree->Branch("St",&St,"St/F");
        myTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        myTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        myTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        myTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        myTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        myTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        myTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        myTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        myTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        myTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        myTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        myTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        myTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
   //     myTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        /*
        myTree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
        myTree->Branch("Match_MC_St",&Match_MC_St,"Match_MC_St/F");
        myTree->Branch("Match_MC_DeltaR_2L",&Match_MC_DeltaR_2L,"Match_MC_DeltaR_2L/F");
        myTree->Branch("Match_MC_DeltaPhi_2L",&Match_MC_DeltaPhi_2L,"Match_MC_DeltaPhi_2L/F");
        myTree->Branch("Match_MC_invMass_2L",&Match_MC_invMass_2L,"Match_MC_invMass_2L/F");
        myTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMu",&Match_SM_MC_DeltaR_Mu0b0_DiMu,"Match_SM_MC_DeltaR_Mu0b0_DiMu/F");
        myTree->Branch("Match_MC_DeltaR_Mu1b0_DiMu",&Match_MC_DeltaR_Mu1b0_DiMu,"Match_MC_DeltaR_Mu1b0_DiMu/F");
        myTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElec",&Match_SM_MC_DeltaR_Elec0b0_DiElec,"Match_SM_MC_DeltaR_Elec0b0_DiElec/F");
        myTree->Branch("Match_MC_DeltaR_Elec1b0_DiElec",&Match_MC_DeltaR_Elec1b0_DiElec,"Match_MC_DeltaR_Elec1b0_DiElec/F");
        myTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMuEl",&Match_SM_MC_DeltaR_Mu0b0_DiMuEl,"Match_SM_MC_DeltaR_Mu0b0_DiMuEl/F");
        myTree->Branch("Match_MC_DeltaR_Mu1b0_DiElMu",&Match_MC_DeltaR_Mu1b0_DiElMu,"Match_MC_DeltaR_Mu1b0_DiElMu/F");
        myTree->Branch("Match_MC_DeltaR_Elec1b0_DiMuEl",&Match_MC_DeltaR_Elec1b0_DiMuEl,"Match_MC_DeltaR_Elec1b0_DiMuEl/F");
        myTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElMu",&Match_SM_MC_DeltaR_Elec0b0_DiElMu,"Match_SM_MC_DeltaR_Elec0b0_DiElMu/F");
        myTree->Branch("Match_FCNC_MC_Mass_WJetPair",&Match_FCNC_MC_Mass_WJetPair,"Match_FCNC_MC_Mass_WJetPair/F");
        myTree->Branch("Match_FCNC_MC_Mass_JetPair",&Match_FCNC_MC_Mass_JetPair,"Match_FCNC_MC_Mass_JetPair/F");
        */
        
        SSLeptonTree->Branch("Ht",&Ht,"Ht/F");
        SSLeptonTree->Branch("St",&St,"St/F");
        SSLeptonTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        SSLeptonTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        SSLeptonTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        SSLeptonTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        SSLeptonTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        SSLeptonTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        SSLeptonTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        SSLeptonTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        SSLeptonTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        SSLeptonTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        SSLeptonTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        SSLeptonTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        SSLeptonTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        /* SSLeptonTree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
        SSLeptonTree->Branch("Match_MC_St",&Match_MC_St,"Match_MC_St/F");
        SSLeptonTree->Branch("Match_MC_DeltaR_2L",&Match_MC_DeltaR_2L,"Match_MC_DeltaR_2L/F");
        SSLeptonTree->Branch("Match_MC_DeltaPhi_2L",&Match_MC_DeltaPhi_2L,"Match_MC_DeltaPhi_2L/F");
        SSLeptonTree->Branch("Match_MC_invMass_2L",&Match_MC_invMass_2L,"Match_MC_invMass_2L/F");
        SSLeptonTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMu",&Match_SM_MC_DeltaR_Mu0b0_DiMu,"Match_SM_MC_DeltaR_Mu0b0_DiMu/F");
        SSLeptonTree->Branch("Match_MC_DeltaR_Mu1b0_DiMu",&Match_MC_DeltaR_Mu1b0_DiMu,"Match_MC_DeltaR_Mu1b0_DiMu/F");
        SSLeptonTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElec",&Match_SM_MC_DeltaR_Elec0b0_DiElec,"Match_SM_MC_DeltaR_Elec0b0_DiElec/F");
        SSLeptonTree->Branch("Match_MC_DeltaR_Elec1b0_DiElec",&Match_MC_DeltaR_Elec1b0_DiElec,"Match_MC_DeltaR_Elec1b0_DiElec/F");
        SSLeptonTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMuEl",&Match_SM_MC_DeltaR_Mu0b0_DiMuEl,"Match_SM_MC_DeltaR_Mu0b0_DiMuEl/F");
        SSLeptonTree->Branch("Match_MC_DeltaR_Mu1b0_DiElMu",&Match_MC_DeltaR_Mu1b0_DiElMu,"Match_MC_DeltaR_Mu1b0_DiElMu/F");
        SSLeptonTree->Branch("Match_MC_DeltaR_Elec1b0_DiMuEl",&Match_MC_DeltaR_Elec1b0_DiMuEl,"Match_MC_DeltaR_Elec1b0_DiMuEl/F");
        SSLeptonTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElMu",&Match_SM_MC_DeltaR_Elec0b0_DiElMu,"Match_SM_MC_DeltaR_Elec0b0_DiElMu/F");
        SSLeptonTree->Branch("Match_FCNC_MC_Mass_WJetPair",&Match_FCNC_MC_Mass_WJetPair,"Match_FCNC_MC_Mass_WJetPair/F");
        SSLeptonTree->Branch("Match_FCNC_MC_Mass_JetPair",&Match_FCNC_MC_Mass_JetPair,"Match_FCNC_MC_Mass_JetPair/F");
    //    SSLeptonTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        */
        OSLeptonTree->Branch("Ht",&Ht,"Ht/F");
        OSLeptonTree->Branch("St",&St,"St/F");
        OSLeptonTree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        OSLeptonTree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        OSLeptonTree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        OSLeptonTree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        OSLeptonTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        OSLeptonTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        OSLeptonTree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        OSLeptonTree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        OSLeptonTree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        OSLeptonTree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        OSLeptonTree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        OSLeptonTree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        OSLeptonTree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        /*OSLeptonTree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
        OSLeptonTree->Branch("Match_MC_St",&Match_MC_St,"Match_MC_St/F");
        OSLeptonTree->Branch("Match_MC_DeltaR_2L",&Match_MC_DeltaR_2L,"Match_MC_DeltaR_2L/F");
        OSLeptonTree->Branch("Match_MC_DeltaPhi_2L",&Match_MC_DeltaPhi_2L,"Match_MC_DeltaPhi_2L/F");
        OSLeptonTree->Branch("Match_MC_invMass_2L",&Match_MC_invMass_2L,"Match_MC_invMass_2L/F");
        OSLeptonTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMu",&Match_SM_MC_DeltaR_Mu0b0_DiMu,"Match_SM_MC_DeltaR_Mu0b0_DiMu/F");
        OSLeptonTree->Branch("Match_MC_DeltaR_Mu1b0_DiMu",&Match_MC_DeltaR_Mu1b0_DiMu,"Match_MC_DeltaR_Mu1b0_DiMu/F");
        OSLeptonTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElec",&Match_SM_MC_DeltaR_Elec0b0_DiElec,"Match_SM_MC_DeltaR_Elec0b0_DiElec/F");
        OSLeptonTree->Branch("Match_MC_DeltaR_Elec1b0_DiElec",&Match_MC_DeltaR_Elec1b0_DiElec,"Match_MC_DeltaR_Elec1b0_DiElec/F");
        OSLeptonTree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMuEl",&Match_SM_MC_DeltaR_Mu0b0_DiMuEl,"Match_SM_MC_DeltaR_Mu0b0_DiMuEl/F");
        OSLeptonTree->Branch("Match_MC_DeltaR_Mu1b0_DiElMu",&Match_MC_DeltaR_Mu1b0_DiElMu,"Match_MC_DeltaR_Mu1b0_DiElMu/F");
        OSLeptonTree->Branch("Match_MC_DeltaR_Elec1b0_DiMuEl",&Match_MC_DeltaR_Elec1b0_DiMuEl,"Match_MC_DeltaR_Elec1b0_DiMuEl/F");
        OSLeptonTree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElMu",&Match_SM_MC_DeltaR_Elec0b0_DiElMu,"Match_SM_MC_DeltaR_Elec0b0_DiElMu/F");
        OSLeptonTree->Branch("Match_FCNC_MC_Mass_WJetPair",&Match_FCNC_MC_Mass_WJetPair,"Match_FCNC_MC_Mass_WJetPair/F");
        OSLeptonTree->Branch("Match_FCNC_MC_Mass_JetPair",&Match_FCNC_MC_Mass_JetPair,"Match_FCNC_MC_Mass_JetPair/F");
     //   OSLeptonTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        */
        ///// MVA Tree /////
        
        
        ////////MC Particle information ////
        myTree->Branch("Jet_matchMC_pdgId",&Jet_matchMC_pdgId,"Jet_matchMC_pdgId[nJets]/D");
        myTree->Branch("Jet_matchMC_motherpdgId",&Jet_matchMC_motherpdgId,"Jet_matchMC_motherpdgId[nJets]/D");
        myTree->Branch("Jet_matchMC_grannypdgId",&Jet_matchMC_grannypdgId,"Jet_matchMC_grannypdgId[nJets]/D");



        ////////////////////////////////////
        ///  Loop on events
        ////////////////////////////////////
        //some bookkeeping variables
      //  Double_t scaleFactor = 1.;
        bool InBarrel = false;
        bool InEndcap = false;
        vector<double> NSSL_Barrel_Vec;
        vector<double> NSSL_EndCap_Vec;
        vector<double> NOSL_Barrel_Vec;
        vector<double> NOSL_EndCap_Vec;
        vector<double> Charge_misId_BVec;
        vector<double> Charge_misId_EVec;
        //                            double Charge_misId_BVec[5];
        //                            double Charge_misId_EVec[5];
        double Nb_2SSl_B1, Nb_2SSl_B2, Nb_2SSl_B3,Nb_2SSl_E1, Nb_2SSl_E2, Nb_2SSl_E3,Nb_2OSl_B1, Nb_2OSl_B2, Nb_2OSl_B3, Nb_2OSl_E1, Nb_2OSl_E2, Nb_2OSl_E3;
        Nb_2SSl_B1= Nb_2SSl_B2= Nb_2SSl_B3=Nb_2SSl_E1=Nb_2SSl_E2= Nb_2SSl_E3=Nb_2OSl_B1= Nb_2OSl_B2= Nb_2OSl_B3= Nb_2OSl_E1= Nb_2OSl_E2= Nb_2OSl_E3 =0;
        double Charge_misId_Ratio , Charge_misId_B1, Charge_misId_B2, Charge_misId_B3 , Charge_misId_E1, Charge_misId_E2, Charge_misId_E3 ;
        Charge_misId_Ratio = Charge_misId_B1 = Charge_misId_B2 = Charge_misId_B3 = Charge_misId_E1= Charge_misId_E2= Charge_misId_E3 =  0;

        nbEvents = 0;
        int nbEvents_eEqLumi = 0;
        int previousRun = -1;
        int currentRun;
        bool Btagged = false;
        nofPosWeights = 0;
        nofNegWeights = 0;
        int nbEvents_0 = 0;
        int nbEvents_1PU = 0;
        int nbEvents_4BTag = 0;
        int nbEvents_3Trig=0;
        int nbEvents_2GPV = 0;
        int nbEvents_5LepSF =0;
        int nbEvents_6 = 0;
        int nbEvents_7 = 0;
        int nbEvents_8 = 0;
        int nbEvents_9 = 0;
        int nbEvents_10 = 0;
        int nbEvents_11 = 0;
        int nb_2l_Jets = 0;
        int nb_2l_CSVLbJets = 0;
        int nb_2l_CSVMbJets = 0;
        int nb_2l_CSVTbJets = 0;
        
        //// Met Filters
        
        bool   HBHEnoise = false;
        bool   HBHEIso = false;
        bool   CSCTight = false;
        bool badchan = false;
        bool badmu = false;
        bool EcalDead = false;
        bool passedMET = false;
        int nbofElecbeforeChargeCut = 0, nbofElecAfterChargeCut = 0 , nbofEventsbeforeChargeCut = 0 , nbofEventsAfterChargeCut = 0;
        
        int itrigger = -1;
        
        
        if (Apply_HLT_Triggers) {
            //DilepTrigger->bookTriggers(isData);
            if(Elec_Elec){DiElecTrigger->bookTriggers(isData);}
            if (Mu_Mu){DiMuTrigger->bookTriggers(isData);}
            if (Elec_Mu){DiElMuTrigger->bookTriggers(isData);}
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
        vector<TRootPFJet*> selectedOrgiJets;
        vector<TRootPFJet*> selectedJets;
        vector < TRootMuon* > selectedMuons;
        vector < TRootMuon* > selectedOrgiMuons;
        vector < TRootElectron* > selectedElectrons;
        vector < TRootMuon* > selectedOrgiLooseMuons;
        vector < TRootMuon* > selectedLooseMuons;
        vector < TRootElectron* > selectedLooseElectrons;
        vector<TRootJet*> selectedBCSVLJets;
        vector<TRootJet*> selectedBCSVMJets;
        vector<TRootJet*> selectedBCSVTJets;
        vector<TRootJet*> selectedNonBCSVLJets;
        vector<TRootJet*> selectedNonBCSVMJets;
        vector<TRootJet*> selectedNonBCSVTJets;
        vector<TRootJet*> JetsExcludingHighestCSVLb;
        vector<TRootPFJet*> selectedCSVLLJets;
        vector<TRootPFJet*> selectedCSVMLJets;
        vector<TRootPFJet*> selectedCSVTLJets;
        vector<TLorentzVector> selectedMuonsTLV, selectedElectronsTLV, metsTLV, selectedJetsTLV, selectedBJetsTLV, selectedLightJetsTLV, selectedLeptonsTLV, JetsExcludingHighestCSVLbTLV;
        vector<TLorentzVector> mcParticlesTLV, mcMuonsTLV, mcPartonsTLV;
        vector<TRootMCParticle*> mcParticlesMatching_;
        vector<int> mcMuonIndex, mcPartonIndex;
        JetPartonMatching muonMatching, jetMatching;
        float rho;
       // vector<int> electroncharge;
        
        
        /////// ************* ///////////
        /// start looping on events ///
        ////////************** /////////
        
        int event_start = startEvent;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();
        cout <<"Number of events = "<<  ending  <<endl;
        SampleLumiWeight = lumiWeight;
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
            nbEvents_eEqLumi++;
            histo1D["h_Nb_Events_Lumi"]->Fill(1);
            vertex.clear();
            init_muons.clear();
            init_jets.clear();
            init_jets_corrected.clear();
            genjets.clear();
            mets.clear();
           // nCuts = 0;
            passedMET = false;
            HBHEnoise = false;
            HBHEIso = false;
            CSCTight = false;
            badchan = false;
            badmu = false;
            EcalDead = false;
          //  cout << "the nCut Value just after start looping over events before any cuts is =  " << nCuts << endl;
            
            
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
            
            ///// define met Filters ////
            HBHEnoise = event->getHBHENoiseFilter();
            HBHEIso = event->getHBHENoiseIsoFilter();
            CSCTight = event->getglobalTightHalo2016Filter();
            EcalDead = event->getEcalDeadCellTriggerPrimitiveFilter();
            badchan   = event-> getBadChCandFilter();
            badmu = event-> getBadPFMuonFilter();

            /////////////////////////////////////
            //  fix negative weights for amc@nlo///
            /////////////////////////////////////
         //   cout << "the number of events before calculating negative weights is =  " << datasets[d]->NofEvtsToRunOver() << endl;
            if(debug) cout << "amc fixing" << endl;
            double hasNegWeight = false;
            double mc_baseweight = 1;
            
            //if(!isData && !isSignal && (event->getWeight(1001) != -9999.))
            if(!isData && (event->getWeight(1001) != -9999.))
            {
                mc_baseweight =  event->getWeight(1001)/abs(event->originalXWGTUP());
                //mc_scaleupweight = event->getWeight(1005)/abs(event->originalXWGTUP());
                //mc_scaledownweight = event->getWeight(1009)/abs(event->originalXWGTUP());
                if(mc_baseweight >= 0)
                {
                    nofPosWeights++;
                  //  histo1D["weightIndex"]->Fill(1.,1.);
                    
                }
                else
                {
                    if(nlo) hasNegWeight = true;
                    nofNegWeights++;
                  //  histo1D["weightIndex"]->Fill(-1.,1.);
                }
            }
            //if( !isData && !isSignal && (event->getWeight(1) != -9999. ))
            if( !isData && (event->getWeight(1) != -9999. ))
            {
                mc_baseweight =  event->getWeight(1)/abs(event->originalXWGTUP());
                //mc_scaleupweight = event->getWeight(5)/abs(event->originalXWGTUP());
                //mc_scaledownweight = event->getWeight(9)/abs(event->originalXWGTUP());
                if(mc_baseweight >= 0)
                {
                    nofPosWeights++;
                   // histo1D["weightIndex"]->Fill(2.,1.);
                    
                }
                else
                {
                    if(nlo) hasNegWeight = true;
                    nofNegWeights++;
                  //  histo1D["weightIndex"]->Fill(-2.,1.);
                }
                
                
            }
           // if(!isData && !isSignal)
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
                
                
            }
           if(debug) cout << "the number of events after calculating negative weights is =  " << datasets[d]->NofEvtsToRunOver() << endl;
         
           //// Trigger.checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
            if (Apply_HLT_Triggers)
            {
                
                if (Elec_Elec)
                {
                    DiElecTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
                    itrigger = DiElecTrigger->checkIfFired();
                }
                
                if (Mu_Mu) {
                    DiMuTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
                    itrigger = DiMuTrigger->checkIfFired();

                }
                if (Elec_Mu) {
                    DiElMuTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
                    itrigger = DiElMuTrigger->checkIfFired();

                }
            }
            
            ////////////////////////////
            ///// JES - JER smearing     ////
            //////////////////////////
            if(applyJER && !isData)
            {
                if(doJERShift == 1)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus", false);
                else if(doJERShift == 2)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus", false);
                else
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
            }else JERon = -1;
            
            
            if(applyJES && !isData)
            {
                JESon = doJESShift;
                if (doJESShift == 1)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
                else if (doJESShift == 2)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
                
                jetTools->correctJets(init_jets,event->fixedGridRhoFastjetAll() ,false);
            } else JESon = -1;
            
            UncorrJES_met_Pt = mets[0]->Pt();
            
            if(applyJES) // jer doesn't need to be applied  --> smeared type-1 corrected MET,  NOW only yes --> Type 1 corrected MET
            {
                
                jetTools->correctMETTypeOne(init_jets, mets[0], isData);
                //  METon = 1;
                //  if JES applied: replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which JEC is applied
                //  if JER applied:  replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which smearing is applied.
                // type 1 correction / sleard pmet correction
                
            }
            
            corrJES_met_Pt = mets[0]->Pt();
            //corrJES_met_Pt = sqrt(UncorrJES_met_Pt*UncorrJES_met_Pt + corrected_met_py*corrected_met_py);

            
            /////////////////////////
            ///  EVENT SELECTION  ///
            /////////////////////////
            selectedElectrons.clear();
            selectedMuons.clear();
            selectedOrgiMuons.clear();
            selectedOrgiJets.clear();
            mcParticles.clear();
            selectedLooseMuons.clear();
            selectedOrgiLooseMuons.clear();
            selectedLooseElectrons.clear();
            bool isBadPfMuon = false;
            bool isDuplicatedMuon = false;
            
            
            ////// **** Removing Bad muons from initmuon collection **** ////
            if (Mu_Mu || Elec_Mu)
            {
                for(int iMu = 0 ; iMu < init_muons.size(); iMu++)
                {
                    if(init_muons[iMu]->Pt() > 10.0)
                    {
                        isBadPfMuon = init_muons[iMu]->isBad80X();
                        isDuplicatedMuon = init_muons[iMu]->isClone80X();
                        
//                        if (isBadPfMuon || isDuplicatedMuon)
//                        {
//                            erase(init_muons[iMu]);
//                        }
                    }
                }

            }
            
            ///// *** use the "selective" method for Electron charge ** //
            bool is3ChargeAgr = false;
            if (Elec_Elec)
            {
                int nbElecWithNo3ChargeAgr = 0;
                nbofEventsbeforeChargeCut++;
                histo1D["h_SelctiveElecCharge_CutFlow"]->Fill(0.,lumiWeight);
                for (unsigned int iElec = 0 ; iElec<init_electrons.size(); iElec++)
                {
                    
                    nbofElecbeforeChargeCut++;
                    is3ChargeAgr = init_electrons[iElec]->isGsfCtfScPixChargeConsistent();
                    // cout << " the electron with index =  " << iElec << " has Pt =  " << init_electrons[iElec]->Pt() << " is3ChargeAgr  is =  " << is3ChargeAgr << endl;
                    if (!is3ChargeAgr)
                    {
                        // cout << " this event will be rmoved !!!  " << endl;
                        nbElecWithNo3ChargeAgr++;
                        // break;
                    }
                }
                if (nbElecWithNo3ChargeAgr >0 )
                {
                    continue;
                }
                nbofEventsAfterChargeCut ++;
                histo1D["h_EventsWith3ChargeAgr"]->Fill(nbofEventsAfterChargeCut);
                histo1D["h_SelctiveElecCharge_CutFlow"]->Fill(1.,lumiWeight);
            }
            
          
            //Declare selection instance
            
            Run2Selection selection(init_jets, init_muons, init_electrons, mets, rho); ///rho = event->fixedGridRhoFastjetAll(); and it is used in electron isolation
            
            
            //// choose good primary vertex
            bool isGoodPV = selection.isPVSelected(vertex, 4 , 24. ,2.); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
            
            // Jets Selection
            selectedOrgiJets = selection.GetSelectedJets(jet_pt_cut, jet_eta_cut, true , "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
            
            
            /// --- Muons Selection -- ///
            selectedOrgiMuons = selection.GetSelectedMuons(mu_pt_cut , mu_eta_cut , mu_iso_cut ,"Tight","Summer16");// GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
            selectedOrgiLooseMuons = selection.GetSelectedMuons(mu_pt_cut , mu_eta_cut , mu_iso_cut ,"Loose","Summer16");
            
            //// --- Electron Selection --- ///
            selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut , "Tight" , "Spring16_80X", true, true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased , bool applyVID) // VID (i.e., a boolean which states whether or not the electron satisfies the officially recommended cuts)
            
            selectedLooseElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut , "Loose" , "Spring16_80X", true, true);
            
           if(debug)cout << "the nb of selected electrons is " << selectedElectrons.size()<<endl;
           if(debug) cout << "the number of events before jet cleaning is =  " << datasets[d]->NofEvtsToRunOver() << endl;

            //// sorting objects in the the event according to Pt
            
            sort(selectedJets.begin(), selectedJets.end(),HighestPt());
            sort(selectedOrgiMuons.begin(), selectedOrgiMuons.end(), HighestPt());
            sort(selectedOrgiLooseMuons.begin(), selectedOrgiLooseMuons.end(), HighestPt());
            sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt());
            sort(selectedLooseElectrons.begin(), selectedLooseElectrons.end(), HighestPt());
            
            ///// Rmoving events where the charge of the electron not calculate by selective method
//            if (Elec_Elec)
//            {
//                int nbElecWithNo3ChargeAgr = 0;
//                nbofEventsbeforeChargeCut++;
//                for (unsigned int iElec = 0 ; iElec<selectedElectrons.size(); iElec++)
//                {
//                    
//                    nbofElecbeforeChargeCut++;
//                    is3ChargeAgr = selectedElectrons[iElec]->isGsfCtfScPixChargeConsistent();
//                    // cout << " the electron with index =  " << iElec << " has Pt =  " << init_electrons[iElec]->Pt() << " is3ChargeAgr  is =  " << is3ChargeAgr << endl;
//                     if (!is3ChargeAgr)
//                     {
//                     // cout << " this event will be rmoved !!!  " << endl;
//                         nbElecWithNo3ChargeAgr++;
//                    // break;
//                    }
//                }
//                if (nbElecWithNo3ChargeAgr >0 )
//                {
//                    continue;
//                }
//             nbofEventsAfterChargeCut ++;
//            }
//            
            ///// ----- Removing bad and duplicated Muons ---- /////
            if ((Mu_Mu ||Elec_Mu) && RemoveBadMuon)
            {
                bool toBeErased = false;
                for (unsigned int imu=0; imu < selectedOrgiMuons.size(); imu++)
                {
                    bool toBeErased = false;
                    isBadPfMuon= selectedOrgiMuons[imu]->isBad80X();
                    isDuplicatedMuon = selectedOrgiMuons[imu]->isClone80X();
                    if (isBadPfMuon || isDuplicatedMuon)
                    {
                        toBeErased = true;
                        break;
                    }
                    if (!toBeErased)
                    {
                        selectedMuons.push_back(selectedOrgiMuons[imu]);
                    }
                }
                for (unsigned int imu=0; imu < selectedOrgiLooseMuons.size(); imu++)
                {
                    bool toBeErased = false;
                    isBadPfMuon= selectedOrgiLooseMuons[imu]->isBad80X();
                    isDuplicatedMuon = selectedOrgiLooseMuons[imu]->isClone80X();
                    if (isBadPfMuon || isDuplicatedMuon)
                    {
                        toBeErased = true;
                        break;
                    }
                    if (!toBeErased)
                    {
                        selectedLooseMuons.push_back(selectedOrgiLooseMuons[imu]);
                    }
                }
                
                if(debug)
                {
                    if( selectedOrgiMuons.size() != selectedMuons.size()) cout << "--> original Muon collection size = " << selectedOrgiMuons.size()  << "  And after removing Bad Muons =  " << selectedMuons.size() << endl;
                    else cout << "--> no change" << endl;
                    if( selectedOrgiLooseMuons.size() != selectedLooseMuons.size()) cout << "--> original Loose Muon collection size = " << selectedOrgiLooseMuons.size()  << "  And after removing Bad Muons =  " << selectedLooseMuons.size() << endl;
                    else cout << "--> no change" << endl;
                }

                

            }
            sort(selectedMuons.begin(), selectedMuons.end(), HighestPt());
            sort(selectedLooseMuons.begin(), selectedLooseMuons.end(), HighestPt());
            
            
            
            ///// --- Jet Cleaning -- ////
            
            /////// New mthod recommended from Lana & Liese
            
            selectedJets.clear();
            
            if (Apply_JetCleaning)
            {
                if(debug) cout << "Applying jet cleaning " << endl;
                for (int iOrgiJets=0 ; iOrgiJets < selectedOrgiJets.size() ; iOrgiJets++)
                {
                    bool toBeErased = false;
                    for (int iMuon = 0 ; iMuon < selectedMuons.size(); iMuon++)
                    {
                        if (selectedOrgiJets[iOrgiJets]->DeltaR(*selectedMuons[iMuon]) < 0.4)
                        {
                            toBeErased = true;
                            break;
                        }
                    }
                    for (int iElecton = 0 ; iElecton < selectedElectrons.size(); iElecton++)
                    {
                        if (selectedOrgiJets[iOrgiJets]->DeltaR(*selectedElectrons[iElecton]) < 0.4)
                        {
                            toBeErased = true;
                            break;
                        }
                    }
                    if(!toBeErased)
                    {
                        selectedJets.push_back(selectedOrgiJets[iOrgiJets]);
                    }
                    
                }
                if(debug)
                {
                    if( selectedOrgiJets.size() != selectedJets.size()) cout << "--> original jet collection size = " << selectedOrgiJets.size()  << "  And after cleaning =  " << selectedJets.size() << endl;
                    else cout << "--> no change" << endl;
                }
                
            }
    
           if(debug)cout << "the number of events after jet cleaning is =  " << datasets[d]->NofEvtsToRunOver() << endl;
        

            
            /////---- bTagging & c-tagging----\\\\\\
            
            //cout << "the size of BCSVLJets before filling and before clear = " << selectedBCSVLJets.size() << endl;
            selectedBCSVLJets.clear();
            selectedBCSVMJets.clear();
            selectedBCSVTJets.clear();
            selectedNonBCSVLJets.clear();
            selectedNonBCSVMJets.clear();
            selectedNonBCSVTJets.clear();
            selectedCSVLLJets.clear();
            selectedCSVMLJets.clear();
            selectedCSVTLJets.clear();
           if(debug)cout << "the size of BCSVLJets before filling and before any selection = " << selectedBCSVLJets.size() << endl;
            
            for (unsigned ibTagjet = 0; ibTagjet < selectedJets.size(); ibTagjet++)
            {
                TRootJet* tempJet = (TRootJet*) selectedJets[ibTagjet];
                if (tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Loose)
                {
                    selectedBCSVLJets.push_back(tempJet);
                    if (tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Medium)
                    {
                        selectedBCSVMJets.push_back(tempJet);
                        if (tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > CSVv2_workingpointvalue_Tight)
                        {
                            selectedBCSVTJets.push_back(tempJet);
                        } else
                        {
                            selectedNonBCSVTJets.push_back(tempJet);
                        }
                    } else
                    {
                        selectedNonBCSVMJets.push_back(tempJet);
                        selectedNonBCSVTJets.push_back(tempJet);
                    }
                    
                } else
                {
                    selectedNonBCSVMJets.push_back(tempJet);
                    selectedNonBCSVTJets.push_back(tempJet);
                    selectedNonBCSVLJets.push_back(tempJet);
                  
                }
                
            }
            

           // cout << "the nb of BCSVLJets After filling and before any selection = " << selectedBCSVLJets.size() << endl;
            if(verbose > 2) cout << "btagging done" << endl;
            /// sorting bTag Jets
            sort(selectedBCSVLJets.begin(),selectedBCSVLJets.end(),HighestPt());
            sort(selectedBCSVMJets.begin(),selectedBCSVMJets.end(),HighestPt());
            sort(selectedBCSVTJets.begin(),selectedBCSVTJets.end(),HighestPt());
            sort(selectedNonBCSVLJets.begin(),selectedNonBCSVLJets.end(),HighestPt());
            sort(selectedNonBCSVMJets.begin(),selectedNonBCSVMJets.end(),HighestPt());
            sort(selectedNonBCSVTJets.begin(),selectedNonBCSVTJets.end(),HighestPt());
            
            
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
            
            //////****************//////
            ///// btagging SF  /////////
            //////****************//////
            float btagWeight  =  1.;
            double btagWeightShape = 1.;
            Double_t btagWeight_shape_up_lf= 1.;
            Double_t btagWeight_shape_down_lf= 1.;
            Double_t btagWeight_shape_up_hf= 1.;
            Double_t btagWeight_shape_down_hf= 1.;
            Double_t btagWeight_shape_up_hfstats1= 1.;
            Double_t btagWeight_shape_down_hfstats1= 1.;
            Double_t btagWeight_shape_up_hfstats2= 1.;
            Double_t btagWeight_shape_down_hfstats2= 1.;
            Double_t btagWeight_shape_up_lfstats1= 1.;
            Double_t btagWeight_shape_down_lfstats1= 1.;
            Double_t btagWeight_shape_up_lfstats2= 1.;
            Double_t btagWeight_shape_down_lfstats2= 1.;
            Double_t btagWeight_shape_up_cferr1= 1.;
            Double_t btagWeight_shape_down_cferr1= 1.;
            Double_t btagWeight_shape_up_cferr2= 1.;
            Double_t btagWeight_shape_down_cferr2= 1.;
            float bTagEff = 1, bTagEff_LFUp = 1, bTagEff_LFDown = 1, bTagEff_HFUp = 1, bTagEff_HFDown = 1, bTagEff_HFStats1Up = 1,
            bTagEff_HFStats1Down = 1, bTagEff_HFStats2Up = 1, bTagEff_HFStats2Down = 1, bTagEff_LFStats1Up = 1, bTagEff_LFStats1Down = 1,
            bTagEff_LFStats2Up = 1, bTagEff_LFStats2Down = 1, bTagEff_CFErr1Up = 1, bTagEff_CFErr1Down = 1, bTagEff_CFErr2Up = 1, bTagEff_CFErr2Down = 1;
            
            ////// non shape b_tag scaling factor //////
            if(fillBtagHisto && !isData && !btagShape)
            {
                btwt_CSVv2L_mujets_central->FillMCEfficiencyHistos(selectedJets);
                
            }
            else if(!fillBtagHisto && !isData && !btagShape)
            {
                btagWeight =  btwt_CSVv2L_mujets_central->getMCEventWeight(selectedJets);
               
            }
            
            ////// shape b_tag scaling factor //////
            
            if(!isData && btagShape)
            {
                
                double jetpt ;
                double jeteta;
                double jetdisc ;
                int jetpartonflav ;
                bool isBFlav = false;
                bool isLFlav = false;
                bool isCFlav = false;
                
                for(int intJet = 0; intJet < selectedJets.size(); intJet++)
                {
                    jetpt = selectedJets[intJet]->Pt();
                    if(jetpt > 1000.) jetpt = 999.;
                    jeteta = selectedJets[intJet]->Eta();
                    jetdisc = selectedJets[intJet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                    if(jetdisc<0.0) jetdisc = -0.05;
                    if(jetdisc>1.0) jetdisc = 1.0;
                    isBFlav = false;
                    isLFlav = false;
                    isCFlav = false;
                    BTagEntry::JetFlavor jflav;
                    jetpartonflav = std::abs(selectedJets[intJet]->partonFlavour());
                    if(debug) cout<<"parton flavour: "<<jetpartonflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
                    if(jetpartonflav == 5)
                    {
                        jflav = BTagEntry::FLAV_B;
                        isBFlav =true;
                    }
                    else if(jetpartonflav == 4)
                    {
                        jflav = BTagEntry::FLAV_C;
                        isCFlav=true;
                    }
                    else
                    {
                        jflav = BTagEntry::FLAV_UDSG;
                        isLFlav = true;
                    }
                    
                    if( doJESShift == 2 && !isCFlav)        bTagEff = reader_JESUp->eval(jflav, jeteta, jetpt, jetdisc);
                    else if( doJESShift == 1 && !isCFlav) bTagEff = reader_JESDown->eval(jflav, jeteta, jetpt, jetdisc);
                    else bTagEff = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    
                    
                    ///// Other Systematics
                    
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
                    if( !isBFlav ) bTagEff_LFUp = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_LFDown = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_HFUp = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_HFDown = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats1Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats1Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats2Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isBFlav ) bTagEff_HFStats2Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats1Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats1Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats2Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isLFlav ) bTagEff_LFStats2Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr1Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr1Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr2Up = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    if( !isCFlav ) bTagEff_CFErr2Down = b_reader_CSVv2_shape->eval(jflav, jeteta, jetpt, jetdisc);
                    
                    
                    // fill btagweights
                    btagWeightShape *= bTagEff;
                    
                    btagWeight_shape_up_lf *= bTagEff_LFUp;
                    btagWeight_shape_down_lf *= bTagEff_LFDown;
                    btagWeight_shape_up_hf *= bTagEff_HFUp;
                    btagWeight_shape_down_hf *= bTagEff_HFDown;
                    btagWeight_shape_up_hfstats1 *= bTagEff_HFStats1Up;
                    btagWeight_shape_down_hfstats1 *= bTagEff_HFStats1Down;
                    btagWeight_shape_up_hfstats2 *= bTagEff_HFStats2Up;
                    btagWeight_shape_down_hfstats2 *= bTagEff_HFStats2Down;
                    btagWeight_shape_up_lfstats1 *= bTagEff_LFStats1Up;
                    btagWeight_shape_down_lfstats1 *= bTagEff_LFStats1Down;
                    btagWeight_shape_up_lfstats2 *= bTagEff_LFStats2Up;
                    btagWeight_shape_down_lfstats2 *= bTagEff_LFStats2Down;
                    btagWeight_shape_up_cferr1 *= bTagEff_CFErr1Up;
                    btagWeight_shape_down_cferr1 *= bTagEff_CFErr1Down;
                    btagWeight_shape_up_cferr2 *= bTagEff_CFErr2Up;
                    btagWeight_shape_down_cferr2 *= bTagEff_CFErr2Down;
                    
                }
                
            }
            
            ////===========================///////
            //// Applying Scale Factors    //////
            ////==========================//////
            
            // scale factor for the event
            Double_t scaleFactor = 1.;
            Double_t cutstep = 1.;
            Double_t Elec_scaleFactor = 1.;
            Double_t Muon_scaleFactor = 1.;
            Double_t MuonID_scaleFactor = 1.;
            Double_t MuonIso_scaleFactor = 1.;
            Double_t PU_scaleFactor =1.;
            Double_t Lep_scaleFactor = 1.;

            btagSFshape = 1.;
            btagSFshape_down_cferr1 = 1.;
            btagSFshape_down_cferr2 = 1.;
            btagSFshape_down_hf= 1.;
            btagSFshape_down_hfstats1 = 1.;
            btagSFshape_down_hfstats2 = 1.;
            btagSFshape_down_lf = 1.;
            btagSFshape_down_lfstats1 = 1.;
            btagSFshape_down_lfstats2 = 1.;
            btagSFshape_up_cferr1 = 1.;
            btagSFshape_up_cferr2 = 1.;
            btagSFshape_up_hf= 1.;
            btagSFshape_up_hfstats1 = 1.;
            btagSFshape_up_hfstats2 = 1.;
            btagSFshape_up_lf = 1.;
            btagSFshape_up_lfstats1 = 1.;
            btagSFshape_up_lfstats2 = 1.;
            float muon1SF, muon2SF,muonID1SF, muonID2SF,muonIso1SF, muonIso2SF,electron1SF, electron2SF;
            muon1SF =  muon2SF = electron1SF = electron2SF = 1.;
            
            histo1D["h_cutFlow"]->Fill(0., scaleFactor*lumiWeight); //// fill histogram before applying any scaling factors or triggers
             nbEvents_0++; //// add to the frist bin of Count_cut branch before applying any scaling factors or triggers
            //nCuts++;
           // Count_cut[0]=Count_cut[0]+cutstep;
            
            
            /////applyNegWeightCorrection
            
            if(hasNegWeight && applyNegWeightCorrection && !isData) scaleFactor *= -1.;
            
            ////// PU SF
            
            float puWeight = 1;
            float puWeight_UP = 1;
            float puWeight_Down = 1;
            
            // if (verbose>2) cout << " isData =  "<< isData << endl;
            if (ApplyPU_SF)
            {
                
                if(!isData)
                {
                    puWeight = LumiWeights.ITweight((int)event->nTruePU()); // simplest reweighting, just use reconstructed number of PV. faco
                    puWeight_UP = LumiWeights_UP.ITweight((int)event->nTruePU());
                    puWeight_Down = LumiWeights_Down.ITweight((int)event->nTruePU());
                    PU_scaleFactor =puWeight;
                    if (debug) cout << "while puWeight  =  "<< puWeight << endl;
                }
                else if (isData) {PU_scaleFactor =1;}
                 if (debug) cout << "PU_scaleFactor is " << PU_scaleFactor << endl;
                scaleFactor *= PU_scaleFactor;
                
            }
           nbEvents_1PU++; //// add to the frist bin of Count_cut branch After applying PU SF
            if (debug) cout << "puSF " << PU_scaleFactor << endl;
            
            /////////////////////////////////
            //// *** Apply selections *** ///
            ////////////////////////////////
            nbEvents++;
           // cout << "Hello : before any selections =  " << endl;
            bool diElectron = false;
            bool diMuon = false;
            bool diEMu = false;
            bool diMuE = false;
            bool SSdiLepton = false;
            bool OSdiLepton = false;
           // bool Event_passed = false;
            float InvMass_ll = 0.;
            float InvMass_lb = 0.;
            const float Zmass = 91.1876; // GeV SM xsections at 13 TeV https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
            const float Wmass = 80.398; // GeV SM xsections at 13 TeV https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
            vector <float> InvMass_JJ;
            TLorentzVector tempLepton_0;
            TLorentzVector tempLepton_1;
            TLorentzVector Lepton0;
            TLorentzVector Lepton1;
            TLorentzVector SM_bJet;
            int qLepton0 , qLepton1, qElec0, qElec1, qMu0, qMu1;
            qLepton0 = qLepton1 = qElec0 = qElec1 =  qMu0 = qMu1 = 0;
            float sum_jet_PT = 0.;
            float Ht_AllJets = 0.;
            float St_AllJets_Leps = 0.;
            float Sum_Leptons_Pt =0.;
            float Zmass_Window = 0.;
            DeltaR_2L = 0.;
            DeltaPhi_2L = 0.;
            invMass_2L = 0.;
           // invMass_2SSL_Zmass = 0;
          //  invMass_2OSL_Zmass =0;
            DeltaR_Mu0b0_DiMu= DeltaR_Mu1b0_DiMu= DeltaR_Mu0b0_DiMuEl= DeltaR_Mu1b0_DiElMu= DeltaR_Elec0b0_DiElec= DeltaR_Elec1b0_DiElec = DeltaR_Elec0b0_DiElMu = DeltaR_Elec1b0_DiMuEl =0.;
        
            Mass_WJetPair = 0.;
            Mass_JetPair = 0.;
            double massZLepPair = 0.;
            bool SSLinZpeak = false;
            bool OSLinZpeak = false;
            bool inZpeak = false;
            
            ///// ***** Selecting Good Primary Vertex ***/////
            if (!isGoodPV) continue;
            // nCuts++;
            
            // Count_cut[4]=Count_cut[4]+cutstep;
            if(debug)cout << "the nCut Value After good PV is =  " << nCuts << "and the nbEvents (nbEvents_2GPV) =  " << nbEvents_2GPV << endl;
            if(verbose>2) cout << "GoodPV" << endl;
            histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
            
            ///// *** Apply Met Filters *** ///
            if(HBHEnoise && HBHEIso && CSCTight && EcalDead  && isGoodPV && badchan && badmu) passedMET = true;
            PassedMETFilter = passedMET;
            
            if (!PassedMETFilter) continue;
            nbEvents_2GPV++;

            
            ///// ***** Applying Triggers ***/////
            if(Apply_HLT_Triggers) trigged = itrigger; //trigged = treeLoader.EventTrigged(itrigger);
            if (!trigged) continue;
            nbEvents_3Trig++;
            histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
            
           if(debug) cout << "the number of events After Applying trigger is =  " << datasets[d]->NofEvtsToRunOver() << endl;
           if(debug) cout << "trigged = " << trigged << endl;
            
           // nCuts++;
          //  Count_cut[3]=Count_cut[3]+cutstep;
            
            //////--- Filling variable branches for Met  ----- //////
            
            float met_px = mets[0]->Px();
            float met_py = mets[0]->Py();
            
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
            met_Phi = mets[0]->Phi();
            met_Eta = mets[0]->Eta();
            
            if(!isData)
            {
                puSF =  puWeight;
                puSF_UP= puWeight_UP;
                puSF_Down= puWeight_Down;
                
                btagSFshape = btagWeightShape;
                btagSFshape_up_lf = btagWeight_shape_up_lf;
                btagSFshape_down_lf = btagWeight_shape_down_lf;
                btagSFshape_up_hf = btagWeight_shape_up_hf;
                btagSFshape_down_hf = btagWeight_shape_down_hf;
                btagSFshape_up_hfstats1 = btagWeight_shape_up_hfstats1;
                btagSFshape_down_hfstats1 = btagWeight_shape_down_hfstats1;
                btagSFshape_up_hfstats2 = btagWeight_shape_up_hfstats2;
                btagSFshape_down_hfstats2 = btagWeight_shape_down_hfstats2;
                btagSFshape_up_lfstats1 = btagWeight_shape_up_lfstats1;
                btagSFshape_down_lfstats1 = btagWeight_shape_down_lfstats1;
                btagSFshape_up_lfstats2 = btagWeight_shape_up_lfstats2;
                btagSFshape_down_lfstats2 = btagWeight_shape_down_lfstats2;
                btagSFshape_up_cferr1 = btagWeight_shape_up_cferr1;
                btagSFshape_down_cferr1 = btagWeight_shape_down_cferr1;
                btagSFshape_up_cferr2 = btagWeight_shape_up_cferr2;
                btagSFshape_down_cferr2 = btagWeight_shape_down_cferr2;
                
            }else if(isData)
            {
                btagSFshape = 1.;   puSF_Down = 1.; puSF_UP = 1.; puSF = 1.;
                
            }
            
            if (!isData &&!btagShape)
            {
                scaleFactor *= btagWeight;
            }else if (!isData && btagShape)
            {
                scaleFactor *= btagSFshape;
            }
            
            nbEvents_4BTag++;

            ////*** Applying Selection on leptons ***/////
            int numLep = selectedElectrons.size()+ selectedMuons.size();
            int numLooseLep = selectedLooseElectrons.size()+ selectedLooseMuons.size();
            nMuons = 0;
            nLeptons = 0;
            nElectrons = 0;
         if(debug)cout<< "the size of the muon cllection is =  " << selectedMuons.size() << " while  the size of the electron cllection is =  " <<selectedElectrons.size() <<endl;
            if (numLep==2 && numLooseLep == 2)
            {
               // histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
                
                ///Making Lorentz vectors for leptons & met
                selectedMuonsTLV.clear();
                selectedElectronsTLV.clear();
                selectedLeptonsTLV.clear();
                metsTLV.clear();
                
                for(unsigned int iMuon=0; iMuon<selectedMuons.size(); iMuon++)
                {
                    selectedMuonsTLV.push_back( *selectedMuons[iMuon]);  // vector filled with TLV of selected Muons
                }
                for(unsigned int iElectron=0; iElectron<selectedElectrons.size(); iElectron++)
                {
                    selectedElectronsTLV.push_back( *selectedElectrons[iElectron]);
                }
                if(mets.size() > 0) metsTLV.push_back( *mets[0] );
                else cout<<"No MET??"<<endl;

                //////--- Filling variable branches for Muons ----- //////
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
                    
                    ///// Applying muon scaling factors
                    
                    if(ApplyMu_SF && !isData)
                    {
                        
                        MuonIDSF[nMuons] = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsBCDEF+lum_RunsGH);
                        
                        MuonIsoSF[nMuons] =  (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsBCDEF + muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                        MuonTrackSF[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta());
                        
                        MuonIDSF_up[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                        MuonIsoSF_up[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                        
                        MuonIDSF_down[nMuons]  = (muonSFWeightID_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                        MuonIsoSF_down[nMuons] = (muonSFWeightIso_BCDEF->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsBCDEF+muonSFWeightIso_GH->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), -1)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF);
                        
                        
                        MuonTrackSF_up[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*1.01;
                        MuonTrackSF_down[nMuons] = h_muonSFWeightTrack->Eval(selectedMuons[selmu]->Eta())*0.99;
                        
                        Muon_scaleFactor *= MuonIDSF[nMuons] * MuonIsoSF[nMuons] * MuonTrackSF[nMuons];
                    }
                    else
                    {
                        
                        Muon_scaleFactor = 1.;
                        MuonIDSF[nMuons] = 1.;
                        MuonIsoSF[nMuons] = 1.;
                        MuonIDSF_up[nMuons] = 1.;
                        MuonIsoSF_up[nMuons] = 1.;
                        MuonIDSF_down[nMuons] = 1.;
                        MuonIsoSF_down[nMuons] = 1.;
                        MuonTrackSF[nMuons] = 1.;
                        MuonTrackSF_up[nMuons] = 1.;
                        MuonTrackSF_down[nMuons] = 1.;
                    }
                    
                    nMuons++;
                }
               // scaleFactor *= Muon_scaleFactor;
                
                //////--- Filling variable branches for Electrons ----- //////
                
                for (unsigned int selelec = 0 ; selelec < selectedElectrons.size() ; selelec++)
                {
                    

                    pt_electron[nElectrons]=selectedElectrons[selelec]->Pt();
                    phi_electron[nElectrons]=selectedElectrons[selelec]->Phi();
                    eta_electron[nElectrons]=selectedElectrons[selelec]->Eta();
                    eta_superCluster_electron[nElectrons]=selectedElectrons[selelec]->superClusterEta();
                    E_electron[nElectrons]=selectedElectrons[selelec]->E();
                    d0_electron[nElectrons]=selectedElectrons[selelec]->d0();
                    d0BeamSpot_electron[nElectrons]=selectedElectrons[selelec]->d0BeamSpot();
                    chargedHadronIso_electron[nElectrons]=selectedElectrons[selelec]->chargedHadronIso(3);
                    neutralHadronIso_electron[nElectrons]=selectedElectrons[selelec]->neutralHadronIso(3);
                    photonIso_electron[nElectrons]=selectedElectrons[selelec]->photonIso(3);
                    pfIso_electron[nElectrons]=selectedElectrons[selelec]->relPfIso(3,0);
                    charge_electron[nElectrons]=selectedElectrons[selelec]->charge();
                    sigmaIEtaIEta_electron[nElectrons]=selectedElectrons[selelec]->sigmaIEtaIEta();
                    deltaEtaIn_electron[nElectrons]=selectedElectrons[selelec]->deltaEtaIn();
                    deltaPhiIn_electron[nElectrons]=selectedElectrons[selelec]->deltaPhiIn();
                    hadronicOverEm_electron[nElectrons]=selectedElectrons[selelec]->hadronicOverEm();
                    missingHits_electron[nElectrons]=selectedElectrons[selelec]->missingHits();
                    passConversion_electron[nElectrons]=selectedElectrons[selelec]->passConversion();
                    isEBEEGap[nElectrons]=selectedElectrons[selelec]->isEBEEGap();
                    isSelectiveElec[nElectrons] =selectedElectrons[selelec]->isGsfCtfScPixChargeConsistent();
                    //cout<< "electron with index =  " << selelec << "  with Pt =  " << selectedElectrons[selelec]->Pt() << endl;
                    if (ApplyElec_SF && !isData)
                    {
                        ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
                        ElectronSF[nElectrons]= electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),0);
                        ElectronSF_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),1)*electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),1);
                        ElectronSF_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),-1)*electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),-1);
                        
                        ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
                        ElectronSFID[nElectrons] = electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),0);
                        ElectronSFID_up[nElectrons] = electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),1);
                        ElectronSFID_down[nElectrons] = electronSFWeightID->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),-1);
                        
                        ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
                        ElectronSFReco[nElectrons] = electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),0);
                        ElectronSFReco_up[nElectrons] = electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),1);
                        ElectronSFReco_down[nElectrons] = electronSFWeightReco->at(selectedElectrons[selelec]->superClusterEta(),selectedElectrons[selelec]->Pt(),-1);
                        
                    }
                    else {
                        ElectronSF[nElectrons] = 1.; ElectronSF_up[nElectrons] = 1.; ElectronSF_down[nElectrons] = 1.;
                        ElectronSFReco[nElectrons] = 1.; ElectronSFReco_up[nElectrons] = 1.; ElectronSFReco_down[nElectrons] = 1.;
                        ElectronSFID[nElectrons] = 1.; ElectronSFID_up[nElectrons] = 1.; ElectronSFID_down[nElectrons] = 1.;
                    }
                    
                    nElectrons++;
                    
                }
                nLeptons = nMuons + nElectrons;
                
                nbEvents_5LepSF++;
                
                histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
                
                //////--- Filling variable branches for Jets  ----- //////
                nJets = 0;
                
                for(Int_t seljet = 0; seljet < selectedJets.size(); seljet++)
                {
                    
                    pt_jet[nJets]=selectedJets[seljet]->Pt();
                    phi_jet[nJets]=selectedJets[seljet]->Phi();
                    eta_jet[nJets]=selectedJets[seljet]->Eta();
                    E_jet[nJets]=selectedJets[seljet]->E();
                    charge_jet[nJets]=selectedJets[seljet]->charge();
                    bdisc_jet[nJets]=selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                    if(selectedJets.size()>0) {pt_1st_jet =selectedJets[0]->Pt(); bdiscCSVv2_1st_jet = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();}
                    if(selectedJets.size()>1) {pt_2nd_jet =selectedJets[1]->Pt();bdiscCSVv2_2nd_jet = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();}
                    if(selectedJets.size()>2) {pt_3rd_jet =selectedJets[2]->Pt();bdiscCSVv2_3rd_jet = selectedJets[2]->btag_combinedInclusiveSecondaryVertexV2BJetTags();}
                    if(selectedJets.size()>3) {pt_4th_jet =selectedJets[3]->Pt();bdiscCSVv2_3rd_jet = selectedJets[3]->btag_combinedInclusiveSecondaryVertexV2BJetTags();}
                    nJets++;
                    
                }
                MVA_nJets = nJets;
                
                
                //////--- Filling variable branches for b-tagged Jets  ----- //////
                nCSVTbJets = selectedBCSVTJets.size();
                nCSVMbJets = selectedBCSVMJets.size();
                nCSVLbJets = selectedBCSVLJets.size();
                MVA_nCSVLbtagJets = nCSVLbJets;
               // bdiscCSVv2_1st_bjet = bdiscCSVv2_2nd_bjet = 0;
                for (unsigned int selbjet = 0; selbjet < selectedBCSVLJets.size(); selbjet++)
                {
                    if(selectedBCSVLJets.size()>0 ) bdiscCSVv2_1st_bjet = selectedBCSVLJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                    if(selectedBCSVLJets.size()>1 ) bdiscCSVv2_2nd_bjet = selectedBCSVLJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                }

                //////// **** Divideing Analysis iinto 3 channels ee - uu - eu **** ////
                if (Elec_Elec && !Mu_Mu && !Elec_Mu) // Electron-Electron Channel
                {
                    
                    if (selectedElectrons.size()==2 && selectedMuons.size()==0 && selectedElectrons[0]->Pt()>25. && selectedElectrons[1]->Pt()>= 20.)
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
                            invMass_2L = InvMass_ll;
                            electron1SF = electronSFWeightID->at(selectedElectrons[0]->superClusterEta(),selectedElectrons[0]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[0]->superClusterEta(),selectedElectrons[0]->Pt(),0);
                            electron2SF = electronSFWeightID->at(selectedElectrons[1]->superClusterEta(),selectedElectrons[1]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[1]->superClusterEta(),selectedElectrons[1]->Pt(),0);
                            Lep_scaleFactor  = electron1SF * electron2SF;
                            Elec_scaleFactor = electron1SF*electron2SF;
                           // scaleFactor *= Elec_scaleFactor;
                            if(selectedElectrons.size()>0) pt_1st_Electron = selectedElectrons[0]->Pt();
                            if(selectedElectrons.size()>1) pt_2nd_Electron = selectedElectrons[1]->Pt();
                            if(selectedElectrons.size()>0) eta_1st_Electron = selectedElectrons[0]->Eta();
                            if(selectedElectrons.size()>1) eta_2nd_Electron = selectedElectrons[1]->Eta();
                            if(selectedElectrons.size()>0) phi_1st_Electron = selectedElectrons[0]->Phi();
                            if(selectedElectrons.size()>1) phi_2nd_Electron = selectedElectrons[1]->Phi();
                            if(selectedElectrons.size()>0) charge_1st_Electron = selectedElectrons[0]->charge();
                            if(selectedElectrons.size()>1) charge_2nd_Electron = selectedElectrons[1]->charge();
                            MT_lep1 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_lep2 = sqrt(2*selectedElectrons[1]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[1].DeltaPhi( metsTLV[0] )) ) );
                            MT_Elec1 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_Elec2 = sqrt(2*selectedElectrons[1]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[1].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedElectrons[0]);
                            selectedLeptonsTLV.push_back(*selectedElectrons[1]);
                            DeltaPhi_met_lep1 = selectedElectronsTLV[0].DeltaPhi(metsTLV[0]);
                            DeltaPhi_met_lep2 = selectedElectronsTLV[1].DeltaPhi(metsTLV[0]);
                            
                            if (fabs(Zmass - InvMass_ll)<= 10 )
                            {
                                massZLepPair= InvMass_ll;
                                inZpeak = true;
                            }

                            
                        }
                      
                    }
                    
                }
                else if (Mu_Mu && !Elec_Elec && !Elec_Mu) // Muon-Muon channel
                {
                    
                    if (selectedMuons.size()==2 && selectedElectrons.size()==0 && selectedMuons[0]->Pt() >= 20. && selectedMuons[1]->Pt()>= 15.)
                    {
                        
                        tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                        qMu0 = selectedMuons[0]->charge();
                        tempLepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
                        qMu1 = selectedMuons[1]->charge();
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        
                        
                        if (InvMass_ll >12.)
                        {
                            diMuon = true;
                            invMass_2L = InvMass_ll;
                            muon1SF = ((muonSFWeightID_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsBCDEF+lum_RunsGH))*((muonSFWeightIso_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF + muonSFWeightIso_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF))*(h_muonSFWeightTrack->Eval(selectedMuons[0]->Eta()));
                            muon2SF = ((muonSFWeightID_BCDEF->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0)*lum_RunsGH)/(lum_RunsBCDEF+lum_RunsGH))*((muonSFWeightIso_BCDEF->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0)*lum_RunsBCDEF + muonSFWeightIso_GH->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF))*(h_muonSFWeightTrack->Eval(selectedMuons[1]->Eta()));
                           // scaleFactor *= (muon1SF * muon2SF);
                            Lep_scaleFactor *= (muon1SF * muon2SF);
                            if(selectedMuons.size()>0) pt_1st_Muon = selectedMuons[0]->Pt();
                            if(selectedMuons.size()>1) pt_2nd_Muon = selectedMuons[1]->Pt();
                            if(selectedMuons.size()>0) phi_1st_Muon = selectedMuons[0]->Phi();
                            if(selectedMuons.size()>1) phi_2nd_Muon = selectedMuons[1]->Phi();
                            if(selectedMuons.size()>0) eta_1st_Muon = selectedMuons[0]->Eta();
                            if(selectedMuons.size()>1) eta_2nd_Muon = selectedMuons[1]->Eta();
                            if(selectedMuons.size()>0) charge_1st_Muon = selectedMuons[0]->charge();
                            if(selectedMuons.size()>1) charge_2nd_Muon = selectedMuons[1]->charge();
                            MT_lep1 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_lep2 = sqrt(2*selectedMuons[1]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[1].DeltaPhi( metsTLV[0] )) ) );
                            MT_Mu1 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_Mu2 = sqrt(2*selectedMuons[1]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[1].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedMuons[0]);
                            selectedLeptonsTLV.push_back(*selectedMuons[1]);
                            DeltaPhi_met_lep1 = selectedMuonsTLV[0].DeltaPhi( metsTLV[0] );
                            DeltaPhi_met_lep2 = selectedMuonsTLV[1].DeltaPhi( metsTLV[0] );
                            
                            if (fabs(Zmass - InvMass_ll)<= 10 )
                            {
                                massZLepPair= InvMass_ll;
                                inZpeak = true;
                            }
                            
                        }
                        
                        
                    }
                }
                else if (Elec_Mu && !Elec_Elec && !Mu_Mu) //// Electron- Muon channel Emu or muE
                {
                   // cout<< "the size of the muon cllection is =  " << selectedMuons.size() << " while  the size of the electron cllection is =  " <<selectedElectrons.size() <<endl;
                    if (selectedMuons.size() == 1 && selectedElectrons.size() ==1)
                    {
                      //  cout << "Hello I am at EMu Before PT cut " <<endl;
                        if (selectedMuons[0]->Pt() > selectedElectrons[0]->Pt() && selectedMuons[0]->Pt()>= 25 && selectedElectrons[0]->Pt()>=20)
                        {
                          //  cout << "Hello I am at EMu After PT cut mu pt > e pT" <<endl;
                            diMuE = true;
                            tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            qLepton0 = selectedMuons[0]->charge();
                            qLepton1 = selectedElectrons[0]->charge();
                            selectedLeptonsTLV.push_back(*selectedMuons[0]);
                            selectedLeptonsTLV.push_back(*selectedElectrons[0]);
                            MT_lep1 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos(selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_lep2 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos(selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                           
                            DeltaPhi_met_lep1 = selectedMuonsTLV[0].DeltaPhi( metsTLV[0] );
                            DeltaPhi_met_lep2 = selectedElectronsTLV[0].DeltaPhi(metsTLV[0]);
                            
                            if(selectedMuons.size()>0) pt_1st_Muon = selectedMuons[0]->Pt();
                            if(selectedMuons.size()>0) phi_1st_Muon = selectedMuons[0]->Phi();
                            if(selectedMuons.size()>0) eta_1st_Muon = selectedMuons[0]->Eta();
                            if(selectedMuons.size()>0) charge_1st_Muon = selectedMuons[0]->charge();
                            
                            if(selectedElectrons.size()>0) pt_2nd_Electron = selectedElectrons[0]->Pt();
                            if(selectedElectrons.size()>0) eta_2nd_Electron = selectedElectrons[0]->Eta();
                            if(selectedElectrons.size()>0) phi_2nd_Electron = selectedElectrons[0]->Phi();
                            if(selectedElectrons.size()>0) charge_2nd_Electron = selectedElectrons[0]->charge();

                        }
                        if (selectedElectrons[0]->Pt() > selectedMuons[0]->Pt() && selectedElectrons[0]->Pt()>=25 && selectedMuons[0]->Pt()>= 15)
                        {
                           // cout << "Hello I am at EMu After PT cut mu pt < e pT" <<endl;
                            diEMu = true;
                            tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                            qLepton0 = selectedElectrons[0]->charge();
                            qLepton1 = selectedMuons[0]->charge();
                            selectedLeptonsTLV.push_back(*selectedElectrons[0]);
                            selectedLeptonsTLV.push_back(*selectedMuons[0]);
                            MT_lep1 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            MT_lep2 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            DeltaPhi_met_lep1 = selectedElectronsTLV[0].DeltaPhi(metsTLV[0]);
                            DeltaPhi_met_lep2 = selectedMuonsTLV[0].DeltaPhi( metsTLV[0] );
                            
                            if(selectedElectrons.size()>0) pt_1st_Electron = selectedElectrons[0]->Pt();
                            if(selectedElectrons.size()>0) eta_1st_Electron = selectedElectrons[0]->Eta();
                            if(selectedElectrons.size()>0) phi_1st_Electron = selectedElectrons[0]->Phi();
                            if(selectedElectrons.size()>0) charge_1st_Electron = selectedElectrons[0]->charge();
                            if(selectedMuons.size()>0) pt_2nd_Muon = selectedMuons[0]->Pt();
                            if(selectedMuons.size()>0) phi_2nd_Muon = selectedMuons[0]->Phi();
                            if(selectedMuons.size()>0) eta_2nd_Muon = selectedMuons[0]->Eta();
                            if(selectedMuons.size()>0) charge_2nd_Muon = selectedMuons[0]->charge();
                         
                        }
                        muon1SF = ((muonSFWeightID_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF+muonSFWeightID_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsBCDEF+lum_RunsGH))*((muonSFWeightIso_BCDEF->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsBCDEF + muonSFWeightIso_GH->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)*lum_RunsGH)/(lum_RunsGH+lum_RunsBCDEF))*(h_muonSFWeightTrack->Eval(selectedMuons[0]->Eta()));
                        electron1SF = electronSFWeightID->at(selectedElectrons[0]->superClusterEta(),selectedElectrons[0]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[0]->superClusterEta(),selectedElectrons[0]->Pt(),0);
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        Lep_scaleFactor = muon1SF * electron1SF;
                        
                    }
                }
               
                if (diElectron || diMuon || diEMu || diMuE )
                {
                   // if (met_Pt < 30 ) {continue;}
                    scaleFactor *= Lep_scaleFactor;
                   // nCuts++;
                    nbEvents_6++;
                  //  Count_cut[6]=Count_cut[6]+cutstep;
                    histo1D["h_cutFlow"]->Fill(4., scaleFactor*lumiWeight);
                   
                   ////// cout << "the nCut Value After passing 2 leptons is =  " << nCuts << "and the nbEvents (nbEvents_3) =  " << nbEvents_3 << endl;
                    
                    histo2D["2L_Nb_jets_vs_CSVLbjets"]->Fill(selectedJets.size(),selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                    histo2D["2L_Nb_jets_vs_CSVMbjets"]->Fill(selectedJets.size(),selectedBCSVMJets.size(),scaleFactor*lumiWeight);
                    histo2D["2L_Nb_jets_vs_CSVTbjets"]->Fill(selectedJets.size(),selectedBCSVTJets.size(),scaleFactor*lumiWeight);
                    histo2D["h_2L_Lep0_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                    
                    
                    //////////////////////////////////////////////////////////////////////
                    //////////////// Calculating Charge mis-Idetification ///////////////
                    ////////////////////////////////////////////////////////////////////
                    
                    
                   // if (isData && diElectron && fabs(Zmass - invMass_2L ) <= 10)
                   // cout << ""
                    if (diElectron && fabs(Zmass - invMass_2L ) <= 10 ) //&& !ApplyCharge_misID)
                    {
                        int eta_range, Pt_range;
                        eta_range = Pt_range = 0;
                        
                        if (qElec0 == qElec1)
                        {
                            for (unsigned iElec = 0 ; iElec < selectedElectrons.size() ; iElec++)
                            {
                                
                                //  if (selectedElectrons[iElec]->Pt() >=50) cout << " At SSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                               // cout << " At SSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                                
                                if (abs(selectedElectrons[iElec]->Eta()) > 0 && abs(selectedElectrons[iElec]->Eta()) <= 1.479)
                                {
                                    //eta_range = 1;
                                    if (selectedElectrons[iElec]->Pt() < 25 )
                                    {
                                        Nb_2SSl_B1++;
                                        //   cout << " Nb_2SSl_B1 =  " << Nb_2SSl_B1 <<endl;
                                        InBarrel = true;
                                        histo1D["h_NSSL_barrel_cutFlow"]->Fill(1);
                                        Pt_range =1;eta_range = 1;
                                    }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                    {
                                        Nb_2SSl_B2++;
                                        //   cout << " Nb_2SSl_B2 =  " << Nb_2SSl_B2 <<endl;
                                        InBarrel = true;
                                        histo1D["h_NSSL_barrel_cutFlow"]->Fill(2);
                                        Pt_range =2;eta_range = 1;
                                    }else if (selectedElectrons[iElec]->Pt() >= 50 )
                                    {
                                        Nb_2SSl_B3++;
                                        //    cout << " Nb_2SSl_B3 =  " << Nb_2SSl_B3 <<endl;
                                        InBarrel = true;
                                        histo1D["h_NSSL_barrel_cutFlow"]->Fill(3);
                                        Pt_range =3;eta_range = 1;
                                    }
                                }else if (abs(selectedElectrons[iElec]->Eta()) > 1.479 && abs(selectedElectrons[iElec]->Eta()) <=2.5)
                                {
                                    //eta_range = 2;
                                    if (selectedElectrons[iElec]->Pt() < 25 )
                                    {
                                        Nb_2SSl_E1++;
                                        InEndcap = true;
                                        histo1D["h_NSSL_EndCap_cutFlow"]->Fill(1);
                                        eta_range = 2;Pt_range =1;
                                    }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                    {
                                        Nb_2SSl_E2++;
                                        InEndcap = true;
                                        histo1D["h_NSSL_EndCap_cutFlow"]->Fill(2);
                                        eta_range = 2;Pt_range =2;
                                    }else //if (selectedElectrons[iElec]->Pt() >= 50)
                                    {   Nb_2SSl_E3++;
                                        InEndcap = true;
                                        histo1D["h_NSSL_EndCap_cutFlow"]->Fill(3);
                                        eta_range = 2;Pt_range =3;
                                    }
                                }
                                histo2D["h_NSSL_Pt_eta"]->Fill(Pt_range,eta_range);
                            }
                            
                        }
                        if (qElec0 != qElec1)
                        {
                            for (unsigned iElec = 0; iElec < selectedElectrons.size() ; iElec++)
                            {
                                //  histo1D["h_NOSL_barrel_cutFlow"]->Fill(0);
                                //  histo1D["h_NOSL_EndCap_cutFlow"]->Fill(0);
                                // if (selectedElectrons[iElec]->Pt() >=50) cout << " At OSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                                
                              //  cout << " At OSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                                
                                if (abs(selectedElectrons[iElec]->Eta()) > 0 && abs(selectedElectrons[iElec]->Eta()) <= 1.479)
                                {
                                    //eta_range = 1;
                                    if (selectedElectrons[iElec]->Pt() < 25 )
                                    {
                                        Nb_2OSl_B1++;
                                        InBarrel = true;
                                        //    cout << " Nb_2OSl_B1 =  " << Nb_2OSl_B1 <<endl;
                                        histo1D["h_NOSL_barrel_cutFlow"]->Fill(1);
                                        eta_range = 1; Pt_range = 1;
                                    }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                    {
                                        Nb_2OSl_B2++;
                                        InBarrel = true;
                                        //  cout << " Nb_2OSl_B2 =  " << Nb_2OSl_B2 <<endl;
                                        histo1D["h_NOSL_barrel_cutFlow"]->Fill(2);
                                        eta_range = 1; Pt_range = 2;
                                    }else //(selectedElectrons[iElec]->Pt()>= 50)
                                    {
                                        Nb_2OSl_B3++;
                                        InBarrel = true;
                                        //  cout << " Nb_2OSl_B3 =  " << Nb_2OSl_B3 <<endl;
                                        histo1D["h_NOSL_barrel_cutFlow"]->Fill(3);
                                        eta_range = 1; Pt_range = 3;
                                    }
                                }else if (abs(selectedElectrons[iElec]->Eta()) > 1.479 && abs(selectedElectrons[iElec]->Eta()) <=2.5)
                                {
                                    //eta_range = 2;
                                    if (selectedElectrons[iElec]->Pt() < 25 )
                                    {
                                        Nb_2OSl_E1++;
                                        InEndcap = true;
                                        histo1D["h_NOSL_EndCap_cutFlow"]->Fill(1);
                                        eta_range = 2; Pt_range =1;
                                    }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() <= 50 )
                                    {
                                        Nb_2OSl_E2++;
                                        InEndcap = true;
                                        histo1D["h_NOSL_EndCap_cutFlow"]->Fill(2);
                                        eta_range = 2; Pt_range=2;
                                    }else //(selectedElectrons[iElec]->Pt() >= 50 )
                                    {
                                        Nb_2OSl_E3++;
                                        InEndcap = true;
                                        histo1D["h_NOSL_EndCap_cutFlow"]->Fill(3);
                                        eta_range = 2; Pt_range =3;
                                    }
                                }
                                histo2D["h_NOSL_Pt_eta"]->Fill(Pt_range,eta_range);
                            }
                            
                        }
                        
                    }
                    

                   if (debug)cout << " Hello selectedJets.size = " << selectedJets.size() <<endl;
                    
                    if (selectedJets.size()>= 3 && selectedJets[0]->Pt() >=30 && selectedJets[1]->Pt() >=30 && selectedJets[2]->Pt() >=30)
                    {
                      //  nCuts++;
                        nbEvents_7++;
                       // Count_cut[7]=Count_cut[7]+cutstep;
                        histo1D["h_cutFlow"]->Fill(5., scaleFactor*lumiWeight);
                        
                        histo2D["2L_3Jets_Nb_jets_vs_CSVLbjets"]->Fill(selectedJets.size(),selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                        histo2D["2L_3Jets_Nb_jets_vs_CSVMbjets"]->Fill(selectedJets.size(),selectedBCSVMJets.size(),scaleFactor*lumiWeight);
                        histo2D["2L_3Jets_Nb_jets_vs_CSVTbjets"]->Fill(selectedJets.size(),selectedBCSVTJets.size(),scaleFactor*lumiWeight);
                        histo2D["h_2L_3Jets_Lep0_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                       
                   //     cout << " Hello selectedBCSVLJets.size = " << selectedBCSVLJets.size() <<endl;
                        if (selectedBCSVLJets.size()>= 1)
                        {
                           // nCuts++;
                            nbEvents_8++;
                           // Count_cut[8]=Count_cut[8]+cutstep;
                            histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
                            histo1D["h_cutFlow_LeptonCuts"]->Fill(0., scaleFactor*lumiWeight);
                            histo2D["h_2L_3Jets1b_Lep0_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                            
                            SM_bJet.SetPxPyPzE(selectedBCSVLJets[0]->Px(),selectedBCSVLJets[0]->Py(),selectedBCSVLJets[0]->Pz(),selectedBCSVLJets[0]->Energy());
                            
                            // TLorentzVector temp_W_1st_Jet;
                            //TLorentzVector temp_W_2nd_Jet;
                            
                            ///////////////////////////////
                            /// Chicking MC Information ///
                            ///////////////////////////////
                            
                            bool Find_FCNC_W = false;
                            bool Find_FCNC_W_1stJet = false;
                            bool Find_FCNC_W_2ndJet = false;
                            bool Find_FCNC_HdecayJets = false;
                            bool Find_SMNegLep = false;
                            bool Find_SMPosLep = false;
                            bool Find_SMNegElec = false;
                            bool Find_SMPosElec = false;
                            bool Find_SMNegMu = false;
                            bool Find_SMPosMu = false;
                            bool Find_SMPosB = false;
                            bool Find_SMNegB = false;
                            
                            bool Find_FCNCNegLep = false;
                            bool Find_FCNCPosLep = false;
                            bool Find_FCNCNegMu = false;
                            bool Find_FCNCPosMu = false;
                            bool Find_FCNCNegElec = false;
                            bool Find_FCNCPosElec = false;
                            bool Find_SM_top = false;
                            bool Find_SM_Antitop = false;
                            int topId = 6;
                            int atopId = -6;
                            
                            
                            selectedJetsTLV.clear();
                            JetsExcludingHighestCSVLbTLV.clear();
                            int FCNCWJetsCounter =0;
                            int FCNCWJetsCounter_nonb =0;
                            TLorentzVector temp_Match_WJet1;
                            TLorentzVector temp_Match_WJet2;
                            
                            for (unsigned int i = 0; i < selectedJets.size(); i++)
                            {
                                selectedJetsTLV.push_back(*selectedJets[i]);
                            }
                            for (unsigned int ijet = 0 ; ijet < JetsExcludingHighestCSVLb.size(); ijet++)
                            {
                                JetsExcludingHighestCSVLbTLV.push_back(*JetsExcludingHighestCSVLb[ijet]);
                            }
                            
                            for( unsigned int nj = 0; nj < selectedJets.size(); nj++)
                            {
                                Jet_matchMC_pdgId[nj] = 0;
                                Jet_matchMC_motherpdgId[nj] = 0;
                                Jet_matchMC_grannypdgId[nj] = 0;
                            }
                            
                            pair<unsigned int, unsigned int> FCNCWJet1_ = pair<unsigned int,unsigned int>(9999,9999); // First one is jet number, second one is mcParticle number
                            pair<unsigned int, unsigned int> FCNCWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
                            pair<unsigned int, unsigned int> FCNCWJet1_nonb_ = pair<unsigned int,unsigned int>(9999,9999); // First one is jet number, second one is mcParticle number
                            pair<unsigned int, unsigned int> FCNCWJet2_nonb_ = pair<unsigned int,unsigned int>(9999,9999);
                            
                            pair<unsigned int, unsigned int> FCNCLep_ = pair<unsigned int,unsigned int>(9999,9999);
                            pair<unsigned int, unsigned int> SMLep_ = pair<unsigned int,unsigned int>(9999,9999);
                            pair<unsigned int, unsigned int> SMbJet_ = pair<unsigned int,unsigned int>(9999,9999);
                            pair<unsigned int, unsigned int> Non1stbJet = pair<unsigned int,unsigned int>(9999,9999);
                            
                            
                            if (!isData)
                            {
                                treeLoader.LoadMCEvent(ievt, 0,  mcParticles, false);
                                sort(mcParticles.begin(),mcParticles.end(),HighestPt());
                                
                                mcParticlesMatching_.clear();
                                mcParticlesTLV.clear();
                                
                                vector< pair<unsigned int, unsigned int> > PPair; // First one is jet number, second one is mcParticle number
                                PPair.clear();
                                vector< pair<unsigned int, unsigned int> > W_1stJet_Pair; // First one is jet number, second one is mcParticle number
                                W_1stJet_Pair.clear();
                                vector< pair<unsigned int, unsigned int> > W_2ndJet_Pair; // First one is jet number, second one is mcParticle number
                                W_2ndJet_Pair.clear();
                                vector<string > NPair; // First one is jet number, second one is mcParticle number
                                NPair.clear();
                                
                                for (unsigned imc = 0; imc < mcParticles.size(); imc++)
                                {
                                    if ( (mcParticles[imc]->status() > 1 && mcParticles[imc]->status() <= 20) || mcParticles[imc]->status() >= 40 ) continue;  /// Final state particle or particle from hardest process
                                    if (mcParticles[imc]->status() == 1 && mcParticles[imc]->motherType() == -24 && mcParticles[imc]->grannyType() == -6)
                                    {
                                        if (mcParticles[imc]->type() == 13 || mcParticles[imc]->type() == 11)
                                        {
                                            Find_SMNegLep =true;
                                            if (mcParticles[imc]->type() == 13)
                                            {
                                                Find_SMNegMu = true;
                                            }else if (mcParticles[imc]->type() == 11)
                                            {
                                                Find_SMNegElec= true;
                                            }
                                        }
                                    }else if (mcParticles[imc]->status() == 1 && mcParticles[imc]->motherType() == 24 && mcParticles[imc]->grannyType() == 6)
                                    {
                                        if (mcParticles[imc]->type() == -13 || mcParticles[imc]->type() == -11)
                                        {
                                            Find_SMPosLep =true;
                                            if (mcParticles[imc]->type() == -13)
                                            {
                                                Find_SMPosMu = true;
                                            }else if (mcParticles[imc]->type() == -11)
                                            {
                                                Find_SMPosElec = true;
                                            }
                                        }
                                    }
                                    if (mcParticles[imc]->status() == 1 &&(mcParticles[imc]->grannyType() == 25 || mcParticles[imc]->motherType() == -24))
                                    {
                                        if (mcParticles[imc]->type() == 13 || mcParticles[imc]->type() == 11)
                                        {
                                            Find_FCNCNegLep = true;
                                            if (mcParticles[imc]->type() == 13)
                                            {
                                                Find_FCNCNegMu = true;
                                            }else if(mcParticles[imc]->type() == 11)
                                            {
                                                Find_FCNCNegElec = true;
                                            }
                                        }else if (mcParticles[imc]->type() == -13 || mcParticles[imc]->type() == -11)
                                        {
                                            Find_FCNCPosLep = true;
                                            if (mcParticles[imc]->type() == -13)
                                            {
                                                Find_FCNCPosMu = true;
                                            }else if(mcParticles[imc]->type() == -11)
                                            {
                                                Find_FCNCPosElec = true;
                                            }
                                        }
                                    }
                                    if (fabs(mcParticles[imc]->type()) <= 5 || fabs(mcParticles[imc]->type() == 21) )  //light/b quarks, 6 for top should stay hardcoded, OR gluon
                                    {
                                        mcParticlesTLV.push_back(*mcParticles[imc]);  // vector contains TLV of MC particles
                                        mcParticlesMatching_.push_back(mcParticles[imc]);
                                    }
                                    
                                }
                                
                                
                                JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
                                
                                if (matching.getNumberOfAvailableCombinations() != 1) cout << "matching.getNumberOfAvailableCombinations() = "<<matching.getNumberOfAvailableCombinations()<<" .  This should be equal to 1 !!!"<<endl;
                                
                                vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
                                JetPartonPair.clear();
                                
                                for (unsigned int i = 0; i < mcParticlesTLV.size(); i++)
                                {
                                    int matchedJetNumber = matching.getMatchForParton(i, 0);
                                    
                                    if (matchedJetNumber > -1) JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
                                }
                                if(debug)cout<< " the size of JetPartonPair vector is =  " << JetPartonPair.size() << endl;
                                
                                for (unsigned int iseljet = 0 ; iseljet < JetPartonPair.size(); iseljet++)
                                {
                                    
                                    unsigned int particlenb = JetPartonPair[iseljet].first; // place in selectedjets vector
                                    unsigned int partonIDnb = JetPartonPair[iseljet].second; // place in mcParticles_ vector
                                    
                                    Jet_matchMC_pdgId[particlenb]=mcParticlesMatching_[partonIDnb]->type();
                                    //cout << "the Id of corresponding MC particle is =  " << Jet_matchMC_pdgId[particlenb] << endl;
                                    Jet_matchMC_motherpdgId[particlenb]=mcParticlesMatching_[partonIDnb]->motherType();
                                    //cout << "the mother Id of corresponding MC particle is =  " << Jet_matchMC_motherpdgId[particlenb] << endl;
                                    Jet_matchMC_grannypdgId[particlenb]=mcParticlesMatching_[partonIDnb]->grannyType();
                                    //cout << "the granny Id of corresponding MC particle is =  " << Jet_matchMC_grannypdgId[particlenb] << endl;
                                    
                                    if (fabs(mcParticlesMatching_[partonIDnb]->type()) < topId)
                                    {
                                        // if ((Find_FCNCPosLep && Find_SMPosLep && mcParticlesMatching_[partonIDnb]->grannyType() == 25 || mcParticlesMatching_[partonIDnb]->motherType() == -24 && mcParticlesMatching_[partonIDnb]->grannyType() != -6) || (Find_FCNCNegLep &&Find_SMNegLep && mcParticlesMatching_[partonIDnb]->grannyType() == 25 || mcParticlesMatching_[partonIDnb]->motherType() == 24 && mcParticlesMatching_[partonIDnb]->grannyType() != topId ))
                                        if ((mcParticlesMatching_[partonIDnb]->grannyType() == 25 || mcParticlesMatching_[partonIDnb]->motherType() == -24 && mcParticlesMatching_[partonIDnb]->grannyType() != atopId) || (mcParticlesMatching_[partonIDnb]->grannyType() == 25 || mcParticlesMatching_[partonIDnb]->motherType() == 24 && mcParticlesMatching_[partonIDnb]->grannyType() != topId ))
                                        {
                                            cout << "Find particles comes from Higgs decay" << endl;
                                            if (FCNCWJetsCounter <1)
                                            {
                                                FCNCWJet1_ =JetPartonPair[particlenb];
                                                FCNCWJetsCounter++;
                                              //  cout << "Find FCNC W first Jet" << endl;
                                            }
                                            if (FCNCWJetsCounter ==1)
                                            {
                                                FCNCWJet2_ =JetPartonPair[particlenb];
                                                FCNCWJetsCounter++;
                                               // cout << "Find FCNC W second Jet" << endl;
                                            }
                                            
                                        }else if (mcParticlesMatching_[partonIDnb]->type() ==5 && mcParticlesMatching_[partonIDnb]->motherType() == topId)
                                        {
                                            Find_SMPosB = true;
                                            SMbJet_= JetPartonPair[particlenb];
                                            // cout << "Find SM +ve b Jet" << endl;
                                        }else if (mcParticlesMatching_[partonIDnb]->type() == -5 && mcParticlesMatching_[partonIDnb]->motherType() == atopId)
                                        {
                                            Find_SMNegB = true;
                                            SMbJet_ = JetPartonPair[particlenb];
                                            //  cout << "Find SM -ve b Jet" << endl;
                                        }
                                    }
                                    
                                    
                                    
                                }
                                
                                JetPartonMatching matching_non1stbJet = JetPartonMatching(mcParticlesTLV, JetsExcludingHighestCSVLbTLV, 2, true, true, 0.3);		// partons, jets, choose algorithm, use maxDist, use dR, set maxDist=0.3
                                
                                if (matching_non1stbJet.getNumberOfAvailableCombinations() != 1) cout << "matching.getNumberOfAvailableCombinations() = "<<matching_non1stbJet.getNumberOfAvailableCombinations()<<" .  This should be equal to 1 !!!"<<endl;
                                
                                vector< pair<unsigned int, unsigned int> > Non1stbJetPartonPair; // First one is jet number, second one is mcParticle number
                                Non1stbJetPartonPair.clear();
                                
                                for (unsigned int i = 0; i < mcParticlesTLV.size(); i++)
                                {
                                    int matchedJetNumber = matching_non1stbJet.getMatchForParton(i, 0);
                                    
                                    if (matchedJetNumber > -1) Non1stbJetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
                                }
                                if(debug)cout<< " the size of JetPartonPair vector is =  " << Non1stbJetPartonPair.size() << endl;
                                
                                for (unsigned int iseljet = 0 ; iseljet < Non1stbJetPartonPair.size(); iseljet++)
                                {
                                    unsigned int particlenb = Non1stbJetPartonPair[iseljet].first; // place in selectedjets vector
                                    unsigned int partonIDnb = Non1stbJetPartonPair[iseljet].second; // place in mcParticles_ vector
                                    
                                    
                                    if (fabs(mcParticlesMatching_[partonIDnb]->type()) < topId)
                                    {
                                        if ((mcParticlesMatching_[partonIDnb]->grannyType() == 25 || fabs(mcParticlesMatching_[partonIDnb]->motherType()) == 24) && fabs(mcParticlesMatching_[partonIDnb]->grannyType()) != topId)
                                        {
                                            // cout << "Find particles comes from Higgs decay" << endl;
                                            if (FCNCWJetsCounter_nonb <1)
                                            {
                                                FCNCWJet1_nonb_ = Non1stbJetPartonPair[particlenb];
                                                FCNCWJetsCounter_nonb++;
                                                //     cout << "Find FCNC W first Jet from JetsExcludingHighestCSVLb " << endl;
                                                W_1stJet_Pair.push_back(pair<unsigned int,unsigned int> (Non1stbJetPartonPair[iseljet].first,Non1stbJetPartonPair[iseljet].second));
                                                NPair.push_back("FCNCWJet1");
                                            }
                                            if (FCNCWJetsCounter_nonb ==1)
                                            {
                                                FCNCWJet2_nonb_ = Non1stbJetPartonPair[particlenb];
                                                FCNCWJetsCounter_nonb++;
                                                //    cout << "Find FCNC W second Jet from JetsExcludingHighestCSVLb" << endl;
                                                W_2ndJet_Pair.push_back(pair<unsigned int,unsigned int> (Non1stbJetPartonPair[iseljet].first,Non1stbJetPartonPair[iseljet].second));
                                                NPair.push_back("FCNCWJet2");
                                            }
                                            if (FCNCWJetsCounter_nonb == 2) {
                                                Find_FCNC_W = true;
                                                //    cout << "Find_FCNC_W =  " << Find_FCNC_W <<endl;
                                            }
                                            if (FCNCWJetsCounter == 2) {
                                                Find_FCNC_HdecayJets = true;
                                                //    cout << "Find_FCNC_HdecayJets =  " << Find_FCNC_HdecayJets <<endl;
                                                
                                            }
                                            
                                        }
                                    }
                                }
                                temp_Match_WJet1.Clear();
                                temp_Match_WJet2.Clear();
                                if (Find_FCNC_W)
                                {
                                    for (unsigned int iWjet =0 ; iWjet < W_1stJet_Pair.size() ; iWjet++)
                                    {
                                        temp_Match_WJet1.SetPxPyPzE(JetsExcludingHighestCSVLb[W_1stJet_Pair[iWjet].first]->Px(), JetsExcludingHighestCSVLb[W_1stJet_Pair[iWjet].first]->Py(), JetsExcludingHighestCSVLb[W_1stJet_Pair[iWjet].first]->Pz(), JetsExcludingHighestCSVLb[W_1stJet_Pair[iWjet].first]->E());
                                    }
                                    for (unsigned int iWjet = 0  ; iWjet < W_2ndJet_Pair.size(); iWjet++)
                                    {
                                        temp_Match_WJet2.SetPxPyPzE(JetsExcludingHighestCSVLb[W_2ndJet_Pair[iWjet].first]->Px(), JetsExcludingHighestCSVLb[W_2ndJet_Pair[iWjet].first]->Py(), JetsExcludingHighestCSVLb[W_2ndJet_Pair[iWjet].first]->Pz(), JetsExcludingHighestCSVLb[W_2ndJet_Pair[iWjet].first]->E());
                                    }
                                    
                                    
                                }
                               // Match_FCNC_MC_Mass_WJetPair= (temp_Match_WJet1 + temp_Match_WJet2).M();
                                
                                //                            for (unsigned int iPart = 0 ; iPart<PPair.size(); iPart++)
                                //                            {
                                //                                if(NPair[iPart].find("FCNCWJet1")!=string::npos){ temp_Match_WJet1.SetPxPyPzE(JetsExcludingHighestCSVLb[PPair[iPart].first]->Px(), JetsExcludingHighestCSVLb[PPair[iPart].first]->Py(), JetsExcludingHighestCSVLb[PPair[iPart].first]->Pz(), JetsExcludingHighestCSVLb[PPair[iPart].first]->E());
                                //                                }
                                //                                if(NPair[iPart].find("FCNCWJet2")!=string::npos){ temp_Match_WJet2.SetPxPyPzE(JetsExcludingHighestCSVLb[PPair[iPart].first]->Px(), JetsExcludingHighestCSVLb[PPair[iPart].first]->Py(), JetsExcludingHighestCSVLb[PPair[iPart].first]->Pz(), JetsExcludingHighestCSVLb[PPair[iPart].first]->E());
                                //
                                //                                }
                                //                            }
                                //                            Match_FCNC_MC_Mass_WJetPair= (temp_Match_WJet1 + temp_Match_WJet2).M();
                                //                            if (Find_SMNegB && Find_SMNegLep) {
                                //                                Find_SM_Antitop = true;
                                //                            }else if (Find_SMPosB && Find_SMPosLep) {
                                //                                Find_SM_top = true;
                                //                            }
                                
                                
                            }
                            

                            float PairMass = 0.;
                            bool W_RecoMass = false;
                            nJetpair = 0;
                            ///// **** Reconstruction of FCNC WJets **** /////
                            
                                if (JetsExcludingHighestCSVLb.size() >= 2)
                            {
                               // cout <<" selectedJets.size() =  " << selectedJets.size() << endl;
                              //  cout << "JetsExcludingHighestCSVLb.size()  =  " << JetsExcludingHighestCSVLb.size() << endl;
                                //cout<< "the size of JetsExcludingHighestCSVLb  =   " << JetsExcludingHighestCSVLb.size() <<endl;
                                TLorentzVector temp_W_1st_Jet;
                                TLorentzVector temp_W_2nd_Jet;
                                float mjj = 0.;
                                float massDiff_from_Wmass = 0.;
                                float MinimassDiff_from_Wmass = 999;
                                //float WPairMass,PairMass;
                                InvMass_JJ.clear();
                                
                                for (unsigned int xjet =JetsExcludingHighestCSVLb.size()-1; xjet > 0 ; xjet--)
                                {
                                   // cout<< "xjet =  " << xjet <<endl;
                                  
                                    for (unsigned int kjet = 0 ; kjet < xjet; kjet++)
                                    {
                                        //cout<< "xjet =  " << xjet <<endl;
                                       // cout << " kjet =  " << kjet << endl;
                                        temp_W_1st_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[xjet]->Px(),JetsExcludingHighestCSVLb[xjet]->Py(),JetsExcludingHighestCSVLb[xjet]->Pz(),JetsExcludingHighestCSVLb[xjet]->Energy());
                                        temp_W_2nd_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[kjet]->Px(),JetsExcludingHighestCSVLb[kjet]->Py(),JetsExcludingHighestCSVLb[kjet]->Pz(),JetsExcludingHighestCSVLb[kjet]->Energy());
                                        mjj = (temp_W_1st_Jet+ temp_W_2nd_Jet).M();
                                        Mass_JetPair = mjj;
                                       // Mass_JetPair_arr[nJetpair]=mjj;
                                       // nJetpair+1;
                                        InvMass_JJ.push_back(mjj);
//
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
                            
                           
                            Mass_WJetPair = PairMass;
//                            for (unsigned int i = 0 ; i<InvMass_JJ.size(); i++)
//                            {
//                              Mass_JetPair_arr[nJetpair] = InvMass_JJ[i];
//                                nJetpair+1;
//                            }
                            
                           // cout << "mass of W jets pair (outside loop ) =  "<< PairMass << endl;
                            
                            if (diElectron && !diMuon)
                            {
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                InvMass_lb = (tempLepton_0+SM_bJet).M();
                                //cout << " InvMass_lb =  " << InvMass_lb;
                                DeltaR_2L = tempLepton_0.DeltaR(tempLepton_1);
                                DeltaPhi_2L = tempLepton_0.DeltaPhi(tempLepton_1);
                                DeltaR_Elec0b0_DiElec = tempLepton_0.DeltaR(SM_bJet);
                                DeltaR_Elec1b0_DiElec = tempLepton_1.DeltaR(SM_bJet);
                            }
                            if (diMuon && !diElectron)
                            {
                                
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                InvMass_lb = (tempLepton_0+SM_bJet).M();
                                //cout << " InvMass_lb =  " << InvMass_lb;
                                //histo1D["h_2L_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                DeltaR_2L = tempLepton_0.DeltaR(tempLepton_1);
                                DeltaPhi_2L = tempLepton_0.DeltaPhi(tempLepton_1);
                                DeltaR_Mu0b0_DiMu = tempLepton_0.DeltaR(SM_bJet);
                                DeltaR_Mu1b0_DiMu = tempLepton_1.DeltaR(SM_bJet);
                            }
                            
                            if (diEMu || diMuE && !diMuon && !diElectron)
                            {
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                InvMass_lb = (tempLepton_0+SM_bJet).M();
                                DeltaR_2L = tempLepton_0.DeltaR(tempLepton_1);
                                DeltaPhi_2L = tempLepton_0.DeltaPhi(tempLepton_1);
                                if (diEMu)
                                {
                                    DeltaR_Elec0b0_DiElMu = tempLepton_0.DeltaR(SM_bJet);
                                    DeltaR_Mu1b0_DiElMu = tempLepton_1.DeltaR(SM_bJet);
                                }
                                if (diMuE)
                                {
                                    DeltaR_Mu0b0_DiMuEl = tempLepton_0.DeltaR(SM_bJet);
                                    DeltaR_Elec1b0_DiMuEl = tempLepton_1.DeltaR(SM_bJet);
                                }
                            }

                            
                            for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                            {
                                
                                sum_jet_PT+= selectedJets[ijet]->Pt();
                            }
                            Ht_AllJets=sum_jet_PT;
                            St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                           // histo1D["h_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                           // histo1D["h_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                            histo2D["h_2L_3Jets1b_Ht_vs_metPt"]->Fill(Ht_AllJets,met_Pt,scaleFactor*lumiWeight);
                            
                            histo2D["h_2L_2lDeltaPhi_vs_metPt"]->Fill(DeltaPhi_2L,met_Pt,scaleFactor*lumiWeight);
                            Ht=sum_jet_PT;
                            St=Ht_AllJets+Sum_Leptons_Pt;
                            
                            
                            
                            //////////////////////////////////////////////////////////////////////////////
                            /////////// Define the SS and OS region /////////////////////////////////////
                            ////////////////*************************** ////////////////////////////////
                            
                            if (diElectron && qElec0 == qElec1) // in case of Same Sign dielectron
                            {
                                SSdiLepton = true;
                                if(inZpeak)invMass_2SSL_Zmass= massZLepPair;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                
                                
                            }else if (diElectron && qElec0 != qElec1) // in case of Opposite Sign dielectron
                            {
                                OSdiLepton = true;
                              if(inZpeak)invMass_2OSL_Zmass= massZLepPair ;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                
                            }else if (diMuon && !diElectron && qMu0 == qMu1) // in case of Same Sign diMuon
                            {
                              if(inZpeak)invMass_2SSL_Zmass= massZLepPair;
                                SSdiLepton = true;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                
                            }else if (diMuon && !diElectron && qMu0 != qMu1) // in case of Opposite Sign diMuon
                            {
                              if(inZpeak)invMass_2OSL_Zmass= massZLepPair ;
                                OSdiLepton = true;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                
                            }else if (diEMu || diMuE && !diMuon && !diElectron && qLepton0 == qLepton1)
                            {
                                SSdiLepton = true;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                            }else if (diEMu || diMuE && !diMuon && !diElectron && qLepton0 != qLepton1) // in case of Opposite Sign dilepton
                            {
                                OSdiLepton = true;
                                Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                            }

                            histo1D["h_cutFlow_ChargeMisId"]->Fill(0.,scaleFactor*lumiWeight);
                            if (OSdiLepton)
                            {
                                histo1D["h_cutFlow_LeptonCuts"]->Fill(1., scaleFactor*lumiWeight);
                                histo1D["h_cutFlow_ChargeMisId"]->Fill(1.,scaleFactor*lumiWeight);
                                for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                {
                                    sum_jet_PT+= selectedJets[ijet]->Pt();
                                }
                                Ht_AllJets=sum_jet_PT;
                                St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                                OSLNoZmassVetoTree->Fill();
                            }

                            if (SSdiLepton)
                            {
                                histo1D["h_cutFlow_LeptonCuts"]->Fill(2., scaleFactor*lumiWeight);
                                histo1D["h_cutFlow_ChargeMisId"]->Fill(2.,scaleFactor*lumiWeight);
                                for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                {
                                    sum_jet_PT+= selectedJets[ijet]->Pt();
                                }
                                Ht_AllJets=sum_jet_PT;
                                St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;

                                SSLNoZmassVetoTree->Fill();
                            }
                            
                            
                            NoZmassVetoTree ->Fill(); //  fill trees before Zmass window cut
                            
                            
                            if(debug)cout << "the mass of dilepton before Zmass veto =  " <<invMass_2L << endl;
                            if(debug)cout << "the mass of WJets before Zmass veto =  " << Mass_WJetPair << endl;
                        
                           // if ((diElectron && fabs(Zmass - invMass_2L ) >= 15) || (diMuon && fabs(Zmass - invMass_2L ) >= 15) || diEMu || diMuE)
                           // {
                                ///// histograms to estimate Charge mis Id in DY samples
                                //if ((diElectron || diEMu || diMuE) && isDY && ApplyCharge_misID)
                                if (diElectron  && isDY && ApplyCharge_misID)
                                {
                                    ElectronChargeMisId(selectedElectrons);
                                }
                                histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                                histo1D["h_cutFlow_LeptonCuts"]->Fill(3., scaleFactor*lumiWeight);
                                if(debug)cout << "the mass of dilepton After Zmass veto =  " <<InvMass_ll << endl;
                                if(debug)cout << "And the mass difference (Zmass - InvMass_ll)  =  " <<InvMass_ll << endl;
                                if(debug)cout << "the mass of WJets After Zmass veto =  " << Mass_WJetPair << endl;
                                
                                ///// filling trees
                                myTree->Fill();
                                nbEvents_9++;
                            
                                histo1D["h_cutFlow_ChargeMisId"]->Fill(3.,scaleFactor*lumiWeight);
                                if (OSdiLepton && !SSdiLepton)
                                {
                                 //   nCuts++;
                                    nbEvents_10++;
                                    
                                   // Count_cut[9]=Count_cut[9]+cutstep;
                                    
                                    histo1D["h_cutFlow"]->Fill(8., scaleFactor*lumiWeight);
                                    histo1D["h_cutFlow_LeptonCuts"]->Fill(4., scaleFactor*lumiWeight);
                                    histo1D["h_cutFlow_ChargeMisId"]->Fill(4.,scaleFactor*lumiWeight);
                                    histo2D["h_2OSL_3Jets1b_Ht_vs_metPt"]->Fill(Ht_AllJets,met_Pt,scaleFactor*lumiWeight);
                                    histo2D["h_2OSL_2lDeltaPhi_vs_metPt"]->Fill(DeltaPhi_2L,met_Pt,scaleFactor*lumiWeight);
                                    
                                    OSLeptonTree->Fill();
                                }
                                if (SSdiLepton && !OSdiLepton)
                                {
                                 //   nCuts++;
                                    nbEvents_11++;
                                   // Count_cut[10]=Count_cut[10]+cutstep;
                                    histo1D["h_cutFlow_ChargeMisId"]->Fill(5.,scaleFactor*lumiWeight);
                                    histo1D["h_cutFlow"]->Fill(9., scaleFactor*lumiWeight);
                                    histo1D["h_cutFlow_LeptonCuts"]->Fill(5., scaleFactor*lumiWeight);
                                    histo2D["h_2SSL_3Jets1b_Ht_vs_metPt"]->Fill(Ht_AllJets,met_Pt,scaleFactor*lumiWeight);
                                    histo2D["h_2SSL_2lDeltaPhi_vs_metPt"]->Fill(DeltaPhi_2L,met_Pt,scaleFactor*lumiWeight);
                                    SSLeptonTree->Fill();
                                }
                                
                                eventSelected = true;
                                nofSelectedEvents++;
                                
                           // }
                            
                            //////*** Apply DD for charge misId /////
                            if (isData && diElectron && OSdiLepton && ApplyCharge_misID)
                            {
//                                bool twoElecInBarrel = false;
//                                bool twoElecInEndcap = false;
//                                bool oneElecInBarreloneInEndCap = false;
//                                
//                                Lepton0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
//                                if(selectedElectrons.size()>1)Lepton1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
//                                
//                                if (abs(Lepton0.Eta()) < 1.47 && abs(Lepton1.Eta()) < 1.47)
//                                {
//                                    twoElecInBarrel = true;
//                                }
//                                if (abs(Lepton0.Eta()) > 1.47 && abs(Lepton0.Eta()) < 2.5 && Lepton1.Eta() > 1.47 && Lepton1.Eta() < 2.5)
//                                {
//                                    twoElecInEndcap = true;
//                                }
//                                if ((abs(Lepton0.Eta()) < 1.47 && Lepton1.Eta() > 1.47 && Lepton1.Eta() < 2.5) || (abs(Lepton0.Eta()) > 1.47 && abs(Lepton0.Eta()) < 2.5 && Lepton1.Eta() < 1.47) {
//                                    oneElecInBarreloneInEndCap=true;
//                                }
                                
                                vector<int> electroncharge;
                                electroncharge.clear();
                                //// my calculated ratios
                                float Pro_B_Pt25= 0.00025 , Pro_B_Pt25_50 = 0.00018 , Pro_B_Pt_50 = 0.00012;
                                float Pro_E_Pt25= 0.0004 , Pro_E_Pt25_50 = 0.0013 , Pro_E_Pt_50 = 0.0015;
                                
                                
                                //    ////// values for checking only
//                                    float Pro_B_Pt25= 0.5 , Pro_B_Pt25_50 = 0.5 , Pro_B_Pt_50 = 0.5;
//                                    float Pro_E_Pt25= 0.5 , Pro_E_Pt25_50 = 0.5 , Pro_E_Pt_50 = 0.5;

                                /////
                                //int nbofFlips = 0;
                                //bool isFlips = false;
                                /////
                                if (abs(selectedElectrons[0]->Eta()) < 1.47 && abs(selectedElectrons[1]->Eta()) < 1.47)
                                {
                                    int nbofFlips = 0;
                                    bool isFlips = false;
                                    for (unsigned int iElec = 0; iElec < selectedElectrons.size() ; iElec++)
                                    {
                                        electroncharge.push_back(selectedElectrons[iElec]->charge());
                                        ///now applying charge-misid to electrons
                                        //creates a randm number between 0-1.
                                        double rdmnr = gRandom->Uniform();
                                        cout<< "the value of random variable =  " << rdmnr << endl;
                                        
                                        if (selectedElectrons[iElec]->Pt() <25 && rdmnr <= Pro_B_Pt25)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && rdmnr <= Pro_B_Pt25_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >= 50 && rdmnr <= Pro_B_Pt_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
                                        }

                                    }
                                    cout << " 2 electrons are in the barrel : isFlip is =   "<< isFlips << " and the nb of electron charge flips =  " << nbofFlips <<endl;
                                    if(isFlips && nbofFlips ==1 )histo1D["h_ElecChargeFlips_2atBarrel_CutFlow"]->Fill(1);
                                    if(isFlips && nbofFlips == 2 )histo1D["h_ElecChargeFlips_2atBarrel_CutFlow"]->Fill(2);
                                    
                                }else if ((abs(selectedElectrons[0]->Eta()) < 1.47 && abs(selectedElectrons[1]->Eta()) > 1.47 && abs(selectedElectrons[1]->Eta()) <2.5) || (abs(selectedElectrons[1]->Eta()) < 1.47 && abs(selectedElectrons[0]->Eta()) > 1.47 && abs(selectedElectrons[0]->Eta()) <2.5) )
                                {
                                    bool isFlips = false;
                                    int nbofFlips = 0;
                                    for (unsigned int iElec = 0; iElec < selectedElectrons.size() ; iElec++)
                                    {
                                        electroncharge.push_back(selectedElectrons[iElec]->charge());
                                        ///now applying charge-misid to electrons
                                        //creates a randm number between 0-1.
                                        double rdmnr = gRandom->Uniform();
                                        cout<< "the value of random variable =  " << rdmnr << endl;
                                        
                                        if (selectedElectrons[iElec]->Pt() <25 && abs(selectedElectrons[iElec]->Eta()) < 1.47 && rdmnr <= Pro_B_Pt25)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && abs(selectedElectrons[iElec]->Eta()) < 1.47 && rdmnr <= Pro_B_Pt25_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >= 50 && abs(selectedElectrons[iElec]->Eta()) < 1.47 && rdmnr <= Pro_B_Pt_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
                                        }else if (selectedElectrons[iElec]->Pt() <25 && abs(selectedElectrons[iElec]->Eta()) > 1.47 && abs(selectedElectrons[iElec]->Eta()) < 2.5 && rdmnr <= Pro_E_Pt25)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && abs(selectedElectrons[iElec]->Eta()) > 1.47 && abs(selectedElectrons[iElec]->Eta()) < 2.5 && rdmnr <= Pro_E_Pt25_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >= 50 && abs(selectedElectrons[iElec]->Eta()) > 1.47 && abs(selectedElectrons[iElec]->Eta()) < 2.5 && rdmnr <= Pro_E_Pt_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
                                        }
                                        
                                    }
                                    cout << " 1 electrons is in the Endcap and 1 at Barrel : isFlip is =   "<< isFlips << " and the nb of electron charge flips =  " << nbofFlips <<endl;
                                    if(isFlips && nbofFlips == 1 )histo1D["h_ElecChargeFlips_1atBarrelEndCap_CutFlow"]->Fill(1);
                                    if(isFlips && nbofFlips == 2 )histo1D["h_ElecChargeFlips_1atBarrelEndCap_CutFlow"]->Fill(2);
                                }else if ((abs(selectedElectrons[0]->Eta()) > 1.47 && abs(selectedElectrons[0]->Eta()) <2.5) || (abs(selectedElectrons[1]->Eta()) > 1.47 && abs(selectedElectrons[1]->Eta()) <2.5) )
                                {
                                    bool isFlips = false;
                                    int nbofFlips = 0;
                                    for (unsigned int iElec = 0; iElec < selectedElectrons.size() ; iElec++)
                                    {
                                        electroncharge.push_back(selectedElectrons[iElec]->charge());
                                        ///now applying charge-misid to electrons
                                        //creates a randm number between 0-1.
                                        double rdmnr = gRandom->Uniform();
                                        cout<< "the value of random variable =  " << rdmnr << endl;
                                        
                                        if (selectedElectrons[iElec]->Pt() <25 && rdmnr < Pro_E_Pt25)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && rdmnr <= Pro_E_Pt25_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
                                        }else if(selectedElectrons[iElec]->Pt() >= 50 && rdmnr <= Pro_E_Pt_50)
                                        {
                                            electroncharge[iElec] = electroncharge[iElec] * -1;
                                            isFlips = true;
                                            nbofFlips ++;
                                            cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
                                        }
                                        
                                    }
                                    cout << " 2 electrons are in the Endcap : isFlip is =   "<< isFlips << " and the nb of electron charge flips =  " << nbofFlips <<endl;
                                    if(isFlips && nbofFlips == 1 )histo1D["h_ElecChargeFlips_2atEndCap_CutFlow"]->Fill(1);
                                    if(isFlips && nbofFlips == 2 )histo1D["h_ElecChargeFlips_2atEndCap_CutFlow"]->Fill(2);
                                }
                            }

                            
                                
                                //////////
//                                for (unsigned int iElec = 0; iElec < selectedElectrons.size() ; iElec++)
//                                {
//                                    electroncharge.push_back(selectedElectrons[iElec]->charge());
//                                    ///now applying charge-misid to electrons
//                                    //creates a randm number between 0-1.
//                                    double rdmnr = gRandom->Uniform();
//                                    cout<< "the value of random variable =  " << rdmnr << endl;
//                                    if (abs(selectedElectrons[iElec]->Eta()) < 1.47)
//                                    {
//                                        if (selectedElectrons[iElec]->Pt() <25 && rdmnr < Pro_B_Pt25)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
//                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && rdmnr < Pro_B_Pt25_50)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
//                                        }else if(selectedElectrons[iElec]->Pt() >= 50 && rdmnr < Pro_B_Pt_50)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
//                                        }
//                                    }
//                                    if (abs(selectedElectrons[iElec]->Eta()) > 1.47 && abs(selectedElectrons[iElec]->Eta()) <2.5)
//                                    {
//                                        if (selectedElectrons[iElec]->Pt() <25 && rdmnr < Pro_E_Pt25)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at EndCap and with Pt <25 is flipped!"<<endl;
//                                        }else if(selectedElectrons[iElec]->Pt() >=25 && selectedElectrons[iElec]->Pt() <50 && rdmnr < Pro_E_Pt25_50)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at EndCap and with Pt >25 <50 is flipped!"<<endl;
//                                        }else if (selectedElectrons[iElec]->Pt() >= 50 && rdmnr < Pro_E_Pt_50)
//                                        {
//                                            electroncharge[iElec] = electroncharge[iElec] * -1;
//                                            cout<<"charge of the electron at EndCap and with Pt >50 is flipped!"<<endl;
//                                        }
//                                    }
//                                    
//                                }
//                            }
                            
                            ///////////////////////////////////////////
                            

                            
                            
                            
                        }//2L+>=2Jets+>=1CSVLB
                        
                        
                    }//2L+>=3Jets
                    
                }//2L
                
            }//2L
            
            
            
            
            
            if (! eventSelected )
            {
                continue;
            }
            
            if(verbose>3)cout << "filling the tree" << endl;
            //myTree->Fill();
           // nofSelectedEvents++;
           // infoFile << "|" << evt_num << "|"  << "|"  <<Channel << "|"  << endl;
            
            //if (verbose > 0) cout << "  Event " << ievt << " is selected" << endl;
            // histo1D["h_Nb_Events_Lumi"]->Fill(2);
        } // end the loop over events

        
       
        sumW = (int) sumWeights;
        nEv = (int) nEvents;
        nEv_lumi = (int)nbEvents_eEqLumi;
        Count_cut[0] = nbEvents_0;
        Count_cut[1] = nbEvents_1PU;
        Count_cut[2] = nbEvents_2GPV;
        Count_cut[3] = nbEvents_3Trig;
        Count_cut[4] = nbEvents_4BTag;
        Count_cut[5] = nbEvents_5LepSF;
        Count_cut[6] = nbEvents_6;
        Count_cut[7] = nbEvents_7;
        Count_cut[8] = nbEvents_8;
        Count_cut[9] = nbEvents_9;
        Count_cut[10] = nbEvents_10;
        Count_cut[11] = nbEvents_11;
        cout << "*****************************************************************************************************" <<endl;
        cout << "**The nb of events before any selections or applying SF (nbEvents_0) =   " << nbEvents_0 << endl;
        cout << "**The nb of events before any selections & after Applying PU SF (nbEvents_1PU) =  " << nbEvents_1PU << endl;
        cout << "**The nb of events before any selections & after Good Primary Vertex and pass met filters (nbEvents_2GPV) =  " << nbEvents_2GPV << endl;
        cout << "**The nb of events before any selections & after triggers (nbEvents_3Trig) =  " << nbEvents_3Trig << endl;
        cout << "**The nb of events before any selections & after Applying bTag SF (nbEvents_4BTag) =  " << nbEvents_4BTag << endl;
        cout << "**The nb of events after 2 leptons & after Applying Lep SF (nbEvents_5LepSF) =  " << nbEvents_5LepSF << endl;
        cout << "**The nb of events after 2 leptons (ee, uu , eu) (nbEvents_6) =  " << nbEvents_6 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets  (nbEvents_7) =  " << nbEvents_7 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged  (nbEvents_8) =  " << nbEvents_8 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto (nbEvents_9) =  " << nbEvents_9 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto 2OSL Channel (nbEvents_10) =  " << nbEvents_10 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto 2SSL Channel (nbEvents_11) =  " << nbEvents_11 << endl;
        cout << "*****************************************************************************************************" <<endl;
//        cout << "the number of jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_Jets <<endl;
//        cout << "the number of CSVLbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVLbJets <<endl;
//        cout << "the number of CSVMbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVMbJets <<endl;
//        cout << "the number of CSVTbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVTbJets <<endl;
        
        globalTree->Fill();
        
        
        infoFile << "**The nb of events before any selections or applying SF (nbEvents_0) =   " << nbEvents_0 << endl;
        infoFile << "**The nb of events before any selections & after Applying PU SF (nbEvents_1PU) =  " << nbEvents_1PU << endl;
        infoFile << "**The nb of events before any selections & after Good Primary Vertex (nbEvents_2GPV) =  " << nbEvents_2GPV << endl;
        infoFile << "**The nb of events before any selections & after triggers (nbEvents_3Trig) =  " << nbEvents_3Trig << endl;
        infoFile << "**The nb of events before any selections & after Applying bTag SF (nbEvents_4BTag) =  " << nbEvents_4BTag << endl;
        infoFile << "**The nb of events before any selections & after Applying Lep SF (nbEvents_5LepSF) =  " << nbEvents_5LepSF << endl;
        infoFile << "**The nb of events after 2 leptons (nbEvents_6) =  " << nbEvents_6 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets  (nbEvents_7) =  " << nbEvents_7 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged  (nbEvents_8) =  " << nbEvents_8 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto (nbEvents_9) =  " << nbEvents_9 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged +Zmass Veto 2OSL Channel (nbEvents_10) =  " << nbEvents_10 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged +Zmass Veto 2SSL Channel (nbEvents_11) =  " << nbEvents_11 << endl;
        
        //////////////////////
        ///  END OF EVENT  ///
        //////////////////////
        
        cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
        cout << "the nb of events before charge cut =  " << nbofEventsbeforeChargeCut << "  &  the nb of events After charge cut =  " << nbofEventsAfterChargeCut << "  So the nb of removed events is =   " << nbofEventsbeforeChargeCut - nbofEventsAfterChargeCut <<endl;
        if (! isData )
        {
            cout << "Data set " << datasets[d]->Title() << " has " << nofPosWeights << " events with positive weights and " << nofNegWeights << " events with negative weights." << endl;
            cout << "         Pos - neg is " << nofPosWeights - nofNegWeights << ", pos + neg is " << nofPosWeights + nofNegWeights << endl;
            cout << "The sum of the weights is " << ((int)sumWeights) << ", whereas the total number of events is " << ((int)nEvents) << endl;
            
            // Determine scale factor due to negative weights
            nloSF = ((double) (nofPosWeights - nofNegWeights))/((double) (nofPosWeights + nofNegWeights));
            cout << "This corresponds to an event scale factor of " << nloSF  << endl;
        }
        infoFile.close();
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
        
        if(!isData && !btagShape) delete btwt_CSVv2L_mujets_central;
        
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
       // cout << "2D Plot: " << it->first << endl;
        TCanvas *ctemp = new TCanvas();
        ctemp->cd();
        TH2F *temp = it->second;
        string name = it->first;
        temp->Draw();
        delete ctemp;
    }
    
    
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


/////// ****** Charge misidetification ******////////
//// charge-misid to electrons: 0.3% in endcap region, 0.03% in barrel region (Gerrit suggestion)
////now applying charge-misid to electrons: 2.9% in endcap region, 0.55% in barrel region // An 2014
//////selective  charge-misid in barrel (MC 0.017 +- 0.002%) (Data 0.020 +- 0.002%) // EGM 13-001
//////selective  charge-misid in endcap (MC 0.21 +- 0.02%) (Data 0.23 + - 0.02%) // EGM 13-001

void ElectronChargeMisId(vector <TRootElectron*> SelElectron)
{
    vector<int> electroncharge;
    electroncharge.clear();
    //// my calculated ratios
    float Pro_B_Pt25= 0.00025 , Pro_B_Pt25_50 = 0.00018 , Pro_B_Pt_50 = 0.00012;
    float Pro_E_Pt25= 0.0004 , Pro_E_Pt25_50 = 0.0013 , Pro_E_Pt_50 = 0.0015;
    
    
//    ////// values for checking only 
//    float Pro_B_Pt25= 1. , Pro_B_Pt25_50 = 1. , Pro_B_Pt_50 = 1.;
//    float Pro_E_Pt25= 1. , Pro_E_Pt25_50 = 1. , Pro_E_Pt_50 = 1.;

    
    for (unsigned int iElec = 0 ; iElec < SelElectron.size() ; iElec++)
    {
        electroncharge.push_back(SelElectron[iElec]->charge());
        ///now applying charge-misid to electrons
        //creates a randm number between 0-1.
        double rdmnr = gRandom->Uniform();
         cout<< "the value of random variable =  " << rdmnr << endl;
        if (abs(SelElectron[iElec]->Eta()) < 1.47)
        {
            if (SelElectron[iElec]->Pt() <25 && rdmnr < Pro_B_Pt25)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at barrel and with Pt <25 is flipped!"<<endl;
            }else if(SelElectron[iElec]->Pt() >=25 && SelElectron[iElec]->Pt() <50 && rdmnr < Pro_B_Pt25_50)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at barrel and with Pt >25 <50 is flipped!"<<endl;
            }else if(SelElectron[iElec]->Pt() >= 50 && rdmnr < Pro_B_Pt_50)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at barrel and with Pt >50 is flipped!"<<endl;
            }
        }
        if (abs(SelElectron[iElec]->Eta()) > 1.47 && abs(SelElectron[iElec]->Eta()) <2.5)
        {
            if (SelElectron[iElec]->Pt() <25 && rdmnr < Pro_E_Pt25)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at EndCap and with Pt <25 is flipped!"<<endl;
            }else if(SelElectron[iElec]->Pt() >=25 && SelElectron[iElec]->Pt() <50 && rdmnr < Pro_E_Pt25_50)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at EndCap and with Pt >25 <50 is flipped!"<<endl;
            }else if (SelElectron[iElec]->Pt() >= 50 && rdmnr < Pro_E_Pt_50)
            {
                electroncharge[iElec] = electroncharge[iElec] * -1;
                cout<<"charge of the electron at EndCap and with Pt >50 is flipped!"<<endl;
            }
        }
        
        
    }
}

