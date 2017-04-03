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
///int ElectronChargeId(TLorentzVector SelElectron);

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
    
    
    //// *** Working Conditions ////
    bool Elec_Elec, Mu_Mu, Elec_Mu, Mu_Elec, Apply_HLT_Triggers, eventSelected, Fake_Electrons, ApplyCharge_misID, ApplyElec_SF , ApplyMu_SF , ApplyPU_SF, Apply_btag_SF, Apply_JetCleaning, trigged,debug, printTrigger, All_lep;
    Elec_Elec = true;
    ApplyElec_SF = true;
    Mu_Mu =false;
    ApplyMu_SF = false;
    Elec_Mu = false;
    Mu_Elec = false;
    All_lep = false;
    Apply_HLT_Triggers = true;
    printTrigger = false;
    eventSelected= false;
    Fake_Electrons = false;
    ApplyCharge_misID = false;
    ApplyPU_SF = true;
    Apply_JetCleaning = true;
    trigged = false;
    debug = false;
    bool applyJER = true;
    bool applyJES = true;
    bool applyNegWeightCorrection = true;
    bool btagShape = true;
    bool fillBtagHisto = false;
    
    std::string channelpostfix = "";
    string runDate = "Test_correctingChMisId_31Mar17";
    
    /////////////////////
    ///  Configuration
    /////////////////////
    
    
    /// xml file
    string xmlFileName ="";
    string Channel = "";
    if (argc > 1) xmlFileName = (string)argv[1];
    const char *xmlfile = xmlFileName.c_str();
    double dataLumi = 0;
    
    //Setting Lepton Channels
    if(Elec_Elec)
    {
        cout << " --> Using the Electron-Electron channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_ElEl_V4_Samples.xml";
        dataLumi = 36768.376571746; //Runs from B to H
        channelpostfix = "_ElEl_";
        Channel = "Dilepton_ElecElec";
    }
    else if(Mu_Mu)
    {
        cout << " --> Using the Muon-Muon channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_MuMu_V4_Samples.xml";
        dataLumi = 36768.376651703; //Run B+C+D+E+F+G+H
        channelpostfix = "_MuMu_";
        Channel = "Dilepton_MuMu";
    }
    else if(Elec_Mu)
    {
        cout << " --> Using the Electron-Muon channel..." << endl;
        xmlFileName ="config/Run2SameSignDiLepton_80X_ElMu_V4_Samples.xml";
        dataLumi = 2298.292131932;
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
    
    string info_dir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_24/src/TopBrussels/FCNCAnalysis/Information/"+Channel +"/";
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
    
    // Print information to a textfile
    
    
    ////////////////////////////////////
    ///  AnalysisEnvironment
    ////////////////////////////////////
    
    //TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    
    AnalysisEnvironment anaEnv;
    cout << " - Loading environment ..." << endl;
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    ////anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.ElectronCollection = "Electrons_calibratedPatElectrons"; ///     used from tag 80X_v2 onwards
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
    
   LumiReWeighting LumiWeights("../TopTreeAnalysisBase/Calibrations/PileUpReweighting/MCPileup_Summer16.root", "../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2016Data80X_Run271036-284044Cert__Full2016DataSet.root", "pileup", "pileup");
    
    
    ///// lepton scaling factors
    
    //MuonSFWeight (const string &sfFile, const string &dataOverMC, const bool &extendRange, const bool &debug, const bool &printWarning)
    
   MuonSFWeight *muonSFWeightID;
   MuonSFWeight *muonSFWeightIso;
    
    
    if (ApplyMu_SF && Mu_Mu)
    {
        if (dName.find("DataRun2016B")!=string::npos || dName.find("DataRun2016C")!=string::npos || dName.find("DataRun2016D")!=string::npos ||dName.find("DataRun2016E")!=string::npos || dName.find("DataRun2016F")!=string::npos )
        {
            muonSFWeightID = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
            muonSFWeightIso = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonIso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
            
        }
        else if(dName.find("DataRun2016G")!=string::npos || dName.find("DataRun2016H_v2")!=string::npos || dName.find("DataRun2016H_v3")!=string::npos)
        {
            
            muonSFWeightID = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonID_EfficienciesAndSF_GH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
            muonSFWeightIso = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonIso_EfficienciesAndSF_GH.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
            
        }
        else
        {
            muonSFWeightID = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
            //
            muonSFWeightIso = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonIso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
//                MuonSFWeight *muonSFWeightID = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonID_EfficienciesAndSF_BCDEF.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio", true, false, false); // Tight ID
//            //
//                MuonSFWeight *muonSFWeightIso = new MuonSFWeight(pathToCaliDir+"LeptonSF/MuonSF/"+"MuonIso_EfficienciesAndSF_BCDEF.root", "TightISO_TightID_pt_eta/abseta_pt_ratio", true, false, false);  // Tight RelIso
            
        }
    }
    
    

    ////Triggers SF for muons to be added
    

    string electronFile= "egammaEffi.txt_EGM2D_CutBasedTightID.root";
    string electronRecoFile = "egammaEffi.txt_EGM2D_RecoEff.root";
    string elecHistName = "EGamma_SF2D";
    
    ElectronSFWeight* electronSFWeight = new ElectronSFWeight (pathToCaliDir+"LeptonSF/ElectronSF/Moriond17/"+electronFile,elecHistName, true,false, false); // (... , ... , debug, print warning)
    ElectronSFWeight* electronSFWeightReco = new ElectronSFWeight(pathToCaliDir+"LeptonSF/ElectronSF/Moriond17/"+electronRecoFile,elecHistName, true,false, false);

    
    ///// b-tagging scaling factor
    
    BTagCalibration * btagcalib;
    BTagCalibrationReader * btagreader;
    BTagCalibrationReader *reader_csvv2;
    BTagWeightTools *btwt = 0;

    
    
    
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

    string histoDir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_24/src/TopBrussels/FCNCAnalysis/Output_Histos/";
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
    
    titlePlot = "NSSL_barrel_cutFlow"+channelpostfix;
    histo1D["h_NSSL_barrel_cutFlow"] = new TH1F(titlePlot.c_str(), "NSSL_barrel_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NOSL_barrel_cutFlow"+channelpostfix;
    histo1D["h_NOSL_barrel_cutFlow"] = new TH1F(titlePlot.c_str(), "NOSL_barrel_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NSSL_EndCap_cutFlow"+channelpostfix;
    histo1D["h_NSSL_EndCap_cutFlow"] = new TH1F(titlePlot.c_str(), "NSSL_EndCap_cutflow", 6,-0.5,5.5);
    
    titlePlot = "NOSL_EndCap_cutFlow"+channelpostfix;
    histo1D["h_NOSL_EndCap_cutFlow"] = new TH1F(titlePlot.c_str(), "NOSL_EndCap_cutflow", 6,-0.5,5.5);
    
    
    //*** histos for Jets *** //
    titlePlot = "initial_Nb_Jets"+channelpostfix;
    histo1D["h_initial_Nb_Jets"] = new TH1F(titlePlot.c_str(), "Initial nb. of jets",  16, - 0.5, 15.5 );

    
    //////***** CSVLBJet *****////////
//    titlePlot = "initial_Nb_CSVLBJets"+channelpostfix;
//    histo1D["h_initial_Nb_CSVLBJets"] = new TH1F(titlePlot.c_str(), "Initial nb. of CSVLBJets",  16, - 0.5, 15.5 );
//    
//    titlePlot = "2L_2J_1bJ_1st_CSVLBJet_Pt"+channelpostfix;
//    histo1D["h_2L_2J_1bJ_1st_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + cut: 1st CSVLBJet P_{T}",  100, 0, 500 );
//    
//    titlePlot = "2L_2J_1bJ_2nd_CSVLBJet_Pt"+channelpostfix;
//    histo1D["h_2L_2J_1bJ_2nd_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 2nd CSVLBJet P_{T}",  100, 0, 500 );
//    
//    titlePlot = "2L_2J_1bJ_3rd_CSVLBJet_Pt"+channelpostfix;
//    histo1D["h_2L_2J_1bJ_3rd_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 3rd CSVLBJet P_{T}",  100, 0, 500 );
//    
//    titlePlot = "2L_2J_1bJ_4th_CSVLBJet_Pt"+channelpostfix;
//    histo1D["h_2L_2J_1bJ_4th_CSVLBJet_Pt"] = new TH1F(titlePlot.c_str(), "After 2L + >= 2J + >=1bCSVLBJet cut: 4th CSVLBJet P_{T}",  100, 0, 500 );
//    
//    
//    
//    //*** histos for Muons *** //
//    if (Mu_Mu)
//    {
//
//        titlePlot = "2SSL_mMu0b0"+channelpostfix;
//        histo1D["h_2SSL_mMu0b0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: Mass b_{0} Mu_{0}",  100, 0., 500);
//        
//        titlePlot = "2OSL_DeltaR_b0_Mu0"+channelpostfix;
//        histo1D["h_2OSL_DeltaR_b0_Mu0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R b_{0} Mu_{0}",  100, 0, 5);
//        
//        titlePlot = "2OSL_mMu0b0"+channelpostfix;
//        histo1D["h_2OSL_mMu0b0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: Mass b_{0} Mu_{0}",  100, 0., 500);
//        
//    }
//    
//    if (Elec_Elec)
//    {
//        
//        titlePlot = "2SSL_mElec0b0"+channelpostfix;
//        histo1D["h_2SSL_mElec0b0"] = new TH1F(titlePlot.c_str(), "After 2SSL  cuts: Mass b_{0} Elec_{0}",  100, 0., 500);
//        
//        titlePlot = "2OSL_DeltaR_b0_Elec0"+channelpostfix;
//        histo1D["h_2OSL_DeltaR_b0_Elec0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: #Delta R b_{0} Elec_{0}",  100, 0, 5);
//        
//        titlePlot = "2OSL_mElec0b0"+channelpostfix;
//        histo1D["h_2OSL_mElec0b0"] = new TH1F(titlePlot.c_str(), "After 2OSL  cuts: Mass b_{0} Elec_{0}",  100, 0., 500);
//        
//    }
//
//    
//    ///*** histos for mets ***///
//    
//    //        titlePlot = "initial_met"+channelpostfix;
//    //        histo1D["h_initial_met"] = new TH1F(titlePlot.c_str(), " initial missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
//    
//    //        titlePlot = "2L_2J_1bJ_met"+channelpostfix;
//    //        histo1D["h_2L_2J_1bJ_met"] = new TH1F(titlePlot.c_str(), "After 2L + >=2J + 1bJet missing E_{T}; E_{T}^{mis} [GeV]", 100, 0,500);
//    
//    titlePlot = "Ht_AllJets"+channelpostfix;
//    histo1D["h_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts:  H_{T}",  150, 0, 1500 );
//    
//    titlePlot = "St_AllJets+Lep"+channelpostfix;
//    histo1D["h_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts:  S_{T}",  150, 0, 1500 );
//    
//    titlePlot = "2SSL_Ht_AllJets"+channelpostfix;
//    histo1D["h_2SSL_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts 2SSL:  H_{T}",  150, 0, 1500 );
//    
//    titlePlot = "2SSL_St_AllJets+Lep"+channelpostfix;
//    histo1D["h_2SSL_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts 2SSL:  S_{T}",  150, 0, 1500 );
//    
//    titlePlot = "2OSL_Ht_AllJets"+channelpostfix;
//    histo1D["h_2OSL_Ht_AllJets"] = new TH1F(titlePlot.c_str(), "After All Cuts 2OSL:  H_{T}",  150, 0, 1500 );
//    
//    titlePlot = "2OSL_St_AllJets+Lep"+channelpostfix;
//    histo1D["h_2OSL_St_AllJets+Lep"] = new TH1F(titlePlot.c_str(), "After All Cuts 2OSL:  S_{T}",  150, 0, 1500 );
//    
    
    
    
    titlePlot = "Charge_misId_Pt_barrel_Vec"+channelpostfix;
    histo1D["h_Charge_misId_Pt_barrel_Vec"] = new TH1F(titlePlot.c_str()," Charge_misId_Ratio_Barrel ", 21, -0.05, .25);
//    titlePlot = "Charge_misId_Pt_barrel2_Vec"+channelpostfix;
//    histo1D["h_Charge_misId_Pt_barrel2_Vec"] = new TH1F(titlePlot.c_str()," Charge_misId_Ratio_Barrel ", 21, -0.05, .25);
    
    titlePlot = "Charge_misId_Pt_EndCap_Vec"+channelpostfix;
    histo1D["h_Charge_misId_Pt_EndCap_Vec"] = new TH1F(titlePlot.c_str()," Charge_misId_Ratio_EndCap ", 11, -0.005, .105);
    

    
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
    
    
    titlePlot = "N2SSL_barrel_Pt"+channelpostfix;
    histo2D["h_N2SSL_barrel_Pt"] = new TH2F(titlePlot.c_str()," P_{T} vs N2SSL barrel  ",6,-0.5,5.5 , 100, 0, 10);
    
    titlePlot = "N2OSL_barrel_Pt"+channelpostfix;
    histo2D["h_N2OSL_barrel_Pt"] = new TH2F(titlePlot.c_str()," P_{T} vs N2OSL barrel  ",6,-0.5,5.5 , 1000, 0, 10);
    
    titlePlot = "N2SSL_EndCap_Pt"+channelpostfix;
    histo2D["h_N2SSL_EndCap_Pt"] = new TH2F(titlePlot.c_str()," P_{T} vs N2SSL EndCap  ",6,-0.5,5.5 , 1000, 0, 10);
    
    titlePlot = "N2OSL_EndCap_Pt"+channelpostfix;
    histo2D["h_N2OSL_EndCap_Pt"] = new TH2F(titlePlot.c_str()," P_{T} vs N2OSL EndCap  ",6,-0.5,5.5 , 1000, 0, 10);
    
    titlePlot = "NSSL_Pt_eta"+channelpostfix;
    histo2D["h_NSSL_Pt_eta"] = new TH2F(titlePlot.c_str()," NSSL P_{T} vs #eta  ",6,-0.5,5.5 , 6,-0.5,5.5);
    
    titlePlot = "NOSL_Pt_eta"+channelpostfix;
    histo2D["h_NOSL_Pt_eta"] = new TH2F(titlePlot.c_str()," NOSL P_{T} vs #eta  ",6,-0.5,5.5 , 6,-0.5,5.5);


    
    titlePlot = "Charge_misId_Pt_barrel"+channelpostfix;
    histo2D["h_Charge_misId_Pt_barrel"] = new TH2F(titlePlot.c_str()," P_{T} vs Charge_misId_Ratio_Barrel ",6,-0.5,5.5 , 10, 0, 0.02);
    
    titlePlot = "Charge_misId_Pt_EndCap"+channelpostfix;
    histo2D["h_Charge_misId_Pt_EndCap"] = new TH2F(titlePlot.c_str()," P_{T} vs Charge_misId_Ratio_EndCap ",6,-0.5,5.5 , 10, 0, 0.002 );


    
    ////////////////////////////////////
    ///  Loop on datasets
    ////////////////////////////////////
    
    if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
    
    ///////----Start looping over datasets -----/////
    bool nlo = false;
    bool isSignal = false;
    
    
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
        
        
       // if(!isData && !isSignal && !btagShape)
        if(!isData && !btagShape)
        {
            // documentation at http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf
            //	   btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv");
            btagcalib = new BTagCalibration("CSVv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
            btagreader = new BTagCalibrationReader(btagcalib, BTagEntry::OP_LOOSE, "mujets","central");
            if(fillBtagHisto)  // before btag reweighting can be applied, you first have to make the histograms
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
            BTagCalibration calib_csvv2("csvv2", "../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv");
            reader_csvv2 = new BTagCalibrationReader(&calib_csvv2, // calibration instance
                                                     BTagEntry::OP_RESHAPING, // operating point
                                                     "iterativefit", // measurement type
                                                     "central"); // systematics type  --> depending on JES up/Down andother reader is needed
            
            
        }
        

        
        
        //// ***************** /////
        /// output Ntuples /////
        /// ***************** ////
        string rootTreesDir = "/user/sabuzeid/FCNC_Study/CMSSW_8_0_24/src/TopBrussels/FCNCAnalysis/Output_Ntuples/";
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
        Float_t Count_cut[15];
        Int_t nCuts;
        // various weights
        
        Int_t nofPosWeights;
        Int_t nofNegWeights;
        Int_t sumW;
        Double_t nloWeight; // for amc@nlo samples
        Double_t MuonIDSF[10];
        Double_t MuonIsoSF[10];
        Double_t sf_muon[10];
        Double_t sf_electron[10];
        Double_t puSF;
        Double_t btagSF;
        Float_t SampleLumiWeight;
        
        
        //// Variables for met
        Int_t PassedMETFilter;
        Float_t met_Pt;
        Float_t met_Eta;
        Float_t met_Phi;
        Float_t corrected_met_Eta ;
        Float_t corrected_met_Phi ;
        Float_t corrected_met_Pt ;
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
        
        //Kinamatic variables for Muons
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
        Float_t bdiscCSVv2_jet_2;
        Float_t bdiscCSVv2_jet_1;
        
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
        globalTree->Branch("sumW", &sumW, "sumW/I");
        globalTree->Branch("nCuts",&nCuts, "nCuts/I");
        globalTree->Branch("Count_cut",&Count_cut,"Count_cut[nCuts]/F");
        
        // define the output trees may I need to make many trees depend on selection cuts
        /////(Integer variables)
        
        InitialTree->Branch("isData",&isData,"isData/I");
        InitialTree->Branch("run_num",&run_num,"run_num/I");
        InitialTree->Branch("evt_num",&evt_num,"evt_num/I");
        InitialTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        InitialTree->Branch("nvtx",&nvtx,"nvtx/I");
        InitialTree->Branch("npu",&npu,"npu/I");
        InitialTree->Branch("puSF",&puSF,"puSF/D");
        InitialTree->Branch("btagSF",&btagSF,"btagSF/D");
        InitialTree->Branch("nLeptons",&nLeptons, "nLeptons/I");//
        
        ////** Tree before applying Z_Mass Veto ** ///
        
        NoZmassVetoTree->Branch("isData",&isData,"isData/I");
        NoZmassVetoTree->Branch("run_num",&run_num,"run_num/I");
        NoZmassVetoTree->Branch("evt_num",&evt_num,"evt_num/I");
        NoZmassVetoTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        NoZmassVetoTree->Branch("nvtx",&nvtx,"nvtx/I");
        NoZmassVetoTree->Branch("npu",&npu,"npu/I");
        NoZmassVetoTree->Branch("puSF",&puSF,"puSF/D");
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
        NoZmassVetoTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        NoZmassVetoTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        NoZmassVetoTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        NoZmassVetoTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        NoZmassVetoTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        NoZmassVetoTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        NoZmassVetoTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        NoZmassVetoTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        
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
        //NoZmassVetoTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // NoZmassVetoTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        NoZmassVetoTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        NoZmassVetoTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        NoZmassVetoTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        NoZmassVetoTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        NoZmassVetoTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        NoZmassVetoTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        NoZmassVetoTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        NoZmassVetoTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        NoZmassVetoTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        NoZmassVetoTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        NoZmassVetoTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        NoZmassVetoTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        NoZmassVetoTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
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
        NoZmassVetoTree->Branch("bdiscCSVv2_jet_1", &bdiscCSVv2_jet_1, "bdiscCSVv2_jet_1/F");
        NoZmassVetoTree->Branch("bdiscCSVv2_jet_2", &bdiscCSVv2_jet_2, "bdiscCSVv2_jet_2/F");
        
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
        NoZmassVetoTree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        NoZmassVetoTree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        
        NoZmassVetoTree->Branch("MT_lep1", &MT_lep1, "MT_lep1/F");
        NoZmassVetoTree->Branch("MT_lep2", &MT_lep2, "MT_lep2/F");
        NoZmassVetoTree->Branch("MT_Elec1", &MT_Elec1, "MT_Elec1/F");
        NoZmassVetoTree->Branch("MT_Elec2", &MT_Elec2, "MT_Elec2/F");
        NoZmassVetoTree->Branch("MT_Mu1", &MT_Mu1, "MT_Mu1/F");
        NoZmassVetoTree->Branch("MT_Mu2", &MT_Mu2, "MT_Mu2/F");

        NoZmassVetoTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");

        
        /////////////////////////////////////////////////////
        
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
        myTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");

        
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
        SSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
        SSLeptonTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");
        
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
        OSLeptonTree->Branch("nLeptons",&nLeptons, "nLeptons/I");
        OSLeptonTree->Branch("SampleLumiWeight", &SampleLumiWeight, "SampleLumiWeight/F");

        
        // Set branches for different Trees
        ////--> Electrons <----////
        InitialTree ->Branch("nElectrons",&nElectrons, "nElectrons/I");
        InitialTree->Branch("pt_electron",pt_electron,"pt_electron[nElectrons]/D");
        InitialTree->Branch("phi_electron",phi_electron,"phi_electron[nElectrons]/D");
        InitialTree->Branch("eta_electron",eta_electron,"eta_electron[nElectrons]/D");
        InitialTree->Branch("eta_superCluster_electron",eta_superCluster_electron,"eta_superCluster_electron[nElectrons]/D");
        InitialTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        InitialTree->Branch("chargedHadronIso_electron",chargedHadronIso_electron,"chargedHadronIso_electron[nElectrons]/D");
        InitialTree->Branch("neutralHadronIso_electron",neutralHadronIso_electron,"neutralHadronIso_electron[nElectrons]/D");
        InitialTree->Branch("photonIso_electron",photonIso_electron,"photonIso_electron[nElectrons]/D");
        InitialTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        InitialTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        InitialTree->Branch("d0_electron",d0_electron,"d0_electron[nElectrons]/D");
        InitialTree->Branch("d0BeamSpot_electron",d0BeamSpot_electron,"d0BeamSpot_electron[nElectrons]/D");
        InitialTree->Branch("sigmaIEtaIEta_electron",sigmaIEtaIEta_electron,"sigmaIEtaIEta_electron[nElectrons]/D");
        InitialTree->Branch("deltaEtaIn_electron",deltaEtaIn_electron,"deltaEtaIn_electron[nElectrons]/D");
        InitialTree->Branch("deltaPhiIn_electron",deltaPhiIn_electron,"deltaPhiIn_electron[nElectrons]/D");
        InitialTree->Branch("hadronicOverEm_electron",hadronicOverEm_electron,"hadronicOverEm_electron[nElectrons]/D");
        InitialTree->Branch("missingHits_electron",missingHits_electron,"missingHits_electron[nElectrons]/I");
        InitialTree->Branch("passConversion_electron",passConversion_electron,"passConversion_electron[nElectrons]/O)");
        InitialTree->Branch("isId_electron",isId_electron,"isId_electron[nElectrons]/O)");
        InitialTree->Branch("isIso_electron",isIso_electron,"isIso_electron[nElectrons]/O)");
        InitialTree->Branch("isEBEEGap",isEBEEGap,"isEBEEGap[nElectrons]/O)");
        InitialTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
        
        InitialTree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        InitialTree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        InitialTree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        InitialTree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        InitialTree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        InitialTree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        InitialTree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        InitialTree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        InitialTree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        InitialTree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");

        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
       // myTree->Branch("ElectronSF",&ElectronSF,"ElectronSF[nElectrons]/D");
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
        SSLeptonTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
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
        OSLeptonTree->Branch("sf_electron",sf_electron,"sf_electron[nElectrons]/D");
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
        InitialTree->Branch("nMuons",&nMuons, "nMuons/I");
        InitialTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        InitialTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //InitialTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
        // InitialTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
        InitialTree->Branch("pt_muon",pt_muon,"pt_muon[nMuons]/D");
        InitialTree->Branch("phi_muon",phi_muon,"phi_muon[nMuons]/D");
        InitialTree->Branch("eta_muon",eta_muon,"eta_muon[nMuons]/D");
        InitialTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        InitialTree->Branch("chargedHadronIso_muon",chargedHadronIso_muon,"chargedHadronIso_muon[nMuons]/D");
        InitialTree->Branch("neutralHadronIso_muon",neutralHadronIso_muon,"neutralHadronIso_muon[nMuons]/D");
        InitialTree->Branch("photonIso_muon",photonIso_muon,"photonIso_muon[nMuons]/D");
        InitialTree->Branch("isId_muon",isId_muon,"isId_muon[nMuons]/O");
        InitialTree->Branch("isIso_muon",isIso_muon,"isIso_muon[nMuons]/O");
        InitialTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        InitialTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        //        InitialTree->Branch("d0_muon",d0_muon,"d0_muon[nMuons]/D");
        //        InitialTree->Branch("d0BeamSpot_muon",d0BeamSpot_muon,"d0BeamSpot_muon[nMuons]/D");
        InitialTree->Branch("sf_muon",sf_muon,"sf_muon[nMuons]/D");
        InitialTree->Branch("pt_1st_Muon",&pt_1st_Muon,"pt_1st_Muon/F");
        InitialTree->Branch("phi_1st_Muon",&phi_1st_Muon,"phi_1st_Muon/F");
        InitialTree->Branch("eta_1st_Muon",&eta_1st_Muon,"eta_1st_Muon/F");
        InitialTree->Branch("E_1st_Muon",&E_1st_Muon,"E_1st_Muon/F");
        InitialTree->Branch("charge_1st_Muon", &charge_1st_Muon, "charge_1st_Muon/F");
        InitialTree->Branch("pt_2nd_Muon",&pt_2nd_Muon,"pt_2nd_Muon/F");
        InitialTree->Branch("phi_2nd_Muon",&phi_2nd_Muon,"phi_2nd_Muon/F");
        InitialTree->Branch("eta_2nd_Muon",&eta_2nd_Muon,"eta_2nd_Muon/F");
        InitialTree->Branch("E_2nd_Muon",&E_2nd_Muon,"E_2nd_Muon/F");
        InitialTree->Branch("charge_2nd_Muon", &charge_2nd_Muon, "charge_2nd_Muon/F");

        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/D");
        myTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/D");
        //myTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
       // myTree->Branch("MuonTrigSFv3",&MuonTrigSFv3,"MuonTrigSFv3[nMuons]/D");
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
        SSLeptonTree->Branch("MuonIDSF",&MuonIDSF,"MuonIDSF[nMuons]/F");
        SSLeptonTree->Branch("MuonIsoSF",&MuonIsoSF, "MuonIsoSF[nMuons]/F");
        //SSLeptonTree->Branch("MuonTrigSFv2",&MuonTrigSFv2,"MuonTrigSFv2[nMuons]/D");
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
        myTree->Branch("bdiscCSVv2_jet_1", &bdiscCSVv2_jet_1, "bdiscCSVv2_jet_1/F");
        myTree->Branch("bdiscCSVv2_jet_2", &bdiscCSVv2_jet_2, "bdiscCSVv2_jet_2/F");
        
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
        myTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        myTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        myTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        myTree->Branch("corrected_met_Pt", &corrected_met_Pt, "corrected_met_Pt/F");
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

        
        InitialTree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        InitialTree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        InitialTree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        
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

        
        
        
        
        InitialTree->Branch("corrected_met_Pt", &corrected_met_Pt, "corrected_met_Pt/F");
        InitialTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        InitialTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        SSLeptonTree->Branch("corrected_met_Pt", &corrected_met_Pt, "corrected_met_Pt/F");
        SSLeptonTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        SSLeptonTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        OSLeptonTree->Branch("corrected_met_Pt", &corrected_met_Pt, "corrected_met_Pt/F");
        OSLeptonTree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        OSLeptonTree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        ////--> bJets <----////
        myTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        myTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        myTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        myTree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        myTree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        myTree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        
        InitialTree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        InitialTree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        InitialTree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        InitialTree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        InitialTree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        InitialTree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        
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
//        myTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
//        myTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
//        myTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
//        myTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
//        myTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
//        myTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
//        myTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
//        myTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
//        myTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
//        myTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
//        
//        
//        SSLeptonTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
//        SSLeptonTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
//        SSLeptonTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
//        SSLeptonTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
//        SSLeptonTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
//        SSLeptonTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
//        SSLeptonTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
//        SSLeptonTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
//        SSLeptonTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
//        SSLeptonTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
//        
//        OSLeptonTree->Branch("nMCParticles",&nMCParticles,"nMCParticles/I");
//        OSLeptonTree->Branch("mc_status",&mc_status,"mc_status[nMCParticles]/I");
//        OSLeptonTree->Branch("mc_pdgId",&mc_pdgId,"mc_pdgId[nMCParticles]/I");
//        OSLeptonTree->Branch("mc_mother",&mc_mother,"mc_mother[nMCParticles]/I");
//        OSLeptonTree->Branch("mc_granny",&mc_granny,"mc_granny[nMCParticles]/I");
//        OSLeptonTree->Branch("mc_pt",&mc_pt,"mc_pt[nMCParticles]/D");
//        OSLeptonTree->Branch("mc_phi",&mc_phi,"mc_phi[nMCParticles]/D");
//        OSLeptonTree->Branch("mc_eta",&mc_eta,"mc_eta[nMCParticles]/D");
//        OSLeptonTree->Branch("mc_E",&mc_E,"mc_E[nMCParticles]/D");
//        OSLeptonTree->Branch("mc_M",&mc_M,"mc_M[nMCParticles]/D");
        
        
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
        myTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        
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
        SSLeptonTree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
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
        SSLeptonTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        
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
        OSLeptonTree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
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
        OSLeptonTree->Branch("Mass_JetPair_arr",Mass_JetPair_arr,"Mass_JetPair_arr[nJetpair]/D");
        
        ///// MVA Tree /////
        
        MVATree->Branch("Ht",&Ht,"Ht/F");
        MVATree->Branch("St",&St,"St/F");
        MVATree->Branch("DeltaR_2L",&DeltaR_2L,"DeltaR_2L/F");
        MVATree->Branch("DeltaPhi_2L",&DeltaPhi_2L,"DeltaPhi_2L/F");
        MVATree->Branch("invMass_2L",&invMass_2L,"invMass_2L/F");
        MVATree->Branch("DeltaR_Mu0b0_DiMu",&DeltaR_Mu0b0_DiMu,"DeltaR_Mu0b0_DiMu/F");
        MVATree->Branch("DeltaR_Mu1b0_DiMu",&DeltaR_Mu1b0_DiMu,"DeltaR_Mu1b0_DiMu/F");
        MVATree->Branch("DeltaR_Elec0b0_DiElec",&DeltaR_Elec0b0_DiElec,"DeltaR_Elec0b0_DiElec/F");
        MVATree->Branch("DeltaR_Elec1b0_DiElec",&DeltaR_Elec1b0_DiElec,"DeltaR_Elec1b0_DiElec/F");
        MVATree->Branch("DeltaR_Mu0b0_DiMuEl",&DeltaR_Mu0b0_DiMuEl,"DeltaR_Mu0b0_DiMuEl/F");
        MVATree->Branch("DeltaR_Mu1b0_DiElMu",&DeltaR_Mu1b0_DiElMu,"DeltaR_Mu1b0_DiElMu/F");
        MVATree->Branch("DeltaR_Elec1b0_DiMuEl",&DeltaR_Elec1b0_DiMuEl,"DeltaR_Elec1b0_DiMuEl/F");
        MVATree->Branch("DeltaR_Elec0b0_DiElMu",&DeltaR_Elec0b0_DiElMu,"DeltaR_Elec0b0_DiElMu/F");
        MVATree->Branch("Mass_WJetPair",&Mass_WJetPair,"Mass_WJetPair/F");
        MVATree->Branch("Mass_JetPair",&Mass_JetPair,"Mass_JetPair/F");
        
        MVATree->Branch("nCSVLbJets",&nCSVLbJets,"nCSVLbJets/I");
        MVATree->Branch("nCSVMbJets",&nCSVMbJets,"nCSVMbJets/I");
        MVATree->Branch("nCSVTbJets",&nCSVTbJets,"nCSVTbJets/I");
        MVATree->Branch("nNonCSVLbJets",&nNonCSVLbJets,"nNonCSVLbJets/I");
        MVATree->Branch("nNonCSVMbJets",&nNonCSVMbJets,"nNonCSVMbJets/I");
        MVATree->Branch("nNonCSVTbJets",&nNonCSVTbJets,"nNonCSVTbJets/I");
        MVATree->Branch("MVA_nJets", &MVA_nJets, "MVA_nJets/F");
        MVATree->Branch("MVA_nCSVLbtagJets", &MVA_nCSVLbtagJets, "MVA_nCSVLbtagJets/F");

        
        MVATree->Branch("met_Pt", &met_Pt, "met_Pt/F");
        MVATree->Branch("met_Eta", &met_Eta,"met_Eta/F");
        MVATree->Branch("met_Phi", &met_Phi, "met_Phi/F");
        MVATree->Branch("corrected_met_Pt", &corrected_met_Pt, "corrected_met_Pt/F");
        MVATree->Branch("corrected_met_Eta", &corrected_met_Eta,"corrected_met_Eta/F");
        MVATree->Branch("corrected_met_Phi", &corrected_met_Phi, "corrected_met_Phi/F");
        
        MVATree->Branch("pt_1st_Electron",&pt_1st_Electron,"pt_1st_Electron/F");
        MVATree->Branch("pt_2nd_Electron",&pt_2nd_Electron,"pt_2nd_Electron/F");
        MVATree->Branch("phi_1st_Electron",&phi_1st_Electron,"phi_1st_Electron/F");
        MVATree->Branch("phi_2nd_Electron",&phi_2nd_Electron,"phi_2nd_Electron/F");
        MVATree->Branch("eta_1st_Electron",&eta_1st_Electron,"eta_1st_Electron/F");
        MVATree->Branch("eta_2nd_Electron",&eta_2nd_Electron,"eta_2nd_Electron/F");
        MVATree->Branch("E_1st_Electron",&E_1st_Electron,"E_1st_Electron/F");
        MVATree->Branch("E_2nd_Electron",&E_2nd_Electron,"E_2nd_Electron/F");
        MVATree->Branch("charge_1st_Electron", &charge_1st_Electron, "charge_1st_Electron/F");
        MVATree->Branch("charge_2nd_Electron", &charge_2nd_Electron, "charge_2nd_Electron/F");
        MVATree->Branch("Match_MC_Ht",&Match_MC_Ht,"Match_MC_Ht/F");
        MVATree->Branch("Match_MC_St",&Match_MC_St,"Match_MC_St/F");
        MVATree->Branch("Match_MC_DeltaR_2L",&Match_MC_DeltaR_2L,"Match_MC_DeltaR_2L/F");
        MVATree->Branch("Match_MC_DeltaPhi_2L",&Match_MC_DeltaPhi_2L,"Match_MC_DeltaPhi_2L/F");
        MVATree->Branch("Match_MC_invMass_2L",&Match_MC_invMass_2L,"Match_MC_invMass_2L/F");
        MVATree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMu",&Match_SM_MC_DeltaR_Mu0b0_DiMu,"Match_SM_MC_DeltaR_Mu0b0_DiMu/F");
        MVATree->Branch("Match_MC_DeltaR_Mu1b0_DiMu",&Match_MC_DeltaR_Mu1b0_DiMu,"Match_MC_DeltaR_Mu1b0_DiMu/F");
        MVATree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElec",&Match_SM_MC_DeltaR_Elec0b0_DiElec,"Match_SM_MC_DeltaR_Elec0b0_DiElec/F");
        MVATree->Branch("Match_MC_DeltaR_Elec1b0_DiElec",&Match_MC_DeltaR_Elec1b0_DiElec,"Match_MC_DeltaR_Elec1b0_DiElec/F");
        MVATree->Branch("Match_SM_MC_DeltaR_Mu0b0_DiMuEl",&Match_SM_MC_DeltaR_Mu0b0_DiMuEl,"Match_SM_MC_DeltaR_Mu0b0_DiMuEl/F");
        MVATree->Branch("Match_MC_DeltaR_Mu1b0_DiElMu",&Match_MC_DeltaR_Mu1b0_DiElMu,"Match_MC_DeltaR_Mu1b0_DiElMu/F");
        MVATree->Branch("Match_MC_DeltaR_Elec1b0_DiMuEl",&Match_MC_DeltaR_Elec1b0_DiMuEl,"Match_MC_DeltaR_Elec1b0_DiMuEl/F");
        MVATree->Branch("Match_SM_MC_DeltaR_Elec0b0_DiElMu",&Match_SM_MC_DeltaR_Elec0b0_DiElMu,"Match_SM_MC_DeltaR_Elec0b0_DiElMu/F");
        MVATree->Branch("Match_FCNC_MC_Mass_WJetPair",&Match_FCNC_MC_Mass_WJetPair,"Match_FCNC_MC_Mass_WJetPair/F");
        MVATree->Branch("Match_FCNC_MC_Mass_JetPair",&Match_FCNC_MC_Mass_JetPair,"Match_FCNC_MC_Mass_JetPair/F");


        
        
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

        //nEvents[d] = 0;
        nbEvents = 0;
        int previousRun = -1;
        int currentRun;
        bool Btagged = false;
        nofPosWeights = 0;
        nofNegWeights = 0;
        int nbEvents_0 = 0;
        int nbEvents_1PU = 0;
        int nbEvents_2BTag = 0;
        int nbEvents_3Trig=0;
        int nbEvents_4GPV = 0;
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
        ////Filters
        bool passedMET = false;
        bool   HBHEnoise = false;
        bool   HBHEIso = false;
        bool   CSCTight = false;
        bool badchan = false;
        bool badmu = false;
        bool EcalDead = false;
        
        
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
        vector < TRootElectron* > selectedElectrons;
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
        vector<int> electroncharge;
        
        
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
           // nEvents[d]++;
            vertex.clear();
            init_muons.clear();
            init_jets.clear();
            init_jets_corrected.clear();
            genjets.clear();
            mets.clear();
            nCuts = 0;
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
            
            ///// define Filters ////
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
           // cout << "the number of events after calculating negative weights is =  " << datasets[d]->NofEvtsToRunOver() << endl;
         
           //// Trigger.checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
            if (Apply_HLT_Triggers)
            {
                //DilepTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
               // itrigger = DilepTrigger->checkIfFired();
                
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
                jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
            }
            if(applyJES && !isData)
            {
                jetTools->correctJets(init_jets,event->fixedGridRhoFastjetAll() ,false);
            }
            
            /////////////////////////
            ///  EVENT SELECTION  ///
            /////////////////////////
            selectedElectrons.clear();
            selectedMuons.clear();
            selectedOrgiJets.clear();
            mcParticles.clear();
            selectedLooseMuons.clear();
            selectedLooseElectrons.clear();

            
            //Declare selection instance
            
            Run2Selection selection(init_jets, init_muons, init_electrons, mets, rho); ///rho = event->fixedGridRhoFastjetAll(); and it is used in electron isolation
            
            
            //// choose good primary vertex
            bool isGoodPV = selection.isPVSelected(vertex, 4 , 24. ,2.); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
            
            // Jets Selection
            selectedOrgiJets = selection.GetSelectedJets(jet_pt_cut, jet_eta_cut, true , "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
            
            
            /// --- Muons Selection -- ///
            selectedMuons = selection.GetSelectedMuons(mu_pt_cut , mu_eta_cut , mu_iso_cut ,"Tight","Spring15");// GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
            selectedLooseMuons = selection.GetSelectedMuons(mu_pt_cut , mu_eta_cut , mu_iso_cut ,"Loose","Spring15");
            
            //// --- Electron Selection --- ///
            selectedElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut , "Tight" , "Spring16_80X", true, true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
            selectedLooseElectrons = selection.GetSelectedElectrons(el_pt_cut, el_eta_cut , "Loose" , "Spring16_80X", true, true);
            
           if(debug)cout << "the nb of selected electrons is " << selectedElectrons.size()<<endl;
           if(debug) cout << "the number of events before jet cleaning is =  " << datasets[d]->NofEvtsToRunOver() << endl;

            //// sorting objects in the the event according to Pt
            
            sort(selectedJets.begin(), selectedJets.end(),HighestPt());
            sort(selectedMuons.begin(), selectedMuons.end(), HighestPt());
            sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt());
            
            
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
            Double_t Lep_scaleFactor = 1.;
            //Double_t bTag_scaleFactor =1.;
            float muon1SF, muon2SF,muonID1SF, muonID2SF,muonIso1SF, muonIso2SF,electron1SF, electron2SF;
            muon1SF =  muon2SF = electron1SF = electron2SF = 1.;
            
            histo1D["h_cutFlow"]->Fill(0., scaleFactor*lumiWeight); //// fill histogram before applying any scaling factors or triggers
            nCuts++;
            nbEvents_0++; //// add to the frist bin of Count_cut branch before applying any scaling factors or triggers
            
            ////initial histograms after applying scaling factors or triggers
            
           // histo1D["h_initial_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
            
            //histo1D["h_initial_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
            // histo1D["h_initial_met"]->Fill(met_pt, scaleFactor*lumiWeight);
           
            
            
            /////applyNegWeightCorrection
            
            if(hasNegWeight && applyNegWeightCorrection && !isData) scaleFactor *= -1.;
            
            ////// PU SF
            
            float puWeight = 1;
            // if (verbose>2) cout << " isData =  "<< isData << endl;
            if (ApplyPU_SF)
            {
                
                if(!isData && !isSignal)
                {
                    puWeight = LumiWeights.ITweight((int)event->nTruePU()); // simplest reweighting, just use reconstructed number of PV. faco
                    PU_scaleFactor =puWeight;
                    if (debug) cout << "while puWeight  =  "<< puWeight << endl;
                }
                else if (isData) {PU_scaleFactor =1;}
                 if (debug) cout << "PU_scaleFactor is " << PU_scaleFactor << endl;
                scaleFactor *= PU_scaleFactor;
                
            }
           
            nCuts++;
            nbEvents_1PU++; //// add to the frist bin of Count_cut branch After applying PU SF
         
            if (debug) cout << "puSF " << PU_scaleFactor << endl;

            
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
         //   histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
            nCuts++;
            nbEvents_2BTag++;
           ///// cout << "the nCut Value After applying bTag SF is =  " << nCuts << "and the nbEvents (nbEvents_2BTag) =  " << nbEvents_2BTag << endl;
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
            DeltaR_Mu0b0_DiMu= DeltaR_Mu1b0_DiMu= DeltaR_Mu0b0_DiMuEl= DeltaR_Mu1b0_DiElMu= DeltaR_Elec0b0_DiElec= DeltaR_Elec1b0_DiElec = DeltaR_Elec0b0_DiElMu = DeltaR_Elec1b0_DiMuEl =0.;
        
            Mass_WJetPair = 0.;
            Mass_JetPair = 0.;
            
              /// Trigger
            if(Apply_HLT_Triggers) trigged = itrigger; //trigged = treeLoader.EventTrigged(itrigger);
           // if(dName.find("FCNC")!=string::npos) trigged = true;
           if(debug) cout << "the number of events After Applying trigger is =  " << datasets[d]->NofEvtsToRunOver() << endl;
           if(debug) cout << "trigged = " << trigged << endl;
            if (!trigged) continue;
            nCuts++;
            nbEvents_3Trig++;
          /////  cout << "the nCut Value After applying Triggers is =  " << nCuts << "and the nbEvents (nbEvents_2) =  " << nbEvents_2 << endl;
            histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight);
            
            ///// *** Apply Met Filters *** ///
            if(HBHEnoise && HBHEIso && CSCTight && EcalDead  && isGoodPV && badchan && badmu) passedMET = true;
            PassedMETFilter = passedMET;
            
            //if (!isGoodPV) continue;
            if (!isGoodPV && !PassedMETFilter) continue;
            nCuts++;
            nbEvents_4GPV++;
            if(debug)cout << "the nCut Value After good PV is =  " << nCuts << "and the nbEvents (nbEvents_4GPV) =  " << nbEvents_4GPV << endl;
            if(verbose>2) cout << "GoodPV" << endl;
            histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
            
            int numLep = selectedElectrons.size()+ selectedMuons.size();
            int numLooseLep = selectedLooseElectrons.size()+ selectedLooseMuons.size();
            nMuons = 0;
            nLeptons = 0;
            nElectrons = 0;
            //////--- Filling variable branches for Met  ----- //////
            
            float met_px = mets[0]->Px();
            float met_py = mets[0]->Py();
            
            met_Pt = sqrt(met_px*met_px + met_py*met_py);
            met_Phi = mets[0]->Phi();
            met_Eta = mets[0]->Eta();
            
            if(applyJES) // jer doesn't need to be applied ||  applyJER)) --> smeared type-1 corrected MET,  NOW only yes --> Type 1 corrected MET
            {
                jetTools->correctMETTypeOne(init_jets, mets[0], isData);
              //  METon = 1;
                //  if JES applied: replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which JEC is applied
                //  if JER applied:  replaces the vector sum of transverse momenta of particles which can be clustered as jets with the vector sum of the transverse momenta of the jets to which smearing is applied.
                // type 1 correction / sleard pmet correction
                
            }
            float corrected_met_px = mets[0]->Px();
            float corrected_met_py = mets[0]->Py();
            corrected_met_Pt = sqrt(corrected_met_px*corrected_met_px + corrected_met_py*corrected_met_py);
            
            puSF = PU_scaleFactor;
            btagSF = btagWeight;
            
            
           // InitialTree -> Fill();
            nCuts++;
            nbEvents_5LepSF++;
           
            if (numLep==2 && numLooseLep == 2)
            {
                histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
                
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
                        
                        MuonIDSF[nMuons] = muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                        MuonIsoSF[nMuons] =  muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                        sf_muon[nMuons]= muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                        Muon_scaleFactor = muonSFWeightIso->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0)* muonSFWeightID->at(selectedMuons[selmu]->Eta(), selectedMuons[selmu]->Pt(), 0);
                    }
                    else
                    {
                        MuonIDSF[nMuons] = 1.;
                        MuonIsoSF[nMuons] = 1.;
                        sf_muon[nMuons] = 1.;
                        Muon_scaleFactor = 1.;
                        
                    }
                    
                    nMuons++;
                }
                //scaleFactor *= Muon_scaleFactor;
                
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
                    
                    if (ApplyElec_SF && !isData)
                    {
                        sf_electron[nElectrons]= electronSFWeight->at(selectedElectrons[selelec]->Eta(),selectedElectrons[selelec]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[selelec]->Eta(),selectedElectrons[selelec]->Pt(),0);
                        
                    }
                    else sf_electron[nElectrons]= 1.;
                    //electronSF *= sf_electron[nElectrons];
                   // Lep_scaleFactor *= sf_electron[nElectrons];
                    nElectrons++;
                    //Elec_scaleFactor *=sf_electron;
                }
                nLeptons = nMuons + nElectrons;
                
                //////--- Filling variable branches for Jets  ----- //////
                nJets = 0;
                bdiscCSVv2_jet_1 = bdiscCSVv2_jet_2 = 0;
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
                MVA_nJets = nJets;
                if(selectedJets.size()>0 ) bdiscCSVv2_jet_1 = selectedJets[0]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                if(selectedJets.size()>1 ) bdiscCSVv2_jet_2 = selectedJets[1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                
                //////--- Filling variable branches for b-tagged Jets  ----- //////
                nCSVTbJets = selectedBCSVTJets.size();
                nCSVMbJets = selectedBCSVMJets.size();
                nCSVLbJets = selectedBCSVLJets.size();
                MVA_nCSVLbtagJets = nNonCSVLbJets;

                
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
                            electron1SF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                            electron2SF = electronSFWeight->at(selectedElectrons[1]->Eta(),selectedElectrons[1]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[1]->Eta(),selectedElectrons[1]->Pt(),0);
                            Lep_scaleFactor = electron1SF * electron2SF;
                            Elec_scaleFactor = electron1SF*electron2SF;
                           // scaleFactor *= electronSF;
                            if(selectedElectrons.size()>0) pt_1st_Electron = selectedElectrons[0]->Pt();
                            if(selectedElectrons.size()>1) pt_2nd_Electron = selectedElectrons[1]->Pt();
                            if(selectedElectrons.size()>0) eta_1st_Electron = selectedElectrons[0]->Eta();
                            if(selectedElectrons.size()>1) eta_2nd_Electron = selectedElectrons[1]->Eta();
                            if(selectedElectrons.size()>0) phi_1st_Electron = selectedElectrons[0]->Phi();
                            if(selectedElectrons.size()>1) phi_2nd_Electron = selectedElectrons[1]->Phi();
                            if(selectedElectrons.size()>0) phi_1st_Electron = selectedElectrons[0]->E();
                            if(selectedElectrons.size()>1) phi_2nd_Electron = selectedElectrons[1]->E();
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
                        
                        // Zmass_Window = fabs(Zmass - InvMass_ll);
                        if (InvMass_ll >12.) //&& Zmass_Window >= 15.)
                        {
                            diMuon = true;
                            invMass_2L = InvMass_ll;
                            muon1SF = muonSFWeightIso->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)* muonSFWeightID->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                            muon2SF = muonSFWeightIso->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0)* muonSFWeightID->at(selectedMuons[1]->Eta(), selectedMuons[1]->Pt(), 0);
                            Lep_scaleFactor = muon1SF * muon2SF;
                            // scaleFactor *= Muon_scaleFactor;
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
                            
                        }
                        
                        
                    }
                }
                
                else if (Elec_Mu && !Elec_Elec && !Mu_Mu) //// Electron- Muon channel
                {
                    if (selectedMuons.size() == 1 && selectedElectrons.size() ==1)
                    {
                        if (selectedMuons[0]->Pt() > selectedElectrons[0]->Pt() && selectedMuons[0]->Pt()>= 20 && selectedElectrons[0]->Pt()>=15)
                        {
                            diMuE = true;
                            tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            qLepton0 = selectedMuons[0]->charge();
                            qLepton1 = selectedElectrons[0]->charge();
                            MT_lep1 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedMuons[0]);
                            MT_lep2 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedElectrons[0]);

                            
                        }
                        if (selectedElectrons[0]->Pt() > selectedMuons[0]->Pt() && selectedElectrons[0]->Pt()>=25 && selectedMuons[0]->Pt()>= 10)
                        {
                            diEMu = true;
                            tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                            qLepton0 = selectedElectrons[0]->charge();
                            qLepton1 = selectedMuons[0]->charge();
                            MT_lep1 = sqrt(2*selectedElectrons[0]->Pt() * mets[0]->Et() * (1-cos( selectedElectronsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedElectrons[0]);
                            MT_lep2 = sqrt(2*selectedMuons[0]->Pt() * mets[0]->Et() * (1-cos( selectedMuonsTLV[0].DeltaPhi( metsTLV[0] )) ) );
                            selectedLeptonsTLV.push_back(*selectedMuons[0]);
                            
                        }
                        muon1SF = muonSFWeightIso->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0)* muonSFWeightID->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                        electron1SF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0) * electronSFWeightReco->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        Lep_scaleFactor = muon1SF * electron1SF;
                        
                    }
                }
               
                if (diElectron || diMuon || diEMu || diMuE)
                {
                    scaleFactor *= Lep_scaleFactor;
                    nCuts++;
                    nbEvents_6++;
                    histo1D["h_cutFlow"]->Fill(4., scaleFactor*lumiWeight);
                   ////// cout << "the nCut Value After passing 2 leptons is =  " << nCuts << "and the nbEvents (nbEvents_3) =  " << nbEvents_3 << endl;
                    
                    histo2D["2L_Nb_jets_vs_CSVLbjets"]->Fill(selectedJets.size(),selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                    histo2D["2L_Nb_jets_vs_CSVMbjets"]->Fill(selectedJets.size(),selectedBCSVMJets.size(),scaleFactor*lumiWeight);
                    histo2D["2L_Nb_jets_vs_CSVTbjets"]->Fill(selectedJets.size(),selectedBCSVTJets.size(),scaleFactor*lumiWeight);
                    histo2D["h_2L_Lep0_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                    
                    InitialTree -> Fill();
                   if (debug)cout << " Hello selectedJets.size = " << selectedJets.size() <<endl;
                    
                    if (selectedJets.size()>= 3 && selectedJets[0]->Pt() >=30 && selectedJets[1]->Pt() >=30 && selectedJets[2]->Pt() >=30)
                    {
                        nCuts++;
                        nbEvents_7++;
                        
                        histo1D["h_cutFlow"]->Fill(5., scaleFactor*lumiWeight);
                        histo2D["2L_3Jets_Nb_jets_vs_CSVLbjets"]->Fill(selectedJets.size(),selectedBCSVLJets.size(),scaleFactor*lumiWeight);
                        histo2D["2L_3Jets_Nb_jets_vs_CSVMbjets"]->Fill(selectedJets.size(),selectedBCSVMJets.size(),scaleFactor*lumiWeight);
                        histo2D["2L_3Jets_Nb_jets_vs_CSVTbjets"]->Fill(selectedJets.size(),selectedBCSVTJets.size(),scaleFactor*lumiWeight);
                        histo2D["h_2L_3Jets_Lep0_vs_Lep1Pt"]->Fill(tempLepton_0.Pt(),tempLepton_1.Pt(),scaleFactor*lumiWeight);
                       
                   //     cout << " Hello selectedBCSVLJets.size = " << selectedBCSVLJets.size() <<endl;
                        if (selectedBCSVLJets.size()>= 1)
                        {
                            nCuts++;
                            nbEvents_8++;
                            histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
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
                                                cout << "Find FCNC W first Jet" << endl;
                                            }
                                            if (FCNCWJetsCounter ==1)
                                            {
                                                FCNCWJet2_ =JetPartonPair[particlenb];
                                                FCNCWJetsCounter++;
                                                cout << "Find FCNC W second Jet" << endl;
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
                                Match_FCNC_MC_Mass_WJetPair= (temp_Match_WJet1 + temp_Match_WJet2).M();
                                
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
                            for (unsigned int i = 0 ; i<InvMass_JJ.size(); i++)
                            {
                              Mass_JetPair_arr[nJetpair] = InvMass_JJ[i];
                                nJetpair+1;
                            }
                            
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
                            
                           // cout << "the nb of selected electrons is =  " << selectedElectrons.size() <<endl;
                          //  cout << "the nb of selected muons is =  " << selectedMuons.size() <<endl;
                            ////** Determine Charge mis-identification ratio **//

                            invMass_2SSL_Zmass = 0;
                            invMass_2OSL_Zmass =0;
                            if (diElectron && fabs(Zmass - InvMass_ll ) <= 10)
                            {
                                if (qElec0 == qElec1) {
                                    invMass_2SSL_Zmass = InvMass_ll;
                                }else if (qElec0 != qElec1) {
                                    invMass_2OSL_Zmass = InvMass_ll;
                                }
                            }
                           // if (isData && diElectron && fabs(Zmass - InvMass_ll ) >= 10)
                            if (isData && diElectron && fabs(Zmass - InvMass_ll ) <= 10)
                            {
                                int eta_range, Pt_range;
                                eta_range = Pt_range = 0;
                                
                                if (qElec0 == qElec1)
                                {
                                    
                                    for (unsigned iElec = 0 ; iElec < selectedElectrons.size() ; iElec++)
                                    {
                                       // histo1D["h_NSSL_barrel_cutFlow"]->Fill(0);
                                      //  histo1D["h_NSSL_EndCap_cutFlow"]->Fill(0);
                                        if (selectedElectrons[iElec]->Pt() >=50) cout << " At SSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                                        
                                        if (abs(selectedElectrons[iElec]->Eta()) > 0 && abs(selectedElectrons[iElec]->Eta()) <= 1.479)
                                        {
                                            eta_range = 1;
                                            if (selectedElectrons[iElec]->Pt() < 25 )
                                            {
                                                Nb_2SSl_B1++;
                                                cout << " Nb_2SSl_B1 =  " << Nb_2SSl_B1 <<endl;
                                                InBarrel = true;
                                                histo1D["h_NSSL_barrel_cutFlow"]->Fill(1);
                                                Pt_range =1;
                                            }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                            {
                                                Nb_2SSl_B2++;
                                                cout << " Nb_2SSl_B2 =  " << Nb_2SSl_B2 <<endl;
                                                InBarrel = true;
                                                histo1D["h_NSSL_barrel_cutFlow"]->Fill(2);
                                                Pt_range =2;
                                            }else //if (selectedElectrons[iElec]->Pt() >= 50 )
                                            {
                                                Nb_2SSl_B3++;
                                                cout << " Nb_2SSl_B3 =  " << Nb_2SSl_B3 <<endl;
                                                InBarrel = true;
                                                histo1D["h_NSSL_barrel_cutFlow"]->Fill(3);
                                                Pt_range =3;
                                            }
                                        }else if (abs(selectedElectrons[iElec]->Eta()) > 1.479 && abs(selectedElectrons[iElec]->Eta()) <=2.4)
                                        {
                                            eta_range = 2;
                                            if (selectedElectrons[iElec]->Pt() < 25 )
                                            {
                                                Nb_2SSl_E1++;
                                                InEndcap = true;
                                                histo1D["h_NSSL_EndCap_cutFlow"]->Fill(1);
                                                Pt_range =1;
                                            }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                            {
                                                Nb_2SSl_E2++;
                                                InEndcap = true;
                                                histo1D["h_NSSL_EndCap_cutFlow"]->Fill(2);
                                                Pt_range =2;
                                            }else //if (selectedElectrons[iElec]->Pt() >= 50)
                                            {   Nb_2SSl_E3++;
                                                InEndcap = true;
                                                histo1D["h_NSSL_EndCap_cutFlow"]->Fill(3);
                                                Pt_range =3;
                                            }
                                        }
                                        histo2D["h_NSSL_Pt_eta"]->Fill(Pt_range,eta_range);
                                    }
                                   // histo2D["h_NSSL_Pt_eta"]->Fill(Pt_range,eta_range);
                                }
                                if (qElec0 != qElec1)
                                {
                                    for (unsigned iElec = 0; iElec < selectedElectrons.size() ; iElec++)
                                    {
                                      //  histo1D["h_NOSL_barrel_cutFlow"]->Fill(0);
                                      //  histo1D["h_NOSL_EndCap_cutFlow"]->Fill(0);
                                        if (selectedElectrons[iElec]->Pt() >=50) cout << " At OSL the electron pt is =  " << selectedElectrons[iElec]->Pt() << " and with eta =  " << abs(selectedElectrons[iElec]->Eta()) << endl;
                                        
                                        if (abs(selectedElectrons[iElec]->Eta()) > 0 && abs(selectedElectrons[iElec]->Eta()) <= 1.479)
                                        {
                                            eta_range = 1;
                                            if (selectedElectrons[iElec]->Pt() < 25 )
                                            {
                                                Nb_2OSl_B1++;
                                                InBarrel = true;
                                                cout << " Nb_2OSl_B1 =  " << Nb_2OSl_B1 <<endl;
                                                histo1D["h_NOSL_barrel_cutFlow"]->Fill(1);
                                                Pt_range = 1;
                                            }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() < 50 )
                                            {
                                                Nb_2OSl_B2++;
                                                InBarrel = true;
                                                cout << " Nb_2OSl_B2 =  " << Nb_2OSl_B2 <<endl;
                                                histo1D["h_NOSL_barrel_cutFlow"]->Fill(2);
                                                Pt_range = 2;
                                            }else //(selectedElectrons[iElec]->Pt()>= 50)
                                            {
                                                Nb_2OSl_B3++;
                                                InBarrel = true;
                                                cout << " Nb_2OSl_B3 =  " << Nb_2OSl_B3 <<endl;
                                                histo1D["h_NOSL_barrel_cutFlow"]->Fill(3);
                                                Pt_range = 3;
                                            }
                                        }else if (abs(selectedElectrons[iElec]->Eta()) > 1.479 && abs(selectedElectrons[iElec]->Eta()) <=2.4)
                                        {
                                            eta_range = 2;
                                            if (selectedElectrons[iElec]->Pt() < 25 )
                                            {
                                                Nb_2OSl_E1++;
                                                InEndcap = true;
                                                histo1D["h_NOSL_EndCap_cutFlow"]->Fill(1);
                                                Pt_range =1;
                                            }else if (25 <= selectedElectrons[iElec]->Pt() && selectedElectrons[iElec]->Pt() <= 50 )
                                            {
                                                Nb_2OSl_E2++;
                                                InEndcap = true;
                                                histo1D["h_NOSL_EndCap_cutFlow"]->Fill(2);
                                                Pt_range=2;
                                            }else //(selectedElectrons[iElec]->Pt() >= 50 )
                                            {
                                                Nb_2OSl_E3++;
                                                InEndcap = true;
                                                histo1D["h_NOSL_EndCap_cutFlow"]->Fill(3);
                                                Pt_range =3;
                                            }
                                        }
                                       histo2D["h_NOSL_Pt_eta"]->Fill(Pt_range,eta_range);
                                    }
                                   // histo2D["h_NOSL_Pt_eta"]->Fill(Pt_range,eta_range);
                                }
                                
                            }
                            
                            NoZmassVetoTree ->Fill();
                            
                            
                            if(debug)cout << "the mass of dilepton before Zmass veto =  " <<InvMass_ll << endl;
                            if(debug)cout << "the mass of WJets before Zmass veto =  " << Mass_WJetPair << endl;
                            if (fabs(Zmass - InvMass_ll ) >= 15)
                            {
                                histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                                if(debug)cout << "the mass of dilepton After Zmass veto =  " <<InvMass_ll << endl;
                                if(debug)cout << "And the mass difference (Zmass - InvMass_ll)  =  " <<InvMass_ll << endl;
                                if(debug)cout << "the mass of WJets After Zmass veto =  " << Mass_WJetPair << endl;
                                myTree->Fill();
                                nbEvents_9++;
                                MVATree->Fill();
                                ////////////****************************************///////////
                                ///// Deviding into 2 categories: SameSign & Opposite sign ////
                                //////////////****************************************/////////
                                
                                if (diElectron && qElec0 == qElec1) // in case of Same Sign dielectron
                                {
                                    SSdiLepton = true;
                                    
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                    // histo1D["h_2SSL_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
                                  //  histo1D["h_2SSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                    
                                }
                                if (diElectron && qElec0 != qElec1) // in case of Opposite Sign dielectron
                                {
                                    OSdiLepton = true;
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                //    histo1D["h_2OSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                    
                                }
                                if (diMuon && !diElectron && qMu0 == qMu1) // in case of Same Sign diMuon
                                {
                                    SSdiLepton = true;
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                //    histo1D["h_2SSL_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                }
                                if (diMuon && !diElectron && qMu0 != qMu1) // in case of Opposite Sign diMuon
                                {
                                    OSdiLepton = true;
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                 //   histo1D["h_2OSL_mMu0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
                                }
                                if (diEMu || diMuE && !diMuon && !diElectron && qLepton0 == qLepton1)
                                {
                                    SSdiLepton = true;
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                }
                                if (diEMu || diMuE && !diMuon && !diElectron && qLepton0 != qLepton1) // in case of Opposite Sign dilepton
                                {
                                    OSdiLepton = true;
                                    Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
                                }
                                
                                if (SSdiLepton && !OSdiLepton)
                                {
                                    nCuts++;
                                    nbEvents_10++;
                                    ///// cout << "the nCut Value After passing 2 leptons + >= 3Jets +>= 1bjet + 2SSL is =  " << nCuts << "and the nbEvents (nbEvents_6) =  " << nbEvents_6 << endl;
                                    histo1D["h_cutFlow"]->Fill(8., scaleFactor*lumiWeight);
                                    for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                    {
                                        sum_jet_PT+= selectedJets[ijet]->Pt();
                                    }
                                    Ht_AllJets=sum_jet_PT;
                                    St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                                   // histo1D["h_2SSL_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                                   // histo1D["h_2SSL_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                                    histo2D["h_2SSL_3Jets1b_Ht_vs_metPt"]->Fill(Ht_AllJets,met_Pt,scaleFactor*lumiWeight);
                                    histo2D["h_2SSL_2lDeltaPhi_vs_metPt"]->Fill(DeltaPhi_2L,met_Pt,scaleFactor*lumiWeight);
                                    SSLeptonTree->Fill();
                                }
                                if (OSdiLepton && !SSdiLepton)
                                {
                                    nCuts++;
                                    nbEvents_11++;
                                    /////  cout << "the nCut Value After passing 2 leptons + >= 3Jets +>= 1bjet  + 2OSSL is =  " << nCuts << "and the nbEvents (nbEvents_7) =  " << nbEvents_7 << endl;
                                    histo1D["h_cutFlow"]->Fill(9., scaleFactor*lumiWeight);
                                    for (int ijet = 0; ijet<selectedJets.size(); ijet++)
                                    {
                                        sum_jet_PT+= selectedJets[ijet]->Pt();
                                    }
                                    Ht_AllJets=sum_jet_PT;
                                    St_AllJets_Leps= Ht_AllJets+Sum_Leptons_Pt;
                                  //  histo1D["h_2OSL_Ht_AllJets"]->Fill(Ht_AllJets,scaleFactor*lumiWeight);
                                   // histo1D["h_2OSL_St_AllJets+Lep"]->Fill(St_AllJets_Leps,scaleFactor*lumiWeight);
                                    histo2D["h_2OSL_3Jets1b_Ht_vs_metPt"]->Fill(Ht_AllJets,met_Pt,scaleFactor*lumiWeight);
                                    histo2D["h_2OSL_2lDeltaPhi_vs_metPt"]->Fill(DeltaPhi_2L,met_Pt,scaleFactor*lumiWeight);
                                    
                                    OSLeptonTree->Fill();
                                }
                                
                                eventSelected = true;
                                nofSelectedEvents++;
                                
                            }
                            
                            
                            
                            
                        }//2L+>=2Jets+>=1CSVLB
                        
                        
                    }//2L+>=3Jets
                    
                }//2L
                
            }//2L
            
            
            
            
            
//            if (! eventSelected )
//            {
//                continue;
//            }
            
            if(verbose>3)cout << "filling the tree" << endl;
            //myTree->Fill();
           // nofSelectedEvents++;
           // infoFile << "|" << evt_num << "|"  << "|"  <<Channel << "|"  << endl;
            
            //if (verbose > 0) cout << "  Event " << ievt << " is selected" << endl;
        } // end the loop over events
        NSSL_Barrel_Vec.clear();
        NOSL_Barrel_Vec.clear();
        NSSL_EndCap_Vec.clear();
        NOSL_EndCap_Vec.clear();
        Charge_misId_BVec.clear();
        Charge_misId_EVec.clear();
        
       // cout<< " the value of the SF  is =   " << scaleFactor << "and the lumiWeight is =  " << lumiWeight <<endl;
        if (isData)
        {
            if (InBarrel)
            {
                NSSL_Barrel_Vec.push_back(Nb_2SSl_B1);
                NSSL_Barrel_Vec.push_back(Nb_2SSl_B2);
                NSSL_Barrel_Vec.push_back(Nb_2SSl_B3);
                NOSL_Barrel_Vec.push_back(Nb_2OSl_B1);
                NOSL_Barrel_Vec.push_back(Nb_2OSl_B2);
                NOSL_Barrel_Vec.push_back(Nb_2OSl_B3);
                
                
            }
            if (InEndcap)
            {
                NSSL_EndCap_Vec.push_back(Nb_2SSl_E1);
                NSSL_EndCap_Vec.push_back(Nb_2SSl_E2);
                NSSL_EndCap_Vec.push_back(Nb_2SSl_E3);
                NOSL_EndCap_Vec.push_back(Nb_2OSl_E1);
                NOSL_EndCap_Vec.push_back(Nb_2OSl_E2);
                NOSL_EndCap_Vec.push_back(Nb_2OSl_E3);
                
            }
        }
        
        for (unsigned i = 0 ; i< NSSL_Barrel_Vec.size() ; i++ ) {histo2D["h_N2SSL_barrel_Pt"]->Fill(i+1 , NSSL_Barrel_Vec[i]);}
        for (unsigned i = 0 ; i< NOSL_Barrel_Vec.size() ; i++ ) {histo2D["h_N2OSL_barrel_Pt"]->Fill(i+1 , NOSL_Barrel_Vec[i]);}
        for (unsigned i = 0 ; i< NSSL_EndCap_Vec.size() ; i++ ) {histo2D["h_N2SSL_EndCap_Pt"]->Fill(i+1 , NSSL_EndCap_Vec[i]);}
        for (unsigned i = 0 ; i< NOSL_EndCap_Vec.size() ; i++ ) {histo2D["h_N2OSL_EndCap_Pt"]->Fill(i+1 , NOSL_EndCap_Vec[i]);}
        
        if (InBarrel)
        {
            if(Nb_2OSl_B1 != 0){Charge_misId_B1 = Nb_2SSl_B1/(2*Nb_2OSl_B1);Charge_misId_BVec.push_back(Charge_misId_B1);}
            if(Nb_2OSl_B2 != 0){Charge_misId_B2 = Nb_2SSl_B2/(2*Nb_2OSl_B2);Charge_misId_BVec.push_back(Charge_misId_B2);}
            if(Nb_2OSl_B3 != 0){Charge_misId_B3 = Nb_2SSl_B3/(2*Nb_2OSl_B3);Charge_misId_BVec.push_back(Charge_misId_B3);}
            
//            cout << " the charge mis_Id ratio at B1 is =   " << Charge_misId_B1 << endl;
//            cout << " the charge mis_Id ratio at B2 is =   " << Charge_misId_B2 << endl;
//            cout << " the charge mis_Id ratio at B3 is =   " << Charge_misId_B3 << endl;
            
        }
        if (InEndcap)
        {
            if(Nb_2OSl_E1 != 0){Charge_misId_E1 = Nb_2SSl_E1/(2*Nb_2OSl_E1);Charge_misId_EVec.push_back(Charge_misId_E1);}
            if(Nb_2OSl_E2 != 0){Charge_misId_E2 = Nb_2SSl_E2/(2*Nb_2OSl_E2);Charge_misId_EVec.push_back(Charge_misId_E2);}
            if(Nb_2OSl_E3 != 0){Charge_misId_E3 = Nb_2SSl_E3/(2*Nb_2OSl_E3);Charge_misId_EVec.push_back(Charge_misId_E3);}
//            cout << " the charge mis_Id ratio at E1 is =   " << Charge_misId_E1 << endl;
//            cout << " the charge mis_Id ratio at E2 is =   " << Charge_misId_E2 << endl;
//            cout << " the charge mis_Id ratio at e3 is =   " << Charge_misId_E3 << endl;
            
        }
        
        for (unsigned iCh = 0; iCh < Charge_misId_BVec.size(); iCh++)
        {
          //  cout << " iCh =  " << iCh << "  And the charge mis_Id ratio =  " << Charge_misId_BVec[iCh] <<endl;
            //if(Charge_misId_BVec[iCh]>0)
          //  histo1D["h_Charge_misId_Pt_barrel_Vec"]->Fill(Charge_misId_BVec[iCh],scaleFactor*lumiWeight);
            histo1D["h_Charge_misId_Pt_barrel_Vec"]->Fill(Charge_misId_BVec[iCh]);
           // if(Charge_misId_BVec[iCh]>0)
            histo2D["h_Charge_misId_Pt_barrel"]->Fill(iCh,Charge_misId_BVec[iCh]);
        }
        for (unsigned iCh = 0; iCh < Charge_misId_EVec.size(); iCh++)
        {
          //  cout << " iCh =  " << iCh << "  And the charge mis_Id ratio =  " << Charge_misId_BVec[iCh] <<endl;
           // if(Charge_misId_EVec[iCh]>0)
           // histo1D["h_Charge_misId_Pt_barrel_Vec"]->Fill(Charge_misId_EVec[iCh],scaleFactor*lumiWeight);
            histo1D["h_Charge_misId_Pt_barrel_Vec"]->Fill(Charge_misId_EVec[iCh]);
          //  if(Charge_misId_EVec[iCh]>0)
          //  histo2D["h_Charge_misId_Pt_barrel"]->Fill(iCh,Charge_misId_EVec[iCh],scaleFactor*lumiWeight);
            histo2D["h_Charge_misId_Pt_barrel"]->Fill(iCh,Charge_misId_EVec[iCh]);
        }

        
        
        sumW = (int) sumWeights;
        nEv = (int) nEvents;
        Count_cut[0] = nbEvents_0;
        Count_cut[1] = nbEvents_1PU;
        Count_cut[2] = nbEvents_2BTag;
        Count_cut[3] = nbEvents_3Trig;
        Count_cut[4] = nbEvents_4GPV;
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
        cout << "**The nb of events before any selections & after Applying bTag SF (nbEvents_2BTag) =  " << nbEvents_2BTag << endl;
        cout << "**The nb of events before any selections & after triggers (nbEvents_3Trig) =  " << nbEvents_3Trig << endl;
        cout << "**The nb of events before any selections & after Good Primary Vertex (nbEvents_4GPV) =  " << nbEvents_4GPV << endl;
        cout << "**The nb of events before any selections & after Applying Lep SF (nbEvents_5LepSF) =  " << nbEvents_5LepSF << endl;
        cout << "**The nb of events after 2 leptons (nbEvents_6) =  " << nbEvents_6 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets  (nbEvents_7) =  " << nbEvents_7 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged  (nbEvents_8) =  " << nbEvents_8 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto (nbEvents_9) =  " << nbEvents_9 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto 2SSL Channel (nbEvents_10) =  " << nbEvents_10 << endl;
        cout << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto 2OSL Channel (nbEvents_11) =  " << nbEvents_11 << endl;
        cout << "*****************************************************************************************************" <<endl;
//        cout << "the number of jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_Jets <<endl;
//        cout << "the number of CSVLbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVLbJets <<endl;
//        cout << "the number of CSVMbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVMbJets <<endl;
//        cout << "the number of CSVTbtagged jets after 2 leptons cut (nb_2l_Jets) =  "<< nb_2l_CSVTbJets <<endl;
        
        globalTree->Fill();
        
        
        infoFile << "**The nb of events before any selections or applying SF (nbEvents_0) =   " << nbEvents_0 << endl;
        infoFile << "**The nb of events before any selections & after Applying PU SF (nbEvents_1PU) =  " << nbEvents_1PU << endl;
        infoFile << "**The nb of events before any selections & after Applying bTag SF (nbEvents_2BTag) =  " << nbEvents_2BTag << endl;
        infoFile << "**The nb of events before any selections & after triggers (nbEvents_3Trig) =  " << nbEvents_3Trig << endl;
        infoFile << "**The nb of events before any selections & after Good Primary Vertex (nbEvents_4GPV) =  " << nbEvents_4GPV << endl;
        infoFile << "**The nb of events before any selections & after Applying Lep SF (nbEvents_5LepSF) =  " << nbEvents_5LepSF << endl;
        infoFile << "**The nb of events after 2 leptons (nbEvents_6) =  " << nbEvents_6 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets  (nbEvents_7) =  " << nbEvents_7 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged  (nbEvents_8) =  " << nbEvents_8 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged + Zmass Veto (nbEvents_9) =  " << nbEvents_9 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged +Zmass Veto 2SSL Channel (nbEvents_10) =  " << nbEvents_10 << endl;
        infoFile << "**The nb of events after 2 Lep + >= 3 Jets + >= 1 CSVLB Tagged +Zmass Veto 2OSL Channel (nbEvents_11) =  " << nbEvents_11 << endl;
        
        //////////////////////
        ///  END OF EVENT  ///
        //////////////////////
        
        cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
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
//
//electroncharge.clear();
//// cout << " applying Charge misId" << endl;
////  cout << " the electron charge after electroncharge.clear() " << electroncharge.size() << endl;
//if (Elec_Elec && ApplyCharge_misID)
//{
//    for(unsigned int iElec=0; iElec<selectedElectrons.size(); iElec++)
//    {
//        electroncharge.push_back(selectedElectrons[iElec]->charge());
//        if (debug) cout << "iElec =  " << iElec << "and the electron charge = " << selectedElectrons[iElec]->charge() << endl;
//            ///now applying charge-misid to electrons
//            
//            //creates a randm number between 0-1.
//            double rdmnr = gRandom->Uniform();
//            if (debug) cout<< "the value of random variable =  " << rdmnr << endl;
//            // to save some CPU time only check when the random number is below 0.25%....
//            if (!isData)
//        {
//            if(rdmnr<=0.0025)
//            {
//                if(fabs(selectedElectrons[iElec]->Eta())<1.479)
//                {
//                    if(rdmnr<0.00017)
//                    {
//                        electroncharge[iElec] = electroncharge[iElec] * -1;
//                        if (debug)cout<<"charge flipped!"<<endl;
//                    }
//                }
//                else if(fabs(selectedElectrons[iElec]->Eta())<2.1 && fabs(selectedElectrons[iElec]->Eta())>1.479)
//                {
//                    if(rdmnr<0.0021)
//                    {
//                        electroncharge[iElec] = electroncharge[iElec] * -1;
//                        if (debug)cout<<"charge flipped!"<<endl;
//                    }
//                }
//            }
//            
//        }
//            else
//        {
//            if(rdmnr<=0.0025)
//            {
//                if(fabs(selectedElectrons[iElec]->Eta())<1.479)
//                {
//                    if(rdmnr<0.0002)
//                    {
//                        electroncharge[iElec] = electroncharge[iElec] * -1;
//                        if (debug)cout<<"charge flipped!"<<endl;
//                    }
//                }
//                else if(fabs(selectedElectrons[iElec]->Eta())<2.1 && fabs(selectedElectrons[iElec]->Eta())>1.479)
//                {
//                    if(rdmnr<0.0023)
//                    {
//                        electroncharge[iElec] = electroncharge[iElec] * -1;
//                        if (debug)cout<<"charge flipped!"<<endl;
//                    }
//                }
//            }
//            
//        }
//            
//            }
//            
//            }
//            

