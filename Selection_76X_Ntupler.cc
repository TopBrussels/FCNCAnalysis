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



int main (int argc, char *argv[])
{
    clock_t start = clock();
    
    /// Some variables from POG/PAG
    const float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
    const float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
    const float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.
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
    ApplyElec_SF = false;
    ApplyMu_SF = true;
    ApplyPU_SF = true;
    Apply_btag_SF = false;
    Apply_JetCleaning = true;
    trigged = false;
    debug = false;
    
    std::string channelpostfix = "";
    /////////////////////
    ///  Configuration
    /////////////////////
    /// xml file
    string xmlFileName ="config/Run2SameSignDiLepton_76XSamples.xml";
    if (argc > 1) xmlFileName = (string)argv[1];
    const char *xmlfile = xmlFileName.c_str();
    
    //Setting Lepton Channels
    if(Elec_Elec)
    {
        cout << " --> Using the Electron-Electron channel..." << endl;
        channelpostfix = "_ElEl_";
    }
    else if(Mu_Mu)
    {
        cout << " --> Using the Muon-Muon channel..." << endl;
        channelpostfix = "_MuMu_";
    }
    else if(Elec_Mu)
    {
        cout << " --> Using the Electron-Muon channel..." << endl;
        channelpostfix = "_ElMu_";
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
    
    
    
    ////////////////////////////////////
    ///  DETERMINE EVENT SCALEFACTOR  ///
    /////////////////////////////////////
    
    
    ///// --- Define Scaling Factors ----- /////
  
    string pathToCaliDir = "../TopTreeAnalysisBase/Calibrations/";
    //////PU SF
    LumiReWeighting LumiWeights(pathToCaliDir+"PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root",pathToCaliDir+"PileUpReweighting/pileup_2015Data74X_25ns-Run254231-258750Cert/nominal.root","pileup60","pileup");
    
    /////Muon SF
     string muon_SF_File= "Muon_SF_TopEA.root";
    MuonSFWeight *muonSFWeight = new MuonSFWeight (pathToCaliDir+"LeptonSF/"+muon_SF_File,"SF_totErr", false, false); // (... , ... , debug, print warning)
    
//    MuonSFWeight *muonSFWeightID_T = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);
//    MuonSFWeight *muonSFWeightIso_TT = new MuonSFWeight(pathToCaliDir+"LeptonSF/"+"MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", false, false);  // Tight RelIso, Tight ID
//    
    
    
    /////Electron SF
//    string electron_SF_File= "Elec_SF_TopEA.root";
//    ElectronSFWeight *electronSFWeight = new ElectronSFWeight (pathToCaliDir+"LeptonSF/"+electron_SF_File,"GlobalSF", false, false); // (... , ... , debug, print warning)
//    
    //
   
    
    
    ///-- dfeine Triggering ---//
    // /// HLT Triggers will used in the analysis is according to Top Trigger (Run2)
    // //// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII

    
    //Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep);
    Trigger* DilepTrigger = new Trigger(1,0,0,1);

    
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
    
    float oldLuminosity = 2612.180735004 ;  //anaEnv.Luminosity;  // in 1/pb
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
    
    cout << "lumi is " << Luminosity << endl;
    
    stringstream ss;
    ss << JobNum;
    string strJobNum = ss.str();
    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    int nofSelectedEvents = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
    ////////////////////////////////////
    ///  Loop on datasets
    ////////////////////////////////////
    
    if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
    
    ///////----Start looping over datasets -----/////
    
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
        string previousFilename = "";
        int iFile = -1;
        int isData;
        string dataSetName = datasets[d]->Name();
        if (verbose > 1)
            cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
        if (verbose > 1)
            cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;

        double lumiWeight = -99.;
        if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        {
            lumiWeight=1;
            isData = 1;
        }
        else
        {
            lumiWeight = Luminosity/datasets[d]->EquivalentLumi();
            cout << "the weight to apply for each event of this data set is " << "Lumi / (EquivalentLumi) = "  << Luminosity << " /(" << datasets[d]->EquivalentLumi() << ") = " << Luminosity/datasets[d]->EquivalentLumi()  <<  endl;
           
        }
        
        if (verbose >1) {
            cout << "the lumiweight of the dataset : "<< datasets[d]->Name() << " is " << lumiWeight << endl;
        }
        
        ///////////////////////\\\\\\\\\\\\\\\\\\\\\\\
        ///// Create root file contains histograms \\\\
        /////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\
        
        string runDate = "Test76X";
        string runChannel = channelpostfix;
        string rootPath = "OutPutHistos";
        rootPath += runChannel+runDate+"/";
        mkdir(rootPath.c_str(),0777); // create the directory rootPath if it is not exist already
        string rootPathSample = channelpostfix+datasets[d]->Name();
        mkdir((rootPath+rootPathSample+"/").c_str(),0777);
        string rootFileName = rootPath+rootPathSample+"/"+channelpostfix+datasets[d]->Name()+"_"+strJobNum+".root";
        TFile *fout = new TFile(rootFileName.c_str(),"RECREATE");
        
        map <string,TH1F*> histo1D;
        std::string titlePlot = "";
        
        titlePlot = "cutFlow"+channelpostfix;
        histo1D["h_cutFlow"] = new TH1F(titlePlot.c_str(), "cutflow", 16,-0.5,15.5);
        
        //*** histos for PV *** //
//        titlePlot = "initial_Nb_PV"+channelpostfix;
//        histo1D["h_initial_Nb_PV"] = new TH1F(titlePlot.c_str(), "Initial nb. of PV",  16, - 0.5, 15.5 );
//        titlePlot = "trigged_Nb_PV"+channelpostfix;
//        histo1D["h_trigged_Nb_PV"] = new TH1F(titlePlot.c_str(), "trigged nb. of PV",  16, - 0.5, 15.5 );
//        titlePlot = "Nb_PV_AfterPU_Weight"+channelpostfix;
//        histo1D["h_Nb_PV_AfterPU_Weight"] = new TH1F(titlePlot.c_str(), "after PU Weight nb. of PV",  16, - 0.5, 15.5 );
        

        
        
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

        
        //*** histos for Electrons *** //
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
        
        ////-- histograms for some variables to be used for control regions
        
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

        
//        titlePlot = "WJets_Mass"+channelpostfix;
//        histo1D["h_WJets_Mass"] = new TH1F(titlePlot.c_str(), "After 2SSL + >= 2J + >=1bCSVLBJet cut:  WJets_Mass",  100, 0, 300 );
        
//        map <string,TH2F*> histo2D;
//        titlePlot = "Nb_2L_bjets_vs_Nb_jets"+channelpostfix;
//        histo2D["h_2L_Nb_bjets_vs_Nb_jets"] = new TH2F(titlePlot.c_str(),"After 2L #CSVLbjet Vs #Jets",15,-0.5,14.5, 15, -0.5,14.5);
//        
//        titlePlot = "Nb_2L_2J_bjets_vs_Nb_jets"+channelpostfix;
//        histo2D["h_2L_2J_Nb_bjets_vs_Nb_jets"] = new TH2F(titlePlot.c_str(),"After 2L+>=2Jets #CSVLbjet Vs #Jets",15,-0.5,14.5, 15, -0.5,14.5);
//        
//        titlePlot = "Nb_2L_2J_1b_bjets_vs_Nb_jets"+channelpostfix;
//        histo2D["h_2L_2J_1b_Nb_bjets_vs_Nb_jets"] = new TH2F(titlePlot.c_str(),"After 2L+>=2Jets+1b #CSVLbjet Vs #Jets",15,-0.5,14.5, 15, -0.5,14.5);
//        
//        titlePlot = "Nb_2SSL_2J_1b_bjets_vs_Nb_jets"+channelpostfix;
//        histo2D["h_2SSL_2J_1b_Nb_bjets_vs_Nb_jets"] = new TH2F(titlePlot.c_str(),"After 2SSL+>=2Jets+1b #CSVLbjet Vs #Jets",15,-0.5,14.5, 15, -0.5,14.5);
//        
//        titlePlot = "Nb_2OSL_2J_1b_bjets_vs_Nb_jets"+channelpostfix;
//        histo2D["h_2OSL_2J_1b_Nb_bjets_vs_Nb_jets"] = new TH2F(titlePlot.c_str(),"After 2OSL+>=2Jets+1b #CSVLbjet Vs #Jets",15,-0.5,14.5, 15, -0.5,14.5);
//        
//        titlePlot = "2L_Lep0Pt_vs_Lep1Pt"+channelpostfix;
//        histo2D["h_2L_Lep0Pt_vs_Lep1Pt"] = new TH2F(titlePlot.c_str(),"After 2L Lep0 P_{T} Vs Lep1 P_{T}",100,0,200, 100, 0,200);
//        
        
        
        //// ***************** /////
        /// output root trees /////
        /// ***************** ////
        string rootTreespath = "Ntuples";
        rootTreespath+=runChannel+runDate+"/";
        mkdir(rootTreespath.c_str(),0777);
        string roottreename = rootTreespath+"";
        roottreename+= datasets[d]->Name();
        roottreename += channelpostfix;
        roottreename+= strJobNum;
        roottreename+="_tree.root";
        
        cout << "  - Recreate outputfile ... " << roottreename.c_str() << endl;
        
        // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
       // fileout->cd();
        
        //////////////////////////////
        // My tree - variables //
        //////////////////////////////
        // various weights
        Double_t pu_weight;
        
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;
        
//        double metPt;
//        double metPx;
//        double metPy;
        int nElectrons;
        Double_t sf_electron[10];
        
        //////vectors
        std::vector<double> *nb_Jets;
        std::vector<double> *nb_electrons;
        std::vector<double> *nb_Muons;
        std::vector<double> *nb_CSVLbJets;
        std::vector<double> *nb_CSVMbJets;
        std::vector<double> *nb_CSVTbJets;
        
        std::vector<double> *ptMuon;
        std::vector<double> *pxMuon;
        std::vector<double> *pyMuon;
        std::vector<double> *pzMuon;
        std::vector<double> *etaMuon;
        std::vector<double> *eMuon;
        std::vector<double> *qMuon;
        
        std::vector<double> *ptElectron;
        std::vector<double> *pxElectron;
        std::vector<double> *pyElectron;
        std::vector<double> *pzElectron;
        std::vector<double> *etaElectron;
        std::vector<double> *eElectron;
        std::vector<double> *qElectron;
        
        std::vector<double> *ptJet;
        std::vector<double> *pxJet;
        std::vector<double> *pyJet;
        std::vector<double> *pzJet;
        std::vector<double> *eJet;
        std::vector<double> *etaJet;
        std::vector<double> *qJet;
        std::vector<double> *BtagCSVjet;
        std::vector<bool> *BtagCSVL;
        
        
        std::vector<double> *ptCSVLJet;
        std::vector<double> *pxCSVLJet;
        std::vector<double> *pyCSVLJet;
        std::vector<double> *pzCSVLJet;
        std::vector<double> *eCSVLJet;
        std::vector<double> *etaCSVLJet;
        std::vector<double> *qCSVLJet;
        
        std::vector<double> *ptCSVMJet;
        std::vector<double> *pxCSVMJet;
        std::vector<double> *pyCSVMJet;
        std::vector<double> *pzCSVMJet;
        std::vector<double> *eCSVMJet;
        std::vector<double> *etaCSVMJet;
        std::vector<double> *qCSVMJet;
        
        std::vector<double> *ptCSVTJet;
        std::vector<double> *pxCSVTJet;
        std::vector<double> *pyCSVTJet;
        std::vector<double> *pzCSVTJet;
        std::vector<double> *eCSVTJet;
        std::vector<double> *etaCSVTJet;
        std::vector<double> *qCSVTJet;
        
        ///////////////////////////////
        // My trees                  //
        ///////////////////////////////
        TTree *bookkeeping = new TTree("startevents","startevents");
        bookkeeping->Branch("run_num",&run_num,"run_num/I");
        bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
        bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
        bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
        bookkeeping->Branch("npu",&npu,"npu/I");
        
        // define the output trees may I need to make many trees depend on selection cuts
        /////(Integer variables)
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isData",&isData,"isData/I");
        myTree->Branch("run_num",&run_num,"run_num/I");
        myTree->Branch("evt_num",&evt_num,"evt_num/I");
        myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        myTree->Branch("nvtx",&nvtx,"nvtx/I");
        myTree->Branch("npu",&npu,"npu/I");
        
       // TTree* InitialTree = new TTree("Initialtree","Initialtree");
        
        // Set branches for vectors
        
        
        myTree->Branch("nb_Jets","std::vector<double>",&nb_Jets);
        myTree->Branch("nb_electrons","std::vector<double>",&nb_electrons);
        myTree->Branch("nb_Muons","std::vector<double>",&nb_Muons);
        myTree->Branch("nb_CSVLbJets","std::vector<double>",&nb_CSVLbJets);
        myTree->Branch("nb_CSVMbJets","std::vector<double>",&nb_CSVMbJets);
        myTree->Branch("nb_CSVTbJets","std::vector<double>",&nb_CSVTbJets);
        
        myTree->Branch("ptMuon","std::vector<double>",&ptMuon);   // Make a branch with name ptMuon, from type vector(double), on loaction ptMuon
        myTree->Branch("pxMuon","std::vector<double>",&pxMuon);
        myTree->Branch("pyMuon","std::vector<double>",&pyMuon);
        myTree->Branch("pzMuon","std::vector<double>",&pzMuon);
        myTree->Branch("eMuon","std::vector<double>",&eMuon);
        myTree->Branch("etaMuon","std::vector<double>",&etaMuon);
        myTree->Branch("qMuon","std::vector<double>",&qMuon);
        
        myTree->Branch("ptElectron","std::vector<double>",&ptElectron);   // Make a branch with name ptElectron, from type vector(double), on loaction ptElectron
        myTree->Branch("pxElectron","std::vector<double>",&pxElectron);
        myTree->Branch("pyElectron","std::vector<double>",&pyElectron);
        myTree->Branch("pzElectron","std::vector<double>",&pzElectron);
        myTree->Branch("eElectron","std::vector<double>",&eElectron);
        myTree->Branch("etaElectron","std::vector<double>",&etaElectron);
        myTree->Branch("qElectron","std::vector<double>",&qElectron);
        
        myTree->Branch("ptJet","std::vector<double>",&ptJet);
        myTree->Branch("pxJet","std::vector<double>",&pxJet);
        myTree->Branch("pyJet","std::vector<double>",&pyJet);
        myTree->Branch("pzJet","std::vector<double>",&pzJet);
        myTree->Branch("eJet","std::vector<double>",&eJet);
        myTree->Branch("etaJet","std::vector<double>",&etaJet);
        myTree->Branch("qJet","std::vector<double>",&qJet);
        myTree->Branch("BtagCSVjet", "std::vector<double>",&BtagCSVjet);
        myTree->Branch("BtagCSVL","std::vector<bool>",&BtagCSVL);
        
        
        myTree->Branch("ptCSVLJet","std::vector<double>",&ptCSVLJet);
        myTree->Branch("pxCSVLJet","std::vector<double>",&pxCSVLJet);
        myTree->Branch("pyCSVLJet","std::vector<double>",&pyCSVLJet);
        myTree->Branch("pzCSVLJet","std::vector<double>",&pzCSVLJet);
        myTree->Branch("eCSVLJet","std::vector<double>",&eCSVLJet);
        myTree->Branch("etaCSVLJet","std::vector<double>",&etaCSVLJet);
        myTree->Branch("qCSVLJet","std::vector<double>",&qCSVLJet);
        
        
        myTree->Branch("ptCSVMJet","std::vector<double>",&ptCSVMJet);
        myTree->Branch("pxCSVMJet","std::vector<double>",&pxCSVMJet);
        myTree->Branch("pyCSVMJet","std::vector<double>",&pyCSVMJet);
        myTree->Branch("pzCSVMJet","std::vector<double>",&pzCSVMJet);
        myTree->Branch("eCSVMJet","std::vector<double>",&eCSVMJet);
        myTree->Branch("etaCSVMJet","std::vector<double>",&etaCSVMJet);
        myTree->Branch("qCSVMJet","std::vector<double>",&qCSVMJet);
        
        
        myTree->Branch("ptCSVTJet","std::vector<double>",&ptCSVTJet);
        myTree->Branch("pxCSVTJet","std::vector<double>",&pxCSVTJet);
        myTree->Branch("pyCSVTJet","std::vector<double>",&pyCSVTJet);
        myTree->Branch("pzCSVTJet","std::vector<double>",&pzCSVTJet);
        myTree->Branch("eCSVTJet","std::vector<double>",&eCSVTJet);
        myTree->Branch("etaCSVTJet","std::vector<double>",&etaCSVTJet);
        myTree->Branch("qCSVTJet","std::vector<double>",&qCSVTJet);
        
//        myTree -> Branch("metPt", &metPt, "metPt/D");
//        myTree -> Branch("metPx", &metPx, "metPx/D");
//        myTree -> Branch("metPy", &metPy, "metPy/D");
        
        //open files and load
        cout << "LoadEvent" << endl;
        treeLoader.LoadDataset(datasets[d], anaEnv);

        ////////////////////////////////////
        ///  Loop on events
        ////////////////////////////////////
        //some bookkeeping variables
        nEvents[d] = 0;
        int previousRun = -1;
        int currentRun;
        bool Btagged = false;
        
        int itrigger = 0;
        
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
        vector<TRootPFJet*> selectedJets;
        vector < TRootMuon* > selectedMuons;
        vector < TRootElectron* > selectedElectrons;
        vector<TRootJet*> selectedBCSVLJets;
        vector<TRootJet*> selectedBCSVMJets;
        vector<TRootJet*> selectedBCSVTJets;
        vector<TRootJet*> JetsExcludingHighestCSVLb;
        
        if (verbose > 1) cout << "	Loop over events " << endl;
        /////// ************* ///////////
        /// start looping on events ///
        ////////************** /////////
        //for (unsigned int ievt = 0; ievt < 1000; ievt++) // used for test code to run over few numbers of events
        
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
        {
            nEvents[d]++;
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
            JetsExcludingHighestCSVLb.clear();
            
            if (ievt%1000 == 0)
                std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
            ////////////////////
            ///  LOAD EVENT  ///
            ///////////////////
            
            //TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets, mets);
            TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets, debug); //updated by adding 5th argument to check analysis at 3 Feb
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if (previousFilename != currentFilename)
            {
                previousFilename = currentFilename;
                iFile++;
                cout << "File changed!!! => iFile = " << iFile << endl;
            }
            
            currentRun = event->runId();
            
            if (verbose>2) {
                cout << "The current Run ID  =  "<< currentRun << endl;
            }
            
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
//            if(run_num > 10000){//data
//                isData=1;}
            bookkeeping->Fill();

            
           //// Trigger.checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
            if (Apply_HLT_Triggers)
            {
                DilepTrigger->checkAvail(currentRun, datasets, d, &treeLoader, event, printTrigger);
                itrigger = DilepTrigger->checkIfFired();
            }
           
            
            /////////////////////////
            ///  EVENT SELECTION  ///
            /////////////////////////
            
            //Declare selection instance
            
            Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
            
            
            //// choose good primary vertex
            bool isGoodPV = selection.isPVSelected(vertex, 4 , 24. ,2.); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
            
            // Jets Selection
            selectedJets = selection.GetSelectedJets(30., 2.4, true , "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
            
            
            /// --- Muons Selection -- ///
            selectedMuons = selection.GetSelectedMuons(12. , 2.1 , .15 ,"Tight","Spring15");  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
            
            //// --- Electron Selection --- ///
            selectedElectrons = selection.GetSelectedElectrons(15., 2.4 , "Tight" , "Spring15_25ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
            

            //// sorting objects in the the event according to Pt
            
            sort(selectedJets.begin(), selectedJets.end(),HighestPt());
            sort(selectedMuons.begin(), selectedMuons.end(), HighestPt());
            sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt());
            
            ///// --- Jet Cleaning -- ////
            
            
            if(Apply_JetCleaning)
            {
                if(debug) cout << "Applying jet cleaning " << endl;
                int OrigSize = selectedJets.size();
                for (int origJets=0; origJets<selectedJets.size(); origJets++)
                {
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
            
            
            /*for (unsigned int ijet = 0; ijet < selectedJets.size() ; ijet++)
            {
                if (selectedBCSVLJets.size() != 0)
                {
                    if (fabs(selectedJets[ijet]->Pt() - selectedBCSVLJets[0]->Pt()) > .000000001)
                    {
                        JetsExcludingHighestCSVLb.push_back(selectedJets[ijet]);
                        
                    }
                    
                }
            }*/
            
            
            //// met parameters
//            double met_px = mets[0]->Px();
//            double met_py = mets[0]->Py();
//            double met_pt = sqrt(met_px*met_px + met_py*met_py);
            
            
            
            ////===========================///////
            //// Applying Scale Factors    //////
            ////==========================//////
            
            // scale factor for the event
            Double_t scaleFactor = 1.;
            Double_t Elec_scaleFactor = 1.;
            Double_t Muon_scaleFactor = 1.;
            Double_t PU_scaleFactor =1.;
            Double_t bTag_scaleFactor =1.;
            float muon1SF, muon2SF,electron1SF, electron2SF, muon3SF;
            muon1SF =  muon2SF = electron1SF = electron2SF = muon3SF = 0.;
            
            
            histo1D["h_cutFlow"]->Fill(0., scaleFactor*lumiWeight); //// fill histogram before applying any scaling factors or triggers
            
            
            //           //// PU SF
                        if (ApplyPU_SF)
                        {
                            double puWeight = LumiWeights.ITweight((int)event->nTruePU()); // simplest reweighting, just use reconstructed number of PV. faco
                            if(!isData)
                            {
                                PU_scaleFactor =puWeight;
                                scaleFactor *= PU_scaleFactor;
                            }
                            else if (isData) PU_scaleFactor =1;
                            if (verbose>3) cout << "PU_scaleFactor is " << PU_scaleFactor << endl;
                           
                        }
             histo1D["h_cutFlow"]->Fill(1., scaleFactor*lumiWeight); // Applying PUWeight
            
            /// Lepton SF
            if (ApplyMu_SF )
            {
                if (!isData)
                {
                    if (selectedMuons.size()>0) {muon1SF = muonSFWeight->at(selectedMuons[0]->Eta(),selectedMuons[0]->Pt(),0); Muon_scaleFactor = muon1SF;
                        // if (verbose>1) cout << " the leading muon with Pt " << selectedMuons[0]->Pt() <<" and Eta  =  " << selectedMuons[0]->Eta() << "the Mu_Scaling Factor is " <<  muon1SF << endl;
                    }
                    
                    if (selectedMuons.size()>1) {muon2SF = muonSFWeight->at(selectedMuons[1]->Eta(),selectedMuons[1]->Pt(),0); Muon_scaleFactor *= muon2SF;
                        // if (verbose>1) cout << " the 2nd leading muon with Pt " << selectedMuons[1]->Pt() <<" and Eta  =  " << selectedMuons[1]->Eta() << "the Mu_Scaling Factor is " <<  muon2SF << endl;
                    }
                    
                    //if (selectedMuons.size()>2) {muon2SF = muonSFWeight->at(selectedMuons[2]->Eta(),selectedMuons[2]->Pt(),0); Muon_scaleFactor *= muon3SF;}
                }else if (isData)
                {Muon_scaleFactor =1.;}
                scaleFactor *=Muon_scaleFactor;
                
            }
            histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
            
            
            
//                        if (ApplyElec_SF)
//                        {
//                            if (!isData)
//                            {
//                               // if (selectedElectrons.size() >0 ) {electron1SF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0); Elec_scaleFactor = electron1SF;}
//                               // if (selectedElectrons.size() >1 ) {electron2SF = electronSFWeight->at(selectedElectrons[1]->Eta(),selectedElectrons[1]->Pt(),0); Elec_scaleFactor *= electron2SF;}
//                            }else if (isData)
//                            {Elec_scaleFactor =1.;}
//scaleFactor *= Elec_scaleFactor;
//                        }
//            
//


            
            ////initial histograms after applying scaling factors or triggers
            //histo1D["h_cutFlow"]->Fill(2., scaleFactor*lumiWeight);
            histo1D["h_initial_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
            //histo1D["h_initial_Nb_PV"]->Fill(vertex.size(),scaleFactor*lumiWeight);
            //histo1D["h_initial_Nb_Elec"]->Fill(selectedElectrons.size(),scaleFactor*lumiWeight);
            histo1D["h_initial_Nb_Mu"]->Fill(selectedMuons.size(),scaleFactor*lumiWeight);
            histo1D["h_initial_Nb_CSVLBJets"]->Fill(selectedBCSVLJets.size(),scaleFactor*lumiWeight);
            // histo1D["h_initial_met"]->Fill(met_pt, scaleFactor*lumiWeight);

            
            /////////////////////////////////
            //// *** Apply selections *** ///
            ////////////////////////////////
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
            float mjj=0. ;
            float massDiff_from_Wmass = 0.;
            float MinimassDiff_from_Wmass = 999;
            float WPairMass,PairMass;
            double qLepton0 , qLepton1, qElec0, qElec1, qMu0, qMu1;
            double sum_jet_PT = 0.;
            double Ht_AllJets = 0.;
            double St_AllJets_Leps = 0.;
            double Sum_Leptons_Pt =0.;
            double Zmass_Window = 0.;
            

            /// Trigger
            if(Apply_HLT_Triggers) trigged = itrigger; //trigged = treeLoader.EventTrigged(itrigger);
            //else trigged = true;
           // cout << "trigged = " << trigged << endl;
            if (!trigged) {
                continue;
            }
            if (trigged)
            {
                histo1D["h_cutFlow"]->Fill(3., scaleFactor*lumiWeight);
                histo1D["h_trigged_Nb_Jets"]->Fill(selectedJets.size(),scaleFactor*lumiWeight);
                //histo1D["h_trigged_Nb_PV"]->Fill(vertex.size(),scaleFactor*lumiWeight);
                //histo1D["h_initial_Nb_Elec"]->Fill(selectedElectrons.size(),scaleFactor*lumiWeight);
                histo1D["h_trigged_Nb_Mu"]->Fill(selectedMuons.size(),scaleFactor*lumiWeight);
               // if (verbose>1) cout << "the nb. of initial muons after pass triggers is  =  " << init_muons.size() << endl;
               // if (verbose>1) cout << "the nb. of muons after pass triggers is  =  " << selectedMuons.size() << endl;
                
                //histo1D["h_Nb_PV_AfterPU_Weight"]->Fill(vertex.size(),scaleFactor*lumiWeight);
                
                if (!isGoodPV) {
                    continue;
                }
                if (isGoodPV)
                {
                    if(verbose>3) cout << "GoodPV" << endl;
                    histo1D["h_cutFlow"]->Fill(4., scaleFactor*lumiWeight);
                    //eventSelected = true;
                    int numLep = selectedElectrons.size()+ selectedMuons.size();
                    if (numLep==2)
                    {
                       histo1D["h_cutFlow"]->Fill(5., scaleFactor*lumiWeight);
                        
                        if (Elec_Elec && !Mu_Mu && !Elec_Mu) // Electron-Electron Channel
                        {
                            if (selectedElectrons.size()==2 && selectedMuons.size()==0 &&selectedElectrons[0]->Pt()>26. && selectedElectrons[1]->Pt()>= 15.)
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
                                }
                                
                            }
                            
                        }
                        
                        else if (Mu_Mu && !Elec_Elec && !Elec_Mu) // Muon-Muon channel
                        {
                            
                            if (selectedMuons.size()==2 && selectedElectrons.size()==0)
                            {
                                histo1D["h_cutFlow"]->Fill(6., scaleFactor*lumiWeight);
                                if (selectedMuons[0]->Pt() >= 20. && selectedMuons[1]->Pt()>= 12.)
                                {
                                   
                                    histo1D["h_cutFlow"]->Fill(7., scaleFactor*lumiWeight);
                                    tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                                    // qLepton0 = selectedMuons[0]->charge();
                                    qMu0 = selectedMuons[0]->charge();
                                    tempLepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
                                    //qLepton1 = selectedMuons[1]->charge();
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
                                        histo1D["h_cutFlow"]->Fill(8., scaleFactor*lumiWeight);
                                    }
                                    
                                    
                                }
                                
                            }
                        }
                        
                        else if (Elec_Mu) // Electron- Muon channel
                        {
                            if (selectedElectrons.size()==1 && selectedMuons.size()==1)
                            {
                                if (selectedMuons[0]->Pt() > selectedElectrons[0]->Pt()&& selectedMuons[0]->Pt()>= 20. && selectedElectrons[0]->Pt()>=15.)
                                {
                                    tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                                    qLepton0 = selectedMuons[0]->charge();
                                    tempLepton_1.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                                    qLepton1 = selectedElectrons[0]->charge();
                                    
                                }
                                if (selectedElectrons[0]->Pt()>selectedMuons[0]->Pt() && selectedElectrons[0]->Pt()>=26. && selectedMuons[0]->Pt())
                                {
                                    tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                                    qLepton0 = selectedElectrons[0]->charge();
                                    tempLepton_1.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                                    qLepton1 = selectedMuons[0]->charge();
                                }
                                InvMass_ll = (tempLepton_0+tempLepton_1).M();
                                if (InvMass_ll >12.){diEMu = true;}
                                
                            }
                            
                        }
                        
                        if (diElectron || diMuon || diEMu)
                        {
                            histo1D["h_cutFlow"]->Fill(9., scaleFactor*lumiWeight);
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
                                histo1D["h_cutFlow"]->Fill(10., scaleFactor*lumiWeight);
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
                                    histo1D["h_cutFlow"]->Fill(11., scaleFactor*lumiWeight);
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
                                    
                                    //// --- selecting W jets pair ---- ///
                                    
                                   /* if (JetsExcludingHighestCSVLb.size() >=2)
                                    {
                                        TLorentzVector temp_W_1st_Jet;
                                        TLorentzVector temp_W_2nd_Jet;
                                        for (int ijet = 0; ijet < JetsExcludingHighestCSVLb.size(); ijet++)
                                        {
                                            temp_W_1st_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[ijet]->Px(),JetsExcludingHighestCSVLb[ijet]->Py(),JetsExcludingHighestCSVLb[ijet]->Pz(),JetsExcludingHighestCSVLb[ijet]->Energy());
                                            for (int kjet = ijet+1 ; kjet < JetsExcludingHighestCSVLb.size(); kjet++)
                                            {
                                                temp_W_2nd_Jet.SetPxPyPzE(JetsExcludingHighestCSVLb[kjet]->Px(),JetsExcludingHighestCSVLb[kjet]->Py(),JetsExcludingHighestCSVLb[kjet]->Pz(),JetsExcludingHighestCSVLb[kjet]->Energy());
                                                mjj = (temp_W_1st_Jet+ temp_W_2nd_Jet).M();
                                                massDiff_from_Wmass = fabs(Wmass - mjj);
                                                if (massDiff_from_Wmass < MinimassDiff_from_Wmass)
                                                {
                                                    PairMass = mjj;
                                                    MinimassDiff_from_Wmass = massDiff_from_Wmass;
                                                }
                                            }
                                        }
                                    }
                                    WPairMass = PairMass;
                                    histo1D["h_WJets_Mass"]->Fill(WPairMass,scaleFactor*lumiWeight);*/
                                    
                                    eventSelected = true;
                                    
//                                    //if (Elec_Elec && selectedElectrons[0]->charge()== selectedElectrons[1]->charge())// in case of Same Sign dilepton
//                                    if (diElectron && qElec0 == qElec1)
//                                    {
//                                        SSdiLepton = true;
//                                        histo1D["h_2SSL_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
//                                        histo1D["h_2SSL_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
//                                        histo1D["h_2SSL_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
//                                        Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
//                                        histo1D["h_2SSL_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
//                                        histo1D["h_2SSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
//                                        histo1D["h_2SSL_DeltaR_b0_Elec0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
//                                    }
//                                    if (diElectron && qElec0 != qElec1) // in case of Opposite Sign dilepton
//                                    {
//                                        OSdiLepton = true;
//                                        histo1D["h_2OSL_Mll_2Elec"]->Fill(InvMass_ll,scaleFactor*lumiWeight);
//                                        histo1D["h_2OSL_DeltaPhi_2Elec"]->Fill(tempLepton_0.DeltaPhi(tempLepton_1),scaleFactor*lumiWeight);
//                                        histo1D["h_2OSL_DeltaR_2Elec"]->Fill(tempLepton_0.DeltaR(tempLepton_1),scaleFactor*lumiWeight);
//                                        Sum_Leptons_Pt = tempLepton_0.Pt()+tempLepton_1.Pt();
//                                        histo1D["h_2OSL_SumPT_2Elec"]->Fill(Sum_Leptons_Pt,scaleFactor*lumiWeight);
//                                        //InvMass_lb = (tempLepton_0+SM_bJet).M();
//                                       // cout << " InvMass_lb =  " << InvMass_lb;
//                                        histo1D["h_2OSL_mElec0b0"]->Fill((tempLepton_0+SM_bJet).M(),scaleFactor*lumiWeight);
//                                        histo1D["h_2OSL_DeltaR_b0_Elec0"]->Fill(tempLepton_0.DeltaR(SM_bJet),scaleFactor*lumiWeight);
//                                        
//                                    }
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
                                        histo1D["h_cutFlow"]->Fill(12., scaleFactor*lumiWeight);
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
                                    }
                                    if (OSdiLepton && !SSdiLepton)
                                    {
                                        histo1D["h_cutFlow"]->Fill(13., scaleFactor*lumiWeight);
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
                                    }
                                    
                                    
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
                    
                }//GoodPV
            }//trigged
            
            
            
            
            
            if (! eventSelected )
            {
                continue;
            }
            if(verbose>3)cout << "filling the tree" << endl;
            nb_Jets= new std::vector<double>;
            nb_electrons = new std::vector<double>;
            nb_Muons = new std::vector<double>;
            nb_CSVLbJets = new std::vector<double>;
            nb_CSVMbJets = new std::vector<double>;
            nb_CSVTbJets= new std::vector<double>;
            
            ptMuon = new std::vector<double>;
            pxMuon = new std::vector<double>;
            pyMuon = new std::vector<double>;
            pzMuon = new std::vector<double>;
            etaMuon = new std::vector<double>;
            eMuon = new std::vector<double>;
            qMuon = new std::vector<double>;
            
            ptElectron = new std::vector<double>;
            pxElectron = new std::vector<double>;
            pyElectron = new std::vector<double>;
            pzElectron = new std::vector<double>;
            etaElectron = new std::vector<double>;
            eElectron = new std::vector<double>;
            qElectron = new std::vector<double>;
            
            
            ptJet = new std::vector<double>;
            pxJet = new std::vector<double>;
            pyJet = new std::vector<double>;
            pzJet = new std::vector<double>;
            eJet = new std::vector<double>;
            etaJet = new std::vector<double>;
            qJet = new std::vector<double>;
            BtagCSVjet = new std::vector<double>;
            BtagCSVL = new std::vector<bool>;
            
            ptCSVLJet = new std::vector<double>;
            pxCSVLJet = new std::vector<double>;
            pyCSVLJet = new std::vector<double>;
            pzCSVLJet = new std::vector<double>;
            eCSVLJet = new std::vector<double>;
            etaCSVLJet = new std::vector<double>;
            qCSVLJet = new std::vector<double>;
            
            ptCSVMJet = new std::vector<double>;
            pxCSVMJet = new std::vector<double>;
            pyCSVMJet = new std::vector<double>;
            pzCSVMJet = new std::vector<double>;
            eCSVMJet = new std::vector<double>;
            etaCSVMJet = new std::vector<double>;
            qCSVMJet = new std::vector<double>;
            
            ptCSVTJet = new std::vector<double>;
            pxCSVTJet = new std::vector<double>;
            pyCSVTJet = new std::vector<double>;
            pzCSVTJet = new std::vector<double>;
            eCSVTJet = new std::vector<double>;
            etaCSVTJet = new std::vector<double>;
            qCSVTJet = new std::vector<double>;
            
//            metPt = met_pt;
//            metPx = met_px;
//            metPy = met_py;
            
            
            Double_t tempJet_Size = 0.;
            Double_t tempElec_Size = 0.;
            Double_t tempMu_Size = 0.;
            Double_t tempCSVLBJet_Size = 0.;
            Double_t tempCSVMBJet_Size = 0.;
            Double_t tempCSVTBJet_Size = 0.;
            
            for (unsigned int i = 0; i < selectedElectrons.size(); i++)
            {
                TRootElectron* tempElectron = (TRootElectron*) selectedElectrons[i];
                ptElectron->push_back(tempElectron->Pt());
                pxElectron->push_back(tempElectron->Px());
                pyElectron->push_back(tempElectron->Py());
                pzElectron->push_back(tempElectron->Pz());
                eElectron->push_back(tempElectron->Energy());
                etaElectron->push_back(tempElectron->Eta());
                qElectron->push_back(tempElectron->charge());
                tempElec_Size++;
            }
            nb_electrons->push_back(tempElec_Size);
            
            
            for (unsigned int i = 0; i < selectedMuons.size(); i++)
            {
                TRootMuon* tempMuon = (TRootMuon*) selectedMuons[i];
                ptMuon->push_back(tempMuon->Pt());
                pxMuon->push_back(tempMuon->Px());
                pyMuon->push_back(tempMuon->Py());
                pzMuon->push_back(tempMuon->Pz());
                eMuon->push_back(tempMuon->Energy());
                etaMuon->push_back(tempMuon->Eta());
                qMuon->push_back(tempMuon->charge());
                tempMu_Size++;
            }
            nb_Muons->push_back(tempMu_Size);
            
            for (unsigned int i =0; i < selectedJets.size(); i ++)
            {
                TRootJet* tempJet = (TRootJet*) selectedJets[i];
                ptJet->push_back(tempJet->Pt());
                pxJet->push_back(tempJet->Px());
                pyJet->push_back(tempJet->Py());
                pzJet->push_back(tempJet->Pz());
                eJet->push_back(tempJet->Energy());
                etaJet->push_back(tempJet->Eta());
                qJet->push_back(tempJet->charge());
                BtagCSVjet->push_back(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags());
                //  BtagCSVL->push_back(BtagBooleans[i]);
                tempJet_Size++;
            }
            nb_Jets->push_back(tempJet_Size);
            
            for (unsigned int i = 0; i< selectedBCSVLJets.size(); i++)
            {
                TRootJet* tempJet = (TRootJet*) selectedBCSVLJets[i];
                ptCSVLJet->push_back(tempJet->Pt());
                pxCSVLJet->push_back(tempJet->Px());
                pyCSVLJet->push_back(tempJet->Py());
                pzCSVLJet->push_back(tempJet->Pz());
                eCSVLJet->push_back(tempJet->Energy());
                etaCSVLJet->push_back(tempJet->Eta());
                qCSVLJet->push_back(tempJet->charge());
                tempCSVLBJet_Size++;
            }
            nb_CSVLbJets->push_back(tempCSVLBJet_Size);
            
            for (unsigned int i = 0; i< selectedBCSVMJets.size(); i++)
            {
                TRootJet* tempJet = (TRootJet*) selectedBCSVMJets[i];
                ptCSVMJet->push_back(tempJet->Pt());
                pxCSVMJet->push_back(tempJet->Px());
                pyCSVMJet->push_back(tempJet->Py());
                pzCSVMJet->push_back(tempJet->Pz());
                eCSVMJet->push_back(tempJet->Energy());
                etaCSVMJet->push_back(tempJet->Eta());
                qCSVMJet->push_back(tempJet->charge());
                tempCSVMBJet_Size++;
            }
            nb_CSVMbJets->push_back(tempCSVMBJet_Size);
            
            for (unsigned int i = 0; i< selectedBCSVTJets.size(); i++)
            {
                TRootJet* tempJet = (TRootJet*) selectedBCSVTJets[i];
                ptCSVTJet->push_back(tempJet->Pt());
                pxCSVTJet->push_back(tempJet->Px());
                pyCSVTJet->push_back(tempJet->Py());
                pzCSVTJet->push_back(tempJet->Pz());
                eCSVTJet->push_back(tempJet->Energy());
                etaCSVTJet->push_back(tempJet->Eta());
                qCSVTJet->push_back(tempJet->charge());
                tempCSVTBJet_Size++;
            }
            nb_CSVTbJets->push_back(tempCSVTBJet_Size);
            
            nofSelectedEvents++;
            myTree->Fill();
            
            delete pxElectron;
            delete pyElectron;
            delete pzElectron;
            delete etaElectron;
            delete eElectron;
            delete qElectron;
            
            delete ptMuon;
            delete pxMuon;
            delete pyMuon;
            delete pzMuon;
            delete etaMuon;
            delete eMuon;
            delete qMuon;
            
            delete ptJet;
            delete pxJet;
            delete pyJet;
            delete pzJet;
            delete eJet;
            delete etaJet;
            delete qJet;
            delete BtagCSVjet;
            delete BtagCSVL;
            delete ptCSVLJet;
            delete pxCSVLJet;
            delete pyCSVLJet;
            delete pzCSVLJet;
            delete eCSVLJet;
            delete etaCSVLJet;
            delete qCSVLJet;
            
            delete ptCSVMJet;
            delete pxCSVMJet;
            delete pyCSVMJet;
            delete pzCSVMJet;
            delete eCSVMJet;
            delete etaCSVMJet;
            delete qCSVMJet;
            
            delete ptCSVTJet;
            delete pxCSVTJet;
            delete pyCSVTJet;
            delete pzCSVTJet;
            delete eCSVTJet;
            delete etaCSVTJet;
            delete qCSVTJet;
        
            if (verbose > 2)
                cout << "  Event " << ievt << " is selected" << endl;
        } // end the loop over events
        
        
        //////////////////////
        ///  END OF EVENT  ///
        //////////////////////
        
        cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
        myTree->Write();
        fileout->cd();
        fileout->Write();
        fileout->Close();
        delete fileout;  //   delete myTree;
        
        ///*****************///
        ///   Write plots   ///
        ///*****************///
        string pathPNG = "OutPutHistos/";
        mkdir(pathPNG.c_str(),0777);
        mkdir((pathPNG+"1DPlot/").c_str(),0777); // 0777 if it doesn't exist already, make it
        fout->cd();
        for (map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
        {
            if(verbose>3) cout << "1D Plot: " << it->first << endl;
            TH1F *temp = it->second;
            string name = it->first;
            temp->Draw();
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
        
        ///////////////////
        /// CLEANING
        /////////////////
        //delete fileout;  //   delete myTree;
        delete fout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }///end the loop over the datasets
    
    
   // delete tcAnaEnv;
    
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;

   
}