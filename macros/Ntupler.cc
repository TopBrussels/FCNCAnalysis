///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntupler.cc: This macro is intended to be an example analysis macro which works out of the box.           /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.             /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "../../TopTreeProducer/interface/TRootRun.h"
#include "../../TopTreeProducer/interface/TRootEvent.h"
#include "../../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../../TopTreeAnalysis/macros/Style.C"
#include "../../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

	

int main (int argc, char *argv[])
{
    
    
    clock_t start = clock();
    
         
    ///////////////////////////////////////////////////////////
    //   different options for executing this macro	     //
    ///////////////////////////////////////////////////////////

    std::string tempxmlName;
    std::string channelName = "undefined";
    bool information = true; 
    bool warnings = true; 
    bool debug = false; 
    int nb_Leptons = 0; 
    bool foundxml= false;

    for(int iarg = 0; iarg < argc && argc>1 ; iarg++)
    {
	    std::string argval=argv[iarg];
    	    
	    if(argval=="--help" || argval =="--h")
    	    {
    		    cout << "--NoWarnings: put warnings off " << endl; 
    		    cout << "--NoInfo: put information off " << endl; 
    		    cout << "--debug: put debug output on" << endl; 
    		    cout << "--xml myxml.xml: change Xml file" << endl;
    		    cout << " " << endl; 
    		    cout << "--1L3B: use the 1 lepton + 3 b-tags channel" << endl; 
    		    cout << "--SSdilepton: use the same sign dilepton channel" << endl;
    		    cout << "--OSdilepton: use the opposite sign dilepton channel" << endl;
    		    cout << "--3L: use the 3 lepton channel (exactly 3)" << endl;
    		    cout << "--45: at least 4 leptons " << endl; 
    		    
    		    return 0;
	    }
	    if (argval=="--xml") {
    		    iarg++;
    		    tempxmlName = argv[iarg];
		    foundxml = true; 
		    if(debug) cout << "foundxml = true" << endl; 
    		    
    	    }
    	    if (argval=="--NoInfo") {
    		    iarg++;
    		    information = false;  
    	    }
    	    if (argval=="--NoWarnings") {
    		    iarg++;
    		    warnings = false;  
    	    }
    	    if (argval=="--debug") {
    		    iarg++;
    		    debug = true;  
    	    }
    	    		    
    	    if (argval=="--1L3B") {
    		    nb_Leptons = 1; 
		    channelName = "1L3B";
    		    if(!foundxml) tempxmlName = "../config/FCNC_1L3B_config.xml";
	    }
    	    if (argval=="--SSdilepton") {
    		    nb_Leptons = 2; 
		    channelName = "SSdilepton";
    		    if(!foundxml) tempxmlName = "../config/FCNC_SSdilepton_config.xml";
	    }
    	    if (argval=="--OSdilepton") {
    		    nb_Leptons = 2; 
		    channelName = "OSdilepton";
    		    if(!foundxml) tempxmlName = "../config/FCNC_OSdilepton_config.xml";
	    }
    	    if (argval=="--3L") {
    		    nb_Leptons = 3;
		    channelName = "3L";
    		    if(!foundxml) tempxmlName = "../config/FCNC_3L_config.xml";
	    }
    	    if (argval=="--45") {
    		    nb_Leptons = 3;
		    channelName = "45";
    		    if(!foundxml) tempxmlName = "../config/FCNC_45_config.xml";
	    }
    	    
    	    

    } 
    //put in a warning 
    if(channelName.find("undefined")!=string::npos && warnings) std::cout << "[WARNING]      No channel was defined" << endl; 
    if(nb_Leptons == 0 && warnings) std::cout << "[WARNING]      No nb of leptons was defined, default setting is 0" << endl; 





    //SetStyle if needed
    //setTDRStyle();
    setMyStyle();
    
    /////////////////////
    // Configuration
    /////////////////////
    
    //xml file
    string xmlFileName = tempxmlName;
    
 
    const char *xmlfile = xmlFileName.c_str();
    const char *channel = channelName.c_str();
    if(information)
    {
    cout << "********************************************************" << endl;
    cout << "used config file: " << xmlfile << endl;
    cout << "used channel: " << channel << endl; 
    cout << "********************************************************" << endl;
    }
    //Configuration output format
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
    ////////////////////////////////////
    /// AnalysisEnvironment
    ////////////////////////////////////
    
    AnalysisEnvironment anaEnv;
    if(debug) std::cout << "Loading the analysisenvironment" << endl; 
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    if(debug) std::cout << "done creating AnalysisEnvironmentLoader" << endl; 

    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    int verbose = anaEnv.Verbose;
    float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
    
    if(debug)   cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
    
    /////////////////////
    // Load Datasets
    /////////////////////
    
    TTreeLoader treeLoader;
    if(debug)    cout << " - Load datasets ..." << endl;
    vector < Dataset* > datasets;
    
    treeLoader.LoadDatasets (datasets, xmlfile);
    for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    
    float Luminosity = oldLuminosity;
    
    cout << "********************************************************" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
        
        string dataSetName = datasets[d]->Name();
        cout << "datasets: " << dataSetName << endl;
    }
    cout << "********************************************************" << endl;
    
    
    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    Double_t *nEvents_Selected = new Double_t[datasets.size()];
    
    ////////////////////////
    // PileUp Reweighting //

    ////////////////////////
    
    //cout << Luminosity << endl;
    
    LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;
    
    LumiWeights = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root","../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
    LumiWeightsUp = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
    LumiWeightsDown = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");
    
   
    if(debug)    cout << " Initialized LumiReWeighting stuff" << endl;
    
    ////////////////////////////////////
    //	Loop on datasets
    ////////////////////////////////////
    
    if(debug) cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
        if (debug)
            std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
//      make root tree file name
        string roottreename = "../ntuples/";
	roottreename+= channelName; 
	roottreename+= "_"; 
	roottreename+= datasets[d]->Name();
        roottreename+="_tree.root";
	//        cout << "creating tree in file " << roottreename << endl;
        
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        /////////////////////
        // My tree         //
        ////////////////////
        Int_t nElectrons;
        Double_t pX_electron[10];
        Double_t pY_electron[10];
        Double_t pZ_electron[10];
        Double_t E_electron[10];
	Double_t pfIso_electron[10];
        Int_t charge_electron[10];
        
        Int_t nMuons;
        Double_t pX_muon[10];
        Double_t pY_muon[10];
        Double_t pZ_muon[10];
	Double_t E_muon[10];
	Double_t pfIso_muon[10];
        Int_t charge_muon[10];
        
        Int_t nJets;
        Double_t pX_jet[10];
        Double_t pY_jet[10];
        Double_t pZ_jet[10];
        Double_t E_jet[10];
	
	Int_t nBJets;
        Double_t pX_Bjet[10];
        Double_t pY_Bjet[10];
        Double_t pZ_Bjet[10];
        Double_t E_Bjet[10];
	
	Int_t nLJets;
        Double_t pX_Ljet[10];
        Double_t pY_Ljet[10];
        Double_t pZ_Ljet[10];
	Double_t Phi_Ljet[10];
	Double_t E_Ljet[10];

        
        Double_t missingEt;
	Double_t missingEt_Phi;
	Double_t missingEt_Theta;
	Double_t missingEt_pX;
	Double_t missingEt_pY;
	Double_t missingEt_pZ;
	

        //45 channel variables
	Double_t InvMass_4lept_Zdecay;
	Double_t InvMass_FCNC_top_Zdecay; 
	Double_t InvMass_SM_lb; 
	Double_t InvMass_SM_W_lv;
	Double_t InvMass_SM_W_qq;
	Double_t InvMass_SM_W;
	Double_t InvMass_SM_top_blv;
	Double_t InvMass_SM_top_bqq;
	Double_t InvMass_SM_top;
	Double_t TrMass_W; 
	Double_t TrMass_W_qq; 
	Double_t TrMass_W_lv; 
	
	Double_t Phi_Higgs; 
	Double_t Eta_Higgs; 
        
	Double_t Bdiscr;
	//3L channel variables
	Double_t InvMass_Z; 
	Double_t InvMass_H;
	Double_t Z_candidate_pT;
	Double_t Z_candidate_Eta;
	Double_t H_candidate_pT;
	Double_t H_candidate_Eta;
	
	Double_t pT_FCNC_top_tcZ;
	Double_t Eta_FCNC_top_tcZ; 
	Double_t InvMass_FCNC_top_tcZ; 
	
	Double_t pT_FCNC_top_candidate;
	Double_t Eta_FCNC_top_candidate; 
	Double_t InvMass_FCNC_top_candidate;
	
	Double_t pT_FCNC_top_tcH_ZZ_llqq;
	Double_t Eta_FCNC_top_tcH_ZZ_llqq; 
	Double_t InvMass_FCNC_top_tcH_ZZ_llqq; 
	
	
	Double_t Bjet_Eta; 
	Double_t Bjet_Phi; 
	Double_t Bjet_Px; 
	Double_t Bjet_Pt;
	Double_t Bjet_Py;
	Double_t Bjet_Pz;
	Double_t FCNC_ll_Eta; 
	Double_t FCNC_ll_Phi; 
	Double_t FCNC_ll_Px; 
	Double_t FCNC_ll_Pt;
	Double_t FCNC_ll_Py;
	Double_t FCNC_ll_Pz;
	
	Double_t InvMass_FCNC_ll; 
	Double_t DeltaR_SMlb_FCNCll; 
	Double_t DeltaPhi_SMlb_FCNCll;
	
	
	
	
	Int_t nEvents_Tree; 
        Int_t isdata;
        // various weights
        Double_t pu_weight;
        
        
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");
	myTree->Branch("nEvents_Tree",&nEvents_Tree,"nEvents_Tree/I");
        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        myTree->Branch("pX_electron",pX_electron,"pX_electron[nElectrons]/D");
        myTree->Branch("pY_electron",pY_electron,"pY_electron[nElectrons]/D");
        myTree->Branch("pZ_electron",pZ_electron,"pZ_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pX_muon",pX_muon,"pX_muon[nMuons]/D");
        myTree->Branch("pY_muon",pY_muon,"pY_muon[nMuons]/D");
        myTree->Branch("pZ_muon",pZ_muon,"pZ_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        
        myTree->Branch("nJets",&nJets, "nJets/I");
        myTree->Branch("pX_jet",pX_jet,"pX_jet[nJets]/D");
        myTree->Branch("pY_jet",pY_jet,"pY_jet[nJets]/D");
        myTree->Branch("pZ_jet",pZ_jet,"pZ_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
	
	myTree->Branch("nBJets",&nJets, "nBJets/I");
        myTree->Branch("pX_Bjet",pX_Bjet,"pX_Bjet[nBJets]/D");
        myTree->Branch("pY_Bjet",pY_Bjet,"pY_Bjet[nBJets]/D");
        myTree->Branch("pZ_Bjet",pZ_Bjet,"pZ_Bjet[nBJets]/D");
        myTree->Branch("E_Bjet",E_Bjet,"E_Bjet[nBJets]/D");
	
	myTree->Branch("nLJets",&nJets, "nLJets/I");
        myTree->Branch("pX_Ljet",pX_Ljet,"pX_Ljet[nLJets]/D");
        myTree->Branch("pY_Ljet",pY_Ljet,"pY_Ljet[nLJets]/D");
        myTree->Branch("pZ_Ljet",pZ_Ljet,"pZ_Ljet[nLJets]/D");
        myTree->Branch("E_Ljet",E_Ljet,"E_Ljet[nLJets]/D");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
	myTree->Branch("missingEt_Phi",&missingEt_Phi,"missingEt_Phi/D");
	myTree->Branch("missingEt_Theta",&missingEt_Theta,"missingEt_Theta/D");
	myTree->Branch("missingEt_pX",&missingEt_pX,"missingEt_pX/D");
	myTree->Branch("missingEt_pY",&missingEt_pY,"missingEt_pY/D");

        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
	
	if(channelName.find("45")!=string::npos)
	{
		//45 channel variables
		myTree->Branch("InvMass_4lept_Zdecay",&InvMass_4lept_Zdecay,"InvMass_4lept_Zdecay/D");
		myTree->Branch("InvMass_FCNC_top_Zdecay",&InvMass_FCNC_top_Zdecay,"InvMass_FCNC_top_Zdecay/D");
		myTree->Branch("InvMass_SM_lb",&InvMass_SM_lb,"InvMass_SM_lb/D");
		myTree->Branch("InvMass_SM_W_lv",&InvMass_SM_W_lv,"InvMass_SM_W_lv/D");
		myTree->Branch("InvMass_SM_W_qq",&InvMass_SM_W_qq,"InvMass_SM_W_qq/D");
		myTree->Branch("InvMass_SM_W",&InvMass_SM_W,"InvMass_SM_W/D");
		myTree->Branch("InvMass_SM_top_blv",&InvMass_SM_top_blv,"InvMass_SM_top_blv/D");
		myTree->Branch("InvMass_SM_top_bqq",&InvMass_SM_top_bqq,"InvMass_SM_top_bqq/D");
		myTree->Branch("InvMass_SM_top",&InvMass_SM_top,"InvMass_SM_top/D");
		myTree->Branch("TrMass_W",&TrMass_W,"TrMass_W/D");
		myTree->Branch("TrMass_W_qq",&TrMass_W_qq,"TrMass_W_qq/D");
		myTree->Branch("TrMass_W_lv",&TrMass_W_lv,"TrMass_W_lv/D");

		myTree->Branch("Phi_Higgs",&Phi_Higgs,"Phi_Higgs/D"); 
		myTree->Branch("Eta_Higgs",&Eta_Higgs,"Eta_Higgs/D");
		
		myTree->Branch("Bdiscr",&Bdiscr,"Bdiscr/D");
	} 
       	if(channelName.find("3L")!=string::npos)
	{
        	//3L channel variables
        	myTree->Branch("InvMass_Z",&InvMass_Z,"InvMass_Z/D"); 
		myTree->Branch("InvMass_H",&InvMass_H,"InvMass_H/D");
		myTree->Branch("Z_candidate_pT",&Z_candidate_pT,"Z_candidate_pT/D");
		myTree->Branch("Z_candidate_Eta",&Z_candidate_Eta,"Z_candidate_Eta/D");
		myTree->Branch("H_candidate_pT",&H_candidate_pT,"H_candidate_pT/D");
		myTree->Branch("H_candidate_Eta",&H_candidate_Eta,"H_candidate_Eta/D");
		myTree->Branch("pT_FCNC_top_tcZ",&pT_FCNC_top_tcZ,"pT_FCNC_top_tcZ/D");
		myTree->Branch("Eta_FCNC_top_tcZ",&Eta_FCNC_top_tcZ,"Eta_FCNC_top_tcZ/D");
		myTree->Branch("InvMass_FCNC_top_tcZ",&InvMass_FCNC_top_tcZ,"InvMass_FCNC_top_tcZ/D");
		myTree->Branch("pT_FCNC_top_candidate",&pT_FCNC_top_candidate,"pT_FCNC_top_candidate/D");
		myTree->Branch("Eta_FCNC_top_candidate",&Eta_FCNC_top_candidate,"Eta_FCNC_top_candidate/D");
		myTree->Branch("InvMass_FCNC_top_candidate",&InvMass_FCNC_top_candidate,"InvMass_FCNC_top_candidate/D");
		myTree->Branch("pT_FCNC_top_tcH_ZZ_llqq",&pT_FCNC_top_tcH_ZZ_llqq,"pT_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("Eta_FCNC_top_tcH_ZZ_llqq",&Eta_FCNC_top_tcH_ZZ_llqq,"Eta_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("InvMass_FCNC_top_tcH_ZZ_llqq",&InvMass_FCNC_top_tcH_ZZ_llqq,"InvMass_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("Bjet_Eta",&Bjet_Eta, "Bjet_Eta/D"); 
		myTree->Branch("Bjet_Phi",&Bjet_Phi, "Bjet_Phi/D"); 
		myTree->Branch("Bjet_Px",&Bjet_Px, "Bjet_Px/D"); 
		myTree->Branch("Bjet_Pt",&Bjet_Pt, "Bjet_Pt/D");
		myTree->Branch("Bjet_Py",&Bjet_Py, "Bjet_Py/D");
		myTree->Branch("Bjet_Pz",&Bjet_Pz, "Bjet_Pz/D");
		myTree->Branch("FCNC_ll_Eta",&FCNC_ll_Eta, "FCNC_ll_Eta/D"); 
		myTree->Branch("FCNC_ll_Phi",&FCNC_ll_Phi, "FCNC_ll_Phi/D"); 
		myTree->Branch("FCNC_ll_Px",&FCNC_ll_Px, "FCNC_ll_Px/D"); 
		myTree->Branch("FCNC_ll_Pt",&FCNC_ll_Pt, "FCNC_ll_Pt/D");
		myTree->Branch("FCNC_ll_Py",&FCNC_ll_Py, "FCNC_ll_Py/D");
		myTree->Branch("FCNC_ll_Pz",&FCNC_ll_Pz, "FCNC_ll_Pz/D");
		myTree->Branch("InvMass_SM_lb",&InvMass_SM_lb, "InvMass_SM_lb/D"); 
		myTree->Branch("InvMass_FCNC_ll",&InvMass_FCNC_ll, "InvMass_FCNC_ll/D"); 
		myTree->Branch("DeltaR_SMlb_FCNCll",&DeltaR_SMlb_FCNCll, "DeltaR_SMlb_FCNCll/D"); 
		myTree->Branch("DeltaPhi_SMlb_FCNCll",&DeltaPhi_SMlb_FCNCll, "DeltaPhi_SMlb_FCNCll/D");
		myTree->Branch("Bdiscr",&Bdiscr,"Bdiscr/D");
	
	}
       
	//        myTree->Print();

        TH1F * EventSummary = new TH1F("EventSummary","EventSummary",2,0,2);
        TH1F * Xsection = new TH1F("Xsection","Xsection",2,0,2);
        
        //open files and load
	if(debug)       cout<<"LoadEvent"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);
	if(debug)       cout<<"LoadEvent"<<endl;
        
        
        
        /////////////////////////////////////
        /// Initialize JEC factors
        /////////////////////////////////////
   	    
        vector<JetCorrectorParameters> vCorrParam;
        
        /*JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
         JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
         JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
         
         //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
         vCorrParam.push_back(*L1JetPar);
         vCorrParam.push_back(*L2JetPar);
         vCorrParam.push_back(*L3JetPar);
         vector<TRootJet*> selectedLightJets
         if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // DATA!
         JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
         vCorrParam.push_back(*ResJetCorPar);
         }*/
        
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
        
        // true means redo also the L1
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
        
        
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
        
        nEvents[d] = 0;
	nEvents_Selected[d] = 0; 
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
	
	// store number of events in ntuple 
        nEvents_Tree = datasets[d]->NofEvtsToRunOver(); 
	EventSummary->SetBinContent(1, datasets[d]->NofEvtsToRunOver());
	
	
	//define all cross sections
	string datasetName = datasets[d]->Name(); 
		
	
	if(datasetName.find("GluGluHiggs4lep")!=string::npos) {Xsection->SetBinContent(1,0.005109893);}
	if(datasetName.find("VBHiggs4lep")!=string::npos) {Xsection->SetBinContent(1,0.00414341777);}
	
	if(datasetName.find("WW_To2L2Nu")!=string::npos) {Xsection->SetBinContent(1,5.757);}
	if(datasetName.find("WZ_To2L2Q")!=string::npos) {Xsection->SetBinContent(1,2.267);}
	if(datasetName.find("WZ_To3LNu")!=string::npos) {Xsection->SetBinContent(1,1.087);}
	if(datasetName.find("ZZ_To2L2Nu")!=string::npos) {Xsection->SetBinContent(1,0.713);}
	if(datasetName.find("ZZ_To2L2Q")!=string::npos) {Xsection->SetBinContent(1,2.492);}
	if(datasetName.find("ZZ_To4L")!=string::npos) {Xsection->SetBinContent(1,0.18);}
	
	if(datasetName.find("TBZ_ToLL_4F")!=string::npos) {Xsection->SetBinContent(1,0.0114);}
	if(datasetName.find("ttH")!=string::npos) {Xsection->SetBinContent(1,0.1293);}
	if(datasetName.find("TTZ")!=string::npos) {Xsection->SetBinContent(1,0.172);}
	if(datasetName.find("TTW")!=string::npos) {Xsection->SetBinContent(1,0.2148);}
		
	if(datasetName.find("ST_T_s-ch")!=string::npos) {Xsection->SetBinContent(1,3.79);}
        if(datasetName.find("ST_TBar_s-ch")!=string::npos) {Xsection->SetBinContent(1,1.76);}
	if(datasetName.find("ST_T_tW-ch")!=string::npos) {Xsection->SetBinContent(1,11.1);}
        if(datasetName.find("ST_TBar_tW-ch")!=string::npos) {Xsection->SetBinContent(1,11.1);}
	
	if(datasetName.find("TT_SemiLeptMGDecays")!=string::npos) {Xsection->SetBinContent(1,110.26;}
	if(datasetName.find("TT_FullLeptMGDecays")!=string::npos) {Xsection->SetBinContent(1,26.42);}
		
	if(datasetName.find("Z_1Jets")!=string::npos) {Xsection->SetBinContent(1,671.83);}
	if(datasetName.find("Z_2Jets")!=string::npos) {Xsection->SetBinContent(1,216.76);}
	if(datasetName.find("Z_3Jets")!=string::npos) {Xsection->SetBinContent(1,61.20);}
	if(datasetName.find("Z_4Jets")!=string::npos) {Xsection->SetBinContent(1,27.59);}
		
	
	if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets")!=string::npos) {Xsection->SetBinContent(1,0.090636);	}
	if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL")!=string::npos) {Xsection->SetBinContent(1,0.022659);	}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL")!=string::npos) {	Xsection->SetBinContent(1, 0.0005135);	}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL")!=string::npos) {Xsection->SetBinContent(1, 0.0018609);}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL")!=string::npos) {	Xsection->SetBinContent(1, 0.00067929);	}		
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL")!=string::npos) {Xsection->SetBinContent(1, 0.00016516);	}
	if(datasetName.find("TTJetsTocZbW")!=string::npos) {Xsection->SetBinContent(1, 0.1575);	}
	
        
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++) // event loop
            //for (unsigned int ievt = 0; ievt < 20000; ievt++)
        {
            
            vector < TRootVertex* > vertex;
            vector < TRootMuon* > init_muons;
            vector < TRootElectron* > init_electrons;
            vector < TRootJet* > init_jets_corrected;
            vector < TRootJet* > init_jets;
            vector < TRootMET* > mets;
            vector < TRootGenJet* > genjets;
            
            nEvents[d]++;
            
            if(ievt%1000 == 0)
                std::cout<<"Processing the "<<ievt<<"th event (" <<
		((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100 <<"%)" << " +> "<< (nEvents_Selected[d]/(double)datasets[d]->NofEvtsToRunOver())*100 << "% selected of the total events" << flush<<"\r";
            
            ////////////////
            // LOAD EVENT //
            ////////////////
            
            TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
            isdata=0;
            if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
                genjets = treeLoader.LoadGenJet(ievt,false);
                sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
            }
            else{
                isdata=1;
            }
            
            
	    
            /////////////////////////////////
            // DETERMINE EVENT SCALEFACTOR //
            /////////////////////////////////
            
            // scale factor for the event
            float scaleFactor = 1.;
            
            // PU reweighting
            
            double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );

            if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
                lumiWeight=1;
            
            // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
            // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
            
            pu_weight=lumiWeight;
            
            scaleFactor = scaleFactor*lumiWeight;
            
            ///////////////////
            // TRIGGER SETUP //
            ///////////////////
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename){
                previousFilename = currentFilename;
                iFile++;
               // cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << endl;
            }
            
            int currentRun = event->runId();
            
            if(previousRun != currentRun)
                previousRun = currentRun;
            
	    //triggering only for data, I only have events that either have at least 2 muons or 2 electrons so use double, 
	    // otherwise include EMU trigger as well
            /*bool triggered = false; 
	    int trigger1; 
	    int trigger2; 
	    if(isdata == 1)
	    {
	    	cout << "isdata with currentRun " << currentRun << " iFile " << iFile << endl; 
		trigger1 = treeLoader.iTrigger("HLT_Mu17_Mu8_v17",currentRun,iFile); //double muon
		cout << "trigger1" << endl; 
		trigger2 = treeLoader.iTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",currentRun,iFile); //double electron
		cout << "trigger2" << endl; 
		
		triggered = (treeLoader.EventTrigged(trigger1) || treeLoader.EventTrigged(trigger2));
		
	 	cout << "triggered is true" << endl; 
		
	    }
            else triggered = true; 
	    
	    
	    if(!triggered) continue; //if the events isn't triggered go to the next event
	    */
	    /////////////////////////////////////////////////////////////////////////////
            // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
            /////////////////////////////////////////////////////////////////////////////
            
            if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
                
                jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
            
            /////////////////////
            // EVENT SELECTION //
            /////////////////////
            
            //Declare selection instance
            Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
/*          selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); // standard TOP jet selection
            selection.setMuonCuts(5,2.5,0.4,0.2,0.3,1,0.5,5,0); // standard mu selection but with looser iso
            selection.setElectronCuts(10,2.5,0.4,0.02,0.5,0.3,0); // standard ele selection but with looser iso
*/            
	    //define selection cuts --> have to be validated!!!
	    // From the class Selection the following functions are used: 
	    //      void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon) 
	    //      void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid)  
	    //      void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
	    //void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits)
		
	    selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1); 
	    selection.setDiMuonCuts(10.,2.5,0.2,0.04);
	    selection.setDiElectronCuts(15.0,2.4,0.15,0.04,0.5,1,0.3,1); 
	    	
	    
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
            
            if(!isGoodPV)
                continue;
            
            missingEt=mets[0]->Pt();
	    missingEt_pX=mets[0]->Px();
	    missingEt_pY=mets[0]->Py();
	    missingEt_Phi = mets[0]->Phi(); 
	    missingEt_Theta = mets[0]->Theta();
	    
            
            vector<TRootJet*> selectedJets= selection.GetSelectedJets(true);
            
//            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
            vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	    
//            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(selectedJets);
	    vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
	    
	    vector<TRootJet*> selectedBJets_CSVM; // B-jets at the Tight working point
	    vector<TRootJet*> selectedLightJets; // light-Jets, to be filled afer b-tagging
	     
            
            nElectrons=0;
            for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
                pX_electron[nElectrons]=selectedElectrons[iele]->Px();
                pY_electron[nElectrons]=selectedElectrons[iele]->Py();
                pZ_electron[nElectrons]=selectedElectrons[iele]->Pz();
                E_electron[nElectrons]=selectedElectrons[iele]->E();
                Double_t isocorr=0;
                
                // get isolation out, start by getting pu corrections
                if(selectedElectrons[iele]->puChargedHadronIso()>0){
                    isocorr = selectedElectrons[iele]->puChargedHadronIso();
                    
                }
                else{
                    // go through loads of pain to get rho correction, no function available. code below taken from TRootElectron selector in TopTreeAnalysisBase/*/Selector.cc
                    double EffectiveArea = 0.;
                    
                    // HCP 2012 updated for electron conesize = 0.3, taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 0.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 1.0 ) EffectiveArea = 0.130;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 1.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 1.479 ) EffectiveArea = 0.137;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 1.479 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.0 ) EffectiveArea = 0.067;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.2 ) EffectiveArea = 0.089;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.2 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.3 ) EffectiveArea = 0.107;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.3 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.4 ) EffectiveArea = 0.110;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.4) EffectiveArea = 0.138;
                    isocorr = event->kt6PFJets_rho()*EffectiveArea;
                }
                
                pfIso_electron[nElectrons]=(selectedElectrons[iele]->chargedHadronIso() + max( selectedElectrons[iele]->neutralHadronIso() + selectedElectrons[iele]->photonIso()  - isocorr, 0.) )/ selectedElectrons[iele]->Pt();
                charge_electron[nElectrons]=selectedElectrons[iele]->charge();
                nElectrons++;
            }
            nMuons=0;
            for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
                pX_muon[nMuons]=selectedMuons[imuo]->Px();
                pY_muon[nMuons]=selectedMuons[imuo]->Py();
                pZ_muon[nMuons]=selectedMuons[imuo]->Pz();
                E_muon[nMuons]=selectedMuons[imuo]->E();
                pfIso_muon[nMuons]=(selectedMuons[imuo]->chargedHadronIso() + max( 0.0, selectedMuons[imuo]->neutralHadronIso() + selectedMuons[imuo]->photonIso() - 0.5*selectedMuons[imuo]->puChargedHadronIso() ) ) / selectedMuons[imuo]->Pt(); // dBeta corrected


                charge_muon[nMuons]=selectedMuons[imuo]->charge();
                nMuons++;
            }
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<10; ijet++){
                pX_jet[nJets]=selectedJets[ijet]->Px();
                pY_jet[nJets]=selectedJets[ijet]->Py();
                pZ_jet[nJets]=selectedJets[ijet]->Pz();
                E_jet[nJets]=selectedJets[ijet]->E();
                nJets++;
            }
	    
	    
	    for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++){
		if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > .679)
		{
			selectedBJets_CSVM.push_back(selectedJets[iJet]);
		}
		else selectedLightJets.push_back(selectedJets[iJet]);
		
			
	     } 
	     
	    nBJets=0;
            for(int ijet=0; ijet<selectedBJets_CSVM.size() && nBJets<10; ijet++){
                pX_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Px();
                pY_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Py();
                pZ_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Pz();
                E_Bjet[nBJets]=selectedBJets_CSVM[ijet]->E();
                nBJets++;
            }
	    
	    nLJets=0;
            for(int ijet=0; ijet<selectedLightJets.size() && nLJets<10; ijet++){
                pX_Ljet[nLJets]=selectedLightJets[ijet]->Px();
                pY_Ljet[nLJets]=selectedLightJets[ijet]->Py();
                pZ_Ljet[nLJets]=selectedLightJets[ijet]->Pz();
                E_Ljet[nLJets]=selectedLightJets[ijet]->E();
                nLJets++;
            }
	    
	    
	    
	    vector<pair<int,int> >  HighestPtLept;
	    HighestPtLept.clear();
	
	    
	    for(int leptonIt = 0; leptonIt < selectedElectrons.size(); leptonIt++)
	    {
	    	    HighestPtLept.push_back(make_pair(leptonIt, selectedElectrons[leptonIt]->Pt()));
	    }
	    for(int k =0; k < selectedMuons.size() ; k++)
	    {
	    	    HighestPtLept.push_back(make_pair(20+k, selectedMuons[k]->Pt()));
	    	    
	    }		
	    
	    sort(HighestPtLept.begin(), HighestPtLept.end(),cmp_big_first);
	    if(debug)
	    {
	    	cout << "****************New event*****************" << endl; 
	    	for(int i = 0; i<HighestPtLept.size(); i++)
	    	{
	         	pair<int,int> aPair = HighestPtLept[i];
	         	if(aPair.first>19) 
		 	{
		 		int number = aPair.first-20;
				cout << "selectedMuons[number]->Pt: " << selectedMuons[number]->Pt() << endl; 
		 	}
		 	else
		 	{
		 		int number = aPair.first;
				cout << "selectedElectrons[number]->Pt: " << selectedElectrons[number]->Pt() << endl; 
		 
		 	}
	    	}
	    }
	    
	    

	    

	    
	    
	    if(channelName.find("45")!=string::npos && (nElectrons+nMuons>3))
	    { 		
	    	
		
		 myTree->Fill(); 
		 nEvents_Selected[d]++;
		

	    } // > 3 leptons
	    
	    if(channelName.find("3L")!=string::npos && (nElectrons+nMuons == 3))
	    {
	    	
		
	    	myTree->Fill();
		nEvents_Selected[d]++;
		if(debug) cout << "filled tree for 3l channel" << endl; 
		
	    }
	    
	    
	    
	    






        }			//loop on events
        
	cout<<endl;
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        //myTree->Write();
	configTree->Write("", TObject::kOverwrite);
	myTree->Write("", TObject::kOverwrite);
        fileout->Write();
        fileout->Close();
	//        delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    
    //Once everything is filled 
    //Once everything is filled ...
        if (debug)
            cout << " We ran over all the data ;-)" << endl;
    
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    //    if (debug)
    //        cout << "Treating the special plots." << endl;
    
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
