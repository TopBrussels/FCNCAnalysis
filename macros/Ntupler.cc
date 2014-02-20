///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntuple.cc: This macro is intended to be an example analysis macro which works out of the box.           /////
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

//for ordering the leptons
bool cmp_big_first (pair<int,int> i, pair<int,int> j) { return i.second > j.second; }

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
        string roottreename = "../data/ntuples/";
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
	Double_t missingEt_pX;
	Double_t missingEt_pY;
	
	Double_t InvMass_4lept_HighPt; 
	Double_t InvMass_4lept_Zdecay;
	Double_t Reco_FCNC_top_Zdecay; 
	Double_t Reco_FCNC_top_HighPt;
	//Double_t Reco_SM_top; 
	//Double_t Reco_SM_W; 
	
	Double_t Phi_Higgs; 
	Double_t Eta_Higgs; 
        
        Int_t isdata;
        // various weights
        Double_t pu_weight;
        
        
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");
        
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
	myTree->Branch("missingEt_pX",&missingEt_pX,"missingEt_pX/D");
	myTree->Branch("missingEt_pY",&missingEt_pY,"missingEt_pY/D");
        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
	
	myTree->Branch("InvMass_4lept_HighPt",&InvMass_4lept_HighPt,"InvMass_4lept_HighPt/D");
        myTree->Branch("InvMass_4lept_Zdecay",&InvMass_4lept_Zdecay,"InvMass_4lept_Zdecay/D");
	myTree->Branch("Reco_FCNC_top_Zdecay",&Reco_FCNC_top_Zdecay,"Reco_FCNC_top_Zdecay/D");
	myTree->Branch("Reco_FCNC_top_HighPt",&Reco_FCNC_top_HighPt,"Reco_FCNC_top_HighPT/D");
	//myTree->Branch("Reco_SM_top",&Reco_SM_top,"Reco_SM_top/D");
	//myTree->Branch("Reco_SM_W",&Reco_SM_W,"Reco_SM_W/D");
	
	myTree->Branch("Phi_Higgs",&Phi_Higgs,"Phi_Higgs/D"); 
	myTree->Branch("Eta_Higgs",&Eta_Higgs,"Eta_Higgs/D"); 
       
	//        myTree->Print();

        
        
        
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
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
        
        
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
                std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << " +> # selected" << flush<<"\r";
            
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
	    
	    
	    
	    
	   
	    
	    
	    TLorentzVector Reco_FCNC_top_combi; 
	    TLorentzVector leptonpair_1;
	    TLorentzVector leptonpair_2;
	    TLorentzVector leptonpair_3;
	    TLorentzVector leptonfour;
	    TLorentzVector lepton_0;
	    TLorentzVector lepton_1;
	    TLorentzVector lepton_2;
	    TLorentzVector lepton_3;
	    TLorentzVector lepton_4;
	    lepton_0.Clear(); 
	    lepton_1.Clear();
	    lepton_2.Clear();
	    lepton_3.Clear();
	    lepton_4.Clear();
	    leptonpair_1.Clear();
	    leptonpair_2.Clear();
	    leptonpair_3.Clear();
	    leptonfour.Clear();
	    
	    
	
	    
	    InvMass_4lept_HighPt =-10;
	    if(nElectrons+nMuons>3)
	    { 
	    	for(int j = 0; j <HighestPtLept.size() ; j++)
	    	{
	       	    pair<int,int> aPair = HighestPtLept[j]; 
	       	    if(aPair.first>19)
	       	    {
	    	    	int number = aPair.first-20;
	    	    	if(j == 0) lepton_0.SetPxPyPzE(selectedMuons[number]->Px(),selectedMuons[number]->Py(),selectedMuons[number]->Pz(),selectedMuons[number]->Energy());
	    	    	if(j == 1) lepton_1.SetPxPyPzE(selectedMuons[number]->Px(),selectedMuons[number]->Py(),selectedMuons[number]->Pz(),selectedMuons[number]->Energy());
	    	    	if(j == 2) lepton_2.SetPxPyPzE(selectedMuons[number]->Px(),selectedMuons[number]->Py(),selectedMuons[number]->Pz(),selectedMuons[number]->Energy());
	    	    	if(j == 3) lepton_3.SetPxPyPzE(selectedMuons[number]->Px(),selectedMuons[number]->Py(),selectedMuons[number]->Pz(),selectedMuons[number]->Energy());
	       
	       	    }
	       	    else
	       	    {
	    	    	int number = aPair.first;
	    	    	if(j == 0) lepton_0.SetPxPyPzE(selectedElectrons[number]->Px(),selectedElectrons[number]->Py(),selectedElectrons[number]->Pz(),selectedElectrons[number]->Energy());
	    	    	if(j == 1) lepton_1.SetPxPyPzE(selectedElectrons[number]->Px(),selectedElectrons[number]->Py(),selectedElectrons[number]->Pz(),selectedElectrons[number]->Energy());
	    	    	if(j == 2) lepton_2.SetPxPyPzE(selectedElectrons[number]->Px(),selectedElectrons[number]->Py(),selectedElectrons[number]->Pz(),selectedElectrons[number]->Energy());
	    	    	if(j == 3) lepton_3.SetPxPyPzE(selectedElectrons[number]->Px(),selectedElectrons[number]->Py(),selectedElectrons[number]->Pz(),selectedElectrons[number]->Energy());
	       
	       
	       	    }    
	    
	    	}
	    
	    	leptonpair_1 = lepton_0 + lepton_1; 
	    	leptonpair_2 = lepton_2 + lepton_3;
	    	leptonfour = leptonpair_1 + leptonpair_2;
	    	if(debug) cout << "leptonfour.M()= " << leptonfour.M() << endl; 
	    	InvMass_4lept_HighPt = leptonfour.M();
		 
		
		
	    }
	    
	    
	    Double_t DeltaR_LightJet_Higgs = 1000; 
	    Reco_FCNC_top_combi.Clear();
	    Reco_FCNC_top_HighPt= -10;
	    if(nElectrons+nMuons>3 & nLJets>0)
	    { 
	    	if(debug) cout << "In nElectrons + nMuons > 3" << endl; 
		for( int i = 0; i <selectedLightJets.size(); i++)
		{
			TLorentzVector LightJet;
			LightJet.SetPxPyPzE(selectedLightJets[i]->Px(),selectedLightJets[i]->Py(),selectedLightJets[i]->Pz(), selectedLightJets[i]->Energy());
			Double_t Phi = selectedLightJets[i]->Phi(); 
			Double_t Eta = selectedLightJets[i]->Phi();
			
			Double_t Phi_H = leptonfour.Phi(); 
			Double_t Theta_H = leptonfour.Theta(); 
			Double_t Eta_H = -TMath::Log(TMath::Tan(Theta_H/2)); 
			
			Double_t Delta_Eta = fabs(Eta_H - Eta); 
			Double_t Delta_Phi = fabs(Phi_H - Phi); 
			Double_t Delta_R = TMath::Sqrt(Delta_Eta*Delta_Eta + Delta_Phi*Delta_Phi);
			if(Delta_R < DeltaR_LightJet_Higgs) 
			{
				DeltaR_LightJet_Higgs = Delta_R; 
				Reco_FCNC_top_combi = leptonfour + LightJet; 
				Reco_FCNC_top_HighPt = Reco_FCNC_top_combi.M();
			}
				
			  
	    	}
	    }
	    if(debug) cout << "Out nElectrons + nMuons > 3" << endl; 	   
	    
	   
	    lepton_0.Clear(); 
	    lepton_1.Clear();
	    lepton_2.Clear();
	    lepton_3.Clear();
	    lepton_4.Clear();
	    leptonpair_1.Clear();
	    leptonpair_2.Clear();
	    leptonfour.Clear();
	    
	    InvMass_4lept_Zdecay =-10;
	    Double_t leptonPhi = -100; 
	    Double_t leptonPt = -10;
	    if(nElectrons+nMuons>3)
	    { 		
	    	    if(debug) cout << "In nElectrons + nMuons > 3" << endl; 
		    bool leptons = false; 
		    bool chargereq1=false; 
		    bool chargereq2=false;
		    bool chargereq3=false;
		    bool chargereq4=false;
		    bool chargereq5=false;
		    bool chargereq6=false;
		    bool chargereq7=false;
		    bool chargereq8=false;
		    bool chargereq9=false;
		    bool chargereq10=false;
		    
		    
		    //select ZZ events, muon flavours come in pairs and OS
		    if( nMuons == 2 & nElectrons == 2){ 
		         if(selectedMuons[0]->charge()!=selectedMuons[1]->charge())   chargereq1 = true; 
			 if(selectedElectrons[0]->charge()!=selectedElectrons[1]->charge())   chargereq2 = true;
			
			 if(debug) cout << "In nMuons == 2 & nElectrons == 2 " << endl; 
			 if(chargereq1 && chargereq2  && !leptons)
			 {
			 	leptons=true; 
		 
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 
		 	 	
		 	 
		         }
			 
			 leptonpair_1 = lepton_0 + lepton_1; 
		 	 leptonpair_2 = lepton_2 + lepton_3;
			 leptonfour = leptonpair_1 + leptonpair_2;
			 
			 
		    } 
		    if( nMuons == 2 & nElectrons == 3)
		    { 
		    	 if(debug) cout << "In nMuons == 2 & nElectrons == 3 " << endl; 
		         if(selectedMuons[0]->charge()!=selectedMuons[1]->charge())   chargereq1 = true; 
			 if(selectedElectrons[0]->charge()!=selectedElectrons[1]->charge())   chargereq2 = true;
			 if(selectedElectrons[0]->charge()!=selectedElectrons[2]->charge())   chargereq3 = true;  
		 	 if(selectedElectrons[2]->charge()!=selectedElectrons[1]->charge())   chargereq4 = true; 
			 
			 if(chargereq1 && chargereq2 && !leptons) // Mu0 and Mu1    El0 and El1
			 {
			 	leptons=true; 
		                 
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	 
				lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	 
				lepton_2.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	 
				lepton_3.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 		 
			        lepton_4.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	 	
		 	 
		         }
			 if(chargereq1 && chargereq3 && !leptons) // Mu0 and Mu1    El0 and El2
			 {
			 	leptons=true;    
		 		
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 
			 }
			 if(chargereq1 && chargereq4 && !leptons) // Mu0 and Mu1    El2 and El1
			 {
			 	leptons=true; 
		 		
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 
			 }
			  
			 leptonpair_1 = lepton_0 + lepton_1; 
		 	 leptonpair_2 = lepton_2 + lepton_3;
			 leptonPhi = lepton_4.Phi(); 
			 if(debug) cout << "leptonPhi = " << leptonPhi << endl; 
			 leptonPt = lepton_4.Pt(); 
			 if(debug) cout << "leptonPt = " << leptonPt << endl;
			 leptonfour = leptonpair_1 + leptonpair_2;
			 
			 
		 } 
		 if( nMuons == 3 & nElectrons == 2)
		 {
		 	 if(debug) cout << "In nMuons == 3 & nElectrons == 2 " << endl; 
			 if(selectedElectrons[0]->charge()!=selectedElectrons[1]->charge())   chargereq1 = true; 
			 if(selectedMuons[0]->charge()!=selectedMuons[1]->charge())   chargereq2 = true;
			 if(selectedMuons[0]->charge()!=selectedMuons[2]->charge())   chargereq3 = true;  
		 	 if(selectedMuons[2]->charge()!=selectedMuons[1]->charge())   chargereq4 = true; 
			 
			 if(chargereq1 && chargereq2 && !leptons) //El0 and El1 Mu0 and Mu1
			 {
			 	leptons=true; 
		 
		 	 	lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 		lepton_4.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 
		 	 	
		 	 
		         }
			 else if(chargereq1 && chargereq3 && !leptons) //El0 and El1 Mu0 and Mu2
			 {
			 	leptons=true; 
		 
		 	 	lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 
			 }
			 else if(chargereq1 && chargereq4 && !leptons) //El0 and El1 Mu2 and Mu1
			 {
			 	leptons=true; 
		 
		 	 	lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 
			 }
			 leptonpair_1 = lepton_0 + lepton_1; 
		 	 leptonpair_2 = lepton_2 + lepton_3;
		 	 leptonfour = leptonpair_1 + leptonpair_2;
		 	 leptonPhi = lepton_4.Phi(); 
			 leptonPt = lepton_4.Pt();
		 
		 }
		 if( nMuons == 5 )
		 {
		 	 if(debug) cout << "In nMuons == 5 " << endl; 
			 if(selectedMuons[0]->charge()!=selectedMuons[1]->charge())   chargereq1 = true; 
			 if(selectedMuons[0]->charge()!=selectedMuons[2]->charge())   chargereq2 = true;
			 if(selectedMuons[0]->charge()!=selectedMuons[3]->charge())   chargereq3 = true;
			 if(selectedMuons[0]->charge()!=selectedMuons[4]->charge())   chargereq4 = true;  
		 	 if(selectedMuons[2]->charge()!=selectedMuons[1]->charge())   chargereq5 = true;
			 if(selectedMuons[3]->charge()!=selectedMuons[1]->charge())   chargereq6 = true;
			 if(selectedMuons[4]->charge()!=selectedMuons[1]->charge())   chargereq7 = true;
			 if(selectedMuons[2]->charge()!=selectedMuons[3]->charge())   chargereq8 = true; 
			 if(selectedMuons[2]->charge()!=selectedMuons[4]->charge())   chargereq9 = true;   
			 if(selectedMuons[4]->charge()!=selectedMuons[3]->charge())   chargereq10 = true; 
		 	 
			 if(chargereq1 && chargereq8  && !leptons) //0 and 1     2 and 3
			 {
			 	leptons=true; 
		 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 		lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
				
		 	 }
			 if(chargereq1 && chargereq9  && !leptons) //0 and 1     2 and 4
			 {
			 	leptons=true; 
		 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 		lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 }
			 if(chargereq1 && chargereq10  && !leptons) //0 and 1	  3 and 4
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		    	 }
			 if(chargereq2 && chargereq6  && !leptons) //0 and 2	 1 and 3
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		    	 }
			 if(chargereq2 && chargereq7  && !leptons) //0 and 2	 1 and 4
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		    	 }
			 if(chargereq2 && chargereq10  && !leptons) //0 and 2	  3 and 4
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		    	 }
			 if(chargereq3 && chargereq5  && !leptons) //0 and 3	 1 and 2
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		    	 }
			 if(chargereq3 && chargereq7  && !leptons) //0 and 3	 1 and 4
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		    	 }
			 if(chargereq3 && chargereq9  && !leptons) //0 and 3	 4 and 2
			 {
			 	leptons=true; 
		 	 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	        lepton_1.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
				lepton_4.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		    	 }
			 if(chargereq4 && chargereq5  && !leptons) //0 and 4	  1 and 2
			 {
			 	 leptons=true; 
		 	 
		 	 	 lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	 lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	 lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	 lepton_1.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
			 	 lepton_4.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		    	 }
			 if(chargereq4 && chargereq6  && !leptons) //0 and 4	  1 and 3
			 {
			 	 leptons=true; 
		 	 
		 	 	 lepton_3.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	 lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 	 lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	 lepton_1.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
			 	 lepton_4.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		    	 }
			 if(chargereq4 && chargereq8  && !leptons) //0 and 4	  3 and 2
			 {
			 	 leptons=true; 
		 	 
		 	 	 lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	 lepton_2.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	 lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	 lepton_1.SetPxPyPzE(selectedMuons[4]->Px(),selectedMuons[4]->Py(),selectedMuons[4]->Pz(),selectedMuons[4]->Energy());
			 	 lepton_4.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		    	 }
			 
		 	 leptonpair_1 = lepton_0 + lepton_1; 
		 	 leptonpair_2 = lepton_2 + lepton_3;
		 	 leptonfour = leptonpair_1 + leptonpair_2;
		 	 leptonPhi = lepton_4.Phi(); 
			 leptonPt = lepton_4.Pt();
		 	
		 	 
		 }
		  if( nElectrons == 5)
		  {
		 	if(debug) cout << "nElectrons == 5 " << endl; 
			if(selectedElectrons[0]->charge()!=selectedElectrons[1]->charge())   chargereq1 = true; 
			if(selectedElectrons[0]->charge()!=selectedElectrons[2]->charge())   chargereq2 = true;
			if(selectedElectrons[0]->charge()!=selectedElectrons[3]->charge())   chargereq3 = true;
			if(selectedElectrons[0]->charge()!=selectedElectrons[4]->charge())   chargereq4 = true;  
		 	if(selectedElectrons[2]->charge()!=selectedElectrons[1]->charge())   chargereq5 = true;
			if(selectedElectrons[3]->charge()!=selectedElectrons[1]->charge())   chargereq6 = true;
			if(selectedElectrons[4]->charge()!=selectedElectrons[1]->charge())   chargereq7 = true;
			if(selectedElectrons[2]->charge()!=selectedElectrons[3]->charge())   chargereq8 = true; 
			if(selectedElectrons[2]->charge()!=selectedElectrons[4]->charge())   chargereq9 = true;   
			if(selectedElectrons[4]->charge()!=selectedElectrons[3]->charge())   chargereq10 = true; 
		 	
			if(chargereq1 && chargereq8  && !leptons) //0 and 1	2 and 3
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		    	}
			if(chargereq1 && chargereq9  && !leptons) //0 and 1	2 and 4
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		    	}
			if(chargereq1 && chargereq10  && !leptons) //0 and 1	 3 and 4
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			
		    	}
			if(chargereq2 && chargereq6  && !leptons) //0 and 2	1 and 3
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
			
		    	}
			if(chargereq2 && chargereq7  && !leptons) //0 and 2	1 and 4
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		    	}
			if(chargereq2 && chargereq10  && !leptons) //0 and 2	 3 and 4
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		    	}
			if(chargereq3 && chargereq5  && !leptons) //0 and 3	1 and 2
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		    	}
			if(chargereq3 && chargereq7  && !leptons) //0 and 3	1 and 4
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		    	}
			if(chargereq3 && chargereq9  && !leptons) //0 and 3	4 and 2
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		    	}
			if(chargereq4 && chargereq5  && !leptons) //0 and 4	1 and 2
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		    	}
			if(chargereq4 && chargereq6  && !leptons) //0 and 4	1 and 3
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		    	}
			if(chargereq4 && chargereq8  && !leptons) //0 and 4	3 and 2
			{
			       leptons=true; 
		 	
		 	       lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
		 	       lepton_2.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
		 	       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
		 	       lepton_1.SetPxPyPzE(selectedElectrons[4]->Px(),selectedElectrons[4]->Py(),selectedElectrons[4]->Pz(),selectedElectrons[4]->Energy());
			       lepton_4.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
		    	}
			
		 	leptonpair_1 = lepton_0 + lepton_1; 
		 	leptonpair_2 = lepton_2 + lepton_3;
		 	leptonfour = leptonpair_1 + leptonpair_2;
		 	leptonPhi = lepton_4.Phi(); 
			leptonPt = lepton_4.Pt();
		 	
		 }
		 if( nMuons == 4)
		 {
		 	 if(debug) cout << "In nMuons == 4 " << endl; 
			 if(selectedMuons[0]->charge()!=selectedMuons[1]->charge())   chargereq1 = true; 
			 if(selectedMuons[0]->charge()!=selectedMuons[2]->charge())   chargereq2 = true;
			 if(selectedMuons[0]->charge()!=selectedMuons[3]->charge())   chargereq3 = true;
			 if(selectedMuons[2]->charge()!=selectedMuons[1]->charge())   chargereq4 = true;
			 if(selectedMuons[3]->charge()!=selectedMuons[1]->charge())   chargereq5 = true;
			 if(selectedMuons[3]->charge()!=selectedMuons[2]->charge())   chargereq6 = true;
			 
			 if(chargereq1 && chargereq6  && !leptons) //0 and 1     2 and 3
			 {
			 	leptons=true; 
		 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 }
			 if(chargereq2 && chargereq5  && !leptons) //0 and 2     1 and 3
			 {
			 	leptons=true; 
		 
		 	 	lepton_3.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_1.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 }
			 if(chargereq3 && chargereq4  && !leptons) //0 and 3     2 and 1
			 {
			 	leptons=true; 
		 
		 	 	lepton_1.SetPxPyPzE(selectedMuons[3]->Px(),selectedMuons[3]->Py(),selectedMuons[3]->Pz(),selectedMuons[3]->Energy());
		 	 	lepton_3.SetPxPyPzE(selectedMuons[2]->Px(),selectedMuons[2]->Py(),selectedMuons[2]->Pz(),selectedMuons[2]->Energy());
		 	 	lepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
		 	 	lepton_2.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
		 	 }
			 
			 
		 	 leptonpair_1 = lepton_0 + lepton_1; 
		 	 leptonpair_2 = lepton_2 + lepton_3;
			 leptonfour = leptonpair_1 + leptonpair_2;
			 
		  
		
		}
		if(nElectrons == 4)
		{
			if(debug) cout << "In nElectrons == 4 " << endl; 
			if(selectedElectrons[0]->charge()!=selectedElectrons[1]->charge())   chargereq1 = true; 
		        if(selectedElectrons[0]->charge()!=selectedElectrons[2]->charge())   chargereq2 = true;
		        if(selectedElectrons[0]->charge()!=selectedElectrons[3]->charge())   chargereq3 = true;
		        if(selectedElectrons[2]->charge()!=selectedElectrons[1]->charge())   chargereq4 = true;
		        if(selectedElectrons[3]->charge()!=selectedElectrons[1]->charge())   chargereq5 = true;
		        if(selectedElectrons[3]->charge()!=selectedElectrons[2]->charge())   chargereq6 = true;
		        
		        if(chargereq1 && chargereq6  && !leptons) //0 and 1	2 and 3
		        {
		               leptons=true; 
		
			       lepton_3.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_2.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
			       lepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			}
		        if(chargereq2 && chargereq5  && !leptons) //0 and 2	1 and 3
		        {
		               leptons=true; 
		
			       lepton_3.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_1.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
			       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			}
		        if(chargereq3 && chargereq4  && !leptons) //0 and 3	2 and 1
		        {
		               leptons=true; 
		
			       lepton_1.SetPxPyPzE(selectedElectrons[3]->Px(),selectedElectrons[3]->Py(),selectedElectrons[3]->Pz(),selectedElectrons[3]->Energy());
			       lepton_3.SetPxPyPzE(selectedElectrons[2]->Px(),selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(),selectedElectrons[2]->Energy());
			       lepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
			       lepton_2.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
			}
		        
		        
			leptonpair_1 = lepton_0 + lepton_1; 
			leptonpair_2 = lepton_2 + lepton_3;
			leptonfour = leptonpair_1 + leptonpair_2;
		
		}
		 

	    	if(debug) cout << "Zdecay leptonfour.M()= " << leptonfour.M() << endl; 
	    	InvMass_4lept_Zdecay = leptonfour.M();
		
		
	    }
	    
	    
	    vector<pair<int,bool> >  LightJets_Paired;
	    LightJets_Paired.clear();
	
	    
	    for(int It = 0; It < selectedLightJets.size(); It++)
	    {
	    	    LightJets_Paired.push_back(make_pair(It, false));
	    }
	    
	    DeltaR_LightJet_Higgs = 1000; 
	    Reco_FCNC_top_combi.Clear();
	    Reco_FCNC_top_Zdecay= -10;
	    if(nElectrons+nMuons>3 & nLJets>0)
	    { 
	    	for( int i = 0; i <selectedLightJets.size(); i++)
		{
			pair<int,int> aPair = LightJets_Paired[i];
			pair<int,int> pPair; 
			if(i == 0) pPair = LightJets_Paired[i];
			else pPair = LightJets_Paired[i-1];
			
			TLorentzVector LightJet;
			LightJet.SetPxPyPzE(selectedLightJets[i]->Px(),selectedLightJets[i]->Py(),selectedLightJets[i]->Pz(), selectedLightJets[i]->Energy());
			Double_t Phi = selectedLightJets[i]->Phi(); 
			Double_t Eta = selectedLightJets[i]->Phi();
			
			Double_t Phi_H = leptonfour.Phi(); 
			Double_t Theta_H = leptonfour.Theta(); 
			Double_t Eta_H = -TMath::Log(TMath::Tan(Theta_H/2)); 
			
			Double_t Delta_Eta = fabs(Eta_H - Eta); 
			Double_t Delta_Phi = fabs(Phi_H - Phi); 
			Double_t Delta_R = TMath::Sqrt(Delta_Eta*Delta_Eta + Delta_Phi*Delta_Phi);
			if(Delta_R < DeltaR_LightJet_Higgs) 
			{
				DeltaR_LightJet_Higgs = Delta_R; 
				Reco_FCNC_top_combi = leptonfour + LightJet; 
				Reco_FCNC_top_Zdecay = Reco_FCNC_top_combi.M();
				Phi_Higgs = Phi_H; 
				Eta_Higgs = Eta_H;
				
				aPair.second = true; 
				if(i != 0) pPair.second = false; 
			}
				
			  
	    	}
	    }
	    
	    /*
	    Reco_SM_W = -10;
	    bool SM_W = false; 
	    if(nElectrons+nMuons > 4)
	    {
	    	Double_t angle = fabs(leptonPhi - missingEt_Phi); 
	    	Reco_SM_W = TMath::Sqrt(2*leptonPt*missingEt*(1 - cos(angle)));
		SM_W = true; 
	    }
	    else if(nElectrons+nMuons > 3 && nLJets > 2)
	    {
	    	float counter = 0; 
		for(int iJet = 0; iJet < selectedLightJets.size(); iJet++)
		{
		    pair<int,int> aPair = LightJets_Paired[iJet];
		    
		    if(!aPair.second && counter < 2)
		    {	
		    	counter++;
			
		    	TLorentzVector LightJet;
		        LightJet.SetPxPyPzE(selectedLightJets[iJet]->Px(),selectedLightJets[iJet]->Py(),selectedLightJets[iJet]->Pz(),selectedLightJets[iJet]->Energy());
			leptonpair_3 += LightJet;
		    
		    }
		}
	    	Reco_SM_W = leptonpair_3.M();
		SM_W = true; 
	    
	    }
	    
	    Reco_SM_top = -10;
	    if(SM_W && selectedBJets_CSVM.size() >0)
	    {
	    	Reco_SM_top = Reco_SM_W
	    
	    
	    }
	    */	    
	    	    
	    
	    
	    

            if(channelName.find("45")!=string::npos )
	    {
	    	if(nElectrons+nMuons > nb_Leptons ){
                	myTree->Fill();
                	if(debug) cout << "found " << nMuons << " muons and " << nElectrons << " electrons!" << endl;
            	}
	    }
	    else{
            	if(nElectrons+nMuons== nb_Leptons ){
                	myTree->Fill();
                	if(debug) cout << "found " << nMuons << " muons and " << nElectrons << " electrons!" << endl;
            	}
            }
        }			//loop on events
        
	cout<<endl;
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        myTree->Write();
        fileout->Write();
        fileout->Close();
	//        delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    
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