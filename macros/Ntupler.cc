///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntuple.cc: This macro is intended to be an example analysis macro which works out of the box. /////
///// It should serve as the first port of call for new users of the TopBrussels framework. /////
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
    // different options for executing this macro //
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
    if(channelName.find("undefined")!=string::npos && warnings) std::cout << "[WARNING] No channel was defined" << endl;
    if(nb_Leptons == 0 && warnings) std::cout << "[WARNING] No nb of leptons was defined, default setting is 0" << endl;





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
    
    if(debug) cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
    
    /////////////////////
    // Load Datasets
    /////////////////////
    
    TTreeLoader treeLoader;
    if(debug) cout << " - Load datasets ..." << endl;
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
    
   
    if(debug) cout << " Initialized LumiReWeighting stuff" << endl;
    
    ////////////////////////////////////
    // Loop on datasets
    ////////////////////////////////////
    
    if(debug) cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        cout << " Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
        if (debug)
            std::cout<<" -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
// make root tree file name
        string roottreename = "../data/ntuples/";
roottreename+= channelName;
roottreename+= "_";
roottreename+= datasets[d]->Name();
        roottreename+="_tree.root";
// cout << "creating tree in file " << roottreename << endl;
        
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        /////////////////////
        // My tree //
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
        Double_t CSV_jet[10];
        
        Double_t missingEt;
        
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
	myTree->Branch("CSV_jet",CSV_jet,"CSV_jet[nJets]/D");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
        
// myTree->Print();

        
        
        
        //open files and load
if(debug) cout<<"LoadEvent"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);
if(debug) cout<<"LoadEvent"<<endl;
        
        
        
        /////////////////////////////////////
        /// Initialize JEC factors
        /////////////////////////////////////
   
        vector<JetCorrectorParameters> vCorrParam;
        
        /*JetCorrectorParameters *L3JetPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
JetCorrectorParameters *L2JetPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
JetCorrectorParameters *L1JetPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
// Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
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
        // Loop on events
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
                std::cout<<"Processing the "<<ievt<<"th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100 << "%)" << " +> # selected" << flush<<"\r";
            
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
/* selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); // standard TOP jet selection
selection.setMuonCuts(5,2.5,0.4,0.2,0.3,1,0.5,5,0); // standard mu selection but with looser iso
selection.setElectronCuts(10,2.5,0.4,0.02,0.5,0.3,0); // standard ele selection but with looser iso
*/
//define selection cuts --> have to be validated!!!
// From the class Selection the following functions are used:
// void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon)
// void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid)
// void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso)
//void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits)

selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1);
selection.setDiMuonCuts(10.,2.5,0.2,0.04);
selection.setDiElectronCuts(15.0,2.4,0.15,0.04,0.5,1,0.3,1);


            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
            
            if(!isGoodPV)
                continue;
            
            missingEt=mets[0]->Pt();
            
            vector<TRootJet*> selectedJets= selection.GetSelectedJets(true);
            
// vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
            vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();

// vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(selectedJets);
vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
            
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
                
                pfIso_electron[nElectrons]=(selectedElectrons[iele]->chargedHadronIso() + max( selectedElectrons[iele]->neutralHadronIso() + selectedElectrons[iele]->photonIso() - isocorr, 0.) )/ selectedElectrons[iele]->Pt();
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
                CSV_jet[nJets]=selectedJets[ijet]->btag_combinedSecondaryVertexBJetTags();
                nJets++;
            }

            if(channelName.find("45")!=string::npos && warnings)
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
        }	//loop on events
      
      
        if(channelName.find("1L3B")!=string::npos){
                 myTree->Fill();
                 if(debug) cout << "Filling 1L3B Tree" << endl;
	}
	
cout<<endl;
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        myTree->Write();
        fileout->Write();
        fileout->Close();
// delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }	//loop on datasets
    
    //Once everything is filled ...
        if (debug)
            cout << " We ran over all the data ;-)" << endl;
    
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    // if (debug)
    // cout << "Treating the special plots." << endl;
    
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << " End of the program !! " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
