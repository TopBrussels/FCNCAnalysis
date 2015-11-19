////////////////////////////////////////////////////////////
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

//used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

////---MultiSample plots ---//

map<string,MultiSamplePlot*> MSPlot;

/////
using namespace std;
using namespace reweight;
using namespace TopTree;

/// Some variables from POG/PAG
float workingpointvalue_Loose = 0.605;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Medium = 0.890;//working points updated to 2015 BTV-POG recommendations.
float workingpointvalue_Tight = 0.970;//working points updated to 2015 BTV-POG recommendations.


int main (int argc, char *argv[])
{
    
    //  string rootFileName = "testAnalyser_output.root";
    
    clock_t start = clock();
    
    
    /////////////////////
    ///  Configuration
    /////////////////////
    
    bool eventSelected = false;
    int nofSelectedEvents = 0;
    bool Fake_Electrons = false;
    bool Charge_misID = false;
    std::string sWPMuon = "Tight";
    std::string sWPElectron = "Tight";
    
    /// xml file
    string xmlFileName ="config/Run2SameSignDiLepton_samples.xml";
    
    if (argc > 1)
        xmlFileName = (string)argv[1];
    
    const char *xmlfile = xmlFileName.c_str();
    
    cout << " - Using config file " << xmlfile << endl;
    
    //Configuration output format
    //Configuration output format
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
//    TTree *configTree = new TTree("configTree","configuration Tree");
//    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
//    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
//    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
//    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
    
    ////////////////////////////////////
    ///  AnalysisEnvironment
    ////////////////////////////////////
    
    AnalysisEnvironment anaEnv;
    cout << " - Loading environment ..." << endl;
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    
    cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
    
    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    int verbose = anaEnv.Verbose;
    float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
    cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
    
    ////////////
    /// object selection and identification
    //////////////////////
    
    int PVertexNdofCut = 4; // anaEnv.PVertexNdofCut;
    int PVertexZCut =24;// anaEnv.PVertexZCut;
    int PVertexRhoCut = 2; // anaEnv.PVertexRhoCut;
    
    /////---- Muon selection --- ///
    int MuonPtCut = 10;  //anaEnv.MuonPtCutSR;
    int MuonEtaCut = 2.4;  //anaEnv.MuonEtaCutSR;
    int MuonRelIsoCut = 0.4; //anaEnv.MuonRelIsoCutSR;
    std::string WPMuon = sWPMuon; // https://indico.cern.ch/event/450085/contribution/4/attachments/1169767/1688138/20151013_MOC.pdf
    std::string CampaignMuon = "Spring15";
    
    ///// ----- Electron Selection ---- ///
    int ElectronPtCut = 12.; //anaEnv.ElectronPtCut;
    int ElectronEtaCut = 2.4; //anaEnv.ElectronEtaCut;
    std::string CampaignElectron = "Spring15_25ns";
    std::string WPElectron = sWPElectron; // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    int cutsBasedElectron = 1;
    
    ///// ------ Jet Selection  ---- ////
    int JetsPtCut = 40; //anaEnv.JetsPtCutSR;
    int applyJetID = anaEnv.applyJetID;
    int JetsEtaCut = 2.4;
    std::string WPJet = "Tight"; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
    
    
    /////////////////////
    ///  Load Datasets
    /////////////////////
    
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Loading datasets ..." << endl;
    
    treeLoader.LoadDatasets(datasets, xmlfile);
    cout << "Number of datasets: " << datasets.size() << endl;
    
    for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    
    float Luminosity =oldLuminosity;
    
    
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
        string dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")==0 || dataSetName.find("data")==0 || dataSetName.find("DATA")==0)
        {
           Luminosity = datasets[d]->EquivalentLumi();
            cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        }
    }
    //  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
     // the previous doesn't work
    cout << "lumi is " << Luminosity << endl;
    
    // set rootfile to store controlplots
    //TFile *fout = new TFile(rootFileName.c_str(), "RECREATE");
    
    
    
    
    //Global variable
    TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    
//    TTreeLoader treeLoader;
//    vector < Dataset* > datasets;
//    
//    
//    cout << " - Loading datasets ..." << endl;
//    treeLoader.LoadDatasets(datasets, xmlfile);
//    cout << "Number of datasets: " << datasets.size() << endl;
//    for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
//    
//    float Luminosity = oldLuminosity;
//    
//    
//    
//    for (unsigned int d = 0; d < datasets.size (); d++)
//    {
//        if (Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
//        string dataSetName = datasets[d]->Name();
//    }
//    if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
//    
//    
//    
//    
//    //Global variable
//    TRootEvent* event = 0;
//    
//    //nof selected events
//    double NEvtsData = 0;
//    Double_t *nEvents = new Double_t[datasets.size()];
    
    
    ///////////////////////////////////
    /// MultisamplePlots            ///
    //////////////////////////////////
    MSPlot["cutFlow"]              = new MultiSamplePlot(datasets, "cutFlow", 15, -0.5, 14.5, "cutFlow");
    MSPlot["BasecutFlow"]              = new MultiSamplePlot(datasets, "BasecutFlow", 15, -0.5, 14.5, "Baseline cutFlow");
    
    MSPlot["init_NbOfVertices"]    = new MultiSamplePlot(datasets, "init_NbOfVertices", 60, -0.5, 59.5, "initial nb of vertices");
    MSPlot["init_NbOfJets"]	  = new MultiSamplePlot(datasets, "init_NbOfJets", 16, -0.5, 15.5, "initial nb of jets");
    MSPlot["init_NbOfCSVLJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVLJets", 16, -0.5, 15.5, "initial nb of CSV loose jets");
    MSPlot["init_NbOfCSVMJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVMJets", 16, -0.5, 15.5, "initial nb of CSV medium jets");
    MSPlot["init_NbOfCSVTJets"]    = new MultiSamplePlot(datasets, "init_NbOfCSVTJets", 16, -0.5, 15.5, "initial nb of CSV tight jets");
    MSPlot["init_NbOfMuons"]       = new MultiSamplePlot(datasets, "init_NbOfMuons", 16, -0.5, 15.5,"initial nb of muons");
    MSPlot["init_NbOfElectrons"]   = new MultiSamplePlot(datasets, "init_NbOfElectrons", 16, -0.5, 15.5, "initial nb of electrons");
    MSPlot["init_NbOfLeptons"]     = new MultiSamplePlot(datasets, "init_NbOfLeptons", 16, -0.5, 15.5, "initial nb of leptons");
    
    MSPlot["2L_NbOfVertices"]      = new MultiSamplePlot(datasets, "2L_NbOfVertices", 60, -0.5, 59.5, "After 2 lepton req.:  nb of vertices");
    MSPlot["2L_NbOfJets"]          = new MultiSamplePlot(datasets, "2L_NbOfJets", 16, -0.5, 15.5, "After 2 lepton req.:  nb of jets");
    MSPlot["2L_NbOfCSVLJets"]      = new MultiSamplePlot(datasets, "2L_NbOfCSVLJets", 16, -0.5, 15.5, "After 2 lepton req.:  nb of CSV loose jets");
    MSPlot["2L_NbOfCSVMJets"]      = new MultiSamplePlot(datasets, "2L_NbOfCSVMJets", 16, -0.5, 15.5, "After 2 lepton req.:  nb of CSV medium jets");
    MSPlot["2L_NbOfCSVTJets"]      = new MultiSamplePlot(datasets, "2L_NbOfCSVTJets", 16, -0.5, 15.5, "After 2 lepton req.:  nb of CSV tight jets");
    MSPlot["2L_NbOfMuons"]         = new MultiSamplePlot(datasets, "2L_NbOfMuons", 16, -0.5, 15.5,"After 2 lepton req.:  nb of muons");
    MSPlot["2L_NbOfElectrons"]     = new MultiSamplePlot(datasets, "2L_NbOfElectrons", 16, -0.5, 15.5, "After 2 lepton req.:  nb of electrons");
    MSPlot["2L_NbOfLeptons"]       = new MultiSamplePlot(datasets, "2L_NbOfLeptons", 16, -0.5, 15.5, "After 2 lepton req.:  nb of leptons");
    
    MSPlot["2SSL_NbOfVertices"]    = new MultiSamplePlot(datasets, "2SSL_NbOfVertices", 60, -0.5, 59.5, "After 2SSL req.:  nb of vertices");
    MSPlot["2SSL_NbOfJets"]        = new MultiSamplePlot(datasets, "2SSL_NbOfJets", 16, -0.5, 15.5, "After 2SSL req.:  nb of jets");
    MSPlot["2SSL_NbOfCSVLJets"]    = new MultiSamplePlot(datasets, "2SSL_NbOfCSVLJets", 16, -0.5, 15.5, "After 2SSL req.:  nb of CSV loose jets");
    MSPlot["2SSL_NbOfCSVMJets"]    = new MultiSamplePlot(datasets, "2SSL_NbOfCSVMJets", 16, -0.5, 15.5, "After 2SSL req.:  nb of CSV medium jets");
    MSPlot["2SSL_NbOfCSVTJets"]    = new MultiSamplePlot(datasets, "2SSL_NbOfCSVTJets", 16, -0.5, 15.5, "After 2SSL req.:  nb of CSV tight jets");
    MSPlot["2SSL_NbOfMuons"]       = new MultiSamplePlot(datasets, "2SSL_NbOfMuons", 16, -0.5, 15.5,"After 2SSL req.:  nb of muons");
    MSPlot["2SSL_NbOfElectrons"]   = new MultiSamplePlot(datasets, "2SSL_NbOfElectrons", 16, -0.5, 15.5, "After 2SSL req.:  nb of electrons");
    MSPlot["2SSL_NbOfLeptons"]     = new MultiSamplePlot(datasets, "2SSL_NbOfLeptons", 16, -0.5, 15.5, "After 2SSL req.:  nb of leptons");

    
    MSPlot["2J_NbOfVertices"]      = new MultiSamplePlot(datasets, "2J_NbOfVertices", 60, -0.5, 59.5, "After 2 jet req.:  nb of vertices");
    MSPlot["2J_NbOfJets"]          = new MultiSamplePlot(datasets, "2J_NbOfJets", 16, -0.5, 15.5, "After 2 jet req.:  nb of jets");
    MSPlot["2J_NbOfCSVLJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVLJets", 16, -0.5, 15.5, "After 2 jet req.:  nb of CSV loose jets");
    MSPlot["2J_NbOfCSVMJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVMJets", 16, -0.5, 15.5, "After 2 jet req.:  nb of CSV medium jets");
    MSPlot["2J_NbOfCSVTJets"]      = new MultiSamplePlot(datasets, "2J_NbOfCSVTJets", 16, -0.5, 15.5, "After 2 jet req.:  nb of CSV tight jets");
    MSPlot["2J_NbOfMuons"]         = new MultiSamplePlot(datasets, "2J_NbOfMuons", 16, -0.5, 15.5,"After 2 jet req.:  nb of muons");
    MSPlot["2J_NbOfElectrons"]     = new MultiSamplePlot(datasets, "2J_NbOfElectrons", 16, -0.5, 15.5, "After 2 jet req.:  nb of electrons");
    MSPlot["2J_NbOfLeptons"]       = new MultiSamplePlot(datasets, "2J_NbOfLeptons", 16, -0.5, 15.5, "After 2 jet req.:  nb of leptons");
    
    MSPlot["1BJ_NbOfVertices"]     = new MultiSamplePlot(datasets, "1BJ_NbOfVertices", 60, -0.5, 59.5, "After 1 Bjet req.:  nb of vertices");
    MSPlot["1BJ_NbOfJets"]         = new MultiSamplePlot(datasets, "1BJ_NbOfJets", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of jets");
    MSPlot["1BJ_NbOfCSVLJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVLJets", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of CSV loose jets");
    MSPlot["1BJ_NbOfCSVMJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVMJets", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of CSV medium jets");
    MSPlot["1BJ_NbOfCSVTJets"]     = new MultiSamplePlot(datasets, "1BJ_NbOfCSVTJets", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of CSV tight jets");
    MSPlot["1BJ_NbOfMuons"]        = new MultiSamplePlot(datasets, "1BJ_NbOfMuons", 16, -0.5, 15.5,"After 1 Bjet req.:  nb of muons");
    MSPlot["1BJ_NbOfElectrons"]    = new MultiSamplePlot(datasets, "1BJ_NbOfElectrons", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of electrons");
    MSPlot["1BJ_NbOfLeptons"]      = new MultiSamplePlot(datasets, "1BJ_NbOfLeptons", 16, -0.5, 15.5, "After 1 Bjet req.:  nb of leptons");
    
    
    
    ////////////////////////////////////
    ///  Selection table
    ////////////////////////////////////
    
    vector<string> CutsSelecTable;
    CutsSelecTable.push_back(string("preselected"));
    CutsSelecTable.push_back(string("trigged"));
    CutsSelecTable.push_back(string("Good PV"));
    CutsSelecTable.push_back(string("2 selected leptons"));
    CutsSelecTable.push_back(string("At least 2 jets"));
    CutsSelecTable.push_back(string("At least 1 CSV loose jet"));
    
    SelectionTable selecTable(CutsSelecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    
    if (verbose > 0)
        cout << " - SelectionTable instantiated ..." << endl;
    
    
    ////////////////////////////////////
    ///  Loop on datasets
    ////////////////////////////////////
    
    if (verbose > 0)
        cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
    
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        int isdata;
        
        if (verbose > 1)
            cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
        if (verbose > 1)
            cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
        
        // make root tree file name
        string roottreename = "Ntuples/";
        roottreename+= datasets[d]->Name();
        roottreename+="_tree.root";
        
        cout << "  - Recreate outputfile ... " << roottreename.c_str() << endl;
        // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        
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
        
        //doubles
        
        
        //vectors
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
        
        
        // define the output tree (Integer variables)
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");
        myTree->Branch("run_num",&run_num,"run_num/I");
        myTree->Branch("evt_num",&evt_num,"evt_num/I");
        myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
        myTree->Branch("nvtx",&nvtx,"nvtx/I");
        myTree->Branch("npu",&npu,"npu/I");
        
        //Set branches for doubles
        ///   myTree -> Branch("metPt", &metPt, "metPt/D");
        //    myTree -> Branch("metPx", &metPx, "metPx/D");
        //    myTree -> Branch("metPy", &metPy, "metPy/D");
        
        
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
        
        
        //open files and load
        cout << "LoadEvent" << endl;
        treeLoader.LoadDataset(datasets[d], anaEnv);
        
        nofSelectedEvents = 0;
        
        ////////////////////////////////////
        ///  Loop on events
        ////////////////////////////////////
        //some bookkeeping variables
        nEvents[d] = 0;
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
        if (verbose > 1)
            cout << "	Loop over events " << endl;
        
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
        {
            vector < TRootVertex* > vertex;
            vector < TRootMuon* > init_muons;
            vector < TRootElectron* > init_electrons;
            vector < TRootJet* > init_jets_corrected;
            vector < TRootJet* > init_jets;
            vector < TRootMET* > mets;
            //vector < TRootGenJet* > genjets;
            
            nEvents[d]++;
            
            if (ievt%500 == 0)
                std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
            
            
            
            ////////////////////
            ///  LOAD EVENT  ///
            ////////////////////
            
            TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
            
            
            ////////////////////////////
            ///  Include trigger set up here when using data
            ////////////////////////////
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if (previousFilename != currentFilename)
            {
                previousFilename = currentFilename;
                iFile++;
                cout << "File changed!!! => iFile = " << iFile << endl;
            }
            
            int currentRun = event->runId();
            
            if (previousRun != currentRun)
                previousRun = currentRun;
            
            
            run_num=event->runId();
            evt_num=event->eventId();
            lumi_num=event->lumiBlockId();
            nvtx=vertex.size();
            npu=(int)event->nTruePU();
            if( run_num > 10000){//data
                isdata=1;
            }
            bookkeeping->Fill();
            
            
            ////////////////////////////////////
            ///  DETERMINE EVENT SCALEFACTOR  ///
            /////////////////////////////////////
            
            // scale factor for the event
            float scaleFactor = 1.;
            
            
            // PU reweighting
            
            // old method
            //cout << "scalefactor " << scaleFactor << endl;
            double lumiWeight = 1; //LumiWeights.ITweight( (int)event->nTruePU() ); // currently no pile-up reweighting applied
            
            if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
                lumiWeight=1;
            
            // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
            // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
            
            scaleFactor = scaleFactor*lumiWeight;
            
            /////////////////////////
            ///  EVENT SELECTION  ///
            /////////////////////////
            
            //Declare selection instance
            Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
            
            bool isGoodPV = selection.isPVSelected(vertex, PVertexNdofCut, PVertexZCut,PVertexRhoCut); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
            
            vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(JetsPtCut, JetsEtaCut, applyJetID, WPJet);  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(MuonPtCut, MuonEtaCut, MuonRelIsoCut,WPMuon,CampaignMuon);  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(ElectronPtCut, ElectronEtaCut, WPElectron, CampaignElectron, cutsBasedElectron);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
            
            
            sort(selectedJets.begin(), selectedJets.end(),HighestPt());
            sort(selectedMuons.begin(), selectedMuons.end(), HighestPt());
            sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt());
            
            /////---- bTagging ----\\\\\\
           // vector<bool> BtagBooleans;
           // BtagBooleans.clear();
            vector<TRootJet*> selectedBCSVLJets;
            vector<TRootJet*> selectedBCSVMJets;
            vector<TRootJet*> selectedBCSVTJets;
            
            for(unsigned int i = 0; i < selectedJets.size() ; i++)
            {
                bool Btagged = false;
                TRootJet* tempJet = (TRootJet*) selectedJets[i];
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Loose)//loose WP
                {
                    Btagged = true;
                    selectedBCSVLJets.push_back(tempJet);
                }
               // BtagBooleans.push_back(Btagged);
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Medium)//medium WP
                {
                    selectedBCSVMJets.push_back(tempJet);
                }
                if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue_Tight)//tight WP
                {
                    selectedBCSVTJets.push_back(tempJet);
                }
            }
            if(verbose > 3) cout << "btagging done" << endl;
            
            
            ////--- Initinal MultiSample plots
            selecTable.Fill(d,0,scaleFactor*Luminosity);
            MSPlot["cutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor );
            MSPlot["init_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["init_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
            
            
            ///////////************\\\\\\\\\\\
            // Start Applying selection cuts \\
            ///////////////////////////////////
            eventSelected = false;
            bool diElectron = false;
            bool diMuon = false;
            bool diEMu = false;
            double InvMass_ll = 0.;
            double Zmass = 91.1876; // SM xsections at 13 TeV https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
            TLorentzVector tempLepton_0;
            TLorentzVector tempLepton_1;
            
            
            selecTable.Fill(d,0,scaleFactor);
            
            /// At the moment do not use trigger
            selecTable.Fill(d,1,scaleFactor);
            MSPlot["cutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor );
            
            if (isGoodPV)
            {
                if(verbose>3) cout << "GoodPV" << endl;
                selecTable.Fill(d,2,scaleFactor*Luminosity);
                MSPlot["cutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor );
                if (selectedMuons.size() + selectedElectrons.size()== 2) //Asking for 2 leptons
                {
                    if(verbose>3) cout << "2 Leptons "<< endl;
                    selecTable.Fill(d,3,scaleFactor);
                    
                    if (selectedElectrons.size()==2 && selectedElectrons[0]->charge()== selectedElectrons[1]->charge() && selectedElectrons[0]->Pt()>26. && selectedElectrons[1]->Pt()>= 15.) // Asking for 2 same sign electrons
                    {
                        tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                        tempLepton_1.SetPxPyPzE(selectedElectrons[1]->Px(),selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(),selectedElectrons[1]->Energy());
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        if (InvMass_ll >12. && fabs(Zmass-InvMass_ll)>15.)
                        {
                            diElectron = true;
                        }
                        
                    }//same sign dielectron channel
                    
                    if (selectedMuons.size()==2 && selectedMuons[0]->charge()==selectedMuons[1]->charge() && selectedMuons[0]->Pt() > 20. && selectedMuons[1]->Pt()>=11. )// Asking for 2 same sign Muons
                    {
                        tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                        tempLepton_1.SetPxPyPzE(selectedMuons[1]->Px(),selectedMuons[1]->Py(),selectedMuons[1]->Pz(),selectedMuons[1]->Energy());
                        InvMass_ll = (tempLepton_0+tempLepton_1).M();
                        if (InvMass_ll >12. && fabs(Zmass-InvMass_ll)>15.)
                        {
                            diMuon = true;
                        }
                        
                    }//same sign diMuon channel
                    
                    if (selectedMuons.size()==1 && selectedElectrons.size()==1 && selectedMuons[0]->charge()==selectedElectrons[0]->charge())// Asking for 2 same sign electron- muon pair
                    {
                        if (selectedMuons[0]->Pt() > selectedElectrons[0]->Pt()&& selectedMuons[0]->Pt()>= 20. && selectedElectrons[0]->Pt()>=15.)
                        {
                            tempLepton_0.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            
                        }
                        if (selectedElectrons[0]->Pt()>selectedMuons[0]->Pt() && selectedElectrons[0]->Pt()>=26. && selectedMuons[0]->Pt())
                        {
                            tempLepton_0.SetPxPyPzE(selectedElectrons[0]->Px(),selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(),selectedElectrons[0]->Energy());
                            tempLepton_1.SetPxPyPzE(selectedMuons[0]->Px(),selectedMuons[0]->Py(),selectedMuons[0]->Pz(),selectedMuons[0]->Energy());
                        }
                        
                            InvMass_ll = (tempLepton_0+tempLepton_1).M();
                            if (InvMass_ll >12.)
                            {
                                diEMu = true;
                            }
                       
                        
                    }// Same Sign electron-muon channel
                    
                    selecTable.Fill(d,3,scaleFactor*Luminosity);
                    MSPlot["cutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["BasecutFlow"]->Fill(0, datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["2L_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["2L_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                    
                    if (diMuon || diElectron || diEMu)
                    {
                        if(verbose>3) cout << " 2SSL "<< endl;
                        selecTable.Fill(d,4,scaleFactor*Luminosity);
                        MSPlot["cutFlow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor );
                        MSPlot["BasecutFlow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor );
                        MSPlot["2SSL_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["2SSL_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                        
                        
                        if (selectedJets.size() >= 2)//at least 2 jets
                        {
                            if(verbose>3) cout << " at least 2 jets " << endl;
                            selecTable.Fill(d,5,scaleFactor*Luminosity);
                            MSPlot["cutFlow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor );
                            MSPlot["BasecutFlow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor );
                            MSPlot["2J_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["2J_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
                            eventSelected = true;
                            
//                            if (selectedBCSVLJets.size()>=1)
//                            {
//                                if(verbose>3) cout << " at least 1 CSVLBjet " << endl;
//                                selecTable.Fill(d,6,scaleFactor*Luminosity);
//                                MSPlot["cutFlow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor );
//                                MSPlot["BasecutFlow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor );
//                                MSPlot["1BJ_NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfCSVLJets"]->Fill(selectedBCSVLJets.size(), datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfCSVMJets"]->Fill(selectedBCSVMJets.size(), datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfCSVTJets"]->Fill(selectedBCSVTJets.size(), datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfMuons"]->Fill(selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfElectrons"]->Fill(selectedElectrons.size(),datasets[d], true, Luminosity*scaleFactor);
//                                MSPlot["1BJ_NbOfLeptons"]->Fill(selectedElectrons.size()+selectedMuons.size(),datasets[d], true, Luminosity*scaleFactor);
//                                
//                                
//                                
//                            }//2SSL+2Jets+>=1CSVLBjet

                            
                        }//2SSL+2Jets
                    }//2SSL
                    
                    
                }//2L
                
            }//isGoodPV
            
            
            
            if (! eventSelected )
            {
                continue;
            }
            
            if(verbose>3) cout << "filling the tree" << endl;
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
            }
            nb_Muons->push_back(selectedMuons.size());
            
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
            }
            nb_Jets->push_back(selectedJets.size());
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
                tempJet_Size++;
            }
            nb_CSVLbJets->push_back(tempJet_Size);
                               
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
            }
            nb_CSVMbJets->push_back(selectedBCSVMJets.size());
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
            }
            nb_CSVTbJets->push_back(selectedBCSVTJets.size());
            
            
            
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
            
            
            //////////////////////
            ///  END OF EVENT  ///
            //////////////////////
            
        }  //loop on events
        
        cout << endl;
        cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
        
        selecTable.TableCalculator(false, true, true, true, true);
        string selectionTable = "SelectionTables/SelectionTable_"+datasets[d]->Name() +".tex";
        //selecTableSemiMu.Write(selectiontableMu.c_str());
        selecTable.Write(selectionTable.c_str(), true, true, true, true, true, true, false);
        
        myTree->Write();
        fileout->Write();
        fileout->Close();
        delete fileout; 
        //   delete myTree;
        
        
        
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }  // Loop on datasets
    
    
    selecTable.TableCalculator(false, true, true, true, true); 
    string selectionTableAll = "SelectionTables/SelectionTable_allSamples.tex";
    selecTable.Write(selectionTableAll.c_str(), true, true, true, true, true, true, false);
    
    //  fout->Close();
    
    //  delete fout;
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
    
}


