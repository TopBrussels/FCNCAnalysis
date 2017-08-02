///// NB: the command line to run this code is ./Ntuples_Processor channel tree debug ApplyElecSF ApplyMuSF ApplyPUSF ApplybTagSF ApplyGlobalSF ApplyAMC e.g is
////// ./Ntuples_Processor El_El 2lep 0 1 0 1 1 1 1
#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <sstream>
#include <TFile.h>

//user code
// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
//#include "../macros/Style.C"



using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TFile*> globalFileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,TTree*> globalttree;
map<string,MultiSamplePlot*> MSPlot;


// functions prototype
std::string intToStr (int number);
//void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG, string DecayChan, string TreeName);

void  MSPCreator (string pathPNG , string TreeName);

void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );


// faco TO BE CHANGED
Bool_t debug = false;
double DataLumi = -1;
Bool_t debug_plot = false;
bool DileptonElMu = false;
bool DileptonMuMu = false;
bool DileptonElEl = false;
string channelpostfix = "";
//applying all appropriate scale factors for individual objects if the bool is set to true
Bool_t applyElectronSF = false; 
Bool_t applyMuonSF = false; 
Bool_t applyPUSF = false;
Bool_t applyGlobalSF = false;
Bool_t applyAMC = false;
Bool_t applyBtagSF = false;
bool elecPlot = false;
bool muPlot = false;




int main(int argc, char* argv[])
{
    if (debug){
        cout << "argc = " << argc << endl;
        for(int i = 0; i < argc; i++)
        {
            cout << "argv[" << i << "] = " << argv[i] << endl;
        }
    }
    if(argc < 10) cout << " ERROR: 9 arguments expected" << endl;
    
    
    //Placing arguments in properly typed variables
    applyElectronSF = false;
    applyMuonSF = false;
    applyPUSF = false;
    applyBtagSF = false;
    applyGlobalSF = false;
    applyAMC = false;
    
    const string channel = argv[1];
    const string tree = argv[2];
    debug = strtol(argv[3],NULL,10);
    applyElectronSF = strtol(argv[4],NULL,10);
    applyMuonSF = strtol(argv[5],NULL,10);
    applyPUSF = strtol(argv[6],NULL,10);
    applyBtagSF = strtol(argv[7],NULL,10);
    applyGlobalSF = strtol(argv[8],NULL,10);
    applyAMC = strtol(argv[9],NULL,10);
    
    string runDate = "Test_DDsingleElec_23July17"; // adapted to the date of the run and according to mereged ntuples
    string xmlFileName;
    string CraneenPath;
    string decayChan;
    string Treename;
    CraneenPath = "Merged_Ntuples/";
    
    if(channel=="Mu_Mu")
    {
        cout << " --> Using the DiMuon channel..." << endl;
        channelpostfix = "_Mu_Mu_";
        xmlFileName = "config/Run2SameSignDiLepton_80X_MuMu_V10_Samples_Plotting.xml"; // to be adapted when have different cfg.xml files according to different dilepton channel
        DileptonMuMu = true;
        DataLumi  = 35826.306833965; // 80X_V4 B->H  //  pb-1
        CraneenPath += "MuMu";
        decayChan = "MuMu";
    }
    else if(channel=="El_El")
    {
        cout << " --> Using the DiElectron channel..." << endl;
        channelpostfix = "_El_El_";
        xmlFileName = "config/Run2SameSignDiLepton_80X_ElEl_V10_Samples_Plotting.xml"; // to be adapted when have different cfg.xml files according to different dilepton channel
        DileptonMuMu = true;
        DataLumi  =  35826.306833965; //RunC //11716.571246596;//  pb-1
        CraneenPath += "ElecElec";
        decayChan = "ElEl";
    }
    else if(channel=="El_Mu")
    {
        cout << " --> Using the DiElMu channel..." << endl;
        channelpostfix = "_El_Mu_";
        xmlFileName = "config/Run2SameSignDiLepton_80X_ElMu_V10_Samples_Plotting.xml"; // to be adapted when have different cfg.xml files according to different dilepton channel
        DileptonElMu = true;
        DataLumi  = 35826.306833965;//  pb-1
        CraneenPath += "ElecMu";
        decayChan = "EMu";
    }
    else
    {
        cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
        exit(1);
    }
    
    if (tree == "2lep") {
        Treename  = "2lep";
    }else if (tree == "2SSlep") {
        Treename = "2SSlep";
    }else if (tree == "2OSlep") {
        Treename = "2OSlep";
    }else if (tree == "NoZCut") {
        Treename = "NoZCut";
        //Treename = "SSLNoZCut";
        //Treename = "OSLNoZCut";
    }else if (tree == "Init") {
        Treename = "Init";
    }

    cout << "Treename =  "<< Treename <<endl;
    
    CraneenPath += "/"+runDate+"/";
    
    if(debug) cout << "input trees are --->   "<< CraneenPath << endl;
    
    string pathPNG = "Output_MSPlots/"+channel+"/";
    mkdir(pathPNG.c_str(),0777);
    pathPNG += ""+channel+"_"+runDate;
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;
    
    cout << "xmlFileName is " << xmlFileName << endl;
    cout << "NtupleFolder is " << CraneenPath << endl;
    if(debug) cout << "debugging on" << endl;
    if(applyElectronSF) cout << "Electron SF on " << endl;
    if(applyMuonSF) cout << "Muon SF on " << endl;
    if(applyPUSF) cout << "PU SF on" <<endl;
    if(applyBtagSF) cout << "BtagSF on" << endl;
    if(applyGlobalSF) cout << "applying SF " << endl;
    if(applyAMC) cout << "nlo reweight on " << endl;
    
   // if(!applyElectronSF && !applyMuonSF && !applyPUSF && !applyAMC) applyGlobalSF = false;
    else cout << "not applying SF " << endl;

    
   //// calling datasetPlotter to create MSPplots
    
    
//    /// General Variables
//    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//
//    /////////////////////////////////////////////////
//    /////****** Jet kinmatics variables *********////////
//    ////////////////////////////////////////////////
//    DatasetPlotter(100, 0, 500, "pt_jet[nJets]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(50, -3.5, 3.5, "eta_jet[nJets]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, -4.0, 4.0, "phi_jet[nJets]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(11, -0.5, 10.5, "nJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(11, -0.5, 10.5, "nCSVLbJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(11, -0.5, 10.5, "nCSVMbJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(11, -0.5, 10.5, "nCSVTbJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    /// *** float variables *** ////
//    DatasetPlotter(11, -0.5, 10.5, "MVA_nJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(11, -0.5, 10.5, "MVA_nCSVLbtagJets", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0, 300, "pt_1st_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0, 300, "pt_2nd_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(100, 0, 300, "pt_3rd_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(100, 0, 300, "pt_4th_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//
    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_1st_bjet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_2nd_bjet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_1st_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_2nd_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_3rd_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(12, -0.1, 1.1, "bdiscCSVv2_4th_jet", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);


    

    
    ///////////////////////////////////
    ////*** DiMuon Varaiables ***/////
    /////////////////////////////////
     muPlot = true;
 //   DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0, 300, "pt_muon[nMuons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(50, -3.5, 3.5, "eta_muon[nMuons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, -4.0, 4.0, "phi_muon[nMuons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    /// *** float variables *** ////
//    DatasetPlotter(100, 0.0, 300 , "pt_1st_Muon", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0.0, 300 , "pt_2nd_Muon", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0.0, 5.0 , "DeltaR_Mu0b0_DiMu", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0.0, 5.0 , "DeltaR_Mu1b0_DiMu", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
  //  DatasetPlotter(5, -2.5, 2.5 , "charge_1st_Muon", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
  //  DatasetPlotter(5, -2.5, 2.5 , "charge_2nd_Muon", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//
    
    ///////////////////////////////////
    ///** DiElectron Varaiables **////
    /////////////////////////////////
    elecPlot = true;
//    DatasetPlotter(11, -0.5, 10.5, "nElectrons", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0, 300, "pt_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(50, -3.5, 3.5, "eta_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, -4.0, 4.0, "phi_electron[nElectrons]", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
////
//    /// *** float variables *** ////
 //   DatasetPlotter(100, 0, 200, "pt_1st_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
 //   DatasetPlotter(100, 0, 200, "pt_2nd_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
 //   DatasetPlotter(50, 0.0, 5.0 , "DeltaR_Elec0b0_DiElec", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(50, 0.0, 5.0 , "DeltaR_Elec1b0_DiElec", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(5, -2.5, 2.5 , "charge_1st_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(5, -2.5, 2.5 , "charge_2nd_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
////
    ///////////////////////////////////
    ///** DiElecMu Varaiables **////
    /////////////////////////////////
    DatasetPlotter(100, 0, 200, "pt_1st_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(100, 0, 200, "pt_2nd_Electron", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 5.0 , "DeltaR_Elec0b0_DiElec", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 5.0 , "DeltaR_Elec1b0_DiElec", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);

////    //////////////////////////////
////    //// ** Met Variables ** ////
//    /////////////////////////////
//    DatasetPlotter(100, 0, 500, "met_Pt", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(50, -3.5, 3.5, "met_Eta", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, -4.0, 4.0, "met_Phi", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
////
////    //////////////////////////
//// /////////    ** MVA Variables ** ////
////    /////////////////////////
////
    DatasetPlotter(50, 0.0, 5.0 , "DeltaR_2L", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(100, -4.0, 4.0, "DeltaPhi_2L", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(100, 0.0, 300.0, "invMass_2L", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);

    DatasetPlotter(50, 0.0, 500.0, "Mass_JetPair", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 500.0, "Mass_WJetPair", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 1000.0, "Ht", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 1000.0, "St", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 500.0, "MT_lep1", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, 0.0, 500.0, "MT_lep2", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, -4.0, 4.0, "DeltaPhi_met_lep1", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
    DatasetPlotter(50, -4.0, 4.0, "DeltaPhi_met_lep2", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//
//////
////
    ////// variables for NoZcut tree ////
//    DatasetPlotter(100, 0.0, 200.0, "invMass_2SSL_Zmass", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0.0, 200.0, "invMass_2OSL_Zmass", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
//    DatasetPlotter(100, 0.0, 200.0, "invMass_2L", xmlFileName,CraneenPath,pathPNG, decayChan, Treename);
////


    
  MSPCreator (pathPNG , Treename);
    
    
}//main

//void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string PathPNG)
//void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath)
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string PathPNG, string DecayChan, string TreeName)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;
    
    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    if (debug) cout << "will start loading from xml file ..." << endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
    if (debug) cout << "finished loading from xml file ..." << endl;
    
    //***************************************************CREATING PLOTS****************************************************
    //  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
    //  outfile->cd();
    string plotname = sVarofinterest;   ///// Non Jet Split plot
    // make for loop here!!!
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
    MSPlot[plotname.c_str()]->setChannel(true, DecayChan); // to add the decay channel name on the plot
    
    
    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath , slumi;
    
    int nEntries;
   // float ScaleFactor, NormFactor, Luminosity;
  //  int varofInterest;  // used for plotting integer variables
    float varofInterest; // used for plotting float variables
    
    double varofInterest_double [20];
    
    vector<string> v;
    // to avoid modifying original string
    // first duplicate the original string and return a char pointer then free the memory
    
    char delim[] = " []";
    char * dup = strdup(sVarofinterest.c_str());
    char * token = strtok(dup, delim);
    while(token != NULL)
    {
        v.push_back(string(token));
        // the call is treated as a subsequent calls to strtok:
        // the function continues from where it left in previous invocation
        token = strtok(NULL, delim);
    }
    free(dup);
    
    if (debug) cout << v[0] << "  " << v[1] << endl;
    
    // get the desired directory
    
    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;
        filepath = TreePath + dataSetName + ".root";
        if (debug) cout<<"filepath: "<<filepath<<endl;
        
        // get the tree corresponding to the final state of interest
        string stree = "tree";
        string sSSLeptonTree = "SSLeptonTree";
        string sOSLeptonTree = "OSLeptonTree";
        string sglobaltree = "globaltree";
        string sNoZmassVetotree = "NoZmassVetoTree";
        string sSSLNoZmassVetotree = "SSLNoZmassVetoTree";
        string sOSLNoZmassVetotree = "OSLNoZmassVetoTree";
        string sInitialTree = "Initialtree";
        
        
        
        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
        string TTreename ;
        
        if (TreeName == "2lep") {
            TTreename = stree;
            cout << "the current tree name is   =   " << TTreename << endl;
        }else if (TreeName == "2SSlep") {
             TTreename = sSSLeptonTree;
            cout << "the current tree name is   =   " << TTreename << endl;
        }else if (TreeName == "2OSlep") {
             TTreename = sOSLeptonTree;
            cout << "the current tree name is   =   " << TTreename << endl;
        }else if (TreeName == "NoZCut") {
           // TTreename = sNoZmassVetotree;
            TTreename = sSSLNoZmassVetotree;
           // TTreename = sOSLNoZmassVetotree;
            cout << "the current tree name is   =   " << TTreename << endl;
        }else if (TreeName == "Init") {
            TTreename = sInitialTree;
            cout << "the current tree name is   =   " << TTreename << endl;
        }else { cout << "No tree is specified or wrong tree name: it should be 2lep or 2SSlep or 2OSlep  " <<endl;
            cout << "the current tree name is   =   " << TTreename << endl;
        
        }
        
        ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
        nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;
        
        globalFileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
        string globalTTreename = sglobaltree;
        
        if(debug) cout << "globalTTreename: " << globalTTreename.c_str() << endl;
        if(debug) cout << "globalFileObj " << globalFileObj[dataSetName.c_str()] << endl;
        globalttree[dataSetName.c_str()] = (TTree*)globalFileObj[dataSetName.c_str()]->Get(globalTTreename.c_str()); //get ttre for each dataset
        if(debug) cout << "globalttree " << globalttree[dataSetName.c_str()]<< endl;
        int globalnEntries = (int)globalttree[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries gt: "<<globalnEntries<<endl;
        
        
        
        // bo logic to set the right branch address depending on the string given as argument of the datasetplotter  (int or double[n] )
        if (v.size() == 2){
            ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&varofInterest);
            ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),varofInterest_double);
        }
        
        else if (v.size() == 1){
            if (debug)	cout << "v.size is to 1" << " and v[0] is " << v[0] << endl ;
            ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);//&varofInterest // faco To be fixed!
            
        }
        else {
            cout << "Vector of string does not have the good size!!!" << endl;
        }
        // be logic to set the right branch address depending on the string given as argument of the datasetplotter
        
        bool isData= false;
        bool isAMC = false;
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
        if(debug) cout << "isData? " << isData << endl;
        if(dataSetName.find("amc")!=string::npos) isAMC =true;
        //      cout << "isAMC? " << isAMC << endl;
        
        ///////////////////////////////////////////
        //  Event Scale Factor
        ///////////////////////////////////////////
        
        // bo of event SF
        // -----------
        
        //applying all appropriate scale factors for individual objects if the bool is set to true
        
        ///////////////////////////////////
        // determine event scalefactor ///
        //////////////////////////////////
        
        if(applyGlobalSF) cout << "    ====== Applying scale factors (not for data)=========" << endl;
        
        // get the SF from the corresponding branch
        Double_t puSF = 1. ;
        ttree[dataSetName.c_str()]->SetBranchAddress("puSF",&puSF);
        
        Double_t muonID[10];
        ttree[dataSetName.c_str()]->SetBranchAddress("MuonIDSF", &muonID);
        
        Double_t muonIso[10];
        ttree[dataSetName.c_str()]->SetBranchAddress("MuonIsoSF", &muonIso);
        
        Double_t muonTrackSF[10];
        ttree[dataSetName.c_str()]->SetBranchAddress("MuonTrackSF", &muonTrackSF);
        
        Double_t ElectronSF[10];
        ttree[dataSetName.c_str()]->SetBranchAddress("ElectronSF", &ElectronSF);
        
        Int_t nEl;
        ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons",&nEl);
        
        Int_t nMu;
        ttree[dataSetName.c_str()]->SetBranchAddress("nMuons",&nMu);
        
        Double_t BSF;
        ttree[dataSetName.c_str()]->SetBranchAddress("btagSF",&BSF);
        
        Double_t BShapeSF;
        ttree[dataSetName.c_str()]->SetBranchAddress("btagSFshape",&BShapeSF);

        
        
        ///// Global SF : SF that is determine for whole dataset
        Int_t nPosW;
        globalttree[dataSetName.c_str()]->SetBranchAddress("nofPosWeights",&nPosW);
        
        Int_t nNegW;
        globalttree[dataSetName.c_str()]->SetBranchAddress("nofNegWeights",&nNegW);
        
        Int_t nEvents;
        globalttree[dataSetName.c_str()]->SetBranchAddress("nEv",&nEvents);
        
        Int_t SumW;
        globalttree[dataSetName.c_str()]->SetBranchAddress("sumW",&SumW);
        
        
        if(debug) cout << "done setting SF addresses " << endl;
        
      
        
        /////////
        
        double globalScaleFactor= 1.;
        double nloSF = 1.;
        int nPos = 0;
        int nNeg = 0;
        int Ev = 0;
        int Weights = 0;
        if(applyAMC && isAMC && !isData)
        {
            
            for (int k = 0; k<globalnEntries; k++)
            {
                globalttree[(dataSetName).c_str()]->GetEntry(k);
                nPos += nPosW;
                nNeg += nNegW;
                Ev += nEvents;
                Weights += SumW;
                // cout << "nPos " << nPos << " vs " << nPosW << " nNeg " << nNeg << " vs " << nNegW << " + " << nPos + nNeg << " - " << nPos - nNeg  << endl;
                // cout << "nEvents " << nEvents << " vs " << Ev << " sumWeights " << SumW << " vs " << Weights << endl;
            }
            //          if(!isData) nloSF *= (double) Weights/(double) Ev; //
            if(!isData) nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
            cout << "                 1/nloSF: " << 1./nloSF << endl;
        }
        
        for (int j = 0; j<nEntries; j++)
        {
            ttree[(dataSetName).c_str()]->GetEntry(j);
            //          cout << "nEl " << nEl << " nMu " << nMu << endl;
           
            
            // here different scaling factors to be defind --> have to be branches on the tree (Scaling factors trees)
            globalScaleFactor= 1.;
            if(v.size() == 1 && sVarofinterest.find("nMuons")!=string::npos) {varofInterest = nMu;}
            if(v.size() == 1 && sVarofinterest.find("nElectrons")!=string::npos) {varofInterest = nEl;}
            
            if (applyGlobalSF && !isData)
            {
                //cout << "Applying Scale factor (applied for non Data)  " << endl;
                
                
                if(applyPUSF)
                {
                    globalScaleFactor *= puSF;
                    if (debug)
                    {
                        cout << "puSF is " << puSF << endl;
                       // cout << "the globalScaleFactor is " << globalScaleFactor << endl;
                    }
                    
                }
                
                if(applyMuonSF)
                {
                    for(unsigned int i = 0; i < nMu ; i ++)
                    {
                        
                        globalScaleFactor *= muonID[i] *  muonIso[i] * muonTrackSF[i] ;
                        if(debug)
                        {
                         cout << "MuonID SF is : " << muonID[i] << "   MuonIso SF is :  " << muonIso[i] <<"   MuonTrack SF is :  " << muonTrackSF[i] <<endl;
                        cout << "the globalScaleFactor is " << globalScaleFactor << endl;
                        }
                    }
                }
                if(applyElectronSF)
                {
                    for(unsigned int i = 0; i < nEl ; i ++)
                    {
                        
                        globalScaleFactor *= ElectronSF[i]  ;
                        if(debug)
                        {
                            cout << "electron SF is : " << ElectronSF[i] << endl;
                            cout << "the globalScaleFactor is " << globalScaleFactor << endl;
                        }
                    }
                }

                if(applyBtagSF)
                {
                   // globalScaleFactor *= BSF;
                    globalScaleFactor *= BShapeSF;
                }
                
            }

            if(applyAMC && !isData) globalScaleFactor =globalScaleFactor/ nloSF ;
            
            
            //////////////////////////////////////////////////////////////////
            
            // make MS plot for single value
            if (v.size() == 1)
            {
                if (isData)
                {
                    // for data, fill once per event, weighted with the event scale factor only ???? what??
                    //		  MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, globalScaleFactor);
                    MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, 1); // false because it is data and at this case golbalSF is 1
                }
                else
                {
                    // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
                    MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*DataLumi); //if true the code will take Eqlumi from xml file
                }
            }

        
            if (v.size() == 2)
            {
                
                // bo of loop over the number of object per entry
                if(elecPlot) varofInterest = nEl;
                if(muPlot) varofInterest = nMu;
                for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
                {
                    //		if (debug) cout << "varofInterest is " << varofInterest << endl;
                    if (isData)
                    {
                        // for data, fill once per event, weighted with the event scale factor
                        MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], false, 1);
                    }
                    else
                    {
                        // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
                        MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, globalScaleFactor*DataLumi);
                        //		  MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*Luminosity);
                    }
                    
                }
                
            }
        
        
        
        }//nEntries
    
        TCanvas *canv = new TCanvas(("canv_"+v[0]+dataSetName).c_str(),("canv_"+v[0]+dataSetName).c_str());
        
        string writename = "";
        if(isData)
        {
            writename = "data_nominal";
        }
        else
        {
            writename = dataSetName +"_nominal";
        }


    
}// loop over datasets
    
    if (debug){
        cout << "before cleaning" << endl;
        if (v.size() == 2){
            cout << " v[0] is " << v[0] << " and v[1] is " << v[1] << endl;
        }
        else if (v.size() == 1){
            cout << " v[0] is " << v[0] << endl;
        }
    }
    
    
    // clearing vector
    v.clear();
    if (debug){
        cout << "after cleaning" << endl ;
        cout << "v.size() is " << v.size() << endl;
    }
}; //// datasetplotter


// function that writes all the MSPlots created in a root file
void MSPCreator (string pathPNG , string TreeName)
{
    // cout << pathPNG.c_str() << endl;
    
    TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
    outfile->cd();
    
    
    // Loop over all the MSPlots
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string treename = "";
        if (TreeName == "2lep") {
            treename = "2l_";
        }else if (TreeName == "2SSlep") {
            treename = "2SSl_";
        }else if (TreeName == "2OSlep") {
            treename = "2OSl_";
        }else if (TreeName == "NoZCut") {
           // treename = "NoZCut_2L";
            treename = "NoZCut_2SSL_";
           // treename = "NoZCut_2OSL_";
        }else if (TreeName == "Init") {
            treename = "Init_";
        }
        
        string name = treename+"MyMSP_" + it->first;
       // string name = "2l_MyMSP_" + it->first;
       
        cout << " name " << name << endl;
       
        MultiSamplePlot *temp = it->second;
      
        if (debug){
            cout << "Saving the MSP" << endl;
            cout << " and it->first is " << it->first << endl;
        }
        temp->Draw("MyMSP_", 1, false, false, false, 10);
      //  cout<< "hello i am here 5" <<endl;
        //      name += "_test";
        if(!applyGlobalSF) name += "_noSF";
        if(!applyPUSF) name += "_noPUSF";
        if(!applyElectronSF) name += "_noElSF";
        if(!applyMuonSF) name+= "_noMuSF";
        if(!applyAMC) name+= "_noAMCcor";
        cout << "name " << name << endl;
        temp->Write(outfile, name, true ,pathPNG.c_str() , "png");
      
        //      vector<string> temp_histo = it->GetTH1FNames();
        //      for (int i_hist=0; i_hist < temp_histo.size();i_hist++  ){
        //	cout << "hist is" << temp_histo[i_hist] << endl;
        //	cout << "integral is " << it->GetTH1F(temp_histo[i_hist].GetSum()) << endl;
        //      }
    }
    
    outfile->Write("kOverwrite");
}





