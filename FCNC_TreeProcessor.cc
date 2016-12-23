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
#include <ctime>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

//includes for MVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"


using namespace std;
using namespace TopTree;
//using namespace KINFIT;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlot_nPV;





// functions prototype
std::string intToStr (int number);
void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug);

inline bool FileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{

    int baseline_jets                 = strtol(argv[1], NULL,10);
    int baseline_bjets             = strtol(argv[2], NULL,10);
    string channel            = argv[3];
    string date            = argv[4];
    bool PVreweighing = strtol(argv[5], NULL,10);
    bool debug         =strtol(argv[6], NULL,10);
    
    bool doInclusive = false;
    string category;
    if(baseline_bjets == 0 && baseline_jets == 0)
    {
        doInclusive = true;
        category = "Inclusive";
    }
    else
    {
        category = "b"+intToStr(baseline_bjets)+"j"+intToStr(baseline_jets);
    }    

    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << "Category: " << category << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();


    cout << " ... Making the TreeProcessor .xml files " << endl;
    system("python scripts/MakeXMLforTreeProcessor.py");

    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    TString TreePath = "Merged/Ntuples" + channel + "/Ntuples" + date;
    if(!FileExists(string(TreePath+"/FCNC_1L3B__Run2_TopTree_Study_Data.root")))
    {
        system(("hadd "+TreePath+"/FCNC_1L3B__Run2_TopTree_Study_Data.root "+TreePath+"/FCNC_1L3B__Run2_TopTree_Study_Data_*.root").Data());
    }


  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
    //***************************************************LOADING DATASETS****************************************************
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  	string dataSetName, filepath;
    float Luminosity = 0;

    //***************************************************CREATING PLOT****************************************************
    //Format of MSPlots: MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float Min, float Max, string XaxisLabel, string YaxisLabel, string Text, string Units)

    MSPlot["NPV"] = new MultiSamplePlot(datasets, "NPV", 51, -0.5, 50.5, "Number of PV","Events", category); 
    MSPlot["NPU"] = new MultiSamplePlot(datasets, "NPU", 51, -0.5, 50.5, "Number of PU","Events", category); 
    MSPlot["NCSVv2Ljets"] = new MultiSamplePlot(datasets, "NCSVv2Ljets", 11, -0.5, 10.5, "Number of CSVv2 L jets","Events", category); 
    MSPlot["NCSVv2Mjets"] = new MultiSamplePlot(datasets, "NCSVv2Mjets", 11, -0.5, 10.5, "Number of CSVv2 M jets","Events", category); 
    MSPlot["NCSVv2Tjets"] = new MultiSamplePlot(datasets, "NCSVv2Tjets", 11, -0.5, 10.5, "Number of CSVv2 T jets","Events", category); 
    MSPlot["Njets"] = new MultiSamplePlot(datasets, "Njets", 11, -0.5, 10.5, "Number jets","Events", category); 
    MSPlot["LeptonPt"] = new MultiSamplePlot(datasets, "LeptonPt", 50, 20., 300., "Lepton Pt","Events", category, "GeV"); 
    MSPlot["LeptonEta"] = new MultiSamplePlot(datasets, "LeptonEta", 50, -2.5, 2.5, "Lepton eta","Events", category); 
    MSPlot["LeptonPhi"] = new MultiSamplePlot(datasets, "LeptonPhi", 50, -3.2, 3.2, "Lepton phi","Events", category); 
    MSPlot["JetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 20., 300., "Jet Pt","Events", category, "GeV"); 
    MSPlot["JetEta"] = new MultiSamplePlot(datasets, "JetEta", 50, -2.5, 2.5, "Jet eta","Events", category); 
    MSPlot["JetPhi"] = new MultiSamplePlot(datasets, "JetPhi", 50, -3.2, 3.2, "Jet phi","Events", category); 
    MSPlot["JetCSVv2"] = new MultiSamplePlot(datasets, "JetCSVv2", 50, 0., 1., "Jet CSVv2","Events", category); 
    MSPlot["JetcMVAv2"] = new MultiSamplePlot(datasets, "JetcMVAv2", 50, -1., 1., "Jet cMVAv2","Events", category); 

    MSPlot["Hmass_TOPTOPLEPHAD"] = new MultiSamplePlot(datasets, "Hmass_TOPTOPLEPHAD", 50, 60., 180., "H mass","Events", category, "GeV"); 
    MSPlot["Hmass_TOPTOPLEPHBB"] = new MultiSamplePlot(datasets, "Hmass_TOPTOPLEPHBB", 50, 60., 180., "H mass","Events", category, "GeV"); 
    MSPlot["Hmass_TOPHLEPBB_hct"] = new MultiSamplePlot(datasets, "Hmass_TOPHLEPBB_hct", 50, 60., 180., "H mass","Events", category, "GeV"); 
    MSPlot["Hmass_TOPHLEPBB_hut"] = new MultiSamplePlot(datasets, "Hmass_TOPHLEPBB_hut", 50, 60., 180., "H mass","Events", category, "GeV"); 
    MSPlot["HadTopMass_TOPTOPLEPHAD"] = new MultiSamplePlot(datasets, "HadTopMass_TOPTOPLEPHAD", 50, 60., 180., "HadTopMass mass","Events", category, "GeV");
    MSPlot["HadTopMass_TOPTOPLEPHBB"] = new MultiSamplePlot(datasets, "HadTopMass_TOPTOPLEPHBB", 50, 60., 180., "HadTopMass mass","Events", category, "GeV");
    MSPlot["HadTopMass_TOPHLEPBB_hct"] = new MultiSamplePlot(datasets, "HadTopMass_TOPHLEPBB_hct", 50, 60., 180., "HadTopMass mass","Events", category, "GeV");
    MSPlot["HadTopMass_TOPHLEPBB_hut"] = new MultiSamplePlot(datasets, "HadTopMass_TOPHLEPBB_hut", 50, 60., 180., "HadTopMass mass","Events", category, "GeV");

    MSPlot["MC_TopPt"] = new MultiSamplePlot(datasets, "MC_TopPt", 60, 80., 250., "Top Pt (Gen)","Events", category, "GeV");
    MSPlot["MC_AntiTopPt"] = new MultiSamplePlot(datasets, "MC_AntiTopPt", 60, 80., 250., "Top Pt (Gen)","Events", category, "GeV");

  
 
    //***************************************************GETTING LUMI FROM DATA IN XML****************************************************
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	  {
		    dataSetName = datasets[d]->Name();
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        Luminosity = datasets[d]->EquivalentLumi();
        }
    }
    if(Luminosity == 0)
    {
            cout << "Luminosity is 0. Please check the data-luminosity in your xml file. Exiting program..." << endl;
            return 1;
    }

  	//***********************************************RUNNING OVER DATASETS**********************************************
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	  {
		    dataSetName = datasets[d]->Name();
		    cout<<"Dataset:  :"<<dataSetName<<endl;
		    filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		    if (debug)
		    {
		        cout<<"filepath: "<<filepath<<endl;
            cout <<"Equivalent luminosity of the dataset is: " << datasets[d]->EquivalentLumi() << endl;
		    }
	
        reweight::LumiReWeighting W_nPV;
        if(PVreweighing)//Before you can apply this, you need to make the nPV distributions first by running this macro once.
        {
            string pathPlot = "MSPlots/";
            mkdir(pathPlot.c_str(),0777);
            pathPlot += "MSPlots";
            pathPlot += channel;
            pathPlot += "/";
            mkdir(pathPlot.c_str(),0777);
            pathPlot += date;
            pathPlot += "/";
            mkdir(pathPlot.c_str(),0777);
            pathPlot += category;
            pathPlot += "/";
            mkdir(pathPlot.c_str(),0777);            
            pathPlot += "Output_NPV.root";            

            if(!FileExists(pathPlot))
            {
                MakeNPV_Distributions(baseline_jets, baseline_bjets, channel, date, debug);
            }
            
            W_nPV = reweight::LumiReWeighting( pathPlot.c_str(), pathPlot.c_str(), ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_Data");    
        }

		    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		                
		                
        bool isData= false;
		    bool isAMC = false;
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        if(debug) cout << "Data found" << endl;
		        isData =true;
	      }
        else if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos || dataSetName.find("amc") !=std::string::npos) isAMC = true;


  	    //***********************************************IMPORTING VARIABLES**********************************************
		    string TTreename = "ObjectVarsTree";	
		    string TTreename_info = "NtupleInfoTree";	
		    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		    ttree[(dataSetName+TTreename_info).c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_info.c_str()); //get ntuple creation information

        int nEntries;

		    nEntries = ttree[dataSetName.c_str()]->GetEntries();
		    cout<<"                 nEntries: "<<nEntries<<endl;
		
        //----------------------------------------------//
        //Import the working points for b-tagging used to create the ntuples
        //----------------------------------------------//
	      Double_t CSVv2_workingpointvalue_Loose;
	      Double_t CSVv2_workingpointvalue_Medium;
	      Double_t CSVv2_workingpointvalue_Tight;
	      Double_t cMVA_workingpointvalue_Loose;
	      Double_t cMVA_workingpointvalue_Medium;
	      Double_t cMVA_workingpointvalue_Tight;
        
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("CSVv2_workingpointvalue_Loose",&CSVv2_workingpointvalue_Loose);
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("CSVv2_workingpointvalue_Medium",&CSVv2_workingpointvalue_Medium);
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("CSVv2_workingpointvalue_Tight",&CSVv2_workingpointvalue_Tight);
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("cMVA_workingpointvalue_Loose",&CSVv2_workingpointvalue_Loose);
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("cMVA_workingpointvalue_Medium",&CSVv2_workingpointvalue_Medium);
        ttree[(dataSetName + TTreename_info).c_str()]->SetBranchAddress("cMVA_workingpointvalue_Tight",&CSVv2_workingpointvalue_Tight);
        ttree[(dataSetName + TTreename_info).c_str()]->GetEntry(0);

        //----------------------------------------------//
        //Get The object variables + weights
        //----------------------------------------------//
        //Weights
        Double_t W_puSF;
        Double_t W_fleptonSF;
        Double_t W_btagWeight_CSVv2M_mujets_central;
        Double_t W_btagWeight_CSVv2M_mujets_up;
        Double_t W_btagWeight_CSVv2M_mujets_down;
        Double_t W_btagWeight_shape;
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
        Double_t W_MuonTrigSF_Runs273158to274093 = 1;//Used in calculation for W_MuonTrigSF
        Double_t W_MuonTrigSF_Runs274094to276097 = 1;//Used in calculation for W_MuonTrigSF
        Double_t W_ElectronIDSF; //One of the 2 components for the total electron SF
        Double_t W_ElectronRecoSF; //One of the 2 components for the total electron SF
        Double_t W_TopPtReweighing;
      
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;

	      // variables for electrons
        Double_t pt_electron;
        Double_t phi_electron;
        Double_t eta_electron;
        Double_t eta_superCluster_electron;
        Double_t E_electron;
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
        Double_t pt_muon;
        Double_t phi_muon;
        Double_t eta_muon;
        Double_t E_muon;
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
	      Int_t TOPTOPLEPHAD_JetIdx_LepTop = -99;
	      Int_t TOPTOPLEPHAD_JetIdx_HadTop = -99;
	      Int_t TOPTOPLEPHAD_JetIdx_W1 = -99;
	      Int_t TOPTOPLEPHAD_JetIdx_W2 = -99;
	      Int_t TOPTOPLEPHBB_JetIdx_LepTop = -99;
	      Int_t TOPTOPLEPHBB_JetIdx_HadTop = -99;
	      Int_t TOPTOPLEPHBB_JetIdx_H1 = -99;
	      Int_t TOPTOPLEPHBB_JetIdx_H2 = -99;
	      Int_t TOPHLEPBB_JetIdx_LepTop_hut = -99;
	      Int_t TOPHLEPBB_JetIdx_H1_hut = -99;
	      Int_t TOPHLEPBB_JetIdx_H2_hut = -99;
	      Int_t TOPHLEPBB_JetIdx_LepTop_hct = -99;
	      Int_t TOPHLEPBB_JetIdx_H1_hct = -99;
	      Int_t TOPHLEPBB_JetIdx_H2_hct = -99;
        Double_t MVA_TOPTOPLEPHAD = -999.;
        Double_t MVA_TOPTOPLEPHBB = -999.;
        Double_t MVA_TOPHLEPBB_hut = -999.;
        Double_t MVA_TOPHLEPBB_hct = -999.;

        Double_t MC_TopPt;
        Double_t MC_AntiTopPt;
        
        // Weights
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_central",&W_btagWeight_CSVv2M_mujets_central); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_up",&W_btagWeight_CSVv2M_mujets_up);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_down",&W_btagWeight_CSVv2M_mujets_down); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape",&W_btagWeight_shape); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_nloWeight",&W_nloWeight); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight1",&W_weight1);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight2",&W_weight2);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight3",&W_weight3); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight4",&W_weight4);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight5",&W_weight5); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight6",&W_weight6);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight7",&W_weight7); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight8",&W_weight8);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_MuonIDSF",&W_MuonIDSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_MuonIsoSF",&W_MuonIsoSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_MuonTrigSF",&W_MuonTrigSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_MuonTrigSF_Runs273158to274093",&W_MuonTrigSF_Runs273158to274093);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_MuonTrigSF_Runs274094to276097",&W_MuonTrigSF_Runs274094to276097);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_ElectronIDSF",&W_ElectronIDSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_ElectronRecoSF",&W_ElectronRecoSF);  
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_TopPtReweighing",&W_TopPtReweighing);  

        ttree[(dataSetName).c_str()]->SetBranchAddress("I_run_num",&run_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_evt_num",&evt_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_lumi_num",&lumi_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_npu",&npu);


        // electrons
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&pt_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&phi_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&eta_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_superCluster_electron",&eta_superCluster_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&E_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("chargedHadronIso_electron",&chargedHadronIso_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("neutralHadronIso_electron",&neutralHadronIso_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("photonIso_electron",&photonIso_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pfIso_electron",&pfIso_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&charge_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("d0_electron",&d0_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("d0BeamSpot_electron",&d0BeamSpot_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("sigmaIEtaIEta_electron",&sigmaIEtaIEta_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("deltaEtaIn_electron",&deltaEtaIn_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("deltaPhiIn_electron",&deltaPhiIn_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("hadronicOverEm_electron",&hadronicOverEm_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_missingHits_electron",&missingHits_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_passConversion_electron",&passConversion_electron);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_isEBEEGap",&isEBEEGap);
      
        // muons
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&pt_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&phi_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&eta_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&E_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("chargedHadronIso_muon",&chargedHadronIso_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("neutralHadronIso_muon",&neutralHadronIso_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("photonIso_muon",&photonIso_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pfIso_muon",&pfIso_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&charge_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("d0_muon",&d0_muon);
        ttree[(dataSetName).c_str()]->SetBranchAddress("d0BeamSpot_muon",&d0BeamSpot_muon);

        //SelectedLepton
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_lepton",&pt_lepton);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_lepton",&phi_lepton);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_lepton",&eta_lepton);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_lepton",&E_lepton);
        
        // jets
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&nJets);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVL",&nJets_CSVL);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&nJets_CSVM);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVT",&nJets_CSVT);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_cMVAL",&nJets_cMVAL);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_cMVAM",&nJets_cMVAM);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_cMVAT",&nJets_cMVAT);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_jet",&charge_jet);	    
        ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&incl_charge_jet);	    
        ttree[(dataSetName).c_str()]->SetBranchAddress("CSVv2",&CSVv2);
        ttree[(dataSetName).c_str()]->SetBranchAddress("cMVA",&cMVA);
        ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&cdiscCvsL_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&cdiscCvsB_jet);
        ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&jet_matchedMC_pdgID);
        ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&jet_matchedMC_motherpdgID);
        ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&jet_matchedMC_grannypdgID);
       
        // met 
        ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px", &met_Px); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py", &met_Py); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("met_Pt", &met_Pt); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("met_Eta", &met_Eta); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("met_Phi", &met_Phi); 

        // Jet-indices associated to the jet-assignment in the bMVA method
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_LepTop",&TOPTOPLEPHAD_JetIdx_LepTop);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_HadTop",&TOPTOPLEPHAD_JetIdx_HadTop);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_W1",&TOPTOPLEPHAD_JetIdx_W1);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_W2",&TOPTOPLEPHAD_JetIdx_W2);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_LepTop",&TOPTOPLEPHBB_JetIdx_LepTop);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_HadTop",&TOPTOPLEPHBB_JetIdx_HadTop);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_H1",&TOPTOPLEPHBB_JetIdx_H1);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_H2",&TOPTOPLEPHBB_JetIdx_H2);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_LepTop_hut",&TOPHLEPBB_JetIdx_LepTop_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_H1_hut",&TOPHLEPBB_JetIdx_H1_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_H2_hut",&TOPHLEPBB_JetIdx_H2_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_LepTop_hct",&TOPHLEPBB_JetIdx_LepTop_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_H1_hct",&TOPHLEPBB_JetIdx_H1_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPHLEPBB_JetIdx_H2_hct",&TOPHLEPBB_JetIdx_H2_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hut",&MVA_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hct",&MVA_TOPHLEPBB_hct);
                  
        //MC variables (affected by TopPtReweighing
        ttree[(dataSetName).c_str()]->SetBranchAddress("MC_TopPt",&MC_TopPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MC_AntiTopPt",&MC_AntiTopPt);

        double nloSF = 1.;
        int nPos = 0; 
        int nNeg = 0;
        if(isAMC && !isData)
        {
            for (int k = 0; k<nEntries; k++)
            {
                ttree[dataSetName.c_str()]->GetEntry(k);
                if( W_nloWeight > 0) nPos++;
                else if( W_nloWeight < 0) nNeg ++;
            }
            nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
        }		

        Double_t average_TopPtWeight = 0;
        if(dataSetName.find("TTJets") != std::string::npos)
        {
            for (int k = 0; k<nEntries; k++)
            {
                ttree[dataSetName.c_str()]->GetEntry(k);
                average_TopPtWeight = average_TopPtWeight + W_TopPtReweighing;
            }
            average_TopPtWeight = average_TopPtWeight/nEntries;
        }
		
  	    //***********************************************RUNNING OVER EVENTS**********************************************
		    for (int j = 0; j<nEntries; j++)
		    {
		              
            if(debug)
            {
                if(!isData)cin.get();
                cout << " " << endl;
                cout << "------------NEW EVENT: " << j << " --------------" << endl;
            }
			      ttree[dataSetName.c_str()]->GetEntry(j);
		        if(!doInclusive)
		        {
		            if(nJets_CSVM != baseline_bjets)  continue;

		            if(baseline_jets == 3 && nJets != baseline_jets) continue;
		            else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		        }

            float ScaleFactor = 1.; // event scale factor
            //Safety triggers in case there are strange things happening in the event weights
            if(!PVreweighing)
            {
                if(W_puSF <= 0|| W_fleptonSF <= 0 || W_btagWeight_CSVv2M_mujets_central <= 0 || nloSF <= 0 || Luminosity <= 0 )
                {
                      cout << "----- Event " << j << " has a negative weight. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << "; Luminosity=" << Luminosity << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_puSF != W_puSF|| W_fleptonSF != W_fleptonSF || W_btagWeight_CSVv2M_mujets_central != W_btagWeight_CSVv2M_mujets_central || nloSF != nloSF)
                {
                      cout << "----- Event " << j << " has a Nan weight. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_puSF >= 20|| W_fleptonSF >= 20 || W_btagWeight_CSVv2M_mujets_central >= 20 || nloSF >= 20)
                {
                      cout << "----- Event " << j << " has a weight larger than 20. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
            }
            else
            {
                if(W_fleptonSF <= 0 || W_btagWeight_CSVv2M_mujets_central <= 0 || nloSF <= 0 || Luminosity <= 0 )
                {
                      cout << "----- Event " << j << " has a negative weight. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << "; Luminosity=" << Luminosity << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_fleptonSF != W_fleptonSF || W_btagWeight_CSVv2M_mujets_central != W_btagWeight_CSVv2M_mujets_central || nloSF != nloSF)
                {
                      cout << "----- Event " << j << " has a Nan weight. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_fleptonSF >= 20 || W_btagWeight_CSVv2M_mujets_central >= 20 || nloSF >= 20)
                {
                      cout << "----- Event " << j << " has a weight larger than 20. Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
			      }
			      if(!isData)
			      {
			          
                double W_puSF_applied = 1.;
			          if(!PVreweighing) W_puSF_applied = W_puSF;
			          else
			          {
			              W_puSF_applied = W_nPV.ITweight( (int)nvtx );
			          }
			          ScaleFactor = ScaleFactor * W_puSF_applied;
			          ScaleFactor = ScaleFactor * W_fleptonSF;
			          ScaleFactor = ScaleFactor * W_btagWeight_CSVv2M_mujets_central;
			          ScaleFactor = ScaleFactor * nloSF;
//                ScaleFactor = ScaleFactor * W_btagWeight_shape;
                if(dataSetName.find("TTJets") != std::string::npos) ScaleFactor = ScaleFactor * W_TopPtReweighing/average_TopPtWeight;
			      }
			      if(debug && !isData)
			      {
                cout << "----- Event " << j << "Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_CSVv2M_mujets_central=" << W_btagWeight_CSVv2M_mujets_central << "; nloSF=" << nloSF << endl;
                cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
			          cout << "   SCALE FACTOR is: " << ScaleFactor << endl;
			      }
			      
			      //**********************EVENT RECONSTRUCTIONS*******************************
			      TLorentzVector Lepton;
			      TLorentzVector Hb1_TOPHLEPBB_hut, Hb2_TOPHLEPBB_hut, LepTopB_TOPHLEPBB_hut;
			      TLorentzVector Hb1_TOPHLEPBB_hct, Hb2_TOPHLEPBB_hct, LepTopB_TOPHLEPBB_hct;
			      TLorentzVector Hb1_TOPTOPLEPHBB, Hb2_TOPTOPLEPHBB, LepTopB_TOPTOPLEPHBB, HadTopB_TOPTOPLEPHBB;
			      TLorentzVector Hb1_TOPTOPLEPHAD, Hb2_TOPTOPLEPHAD, LepTopB_TOPTOPLEPHAD, HadTopB_TOPTOPLEPHAD;
			      TLorentzVector Higgs_TOPHLEPBB_hut, Higgs_TOPHLEPBB_hct, Higgs_TOPTOPLEPHBB, Higgs_TOPTOPLEPHAD;
			      TLorentzVector HadTop_TOPTOPLEPHBB, HadTop_TOPTOPLEPHAD;
			          
			      Lepton.SetPtEtaPhiE(pt_lepton, eta_lepton, phi_lepton, E_lepton);
			          
			      Hb1_TOPHLEPBB_hct.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_H1_hct],eta_jet[TOPHLEPBB_JetIdx_H1_hct],phi_jet[TOPHLEPBB_JetIdx_H1_hct],E_jet[TOPHLEPBB_JetIdx_H1_hct]);
			      Hb2_TOPHLEPBB_hct.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_H2_hct],eta_jet[TOPHLEPBB_JetIdx_H2_hct],phi_jet[TOPHLEPBB_JetIdx_H2_hct],E_jet[TOPHLEPBB_JetIdx_H2_hct]);
			      LepTopB_TOPHLEPBB_hct.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_LepTop_hct],eta_jet[TOPHLEPBB_JetIdx_LepTop_hct],phi_jet[TOPHLEPBB_JetIdx_LepTop_hct],E_jet[TOPHLEPBB_JetIdx_LepTop_hct]);
            Higgs_TOPHLEPBB_hct = Hb1_TOPHLEPBB_hct + Hb2_TOPHLEPBB_hct;

			      Hb1_TOPHLEPBB_hut.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_H1_hut],eta_jet[TOPHLEPBB_JetIdx_H1_hut],phi_jet[TOPHLEPBB_JetIdx_H1_hut],E_jet[TOPHLEPBB_JetIdx_H1_hut]);
			      Hb2_TOPHLEPBB_hut.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_H2_hut],eta_jet[TOPHLEPBB_JetIdx_H2_hut],phi_jet[TOPHLEPBB_JetIdx_H2_hut],E_jet[TOPHLEPBB_JetIdx_H2_hut]);
			      LepTopB_TOPHLEPBB_hut.SetPtEtaPhiE(pt_jet[TOPHLEPBB_JetIdx_LepTop_hut],eta_jet[TOPHLEPBB_JetIdx_LepTop_hut],phi_jet[TOPHLEPBB_JetIdx_LepTop_hut],E_jet[TOPHLEPBB_JetIdx_LepTop_hut]);
            Higgs_TOPHLEPBB_hut = Hb1_TOPHLEPBB_hut + Hb2_TOPHLEPBB_hut;

			      Hb1_TOPTOPLEPHBB.SetPtEtaPhiE(pt_jet[TOPTOPLEPHBB_JetIdx_H1],eta_jet[TOPTOPLEPHBB_JetIdx_H1],phi_jet[TOPTOPLEPHBB_JetIdx_H1],E_jet[TOPTOPLEPHBB_JetIdx_H1]);
			      Hb2_TOPTOPLEPHBB.SetPtEtaPhiE(pt_jet[TOPTOPLEPHBB_JetIdx_H2],eta_jet[TOPTOPLEPHBB_JetIdx_H2],phi_jet[TOPTOPLEPHBB_JetIdx_H2],E_jet[TOPTOPLEPHBB_JetIdx_H2]);
			      LepTopB_TOPTOPLEPHBB.SetPtEtaPhiE(pt_jet[TOPTOPLEPHBB_JetIdx_LepTop],eta_jet[TOPTOPLEPHBB_JetIdx_LepTop],phi_jet[TOPTOPLEPHBB_JetIdx_LepTop],E_jet[TOPTOPLEPHBB_JetIdx_LepTop]);
			      HadTopB_TOPTOPLEPHBB.SetPtEtaPhiE(pt_jet[TOPTOPLEPHBB_JetIdx_HadTop],eta_jet[TOPTOPLEPHBB_JetIdx_HadTop],phi_jet[TOPTOPLEPHBB_JetIdx_HadTop],E_jet[TOPTOPLEPHBB_JetIdx_HadTop]);
            Higgs_TOPTOPLEPHBB = Hb1_TOPTOPLEPHBB + Hb2_TOPTOPLEPHBB;
            HadTop_TOPTOPLEPHBB  = Higgs_TOPTOPLEPHBB + HadTopB_TOPTOPLEPHBB;

			      Hb1_TOPTOPLEPHAD.SetPtEtaPhiE(pt_jet[TOPTOPLEPHAD_JetIdx_W1],eta_jet[TOPTOPLEPHAD_JetIdx_W1],phi_jet[TOPTOPLEPHAD_JetIdx_W1],E_jet[TOPTOPLEPHAD_JetIdx_W1]);
			      Hb2_TOPTOPLEPHAD.SetPtEtaPhiE(pt_jet[TOPTOPLEPHAD_JetIdx_W2],eta_jet[TOPTOPLEPHAD_JetIdx_W2],phi_jet[TOPTOPLEPHAD_JetIdx_W2],E_jet[TOPTOPLEPHAD_JetIdx_W2]);
			      LepTopB_TOPTOPLEPHAD.SetPtEtaPhiE(pt_jet[TOPTOPLEPHAD_JetIdx_LepTop],eta_jet[TOPTOPLEPHAD_JetIdx_LepTop],phi_jet[TOPTOPLEPHAD_JetIdx_LepTop],E_jet[TOPTOPLEPHAD_JetIdx_LepTop]);
			      HadTopB_TOPTOPLEPHAD.SetPtEtaPhiE(pt_jet[TOPTOPLEPHAD_JetIdx_HadTop],eta_jet[TOPTOPLEPHAD_JetIdx_HadTop],phi_jet[TOPTOPLEPHAD_JetIdx_HadTop],E_jet[TOPTOPLEPHAD_JetIdx_HadTop]);
            Higgs_TOPTOPLEPHAD = Hb1_TOPTOPLEPHAD + Hb2_TOPTOPLEPHAD;
            HadTop_TOPTOPLEPHAD  = Higgs_TOPTOPLEPHAD + HadTopB_TOPTOPLEPHAD;
/*
if(!isData)
{
    cout << "----- TOPTOPLEPHAD -----" << endl;
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

    cout << "HadTopMass (TOPTOPLEPHAD): " << HadTop_TOPTOPLEPHAD.M() << endl;
    cout << "HiggsMass (TOPTOPLEPHAD): " << Higgs_TOPTOPLEPHAD.M() << endl;
}
*/
  	        //***********************************************FILLING PLOTS**********************************************
				    MSPlot["NCSVv2Ljets"]->Fill(nJets_CSVL, datasets[d], true, Luminosity * ScaleFactor);
				    MSPlot["NCSVv2Mjets"]->Fill(nJets_CSVM, datasets[d], true, Luminosity * ScaleFactor);
				    MSPlot["NCSVv2Tjets"]->Fill(nJets_CSVT, datasets[d], true, Luminosity * ScaleFactor);
				    MSPlot["Njets"]->Fill(nJets, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["LeptonPt"]->Fill(pt_lepton, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["LeptonEta"]->Fill(eta_lepton, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["LeptonPhi"]->Fill(phi_lepton, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["NPV"]->Fill(nvtx, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["NPU"]->Fill(npu, datasets[d], true, Luminosity * ScaleFactor);
            for(int i_Jet = 0; i_Jet < nJets; i_Jet++)
            {
                MSPlot["JetPt"]->Fill(pt_jet[i_Jet], datasets[d], true, Luminosity * ScaleFactor);
                MSPlot["JetEta"]->Fill(eta_jet[i_Jet], datasets[d], true, Luminosity * ScaleFactor);
                MSPlot["JetPhi"]->Fill(phi_jet[i_Jet], datasets[d], true, Luminosity * ScaleFactor);
                MSPlot["JetCSVv2"]->Fill(CSVv2[i_Jet], datasets[d], true, Luminosity * ScaleFactor);
                MSPlot["JetcMVAv2"]->Fill(cMVA[i_Jet], datasets[d], true, Luminosity * ScaleFactor);
            }			                
            MSPlot["Hmass_TOPTOPLEPHAD"]->Fill(Higgs_TOPTOPLEPHAD.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["Hmass_TOPTOPLEPHBB"]->Fill(Higgs_TOPTOPLEPHBB.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["Hmass_TOPHLEPBB_hct"]->Fill(Higgs_TOPHLEPBB_hct.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["Hmass_TOPHLEPBB_hut"]->Fill(Higgs_TOPHLEPBB_hut.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["HadTopMass_TOPTOPLEPHAD"]->Fill(HadTop_TOPTOPLEPHAD.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["HadTopMass_TOPTOPLEPHBB"]->Fill(HadTop_TOPTOPLEPHBB.M(), datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["MC_TopPt"]->Fill(MC_TopPt, datasets[d], true, Luminosity * ScaleFactor);
            MSPlot["MC_AntiTopPt"]->Fill(MC_AntiTopPt, datasets[d], true, Luminosity * ScaleFactor);
			                
		  }//for-loop events
		              
    }//for-loop datasets
               



  
  cout << "MSPlot size: " << MSPlot.size() << endl;      




  string pathPNG = "MSPlots/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MSPlots";
  pathPNG += channel;
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/";
  pathPNG += date;
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/";
  pathPNG += category;
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  outfile->cd();


  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
     	string name = it->first;
     	MultiSamplePlot *temp = it->second;
      if (debug)
      {
          cout << "Saving the MSP" << endl;
          cout << " and it->first is " << it->first << endl;
          cout << " Luminosity is " << Luminosity << endl;
      }
//      temp->setDataLumi(Luminosity);
      	//    temp->Draw(name,RatioType, addRatioErrorBand, addErrorBand, ErrorBandAroundTotalInput, scaleNPSignal);              
      	//MultiSamplePlot options
          /*
          bool showEntriesLegend = false; //to show number of (weighted) events of the samples in the legend
          bool setCMSPrelim = false; //if true, will display "CMS Preliminary", otherwise "CMS"
          int RatioType = 0; //0: no ratio plot, 1: ratio = data/MC, 2: ratio = (data-MC)/MC
          bool addErrorBand = false; //display an error band around the stacked SM MC on the main canvas
          bool addRatioErrorBand = false; //display an error band on the ratio plot below the main canvas
          bool ErrorBandAroundTotalInput = false; //see dedicated discussion below.
          string errorbandfile = "ErrorBands/ErrorBandFile_15Jul15.root";  //a root file containing systematically shifted distributions to create error bands around the stacked SM MC. See dedicated discussion below.
          bool dosystfile = false; //see dedicated discussion below.
          int scaleNPSignal = 20; //determines the factor with which the new physics signal samples are scaled, only on the canvas (note that the TH1F histogram in the MSPlot output root file itself is not scaled with this factor!)
          bool savePNG = false; //automatically save png files of MSPlots.
          */
      cout << "Drawing MSP: " << it->first << endl;
      temp->Draw("MyMSP_"+it->first, 1, false, false, false, 1);
      temp->Write(outfile, it->first, true,pathPNG, "png");
	}

  	outfile->Write("kOverwrite");

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

}

// function that converts an int into a string
std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}

void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug)
{
    bool doInclusive = false;
    string category;
    if(baseline_bjets == 0 && baseline_jets == 0)
    {
        doInclusive = true;
        category = "Inclusive";
    }
    else
    {
        category = "b"+intToStr(baseline_bjets)+"j"+intToStr(baseline_jets);
    }    

    cout << ".. ..Making nPV_unw distributions for all samples.. .." << endl;


    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    TString TreePath = "Merged/Ntuples" + channel + "/Ntuples" + date;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
    //***************************************************LOADING DATASETS****************************************************
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  	string dataSetName, filepath;

    //***************************************************CREATING PLOT****************************************************
    //Format of MSPlots: MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float Min, float Max, string XaxisLabel, string YaxisLabel, string Text, string Units)

    MSPlot_nPV["NPV_unw"] = new MultiSamplePlot(datasets, "NPV_unw", 51, -0.5, 50.5, "Number of PV","Events", category); 

  
 

  	//***********************************************RUNNING OVER DATASETS**********************************************
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	  {
		    dataSetName = datasets[d]->Name();
		    cout<<".. ..Dataset:  :"<<dataSetName<<endl;
		    filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		                

  	    //***********************************************IMPORTING VARIABLES**********************************************
		    string TTreename = "ObjectVarsTree";	
		    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset

        int nEntries;

		    nEntries = ttree[dataSetName.c_str()]->GetEntries();
		    cout<<"                 nEntries: "<<nEntries<<endl;
        Int_t nvtx;
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
		
  	    //***********************************************RUNNING OVER EVENTS**********************************************
		    for (int j = 0; j<nEntries; j++)
		    {
			      ttree[dataSetName.c_str()]->GetEntry(j);
            MSPlot_nPV["NPV_unw"]->Fill(nvtx, datasets[d], true, 1.);
			                
		  }//for-loop events
		              
    }//for-loop datasets
  
  cout << "MSPlot size: " << MSPlot_nPV.size() << endl;      

  string pathPNG = "MSPlots/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MSPlots";
  pathPNG += channel;
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/";
  pathPNG += date;
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/";
  pathPNG += category;
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  TFile *outfile = new TFile((pathPNG+"/Output_NPV.root").c_str(),"recreate");
  outfile->cd();


  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot_nPV.begin(); it != MSPlot_nPV.end(); it++)
  {
     	string name = it->first;
     	MultiSamplePlot *temp = it->second;
      cout << "Drawing MSP: " << it->first << endl;
      temp->Draw("MyMSP_"+it->first, 1, false, false, false, 1);
      temp->Write(outfile, it->first, true,pathPNG, "png");
	}

  	outfile->Write("kOverwrite");
    cout << ".. .. Finished making nPV_unw distributions for all samples.. .." << endl;

}





