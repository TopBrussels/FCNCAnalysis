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

//includes for MVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

//includes for Kinematic fitting
#include "FCNCAnalysis_76X/TopKinFit/kinfit.h"


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


float Luminosity = 0.; // pb-1 Muon  = 2628.727204156, Electron = 2094.087
std::string channel = "_El";
std::string date = "_26_7_2016";
int maxNumbObjToPlot = 5;
Bool_t debug = false;
bool train_mva = true;
float workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Medium = 0.800;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Tight = 0.935;//working points updated to 2016 BTV-POG recommendations.

std::string xmlNom = "config/FullMcBkgdSamplesV8_TreeProcessor.xml";
TString TreePath = "Merged/Ntuples" + channel + "/Ntuples" + date;


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string NTupleName);
void MVA_JetCombTraining(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP);
void MVA_JetCombComputer(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP);
void MVA_Analysis(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP, bool doEventTraining);
void MSPCreator ();



int main(int argc, char *argv[])
{

    string TopKinFitMethod              = argv[1]; //Valid options are TThypo, SThypo and SMttHypo
    int baseline_bjets             = strtol(argv[2], NULL,10);
    int baseline_nonbjets                 = strtol(argv[3], NULL,10);
    string btag_WP                    = argv[4]; // Valid options are Loose, Medium and Tight
    string SignalSample                    = argv[5]; // Valid options are ST_tcH, ST_tuH, TT_tcH, TT_tuH
    bool doJetCombTraining                    = strtol(argv[6],NULL,10);
    bool doJetCombComputer                    = strtol(argv[7],NULL,10);
    bool doTraining                    = strtol(argv[8],NULL,10);

    cout << "--------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << "MACRO command line arguments: " << TopKinFitMethod << " " << baseline_bjets << " " << baseline_nonbjets << " " << btag_WP << " " << SignalSample << " " << doJetCombTraining << " " << doJetCombComputer << " " << endl;
    cout << "--------------------------------------------------------------------------" << endl;


    float WP_forBaseline = workingpointvalue_Medium;
    if(btag_WP == "Loose") WP_forBaseline = workingpointvalue_Loose;
    else if(btag_WP == "Tight") WP_forBaseline = workingpointvalue_Tight;

    string SignalSampleTraining;
    if(SignalSample == "ST_tcH") SignalSampleTraining = "NP_overlay_ST_tHToBB_1L_Kappa_hct"; 
    else if(SignalSample == "ST_tuH") SignalSampleTraining = "NP_overlay_ST_tHToBB_1L_Kappa_hut"; 
    else if(SignalSample == "TT_tcH") SignalSampleTraining = "NP_overlay_TTtoTHToBB-1L-Kappa-hct"; 
    else if(SignalSample == "TT_tuH") SignalSampleTraining = "NP_overlay_TTtoTHToBB-1L-Kappa-hut"; 

    clock_t start = clock();

    // calling datasetPlotter to create MSPplots

/*    DatasetPlotter(11, -0.5, 10.5, "I_nJets", "ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVL", "ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVM", "ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVT", "ObjectVarsTree");
    DatasetPlotter(40, 0, 400, "pt_muon", "ObjectVarsTree");
    DatasetPlotter(40, 0, 400, "pt_electron", "ObjectVarsTree");
    DatasetPlotter(35, -0.5, 34.5, "I_nvtx", "ObjectVarsTree");
    DatasetPlotter(70, 0, 700, "pt_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(50, -3.15, 3.15, "eta_Jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(30, -3.15, 3.15, "phi_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(21, -10.5, 10.5, "charge_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(25, 0, 1, "bdisc_Jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_pdgID[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_motherpdgID[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_grannypdgID[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(50,-1, 1, "incl_charge_jet[I_nJets]", "ObjectVarsTree");
    DatasetPlotter(40, 0., 500., "MTlepmet", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MLepTop_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW_GenMatch", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 5., "dR_lepJet_min", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW", "AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "Mbb", "AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW", "AdvancedVarsTree");
*/

  int nToys = 500;
  std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kfit_SMttHypo = new KINFIT::kfit();
  KINFIT::kfit *kfit_TThypo = new KINFIT::kfit();
  KINFIT::kfit *kfit_SThypo = new KINFIT::kfit();

  kfit_SMttHypo->Init(TOPTOPLEPHAD);
  kfit_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
  kfit_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
  kfit_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
  kfit_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
  kfit_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
  kfit_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
  kfit_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
  kfit_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
  kfit_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
  kfit_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
  kfit_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
  kfit_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
  kfit_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
  kfit_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
  kfit_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
  kfit_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
  kfit_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
  kfit_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
  kfit_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
  kfit_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
  kfit_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
  kfit_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
  kfit_SMttHypo->SetNToy(nToys);

  kfit_TThypo->Init(TOPTOPLEPHBB);
  kfit_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
  kfit_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
  kfit_TThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
  kfit_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
  kfit_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
  kfit_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
  kfit_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
  kfit_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
  kfit_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
  kfit_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
  kfit_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
  kfit_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
  kfit_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
  kfit_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
  kfit_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
  kfit_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
  kfit_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
  kfit_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
  kfit_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
  kfit_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
  kfit_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
  kfit_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
  kfit_TThypo->SetNToy(nToys);

  kfit_SThypo->Init(TOPHLEPBB);
  kfit_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
  kfit_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
  kfit_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
  kfit_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
  kfit_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
  kfit_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
  kfit_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
  kfit_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
  kfit_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
  kfit_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
  kfit_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
  kfit_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
  kfit_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
  kfit_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
  kfit_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
  kfit_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
  kfit_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
  kfit_SThypo->SetNToy(nToys);

  vector<std::string> MVAvars_TThypo;
  vector<std::string> MVAvars_SThypo;

  MVAvars_TThypo.push_back("Hmass");
  MVAvars_TThypo.push_back("LepTopmass");
  MVAvars_TThypo.push_back("HadTopmass");
  MVAvars_TThypo.push_back("DR_H_HadTop");
  MVAvars_TThypo.push_back("DR_H_LepTop");
  MVAvars_TThypo.push_back("LepTopPt");
  MVAvars_TThypo.push_back("HadTopPt");

  MVAvars_SThypo.push_back("Hmass");
  MVAvars_SThypo.push_back("LepTopmass");
  MVAvars_SThypo.push_back("DR_H_LepTop");
  MVAvars_SThypo.push_back("LepTopPt");




  if(TopKinFitMethod == "SThypo")
  {
      if(doJetCombTraining) MVA_JetCombTraining(4, SignalSampleTraining, TopKinFitMethod, kfit_SThypo, MVAvars_SThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doJetCombComputer) MVA_JetCombComputer(4, SignalSampleTraining, TopKinFitMethod, kfit_SThypo, MVAvars_SThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doTraining) MVA_Analysis(4, SignalSampleTraining, TopKinFitMethod, kfit_SThypo, MVAvars_SThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline,doTraining);
  }
  else if(TopKinFitMethod == "TThypo")
  {
      if(doJetCombTraining) MVA_JetCombTraining(4, SignalSampleTraining, TopKinFitMethod, kfit_TThypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doJetCombComputer) MVA_JetCombComputer(4, SignalSampleTraining, TopKinFitMethod, kfit_TThypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doTraining) MVA_Analysis(4, SignalSampleTraining, TopKinFitMethod, kfit_TThypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline,doTraining);
  }
  else if(TopKinFitMethod == "SMttHypo")
  {
      if(doJetCombTraining) MVA_JetCombTraining(4, SignalSampleTraining, TopKinFitMethod, kfit_SMttHypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doJetCombComputer) MVA_JetCombComputer(4, SignalSampleTraining, TopKinFitMethod, kfit_SMttHypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline);
      if(doTraining) MVA_Analysis(4, SignalSampleTraining, TopKinFitMethod, kfit_SMttHypo, MVAvars_TThypo, baseline_bjets, baseline_nonbjets, WP_forBaseline,doTraining);
  }

	MSPCreator ();

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

}




void MVA_JetCombTraining(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP) 
{

  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);
  string pathMVA_ = "MVA/TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);

  std::string TrainMethod = KinFitMethod+"_"+intToStr(MVAvars.size())+"Vars"+"_"+intToStr(nBaselinebjets_)+"B"+intToStr(nBaselinenonbjets_)+"J"+channel;

  MVATrainer* Eventtrainer_PartReco = 0;
  MVATrainer* Eventtrainer_FullReco = 0;
  Eventtrainer_PartReco = new MVATrainer("BDT","TrainedJetCombMVA_PartReco_"+TrainMethod, pathMVA_+"TrainedJetCombMVA_PartReco_"+TrainMethod+"_"+SignalName+".root");
  Eventtrainer_FullReco = new MVATrainer("BDT","TrainedJetCombMVA_FullReco_"+TrainMethod, pathMVA_+"TrainedJetCombMVA_FullReco_"+TrainMethod+"_"+SignalName+".root");

  
  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_PartReco->bookInputVar(MVAvars[N_var]);
          	Eventtrainer_FullReco->bookInputVar(MVAvars[N_var]);
      }
  }
  if(KinFitMethod == "SThypo")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_PartReco->bookInputVar(MVAvars[N_var]);
          	Eventtrainer_FullReco->bookInputVar(MVAvars[N_var]);
      }
  }

  ////////////////////////////////////////////////////////////
  // Load Datasets
  //////////////////////////////////////////////////////////////////////
 	const char *xmlfile = xmlNom.c_str();
 	cout << "Using config file: " << xmlfile << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
	
		dataSetName = datasets[d]->Name();
    if(SignalName != dataSetName) continue;

    bool SingleTop = false;
    if(dataSetName.find("NP_overlay_ST")!=string::npos) SingleTop = true;

		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "ObjectVarsTree";
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
    /////////////////////////////////////////
    // Get object variables
    ////////////////////////////////////////
    int NumberOfJets;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
    
    double InclJetCharge[20];
    double bDiscJet[20];
    double CvsBJet[20];
    double CvsLJet[20];
    double pdgID[20];
    double MotherpdgID[20];
    double GrandMotherpdgID[20];
    double lepCharge;
    double lepPt;
    double lepEta;
    double lepPhi;
    double lepE;

    double pt_jet[20];
    double phi_jet[20];
    double eta_jet[20];
    double E_jet[20];

    double met_px;
    double met_py;


    ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
    if(channel == "_Mu")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
    }
    if(channel == "_El")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
    }
		  


		//////////////////////////////////////////////////////////
		// Running on events
		//////////////////////////////////////////////////////////
		for (int j = 0; j<int(nEntries/skipEvents); j++)
		{
	        std::vector<float> BJetPt;
	        std::vector<float> BJetEta;
	        std::vector<float> BJetPhi;
	        std::vector<float> BJetE;

	        std::vector<float> NonBJetPt;
	        std::vector<float> NonBJetEta;
	        std::vector<float> NonBJetPhi;
	        std::vector<float> NonBJetE;

	        std::vector<float> LeptonPt;
	        std::vector<float> LeptonEta;
	        std::vector<float> LeptonPhi;
	        std::vector<float> LeptonE;
	        
	        vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	        vector <int> MapIndex_NonBindex;
	        
			    ttree[(dataSetName).c_str()]->GetEntry(j);
          for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
          {
          
 
                if(bDiscJet[i_Jet]  > btagWP)
                {
                    BJetPt.push_back(pt_jet[i_Jet]);
                    BJetEta.push_back(eta_jet[i_Jet]);
                    BJetPhi.push_back(phi_jet[i_Jet]);
                    BJetE.push_back(E_jet[i_Jet]);
                    
                    MapIndex_Bindex.push_back(i_Jet);
                }
                else
                {
                    NonBJetPt.push_back(pt_jet[i_Jet]);
                    NonBJetEta.push_back(eta_jet[i_Jet]);
                    NonBJetPhi.push_back(phi_jet[i_Jet]);
                    NonBJetE.push_back(E_jet[i_Jet]);

                    MapIndex_NonBindex.push_back(i_Jet);
                }
          }
          
          LeptonPt.push_back(lepPt);
          LeptonEta.push_back(lepEta);
          LeptonPhi.push_back(lepPhi);
          LeptonE.push_back(lepE);



          kfit->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
          kfit->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
          if(channel == "_El") kfit->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          else if(channel == "_Mu") kfit->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          kfit->SetMet(met_px,met_py);

          if(BJetPt.size()>=nBaselinebjets_ && NonBJetPt.size()>=nBaselinenonbjets_)// BASELINE SELECTION
          {
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              if(KinFitMethod == "SMttHypo")
              {
                  if(BJetPt.size()>=2 && NonBJetPt.size()>=2)//Additional selection to baseline selection -- This is the selection defined according to the kinfit method
                  {
                      kfit->Run();
                      int nPerm_SMttHypo = kfit->GetNPerm();


                      for(int ip=0;ip<nPerm_SMttHypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);

                          int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kfit->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kfit->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kfit->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];


                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJETHAD_, NONBJET1_, NONBJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;

                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          BJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETHAD],eta_jet[IndexAllJetColl_BJETHAD],phi_jet[IndexAllJetColl_BJETHAD],E_jet[IndexAllJetColl_BJETHAD]);
                          NONBJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET1],eta_jet[IndexAllJetColl_NONBJET1],phi_jet[IndexAllJetColl_NONBJET1],E_jet[IndexAllJetColl_NONBJET1]);
                          NONBJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET2],eta_jet[IndexAllJetColl_NONBJET2],phi_jet[IndexAllJetColl_NONBJET2],E_jet[IndexAllJetColl_NONBJET2]);

                          Higgs_ = NONBJET1_+NONBJET2_;
                          HadTop_ = Higgs_+BJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;

                          double Hmass = -999.;
                          double LepTopmass = -999.;
                          double HadTopmass = -999.;
                          double DR_H_HadTop = -999.;
                          double DR_H_LepTop = -999.;
                          double LepTopPt = -999.;
                          double HadTopPt = -999.;
                          if(chi2>10E+9) //No neutrino reconstructed
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              HadTopmass = HadTop_.M();
                              DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              DR_H_LepTop = (LEP_+BJETLEP_).DeltaPhi(Nu_);
                              LepTopPt = LepTop_.Pt();
                              HadTopPt = HadTop_.Pt();
                          }
                          else
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = LepTop_.M();
                              HadTopmass = HadTop_.M();
                              DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              LepTopPt = LepTop_.Pt();
                              HadTopPt = HadTop_.Pt();
                          }

                          if(SingleTop)
                          {
                              if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                              }
                              else
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                              }
                          }//SingleTop
                          else
                          {
                              if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                              }
                              else
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                              }
                          }// Not SingleTop
                      }//jet permutation loop
                  }//TOPTOPLEPHAD selection
              }//SMtthypo
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "TThypo")
              {
                  if(BJetPt.size()>=3 && NonBJetPt.size()>=1)//Additional selection to baseline selection -- This is the selection defined according to the kinfit method
                  {
                      kfit->Run();
                      int nPerm_TThypo = kfit->GetNPerm();


                      for(int ip=0;ip<nPerm_TThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
                          if(chi2<10E+9)  continue;

		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kfit->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, NONBJETHAD_, BJET1_, BJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;

                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          HadTop_ = Higgs_+NONBJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double Hmass = -999.;
                          double LepTopmass = -999.;
                          double HadTopmass = -999.;
                          double DR_H_HadTop = -999.;
                          double DR_H_LepTop = -999.;
                          double LepTopPt = -999.;
                          double HadTopPt = -999.;
                          if(chi2>10E+9) //No neutrino reconstructed
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              HadTopmass = HadTop_.M();
                              DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              DR_H_LepTop = (LEP_+BJETLEP_).DeltaPhi(Nu_);
                              LepTopPt = LepTop_.Pt();
                              HadTopPt = HadTop_.Pt();
                          }
                          else
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = LepTop_.M();
                              HadTopmass = HadTop_.M();
                              DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              LepTopPt = LepTop_.Pt();
                              HadTopPt = HadTop_.Pt();
                          }

                          if(SingleTop)
                          {
                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                              }
                              else
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                              }
                          }//SingleTop
                          else
                          {
                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("S","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("S","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("S","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("S","HadTopPt", HadTopPt);
                                  }
                              }
                              else
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_PartReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_PartReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_PartReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("B","HadTopmass", HadTopmass);
                                      Eventtrainer_FullReco->Fill("B","DR_H_HadTop", DR_H_HadTop);
                                      Eventtrainer_FullReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("B","LepTopPt", LepTopPt);
                                      Eventtrainer_FullReco->Fill("B","HadTopPt", HadTopPt);
                                  }
                              }
                          }// Not SingleTop
                      }//jet permutation loop
                  }//TOPTOPLEPHBB selection
              }//TThypo
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "SThypo")
              {
                  if(BJetPt.size()>=3)//Additional selection to baseline selection -- This is the selection defined according to the kinfit method
                  {
                  

                      kfit->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                      kfit->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                      if(channel == "_El") kfit->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                      else if(channel == "_Mu") kfit->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                      kfit->SetMet(met_px,met_py);
                                
                                
                      kfit->Run();
                      int nPerm_SThypo = kfit->GetNPerm();
                                
                                
                      for(int ip=0;ip<nPerm_SThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
                          if(chi2<10E+9)  continue;
		                       
		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPHLEPBB,ip)];
                          
                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                          TLorentzVector Higgs_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double Hmass = -999.;
                          double LepTopmass = -999.;
                          double DR_H_LepTop = -999.;
                          double LepTopPt = -999.;
                          if(chi2>10E+9) //No neutrino reconstructed
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              LepTopPt = LepTop_.Pt();
                          }
                          else
                          {
                              Hmass = Higgs_.M();
                              LepTopmass = LepTop_.M();
                              DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              LepTopPt = LepTop_.Pt();
                          }

                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("S","LepTopPt", LepTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("S","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("S","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("S","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("S","LepTopPt", LepTopPt);
                                  }
                              }
                              else
                              {
                                  if(chi2>10E+9) //No neutrino reconstructed
                                  {
                                      Eventtrainer_PartReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_PartReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_PartReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_PartReco->Fill("B","LepTopPt", LepTopPt);
                                  }
                                  else
                                  {
                                      Eventtrainer_FullReco->Fill("B","Hmass", Hmass);
                                      Eventtrainer_FullReco->Fill("B","LepTopmass", LepTopmass);
                                      Eventtrainer_FullReco->Fill("B","DR_H_LepTop", DR_H_LepTop);
                                      Eventtrainer_FullReco->Fill("B","LepTopPt", LepTopPt);
                                  }
                              }
                      }//jet permutation loop
                  }//TOPHLEPBB selection
              }//SThypo
          }//Baseline selection
		}//for-loop events
  }//for-loop datasets


  Eventtrainer_PartReco->TrainMVA("Block","",0,0,"",0,0,"test",false);
      

  delete Eventtrainer_PartReco;
}





void MVA_JetCombComputer(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP)
{

  ///////////////////////////////////////////////////////////////////
  // Initializing MVA
  ///////////////////////////////////////////////////////////////////
  cout << "Initializing MVA for correct Jet combination" << endl;
  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);
  string pathMVA_ = "MVA/TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);

  std::string TrainMethod = KinFitMethod+"_"+intToStr(MVAvars.size())+"Vars"+"_"+intToStr(nBaselinebjets_)+"B"+intToStr(nBaselinenonbjets_)+"J"+channel;

  MVAComputer* Eventcomputer_FullReco_ =0;   
  MVAComputer* Eventcomputer_PartReco_ =0;   

  Eventcomputer_FullReco_ = new MVAComputer("BDT",pathMVA_+"TrainedJetCombMVA_FullReco_"+TrainMethod+"_"+SignalName+".root","TrainedJetCombMVA_PartReco_"+TrainMethod,MVAvars, "test");
  Eventcomputer_PartReco_ = new MVAComputer("BDT",pathMVA_+"TrainedJetCombMVA_PartReco_"+TrainMethod+"_"+SignalName+".root","TrainedJetCombMVA_PartReco_"+TrainMethod,MVAvars, "test");

  ////////////////////////////////////////////////////////////
  // Load Datasets
  //////////////////////////////////////////////////////////////////////
 	const char *xmlfile = xmlNom.c_str();
 	cout << "Using config file: " << xmlfile << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
	                  int NumberOfEvents = 0;
	                  
	                  //TopTopLepHad initializations
	                  int NSelectionPassedEvents_SMttHypo = 0;
	                  int NMCIdentifiedEvents_SMttHypo_STsignal = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_SMttHypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTbackground_TKF = 0;

	                  //TopTopLepHbb initializations
	                  int NSelectionPassedEvents_TThypo = 0;
	                  int NMCIdentifiedEvents_TThypo_STsignal = 0;
	                  int NMCIdentifiedEvents_TThypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_TThypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_TThypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_TThypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_TThypo_TTbackground_TKF = 0;

	                  //TopHLepbb initializations
	                  int NSelectionPassedEvents_SThypo = 0;
	                  int NMCIdentifiedEvents_SThypo_STsignal = 0;
	                  int NMCIdentifiedEvents_SThypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_SThypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_SThypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SThypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SThypo_TTbackground_TKF = 0;


	              
                    int nMCMatchedPassedEvents_SMttHypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_SMttHypo_STsignal = 0;
                    int nMCMatchedPassedEvents_SMttHypo_TTsignal = 0;
                    int nMCMatchedPassedEvents_SThypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_SThypo_STsignal = 0;
                    int nMCMatchedPassedEvents_SThypo_TTsignal = 0;
                    int nMCMatchedPassedEvents_TThypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_TThypo_STsignal = 0;
                    int nMCMatchedPassedEvents_TThypo_TTsignal = 0;

              
	
		dataSetName = datasets[d]->Name();

    bool SingleTop = false;
    if(dataSetName.find("NP_overlay_ST")!=string::npos) SingleTop = true;

		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "ObjectVarsTree";
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
    /////////////////////////////////////////
    // Get object variables
    ////////////////////////////////////////
    int NumberOfJets;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
    
    double InclJetCharge[20];
    double bDiscJet[20];
    double CvsBJet[20];
    double CvsLJet[20];
    double pdgID[20];
    double MotherpdgID[20];
    double GrandMotherpdgID[20];
    double lepCharge;
    double lepPt;
    double lepEta;
    double lepPhi;
    double lepE;

    double pt_jet[20];
    double phi_jet[20];
    double eta_jet[20];
    double E_jet[20];

    double met_px;
    double met_py;


    ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
    if(channel == "_Mu")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
    }
    if(channel == "_El")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
    }
		  


		//////////////////////////////////////////////////////////
		// Running on events
		//////////////////////////////////////////////////////////
		for (int j = int(nEntries/skipEvents)+1; j<2*int(nEntries/skipEvents); j++)
		{
		
		      NumberOfEvents++;
		
          int SMTT_Matched = 0;
          int Signal_Matched = 0;
	        
	        std::vector<float> BJetPt;
	        std::vector<float> BJetEta;
	        std::vector<float> BJetPhi;
	        std::vector<float> BJetE;

	        std::vector<float> NonBJetPt;
	        std::vector<float> NonBJetEta;
	        std::vector<float> NonBJetPhi;
	        std::vector<float> NonBJetE;

	        std::vector<float> LeptonPt;
	        std::vector<float> LeptonEta;
	        std::vector<float> LeptonPhi;
	        std::vector<float> LeptonE;
	        
	        vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	        vector <int> MapIndex_NonBindex;

	        
			    ttree[(dataSetName).c_str()]->GetEntry(j);
          for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
          {
          
                if( (MotherpdgID[i_Jet] == 25) || (fabs(MotherpdgID[i_Jet]) == 6) ) Signal_Matched++;
                if( (MotherpdgID[i_Jet] == 24) || (fabs(MotherpdgID[i_Jet]) == 6) ) SMTT_Matched++;
 
                if(bDiscJet[i_Jet]  > btagWP)
                {
                    BJetPt.push_back(pt_jet[i_Jet]);
                    BJetEta.push_back(eta_jet[i_Jet]);
                    BJetPhi.push_back(phi_jet[i_Jet]);
                    BJetE.push_back(E_jet[i_Jet]);
                    
                    MapIndex_Bindex.push_back(i_Jet);
                }
                else
                {
                    NonBJetPt.push_back(pt_jet[i_Jet]);
                    NonBJetEta.push_back(eta_jet[i_Jet]);
                    NonBJetPhi.push_back(phi_jet[i_Jet]);
                    NonBJetE.push_back(E_jet[i_Jet]);

                    MapIndex_NonBindex.push_back(i_Jet);
                }
          }
          
          LeptonPt.push_back(lepPt);
          LeptonEta.push_back(lepEta);
          LeptonPhi.push_back(lepPhi);
          LeptonE.push_back(lepE);

          kfit->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
          kfit->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
          if(channel == "_El") kfit->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          else if(channel == "_Mu") kfit->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          kfit->SetMet(met_px,met_py);


          if(BJetPt.size()>=nBaselinebjets_ && NonBJetPt.size()>=nBaselinenonbjets_)// BASELINE SELECTION
          {
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              if(KinFitMethod == "SMttHypo")
              {
                  if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
                  {
                  
	                    NSelectionPassedEvents_SMttHypo++;
	                    if(SMTT_Matched == 4) nMCMatchedPassedEvents_SMttHypo_TTbackground++;
	                    if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_SMttHypo_TTsignal++;
	                    if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SMttHypo_STsignal++;

                      kfit->Run();
                      int nPerm_SMttHypo = kfit->GetNPerm();
                                
                      double LowestDisc_SMttHypo = kfit->GetDisc(); //The minimum of likelihood == the best jet-combination
                      double BDTscore = -9999.;
                      int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                      int HighestBDT_IndexAllJetColl_BJETHAD = -1;
                      int HighestBDT_IndexAllJetColl_NONBJET1 = -1;
                      int HighestBDT_IndexAllJetColl_NONBJET2 = -1;
                                
                      for(int ip=0;ip<nPerm_SMttHypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);

		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kfit->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kfit->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kfit->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];
		                      
		                      //Counting whether the TopKinFit method can match to the correct signal
                          if(chi2 == LowestDisc_SMttHypo)
                          {
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                                              {
	                                                NMCIdentifiedEvents_SMttHypo_STsignal_TKF++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                                              {
	                                              NMCIdentifiedEvents_SMttHypo_TTsignal_TKF++;
	                                            }
                                              if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 24) && (MotherpdgID[IndexAllJetColl_NONBJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) == 5) )
                                              {
	                                              NMCIdentifiedEvents_SMttHypo_TTbackground_TKF++;
	                                            }
                                          }                      
                          } 

                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJETHAD_, NONBJET1_, NONBJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          BJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETHAD],eta_jet[IndexAllJetColl_BJETHAD],phi_jet[IndexAllJetColl_BJETHAD],E_jet[IndexAllJetColl_BJETHAD]);
                          NONBJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET1],eta_jet[IndexAllJetColl_NONBJET1],phi_jet[IndexAllJetColl_NONBJET1],E_jet[IndexAllJetColl_NONBJET1]);
                          NONBJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET2],eta_jet[IndexAllJetColl_NONBJET2],phi_jet[IndexAllJetColl_NONBJET2],E_jet[IndexAllJetColl_NONBJET2]);
                          
                          Higgs_ = NONBJET1_+NONBJET2_;
                          HadTop_ = Higgs_+BJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;


                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_FullReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                              HighestBDT_IndexAllJetColl_BJETHAD = IndexAllJetColl_BJETHAD;
                              HighestBDT_IndexAllJetColl_NONBJET1 = IndexAllJetColl_NONBJET1;
                              HighestBDT_IndexAllJetColl_NONBJET2 = IndexAllJetColl_NONBJET2;
                          }
                      }//jet permutation loop

                      //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) )
                                              {
	                                                NMCIdentifiedEvents_SMttHypo_STsignal++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) != 5) )
                                              {
	                                              NMCIdentifiedEvents_SMttHypo_TTsignal++;
	                                            }
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 5) )
                                              {
	                                              NMCIdentifiedEvents_SMttHypo_TTbackground++;
	                                            }
                                          }

                  }//TOPTOPLEPHAD selection
              }
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "TThypo")
              {
                  if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
                  {
	                    NSelectionPassedEvents_TThypo++;
	                    if(SMTT_Matched == 4) nMCMatchedPassedEvents_TThypo_TTbackground++;
	                    if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_TThypo_TTsignal++;
	                    if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_TThypo_STsignal++;
                                
                      kfit->Run();
                      int nPerm_TThypo = kfit->GetNPerm();
                                

                      double LowestDisc_TThypo = kfit->GetDisc(); //The minimum of likelihood == the best jet-combination
                      double BDTscore = -9999.;
                      int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                      int HighestBDT_IndexAllJetColl_NONBJETHAD = -1;
                      int HighestBDT_IndexAllJetColl_BJET1 = -1;
                      int HighestBDT_IndexAllJetColl_BJET2 = -1;

                                
                      for(int ip=0;ip<nPerm_TThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
		                       
		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kfit->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                          if(chi2 == LowestDisc_TThypo)
                          {
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                                              {
	                                                NMCIdentifiedEvents_TThypo_STsignal_TKF++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                                              {
	                                              NMCIdentifiedEvents_TThypo_TTsignal_TKF++;
	                                            }
                                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 24) && (MotherpdgID[IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) == 5) )
                                              {
	                                              NMCIdentifiedEvents_TThypo_TTbackground_TKF++;
	                                            }
                                          }
                          }

                          
                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, NONBJETHAD_, BJET1_, BJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          HadTop_ = Higgs_+NONBJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_FullReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                                  
                              if(BDTscore_tmp > BDTscore)
                              {
                                  BDTscore = BDTscore_tmp;
                                  HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                                  HighestBDT_IndexAllJetColl_NONBJETHAD = IndexAllJetColl_NONBJETHAD;
                                  HighestBDT_IndexAllJetColl_BJET1 = IndexAllJetColl_BJET1;
                                  HighestBDT_IndexAllJetColl_BJET2 = IndexAllJetColl_BJET2;
                              }


                      }

                      //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) )
                                              {
	                                                NMCIdentifiedEvents_TThypo_STsignal++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) != 5) )
                                              {
	                                              NMCIdentifiedEvents_TThypo_TTsignal++;
	                                            }
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 5) )
                                              {
	                                              NMCIdentifiedEvents_TThypo_TTbackground++;
	                                            }
                                          }

                  }//TOPTOPLEPHBB selection
              }
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "SThypo")
              {
                  if(BJetPt.size()>=3)
                  {
	                    NSelectionPassedEvents_SThypo++;
	                    if(SMTT_Matched >= 3) nMCMatchedPassedEvents_SThypo_TTbackground++;
	                    if(Signal_Matched >= 3 && !SingleTop) nMCMatchedPassedEvents_SThypo_TTsignal++;
	                    if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SThypo_STsignal++;
                  
                      kfit->Run();
                      int nPerm_SThypo = kfit->GetNPerm();
                                
                      double LowestDisc_SThypo = kfit->GetDisc(); //The minimum of likelihood == the best jet-combination
                      double BDTscore = -9999.;
                      int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                      int HighestBDT_IndexAllJetColl_BJET1 = -1;
                      int HighestBDT_IndexAllJetColl_BJET2 = -1;
                                
                      for(int ip=0;ip<nPerm_SThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
		                       
		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPHLEPBB,ip)];

                          if(chi2 == LowestDisc_SThypo)
                          {
                                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5))
                                              {
	                                              if(SingleTop) NMCIdentifiedEvents_SThypo_STsignal_TKF++;
	                                              else NMCIdentifiedEvents_SThypo_TTsignal_TKF++;
	                                            }
                                              if( (MotherpdgID[IndexAllJetColl_BJET1] == 24) && (MotherpdgID[IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5))
                                              {
	                                              NMCIdentifiedEvents_SThypo_TTbackground_TKF++;
	                                            }
                          }
                          
                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                          TLorentzVector Higgs_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);

                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);

                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                              if(BDTscore_tmp > BDTscore)
                              {
                                  BDTscore = BDTscore_tmp;
                                  HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                                  HighestBDT_IndexAllJetColl_BJET1 = IndexAllJetColl_BJET1;
                                  HighestBDT_IndexAllJetColl_BJET2 = IndexAllJetColl_BJET2;
                              }
                      }

                      //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5))
                                              {
	                                              if(SingleTop) NMCIdentifiedEvents_SThypo_STsignal++;
	                                              else NMCIdentifiedEvents_SThypo_TTsignal++;
	                                            }
                                              if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5))
                                              {
	                                              NMCIdentifiedEvents_SThypo_TTbackground++;
	                                            }


                  }//TOPHLEPBB selection
              }
        }
		}//for-loop events
		              

    if(KinFitMethod == "SMttHypo")
    {
		    cout << "************ TOPTOPLEPHAD ************" << endl;
        
        //cout efficiencies for TopKinFit method
        cout << " TOPTOPLEPHAD & " << 100*double(NSelectionPassedEvents_SMttHypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_STsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal_TKF)/double(nMCMatchedPassedEvents_SMttHypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_SMttHypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTbackground)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_SMttHypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }    


        //cout efficiencies for MVA method
        cout << " MVA & " << 100*double(NSelectionPassedEvents_SMttHypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_STsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(nMCMatchedPassedEvents_SMttHypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(nMCMatchedPassedEvents_SMttHypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTbackground)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(nMCMatchedPassedEvents_SMttHypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "TThypo")
    {
		    cout << "************ TOPTOPLEPHBB ************" << endl;

        //cout TopKinFit method
        cout << " TOPTOPLEPHBB & " << 100*double(NSelectionPassedEvents_TThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_STsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal_TKF)/double(nMCMatchedPassedEvents_TThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_TThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTbackground)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_TThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }    


        //cout MVA method
        cout << " MVA & " << 100*double(NSelectionPassedEvents_TThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_STsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(nMCMatchedPassedEvents_TThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(nMCMatchedPassedEvents_TThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTbackground)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(nMCMatchedPassedEvents_TThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "SThypo")
    {
		    cout << "************ TOPHLEPBB ************" << endl;

        //cout TopKinFit method efficiencies
        cout << " TOPHLEPBB & " << 100*double(NSelectionPassedEvents_SThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_STsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal_TKF)/double(nMCMatchedPassedEvents_SThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_SThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTbackground)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_SThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }


        //cout MVA method efficiencies
        cout << " MVA & " << 100*double(NSelectionPassedEvents_SThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_STsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(nMCMatchedPassedEvents_SThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(nMCMatchedPassedEvents_SThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTbackground)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(nMCMatchedPassedEvents_SThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }
    }
    cout << " " << endl;		              
  }//for-loop datasets



  delete Eventcomputer_FullReco_;
  delete Eventcomputer_PartReco_;
}







void MVA_Analysis(int skipEvents, std::string SignalName, std::string KinFitMethod, KINFIT::kfit *kfit, vector<std::string> MVAvars, int nBaselinebjets_, int nBaselinenonbjets_, float btagWP, bool doEventTraining)
{

  ///////////////////////////////////////////////////////////////////
  // Initializing MVA
  ///////////////////////////////////////////////////////////////////
  cout << "Initializing MVA for correct Jet combination" << endl;
  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);
  string pathMVA_ = "MVA/TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);

  std::string TrainMethod = KinFitMethod+"_"+intToStr(MVAvars.size())+"Vars"+"_"+intToStr(nBaselinebjets_)+"B"+intToStr(nBaselinenonbjets_)+"J"+channel;

  MVAComputer* Eventcomputer_FullReco_ =0;   
  MVAComputer* Eventcomputer_PartReco_ =0;   

  Eventcomputer_FullReco_ = new MVAComputer("BDT",pathMVA_+"TrainedJetCombMVA_FullReco_"+TrainMethod+"_"+SignalName+".root","TrainedJetCombMVA_PartReco_"+TrainMethod,MVAvars, "test");
  Eventcomputer_PartReco_ = new MVAComputer("BDT",pathMVA_+"TrainedJetCombMVA_PartReco_"+TrainMethod+"_"+SignalName+".root","TrainedJetCombMVA_PartReco_"+TrainMethod,MVAvars, "test");

  MVATrainer*  Eventtrainer_ = 0;
//  MVAComputer* Eventcomputer_ =0;   

  vector<std::string> MVAvars_;

      MVAvars_.push_back("SumCharge_Hjets");
      MVAvars_.push_back("SumCharge_TopJets");
      MVAvars_.push_back("SumCharge_FCNHJetLep");
      MVAvars_.push_back("CvsL_Hjet1");
      MVAvars_.push_back("CvsL_Hjet2");
      MVAvars_.push_back("CvsL_SMb");
      MVAvars_.push_back("CvsL_FCNHjet");
//      MVAvars_.push_back("CvsB_Hjet1");
//      MVAvars_.push_back("CvsB_Hjet2");
//      MVAvars_.push_back("CvsB_SMb");
//      MVAvars_.push_back("CvsB_FCNHjet");
      MVAvars_.push_back("Hmass");
      MVAvars_.push_back("TransvLepTopmass");
//      MVAvars_.push_back("HadTopmass");
      MVAvars_.push_back("DR_H_HadTop");
      MVAvars_.push_back("DPhi_H_LepTop");
      MVAvars_.push_back("LepTopPt");
//      MVAvars_.push_back("HadTopPt");
      MVAvars_.push_back("JetCombBDT");
  

  if(doEventTraining)
  {
      Eventtrainer_ = new MVATrainer("BDT","EventMVA_"+TrainMethod, pathMVA_+"EventMVA_"+TrainMethod+"_"+SignalName+".root");
        for(unsigned int N_var = 0; N_var < MVAvars_.size(); N_var++)
      {
          	Eventtrainer_->bookInputVar(MVAvars_[N_var]);
      }
  }

//  else Eventcomputer_ = new MVAComputer("BDT",pathMVA_+"EventMVA_"+TrainMethod+"_"+SignalName+".root","EventMVA"+TrainMethod,MVAvars_, "test");



  ////////////////////////////////////////////////////////////
  // Load Datasets
  //////////////////////////////////////////////////////////////////////
 	const char *xmlfile = xmlNom.c_str();
 	cout << "Using config file: " << xmlfile << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
		dataSetName = datasets[d]->Name();

    bool SignalSample = false;
    if(dataSetName.find("NP_overlay")!=string::npos) SignalSample = true;

		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "ObjectVarsTree";
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
    /////////////////////////////////////////
    // Get object variables
    ////////////////////////////////////////
    int NumberOfJets;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
    
    double InclJetCharge[20];
    double bDiscJet[20];
    double CvsBJet[20];
    double CvsLJet[20];
    double pdgID[20];
    double MotherpdgID[20];
    double GrandMotherpdgID[20];
    double lepCharge;
    double lepPt;
    double lepEta;
    double lepPhi;
    double lepE;

    double pt_jet[20];
    double phi_jet[20];
    double eta_jet[20];
    double E_jet[20];

    double met_px;
    double met_py;


    ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
    if(channel == "_Mu")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
    }
    if(channel == "_El")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
    }
		  


		//////////////////////////////////////////////////////////
		// Running on events
		//////////////////////////////////////////////////////////
		int eventStart_;
		int eventEnd_;
		
		if(doEventTraining)
		{
		    eventStart_ = 0;
		    eventEnd_ = int(nEntries/skipEvents);
		}
		else
		{
		    eventStart_ = int(nEntries/skipEvents)+1;
		    eventEnd_ = 2*int(nEntries/skipEvents);
		}
		for (int j = 0; j<100000; j++)
		{
		
          int SMTT_Matched = 0;
          int Signal_Matched = 0;
	        
	        std::vector<float> BJetPt;
	        std::vector<float> BJetEta;
	        std::vector<float> BJetPhi;
	        std::vector<float> BJetE;

	        std::vector<float> NonBJetPt;
	        std::vector<float> NonBJetEta;
	        std::vector<float> NonBJetPhi;
	        std::vector<float> NonBJetE;

	        std::vector<float> LeptonPt;
	        std::vector<float> LeptonEta;
	        std::vector<float> LeptonPhi;
	        std::vector<float> LeptonE;
	        
	        vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	        vector <int> MapIndex_NonBindex;

	        
			    ttree[(dataSetName).c_str()]->GetEntry(j);
          for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
          {
                if(bDiscJet[i_Jet]  > btagWP)
                {
                    BJetPt.push_back(pt_jet[i_Jet]);
                    BJetEta.push_back(eta_jet[i_Jet]);
                    BJetPhi.push_back(phi_jet[i_Jet]);
                    BJetE.push_back(E_jet[i_Jet]);
                    
                    MapIndex_Bindex.push_back(i_Jet);
                }
                else
                {
                    NonBJetPt.push_back(pt_jet[i_Jet]);
                    NonBJetEta.push_back(eta_jet[i_Jet]);
                    NonBJetPhi.push_back(phi_jet[i_Jet]);
                    NonBJetE.push_back(E_jet[i_Jet]);

                    MapIndex_NonBindex.push_back(i_Jet);
                }
          }
          
          LeptonPt.push_back(lepPt);
          LeptonEta.push_back(lepEta);
          LeptonPhi.push_back(lepPhi);
          LeptonE.push_back(lepE);

          kfit->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
          kfit->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
          if(channel == "_El") kfit->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          else if(channel == "_Mu") kfit->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
          kfit->SetMet(met_px,met_py);


          double HighestBDT_SumCharge_Hjets = -9999;
          double HighestBDT_SumCharge_TopJets = -9999;
          double HighestBDT_SumCharge_FCNHJetLep = -9999;
          double HighestBDT_CvsL_Hjet1 = -9999;
          double HighestBDT_CvsL_Hjet2 = -9999;
          double HighestBDT_CvsL_SMb = -9999;
          double HighestBDT_CvsL_FCNHjet = -9999;
//          double HighestBDT_CvsB_Hjet1 = -9999;
//          double HighestBDT_CvsB_Hjet2 = -9999;
//          double HighestBDT_CvsB_SMb = -9999;
//          double HighestBDT_CvsB_FCNHjet = -9999;
          double HighestBDT_Hmass = -9999;
          double HighestBDT_TransvLepTopmass = -9999;
//          double HighestBDT_HadTopmass = -9999;
          double HighestBDT_DR_H_HadTop = -9999;
          double HighestBDT_DPhi_H_LepTop = -9999;
          double HighestBDT_LepTopPt = -9999;
//          double HighestBDT_HadTopPt = -9999;


          if(BJetPt.size()>=nBaselinebjets_ && NonBJetPt.size()>=nBaselinenonbjets_)// BASELINE SELECTION
          {
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              if(KinFitMethod == "SMttHypo")
              {
                  if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
                  {
                      kfit->Run();
                      int nPerm_SMttHypo = kfit->GetNPerm();
                                
                      double BDTscore = -9999.;
                      int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                      int HighestBDT_IndexAllJetColl_BJETHAD = -1;
                      int HighestBDT_IndexAllJetColl_NONBJET1 = -1;
                      int HighestBDT_IndexAllJetColl_NONBJET2 = -1;



                                
                      for(int ip=0;ip<nPerm_SMttHypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);

		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kfit->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kfit->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                      int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kfit->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];

                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJETHAD_, NONBJET1_, NONBJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          BJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETHAD],eta_jet[IndexAllJetColl_BJETHAD],phi_jet[IndexAllJetColl_BJETHAD],E_jet[IndexAllJetColl_BJETHAD]);
                          NONBJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET1],eta_jet[IndexAllJetColl_NONBJET1],phi_jet[IndexAllJetColl_NONBJET1],E_jet[IndexAllJetColl_NONBJET1]);
                          NONBJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET2],eta_jet[IndexAllJetColl_NONBJET2],phi_jet[IndexAllJetColl_NONBJET2],E_jet[IndexAllJetColl_NONBJET2]);
                          
                          Higgs_ = NONBJET1_+NONBJET2_;
                          HadTop_ = Higgs_+BJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;

                          HighestBDT_SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
                          HighestBDT_SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
                          HighestBDT_SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
                          HighestBDT_CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
                          HighestBDT_CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
                          HighestBDT_CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                          HighestBDT_CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
//                          HighestBDT_CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
//                          HighestBDT_CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
//                          HighestBDT_CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                          HighestBDT_CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                          HighestBDT_Hmass = Higgs_.M();
                          HighestBDT_TransvLepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
//                          HighestBDT_HadTopmass = HadTop_.M();
                          HighestBDT_DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                          HighestBDT_DPhi_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                          HighestBDT_LepTopPt = LepTop_.Pt();
//                          HighestBDT_HadTopPt = HadTop_.Pt();


                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_FullReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
                              HighestBDT_SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
                              HighestBDT_SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
                              HighestBDT_CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
                              HighestBDT_CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
                              HighestBDT_CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                              HighestBDT_CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
//                              HighestBDT_CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
//                              HighestBDT_CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
//                              HighestBDT_CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                              HighestBDT_CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                              HighestBDT_Hmass = Higgs_.M();
                              HighestBDT_TransvLepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
//                              HighestBDT_HadTopmass = HadTop_.M();
                              HighestBDT_DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              HighestBDT_DPhi_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              HighestBDT_LepTopPt = LepTop_.Pt();
//                              HighestBDT_HadTopPt = HadTop_.Pt();
                          }
                      }//jet permutation loop
                      
                      /////////////////////////////
                      // Do the event training
                      /////////////////////////////
                      if(SignalSample)
                      {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("S","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("S","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("S","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("S","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("S","JetCombBDT", BDTscore);
                      }
                      else
                      {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("B","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("B","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("B","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("B","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("B","JetCombBDT", BDTscore);
                      }
                      
                      
                  }//TOPTOPLEPHAD selection
              }
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "TThypo")
              {
                  if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
                  {
                      kfit->Run();
                      int nPerm_TThypo = kfit->GetNPerm();
                      double BDTscore = -9999.;

                      for(int ip=0;ip<nPerm_TThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
		                       
		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kfit->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                          
                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, NONBJETHAD_, BJET1_, BJET2_;
                          TLorentzVector Higgs_, HadTop_, LepTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          HadTop_ = Higgs_+NONBJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_FullReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double HadTopmass = HadTop_.M();
                              double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();
                              double HadTopPt = HadTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                              Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                                  
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
                              HighestBDT_SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
                              HighestBDT_SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
                              HighestBDT_CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                              HighestBDT_CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                              HighestBDT_CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                              HighestBDT_CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
//                              HighestBDT_CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
//                              HighestBDT_CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
//                              HighestBDT_CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                              HighestBDT_CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                              HighestBDT_Hmass = Higgs_.M();
                              HighestBDT_TransvLepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
//                              HighestBDT_HadTopmass = HadTop_.M();
                              HighestBDT_DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              HighestBDT_DPhi_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              HighestBDT_LepTopPt = LepTop_.Pt();
//                              HighestBDT_HadTopPt = HadTop_.Pt();
                          }
                      }//jet permutation loop
                      
                      /////////////////////////////
                      // Do the event training
                      /////////////////////////////
                      if(SignalSample)
                      {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("S","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("S","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("S","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("S","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("S","JetCombBDT", BDTscore);
                      }
                      else
                      {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("B","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("B","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("B","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("B","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("B","JetCombBDT", BDTscore);
                      }

                  }//TOPTOPLEPHBB selection
              }
              /////////////////////////////////////////////////////////////////////////////
              // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
              /////////////////////////////////////////////////////////////////////////////
              else if(KinFitMethod == "SThypo")
              {
                  if(BJetPt.size()>=3)
                  {
                      kfit->Run();
                      int nPerm_SThypo = kfit->GetNPerm();
                                
                      double LowestDisc_SThypo = kfit->GetDisc(); //The minimum of likelihood == the best jet-combination
                      double BDTscore = -9999.;
                                
                      for(int ip=0;ip<nPerm_SThypo;ip++)
                      {
                          double chi2 = kfit->GetDisc(ip);
		                       
		                      int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kfit->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET1 = MapIndex_Bindex[kfit->GetIndex(BJET1_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_BJET2 = MapIndex_Bindex[kfit->GetIndex(BJET2_TOPHLEPBB,ip)];
		                      int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[0];//Assign the highest pT non-b jet as thee BJETHAD candidate
                          
                          double nuPx = kfit->GetNuPx(ip,0);
                          double nuPy = kfit->GetNuPy(ip,0);
                          double nuPz = kfit->GetNuPz(ip,0);
                          TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_, NONBJETHAD_;
                          TLorentzVector Higgs_, LepTop_, HadTop_;
                          
                          Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                          LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                          BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                          NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                          BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                          BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                          
                          Higgs_ = BJET1_+BJET2_;
                          HadTop_ = Higgs_+NONBJETHAD_;
                          LepTop_ = Nu_ + LEP_ + BJETLEP_;
                          
                          double BDTscore_tmp = -9999.;
                          if(chi2<10E+9)
                          {
                              double Hmass = Higgs_.M();
                              double LepTopmass = LepTop_.M();
                              double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                              double LepTopPt = LepTop_.Pt();

                              Eventcomputer_FullReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_FullReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_FullReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_FullReco_->FillVar("LepTopPt", LepTopPt);

                              std::map<std::string,Float_t> MVAVals = Eventcomputer_FullReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                   BDTscore_tmp = it->second;
                              }
                          }
                          else//MVA_PartReco
                          {

                              double Hmass = Higgs_.M();
                              double LepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
                              double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              double LepTopPt = LepTop_.Pt();

                              Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                              Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                              Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                              Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);

                              std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  BDTscore_tmp = it->second;
                              }
                          }                              
                                  
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
                              HighestBDT_SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
                              HighestBDT_SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
                              HighestBDT_CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                              HighestBDT_CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                              HighestBDT_CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                              HighestBDT_CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
//                              HighestBDT_CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
//                              HighestBDT_CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
//                              HighestBDT_CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                              HighestBDT_CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                              HighestBDT_Hmass = Higgs_.M();
                              HighestBDT_TransvLepTopmass = sqrt(2*(LEP_+BJETLEP_).Pt() * Nu_.Pt() * (1-cos( (LEP_+BJETLEP_).DeltaPhi( Nu_ )) ) );
//                              HighestBDT_HadTopmass = HadTop_.M();
                              HighestBDT_DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                              HighestBDT_DPhi_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                              HighestBDT_LepTopPt = LepTop_.Pt();
//                              HighestBDT_HadTopPt = HadTop_.Pt();
                          }
                      }//jet permutation loop
                      
                      /////////////////////////////
                      // Do the event training
                      /////////////////////////////
                      if(SignalSample)
                      {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("S","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("S","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("S","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("S","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("S","JetCombBDT", BDTscore);
                      }
                      else
                      {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", HighestBDT_SumCharge_Hjets);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", HighestBDT_SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", HighestBDT_SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", HighestBDT_CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", HighestBDT_CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", HighestBDT_CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", HighestBDT_CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", HighestBDT_CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", HighestBDT_CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", HighestBDT_CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", HighestBDT_CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", HighestBDT_Hmass);
                              Eventtrainer_->Fill("B","TransvLepTopmass", HighestBDT_TransvLepTopmass);
//                              Eventtrainer_->Fill("B","HadTopmass", HighestBDT_HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", HighestBDT_DR_H_HadTop);
                              Eventtrainer_->Fill("B","DPhi_H_LepTop", HighestBDT_DPhi_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", HighestBDT_LepTopPt);
//                              Eventtrainer_->Fill("B","HadTopPt", HighestBDT_HadTopPt);
                              Eventtrainer_->Fill("B","JetCombBDT", BDTscore);
                      }

                  }//TOPHLEPBB selection
              }
        }
		}//for-loop events
  }//for-loop datasets


  Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test",false);


  delete Eventtrainer_;
  delete Eventcomputer_FullReco_;
  delete Eventcomputer_PartReco_;
}





void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string NTupleName)
{
  	cout<<""<<endl;
  	cout<<"RUNNING NOMINAL DATASETS"<<endl;
  	cout<<""<<endl;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
  	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  
    //***************************************************CREATING PLOT****************************************************
  	string plotname = s_varofInterest;
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, s_varofInterest.c_str()); 

  
  	//***********************************************OPEN FILES & GET NTUPLES**********************************************
  	string dataSetName, filepath;
  	int nEntries;
  	float ScaleFactor, NormFactor;
  	int varofInterest;
  	Double_t d_varofInterest;
 	  int n_object = 0;
  	double v_d_varofInterest_double [20];
 
  	vector<string> v;
  	// to avoid modifying original string
  	// first duplicate the original string and return a char pointer then free the memory
  
  	char delim[] = "[]";
  	char * dup = strdup(s_varofInterest.c_str());
  	char * token = strtok(dup, delim);//split string of variable according to the delim
  	while(token != NULL){
    	v.push_back(string(token));
    	// the call is treated as a subsequent calls to strtok:
    	// the function continues from where it left in previous invocation
    	token = strtok(NULL, delim);
  	}
  	free(dup);


     if (v.size() == 2)//Meaning we have a variable of the form "var[n_obj]", which is an array of values for the variable 'var'
     {

                //If plotting a variable which consists of several values (e.g.: jet_pt contains the pt of all jets), make also plots for the individual values (e.g.: plot the pt for every jet separately). For now, only done for 5 first objects
                for(int iToPlot = 1; iToPlot <= maxNumbObjToPlot; iToPlot++)
                {
                          string conv_str;
                          ostringstream conv;   // stream used for the conversion
                          conv << (iToPlot);      // insert the textual representation of 'Number' in the characters in the stream
                          conv_str = "_"+conv.str(); // set 'Result' to the contents of the stream

                          MSPlot[(v[0]+conv_str).c_str()] = new MultiSamplePlot(datasets, (v[0]+conv_str).c_str(), nBins, plotLow, plotHigh, (v[0]+conv_str).c_str()); 
                }                

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = NTupleName;	
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;
		                
		                
                   // bo logic to set the right branch address depending on the string given as argument of the datasetplotter
	                 ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),v_d_varofInterest_double); //v[0] is the string of the variable you want to plot. This variable should be an array of values, according to the number of objects
	                 ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&n_object); // v[1] is the string of the variable between [] in the string. This should correspond to the number of objects


		              bool isData= false;
		              bool isAMC = true;
		              if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		              {
		                if(debug) cout << "Data found" << endl;
		                isData =true;
	                }
                  if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos || dataSetName.find("amc") !=std::string::npos) isAMC = true;


                  ////////////////////////////////////////////////////////////
                  // Tree for reweighting
                  ////////////////////////////////////////////////////////////		  
		              string TTreename_Weights = "Weights";	
		              string TTreename_NtupleInfo = "NtupleInfoTree";	
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttree for each dataset
		              ttree[(dataSetName + "NtupleInfoTree").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_NtupleInfo.c_str()); //get ttree for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->GetEntry(0);
                  Luminosity = luminosity_;
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);

                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(isAMC && !isData)
                    {
                       
                        for (int k = 0; k<nEntries; k++)
                        {
                           ttree[(dataSetName + "weights").c_str()]->GetEntry(k);
                           if( nloweight > 0) nPos++;
                           else if( nloweight < 0) nNeg ++;
                           Ev ++; 
                         }
                         
                         nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
                     }		
		
		
		              Int_t ThreeJets = 0;
		              Int_t MoreThanThreeJets = 0;
		              //////////////////////////////////////////////////////////
		              // Making MS plots
		              //////////////////////////////////////////////////////////
		              for (int j = 0; j<nEntries; j++)
		              {
                  		ScaleFactor = 1.; // event scale factor
			                ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                if(ScaleFactor < 0) ScaleFactor = 0;
			                ttree[dataSetName.c_str()]->GetEntry(j);
			                
			                if(n_object == 3) ThreeJets++;
			                else MoreThanThreeJets++; 

                      for(int i_obj = 0; i_obj < n_object;  i_obj++)
                      {
                          string conversion_str;
                          ostringstream convert;   // stream used for the conversion
                          convert << (1+i_obj);      // insert the textual representation of 'Number' in the characters in the stream
                          conversion_str = "_"+convert.str(); // set 'Result' to the contents of the stream

			                    if(debug) cout << "varofInterest is " << v_d_varofInterest_double[i_obj] << endl;
			                    if(isData)
			                    {// for data, fill once per event, weighted with the event scale factor
				                    MSPlot[plotname.c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, 1.);
				                    if(i_obj< maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, 1.);//Fill MSPlot for first 5 variables
			                    }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, ScaleFactor*Luminosity);
				                    if(i_obj<maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true,  ScaleFactor*Luminosity);//Fill MSPlot for first 5 variables
			                    }
			                }
			                
			                
		              }
		              
		              
		              cout << dataSetName << ": Number of Events with exactly 3 jets: " << ThreeJets << endl;
		              cout << dataSetName << ": Number of Events with more than 3 jets: " << MoreThanThreeJets << endl;
		              
		          }//for-loop datasets
               

      }//end statement on variable-plotting consisting of array		   (v.size()==2)
     else if (v.size() == 1)
     {
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = NTupleName;	
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;



                  bool isInteger = false;
	                if (v[0].compare(0,2,"I_") == 0)//these are the variables that are an integer
	                {
	                  ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);
	                  isInteger = true;
	                }
	                else //The others are doubles
	                {
	                  ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&d_varofInterest);
	                }

          		    bool isData= false;
		              bool isAMC = false;
		              if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		              {
		                if(debug) cout << "Data found" << endl;
		                isData =true;
	                }
                  if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos || dataSetName.find("amc") !=std::string::npos) isAMC = true;


                  ////////////////////////////////////////////////////////////
                  // Tree for reweighting
                  ////////////////////////////////////////////////////////////		  
		              string TTreename_Weights = "Weights";	
		              string TTreename_NtupleInfo = "NtupleInfoTree";	
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttree for each dataset
		              ttree[(dataSetName + "NtupleInfoTree").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_NtupleInfo.c_str()); //get ttree for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->GetEntry(0);
                  Luminosity = luminosity_;
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);

                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(isAMC && !isData)
                    {
                       
                        for (int k = 0; k<nEntries; k++)
                        {
                           ttree[(dataSetName + "weights").c_str()]->GetEntry(k);
                           if( nloweight > 0) nPos++;
                           else if( nloweight < 0) nNeg ++;
                           Ev ++; 
                         }
                         
                         nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
                     }		
		
		              //////////////////////////////////////////////////////////
		              // Making MS plots
		              //////////////////////////////////////////////////////////
		              for (int j = 0; j<nEntries; j++)
		              {
                  		ScaleFactor = 1.; // event scale factor
			                ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                ttree[dataSetName.c_str()]->GetEntry(j);

			                if(isInteger)
			                {
			                    if(debug) cout << "varofInterest is " << varofInterest << endl;
			                    if(isData)
			                    {// for data, fill once per event, weighted with the event scale factor
				                    MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, 1.);
				                  }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
			                    }
			                }
			                else
			                {
			                		if(debug) cout << "varofInterest is " << d_varofInterest << endl;
			                    if(isData)
			                    {// for data, fill once per event, weighted with the event scale factor
				                    MSPlot[plotname.c_str()]->Fill(d_varofInterest, datasets[d], true, 1.);
				                  }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(d_varofInterest, datasets[d], true, ScaleFactor*Luminosity);
			                    }
			                }
			                		if(debug) cout << "Event " << j << endl;

			            }
			                
			                
		           }//for-loop datasets


     }//end of statement if variable is not an array of values
     else 
     {
	      cout << "Vector of string does not have the good size!!!" << endl;
      }
      
      




//	treeLoader.UnLoadDataset();
  	// clearing vector
  	v.clear();
  	if (debug){
    	cout << "after cleaning" << endl ;
    	cout << "v.size() is " << v.size() << endl;
  	}
  
cout << "MSPlot size: " << MSPlot.size() << endl;      


};






// function that writes all the MSPlots created in a root file
void MSPCreator ()
{

  	string pathPNG = "MSPlots/";
  	mkdir(pathPNG.c_str(),0777);
  	pathPNG += "MSPlots";
  	pathPNG += channel;
  	mkdir(pathPNG.c_str(),0777);
  	pathPNG += "/";
  	pathPNG += date;
  	mkdir(pathPNG.c_str(),0777);
  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  	TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  	outfile->cd();


  	// Loop over all the MSPlots
  	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      	string name = it->first;
      	MultiSamplePlot *temp = it->second;
      	//if (debug){
			cout << "Saving the MSP" << endl;
			cout << " and it->first is " << it->first << endl;
			cout << " Luminosity is " << Luminosity << endl;
      	//}
      	temp->setDataLumi(Luminosity);
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
        cout << "Drawing MSP: " << it->second << endl;
		    temp->Draw("MyMSP_"+it->first, 1, false, false, false, 1);
      	temp->Write(outfile, it->first, true,pathPNG, "png");
	}

  	outfile->Write("kOverwrite");
}

// function that converts an int into a string
std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}

