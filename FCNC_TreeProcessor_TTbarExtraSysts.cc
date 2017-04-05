#include "TStyle.h"
#include "TPaveText.h"
#include <TLatex.h>

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
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/LHAPDF.h"
#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/Reweighting.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"


using namespace std;
using namespace TopTree;
//using namespace LHAPDF;
//using namespace KINFIT;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlot_nPV;

void GetEnvelopeScale();
void GetEnvelopePDF();
double maximumValue(vector<double> array);
double minimumValue(vector<double> array);
double OptimalCut_CombTraining(string category, string coupling);

bool applySumWeights_scaleEnvelope = true;

// functions prototype
string intToStr (int number);
void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug);
bool applyPDFs = true;

inline bool FileExists (const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{

    if(argc < 7)
    {
        cerr << "INVALID number of arguments. The necessary arguments are: " << endl;
        cout << "    int baseline_bjets             = strtol(argv[1], NULL,10);" << endl;
        cout << "    int baseline_jets                 = strtol(argv[2], NULL,10);" << endl;
        cout << "    string SignalSample            = argv[3];" << endl;
        cout << "    string channel            = argv[4];" << endl;
        cout << "    string date            = argv[5];" << endl;
        cout << "    bool PVreweighing = strtol(argv[6], NULL,10);" << endl;
        cout << "    bool doControlPlots  = strtol(argv[7], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[8], NULL,10);" << endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    string SignalSample  = argv[3];//Valid arguments are: hut,  hct
    string channel            = argv[4];
    string date            = argv[5];
    bool PVreweighing = strtol(argv[6], NULL,10);
    bool doControlPlots = strtol(argv[7], NULL,10);
    bool debug         =strtol(argv[8], NULL,10);

    bool split_ttbar = true;   
    
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
    string TrainingName = "Training_" + SignalSample + channel + "_" +  category;//Example: Training_SThut_El_b3j3

    vector<string> WhatSysts;
    
    WhatSysts.push_back("weight1");   
    WhatSysts.push_back("weight2");   
    WhatSysts.push_back("weight3");   
    WhatSysts.push_back("weight4");   
    WhatSysts.push_back("weight6");   
    WhatSysts.push_back("weight8");   
    WhatSysts.push_back("_hdampup");   
    WhatSysts.push_back("_hdampdown");   
    WhatSysts.push_back("_isrup");
    WhatSysts.push_back("_isrdown");   
    WhatSysts.push_back("_fsrup");   
    WhatSysts.push_back("_fsrdown");   
    WhatSysts.push_back("_UEup");
    WhatSysts.push_back("_UEdown");   
    WhatSysts.push_back("");   


    map<string,string> namingConventionFit;

    namingConventionFit["_isrup"] = "_isrup";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["_isrdown"] = "_isrdown";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["_fsrup"] = "_fsrup";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["_fsrdown"] = "_fsrdown";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight1"] = "_weight1";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight2"] = "_weight2";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight3"] = "_weight3";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight4"] = "_weight4";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight6"] = "_weight6";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["weight8"] = "_weight8";   //Will contain later on the envelope from weight*, isrup, isrdown, fsrup and fsrdown
    namingConventionFit["_UEup"] = "_UEUp";   
    namingConventionFit["_UEdown"] = "_UEDown";   
    namingConventionFit["_hdampup"] = "_hdampUp";   
    namingConventionFit["_hdampdown"] = "_hdampDown";   
    namingConventionFit[""] = "";


    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << " - Category: " << category << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();


    string xmlNom = "config/FullMcBkgdSamples_ExtraSystematicSamples_TreeProcessor.xml";

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
  	vector < Dataset* > datasets_splittedTTbar; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  	string dataSetName, filepath;
    float Luminosity = 0;

    ///////////////////////////////////////////////////////////////////
    //// S p l i t t i n g   T T b a r ////////////////////////////////
    Dataset* ttbar_ll = new Dataset();
    Dataset* ttbar_cc = new Dataset();
    Dataset* ttbar_bb = new Dataset();

    
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
	  {
		    dataSetName = datasets[d]->Name();
		    datasets_splittedTTbar.push_back(datasets[d]);
		    
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        Luminosity = datasets[d]->EquivalentLumi();
        }
        else if (dataSetName.find("TTJets")!=string::npos && split_ttbar)
        {
            
            ttbar_ll->SetName("TTJets_ll");
            ttbar_cc->SetName("TTJets_cc");
            ttbar_bb->SetName("TTJets_bb");

            ttbar_ll->SetTitle("tt+lf");
            ttbar_cc->SetTitle("tt+cc");
            ttbar_bb->SetTitle("tt+bb");


            ttbar_ll->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            ttbar_cc->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            ttbar_bb->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            
            ttbar_ll->SetColor(kAzure+7);
            ttbar_cc->SetColor(kAzure+5);
            ttbar_bb->SetColor(kAzure+3);
            
            datasets_splittedTTbar.pop_back();
            datasets_splittedTTbar.push_back(ttbar_bb);
            datasets_splittedTTbar.push_back(ttbar_cc);
            datasets_splittedTTbar.push_back(ttbar_ll);

            cout << " - split TTBar dataset into ..."  << ttbar_ll->Name() << ", " << ttbar_cc->Name() << " and " << ttbar_ll->Name()  << endl;

        }     
    }

    if(Luminosity == 0)
    {
            cout << "Luminosity is 0. Please check the data-luminosity in your xml file. Exiting program..." << endl;
            return 1;
    }

    //Storing the datasetNames in a vector for which the variables are plotted
    //This will be used later on in the tool to plot the error bands, so do not store the Data name and NewPhysics names
    vector <string> datasetnames_backgrounds;
	  for (int d = 0; d < datasets_splittedTTbar.size(); d++)   //Loop through datasets
	  {
          string n = datasets_splittedTTbar[d]->Name();
          if(n.find("Data")!=string::npos || n.find("NP_")!=string::npos) continue;
          datasetnames_backgrounds.push_back(n);
    }
    //***************************************************CREATING PLOT****************************************************
    //Format of MSPlots: MultiSamplePlot(vector<Dataset*> datasets, string PlotName, int Nbins, float Min, float Max, string XaxisLabel, string YaxisLabel, string Text, string Units)

    string xaxislabelcoupling = " #kappa_{";
    if(SignalSample == "hut") xaxislabelcoupling +=  "Hut}";
    else if(SignalSample == "hct") xaxislabelcoupling +=  "Hct}";
    
    for(int iSyst = 0; iSyst<WhatSysts.size();iSyst++)
    {
        histo1D[("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);

    }
    for(int iCount = 0; iCount < 101; iCount++)
    {
        histo1D[("maxSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("maxSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str(),("maxSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("maxSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str(),("maxSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("maxSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str(),("maxSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);

        histo1D[("ST_ttbb_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("ST_ttbb_pdfweight"+intToStr(iCount)).c_str(),("ST_ttbb_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("ST_ttcc_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("ST_ttcc_pdfweight"+intToStr(iCount)).c_str(),("ST_ttcc_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("ST_ttlf_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("ST_ttlf_pdfweight"+intToStr(iCount)).c_str(),("ST_ttlf_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);

        histo1D[("TT_ttbb_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("TT_ttbb_pdfweight"+intToStr(iCount)).c_str(),("TT_ttbb_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("TT_ttcc_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("TT_ttcc_pdfweight"+intToStr(iCount)).c_str(),("TT_ttcc_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("TT_ttlf_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("TT_ttlf_pdfweight"+intToStr(iCount)).c_str(),("TT_ttlf_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);

        histo1D[("combSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_ttbb_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_ttcc_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_ttlf_pdfweight"+intToStr(iCount)).c_str(),20,-1,1);

        histo1D[("combSTandTT_cutCount_ttbb_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_cutCount_ttbb_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_cutCount_ttbb_pdfweight"+intToStr(iCount)).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttcc_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_cutCount_ttcc_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_cutCount_ttcc_pdfweight"+intToStr(iCount)).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttlf_pdfweight"+intToStr(iCount)).c_str()] = new TH1F(("combSTandTT_cutCount_ttlf_pdfweight"+intToStr(iCount)).c_str(),("combSTandTT_cutCount_ttlf_pdfweight"+intToStr(iCount)).c_str(),1,-1,1);
        
    }  
 


  	//***********************************************RUNNING OVER DATASETS**********************************************
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	  {

        dataSetName = datasets[d]->Name();
        bool isData= false;
		    bool isAMC = false;
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        if(debug) cout << "Data found" << endl;
		        isData =true;
		        continue;
	      }
	      
        for(int Counter = 0; Counter < WhatSysts.size(); Counter++)
        {
            string postfix = WhatSysts[Counter];
            if(isData && WhatSysts[Counter] != "") continue;
            if(!isData) postfix = WhatSysts[Counter];
            if(WhatSysts[Counter].find("weight")!=string::npos) postfix = "";

		        cout<<"Dataset:  :"<<(dataSetName+WhatSysts[Counter]).c_str()<<endl;
		        filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + postfix + ".root";
		        if (debug)
		        {
		            cout<<"filepath: "<<filepath<<endl;
                cout <<"Equivalent luminosity of the dataset is: " << datasets[d]->EquivalentLumi() << endl;
		        }

            //*************Variables to be used for Reading the training********************
            Float_t LepCharge_;
            Float_t MVA_TOPTOPLEPHAD_;
            Float_t MVA_TOPTOPLEPHBB_;
            Float_t MVA_TOPHLEPBB_hut_;
            Float_t MVA_TOPHLEPBB_hct_;
	          Float_t HiggsMass_TOPHLEPBB_hut_;
	          Float_t HiggsMass_TOPHLEPBB_hct_;
	          Float_t HiggsEta_TOPHLEPBB_hut_;
	          Float_t HiggsEta_TOPHLEPBB_hct_;
	          Float_t TopLepMass_TOPHLEPBB_hut_;
	          Float_t TopLepMass_TOPHLEPBB_hct_;
            Float_t TopLepPt_TOPHLEPBB_hut_;
            Float_t TopLepPt_TOPHLEPBB_hct_;
            Float_t TopLepEta_TOPHLEPBB_hut_;
            Float_t TopLepEta_TOPHLEPBB_hct_;
            Float_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_;
            Float_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_;
            Float_t TopLepHiggsDr_TOPHLEPBB_hut_;
            Float_t TopLepHiggsDr_TOPHLEPBB_hct_;
            Float_t HiggsBJet1CSVv2_TOPHLEPBB_hut_;
            Float_t HiggsBJet1CSVv2_TOPHLEPBB_hct_;
            Float_t HiggsBJet2CSVv2_TOPHLEPBB_hut_;
            Float_t HiggsBJet2CSVv2_TOPHLEPBB_hct_;
            Float_t TopLepBJetCSVv2_TOPHLEPBB_hut_;
            Float_t TopLepBJetCSVv2_TOPHLEPBB_hct_;
            Float_t TopHadMass_TOPTOPLEPHAD_;
            Float_t TopLepMass_TOPTOPLEPHAD_;
            Float_t TopLepTopHadDr_TOPTOPLEPHAD_;
            Float_t TopLepBJetCSVv2_TOPTOPLEPHAD_;
            Float_t TopHadBJetCSVv2_TOPTOPLEPHAD_;
            Float_t TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_;
            Float_t TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_;
            Float_t HiggsMass_TOPTOPLEPHBB_;
            Float_t TopLepMass_TOPTOPLEPHBB_;
            Float_t HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_;
            Float_t TopLepHiggsDr_TOPTOPLEPHBB_;
            Float_t HiggsBJet1CSVv2_TOPTOPLEPHBB_;
            Float_t HiggsBJet2CSVv2_TOPTOPLEPHBB_;
            Float_t TopLepBJetCSVv2_TOPTOPLEPHBB_;
            Float_t TopHadNonBJetCSVv2_TOPTOPLEPHBB_;
            

            //********************INITIALIZING MVA READER********************************
	          TMVA::Reader* reader_ST = new TMVA::Reader("!Color:!Silent");
	          TMVA::Reader* reader_TT = new TMVA::Reader("!Color:!Silent");
	          TMVA::Reader* reader_combSTandTT = new TMVA::Reader("!Color:!Silent");

            if(TrainingName.find("hut")!=string::npos && TrainingName.find("j3")!=string::npos)
            {
	              reader_ST->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("LepCharge",&LepCharge_);
	              reader_ST->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hut_);

	              reader_TT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hut_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hut_);

	              reader_combSTandTT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("LepCharge",&LepCharge_);
	              reader_combSTandTT->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hut_);
            }
            else if(TrainingName.find("hct")!=string::npos && TrainingName.find("j3")!=string::npos)
            {
	              reader_ST->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hct_);


	              reader_TT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hct_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hct_);

	              reader_combSTandTT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hct_);
            }  
            else if(TrainingName.find("hut")!=string::npos && TrainingName.find("j4")!=string::npos)
            {
	              reader_ST->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("LepCharge",&LepCharge_);
	              reader_ST->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hut_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);

	              reader_TT->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);


	              reader_combSTandTT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("LepCharge",&LepCharge_);
	              reader_combSTandTT->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hut_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB_);
            }
            else if(TrainingName.find("hct")!=string::npos && TrainingName.find("j4")!=string::npos)
            {
	              reader_ST->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hct_);
	              reader_ST->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_ST->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);

	              reader_TT->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB_);
	              reader_TT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_TT->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);

	              reader_combSTandTT->AddVariable("HiggsMass_TOPHLEPBB",&HiggsMass_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("MVA_TOPHLEPBB",&MVA_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("HiggsEta_TOPHLEPBB",&HiggsEta_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPHLEPBB",&TopLepMass_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopLepPt_TOPHLEPBB",&TopLepPt_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepEta_TOPHLEPBB",&TopLepEta_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPHLEPBB",&TopLepHiggsDr_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",&HiggsBJet1CSVv2_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",&HiggsBJet2CSVv2_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",&TopLepBJetCSVv2_TOPHLEPBB_hct_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_);
	              reader_combSTandTT->AddVariable("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB_);
	              reader_combSTandTT->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB_);
            }
            else
            {
                cerr << "No correct signal selected" << endl;
                return 0;
            }
            
	          std::string weightsFile_ST= "weights/Training_ST" + SignalSample + channel + "_" +  category+"_BDT.weights.xml";
	          std::string weightsFile_TT= "weights/Training_TT" + SignalSample + channel + "_" +  category+"_BDT.weights.xml";
	          std::string weightsFile_combSTandTT= "weights/CombTraining_" + SignalSample + channel + "_" +  category+"_BDT.weights.xml";
	          reader_ST->BookMVA("BDTG method",weightsFile_ST.c_str());
	          reader_TT->BookMVA("BDTG method",weightsFile_TT.c_str());
	          reader_combSTandTT->BookMVA("BDTG method",weightsFile_combSTandTT.c_str());






	
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
		                    
		                    


      	    //***********************************************IMPORTING VARIABLES**********************************************
		        string TTreename = "ObjectVarsTree";	
		        string TTreename_info = "NtupleInfoTree";	
		        ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		        ttree[(dataSetName+TTreename_info).c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_info.c_str()); //get ntuple creation information

            int nEntries;

		        nEntries = ttree[dataSetName.c_str()]->GetEntries();
		        cout<<"                 nEntries: "<<nEntries<<endl;

            //----------------------------------------------//
            //Get The object variables + weights
            //----------------------------------------------//
            //Weights
            Double_t W_puSF;
            Double_t W_fleptonSF;
            Double_t W_btagWeight_shape;
//            Double_t W_weight0; //Nominal weight, not to be applied
            Double_t W_weight1;
            Double_t W_weight2;
            Double_t W_weight3;
            Double_t W_weight4;
//            Double_t W_weight5; //Unphysical weight
            Double_t W_weight6;
//            Double_t W_weight7; //Unphysical weight
            Double_t W_weight8; 
            Double_t W_hdamp_up; 
            Double_t W_hdamp_dw; 
          
            Int_t run_num;
            Int_t evt_num;
            Int_t lumi_num;
            Int_t nvtx;
            Int_t npu;
            Int_t genTTX;


            //variable for  leptons
            Int_t LepCharge;
      
            //variable for jets 
            Int_t nJets;
	          Int_t nJets_CSVM; 

            //Info for PDF calculation
            Int_t id1;
            Int_t id2;
            Float_t x1;
            Float_t x2;
            Float_t q;
	          
	          //JetIndices_correctJetComb
            Double_t MVA_TOPTOPLEPHAD = -999.;
            Double_t MVA_TOPTOPLEPHBB = -999.;
            Double_t MVA_TOPHLEPBB_hut = -999.;
            Double_t MVA_TOPHLEPBB_hct = -999.;
            //Variables for signal/background training
	          Double_t HiggsMass_TOPHLEPBB_hut;
	          Double_t HiggsMass_TOPHLEPBB_hct;
	          Double_t HiggsEta_TOPHLEPBB_hut;
	          Double_t HiggsEta_TOPHLEPBB_hct;
	          Double_t TopLepMass_TOPHLEPBB_hut;
	          Double_t TopLepMass_TOPHLEPBB_hct;
            Double_t TopLepPt_TOPHLEPBB_hut;
            Double_t TopLepPt_TOPHLEPBB_hct;
            Double_t TopLepEta_TOPHLEPBB_hut;
            Double_t TopLepEta_TOPHLEPBB_hct;
            Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
            Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
            Double_t TopLepHiggsDr_TOPHLEPBB_hut;
            Double_t TopLepHiggsDr_TOPHLEPBB_hct;
            Double_t HiggsBJet1CSVv2_TOPHLEPBB_hut;
            Double_t HiggsBJet1CSVv2_TOPHLEPBB_hct;
            Double_t HiggsBJet2CSVv2_TOPHLEPBB_hut;
            Double_t HiggsBJet2CSVv2_TOPHLEPBB_hct;
            Double_t TopLepBJetCSVv2_TOPHLEPBB_hut;
            Double_t TopLepBJetCSVv2_TOPHLEPBB_hct;
            Double_t TopHadMass_TOPTOPLEPHAD;
            Double_t TopLepMass_TOPTOPLEPHAD;
            Double_t TopLepTopHadDr_TOPTOPLEPHAD;
            Double_t TopLepBJetCSVv2_TOPTOPLEPHAD;
            Double_t TopHadBJetCSVv2_TOPTOPLEPHAD;
            Double_t TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
            Double_t TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
            Double_t HiggsMass_TOPTOPLEPHBB;
            Double_t TopLepMass_TOPTOPLEPHBB;
            Double_t HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
            Double_t TopLepHiggsDr_TOPTOPLEPHBB;
            Double_t HiggsBJet1CSVv2_TOPTOPLEPHBB;
            Double_t HiggsBJet2CSVv2_TOPTOPLEPHBB;
            Double_t TopLepBJetCSVv2_TOPTOPLEPHBB;
            Double_t TopHadNonBJetCSVv2_TOPTOPLEPHBB;

            
            // Weights
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape",&W_btagWeight_shape); 
//            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight0",&W_weight0);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight1",&W_weight1);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight2",&W_weight2);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight3",&W_weight3); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight4",&W_weight4);  
//            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight5",&W_weight5); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight6",&W_weight6);  
//            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight7",&W_weight7); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight8",&W_weight8);  
//            ttree[(dataSetName).c_str()]->SetBranchAddress("W_hdamp_up",&W_hdamp_up);  
//            ttree[(dataSetName).c_str()]->SetBranchAddress("W_hdamp_dw",&W_hdamp_dw);  

            ttree[(dataSetName).c_str()]->SetBranchAddress("I_run_num",&run_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_evt_num",&evt_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_lumi_num",&lumi_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_npu",&npu);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_genTTX",&genTTX);



            //SelectedLepton
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_LepCharge",&LepCharge);
            
            // jets
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&nJets);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&nJets_CSVM);

            ttree[(dataSetName).c_str()]->SetBranchAddress("I_id1",&id1);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_id2",&id2);
            ttree[(dataSetName).c_str()]->SetBranchAddress("x1",&x1);
            ttree[(dataSetName).c_str()]->SetBranchAddress("x2",&x2);
            ttree[(dataSetName).c_str()]->SetBranchAddress("q",&q);
           
            // Jet-indices associated to the jet-assignment in the bMVA method
            ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hut",&MVA_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hct",&MVA_TOPHLEPBB_hct);
            //Variables for signal/background training
	          ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPHLEPBB_hut",&HiggsMass_TOPHLEPBB_hut);
	          ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPHLEPBB_hct",&HiggsMass_TOPHLEPBB_hct);
	          ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsEta_TOPHLEPBB_hut",&HiggsEta_TOPHLEPBB_hut);
	          ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsEta_TOPHLEPBB_hct",&HiggsEta_TOPHLEPBB_hct);
	          ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPHLEPBB_hut",&TopLepMass_TOPHLEPBB_hut);
	          ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPHLEPBB_hct",&TopLepMass_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepPt_TOPHLEPBB_hut",&TopLepPt_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepPt_TOPHLEPBB_hct",&TopLepPt_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepEta_TOPHLEPBB_hut",&TopLepEta_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepEta_TOPHLEPBB_hct",&TopLepEta_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPHLEPBB_hut",&TopLepHiggsDr_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPHLEPBB_hct",&TopLepHiggsDr_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB_hut",&HiggsBJet1CSVv2_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB_hct",&HiggsBJet1CSVv2_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB_hut",&HiggsBJet2CSVv2_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB_hct",&HiggsBJet2CSVv2_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB_hut",&TopLepBJetCSVv2_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB_hct",&TopLepBJetCSVv2_TOPHLEPBB_hct);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
            ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);
                      


            
            double EqLumi = datasets[d]->EquivalentLumi();
            int EntryStart = 0;

            if(WhatSysts[Counter] == "_isrup") EqLumi = 70077.8481774;
            else if(WhatSysts[Counter] == "_isrdown") EqLumi = 35755.2406944;
            else if(WhatSysts[Counter] == "_fsrdown") EqLumi = 35577.7940752;
            else if(WhatSysts[Counter] == "_fsrup") EqLumi = 35799.6357122;
            else if(WhatSysts[Counter] == "_UEup")
            {
                EqLumi = 35238.2898913;
                Luminosity = Luminosity * 1.003487;
            }
            else if(WhatSysts[Counter] == "_UEdown")
            {
                EqLumi = 34088.4822545;
                Luminosity = Luminosity * 1.003487;
            }
            else if(WhatSysts[Counter] == "_hdampup")
            {
                EqLumi = 35042.3078773;
                Luminosity = Luminosity * 1.003487;
            }
            else if(WhatSysts[Counter] == "_hdampdown")
            {
                EqLumi =  34857.4276234;
                Luminosity = Luminosity * 1.003487;
            }
            else 
            {
                if(category=="b2j4")
                {
                    EntryStart = (int) 119*nEntries/120;
                    EqLumi = EqLumi/120;
                }
                else if(category=="b2j3")
                {
                    EntryStart = (int) 59*nEntries/60;
                    EqLumi = EqLumi/60;
                }
                else if(category=="b3j4")
                {
                    EntryStart = (int) 9*nEntries/10;
                    EqLumi = EqLumi/10;
                }
                else if(category=="b3j3" || category=="b4j4")
                {
                    EntryStart = (int) nEntries/2;
                    EqLumi = EqLumi/2;
                }
            }
//            EntryStart = nEntries/2+nEntries/3+1;//Manually overwriting the number of events to run over.

            Double_t average_weight1 = 0.;
            Double_t average_weight2 = 0.;
            Double_t average_weight3 = 0.;
            Double_t average_weight4 = 0.;
            Double_t average_weight6 = 0.;
            Double_t average_weight8 = 0.;

            if(dataSetName.find("TTJets") != string::npos)
            {
                int nEventsPassed = 0;
                for (int k = 0; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);

                    average_weight1 = average_weight1 + W_weight1;
                    average_weight2 = average_weight2 + W_weight2;
                    average_weight3 = average_weight3 + W_weight3;
                    average_weight4 = average_weight4 + W_weight4;
                    average_weight6 = average_weight6 + W_weight6;
                    average_weight8 = average_weight8 + W_weight8;
                    nEventsPassed++;
                }
                if(nEventsPassed != 0)
                {
                    average_weight1 = average_weight1/nEventsPassed;
                    average_weight2 = average_weight2/nEventsPassed;
                    average_weight3 = average_weight3/nEventsPassed;
                    average_weight4 = average_weight4/nEventsPassed;
                    average_weight6 = average_weight6/nEventsPassed;
                    average_weight8 = average_weight8/nEventsPassed;
                }
            }


      	    //***********************************************RUNNING OVER EVENTS**********************************************
		        for (int j = EntryStart; j<nEntries; j++)
		        {
		                  
                if(debug)
                {
                    if(!isData) cin.get();
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
                Dataset * Sample = 0;
                if (dataSetName.find("TTJets")!=string::npos && split_ttbar)
                {
                    bool isttbb = (genTTX == 051 || genTTX == 151 || genTTX == 251 ||
		                  genTTX == 052 || genTTX == 152 || genTTX == 252 ||
		                  genTTX == 053 || genTTX == 153 || genTTX == 253 ||
		                  genTTX == 054 || genTTX == 154 || genTTX == 254 ||
		                  genTTX == 055 || genTTX == 155 || genTTX == 255);
                   
                    bool isttcc = (genTTX == 041 || genTTX == 141 || genTTX == 241 ||
		                  genTTX == 042 || genTTX == 142 || genTTX == 242 ||
		                  genTTX == 043 || genTTX == 143 || genTTX == 243 ||
		                  genTTX == 044 || genTTX == 144 || genTTX == 244 ||
		                  genTTX == 045 || genTTX == 145 || genTTX == 245);
                   
                    bool isttlf = (!isttbb && !isttcc);

                    if(isttlf) Sample = ttbar_ll;
                    else if(isttcc) Sample = ttbar_cc;
                    else if(isttbb) Sample = ttbar_bb;
                    
                    if(debug) cout << "   Sample split into " << Sample->Name() << endl;
                }
                else Sample = datasets[d];

                //////////////////////////////////////
                //Applying the scale factors
                ///////////////////////////////////////
                double ScaleFactor = 1.;
                double W_puSF_applied = 1.;

			              if(!PVreweighing) W_puSF_applied = W_puSF;
			              else
			              {
			                  W_puSF_applied = W_nPV.ITweight( (int)nvtx );
			              }

                    if(debug)
                    {
                        //Safety triggers in case there are strange things happening in the event weights
                        if(W_fleptonSF < 0 || W_btagWeight_shape < 0 || Luminosity < 0 || W_puSF_applied < 0)
                        {
                              cout << "----- Event " << j << " has a negative weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; Luminosity=" << Luminosity << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              cout << "----- The event will be skipped....." << endl;
                              continue;
                        }
                        else if(W_fleptonSF != W_fleptonSF || W_btagWeight_shape != W_btagWeight_shape || W_puSF_applied != W_puSF_applied)
                        {
                              cout << "----- Event " << j << " has a Nan weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              cout << "----- The event will be skipped....." << endl;
                              continue;
                        }
                        else if(W_fleptonSF >= 40 || W_btagWeight_shape >= 40 || W_puSF_applied >= 40)
                        {
                              cout << "----- Event " << j << " has a weight larger than 40. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              //cout << "----- The event will be skipped....." << endl;
                              //continue;
                        }
                    }                

	          if( HiggsMass_TOPHLEPBB_hut > 500. ) HiggsMass_TOPHLEPBB_hut = 500.;
	          if( TopLepMass_TOPHLEPBB_hut > 500. ) TopLepMass_TOPHLEPBB_hut = 500.;
	          if( TopLepPt_TOPHLEPBB_hut > 1000. ) TopLepPt_TOPHLEPBB_hut = 1000.;
	          if( HiggsMass_TOPHLEPBB_hct > 500. ) HiggsMass_TOPHLEPBB_hct = 500.;
	          if( TopLepMass_TOPHLEPBB_hct > 500. ) TopLepMass_TOPHLEPBB_hct = 500.;
	          if( TopLepPt_TOPHLEPBB_hct > 1000. ) TopLepPt_TOPHLEPBB_hct = 1000.;
	          if( TopLepMass_TOPTOPLEPHAD > 500.) TopLepMass_TOPTOPLEPHAD = 500.;
	          if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	          if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	          if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
            if( HiggsBJet1CSVv2_TOPHLEPBB_hut < 0.) HiggsBJet1CSVv2_TOPHLEPBB_hut= 0.;
            if( HiggsBJet1CSVv2_TOPHLEPBB_hct < 0.) HiggsBJet1CSVv2_TOPHLEPBB_hct= 0.;
            if( HiggsBJet2CSVv2_TOPHLEPBB_hut < 0.) HiggsBJet2CSVv2_TOPHLEPBB_hut= 0.;
            if( HiggsBJet2CSVv2_TOPHLEPBB_hct < 0.) HiggsBJet2CSVv2_TOPHLEPBB_hct= 0.;
            if( TopLepBJetCSVv2_TOPHLEPBB_hut < 0.) TopLepBJetCSVv2_TOPHLEPBB_hut= 0.;
            if( TopLepBJetCSVv2_TOPHLEPBB_hct < 0.) TopLepBJetCSVv2_TOPHLEPBB_hct= 0.;
            if( TopLepBJetCSVv2_TOPTOPLEPHAD < 0.) TopLepBJetCSVv2_TOPTOPLEPHAD= 0.;
            if( TopHadBJetCSVv2_TOPTOPLEPHAD < 0.) TopHadBJetCSVv2_TOPTOPLEPHAD= 0.;
            if( TopHadWNonBJet1CSVv2_TOPTOPLEPHAD < 0.) TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = 0.;
            if( TopHadWNonBJet2CSVv2_TOPTOPLEPHAD < 0.) TopHadWNonBJet2CSVv2_TOPTOPLEPHAD= 0.;
            if( HiggsBJet1CSVv2_TOPTOPLEPHBB < 0.) HiggsBJet1CSVv2_TOPTOPLEPHBB= 0.;
            if( HiggsBJet2CSVv2_TOPTOPLEPHBB < 0.) HiggsBJet2CSVv2_TOPTOPLEPHBB= 0.;
            if( TopLepBJetCSVv2_TOPTOPLEPHBB < 0.) TopLepBJetCSVv2_TOPTOPLEPHBB= 0.;
            if( TopHadNonBJetCSVv2_TOPTOPLEPHBB < 0.) TopHadNonBJetCSVv2_TOPTOPLEPHBB = 0.;

                

                LepCharge_ = (float) LepCharge;
                MVA_TOPTOPLEPHAD_ = (float) MVA_TOPTOPLEPHAD;
                MVA_TOPTOPLEPHBB_ = (float) MVA_TOPTOPLEPHBB;
                MVA_TOPHLEPBB_hut_ = (float) MVA_TOPHLEPBB_hut;
                MVA_TOPHLEPBB_hct_ = (float) MVA_TOPHLEPBB_hct;
	              HiggsMass_TOPHLEPBB_hut_ = (float) HiggsMass_TOPHLEPBB_hut;
	              HiggsMass_TOPHLEPBB_hct_ = (float) HiggsMass_TOPHLEPBB_hct;
	              HiggsEta_TOPHLEPBB_hut_ = (float) HiggsEta_TOPHLEPBB_hut;
	              HiggsEta_TOPHLEPBB_hct_ = (float) HiggsEta_TOPHLEPBB_hct;
	              TopLepMass_TOPHLEPBB_hut_ = (float) TopLepMass_TOPHLEPBB_hut;
	              TopLepMass_TOPHLEPBB_hct_ = (float) TopLepMass_TOPHLEPBB_hct;
                TopLepPt_TOPHLEPBB_hut_ = (float) TopLepPt_TOPHLEPBB_hut;
                TopLepPt_TOPHLEPBB_hct_ = (float) TopLepPt_TOPHLEPBB_hct;
                TopLepEta_TOPHLEPBB_hut_ = (float) TopLepEta_TOPHLEPBB_hut;
                TopLepEta_TOPHLEPBB_hct_ = (float) TopLepEta_TOPHLEPBB_hct;
                HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut_ = (float) HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
                HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct_ = (float) HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
                TopLepHiggsDr_TOPHLEPBB_hut_ = (float) TopLepHiggsDr_TOPHLEPBB_hut;
                TopLepHiggsDr_TOPHLEPBB_hct_ = (float) TopLepHiggsDr_TOPHLEPBB_hct;
                HiggsBJet1CSVv2_TOPHLEPBB_hut_ = (float) HiggsBJet1CSVv2_TOPHLEPBB_hut;
                HiggsBJet1CSVv2_TOPHLEPBB_hct_ = (float) HiggsBJet1CSVv2_TOPHLEPBB_hct;
                HiggsBJet2CSVv2_TOPHLEPBB_hut_ = (float) HiggsBJet2CSVv2_TOPHLEPBB_hut;
                HiggsBJet2CSVv2_TOPHLEPBB_hct_ = (float) HiggsBJet2CSVv2_TOPHLEPBB_hct;
                TopLepBJetCSVv2_TOPHLEPBB_hut_ = (float) TopLepBJetCSVv2_TOPHLEPBB_hut;
                TopLepBJetCSVv2_TOPHLEPBB_hct_ = (float) TopLepBJetCSVv2_TOPHLEPBB_hct;
                TopHadMass_TOPTOPLEPHAD_ = (float) TopHadMass_TOPTOPLEPHAD;
                TopLepMass_TOPTOPLEPHAD_ = (float) TopLepMass_TOPTOPLEPHAD;
                TopLepTopHadDr_TOPTOPLEPHAD_ = (float) TopLepTopHadDr_TOPTOPLEPHAD;
                TopLepBJetCSVv2_TOPTOPLEPHAD_ = (float) TopLepBJetCSVv2_TOPTOPLEPHAD;
                TopHadBJetCSVv2_TOPTOPLEPHAD_ = (float) TopHadBJetCSVv2_TOPTOPLEPHAD;
                TopHadWNonBJet1CSVv2_TOPTOPLEPHAD_ = (float) TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
                TopHadWNonBJet2CSVv2_TOPTOPLEPHAD_ = (float) TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
                HiggsMass_TOPTOPLEPHBB_ = (float) HiggsMass_TOPTOPLEPHBB;
                TopLepMass_TOPTOPLEPHBB_ = (float) TopLepMass_TOPTOPLEPHBB;
                HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB_ = (float) HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
                TopLepHiggsDr_TOPTOPLEPHBB_ = (float) TopLepHiggsDr_TOPTOPLEPHBB;
                HiggsBJet1CSVv2_TOPTOPLEPHBB_ = (float) HiggsBJet1CSVv2_TOPTOPLEPHBB;
                HiggsBJet2CSVv2_TOPTOPLEPHBB_ = (float) HiggsBJet2CSVv2_TOPTOPLEPHBB;
                TopLepBJetCSVv2_TOPTOPLEPHBB_ = (float) TopLepBJetCSVv2_TOPTOPLEPHBB;
	              TopHadNonBJetCSVv2_TOPTOPLEPHBB_ = (float) TopHadNonBJetCSVv2_TOPTOPLEPHBB;

                double MVAvalue_ST = reader_ST->EvaluateMVA("BDTG method");
                double MVAvalue_TT = reader_TT->EvaluateMVA("BDTG method");
                double MVAvalue_maxSTandTT = max(MVAvalue_ST,MVAvalue_TT);
                double MVAvalue_combSTandTT = reader_combSTandTT->EvaluateMVA("BDTG method");

                    //Nominal scale factor -- scale factors for systematic shifts are calculated below
                    ScaleFactor *= W_puSF_applied;
                    ScaleFactor *= W_fleptonSF;
                    ScaleFactor *= W_btagWeight_shape;
                    std::vector<double> pdfweights;

               if(applySumWeights_scaleEnvelope)
               {
                    if(WhatSysts[Counter]=="weight1") ScaleFactor *= W_weight1/average_weight1;   
                    else if(WhatSysts[Counter]=="weight2") ScaleFactor *= W_weight2/average_weight2;   
                    else if(WhatSysts[Counter]=="weight3") ScaleFactor *= W_weight3/average_weight3;    
                    else if(WhatSysts[Counter]=="weight4") ScaleFactor *= W_weight4/average_weight4;    
                    else if(WhatSysts[Counter]=="weight6") ScaleFactor *= W_weight6/average_weight6;   
                    else if(WhatSysts[Counter]=="weight8") ScaleFactor *= W_weight8/average_weight8;  
               }
               else
               {
                    if(WhatSysts[Counter]=="weight1") ScaleFactor *= W_weight1;   
                    else if(WhatSysts[Counter]=="weight2") ScaleFactor *= W_weight2;   
                    else if(WhatSysts[Counter]=="weight3") ScaleFactor *= W_weight3;    
                    else if(WhatSysts[Counter]=="weight4") ScaleFactor *= W_weight4;    
                    else if(WhatSysts[Counter]=="weight6") ScaleFactor *= W_weight6;   
                    else if(WhatSysts[Counter]=="weight8") ScaleFactor *= W_weight8;  
               }
               if(WhatSysts[Counter]=="" && applyPDFs)
                    {
                        //PDF weights calculation
                        LHAPDF::setVerbosity(0);
                        LHAPDF::PDFSet basepdfSet("NNPDF30_nlo_as_0118");
                        LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_100"); // give the correct name

                        const LHAPDF::PDF* basepdf = basepdfSet.mkPDF(0);
                        for ( size_t i=0; i<newpdfSet.size(); ++i )
                        {
                          const LHAPDF::PDF* newpdf = newpdfSet.mkPDF(i);
                          const double weight = LHAPDF::weightxxQ(id1, id2, x1, x2, q, *basepdf, *newpdf);
                          pdfweights.push_back(weight);
                          delete newpdf;

                                
                                if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("maxSTandTT_ttbb_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("maxSTandTT_ttcc_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("maxSTandTT_ttlf_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);
                                
                                if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("ST_ttbb_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("ST_ttcc_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("ST_ttlf_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi * weight);
                                
                                if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("TT_ttbb_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("TT_ttcc_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("TT_ttlf_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi * weight);
                                
                                if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_ttbb_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_ttcc_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);
                                else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_ttlf_pdfweight"+intToStr(i)).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi * weight);			                

                                //Cut-and-Count
                                if(MVAvalue_combSTandTT > OptimalCut_CombTraining(category, SignalSample))
                                {
                                    if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_cutCount_ttbb_pdfweight"+intToStr(i)).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi * weight);
                                    else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_cutCount_ttcc_pdfweight"+intToStr(i)).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi * weight);
                                    else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_cutCount_ttlf_pdfweight"+intToStr(i)).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi * weight);			                
                                }
                        }
                        delete basepdf;
                    }


            if(isData) ScaleFactor =1.;
		          
      	        //***********************************************FILLING PLOTS**********************************************
                
                bool ScalePlots = (!isData);

                        //-----------------------------------------------------------------------------------------------------------
                        // Fill Plots
                        //-----------------------------------------------------------------------------------------------------------
                        if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor  / EqLumi );
                        
                        if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("ST_ttbb"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("ST_ttcc"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("ST_ttlf"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor  / EqLumi );
                        
                        if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("TT_ttbb"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("TT_ttcc"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("TT_ttlf"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor  / EqLumi );
                        
                        if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_ttbb"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_ttcc"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi );
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_ttlf"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor  / EqLumi );			                
                                
                        //Cut-and-Count
                        if(MVAvalue_combSTandTT > OptimalCut_CombTraining(category, SignalSample))
                        {
                            if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi );
                            else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi );
                            else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[Counter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor  / EqLumi );			                
                        }

		        }//for-loop events
		    }//for-loop systematics            
    }//for-loop datasets
               

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


    GetEnvelopeScale();
    if(applyPDFs) GetEnvelopePDF();

    //Now store histo's to be used for limit setting
    string outname_limitsetting_maxSTandTT = pathPNG+"/inputExtra_MVA";
    if(SignalSample == "hct") outname_limitsetting_maxSTandTT += "HctMAX_"+category+"_"+SignalSample+".root";
    else if(SignalSample == "hut") outname_limitsetting_maxSTandTT += "HutMAX_"+category+"_"+SignalSample+".root";

    TFile *outfile_limitsetting_maxSTandTT = new TFile(outname_limitsetting_maxSTandTT.c_str(),"recreate");
    outfile_limitsetting_maxSTandTT->cd();

    for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
       	string name = it->first;
/*
        if(name.find("weight") != string::npos || name.find("isrup")!= string::npos 
          || name.find("isrdown") != string::npos || name.find("fsrup") != string::npos || name.find("fsrdown") != string::npos
          )
        {
            continue;
        }
        if(name.find("maxSTandTT") == string::npos) continue;
        if(name.find("Up") == string::npos && name.find("Down") == string::npos) continue;
*/

       	TH1F *temp = it->second;

        temp->Write();
	  }
	  outfile_limitsetting_maxSTandTT->Write("kOverwrite");

    string outname_limitsetting_ST = pathPNG+"/inputExtra_MVA";
    if(SignalSample == "hct") outname_limitsetting_ST += "HctST_"+category+"_"+SignalSample+".root";
    else if(SignalSample == "hut") outname_limitsetting_ST += "HutST_"+category+"_"+SignalSample+".root";

    TFile *outfile_limitsetting_ST = new TFile(outname_limitsetting_ST.c_str(),"recreate");
    outfile_limitsetting_ST->cd();

    for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
       	string name = it->first;

        if(name.find("TT") != string::npos || name.find("combSTandTT") != string::npos || name.find("maxSTandTT") != string::npos || name.find("weight") != string::npos || name.find("isrup")!= string::npos 
          || name.find("isrdown") != string::npos || name.find("fsrup") != string::npos || name.find("fsrdown") != string::npos)//Do not save the pictures of the systematics
        {
            continue;
        }
        if(name.find("Up") == string::npos && name.find("Down") == string::npos) continue;


       	TH1F *temp = it->second;

        temp->Write();
	  }
	  outfile_limitsetting_ST->Write("kOverwrite");

    //Now store histo's to be used for limit setting
    string outname_limitsetting_TT = pathPNG+"/inputExtra_MVA";
    if(SignalSample == "hct") outname_limitsetting_TT += "HctTT_"+category+"_"+SignalSample+".root";
    else if(SignalSample == "hut") outname_limitsetting_TT += "HutTT_"+category+"_"+SignalSample+".root";

    TFile *outfile_limitsetting_TT = new TFile(outname_limitsetting_TT.c_str(),"recreate");
    outfile_limitsetting_TT->cd();

    for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
       	string name = it->first;

        if(name.find("weight") != string::npos || name.find("isrup")!= string::npos 
          || name.find("isrdown") != string::npos || name.find("fsrup") != string::npos || name.find("fsrdown") != string::npos)//Do not save the pictures of the systematics
        {
            continue;
        }
        if(name.find("TT") == string::npos) continue;
        if(name.find("Up") == string::npos && name.find("Down") == string::npos) continue;


       	TH1F *temp = it->second;

        temp->Write();
	  }
	  outfile_limitsetting_TT->Write("kOverwrite");

    //Now store histo's to be used for limit setting
    string outname_limitsetting_combSTandTT = pathPNG+"/inputExtra_MVA";
    if(SignalSample == "hct") outname_limitsetting_combSTandTT += "HctComb_"+category+"_"+SignalSample+".root";
    else if(SignalSample == "hut") outname_limitsetting_combSTandTT += "HutComb_"+category+"_"+SignalSample+".root";

    TFile *outfile_limitsetting_combSTandTT = new TFile(outname_limitsetting_combSTandTT.c_str(),"recreate");
    outfile_limitsetting_combSTandTT->cd();

    for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
       	string name = it->first;

        if(name.find("weight") != string::npos || name.find("isrup")!= string::npos 
          || name.find("isrdown") != string::npos || name.find("fsrup") != string::npos || name.find("fsrdown") != string::npos)//Do not save the pictures of the systematics
        {
            continue;
        }
        if(name.find("combSTandTT") == string::npos) continue;
        if(name.find("Up") == string::npos && name.find("Down") == string::npos) continue;


       	TH1F *temp = it->second;

        temp->Write();
	  }
	  outfile_limitsetting_combSTandTT->Write("kOverwrite");


    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;

}

// function that converts an int into a string
string intToStr (int number)
{
  	ostringstream buff;
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
    else if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    else if(channel == "_All") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
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

    MSPlot_nPV["NPV_unw"] = new MultiSamplePlot(datasets, "NPV_unw", 51, -0.5, 50.5, "Number of PV","Events", ""); 

  
 

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


void GetEnvelopeScale()
{


  vector <string> VariablesToEnvelope;
  VariablesToEnvelope.push_back("maxSTandTT");
  VariablesToEnvelope.push_back("ST");
  VariablesToEnvelope.push_back("TT");
  VariablesToEnvelope.push_back("combSTandTT");

  vector <string> SamplesToEnvelope;
  SamplesToEnvelope.push_back("ttbb");
  SamplesToEnvelope.push_back("ttcc");
  SamplesToEnvelope.push_back("ttlf");
  
  for(int iSample = 0; iSample < SamplesToEnvelope.size(); iSample++)
  {  
      for(int i = 0; i < VariablesToEnvelope.size(); i++)
      {
          string histname_ = VariablesToEnvelope[i]+"_"+SamplesToEnvelope[iSample];

          TH1F * nominal_ = (TH1F*) histo1D[(histname_).c_str()]->Clone();
          TH1F * weight1_ = (TH1F*) histo1D[(histname_+"_weight1").c_str()]->Clone();
          TH1F * weight2_ = (TH1F*) histo1D[(histname_+"_weight2").c_str()]->Clone();
          TH1F * weight3_ = (TH1F*) histo1D[(histname_+"_weight3").c_str()]->Clone();
          TH1F * weight4_ = (TH1F*) histo1D[(histname_+"_weight4").c_str()]->Clone();
          TH1F * weight6_ = (TH1F*) histo1D[(histname_+"_weight6").c_str()]->Clone();
          TH1F * weight8_ = (TH1F*) histo1D[(histname_+"_weight8").c_str()]->Clone();
          TH1F * isrup_ = (TH1F*) histo1D[(histname_+"_isrup").c_str()]->Clone();
          TH1F * isrdown_ = (TH1F*) histo1D[(histname_+"_isrdown").c_str()]->Clone();
          TH1F * fsrup_ = (TH1F*) histo1D[(histname_+"_fsrup").c_str()]->Clone();
          TH1F * fsrdown_ = (TH1F*) histo1D[(histname_+"_fsrdown").c_str()]->Clone();

          histo1D[(histname_+"_scaleEnvelopeUp").c_str()] = new TH1F((histname_+"_scaleEnvelopeUp").c_str(),(histname_+"_scaleEnvelopeUp").c_str(),nominal_->GetNbinsX(),-1,1);
          histo1D[(histname_+"_scaleEnvelopeDown").c_str()] = new TH1F((histname_+"_scaleEnvelopeDown").c_str(),(histname_+"_scaleEnvelopeDown").c_str(),nominal_->GetNbinsX(),-1,1);


          for(int binN = 1; binN < nominal_->GetNbinsX(); ++binN)
          {
              float binContentMax = -1, binContentMin = 1000000000000000;

              vector<double> bincontents;
              bincontents.push_back(nominal_->GetBinContent(binN));
              bincontents.push_back(weight1_->GetBinContent(binN));
              bincontents.push_back(weight2_->GetBinContent(binN));
              bincontents.push_back(weight3_->GetBinContent(binN));
              bincontents.push_back(weight4_->GetBinContent(binN));
              bincontents.push_back(weight6_->GetBinContent(binN));
              bincontents.push_back(weight8_->GetBinContent(binN));
              bincontents.push_back(isrup_->GetBinContent(binN));
              bincontents.push_back(isrdown_->GetBinContent(binN));
              double fsrupDiffNom = sqrt(2)/2*(nominal_->GetBinContent(binN) - fsrup_->GetBinContent(binN));//The recommendation is to downscale fsr variation by sqrt(2)/2
              double fsrdwDiffNom = sqrt(2)/2*(nominal_->GetBinContent(binN) - fsrdown_->GetBinContent(binN));
              bincontents.push_back(fsrup_->GetBinContent(binN)+fsrupDiffNom);
              bincontents.push_back(fsrdown_->GetBinContent(binN)+fsrdwDiffNom);

                  if(binContentMin > minimumValue(bincontents)) binContentMin = minimumValue(bincontents);
                  else cout << "ERROR: no minimal bincontent found" << endl;
                  if(binContentMax < maximumValue(bincontents)) binContentMax = maximumValue(bincontents);
                  else cout << "ERROR: no maximal bincontent found" << endl;
              
              histo1D[(histname_+"_scaleEnvelopeUp").c_str()]->SetBinContent(binN, binContentMax);
              histo1D[(histname_+"_scaleEnvelopeDown").c_str()]->SetBinContent(binN, binContentMin);
          }
      }
      
  }

}


void GetEnvelopePDF()
{

  vector <string> VariablesToEnvelope;
  VariablesToEnvelope.push_back("maxSTandTT");
  VariablesToEnvelope.push_back("ST");
  VariablesToEnvelope.push_back("TT");
  VariablesToEnvelope.push_back("combSTandTT");

  vector <string> SamplesToEnvelope;
  SamplesToEnvelope.push_back("ttbb");
  SamplesToEnvelope.push_back("ttcc");
  SamplesToEnvelope.push_back("ttlf");
  
  for(int iSample = 0; iSample < SamplesToEnvelope.size(); iSample++)
  {  
      for(int i = 0; i < VariablesToEnvelope.size(); i++)
      {
          string histname_ = VariablesToEnvelope[i]+"_"+SamplesToEnvelope[iSample];

          TH1F * nominal_ = (TH1F*) histo1D[(histname_).c_str()]->Clone();

          histo1D[(histname_+"_PDFEnvelopeUp").c_str()] = new TH1F((histname_+"_PDFEnvelopeUp").c_str(),(histname_+"_PDFEnvelopeUp").c_str(),nominal_->GetNbinsX(),-1,1);
          histo1D[(histname_+"_PDFEnvelopeDown").c_str()] = new TH1F((histname_+"_PDFEnvelopeDown").c_str(),(histname_+"_PDFEnvelopeDown").c_str(),nominal_->GetNbinsX(),-1,1);

          for(int binN = 1; binN < nominal_->GetNbinsX(); ++binN)
          {
              float binContentMax = -1, binContentMin = 1000000000000000;
              vector<double> bincontents;


              for(int iCount = 0; iCount < 101; iCount++)
              {
                  bincontents.push_back(histo1D[(histname_+"_pdfweight"+intToStr(iCount)).c_str()]->GetBinContent(binN));             
              }

              if(binContentMin > minimumValue(bincontents)) binContentMin = minimumValue(bincontents);
              else binContentMin = nominal_->GetBinContent(binN);
              if(binContentMax < maximumValue(bincontents)) binContentMax = maximumValue(bincontents);
              else binContentMax = nominal_->GetBinContent(binN);

              histo1D[(histname_+"_PDFEnvelopeUp").c_str()]->SetBinContent(binN, binContentMax);
              histo1D[(histname_+"_PDFEnvelopeDown").c_str()]->SetBinContent(binN, binContentMin);
          }//Bins
      }//Variables
  }//Samples

}



double maximumValue(vector<double> array)
{
     int length = array.size();  // establish size of array
     double max = array[0];       // start with max = first element

     for(int i = 1; i<length; i++)
     {
          if(array[i] > max)
                max = array[i];
     }
     return max;                // return highest value in array
}

double minimumValue(vector<double> array)
{
     int length = array.size();  // establish size of array
     double max = array[0];       // start with max = first element

     for(int i = 1; i<length; i++)
     {
          if(array[i] < max)
                max = array[i];
     }
     return max;                // return highest value in array
}

double OptimalCut_CombTraining(string category, string coupling)
{
    double MVA_cutvalue = -1.;
    if(coupling == "hut")
    {
        if(category == "b2j3") MVA_cutvalue = -0.334589;
        else if(category == "b2j4") MVA_cutvalue = -0.445078;
        else if(category == "b3j3") MVA_cutvalue = -0.258389;
        else if(category == "b3j4") MVA_cutvalue = -0.270175;
        else if(category == "b4j4") MVA_cutvalue = -0.46595;
    }
    else if(coupling == "hct")
    {
        if(category == "b2j3") MVA_cutvalue = -0.380768;
        else if(category == "b2j4") MVA_cutvalue = -0.322955;
        else if(category == "b3j3") MVA_cutvalue = -0.333675;
        else if(category == "b3j4") MVA_cutvalue = -0.277198;
        else if(category == "b4j4") MVA_cutvalue = -0.296879;
    }
    
    return MVA_cutvalue;
}
