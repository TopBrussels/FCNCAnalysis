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

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"


using namespace std;
using namespace TopTree;
//using namespace KINFIT;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,double> SystScaleFactor;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlot_nPV;


bool Manual_XML = false; //This boolean controls the luminosity by hand, the nPV reweighing and which xml file to be taken as input
string manualxml = "config/FullMcBkgdSamples_Manual.xml";

// functions prototype
string intToStr (int number);
void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug);
void MakeTotalSystErrorBand_Distributions(string  outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile);
double WeightPrivateSignalSample(Int_t n_jets, string samplename);

inline bool FileExists (const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{

    if(argc < 8)
    {
        cerr << "INVALID number of arguments. The necessary arguments are: " << endl;
        cout << "    int baseline_bjets             = strtol(argv[1], NULL,10);" << endl;
        cout << "    int baseline_jets                 = strtol(argv[2], NULL,10);" << endl;
        cout << "    string SignalSample            = argv[3];" << endl;
        cout << "    string channel            = argv[4];" << endl;
        cout << "    string date            = argv[5];" << endl;
        cout << "    bool PVreweighing = strtol(argv[6], NULL,10);" << endl;
        cout << "    bool doJESSys  = strtol(argv[7], NULL,10);" << endl;
        cout << "    bool doJERSys  = strtol(argv[8], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[9], NULL,10);" << endl;
        cout << "    bool ApplyPostFit         =strtol(argv[10], NULL,10);" << endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    string SignalSample            = argv[3];
    string channel            = argv[4];
    string date            = argv[5];
    bool PVreweighing = strtol(argv[6], NULL,10);
    bool doJESSys  = strtol(argv[7], NULL,10);
    bool doJERSys  = strtol(argv[8], NULL,10);
    bool debug         =strtol(argv[9], NULL,10);
    bool ApplyPostFit         =strtol(argv[10], NULL,10);

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


    vector<string> WhatSysts;
    
    WhatSysts.push_back("iterativefit_lfPlus");   
    WhatSysts.push_back("iterativefit_lfMinus");   
    WhatSysts.push_back("iterativefit_hfPlus");   
    WhatSysts.push_back("iterativefit_hfMinus");   
    WhatSysts.push_back("iterativefit_lfstats1Plus");   
    WhatSysts.push_back("iterativefit_lfstats1Minus");   
    WhatSysts.push_back("iterativefit_lfstats2Plus");   
    WhatSysts.push_back("iterativefit_lfstats2Minus");   
    WhatSysts.push_back("iterativefit_hfstats1Plus");   
    WhatSysts.push_back("iterativefit_hfstats1Minus");   
    WhatSysts.push_back("iterativefit_hfstats2Plus");   
    WhatSysts.push_back("iterativefit_hfstats2Minus");   
    WhatSysts.push_back("iterativefit_cferr1Plus");   
    WhatSysts.push_back("iterativefit_cferr1Minus");   
    WhatSysts.push_back("iterativefit_cferr2Plus");   
    WhatSysts.push_back("iterativefit_cferr2Minus");   
    WhatSysts.push_back("pileupPlus");   
    WhatSysts.push_back("pileupMinus");   
    WhatSysts.push_back("leptonPlus");
    WhatSysts.push_back("leptonMinus");
    WhatSysts.push_back("TopPtPlus");
    WhatSysts.push_back("TopPtMinus");
    if(doInclusive)
    {
        WhatSysts.push_back("noSF");
        WhatSysts.push_back("OnlyTopPtSF");
        WhatSysts.push_back("OnlyBTagSF");
        WhatSysts.push_back("OnlyPUSF");
        WhatSysts.push_back("OnlyLepSF");
        WhatSysts.push_back("OnlyNLOSF");
        WhatSysts.push_back("NoBTagSF");
        WhatSysts.push_back("NoPUSF");
        WhatSysts.push_back("NoLepSF");
        WhatSysts.push_back("NoNLOSF");
    }
    if(doJESSys) WhatSysts.push_back("JESPlus");
    if(doJESSys) WhatSysts.push_back("JESMinus");
    if(doJERSys) WhatSysts.push_back("JERPlus");
    if(doJERSys) WhatSysts.push_back("JERMinus");
    WhatSysts.push_back("");   


    map<string,string> namingConventionFit;

    namingConventionFit["iterativefit_lfPlus"] = "_SfIteraviveFitLfUp";   
    namingConventionFit["iterativefit_lfMinus"] = "_SfIteraviveFitLfDown";   
    namingConventionFit["iterativefit_hfPlus"] = "_SfIteraviveFitHfUp";   
    namingConventionFit["iterativefit_hfMinus"] = "_SfIteraviveFitHfDown";   
    namingConventionFit["iterativefit_lfstats1Plus"] = "_SfIteraviveFitLfstats1Up";   
    namingConventionFit["iterativefit_lfstats1Minus"] = "_SfIteraviveFitLfstats1Down";   
    namingConventionFit["iterativefit_lfstats2Plus"] = "_SfIteraviveFitLfstats2Up";   
    namingConventionFit["iterativefit_lfstats2Minus"] = "_SfIteraviveFitLfstats2Down";   
    namingConventionFit["iterativefit_hfstats1Plus"] = "_SfIteraviveFitHfstats1Up";   
    namingConventionFit["iterativefit_hfstats1Minus"] = "_SfIteraviveFitHfstats1Down";   
    namingConventionFit["iterativefit_hfstats2Plus"] = "_SfIteraviveFitHfstats2Up";   
    namingConventionFit["iterativefit_hfstats2Minus"] = "_SfIteraviveFitHfstats2Down";   
    namingConventionFit["iterativefit_cferr1Plus"] = "_SfIteraviveFitCferr1Up";   
    namingConventionFit["iterativefit_cferr1Minus"] = "_SfIteraviveFitCferr1Down";   
    namingConventionFit["iterativefit_cferr2Plus"] = "_SfIteraviveFitCferr2Up";   
    namingConventionFit["iterativefit_cferr2Minus"] = "_SfIteraviveFitCferr2Down";   
    namingConventionFit["pileupPlus"] = "_SfPileupUp";   
    namingConventionFit["pileupMinus"] = "_SfPileupDown";   
    namingConventionFit["leptonPlus"] = "_SfLeptonUp";
    namingConventionFit["leptonMinus"] = "_SfLeptonDown";
    namingConventionFit["TopPtPlus"] = "_SfTopPtUp";
    namingConventionFit["TopPtMinus"] = "_SfTopPtDown";
    if(doJESSys) namingConventionFit["JESPlus"] = "_JesUp";
    if(doJESSys) namingConventionFit["JESMinus"] = "_JesDown";
    if(doJERSys) namingConventionFit["JERPlus"] = "_JerUp";
    if(doJERSys) namingConventionFit["JERMinus"] = "_JerDown";
    namingConventionFit[""] = "";

    vector<string> WhatSysts_noJECs;
    
    WhatSysts_noJECs.push_back("iterativefit_lfPlus");   
    WhatSysts_noJECs.push_back("iterativefit_lfMinus");   
    WhatSysts_noJECs.push_back("iterativefit_hfPlus");   
    WhatSysts_noJECs.push_back("iterativefit_hfMinus");   
    WhatSysts_noJECs.push_back("iterativefit_lfstats1Plus");   
    WhatSysts_noJECs.push_back("iterativefit_lfstats1Minus");   
    WhatSysts_noJECs.push_back("iterativefit_lfstats2Plus");   
    WhatSysts_noJECs.push_back("iterativefit_lfstats2Minus");   
    WhatSysts_noJECs.push_back("iterativefit_hfstats1Plus");   
    WhatSysts_noJECs.push_back("iterativefit_hfstats1Minus");   
    WhatSysts_noJECs.push_back("iterativefit_hfstats2Plus");   
    WhatSysts_noJECs.push_back("iterativefit_hfstats2Minus");   
    WhatSysts_noJECs.push_back("iterativefit_cferr1Plus");   
    WhatSysts_noJECs.push_back("iterativefit_cferr1Minus");   
    WhatSysts_noJECs.push_back("iterativefit_cferr2Plus");   
    WhatSysts_noJECs.push_back("iterativefit_cferr2Minus");   
    WhatSysts_noJECs.push_back("pileupPlus");   
    WhatSysts_noJECs.push_back("pileupMinus");   
    WhatSysts_noJECs.push_back("leptonPlus");
    WhatSysts_noJECs.push_back("leptonMinus");
    WhatSysts_noJECs.push_back("TopPtPlus");
    WhatSysts_noJECs.push_back("TopPtMinus");
    if(doInclusive)
    {
        WhatSysts_noJECs.push_back("noSF");
        WhatSysts_noJECs.push_back("OnlyTopPtSF");
        WhatSysts_noJECs.push_back("OnlyBTagSF");
        WhatSysts_noJECs.push_back("OnlyPUSF");
        WhatSysts_noJECs.push_back("OnlyLepSF");
        WhatSysts_noJECs.push_back("OnlyNLOSF");
        WhatSysts_noJECs.push_back("NoBTagSF");
        WhatSysts_noJECs.push_back("NoPUSF");
        WhatSysts_noJECs.push_back("NoLepSF");
        WhatSysts_noJECs.push_back("NoNLOSF");
    }

    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << " - Category: " << category << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();


    cout << " ... Making the TreeProcessor .xml files " << endl;
    system("python scripts/MakeXMLforTreeProcessor.py");
    Double_t CorrectionForAllChannel = 1.; 

    string xmlNom;
    if(Manual_XML) xmlNom = manualxml;
    else if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    else if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    else if(channel == "_All")
    {
        xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
        CorrectionForAllChannel = 1.;
        
    }
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
            //creating three new datsets to contain the tt+ll, tt+cc and tt+ bb compenents.
//            vector<string> ttbar_filenames = datasets[d]->Filenames();
//            cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;
            
//            Dataset* ttbar_ll = new Dataset("TTJets_ll","tt+lf" , true, 633, 1, 2, 1, datasets[d]->Xsection(),ttbar_filenames );
//            Dataset* ttbar_cc = new Dataset("TTJets_cc","tt+cc" , true, 633, 1, 2, 1, datasets[d]->Xsection(), ttbar_filenames );
//            Dataset* ttbar_bb = new Dataset("TTJets_bb","tt+bb" , true, 633, 1, 2, 1, datasets[d]->Xsection(), ttbar_filenames );
            
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
    if(Manual_XML) Luminosity = 1.;
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
    histo1D["LeptonCharge_data_obs"] = new TH1F("LeptonCharge_data_obs","LeptonCharge_data_obs",20,-1,1);
    histo1D["HiggsMass_TOPHLEPBB_hut_data_obs"] = new TH1F("HiggsMass_TOPHLEPBB_hut_data_obs","HiggsMass_TOPHLEPBB_hut_data_obs",20,-1,1);
    histo1D["HiggsBJet2CSVv2_TOPHLEPBB_hut_data_obs"] = new TH1F("HiggsBJet2CSVv2_TOPHLEPBB_hut_data_obs","HiggsBJet2CSVv2_TOPHLEPBB_hut_data_obs",20,-1,1);
    histo1D["MVA_TOPHLEPBB_hut_data_obs"] = new TH1F("MVA_TOPHLEPBB_hut_data_obs","MVA_TOPHLEPBB_hut_data_obs",20,-1,1);

    double CSVv2min = 0.;
    if(baseline_jets == baseline_bjets) CSVv2min = 0.8;

    for(int iSyst = 0; iSyst<WhatSysts.size();iSyst++)
    {
    
        MSPlot[("LeptonCharge"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("LeptonCharge"+WhatSysts[iSyst]).c_str(), 3, -1.5, 1.5, "q_{lep}","Entries", ""); 
        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 30, 0., 250, "M_{H}","Entries", "","GeV");
        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 30, CSVv2min, 1., "b_{H2} CSVv2","Entries", "");
        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 30, -1., 1., "bMVA TopHLepbb","Entries", "");
            
        histo1D[("LeptonCharge_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("LeptonCharge_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("LeptonCharge_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("LeptonCharge_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("HiggsMass_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsMass_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsMass_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsMass_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("HiggsBJet2CSVv2_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("HiggsBJet2CSVv2_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("MVA_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("MVA_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("MVA_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("MVA_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

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
	      }
        else if(dataSetName.find("NLO") != string::npos || dataSetName.find("nlo") !=string::npos || dataSetName.find("amc") !=string::npos) isAMC = true;
        
        for(int JecCounter = WhatSysts_noJECs.size(); JecCounter < WhatSysts.size(); JecCounter++)
        {
            string postfix = "";
            if(isData && WhatSysts[JecCounter] != "") continue;
            if(!isData) postfix = WhatSysts[JecCounter];
	      

		        cout<<"Dataset:  :"<<(dataSetName+WhatSysts[JecCounter]).c_str()<<endl;
		        filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + postfix + ".root";
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
                
//                if(Manual_XML) W_nPV = reweight::LumiReWeighting("MSPlots/MSPlots_All/_19_1_2017/Inclusive/Output.root", "MSPlots/MSPlots_All/_19_1_2017/Inclusive/Output.root", ("MultiSamplePlot_Njets/Njets_"+dataSetName).c_str(), "MultiSamplePlot_Njets/Njets_NP_overlay_HiggsMass_TOPHLEPBB_hut_tHToBB_1L_Kappa_hct");
                if(Manual_XML) W_nPV = reweight::LumiReWeighting(pathPlot.c_str(), "MSPlots/MSPlots_All/_12_1_2017/Inclusive/Output_NPV.root", ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_NP_overlay_HiggsMass_TOPHLEPBB_hut_tHToBB_1L_Kappa_hct");
                else W_nPV = reweight::LumiReWeighting( pathPlot.c_str(), pathPlot.c_str(), ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_Data");
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
            Double_t W_puSF_Minus;
            Double_t W_puSF_Plus;
            Double_t W_fleptonSF;
            Double_t W_fleptonSF_Plus;
            Double_t W_fleptonSF_Minus;
            Double_t W_btagWeight_shape;
            Double_t W_btagWeight_shape_up_lf; 
            Double_t W_btagWeight_shape_down_lf; 
            Double_t W_btagWeight_shape_up_hf; 
            Double_t W_btagWeight_shape_down_hf; 
            Double_t W_btagWeight_shape_up_hfstats1; 
            Double_t W_btagWeight_shape_down_hfstats1; 
            Double_t W_btagWeight_shape_up_hfstats2; 
            Double_t W_btagWeight_shape_down_hfstats2; 
            Double_t W_btagWeight_shape_up_lfstats1; 
            Double_t W_btagWeight_shape_down_lfstats1; 
            Double_t W_btagWeight_shape_up_lfstats2; 
            Double_t W_btagWeight_shape_down_lfstats2; 
            Double_t W_btagWeight_shape_up_cferr1; 
            Double_t W_btagWeight_shape_down_cferr1; 
            Double_t W_btagWeight_shape_up_cferr2; 
            Double_t W_btagWeight_shape_down_cferr2; 
            Double_t W_nloWeight;// for amc@nlo samples
            Double_t W_TopPtReweighing;
          
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
	          Int_t nJets_CSVL; 
	          Int_t nJets_CSVM; 
	          Int_t nJets_CSVT;

            //Variables from bMVA method
            Double_t MVA_TOPHLEPBB_hut;
            Double_t HiggsMass_TOPHLEPBB_hut;
            Double_t HiggsBJet2CSVv2_TOPHLEPBB_hut;
            
            // Weights
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Plus",&W_fleptonSF_Plus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Minus",&W_fleptonSF_Minus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Minus",&W_puSF_Minus);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Plus",&W_puSF_Plus);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape",&W_btagWeight_shape); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_lf",&W_btagWeight_shape_up_lf); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_lf",&W_btagWeight_shape_down_lf); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_hf",&W_btagWeight_shape_up_hf); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_hf",&W_btagWeight_shape_down_hf); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_hfstats1",&W_btagWeight_shape_up_hfstats1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_hfstats1",&W_btagWeight_shape_down_hfstats1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_hfstats2",&W_btagWeight_shape_up_hfstats2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_hfstats2",&W_btagWeight_shape_down_hfstats2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_lfstats1",&W_btagWeight_shape_up_lfstats1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_lfstats1",&W_btagWeight_shape_down_lfstats1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_lfstats2",&W_btagWeight_shape_up_lfstats2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_lfstats2",&W_btagWeight_shape_down_lfstats2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_cferr1",&W_btagWeight_shape_up_cferr1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_cferr1",&W_btagWeight_shape_down_cferr1); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_up_cferr2",&W_btagWeight_shape_up_cferr2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape_down_cferr2",&W_btagWeight_shape_down_cferr2); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_nloWeight",&W_nloWeight); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_TopPtReweighing",&W_TopPtReweighing);  

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
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVL",&nJets_CSVL);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&nJets_CSVM);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVT",&nJets_CSVT);

            ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hut",&MVA_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPHLEPBB_hut",&HiggsMass_TOPHLEPBB_hut);
            ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB_hut",&HiggsBJet2CSVv2_TOPHLEPBB_hut);
	
            double nloSF = 1.;
            int nPos = 0; 
            int nNeg = 0;
            int EntryStart = 0;//Manually overwriting the number of events to run over.
                if(category=="b2j4")
                {
                    EntryStart = (int) 29*nEntries/30;
                    CorrectionForAllChannel = 30;
                }
                else if(category=="b2j3")
                {
                    EntryStart = (int) 9*nEntries/10;
                    CorrectionForAllChannel = 10;
                }
                else if(category=="b3j4")
                {
                    EntryStart = (int) nEntries/2;
                    CorrectionForAllChannel = 2;
                }
                else if(category=="b3j3" || category=="b4j4")
                {
                    EntryStart = 0;
                }
                else if(doInclusive) 
                {
                    EntryStart = (int) 79*nEntries/80;
                    CorrectionForAllChannel = 80;
                }

            if(isAMC && !isData)
            {
                for (int k = EntryStart; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);
		                if(!doInclusive)
		                {
		                    if(nJets_CSVM != baseline_bjets)  continue;

		                    if(baseline_jets == 3 && nJets != baseline_jets) continue;
		                    else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		                }
                    if( W_nloWeight > 0) nPos++;
                    else if( W_nloWeight < 0) nNeg ++;
                }
                nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
            }		

            Double_t average_TopPtWeight = 0.;
            if(dataSetName.find("TTJets") != string::npos)
            {
                int nEventsPassed = 0;
                for (int k = EntryStart; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);
		                if(!doInclusive)
		                {
		                    if(nJets_CSVM != baseline_bjets)  continue;

		                    if(baseline_jets == 3 && nJets != baseline_jets) continue;
		                    else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		                }
                    average_TopPtWeight = average_TopPtWeight + W_TopPtReweighing;
                    nEventsPassed++;
                }
                average_TopPtWeight = average_TopPtWeight/nEventsPassed;
            }
            
            if(isData) EntryStart = 0;

            double EqLumi = datasets[d]->EquivalentLumi();
		
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
                double ScaleFactor = CorrectionForAllChannel;
//                ScaleFactor *= WeightPrivateSignalSample(nJets,Sample->Name());//Reweigh privately produced samples such that number of jets agree
                double W_puSF_applied = 1.;
			          if(!isData)
			          {
			              if(!PVreweighing) W_puSF_applied = W_puSF;
			              else
			              {
			                  W_puSF_applied = W_nPV.ITweight( (int)nvtx );
			              }
			              

                    if(debug)
                    {
                        //Safety triggers in case there are strange things happening in the event weights
                        if(W_fleptonSF < 0 || W_btagWeight_shape < 0 || nloSF < 0 || Luminosity < 0 || W_puSF_applied < 0)
                        {
                              cout << "----- Event " << j << " has a negative weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << "; Luminosity=" << Luminosity << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              cout << "----- The event will be skipped....." << endl;
                              continue;
                        }
                        else if(W_fleptonSF != W_fleptonSF || W_btagWeight_shape != W_btagWeight_shape || nloSF != nloSF || W_puSF_applied != W_puSF_applied)
                        {
                              cout << "----- Event " << j << " has a Nan weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              cout << "----- The event will be skipped....." << endl;
                              continue;
                        }
                        else if(W_fleptonSF >= 40 || W_btagWeight_shape >= 40 || nloSF >= 40 || W_puSF_applied >= 40)
                        {
                              cout << "----- Event " << j << " has a weight larger than 40. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              //cout << "----- The event will be skipped....." << endl;
                              //continue;
                        }
                    }
                    


                    //Nominal scale factor -- scale factors for systematic shifts are calculated below
                    ScaleFactor *= W_puSF_applied;
                    ScaleFactor *= W_fleptonSF;
                    ScaleFactor *= W_btagWeight_shape;
                    ScaleFactor *= nloSF * W_nloWeight;
                }
                else ScaleFactor = 1.;    
		          
      	        //***********************************************FILLING PLOTS**********************************************
                bool ScalePlots = true;
                if(isData) ScalePlots = false;

                if(filepath.find("JESMinus") == string::npos && filepath.find("JESPlus") == string::npos  && filepath.find("JERMinus") == string::npos && filepath.find("JERPlus") == string::npos)
                {
                    for(int iSyst_ = 0; iSyst_ < WhatSysts_noJECs.size(); iSyst_++)
                    {
                    
                        //-----------------------------------------------------------------------------------------------------------
                        // Calculate Scale factors
                        //-----------------------------------------------------------------------------------------------------------
                        SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] = CorrectionForAllChannel;
                    
			                  if(!isData)
			                  {
                            if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtMinus")//Apply no TopPt reweighing
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "OnlyTopPtSF")
                            {
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "OnlyBTagSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "OnlyPUSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "OnlyLepSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "OnlyNLOSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoBTagSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoPUSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoLepSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF* W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoNLOSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                            }
                        }//if(!isData)
                        else SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] = 1.;


                        //-----------------------------------------------------------------------------------------------------------
                        // Fill Plots
                        //-----------------------------------------------------------------------------------------------------------
                        MSPlot[("LeptonCharge"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(LepCharge, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str() ]->Fill(HiggsMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str() ]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str() ]->Fill(MVA_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);//Downscaling signal by 0.1
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("LeptonCharge_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("LeptonCharge_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("LeptonCharge_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("LeptonCharge_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("LeptonCharge_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("LeptonCharge_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("LeptonCharge_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(LepCharge,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("MVA_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * CorrectionForAllChannel / EqLumi);

                    }
               }
                if(filepath.find("JESMinus") != string::npos || filepath.find("JESPlus") != string::npos  || filepath.find("JERMinus") != string::npos || filepath.find("JERPlus") != string::npos || isData || WhatSysts[JecCounter] == "")
               {
                        MSPlot[("LeptonCharge"+WhatSysts[JecCounter]).c_str()]->Fill(LepCharge, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str() ]->Fill(HiggsMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str() ]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str() ]->Fill(MVA_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("LeptonCharge_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("LeptonCharge_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("LeptonCharge_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("LeptonCharge_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("LeptonCharge_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("LeptonCharge_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("LeptonCharge_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("LeptonCharge_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(LepCharge,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("HiggsMass_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("HiggsBJet2CSVv2_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("MVA_TOPHLEPBB_hut_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi * 0.1);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("MVA_TOPHLEPBB_hut_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("MVA_TOPHLEPBB_hut_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVA_TOPHLEPBB_hut,Luminosity * ScaleFactor * CorrectionForAllChannel / EqLumi);

               }
               if(Sample->Name().find("Data") != string::npos) histo1D["LeptonCharge_data_obs"]->Fill(LepCharge);
               if(Sample->Name().find("Data") != string::npos) histo1D["HiggsMass_TOPHLEPBB_hut_data_obs"]->Fill(HiggsMass_TOPHLEPBB_hut);
               if(Sample->Name().find("Data") != string::npos) histo1D["HiggsBJet2CSVv2_TOPHLEPBB_hut_data_obs"]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut);
               if(Sample->Name().find("Data") != string::npos) histo1D["MVA_TOPHLEPBB_hut_data_obs"]->Fill(MVA_TOPHLEPBB_hut);
			                
		        }//for-loop events
		    }//for-loop JEC systematic samples              
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
  if(ApplyPostFit)
  {
      pathPNG += "/PostFit";
      mkdir(pathPNG.c_str(),0777);
  }
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  string outfilename = pathPNG+"/Output.root";

  TFile *outfile = new TFile(outfilename.c_str(),"recreate");
  outfile->cd();

  vector<string> NominalVariableNames;
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
     	string name = it->first;
     	MultiSamplePlot *temp = it->second;
      if (debug)
      {
          cout << "Saving the MSP" << endl;
          cout << " and it->first is " << name << endl;
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
      cout << "Drawing MSP: " << name << endl;
      temp->showNumberEntries(false);
//      if(doInclusive) temp->Draw("MyMSP_"+name, 1, false, false, false, 1);
      bool writePng = false;

      if(name.find("Minus") == string::npos && name.find("Plus")== string::npos 
        && name.find("noSF") == string::npos && name.find("OnlyTopPtSF") == string::npos && name.find("OnlyBTagSF") == string::npos && 
        name.find("OnlyPUSF") == string::npos && name.find("OnlyLepSF") == string::npos && name.find("OnlyNLOSF") == string::npos && 
        name.find("NoTopPtSF") == string::npos && name.find("NoBTagSF") == string::npos && name.find("NoPUSF") == string::npos && 
        name.find("NoLepSF") == string::npos && name.find("NoNLOSF") == string::npos)//Do not save the pictures of the systematics
      {
          NominalVariableNames.push_back(name);
          writePng = true;
      }
//      temp->Write(outfile, name, false,pathPNG, "png");
      temp->Write(outfile, name, false,pathPNG, "eps");
	}

  outfile->Write("kOverwrite");
  outfile->Close();
  
  cout << "  - Making total systematic bands " << endl;
  string errorbandfile = "";
  WhatSysts.pop_back();//Delete the last entry (which should be "") for the systematics plotting
  if(!ApplyPostFit)
  {
      errorbandfile = (pathPNG+"/Systematics_BareHistos.root");
      MakeTotalSystErrorBand_Distributions(outfilename, WhatSysts, datasetnames_backgrounds, NominalVariableNames, errorbandfile);
  }
  else
  {
      errorbandfile = "";//Set the name of the error-band file obtained from the post-fit script.
  }


  //Now remake MSPlots with systematic error bands
  TFile *outfile_errorbands = new TFile((pathPNG+"/Output_withErrorBands.root").c_str(),"recreate");
  outfile_errorbands->cd();

  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
     	string name = it->first;

      if(name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }


     	MultiSamplePlot *temp = it->second;
     	
     	temp->setErrorBandFile(errorbandfile);


     	
      if (debug)
      {
          cout << "Saving the MSP" << endl;
          cout << " and it->first is " << name << endl;
          cout << " Luminosity is " << Luminosity << endl;
      }
      cout << "Drawing MSP: " << name << endl;
      temp->showNumberEntries(false);
      temp->setPreliminary(false);
      temp->setChannel(true,category);
      temp->Draw("MyMSP_"+name, 1, true, true, true, 1);
      bool writePng = false;
//      temp->Write(outfile_errorbands, name, true,pathPNG, "png");
      temp->Write(outfile_errorbands, name, true,pathPNG, "eps");
//      temp->Write(outfile_errorbands, name, true,pathPNG, "pdf");
	}
	outfile_errorbands->Write("kOverwrite");


  //Now store histo's to be used for limit setting
  string outname_PreFit_LeptonCharge = pathPNG+"/input_LeptonCharge_"+category+"_"+SignalSample+".root";

  TFile *outfile_PreFit_LeptonCharge = new TFile(outname_PreFit_LeptonCharge.c_str(),"recreate");
  outfile_PreFit_LeptonCharge->cd();

  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;

      if(name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }
      if(name.find("LeptonCharge") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_PreFit_LeptonCharge->Write("kOverwrite");


  //Now store histo's to be used for limit setting
  string outname_PreFit_HiggsMass_TOPHLEPBB_hut = pathPNG+"/input_HiggsMass_TOPHLEPBB_hut_"+category+"_"+SignalSample+".root";

  TFile *outfile_PreFit_HiggsMass_TOPHLEPBB_hut = new TFile(outname_PreFit_HiggsMass_TOPHLEPBB_hut.c_str(),"recreate");
  outfile_PreFit_HiggsMass_TOPHLEPBB_hut->cd();

  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;

      if(name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }
      if(name.find("HiggsMass_TOPHLEPBB_hut") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_PreFit_HiggsMass_TOPHLEPBB_hut->Write("kOverwrite");



  //Now store histo's to be used for limit setting
  string outname_PreFit_HiggsBJet2CSVv2_TOPHLEPBB_hut = pathPNG+"/input_HiggsBJet2CSVv2_TOPHLEPBB_hut_"+category+"_"+SignalSample+".root";

  TFile *outfile_PreFit_HiggsBJet2CSVv2_TOPHLEPBB_hut = new TFile(outname_PreFit_HiggsBJet2CSVv2_TOPHLEPBB_hut.c_str(),"recreate");
  outfile_PreFit_HiggsBJet2CSVv2_TOPHLEPBB_hut->cd();

  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;

      if(name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }
      if(name.find("HiggsBJet2CSVv2_TOPHLEPBB_hut") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_PreFit_HiggsBJet2CSVv2_TOPHLEPBB_hut->Write("kOverwrite");


  //Now store histo's to be used for limit setting
  string outname_PreFit_MVA_TOPHLEPBB_hut = pathPNG+"/input_MVA_TOPHLEPBB_hut_"+category+"_"+SignalSample+".root";

  TFile *outfile_PreFit_MVA_TOPHLEPBB_hut = new TFile(outname_PreFit_MVA_TOPHLEPBB_hut.c_str(),"recreate");
  outfile_PreFit_MVA_TOPHLEPBB_hut->cd();

  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;

      if(name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }
      if(name.find("MVA_TOPHLEPBB_hut") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_PreFit_MVA_TOPHLEPBB_hut->Write("kOverwrite");

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
    if(Manual_XML) xmlNom = manualxml;
    else if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
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

    MSPlot_nPV["NPV_unw"] = new MultiSamplePlot(datasets, "NPV_unw", 51, -0.5, 50.5, "Number of PV","Entries", ""); 

  
 

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


void MakeTotalSystErrorBand_Distributions(string outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile)
{

    TFile *MSPlotFile = new TFile(outfilename.c_str(),"read");

    map<string,MultiSamplePlot*> MSPlot_ErrorBands;
    map<string,TH1F*> histo1D_nominal;
    map<string,TH1F*> histo1D_Up_SamplesAdded;
    map<string,TH1F*> histo1D_Down_SamplesAdded;

    map<string,TH1F*> histo1D_TotalUp;
    map<string,TH1F*> histo1D_TotalDown;

    //Define rate uncertainties
    Double_t LumiUncPlus = 0.062;
    Double_t LumiUncMinus = 0.062;
    Double_t XSecTTJetPlus = 0.055;
    Double_t XSecTTJetMinus = 0.055;
    Double_t XSecOtherPlus = 0.1;
    Double_t XSecOtherMinus = 0.1;

    for(int iVar = 0; iVar < NominalVariableNames.size(); iVar++)
    {

        cout << "  - MakeTotalSystErrorBand_Distributions(): Variable " << NominalVariableNames[iVar] << endl;

        histo1D_nominal[NominalVariableNames[iVar].c_str()] = 0;
        histo1D_TotalDown[NominalVariableNames[iVar].c_str()] =  0;
        histo1D_TotalUp[NominalVariableNames[iVar].c_str()] =  0;


        //Add the nominal samples into 1 histogram
        for(int iDataName = 0; iDataName < datasetNames.size(); iDataName++)
        {
            TDirectory *subdir_nominal = (TDirectory*) MSPlotFile->Get(("MultiSamplePlot_"+NominalVariableNames[iVar]).c_str());
            subdir_nominal->cd();
            string nominalname = (NominalVariableNames[iVar]+"_"+datasetNames[iDataName]+"_");
            
            TH1F *h_tmp =  (TH1F*)subdir_nominal->Get(nominalname.c_str());
//cout << "h_tmp->GetBinContent(h_tmp->GetNbinsX()): " << h_tmp->GetBinContent(h_tmp->GetNbinsX()) << endl;
/*
            //making sure that the overflow is transferred to the last 'visible' bin; analogously for underflow...
            TH1F* h_tmp = (TH1F*) h_tmp_->Clone();
            int Nbins_ = h_tmp->GetNbinsX();
            h_tmp->SetBinContent(Nbins_,h_tmp->GetBinContent(Nbins_)+h_tmp->GetBinContent(Nbins_+1));
            h_tmp->SetBinContent(Nbins_+1,0);
            h_tmp->SetBinContent(1,h_tmp->GetBinContent(0)+h_tmp->GetBinContent(1));
            h_tmp->SetBinContent(0,0);
*/
            TH1F* h_tmp__scaleup = (TH1F*) h_tmp->Clone();//Make a new tmp which will be scaled according to the cross section uncertainty 
            TH1F* h_tmp__scaledown = (TH1F*) h_tmp->Clone();//Make a new tmp which will be scaled according to the cross section uncertainty



            if(datasetNames[iDataName].find("TTJets")!= string::npos)
            {
                h_tmp__scaleup->Scale(1+XSecTTJetPlus);
                h_tmp__scaledown->Scale(1-XSecTTJetMinus);
            }
            else
            {
                h_tmp__scaleup->Scale(1+XSecOtherPlus);
                h_tmp__scaledown->Scale(1-XSecOtherMinus);
            }

            if(iDataName == 0)
            {
                histo1D_nominal[NominalVariableNames[iVar].c_str()] = (TH1F*) h_tmp->Clone();
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()] = (TH1F*) h_tmp->Clone();//PlusPlus in the object name due to convention down below
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()] = (TH1F*) h_tmp->Clone();

                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncPlusPlus").c_str()] = (TH1F*) h_tmp__scaleup->Clone();//PlusPlus in the object name due to convention down below
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncMinusMinus").c_str()] = (TH1F*) h_tmp__scaledown->Clone();
            }
            else
            {
                histo1D_nominal[NominalVariableNames[iVar].c_str()]->Add(h_tmp);
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()]->Add(h_tmp);
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()]->Add(h_tmp);
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncPlusPlus").c_str()]->Add(h_tmp__scaleup);
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncMinusMinus").c_str()]->Add(h_tmp__scaledown);
            }
        }

        //Scale the  Lumi Uncertainties
        histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()]->Scale(1+LumiUncPlus);
        histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()]->Scale(1-LumiUncMinus);
        
            
        //Add the systematic samples for up-variation into 1 histogram and down-variation
        for(int iSys = 0; iSys < systematics.size(); iSys++)
        {
            string varNameSys = NominalVariableNames[iVar]+systematics[iSys];
            TDirectory *subdir_sys = (TDirectory*) MSPlotFile->Get(("MultiSamplePlot_"+varNameSys).c_str());
            subdir_sys->cd();
            cout << "  - MakeTotalSystErrorBand_Distributions(): Making systematic band for " << systematics[iSys] << endl;

            histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()] = 0;
            histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()] = 0;

            for(int iDataName = 0; iDataName < datasetNames.size(); iDataName++)
            {

                TH1F *h_tmp =  (TH1F*)subdir_sys->Get((varNameSys+"_"+datasetNames[iDataName]+"_").c_str());
//cout << "h_tmp->GetBinContent(h_tmp->GetNbinsX()): " << h_tmp->GetBinContent(h_tmp->GetNbinsX()) << endl;
                
                if(systematics[iSys].find("Plus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()] = (TH1F*) h_tmp->Clone();
                    else histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()]->Add(h_tmp);
                }
                else if(systematics[iSys].find("Minus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()] = (TH1F*) h_tmp->Clone();
                    else histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()]->Add(h_tmp);
                }
            }
        }
            

        // Add the rate systematic undertainties to the systematics list, as well as the MC statistical uncertainty
        systematics.push_back("LumiUncPlus");
        systematics.push_back("LumiUncMinus");
        systematics.push_back("XSecUncPlus");
        systematics.push_back("XSecUncMinus");
       
        //Run over all systematics to add their effect in each bin in quadrature.
        int nBins = histo1D_nominal[NominalVariableNames[iVar].c_str()]->GetNbinsX();
        
        //Initialize the total uncertainty histograms
        histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()] = (TH1F*) histo1D_nominal[NominalVariableNames[iVar].c_str()]->Clone(/*(NominalVariableNames[iVar]+"Plus").c_str()*/);
        histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()] = (TH1F*) histo1D_nominal[NominalVariableNames[iVar].c_str()]->Clone(/*(NominalVariableNames[iVar]+"Minus").c_str()*/);

        for(int iBin = 0; iBin < nBins+1; iBin++)
        {
            float bincontent_nominal  = histo1D_nominal[NominalVariableNames[iVar].c_str()]->GetBinContent(iBin);
            float bincontent_up_squared = 0.;
            float bincontent_down_squared = 0.;
            for(int iSys = 0; iSys < systematics.size(); iSys++)
            {
                string varNameSys = NominalVariableNames[iVar]+systematics[iSys];

                float bincontent_Syst_vs_Nom = 0.;
                
                if(systematics[iSys].find("Minus")!= string::npos)  bincontent_Syst_vs_Nom = bincontent_nominal - histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()]->GetBinContent(iBin);
                else if(systematics[iSys].find("Plus")!= string::npos) bincontent_Syst_vs_Nom = bincontent_nominal - histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()]->GetBinContent(iBin);

                //Check whether the variation goes up or down wrt to the nominal and add in quadrature the contents to the relevant histo.
                if(bincontent_Syst_vs_Nom > 0.) bincontent_down_squared += bincontent_Syst_vs_Nom*bincontent_Syst_vs_Nom;
                else if(bincontent_Syst_vs_Nom < 0.) bincontent_up_squared += bincontent_Syst_vs_Nom*bincontent_Syst_vs_Nom;
            }
            
//if(iBin == nBins) cout << "Content of last bin for Plus: " << bincontent_nominal + sqrt(bincontent_up_squared + bincontent_nominal) << endl;
//if(iBin == nBins) cout << "Content of last bin for Minus: " << bincontent_nominal  - sqrt(bincontent_down_squared + bincontent_nominal) << endl;
//if(iBin == nBins) cout << "Content of last bin for Nominal: " << bincontent_nominal << endl;

            histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()]->SetBinContent(iBin,bincontent_nominal + sqrt(bincontent_up_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
            histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()]->SetBinContent(iBin,bincontent_nominal  - sqrt(bincontent_down_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
        }
        
        //Delete the last entries for the systematics that were added in MakeTotalSystErrorBand_Distributions
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
    }


    //Write output histos to a file  
    TFile *fout = new TFile(outputFile.c_str(),"recreate");
    for(int iVar = 0; iVar < NominalVariableNames.size(); iVar++)
    {
        fout->cd();
        TDirectory* subdir = fout->mkdir(("MultiSamplePlot_"+NominalVariableNames[iVar]).c_str());
        subdir->cd();
        
        //Write the histos according to the definitions from MultiSamplePlot to read the systematics
        histo1D_nominal[NominalVariableNames[iVar].c_str()]->Write("Nominal");
        histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()]->Write("Plus");
        histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()]->Write("Minus");
        
        subdir->Write("kOverwrite");
        subdir->Close();
        delete subdir;

    }
    fout->Write("kOverwrite");
}


double WeightPrivateSignalSample(Int_t n_jets, string samplename)
{
    double weight = 1.;
    if(samplename.find("Private")== string::npos) weight = 1.;
    else
    {
        if(samplename.find("")!= string::npos)
        {
            if (n_jets == 0) weight = 0;
            else if (n_jets == 1) weight = 0;
            else if (n_jets == 2) weight = 0;
            else if (n_jets == 3) weight = 1.094071257;
            else if (n_jets == 4) weight = 0.956789131;
            else if (n_jets == 5) weight = 0.836715582;
            else if (n_jets == 6) weight = 0.719082739;
            else if (n_jets == 7) weight = 0.641676573;
            else if (n_jets == 8) weight = 0.549510742;
            else if (n_jets == 9) weight = 0.457332313;
            else if (n_jets == 10) weight = 0.715143217;
        }
    }
    
//cout << "JetWeight: " << weight << endl;
//weight = 1;
    return weight;
}
