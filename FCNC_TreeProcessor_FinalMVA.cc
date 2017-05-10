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
bool PrivateSampleTraining = false;

// functions prototype
string intToStr (int number);
void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug);
void MakeTotalSystErrorBand_Distributions(string directory, string category, string coupling, string outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile);
double WeightPrivateSignalSample(Int_t n_jets, string samplename);
double OptimalCut_CombTraining(string category, string coupling);
double PostFitScaleFactor(string category, string coupling, string samplename);

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
        cout << "    bool doJESSys  = strtol(argv[7], NULL,10);" << endl;
        cout << "    bool doJERSys  = strtol(argv[8], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[9], NULL,10);" << endl;
        cout << "    bool ApplyPostFit         =strtol(argv[10], NULL,10);" << endl;
        cout << "    int khut         =strtod(argv[11], NULL,10);" << endl;
        cout << "    int khct         =strtod(argv[12], NULL,10);" << endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    string SignalSample  = argv[3];//Valid arguments are: hut,  hct, 2D
    string channel            = argv[4];
    string date            = argv[5];
    bool PVreweighing = strtol(argv[6], NULL,10);
    bool doJESSys  = strtol(argv[7], NULL,10);
    bool doJERSys  = strtol(argv[8], NULL,10);
    bool debug         =strtol(argv[9], NULL,10);
    bool ApplyPostFit         =strtol(argv[10], NULL,10);
    int khut         =strtol(argv[11], NULL,10); //Divide this number by 100
    int khct         =strtol(argv[12], NULL,10);

    bool split_ttbar = true;   

   string coupling_hut = "khut0p" + intToStr(khut) + "_";
   if(khut < 10) coupling_hut = "khut0p0" + intToStr(khut) + "_";
   string coupling_hct = "khct0p" + intToStr(khct) + "_";
   if(khct < 10) coupling_hct = "khct0p0" + intToStr(khct) + "_";

    
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
//    string TrainingName = "Training_" + SignalSample + channel + "_" +  category;//Example: Training_SThut_El_b3j3

    string TrainingName = "";
    if(khut == 0 && khct == 0) TrainingName = "Training_" + SignalSample + channel + "_" +  category;
    else TrainingName = "Training_" + SignalSample + "_" + coupling_hut + coupling_hct + channel + "_" +  category;

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
/*    WhatSysts.push_back("noSF");SignalSample
    WhatSysts.push_back("OnlyTopPtSF");
    WhatSysts.push_back("OnlyBTagSF");
    WhatSysts.push_back("OnlyPUSF");
//    WhatSysts.push_back("OnlyLepSF");
//    WhatSysts.push_back("OnlyNLOSF");
    WhatSysts.push_back("NoBTagSF");
    WhatSysts.push_back("NoPUSF");
//    WhatSysts.push_back("NoLepSF");
//    WhatSysts.push_back("NoNLOSF");
*/    if(doJESSys) WhatSysts.push_back("JESPlus");
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
/*    WhatSysts_noJECs.push_back("noSF");
    WhatSysts_noJECs.push_back("OnlyTopPtSF");
    WhatSysts_noJECs.push_back("OnlyBTagSF");
    WhatSysts_noJECs.push_back("OnlyPUSF");
//    WhatSysts_noJECs.push_back("OnlyLepSF");
//    WhatSysts_noJECs.push_back("OnlyNLOSF");
    WhatSysts_noJECs.push_back("NoBTagSF");
    WhatSysts_noJECs.push_back("NoPUSF");
//    WhatSysts_noJECs.push_back("NoLepSF");
//    WhatSysts_noJECs.push_back("NoNLOSF");
*/

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

    int nEvents_forStatisticalUnc_ttbb = 0;
    int nEvents_forStatisticalUnc_ttcc = 0;
    int nEvents_forStatisticalUnc_ttlf = 0;
    int nEvents_forStatisticalUnc_other = 0;
    int nEvents_forStatisticalUnc_Signal_SThct = 0;
    int nEvents_forStatisticalUnc_Signal_TThct = 0;
    int nEvents_forStatisticalUnc_Signal_SThut = 0;
    int nEvents_forStatisticalUnc_Signal_TThut = 0;


    histo1D["maxSTandTT_data_obs"] = new TH1F("maxSTandTT_data_obs","maxSTandTT_data_obs",20,-1,1);
    histo1D["ST_data_obs"] = new TH1F("ST_data_obs","ST_data_obs",20,-1,1);
    histo1D["TT_data_obs"] = new TH1F("TT_data_obs","TT_data_obs",20,-1,1);
    histo1D["combSTandTT_data_obs"] = new TH1F("combSTandTT_data_obs","combSTandTT_data_obs",20,-1,1);
    histo1D["combSTandTT_cutCount_data_obs"] = new TH1F("combSTandTT_cutCount_data_obs","combSTandTT_cutCount_data_obs",1,-1,1);
    
    string xaxislabelcoupling = " #kappa_{";
    if(SignalSample == "hut") xaxislabelcoupling +=  "Hut}";
    else if(SignalSample == "hct") xaxislabelcoupling +=  "Hct}";
    
    for(int iSyst = 0; iSyst<WhatSysts.size();iSyst++)
    {
        MSPlot[("MVA_MaxTT-ST_"+TrainingName+WhatSysts[iSyst]).c_str()] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_MaxTT-ST_"+TrainingName+WhatSysts[iSyst]).c_str(), 20, -1., 1., ("BDT discriminator"+xaxislabelcoupling),"Events", "");
        MSPlot[("MVA_CombTT-ST_"+TrainingName+WhatSysts[iSyst]).c_str()] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_CombTT-ST_"+TrainingName+WhatSysts[iSyst]).c_str(), 20, -1., 1., ("BDT discriminator"+xaxislabelcoupling),"Events", "");

        MSPlot[("MVA_ST"+TrainingName+WhatSysts[iSyst]).c_str()] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_ST"+TrainingName+WhatSysts[iSyst]).c_str(), 20, -1., 1., ("ST BDT disc."+xaxislabelcoupling),"Events", "");
        MSPlot[("MVA_TT"+TrainingName+WhatSysts[iSyst]).c_str()] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TT"+TrainingName+WhatSysts[iSyst]).c_str(), 20, -1., 1., ("TT BDT disc."+xaxislabelcoupling),"Events", "");

        histo1D[("maxSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("maxSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("maxSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("maxSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("ST_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("ST_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("ST_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("ST_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("TT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("TT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("TT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("TT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("combSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);
        histo1D[("combSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),20,-1,1);

        histo1D[("combSTandTT_cutCount_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_sig"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_sig_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_sig_ttbar"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_stop"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_zjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_wjets"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
        histo1D[("combSTandTT_cutCount_other"+namingConventionFit[WhatSysts[iSyst]]).c_str()] = new TH1F(("combSTandTT_cutCount_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),("combSTandTT_cutCount_other"+namingConventionFit[WhatSysts[iSyst]]).c_str(),1,-1,1);
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

        if(dataSetName.find("Private") != string::npos) continue;//Do not read out on private signal samples
        if(khut == 0 && khct == 0 && dataSetName.find("NP_overlay") != string::npos && dataSetName.find(SignalSample.c_str()) == string::npos) continue;//Do not read out on signal samples of the wrong coupling


        double SignalWeight = 1.;
        if(dataSetName.find("NP_overlay") != string::npos && dataSetName.find("hut") != string::npos) SignalWeight = (double(khut)/100)*(double(khut)/100);
        else if(dataSetName.find("NP_overlay") != string::npos && dataSetName.find("hct") != string::npos) SignalWeight = (double(khct)/100)*(double(khct)/100);
        if(khut == 0 && khct == 0 ) SignalWeight = 0.1;
        
        
        for(int JecCounter = WhatSysts_noJECs.size(); JecCounter < WhatSysts.size(); JecCounter++)
        {
            string postfix = "";
            if(isData && WhatSysts[JecCounter] != "") continue;
            if(!isData) postfix = WhatSysts[JecCounter];
	      

		        cout<<"Dataset:  :"<<(dataSetName+WhatSysts[JecCounter]).c_str()<<endl;
		        filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + postfix + ".root";
//if(filepath.find("ttWJERPlus")!=string::npos) continue;
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
	          std::string weightsFile_combSTandTT= "weights/Comb"+TrainingName+"_BDT.weights.xml";
            if(SignalSample == "2D")
            {
                weightsFile_ST = "weights/Training_SThut" + channel + "_" +  category+"_BDT.weights.xml";
                weightsFile_TT= "weights/Training_TThut" + channel + "_" +  category+"_BDT.weights.xml";
            }
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
                
//                if(Manual_XML) W_nPV = reweight::LumiReWeighting("MSPlots/MSPlots_All/_19_1_2017/Inclusive/Output.root", "MSPlots/MSPlots_All/_19_1_2017/Inclusive/Output.root", ("MultiSamplePlot_Njets/Njets_"+dataSetName).c_str(), "MultiSamplePlot_Njets/Njets_NP_overlay_ST_tHToBB_1L_Kappa_hct");
                if(Manual_XML) W_nPV = reweight::LumiReWeighting(pathPlot.c_str(), "MSPlots/MSPlots_All/_12_1_2017/Inclusive/Output_NPV.root", ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_NP_overlay_ST_tHToBB_1L_Kappa_hct");
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
            Double_t W_btagWeight_CSVv2M_mujets_central;
            Double_t W_btagWeight_CSVv2M_mujets_up;
            Double_t W_btagWeight_CSVv2M_mujets_down;
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
/*            Double_t W_weight0;
            Double_t W_weight1;
            Double_t W_weight2;
            Double_t W_weight3;
            Double_t W_weight4;
            Double_t W_weight5;
            Double_t W_weight6;
            Double_t W_weight7;
            Double_t W_weight8; 
            Double_t W_hdamp_up; 
            Double_t W_hdamp_dw; 
*/            Double_t W_TopPtReweighing;
          
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
/*            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight0",&W_weight0);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight1",&W_weight1);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight2",&W_weight2);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight3",&W_weight3); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight4",&W_weight4);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight5",&W_weight5); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight6",&W_weight6);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight7",&W_weight7); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight8",&W_weight8);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_hdamp_up",&W_hdamp_up);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_hdamp_dw",&W_hdamp_dw);  
*/            ttree[(dataSetName).c_str()]->SetBranchAddress("W_TopPtReweighing",&W_TopPtReweighing);  

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
                      

            int EntryStart = 0;
            double Doubling = 1.;
            if(PrivateSampleTraining)
            {
                if(!isData && dataSetName.find("NP_") == string::npos)
                {
                    EntryStart = (int) nEntries/2+1;//Only read in the other half of simulation events on which the BDT is not trained.
                    Doubling = 2;
                }
            }
            else if(!isData)
            {
                EntryStart = (int) nEntries/2+1;
                Doubling = 2;
            }
//            EntryStart = nEntries/2+nEntries/3+1;//Manually overwriting the number of events to run over.
//            Doubling = 1;

            double nloSF = 1.;
            int nPos = 0; 
            int nNeg = 0;
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
                if(nPos + nNeg != 0) nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
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
                if(nEventsPassed != 0)
                {
                    average_TopPtWeight = average_TopPtWeight/nEventsPassed;
                }
            }
            
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

                if(filepath.find("JESMinus") == string::npos && filepath.find("JESPlus") == string::npos  && filepath.find("JERMinus") == string::npos && filepath.find("JERPlus") == string::npos)
                {
                    if(Sample->Name().find("TTJets_bb") != string::npos) nEvents_forStatisticalUnc_ttbb++;
                    else if(Sample->Name().find("TTJets_cc") != string::npos) nEvents_forStatisticalUnc_ttcc++;
                    else if(Sample->Name().find("TTJets_ll") != string::npos) nEvents_forStatisticalUnc_ttlf++;
                    else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos)  nEvents_forStatisticalUnc_other++;
                    else if(Sample->Name().find("hct") != string::npos)
                    {
                            if(Sample->Name().find("NP_overlay_ST") != string::npos) nEvents_forStatisticalUnc_Signal_SThct++;
                            else if(Sample->Name().find("NP_overlay_TT") != string::npos) nEvents_forStatisticalUnc_Signal_TThct++;
                    }
                    else if(Sample->Name().find("hut"))
                    {
                            if(Sample->Name().find("NP_overlay_ST") != string::npos) nEvents_forStatisticalUnc_Signal_SThut++;
                            else if(Sample->Name().find("NP_overlay_TT") != string::npos) nEvents_forStatisticalUnc_Signal_TThut++;
                    }
                }
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

//                    if(debug)
//                    {
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
                        else if(W_fleptonSF >= 50 || W_btagWeight_shape >= 50 || nloSF >= 40 || W_puSF_applied >= 50)
                        {
                              cout << "----- Event " << j << " has a weight larger than 40. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                              cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                              //cout << "----- The event will be skipped....." << endl;
                              //continue;
                        }
//                    }                


                    //Nominal scale factor -- scale factors for systematic shifts are calculated below
                    ScaleFactor *= W_puSF_applied;
                    ScaleFactor *= W_fleptonSF;
                    ScaleFactor *= W_btagWeight_shape;
                    ScaleFactor *= nloSF * W_nloWeight;
                }
                else ScaleFactor = 1.;    
		          
      	        //***********************************************FILLING PLOTS**********************************************
	              if( HiggsMass_TOPHLEPBB_hut > 500. ) HiggsMass_TOPHLEPBB_hut = 500.;
	              if( TopLepMass_TOPHLEPBB_hut > 500. ) TopLepMass_TOPHLEPBB_hut = 500.;
	              if( TopLepPt_TOPHLEPBB_hut > 1000. ) TopLepPt_TOPHLEPBB_hut = 1000.;
	              if( HiggsMass_TOPHLEPBB_hct > 500. ) HiggsMass_TOPHLEPBB_hct = 500.;
	              if( TopLepMass_TOPHLEPBB_hct > 500. ) TopLepMass_TOPHLEPBB_hct = 500.;
	              if( TopLepPt_TOPHLEPBB_hct > 1000. ) TopLepPt_TOPHLEPBB_hct = 1000.;
	              if( TopLepMass_TOPTOPLEPHAD > 500. || TopLepMass_TOPTOPLEPHAD != TopLepMass_TOPTOPLEPHAD) TopLepMass_TOPTOPLEPHAD = 500.;
	              if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
                if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
                

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


                bool ScalePlots = true;
                if(isData) ScalePlots = false;
                if(ApplyPostFit) Doubling *= PostFitScaleFactor(category,SignalSample,Sample->Name()); 

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
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtMinus")//Apply no TopPt reweighing
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
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
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoBTagSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoPUSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoLepSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
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
                        MSPlot[("MVA_MaxTT-ST_"+TrainingName+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVAvalue_maxSTandTT, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling); 
                        MSPlot[("MVA_CombTT-ST_"+TrainingName+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVAvalue_combSTandTT, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling); 

                        MSPlot[("MVA_ST"+TrainingName+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVAvalue_ST, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling); 
                        MSPlot[("MVA_TT"+TrainingName+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVAvalue_TT, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling); 
                        
                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);//Downscaling signal by 0.1
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("maxSTandTT_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("maxSTandTT_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("maxSTandTT_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("maxSTandTT_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("maxSTandTT_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("maxSTandTT_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("maxSTandTT_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("ST_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("ST_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("ST_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("ST_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("ST_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("ST_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("ST_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_ST,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("TT_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("TT_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("TT_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("TT_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("TT_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("TT_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("TT_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_TT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("combSTandTT_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("combSTandTT_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("combSTandTT_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("combSTandTT_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);

                        if(MVAvalue_combSTandTT > OptimalCut_CombTraining(category, SignalSample))
                        {
                            if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                            if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                            else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig_ttbar"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi * SignalWeight);
                            else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("ST-") != string::npos) histo1D[("combSTandTT_cutCount_stop"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("DYJets") != string::npos) histo1D[("combSTandTT_cutCount_zjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("WJets") != string::npos) histo1D[("combSTandTT_cutCount_wjets"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                            else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("combSTandTT_cutCount_other"+namingConventionFit[WhatSysts_noJECs[iSyst_]]).c_str()]->Fill(0.,Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] * Doubling / EqLumi);
                        }
                    }
               }
               if(filepath.find("JESMinus") != string::npos || filepath.find("JESPlus") != string::npos  || filepath.find("JERMinus") != string::npos || filepath.find("JERPlus") != string::npos || isData || WhatSysts[JecCounter] == "")
               {
                        MSPlot[("MVA_MaxTT-ST_"+TrainingName+WhatSysts[JecCounter]).c_str()]->Fill(MVAvalue_maxSTandTT, Sample, ScalePlots, Luminosity * ScaleFactor * Doubling); //Factor 2 to compensate for the fact we're running over half the number of simulated events
                        MSPlot[("MVA_CombTT-ST_"+TrainingName+WhatSysts[JecCounter]).c_str()]->Fill(MVAvalue_combSTandTT, Sample, ScalePlots, Luminosity * ScaleFactor * Doubling); //Factor 2 to compensate for the fact we're running over half the number of simulated events

                        MSPlot[("MVA_ST"+TrainingName+WhatSysts[JecCounter]).c_str()]->Fill(MVAvalue_ST, Sample, ScalePlots, Luminosity * ScaleFactor * Doubling); //Factor 2 to compensate for the fact we're running over half the number of simulated events
                        MSPlot[("MVA_TT"+TrainingName+WhatSysts[JecCounter]).c_str()]->Fill(MVAvalue_TT, Sample, ScalePlots, Luminosity * ScaleFactor * Doubling); //Factor 2 to compensate for the fact we're running over half the number of simulated events

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("maxSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("maxSTandTT_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("maxSTandTT_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("maxSTandTT_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("maxSTandTT_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("maxSTandTT_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("maxSTandTT_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("maxSTandTT_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_maxSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("ST_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("ST_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("ST_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("ST_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("ST_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("ST_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("ST_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("ST_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_ST,Luminosity * ScaleFactor * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("TT_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("TT_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("TT_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("TT_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("TT_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("TT_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("TT_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("TT_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_TT,Luminosity * ScaleFactor * Doubling / EqLumi);

                        if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                        else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("ST-") != string::npos) histo1D[("combSTandTT_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("DYJets") != string::npos) histo1D[("combSTandTT_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("WJets") != string::npos) histo1D[("combSTandTT_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);
                        else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("combSTandTT_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(MVAvalue_combSTandTT,Luminosity * ScaleFactor * Doubling / EqLumi);

                        if(MVAvalue_combSTandTT > OptimalCut_CombTraining(category, SignalSample))
                        {
                            if(Sample->Name().find("NP_") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                            if(Sample->Name().find("NP_overlay_ST") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                            else if(Sample->Name().find("NP_overlay_TT") != string::npos && dataSetName.find(SignalSample.c_str()) != string::npos) histo1D[("combSTandTT_cutCount_sig_ttbar"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi * SignalWeight);
                            else if(Sample->Name().find("TTJets_bb") != string::npos) histo1D[("combSTandTT_cutCount_ttbb"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("TTJets_cc") != string::npos) histo1D[("combSTandTT_cutCount_ttcc"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("TTJets_ll") != string::npos) histo1D[("combSTandTT_cutCount_ttlf"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("ST-") != string::npos) histo1D[("combSTandTT_cutCount_stop"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("DYJets") != string::npos) histo1D[("combSTandTT_cutCount_zjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("WJets") != string::npos) histo1D[("combSTandTT_cutCount_wjets"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                            else if(Sample->Name().find("Data") == string::npos && Sample->Name().find("NP_overlay") == string::npos) histo1D[("combSTandTT_cutCount_other"+namingConventionFit[WhatSysts[JecCounter]]).c_str()]->Fill(0.,Luminosity * ScaleFactor * Doubling / EqLumi);
                        }
               }
               if(Sample->Name().find("Data") != string::npos) histo1D["maxSTandTT_data_obs"]->Fill(MVAvalue_maxSTandTT);
               if(Sample->Name().find("Data") != string::npos) histo1D["ST_data_obs"]->Fill(MVAvalue_ST);
               if(Sample->Name().find("Data") != string::npos) histo1D["TT_data_obs"]->Fill(MVAvalue_TT);
               if(Sample->Name().find("Data") != string::npos) histo1D["combSTandTT_data_obs"]->Fill(MVAvalue_combSTandTT);
               if(Sample->Name().find("Data") != string::npos && MVAvalue_combSTandTT > OptimalCut_CombTraining(category, SignalSample)) histo1D["combSTandTT_cutCount_data_obs"]->Fill(0.);
			                
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


  if(khut ==0 && khct == 0)
  {
      string outfilename = pathPNG+"/OutputMVA_"+SignalSample+".root";

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
    //      temp->Draw("MyMSP_"+name, 1, false, false, false, 1);

          if(name.find("Minus") == string::npos && name.find("Plus")== string::npos 
            && name.find("noSF") == string::npos && name.find("OnlyTopPtSF") == string::npos && name.find("OnlyBTagSF") == string::npos && 
            name.find("OnlyPUSF") == string::npos && name.find("OnlyLepSF") == string::npos && name.find("OnlyNLOSF") == string::npos && 
            name.find("NoTopPtSF") == string::npos && name.find("NoBTagSF") == string::npos && name.find("NoPUSF") == string::npos && 
            name.find("NoLepSF") == string::npos && name.find("NoNLOSF") == string::npos)//Do not save the pictures of the systematics
          {
              NominalVariableNames.push_back(name);
          }
          temp->Write(outfile, name, false,pathPNG, "png");
	    }

      outfile->Write("kOverwrite");
      outfile->Close();
      
      cout << "  - Making total systematic bands " << endl;
      string errorbandfile = "";
      WhatSysts.pop_back();//Delete the last entry (which should be "") for the systematics plotting
      if(!ApplyPostFit && khut == 0 && khct == 0)
      {
          errorbandfile = (pathPNG+"/Systematics_BareHistosMVA"+SignalSample+".root");
          MakeTotalSystErrorBand_Distributions(pathPNG, category, TrainingName, outfilename, WhatSysts, datasetnames_backgrounds, NominalVariableNames, errorbandfile);
      }
      else
      {
          errorbandfile = (pathPNG+"/PostFitSystematics_BareHistosMVA"+SignalSample+".root");//Set the name of the error-band file obtained from the post-fit script.
      }



      //Now remake MSPlots with systematic error bands
      TFile *outfile_errorbands = new TFile((pathPNG+"/Output_withErrorBandsMVA_"+TrainingName+".root").c_str(),"recreate");
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
          
          if(ApplyPostFit && name.find("Comb") == string::npos) continue;

         	MultiSamplePlot *temp = it->second;
         	
         	temp->setErrorBandFile(errorbandfile);


          if(name.find("CategoryRates") != string::npos)
          {
              vector<string> label;
              label.push_back("(nj=3,nb=2)");
              label.push_back("(nj>3,nb=2)");
              label.push_back("(nj=3,nb=3)");
              label.push_back("(nj>3,nb=3)");
              label.push_back("(nj>3,nb=4)");
              temp->setBins(label);

              temp->showNumberEntries(false);
              temp->setChannel(true,category);
              temp->Draw("MyMSP_"+name, 1, true, true, true, 1);
              temp->Write(outfile_errorbands, name, true,pathPNG, "eps");

          }
          else
          {
               	
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
              temp->Write(outfile_errorbands, name, true,pathPNG, "eps");//You can only call 1 format for saving the plots. The second time you want to draw, the THStacks are empty, because the object has been written to the root-file.
        //      temp->Write(outfile_errorbands, name, true,pathPNG, "png");
        //      temp->Write(outfile_errorbands, name, true,pathPNG, "pdf");//.pdf files are corrputed due to the #backslash symbol defined in MultiSamplePlot.cc
         }
	    }
	    outfile_errorbands->Write("kOverwrite");
  }


  //Now store histo's to be used for limit setting
  string outname_limitsetting_maxSTandTT = pathPNG+"/input_MVA";
  if(SignalSample == "hct") outname_limitsetting_maxSTandTT += "HctMAX_"+category+"_"+SignalSample+".root";
  else if(SignalSample == "hut") outname_limitsetting_maxSTandTT += "HutMAX_"+category+"_"+SignalSample+".root";

  TFile *outfile_limitsetting_maxSTandTT = new TFile(outname_limitsetting_maxSTandTT.c_str(),"recreate");
  outfile_limitsetting_maxSTandTT->cd();

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
      if(name.find("maxSTandTT") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_limitsetting_maxSTandTT->Write("kOverwrite");


  //Now store histo's to be used for limit setting
  string outname_limitsetting_ST = pathPNG+"/input_MVA";
  if(SignalSample == "hct") outname_limitsetting_ST += "HctST_"+category+"_"+SignalSample+".root";
  else if(SignalSample == "hut") outname_limitsetting_ST += "HutST_"+category+"_"+SignalSample+".root";

  TFile *outfile_limitsetting_ST = new TFile(outname_limitsetting_ST.c_str(),"recreate");
  outfile_limitsetting_ST->cd();

  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;

      if(name.find("TT") != string::npos || name.find("combSTandTT") != string::npos || name.find("maxSTandTT") != string::npos || name.find("Minus") != string::npos || name.find("Plus")!= string::npos 
        || name.find("noSF") != string::npos || name.find("OnlyTopPtSF") != string::npos || name.find("OnlyBTagSF") != string::npos || 
        name.find("OnlyPUSF") != string::npos || name.find("OnlyLepSF") != string::npos || name.find("OnlyNLOSF") != string::npos || 
        name.find("NoTopPtSF") != string::npos || name.find("NoBTagSF") != string::npos || name.find("NoPUSF") != string::npos || 
        name.find("NoLepSF") != string::npos || name.find("NoNLOSF") != string::npos)//Do not save the pictures of the systematics
      {
          continue;
      }


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_limitsetting_ST->Write("kOverwrite");

  //Now store histo's to be used for limit setting
  string outname_limitsetting_TT = pathPNG+"/input_MVA";
  if(SignalSample == "hct") outname_limitsetting_TT += "HctTT_"+category+"_"+SignalSample+".root";
  else if(SignalSample == "hut") outname_limitsetting_TT += "HutTT_"+category+"_"+SignalSample+".root";

  TFile *outfile_limitsetting_TT = new TFile(outname_limitsetting_TT.c_str(),"recreate");
  outfile_limitsetting_TT->cd();

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
      if(name.find("TT") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_limitsetting_TT->Write("kOverwrite");

  //Now store histo's to be used for limit setting
  string outname_limitsetting_combSTandTT = "";
  if(khut == 0 && khct == 0)
  {
      outname_limitsetting_combSTandTT = pathPNG+"/input_MVA"+coupling_hut+"_"+coupling_hct;
      if(SignalSample == "hct") outname_limitsetting_combSTandTT += "HctComb_"+category+"_"+SignalSample+".root";
      else if(SignalSample == "hut") outname_limitsetting_combSTandTT += "HutComb_"+category+"_"+SignalSample+".root";
  }
  else outname_limitsetting_combSTandTT = pathPNG+"/input_MVA2DComb_"+coupling_hut+"_"+coupling_hct+category+".root";

  TFile *outfile_limitsetting_combSTandTT = new TFile(outname_limitsetting_combSTandTT.c_str(),"recreate");
  outfile_limitsetting_combSTandTT->cd();

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
      if(name.find("combSTandTT") == string::npos) continue;


     	TH1F *temp = it->second;

      temp->Write();
	}
	outfile_limitsetting_combSTandTT->Write("kOverwrite");

  cout << " ---------------------------------------- " << endl;
  cout << " Statistical unc. on MC " << endl;
  cout << " ... ttbb: " << sqrt(nEvents_forStatisticalUnc_ttbb)/nEvents_forStatisticalUnc_ttbb << " for eventcount: " << nEvents_forStatisticalUnc_ttbb << endl;
  cout << " ... ttcc: " << sqrt(nEvents_forStatisticalUnc_ttcc)/nEvents_forStatisticalUnc_ttcc << " for eventcount: " << nEvents_forStatisticalUnc_ttcc << endl;
  cout << " ... ttlf: " << sqrt(nEvents_forStatisticalUnc_ttlf)/nEvents_forStatisticalUnc_ttlf << " for eventcount: " << nEvents_forStatisticalUnc_ttlf << endl;
  cout << " ... other: " << sqrt(nEvents_forStatisticalUnc_other)/nEvents_forStatisticalUnc_other << " for eventcount: " << nEvents_forStatisticalUnc_other << endl;
  cout << " ... SThct: " << sqrt(nEvents_forStatisticalUnc_Signal_SThct)/nEvents_forStatisticalUnc_Signal_SThct << " for eventcount: " << nEvents_forStatisticalUnc_Signal_SThct << endl;
  cout << " ... TThct: " << sqrt(nEvents_forStatisticalUnc_Signal_TThct)/nEvents_forStatisticalUnc_Signal_TThct << " for eventcount: " << nEvents_forStatisticalUnc_Signal_TThct << endl;
  cout << " ... SThut: " << sqrt(nEvents_forStatisticalUnc_Signal_SThut)/nEvents_forStatisticalUnc_Signal_SThut << " for eventcount: " << nEvents_forStatisticalUnc_Signal_SThut << endl;
  cout << " ... TThut: " << sqrt(nEvents_forStatisticalUnc_Signal_TThut)/nEvents_forStatisticalUnc_Signal_TThut << " for eventcount: " << nEvents_forStatisticalUnc_Signal_TThut << endl;
  


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


void MakeTotalSystErrorBand_Distributions(string directory, string category, string coupling, string outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile)
{

    TFile *MSPlotFile = new TFile(outfilename.c_str(),"read");

    map<string,MultiSamplePlot*> MSPlot_ErrorBands;
    map<string,TH1F*> histo1D_nominal;
    map<string,TH1F*> histo1D_Up_SamplesAdded;
    map<string,TH1F*> histo1D_Down_SamplesAdded;

    map<string,TH1F*> histo1D_TotalUp;
    map<string,TH1F*> histo1D_TotalDown;

    //Define rate uncertainties
    Double_t LumiUncPlus = 0.026;
    Double_t LumiUncMinus = 0.026;
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
        
        if(NominalVariableNames[iVar].find("MVA_Max") != string::npos)//Extra TT systematics
        {
            string filenameExtra = directory+"/inputExtra_MVA";
            if(coupling == "hut") filenameExtra += "HutMAX_"+category+"_"+coupling+".root";
            else if(coupling == "hct") filenameExtra += "HctMAX_"+category+"_"+coupling+".root";
            TFile *ExtraTTSysts = new TFile( filenameExtra.c_str(),"read");//Here we read in extra systematics from another file


            TH1F *h_tmp_extrasysts_UEup = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_UEUp");
            h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_UEUp") );
            h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_UEUp") );
            histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"UEPlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_UEup->Clone();
            
            TH1F *h_tmp_extrasysts_scaleEnvelopeup = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_scaleEnvelopeUp");
            h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_scaleEnvelopeUp") );
            h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_scaleEnvelopeUp") );
            histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_scaleEnvelopeup->Clone();

            TH1F *h_tmp_extrasysts_PDFEnvelopeup = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_PDFEnvelopeUp");
            h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_PDFEnvelopeUp") );
            h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_PDFEnvelopeUp") );
            histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_PDFEnvelopeup->Clone();
            
            TH1F *h_tmp_extrasysts_UEdw = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_UEDown");
            h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_UEDown") );
            h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_UEDown") );
            histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"UEMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_UEdw->Clone();

            TH1F *h_tmp_extrasysts_scaleEnvelopedw = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_scaleEnvelopeDown");
            h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_scaleEnvelopeDown") );
            h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_scaleEnvelopeDown") );
            histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_scaleEnvelopedw->Clone();

            TH1F *h_tmp_extrasysts_PDFEnvelopedw = (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttbb_PDFEnvelopeDown");
            h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttcc_PDFEnvelopeDown") );
            h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("maxSTandTT_ttlf_PDFEnvelopeDown") );
            histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_PDFEnvelopedw->Clone();
         }
         else if(NominalVariableNames[iVar].find("MVA_Comb") != string::npos)//Extra TT systematics
         {
              string filenameExtra = directory+"/inputExtra_MVA";
              if(coupling == "hut") filenameExtra += "HutComb_"+category+"_"+coupling+".root";
              else if(coupling == "hct") filenameExtra += "HctComb_"+category+"_"+coupling+".root";
              TFile *ExtraTTSysts = new TFile( filenameExtra.c_str(),"read");//Here we read in extra systematics from another file

              TH1F *h_tmp_extrasysts_UEup = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_UEUp");
              h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_UEUp") );
              h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_UEUp") );
              histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"UEPlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_UEup->Clone();
              
              TH1F *h_tmp_extrasysts_scaleEnvelopeup = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_scaleEnvelopeUp");
              h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_scaleEnvelopeUp") );
              h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_scaleEnvelopeUp") );
              histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_scaleEnvelopeup->Clone();

              TH1F *h_tmp_extrasysts_PDFEnvelopeup = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_PDFEnvelopeUp");
              h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_PDFEnvelopeUp") );
              h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_PDFEnvelopeUp") );
              histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_PDFEnvelopeup->Clone();
              
              TH1F *h_tmp_extrasysts_UEdw = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_UEDown");
              h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_UEDown") );
              h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_UEDown") );
              histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"UEMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_UEdw->Clone();

              TH1F *h_tmp_extrasysts_scaleEnvelopedw = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_scaleEnvelopeDown");
              h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_scaleEnvelopeDown") );
              h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_scaleEnvelopeDown") );
              histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_scaleEnvelopedw->Clone();

              TH1F *h_tmp_extrasysts_PDFEnvelopedw = (TH1F*) ExtraTTSysts->Get("combSTandTT_ttbb_PDFEnvelopeDown");
              h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttcc_PDFEnvelopeDown") );
              h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("combSTandTT_ttlf_PDFEnvelopeDown") );
              histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_PDFEnvelopedw->Clone();
         }
         else if(NominalVariableNames[iVar].find("MVA_ST") != string::npos)//Extra TT systematics
         {
                string filenameExtra = directory+"/inputExtra_MVA";
                if(coupling == "hut") filenameExtra += "HutST_"+category+"_"+coupling+".root";
                else if(coupling == "hct") filenameExtra += "HctST_"+category+"_"+coupling+".root";
                TFile *ExtraTTSysts = new TFile( filenameExtra.c_str(),"read");//Here we read in extra systematics from another file

                TH1F *h_tmp_extrasysts_UEup = (TH1F*) ExtraTTSysts->Get("ST_ttbb_UEUp");
                h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_UEUp") );
                h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_UEUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"UEPlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_UEup->Clone();
                
                TH1F *h_tmp_extrasysts_scaleEnvelopeup = (TH1F*) ExtraTTSysts->Get("ST_ttbb_scaleEnvelopeUp");
                h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_scaleEnvelopeUp") );
                h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_scaleEnvelopeUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_scaleEnvelopeup->Clone();

                TH1F *h_tmp_extrasysts_PDFEnvelopeup = (TH1F*) ExtraTTSysts->Get("ST_ttbb_PDFEnvelopeUp");
                h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_PDFEnvelopeUp") );
                h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_PDFEnvelopeUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_PDFEnvelopeup->Clone();
                
                TH1F *h_tmp_extrasysts_UEdw = (TH1F*) ExtraTTSysts->Get("ST_ttbb_UEDown");
                h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_UEDown") );
                h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_UEDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"UEMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_UEdw->Clone();

                TH1F *h_tmp_extrasysts_scaleEnvelopedw = (TH1F*) ExtraTTSysts->Get("ST_ttbb_scaleEnvelopeDown");
                h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_scaleEnvelopeDown") );
                h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_scaleEnvelopeDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_scaleEnvelopedw->Clone();

                TH1F *h_tmp_extrasysts_PDFEnvelopedw = (TH1F*) ExtraTTSysts->Get("ST_ttbb_PDFEnvelopeDown");
                h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttcc_PDFEnvelopeDown") );
                h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("ST_ttlf_PDFEnvelopeDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_PDFEnvelopedw->Clone();
         }
         else if(NominalVariableNames[iVar].find("MVA_TT") != string::npos)//Extra TT systematics
         {
                string filenameExtra = directory+"/inputExtra_MVA";
                if(coupling == "hut") filenameExtra += "HutTT_"+category+"_"+coupling+".root";
                else if(coupling == "hct") filenameExtra += "HctTT_"+category+"_"+coupling+".root";
                TFile *ExtraTTSysts = new TFile( filenameExtra.c_str(),"read");//Here we read in extra systematics from another file


                TH1F *h_tmp_extrasysts_UEup = (TH1F*) ExtraTTSysts->Get("TT_ttbb_UEUp");
                h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_UEUp") );
                h_tmp_extrasysts_UEup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_UEUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"UEPlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_UEup->Clone();
                
                TH1F *h_tmp_extrasysts_scaleEnvelopeup = (TH1F*) ExtraTTSysts->Get("TT_ttbb_scaleEnvelopeUp");
                h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_scaleEnvelopeUp") );
                h_tmp_extrasysts_scaleEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_scaleEnvelopeUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_scaleEnvelopeup->Clone();

                TH1F *h_tmp_extrasysts_PDFEnvelopeup = (TH1F*) ExtraTTSysts->Get("TT_ttbb_PDFEnvelopeUp");
                h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_PDFEnvelopeUp") );
                h_tmp_extrasysts_PDFEnvelopeup->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_PDFEnvelopeUp") );
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopePlusPlus").c_str()] = (TH1F*) h_tmp_extrasysts_PDFEnvelopeup->Clone();
                
                TH1F *h_tmp_extrasysts_UEdw = (TH1F*) ExtraTTSysts->Get("TT_ttbb_UEDown");
                h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_UEDown") );
                h_tmp_extrasysts_UEdw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_UEDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"UEMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_UEdw->Clone();

                TH1F *h_tmp_extrasysts_scaleEnvelopedw = (TH1F*) ExtraTTSysts->Get("TT_ttbb_scaleEnvelopeDown");
                h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_scaleEnvelopeDown") );
                h_tmp_extrasysts_scaleEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_scaleEnvelopeDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_scaleEnvelopedw->Clone();

                TH1F *h_tmp_extrasysts_PDFEnvelopedw = (TH1F*) ExtraTTSysts->Get("TT_ttbb_PDFEnvelopeDown");
                h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttcc_PDFEnvelopeDown") );
                h_tmp_extrasysts_PDFEnvelopedw->Add( (TH1F*) ExtraTTSysts->Get("TT_ttlf_PDFEnvelopeDown") );
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopeMinusMinus").c_str()] = (TH1F*)h_tmp_extrasysts_PDFEnvelopedw->Clone();
         }


        //Add the nominal samples into 1 histogram
        for(int iDataName = 0; iDataName < datasetNames.size(); iDataName++)
        {
            TDirectory *subdir_nominal = (TDirectory*) MSPlotFile->Get(("MultiSamplePlot_"+NominalVariableNames[iVar]).c_str());
            subdir_nominal->cd();
            string nominalname = (NominalVariableNames[iVar]+"_"+datasetNames[iDataName]+"_");
            
            TH1F *h_tmp =  (TH1F*)subdir_nominal->Get(nominalname.c_str());

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
                
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"UEPlusPlus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopePlusPlus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopePlusPlus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"UEMinusMinus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"scaleEnvelopeMinusMinus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"PDFEnvelopeMinusMinus").c_str()]->Add(h_tmp); //Add the non-ttbar nominal backgrounds to the regular systematics

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
            cout << "  - MakeTotalSystErrorBand_Distributions(): Making systematic band for " << varNameSys << endl;
            TDirectory *subdir_sys = (TDirectory*) MSPlotFile->Get(("MultiSamplePlot_"+varNameSys).c_str());
            subdir_sys->cd();

            histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()] = 0;
            histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()] = 0;

            for(int iDataName = 0; iDataName < datasetNames.size(); iDataName++)
            {

                TH1F *h_tmp =  (TH1F*)subdir_sys->Get((varNameSys+"_"+datasetNames[iDataName]+"_").c_str());

                
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
        systematics.push_back("UEPlus");
        systematics.push_back("UEMinus");
        systematics.push_back("scaleEnvelopePlus");
        systematics.push_back("scaleEnvelopeMinus");
        systematics.push_back("PDFEnvelopePlus");
        systematics.push_back("PDFEnvelopeMinus");
       
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
            
            histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()]->SetBinContent(iBin,bincontent_nominal + sqrt(bincontent_up_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
            histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()]->SetBinContent(iBin,bincontent_nominal  - sqrt(bincontent_down_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
        }
        
        //Delete the last entries for the systematics that were added manually in MakeTotalSystErrorBand_Distributions
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
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
        if(samplename.find("ST_tHToBB_1L_Kappa_hct")!= string::npos)
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

double OptimalCut_CombTraining(string category, string coupling)
{
    double MVA_cutvalue = -1.;
/* // Cut values for S/sqrt(S+b)
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
*/
//    MVA_cutvalue = 0.3;

    //Cut values for ROC-optimal
    if(coupling == "hut")
    {
        if(category == "b2j3") MVA_cutvalue = 0.0385566;
        else if(category == "b2j4") MVA_cutvalue = 0.0134737;
        else if(category == "b3j3") MVA_cutvalue = 0.036568;
        else if(category == "b3j4") MVA_cutvalue = 0.0326008;
        else if(category == "b4j4") MVA_cutvalue = -0.0120028;
    }
    else if(coupling == "hct")
    {
        if(category == "b2j3") MVA_cutvalue = 0.0169624;
        else if(category == "b2j4") MVA_cutvalue = 0.0502082;
        else if(category == "b3j3") MVA_cutvalue = 0.0312207;
        else if(category == "b3j4") MVA_cutvalue = 0.0371161;
        else if(category == "b4j4") MVA_cutvalue = 0.0295859;
    }
    
    return MVA_cutvalue;
}

double PostFitScaleFactor(string category, string coupling, string samplename)
{
    double SF = 1.;
    
    if(coupling == "hct")
    {
         if(category == "b2j3")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.98691;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.66657;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.23614;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.25245;
            else SF = 1.40752;
         }
         else if(category == "b2j4")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.93653;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.63244;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.14704;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.24838;
            else SF = 1.42345;
         }
         else if(category == "b3j3")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.94788;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.63917;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.22019;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.25842;
            else SF = 1.43382;
         }
         else if(category == "b3j4")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.88371;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.59189;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.16575;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.24746;
            else SF = 1.46548;
         }
         else if(category == "b4j4")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.77713;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.53112;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.11378;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.24292;
            else SF = 1.57951;
         }
    }
    else if(coupling == "hut")
    {
         if(category == "b2j3")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.02092;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.69053;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.18725;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.09828;
            else SF = 1.23193;
         }
         else if(category == "b2j4")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.96791;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.64368;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.07841;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.09192;
            else SF = 1.19029;
         }
         else if(category == "b3j3")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.98820;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.64636;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.23138;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.10573;
            else SF = 1.30119;
         }
         else if(category == "b3j4")
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.92554;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.60601;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.16177;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.10060;
            else SF = 1.24162;
         }
/*         else if(category == "b4j4)
         {
            if(samplename.find("TTJets_ll") != string::npos) SF = 0.77713;
            else if(samplename.find("TTJets_cc") != string::npos) SF = 0.53112;
            else if(samplename.find("TTJets_bb") != string::npos) SF = 1.11378;
            else if(samplename.find("NP_overlay") != string::npos) SF = 0.24292;
            else SF = 1.57951;
         }
*/    }

    return SF;
}
