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
        cout << "    string channel            = argv[3];" << endl;
        cout << "    string date            = argv[4];" << endl;
        cout << "    bool PVreweighing = strtol(argv[5], NULL,10);" << endl;
        cout << "    bool doJESSys  = strtol(argv[6], NULL,10);" << endl;
        cout << "    bool doJERSys  = strtol(argv[7], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[6], NULL,10);" << endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    string channel            = argv[3];
    string date            = argv[4];
    bool PVreweighing = strtol(argv[5], NULL,10);
    bool doJESSys  = strtol(argv[6], NULL,10);
    bool doJERSys  = strtol(argv[7], NULL,10);
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
/*    WhatSysts.push_back("noSF");
    WhatSysts.push_back("OnlyTopPtSF");
    WhatSysts.push_back("OnlyBTagSF");
    WhatSysts.push_back("OnlyPUSF");
//    WhatSysts.push_back("OnlyLepSF");
//    WhatSysts.push_back("OnlyNLOSF");
    WhatSysts.push_back("NoTopPtSF");
    WhatSysts.push_back("NoBTagSF");
    WhatSysts.push_back("NoPUSF");
//    WhatSysts.push_back("NoLepSF");
//    WhatSysts.push_back("NoNLOSF");
    if(doJESSys) WhatSysts.push_back("JESPlus");
    if(doJESSys) WhatSysts.push_back("JESMinus");
    if(doJERSys) WhatSysts.push_back("JERPlus");
    if(doJERSys) WhatSysts.push_back("JERMinus");
*/    WhatSysts.push_back("");   

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
    WhatSysts_noJECs.push_back("NoTopPtSF");
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
            
            ttbar_ll->SetColor(kBlue);
            ttbar_cc->SetColor(kBlue-3);
            ttbar_bb->SetColor(kBlue+2);
            
            datasets_splittedTTbar.pop_back();
            datasets_splittedTTbar.push_back(ttbar_ll);
            datasets_splittedTTbar.push_back(ttbar_cc);
            datasets_splittedTTbar.push_back(ttbar_bb);

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

    for(int iSyst = 0; iSyst<WhatSysts.size();iSyst++)
    {
    
    
        MSPlot[("NPV"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("NPV"+WhatSysts[iSyst]).c_str(), 51, -0.5, 50.5, "Number of PV","Events", category); 
        MSPlot[("NCSVv2Ljets"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("NCSVv2Ljets"+WhatSysts[iSyst]).c_str(), 11, -0.5, 10.5, "Number of CSVv2 L jets","Events", category); 
        MSPlot[("NCSVv2Mjets"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("NCSVv2Mjets"+WhatSysts[iSyst]).c_str(), 11, -0.5, 10.5, "Number of CSVv2 M jets","Events", category); 
        MSPlot[("NCSVv2Tjets"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("NCSVv2Tjets"+WhatSysts[iSyst]).c_str(), 11, -0.5, 10.5, "Number of CSVv2 T jets","Events", category); 
        MSPlot[("Njets"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("Njets"+WhatSysts[iSyst]).c_str(), 11, -0.5, 10.5, "Number jets","Events", category); 
        MSPlot[("LeptonPt"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("LeptonPt"+WhatSysts[iSyst]).c_str(), 50, 20., 300., "Lepton Pt","Events", category, "GeV"); 
        MSPlot[("LeptonEta"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("LeptonEta"+WhatSysts[iSyst]).c_str(), 50, -2.5, 2.5, "Lepton eta","Events", category); 
        MSPlot[("LeptonPhi"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("LeptonPhi"+WhatSysts[iSyst]).c_str(), 50, -3.2, 3.2, "Lepton phi","Events", category); 
        MSPlot[("LeptonCharge"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("LeptonCharge"+WhatSysts[iSyst]).c_str(), 3, -1.5, 1.5, "Lepton Charge","Events", category); 
/*        MSPlot[("JetPt"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("JetPt"+WhatSysts[iSyst]).c_str(), 50, 20., 300., "Jet Pt","Events", category, "GeV"); 
        MSPlot[("JetEta"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("JetEta"+WhatSysts[iSyst]).c_str(), 50, -2.5, 2.5, "Jet eta","Events", category); 
        MSPlot[("JetPhi"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("JetPhi"+WhatSysts[iSyst]).c_str(), 50, -3.2, 3.2, "Jet phi","Events", category);
        MSPlot[("JetCSVv2"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("JetCSVv2"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "Jet CSVv2","Events", category);
        MSPlot[("JetcMVAv2"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("JetcMVAv2"+WhatSysts[iSyst]).c_str(), 50, -1., 1., "Jet cMVAv2","Events", category); 

        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
        MSPlot[("HiggsMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
        MSPlot[("HiggsEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (Higgs)","Events", category);
        MSPlot[("HiggsEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (Higgs)","Events", category);
        MSPlot[("TopLepMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepMass_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
        MSPlot[("TopLepMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepMass_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
        MSPlot[("TopLepPt_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepPt_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "Pt(LepTop)","Events", category,"GeV");
        MSPlot[("TopLepPt_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepPt_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "Pt(LepTop)","Events", category,"GeV");
        MSPlot[("TopLepEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepEta_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (LepTop)","Events", category);
        MSPlot[("TopLepEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepEta_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, -5.0, 5.0, "eta (LepTop)","Events", category);
        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "#Delta R(H,LepTop)","Events", category);
        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "#Delta R(H,LepTop)","Events", category);
        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopHadMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopHadMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(HadTop)","Events", category,"GeV");
        MSPlot[("TopLepMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepMass_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
        MSPlot[("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(HadTop,LepTop)","Events", category);
        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepBJetCSVv2_TOPTOPLEPHAF"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("HiggsMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 50., 250, "M(Higgs)","Events", category,"GeV");
        MSPlot[("TopLepMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepMass_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 80., 300., "M(LepTop)","Events", category,"GeV");
        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(Hb1,Hb2)","Events", category);
        MSPlot[("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 5.0, "DR(H,LepTop)","Events", category);
        MSPlot[("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, 0., 1., "CSVv2 disc.","Events", category);
        MSPlot[("MVA_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TOPTOPLEPHBB"+WhatSysts[iSyst]).c_str(), 50, -1., 1., "BDT disc.","Events", category);
        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TOPHLEPBB_hut"+WhatSysts[iSyst]).c_str(), 50, -1., 1., "BDT disc.","Events", category);
        MSPlot[("MVA_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TOPHLEPBB_hct"+WhatSysts[iSyst]).c_str(), 50, -1., 1., "BDT disc.","Events", category);
        MSPlot[("MVA_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str() ] = new MultiSamplePlot(datasets_splittedTTbar, ("MVA_TOPTOPLEPHAD"+WhatSysts[iSyst]).c_str(), 50, -1., 1., "BDT disc.","Events", category);
*/    }
  
 


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
            if(!isData) postfix = WhatSysts[JecCounter];
	      

		        cout<<"Dataset:  :"<<dataSetName<<endl;
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
            Double_t W_MuonTrigSF_Runs273158to274093;//Used in calculation for W_MuonTrigSF
            Double_t W_MuonTrigSF_Runs274094to276097;//Used in calculation for W_MuonTrigSF
            Double_t W_ElectronIDSF; //One of the 2 components for the total electron SF
            Double_t W_ElectronRecoSF; //One of the 2 components for the total electron SF
            Double_t W_TopPtReweighing;
          
            Int_t run_num;
            Int_t evt_num;
            Int_t lumi_num;
            Int_t nvtx;
            Int_t npu;
            Int_t genTTX;

	          // variables for electrons
            Double_t eta_superCluster_electron;
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
            Int_t LepCharge;
      
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

            Double_t MC_TopPt;
            Double_t MC_AntiTopPt;
            
            // Weights
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Plus",&W_fleptonSF_Plus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Minus",&W_fleptonSF_Minus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Minus",&W_puSF_Minus);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Plus",&W_puSF_Plus);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_central",&W_btagWeight_CSVv2M_mujets_central); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_up",&W_btagWeight_CSVv2M_mujets_up);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_CSVv2M_mujets_down",&W_btagWeight_CSVv2M_mujets_down); 
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
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight1",&W_weight1);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight2",&W_weight2);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight3",&W_weight3); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight4",&W_weight4);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight5",&W_weight5); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight6",&W_weight6);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight7",&W_weight7); 
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_weight8",&W_weight8);  
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_TopPtReweighing",&W_TopPtReweighing);  

            ttree[(dataSetName).c_str()]->SetBranchAddress("I_run_num",&run_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_evt_num",&evt_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_lumi_num",&lumi_num);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_npu",&npu);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_genTTX",&genTTX);


            // electrons
            ttree[(dataSetName).c_str()]->SetBranchAddress("eta_superCluster_electron",&eta_superCluster_electron);
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
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_LepCharge",&LepCharge);
            
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
		                if(!doInclusive)
		                {
		                    if(nJets_CSVM != baseline_bjets)  continue;

		                    if(baseline_jets == 3 && nJets != baseline_jets) continue;
		                    else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		                }
                    if( W_nloWeight > 0) nPos++;
                    else if( W_nloWeight < 0) nNeg ++;
                }
                nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
            }		

            Double_t average_TopPtWeight = 0.;
            Double_t average_TopPtWeight_Up = 0.;
            if(dataSetName.find("TTJets") != string::npos)
            {
                int nEventsPassed = 0;
                for (int k = 0; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);
		                if(!doInclusive)
		                {
		                    if(nJets_CSVM != baseline_bjets)  continue;

		                    if(baseline_jets == 3 && nJets != baseline_jets) continue;
		                    else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		                }
		                double TopPtReweighing_Up = 1+ 2*(1-W_TopPtReweighing);

                    average_TopPtWeight = average_TopPtWeight + W_TopPtReweighing;
                    average_TopPtWeight_Up = average_TopPtWeight_Up + TopPtReweighing_Up;
                    nEventsPassed++;
                }
                average_TopPtWeight = average_TopPtWeight/nEventsPassed;
                average_TopPtWeight_Up = average_TopPtWeight_Up/nEventsPassed;
            }
		
      	    //***********************************************RUNNING OVER EVENTS**********************************************
		        for (int j = 0; j<nEntries; j++)
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
                


                    //Nominal scale factor -- scale factors for systematic shifts are calculated below
                    ScaleFactor *= W_puSF_applied;
                    ScaleFactor *= W_fleptonSF;
                    ScaleFactor *= W_btagWeight_shape;
                    ScaleFactor *= nloSF;
                    if(dataSetName.find("TTJets") != string::npos) ScaleFactor *= W_TopPtReweighing/average_TopPtWeight;
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
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos)
                                {
                                    double TopPtReweighing_Up = 1+ 2*(1-W_TopPtReweighing);
                                    SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= TopPtReweighing_Up/average_TopPtWeight_Up;
                                }
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtMinus")//Apply no TopPt reweighing
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
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
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoTopPtSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoBTagSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoPUSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoLepSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoNLOSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                        }//if(!isData)
                        else SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] = 1.;


                        //-----------------------------------------------------------------------------------------------------------
                        // Fill Plots
                        //-----------------------------------------------------------------------------------------------------------
				                MSPlot[("NCSVv2Ljets"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(nJets_CSVL, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
				                MSPlot[("NCSVv2Mjets"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(nJets_CSVM, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
				                MSPlot[("NCSVv2Tjets"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(nJets_CSVT, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
				                MSPlot[("Njets"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(nJets, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("LeptonPt"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(pt_lepton, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("LeptonEta"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(eta_lepton, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("LeptonPhi"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(phi_lepton, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("LeptonCharge"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(LepCharge, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("NPV"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(nvtx, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
/*                        for(int i_Jet = 0; i_Jet < nJets; i_Jet++)
                        {
                            MSPlot[("JetPt"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(pt_jet[i_Jet], Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                            MSPlot[("JetEta"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(eta_jet[i_Jet], Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                            MSPlot[("JetPhi"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(phi_jet[i_Jet], Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                            MSPlot[("JetCSVv2"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(CSVv2[i_Jet], Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                            MSPlot[("JetcMVAv2"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(cMVA[i_Jet], Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        }			                


                        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsMass_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsEta_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsEta_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsEta_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsEta_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepMass_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepMass_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepMass_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepPt_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepPt_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepPt_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepPt_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepEta_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepEta_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepEta_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepEta_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepHiggsDr_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepHiggsDr_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1CSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepBJetCSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepBJetCSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopHadMass_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopHadMass_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepMass_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepMass_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepTopHadDr_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepBJetCSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopHadBJetCSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopHadWNonBJet1CSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopHadWNonBJet2CSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsMass_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsMass_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepMass_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepMass_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepHiggsDr_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet1CSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(HiggsBJet2CSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopLepBJetCSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(TopHadNonBJetCSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("MVA_TOPTOPLEPHBB"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVA_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVA_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("MVA_TOPHLEPBB_hct"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVA_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        MSPlot[("MVA_TOPTOPLEPHAD"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(MVA_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
*/                    }
               }
                if(filepath.find("JESMinus") != string::npos || filepath.find("JESPlus") != string::npos  || filepath.find("JERMinus") != string::npos || filepath.find("JERPlus") != string::npos || isData || WhatSysts[JecCounter] == "")
               {
				                MSPlot[("NCSVv2Ljets"+WhatSysts[JecCounter]).c_str()]->Fill(nJets_CSVL, Sample, ScalePlots, Luminosity * ScaleFactor);
				                MSPlot[("NCSVv2Mjets"+WhatSysts[JecCounter]).c_str()]->Fill(nJets_CSVM, Sample, ScalePlots, Luminosity * ScaleFactor);
				                MSPlot[("NCSVv2Tjets"+WhatSysts[JecCounter]).c_str()]->Fill(nJets_CSVT, Sample, ScalePlots, Luminosity * ScaleFactor);
				                MSPlot[("Njets"+WhatSysts[JecCounter]).c_str()]->Fill(nJets, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("LeptonPt"+WhatSysts[JecCounter]).c_str()]->Fill(pt_lepton, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("LeptonEta"+WhatSysts[JecCounter]).c_str()]->Fill(eta_lepton, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("LeptonPhi"+WhatSysts[JecCounter]).c_str()]->Fill(phi_lepton, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("LeptonCharge"+WhatSysts[JecCounter]).c_str()]->Fill(LepCharge, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("NPV"+WhatSysts[JecCounter]).c_str()]->Fill(nvtx, Sample, ScalePlots, Luminosity * ScaleFactor);
/*                        for(int i_Jet = 0; i_Jet < nJets; i_Jet++)
                        {
                            MSPlot[("JetPt"+WhatSysts[JecCounter]).c_str()]->Fill(pt_jet[i_Jet], Sample, ScalePlots, Luminosity * ScaleFactor);
                            MSPlot[("JetEta"+WhatSysts[JecCounter]).c_str()]->Fill(eta_jet[i_Jet], Sample, ScalePlots, Luminosity * ScaleFactor);
                            MSPlot[("JetPhi"+WhatSysts[JecCounter]).c_str()]->Fill(phi_jet[i_Jet], Sample, ScalePlots, Luminosity * ScaleFactor);
                            MSPlot[("JetCSVv2"+WhatSysts[JecCounter]).c_str()]->Fill(CSVv2[i_Jet], Sample, ScalePlots, Luminosity * ScaleFactor);
                            MSPlot[("JetcMVAv2"+WhatSysts[JecCounter]).c_str()]->Fill(cMVA[i_Jet], Sample, ScalePlots, Luminosity * ScaleFactor);
                        }			                


                        MSPlot[("HiggsMass_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsMass_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsMass_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsEta_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsEta_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsEta_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsEta_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepMass_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepMass_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepMass_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepMass_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepPt_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepPt_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepPt_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepPt_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepEta_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepEta_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepEta_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepEta_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepHiggsDr_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepHiggsDr_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepHiggsDr_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1CSVv2_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1CSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet2CSVv2_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet2CSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepBJetCSVv2_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepBJetCSVv2_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepBJetCSVv2_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopHadMass_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopHadMass_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepMass_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepMass_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepTopHadDr_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepTopHadDr_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepBJetCSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopHadBJetCSVv2_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopHadBJetCSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopHadWNonBJet1CSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(TopHadWNonBJet2CSVv2_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsMass_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsMass_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepMass_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepMass_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepHiggsDr_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepHiggsDr_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet1CSVv2_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet1CSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("HiggsBJet2CSVv2_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(HiggsBJet2CSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopLepBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(TopLepBJetCSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("TopHadNonBJetCSVv2_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(TopHadNonBJetCSVv2_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("MVA_TOPTOPLEPHBB"+WhatSysts[JecCounter]).c_str()]->Fill(MVA_TOPTOPLEPHBB, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("MVA_TOPHLEPBB_hut"+WhatSysts[JecCounter]).c_str()]->Fill(MVA_TOPHLEPBB_hut, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("MVA_TOPHLEPBB_hct"+WhatSysts[JecCounter]).c_str()]->Fill(MVA_TOPHLEPBB_hct, Sample, ScalePlots, Luminosity * ScaleFactor);
                        MSPlot[("MVA_TOPTOPLEPHAD"+WhatSysts[JecCounter]).c_str()]->Fill(MVA_TOPTOPLEPHAD, Sample, ScalePlots, Luminosity * ScaleFactor);
*/
               }
			                
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
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  string outfilename = pathPNG+"/Output.root";

  TFile *outfile = new TFile(outfilename.c_str(),"recreate");
//  outfile->cd();

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
      temp->Draw("MyMSP_"+name, 1, false, false, false, 1);
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
      temp->Write(outfile, name, false,pathPNG, "png");
	}

  outfile->Write("kOverwrite");
  outfile->Close();
  
  cout << "  - Making total systematic bands " << endl;
  string errorbandfile = (pathPNG+"/Systematics_BareHistos.root");
  MakeTotalSystErrorBand_Distributions(outfilename, WhatSysts, datasetnames_backgrounds, NominalVariableNames, errorbandfile);



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
      temp->Draw("MyMSP_"+name, 1, true, true, true, 1);
      bool writePng = false;
      temp->Write(outfile_errorbands, name, true,pathPNG, "png");
	}
	outfile_errorbands->Write("kOverwrite");



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


void MakeTotalSystErrorBand_Distributions(string outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile)
{

    TFile *MSPlotFile = new TFile(outfilename.c_str(),"read");

    map<string,MultiSamplePlot*> MSPlot_ErrorBands;
    map<string,TH1F*> histo1D_nominal;
    map<string,TH1F*> histo1D_Up_SamplesAdded;
    map<string,TH1F*> histo1D_Down_SamplesAdded;

    map<string,TH1F*> histo1D_TotalUp;
    map<string,TH1F*> histo1D_TotalDown;

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
            if(iDataName == 0) histo1D_nominal[NominalVariableNames[iVar].c_str()] = (TH1F*) h_tmp->Clone(/*NominalVariableNames[iVar].c_str()*/);
            else histo1D_nominal[NominalVariableNames[iVar].c_str()]->Add(h_tmp);
        }
            
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
                
                if(systematics[iSys].find("Plus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()] = (TH1F*) h_tmp->Clone(/*(varNameSys+"Plus").c_str()*/);
                    else histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()]->Add(h_tmp);
                }
                else if(systematics[iSys].find("Minus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()] = (TH1F*) h_tmp->Clone(/*(varNameSys+"Minus").c_str()*/);
                    else histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()]->Add(h_tmp);
//cout << "Integral of Minus systematic " << histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()]->GetBinContent(0) << endl;
                }
            }
        }
            
       
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
            

            histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()]->SetBinContent(iBin,bincontent_nominal + sqrt(bincontent_up_squared));
            histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()]->SetBinContent(iBin,bincontent_nominal  - sqrt(bincontent_down_squared));
        }
        
        
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
