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




using namespace std;
using namespace TopTree;



/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;


map<string,TFile*> globalFileObj;
map<string,TTree*> globalttree;
map<string,TFile*> globalFileObj0;
map<string,TTree*> globalttree0;

vector<string> CutsselecTable;
vector<string> CutsselecTabledetail;
vector<string> CutsselecTabledetail3;

// functions prototype
string intToStr (int number);

inline bool FileExists (const string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
  
}


int main(int argc, char *argv[])
{
  
  if(argc < 3)
  {
    cerr << "INVALID number of arguments. The necessary arguments are: " << endl;
    cout << "    string channel            = argv[1];" << endl;
    cout << "    string date            = argv[2];" << endl;
    cout << "    bool debug         =strtol(argv[3], NULL,10);" << endl;
    
    return 1;
  }
  
  
  string decaychannel            = argv[1];
  string date            = argv[2];
  bool debug         =strtol(argv[3], NULL,10);
  
  int dChan = 9;
  if(decaychannel.find("eee")!=string::npos) dChan = 3;
  if(decaychannel.find("eeu")!=string::npos) dChan = 2;
  if(decaychannel.find("uue")!=string::npos) dChan = 1;
  if(decaychannel.find("uuu")!=string::npos) dChan = 0;
  if(decaychannel.find("all")!=string::npos) dChan = 9;
  
  cout << "------------------------------------------------------------------------------------------------" << endl;
  cout << "Begin program" << endl;
  cout << "------------------------------------------------------------------------------------------------" << endl;
  
  
  clock_t start = clock();
  
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////
  // **************** Preparing samples ********************//
  ///////////////////////////////////////////////////////////////////////////////////
  cout << " ... Making the TreeProcessor .xml files " << endl;
  //system("python MakeXMLforTreeProcessor.py");
  
  string xmlNom;
  //  xmlNom = "config/Run2TriLepton_samples_TreeProcessor.xml";
  xmlNom = "config/Run2TriLepton_samples_analy.xml";
  TString TreePath = "NtupleMakerOutput/MergedTuples/all/"+ date;
  cout << "tree path " << TreePath << endl;
  if(!FileExists(string(TreePath+"/data.root")))
  {
    system(("hadd -f "+TreePath+"/data.root "+TreePath+"/data_*.root").Data());
  }
  
  
  const char *xmlfile = xmlNom.c_str();
  cout << "used config file: " << xmlfile << endl;
  
  //***************************************************LOADING DATASETS****************************************************
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  string dataSetName, filepath;
  string dataSetName0, filepath0;
  float Luminosity = 0;
  
  
  MultiSamplePlot* cutFlowPlot = new MultiSamplePlot(datasets, "Cut Flow", 12, -0.5, 11.5, "Applied Cut");
  MultiSamplePlot* detailcutFlowPlot = new MultiSamplePlot(datasets, "Cut Flow Detail", 10, -0.5, 9.5, "Applied Cut");
  MultiSamplePlot* detailcutFlowPlot3 = new MultiSamplePlot(datasets, "Cut Flow Last 3 cuts", 10, -0.5, 9.5, "Applied Cut");
  
  CutsselecTable.clear();
  CutsselecTabledetail.clear();
  
  
  
  CutsselecTabledetail3.clear();
  
  
  
  cout << "defining cuts" << endl;
  vector <string> cutstep_string;
  cutstep_string.clear();
  cutstep_string.push_back("trigger");
  cutstep_string.push_back("3lep");
  cutstep_string.push_back("VetoMu");
  cutstep_string.push_back("VetoEl");
  cutstep_string.push_back("OSSF");
  cutstep_string.push_back("Zmass");
  cutstep_string.push_back(">1jet");
  cutstep_string.push_back(">0CSVL");
  cutstep_string.push_back("mWt");
  cutstep_string.push_back("SMtop");
  cutstep_string.push_back("METfilter");
  
  for(int iC = 0 ; iC < cutstep_string.size(); iC++){
    if(debug) cout << "cut " << iC  << " label " << cutstep_string[iC] <<  endl;
    
    CutsselecTable.push_back(cutstep_string[iC]);
    
  }
  for(int iC = 1 ; iC < cutstep_string.size(); iC++){
    if(debug) cout << "cut " << iC  << " label " << cutstep_string[iC] <<  endl;
    CutsselecTabledetail.push_back(cutstep_string[iC]);
  }
  for(int iC = 6 ; iC < cutstep_string.size(); iC++){
    if(debug) cout << "cut " << iC  << " label " << cutstep_string[iC] <<  endl;
    CutsselecTabledetail3.push_back(cutstep_string[iC]);
  }
  cout << " end defining cuts" << endl;
  
  //***************************************************GETTING LUMI FROM DATA IN XML****************************************************
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    dataSetName = datasets[d]->Name();
    
    
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
    {
      Luminosity = datasets[d]->EquivalentLumi();
      cout << "lumi " << Luminosity << endl;
    }
    
    
  }
  if(Luminosity == 0)
  {
    cout << "Luminosity is 0. Please check the data-luminosity in your xml file. Exiting program..." << endl;
    return 1;
  }
  
  SelectionTable selecTable(CutsselecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  selecTable.SetPrecision(2);
  
  //***********************************************RUNNING OVER DATASETS**********************************************
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    
    dataSetName = datasets[d]->Name();
    
    cout<<"Dataset:  :"<<dataSetName<< " with title " << datasets[d]->Title() <<  endl;
 //   if(dataSetName.find("tZq")==string::npos) continue;
    filepath = TreePath+"/"+dataSetName + ".root";
  //  filepath = "NtupleMakerOutput/Ntuples/Ntuples_170117/FCNC_3L_tZq_80X_1.root";
    
    if (debug)
    {
      cout<<"filepath: "<<filepath<<endl;
      cout <<"Equivalent luminosity of the dataset is: " << datasets[d]->EquivalentLumi() << endl;
    }
    
    
    bool isData= false;
    bool isAMC = false;
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
    {
      if(debug) cout << "Data found" << endl;
      isData =true;
    }
    else if(dataSetName.find("NLO") != string::npos || dataSetName.find("nlo") !=string::npos || dataSetName.find("amc") !=string::npos) isAMC = true;
    
    std::vector<double> cuts;
    
    
    //***********************************************IMPORTING VARIABLES**********************************************
    
    string sglobaltree = "globaltree";
    globalFileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset
    string globalTTreename = sglobaltree;
    globalttree[dataSetName.c_str()] = (TTree*)globalFileObj[dataSetName.c_str()]->Get(globalTTreename.c_str()); //get ttre for each dataset
    if(debug) cout << "globalttree " << globalttree[dataSetName.c_str()]<< endl;
    int globalnEntries = (int)globalttree[dataSetName.c_str()]->GetEntries();
    cout<<"                 nEntries gt: "<<globalnEntries<<endl;
    
    //globalttree[dataSetName.c_str()]->SetMakeClass(1);
    
    
    
    Int_t cutst[15]; //
    Int_t cutst_eee[15];
    Int_t cutst_eeu[15];
    Int_t cutst_uuu[15];
    Int_t cutst_uue[15];
    
    Int_t nCuts_ = 10; //REDEFINE if ncuts change
    
    Int_t nofPosWeights = 0;
    Int_t nofNegWeights = 0;
    
    Int_t nPosW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofPosWeights",&nPosW);
    
    Int_t nNegW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nofNegWeights",&nNegW);
    
    Int_t nEvents;
    globalttree[dataSetName.c_str()]->SetBranchAddress("nEv",&nEvents);
    
    Int_t SumW;
    globalttree[dataSetName.c_str()]->SetBranchAddress("sumW",&SumW);
    
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nCuts",&nCuts_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep",cutst);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_eee",cutst_eee);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_eeu",cutst_eeu);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_uuu",cutst_uuu);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_uue",cutst_uue);
    
    
    
    int nTrigg_ = 0;
    int n3lep_ = 0;
    int nVetoMu_ = 0;
    int nVetoEl_ = 0;
    int nOS_ = 0;
    int nZmass_ = 0;
    int nJet_ = 0;
    int nBJet_ = 0;
    int nMWT_ = 0;
    int nSMtop_ = 0;
    int nMET_ = 0;
    
    int nTrigg_eee_ = 0;
    int n3lep_eee_ = 0;
    int nVetoMu_eee_ = 0;
    int nVetoEl_eee_ = 0;
    int nOS_eee_ = 0;
    int nZmass_eee_ = 0;
    int nJet_eee_ = 0;
    int nBJet_eee_ = 0;
    int nMWT_eee_ = 0;
    int nSMtop_eee_ = 0;
    int nMET_eee_ = 0;
    
    int nTrigg_eeu_ = 0;
    int n3lep_eeu_ = 0;
    int nVetoMu_eeu_ = 0;
    int nVetoEl_eeu_ = 0;
    int nOS_eeu_ = 0;
    int nZmass_eeu_ = 0;
    int nJet_eeu_ = 0;
    int nBJet_eeu_ = 0;
    int nMWT_eeu_ = 0;
    int nSMtop_eeu_ = 0;
    int nMET_eeu_ = 0;
    
    int nTrigg_uuu_ = 0;
    int n3lep_uuu_ = 0;
    int nVetoMu_uuu_ = 0;
    int nVetoEl_uuu_ = 0;
    int nOS_uuu_ = 0;
    int nZmass_uuu_ = 0;
    int nJet_uuu_ = 0;
    int nBJet_uuu_ = 0;
    int nMWT_uuu_ = 0;
    int nSMtop_uuu_ = 0;
    int nMET_uuu_ = 0;
    
    int nTrigg_uue_ = 0;
    int n3lep_uue_ = 0;
    int nVetoMu_uue_ = 0;
    int nVetoEl_uue_ = 0;
    int nOS_uue_ = 0;
    int nZmass_uue_ = 0;
    int nJet_uue_ = 0;
    int nBJet_uue_ = 0;
    int nMWT_uue_ = 0;
    int nSMtop_uue_ = 0;
    int nMET_uue_ = 0;
    
    
    
    
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg", &nTrigg_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep", &n3lep_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu", &nVetoMu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl", &nVetoEl_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS", &nOS_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass",&nZmass_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet", &nJet_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet",&nBJet_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT", &nMWT_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop",&nSMtop_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET",&nMET_);
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg_eee", &nTrigg_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep_eee", &n3lep_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu_eee", &nVetoMu_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl_eee", &nVetoEl_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS_eee", &nOS_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass_eee",&nZmass_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet_eee", &nJet_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet_eee",&nBJet_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT_eee", &nMWT_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop_eee",&nSMtop_eee_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET_eee",&nMET_eee_);
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg_eeu", &nTrigg_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep_eeu", &n3lep_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu_eeu", &nVetoMu_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl_eeu", &nVetoEl_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS_eeu", &nOS_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass_eeu",&nZmass_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet_eeu", &nJet_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet_eeu",&nBJet_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT_eeu", &nMWT_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop_eeu",&nSMtop_eeu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET_eeu",&nMET_eeu_);
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg_uuu", &nTrigg_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep_uuu", &n3lep_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu_uuu", &nVetoMu_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl_uuu", &nVetoEl_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS_uuu", &nOS_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass_uuu",&nZmass_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet_uuu", &nJet_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet_uuu",&nBJet_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT_uuu", &nMWT_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop_uuu",&nSMtop_uuu_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET_uuu",&nMET_uuu_);
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg_uue", &nTrigg_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep_uue", &n3lep_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu_uue", &nVetoMu_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl_uue", &nVetoEl_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS_uue", &nOS_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass_uue",&nZmass_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet_uue", &nJet_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet_uue",&nBJet_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT_uue", &nMWT_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop_uue",&nSMtop_uue_);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET_uue",&nMET_uue_);
    
    
    
    double Xsect = 0.;
    Xsect = datasets[d]->Xsection();
    double globalSF = 1.;
    double nloSF = 1.;
    int nPos = 0;
    int nNeg = 0;
    int Ev = 0;
    int Weights = 0;
    double EquilumiSF = 1.0;
    if(!isData){
      int eventCounter = 0;
      for (int k = 0; k<globalnEntries; k++)
      {
        globalttree[(dataSetName).c_str()]->GetEntry(k);
        eventCounter += nEvents;
        
      }
      EquilumiSF = eventCounter / Xsect;
      cout << "                 equilumi = " <<  eventCounter <<" / " << Xsect <<" = " << EquilumiSF << endl;
    }
    
    
    
    
    
    if(isAMC){
      for (int k = 0; k<globalnEntries; k++)
      {
        globalttree[dataSetName.c_str()]->GetEntry(k);
        nPos += nPosW;
        nNeg += nNegW;
        Ev += nEvents;
        Weights += SumW;
      }
      nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
    }
    
    globalSF = globalSF * nloSF;
    
    //***********************************************RUNNING OVER EVENTS**********************************************
    int nTrigg_total = 0;
    int n3lep_total = 0;
    int nVetoMu_total = 0;
    int nVetoEl_total = 0;
    int nOS_total = 0;
    int nZmass_total = 0;
    int nSMtop_total = 0;
    int nJet_total = 0;
    int nBJet_total = 0;
    int nMWT_total = 0;
    int nMET_total = 0;
    
    int nTrigg_eee_total = 0;
    int n3lep_eee_total = 0;
    int nVetoMu_eee_total = 0;
    int nVetoEl_eee_total = 0;
    int nOS_eee_total = 0;
    int nZmass_eee_total = 0;
    int nSMtop_eee_total = 0;
    int nJet_eee_total = 0;
    int nBJet_eee_total = 0;
    int nMWT_eee_total = 0;
    int nMET_eee_total = 0;
    
    int nTrigg_eeu_total = 0;
    int n3lep_eeu_total = 0;
    int nVetoMu_eeu_total = 0;
    int nVetoEl_eeu_total = 0;
    int nOS_eeu_total = 0;
    int nZmass_eeu_total = 0;
    int nSMtop_eeu_total = 0;
    int nJet_eeu_total = 0;
    int nBJet_eeu_total = 0;
    int nMWT_eeu_total = 0;
    int nMET_eeu_total = 0;
    
    int nTrigg_uuu_total = 0;
    int n3lep_uuu_total = 0;
    int nVetoMu_uuu_total = 0;
    int nVetoEl_uuu_total = 0;
    int nOS_uuu_total = 0;
    int nZmass_uuu_total = 0;
    int nSMtop_uuu_total = 0;
    int nJet_uuu_total = 0;
    int nBJet_uuu_total = 0;
    int nMWT_uuu_total = 0;
    int nMET_uuu_total = 0;
    
    int nTrigg_uue_total = 0;
    int n3lep_uue_total = 0;
    int nVetoMu_uue_total = 0;
    int nVetoEl_uue_total = 0;
    int nOS_uue_total = 0;
    int nZmass_uue_total = 0;
    int nSMtop_uue_total = 0;
    int nJet_uue_total = 0;
    int nBJet_uue_total = 0;
    int nMWT_uue_total = 0;
    int nMET_uue_total = 0;
    
    
    for (int j = 0; j<globalnEntries; j++)
    {
      
      if(debug)
      {
        cin.get();
        cout << " " << endl;
        cout << "------------NEW File in Merged sample: " << j << " --------------" << endl;
      }
      globalttree[dataSetName.c_str()]->GetEntry(j);
      
      
      
      nTrigg_total = nTrigg_total + nTrigg_ ;
      n3lep_total = n3lep_total + n3lep_ ;
      nVetoMu_total = nVetoMu_total + nVetoMu_ ;
      nVetoEl_total = nVetoEl_total + nVetoEl_ ;
      nOS_total = nOS_total + nOS_ ;
      nZmass_total = nZmass_total + nZmass_ ;
      nJet_total = nJet_total + nJet_;
      nBJet_total = nBJet_total + nBJet_;
      nMWT_total = nMWT_total + nMWT_;
      nSMtop_total = nSMtop_total + nSMtop_;
      nMET_total = nMET_total + nMET_;
      
      nTrigg_eee_total = nTrigg_eee_total + nTrigg_eee_ ;
      n3lep_eee_total = n3lep_eee_total + n3lep_eee_ ;
      nVetoMu_eee_total = nVetoMu_eee_total + nVetoMu_eee_ ;
      nVetoEl_eee_total = nVetoEl_eee_total + nVetoEl_eee_ ;
      nOS_eee_total = nOS_eee_total + nOS_eee_ ;
      nZmass_eee_total = nZmass_eee_total + nZmass_eee_ ;
      nJet_eee_total = nJet_eee_total + nJet_eee_;
      nBJet_eee_total = nBJet_eee_total + nBJet_eee_;
      nMWT_eee_total = nMWT_eee_total + nMWT_eee_;
      nSMtop_eee_total = nSMtop_eee_total + nSMtop_eee_;
      nMET_eee_total = nMET_eee_total + nMET_eee_;
      
      nTrigg_eeu_total = nTrigg_eeu_total + nTrigg_eeu_ ;
      n3lep_eeu_total = n3lep_eeu_total + n3lep_eeu_ ;
      nVetoMu_eeu_total = nVetoMu_eeu_total + nVetoMu_eeu_ ;
      nVetoEl_eeu_total = nVetoEl_eeu_total + nVetoEl_eeu_ ;
      nOS_eeu_total = nOS_eeu_total + nOS_eeu_ ;
      nZmass_eeu_total = nZmass_eeu_total + nZmass_eeu_ ;
      nJet_eeu_total = nJet_eeu_total + nJet_eeu_;
      nBJet_eeu_total = nBJet_eeu_total + nBJet_eeu_;
      nMWT_eeu_total = nMWT_eeu_total + nMWT_eeu_;
      nSMtop_eeu_total = nSMtop_eeu_total + nSMtop_eeu_;
      nMET_eeu_total = nMET_eeu_total + nMET_eeu_;
      
      nTrigg_uuu_total = nTrigg_uuu_total + nTrigg_uuu_ ;
      n3lep_uuu_total = n3lep_uuu_total + n3lep_uuu_ ;
      nVetoMu_uuu_total = nVetoMu_uuu_total + nVetoMu_uuu_ ;
      nVetoEl_uuu_total = nVetoEl_uuu_total + nVetoEl_uuu_ ;
      nOS_uuu_total = nOS_uuu_total + nOS_uuu_ ;
      nZmass_uuu_total = nZmass_uuu_total + nZmass_uuu_ ;
      nJet_uuu_total = nJet_uuu_total + nJet_uuu_;
      nBJet_uuu_total = nBJet_uuu_total + nBJet_uuu_;
      nMWT_uuu_total = nMWT_uuu_total + nMWT_uuu_;
      nSMtop_uuu_total = nSMtop_uuu_total + nSMtop_uuu_;
      nMET_uuu_total = nMET_uuu_total + nMET_uuu_;
      
      nTrigg_uue_total = nTrigg_uue_total + nTrigg_uue_ ;
      n3lep_uue_total = n3lep_uue_total + n3lep_uue_ ;
      nVetoMu_uue_total = nVetoMu_uue_total + nVetoMu_uue_ ;
      nVetoEl_uue_total = nVetoEl_uue_total + nVetoEl_uue_ ;
      nOS_uue_total = nOS_uue_total + nOS_uue_ ;
      nZmass_uue_total = nZmass_uue_total + nZmass_uue_ ;
      nJet_uue_total = nJet_uue_total + nJet_uue_;
      nBJet_uue_total = nBJet_uue_total + nBJet_uue_;
      nMWT_uue_total = nMWT_uue_total + nMWT_uue_;
      nSMtop_uue_total = nSMtop_uue_total + nSMtop_uue_;
      nMET_uue_total = nMET_uue_total + nMET_uue_;
      
      
    }//for-loop events
    
    if(decaychannel.find("eee")!=string::npos) dChan = 3;
    if(decaychannel.find("eeu")!=string::npos) dChan = 2;
    if(decaychannel.find("uue")!=string::npos) dChan = 1;
    if(decaychannel.find("uuu")!=string::npos) dChan = 0;
    if(decaychannel.find("all")!=string::npos) dChan = 9;
    if(dChan == 9){
      cutFlowPlot->Fill(0, datasets[d], true, nTrigg_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(1, datasets[d], true, n3lep_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(2, datasets[d], true, nVetoMu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(3, datasets[d], true, nVetoEl_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(4, datasets[d], true, nOS_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(5, datasets[d], true, nZmass_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(6, datasets[d], true, nJet_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(7, datasets[d], true, nBJet_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(8, datasets[d], true, nMWT_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(9, datasets[d], true, nMET_total*Luminosity*globalSF/EquilumiSF);
      
      
      detailcutFlowPlot->Fill(0, datasets[d], true, n3lep_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(1, datasets[d], true, nVetoMu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(2, datasets[d], true, nVetoEl_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(3, datasets[d], true, nOS_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(4, datasets[d], true, nZmass_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(5, datasets[d], true, nJet_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(6, datasets[d], true, nBJet_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(7, datasets[d], true, nMWT_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(8, datasets[d], true, nMET_total*Luminosity*globalSF/EquilumiSF);
      
      
      
      detailcutFlowPlot3->Fill(0, datasets[d], true, nZmass_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(1, datasets[d], true, nJet_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(2, datasets[d], true, nBJet_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(3, datasets[d], true, nMWT_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(4, datasets[d], true, nMET_total*Luminosity*globalSF/EquilumiSF);
      
      selecTable.Fill(d, 0, nTrigg_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 1, n3lep_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 2, nVetoMu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 3, nVetoEl_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 4, nOS_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 5, nZmass_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 6, nJet_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 7, nBJet_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 8, nMWT_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 9, nMET_total *Luminosity*globalSF/EquilumiSF);
      
    }
    else if(dChan == 3) // eee
    {
      cutFlowPlot->Fill(0, datasets[d], true, nTrigg_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(1, datasets[d], true, n3lep_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(2, datasets[d], true, nVetoMu_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(3, datasets[d], true, nVetoEl_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(4, datasets[d], true, nOS_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(5, datasets[d], true, nZmass_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(6, datasets[d], true, nJet_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(7, datasets[d], true, nBJet_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(8, datasets[d], true, nMWT_eee_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(9, datasets[d], true, nMET_eee_total*Luminosity*globalSF/EquilumiSF);
      
      
      detailcutFlowPlot->Fill(0, datasets[d], true, n3lep_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(1, datasets[d], true, nVetoMu_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(2, datasets[d], true, nVetoEl_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(3, datasets[d], true, nOS_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(4, datasets[d], true, nZmass_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(5, datasets[d], true, nJet_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(6, datasets[d], true, nBJet_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(7, datasets[d], true, nMWT_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(8, datasets[d], true, nMET_eee_total*Luminosity*globalSF/EquilumiSF);
      
      
      
      detailcutFlowPlot3->Fill(0, datasets[d], true, nZmass_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(1, datasets[d], true, nJet_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(2, datasets[d], true, nBJet_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(3, datasets[d], true, nMWT_eee_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(4, datasets[d], true, nMET_eee_total*Luminosity*globalSF/EquilumiSF);
      
      selecTable.Fill(d, 0, nTrigg_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 1, n3lep_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 2, nVetoMu_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 3, nVetoEl_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 4, nOS_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 5, nZmass_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 6, nJet_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 7, nBJet_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 8, nMWT_eee_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 9, nMET_eee_total *Luminosity*globalSF/EquilumiSF);
      
    }
    else if(dChan == 2) // eeu
    {
      cutFlowPlot->Fill(0, datasets[d], true, nTrigg_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(1, datasets[d], true, n3lep_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(2, datasets[d], true, nVetoMu_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(3, datasets[d], true, nVetoEl_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(4, datasets[d], true, nOS_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(5, datasets[d], true, nZmass_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(6, datasets[d], true, nJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(7, datasets[d], true, nBJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(8, datasets[d], true, nMWT_eeu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(9, datasets[d], true, nMET_eeu_total*Luminosity*globalSF/EquilumiSF);
      
      
      detailcutFlowPlot->Fill(0, datasets[d], true, n3lep_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(1, datasets[d], true, nVetoMu_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(2, datasets[d], true, nVetoEl_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(3, datasets[d], true, nOS_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(4, datasets[d], true, nZmass_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(5, datasets[d], true, nJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(6, datasets[d], true, nBJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(7, datasets[d], true, nMWT_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(8, datasets[d], true, nMET_eeu_total*Luminosity*globalSF/EquilumiSF);
      
      
      
      detailcutFlowPlot3->Fill(0, datasets[d], true, nZmass_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(1, datasets[d], true, nJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(2, datasets[d], true, nBJet_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(3, datasets[d], true, nMWT_eeu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(4, datasets[d], true, nMET_eeu_total*Luminosity*globalSF/EquilumiSF);
      
      selecTable.Fill(d, 0, nTrigg_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 1, n3lep_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 2, nVetoMu_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 3, nVetoEl_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 4, nOS_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 5, nZmass_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 6, nJet_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 7, nBJet_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 8, nMWT_eeu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 9, nMET_eeu_total *Luminosity*globalSF/EquilumiSF);
      
    }
    else if(dChan==1) // uue
    {
      cutFlowPlot->Fill(0, datasets[d], true, nTrigg_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(1, datasets[d], true, n3lep_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(2, datasets[d], true, nVetoMu_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(3, datasets[d], true, nVetoEl_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(4, datasets[d], true, nOS_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(5, datasets[d], true, nZmass_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(6, datasets[d], true, nJet_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(7, datasets[d], true, nBJet_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(8, datasets[d], true, nMWT_uue_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(9, datasets[d], true, nMET_uue_total*Luminosity*globalSF/EquilumiSF);
      
      
      detailcutFlowPlot->Fill(0, datasets[d], true, n3lep_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(1, datasets[d], true, nVetoMu_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(2, datasets[d], true, nVetoEl_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(3, datasets[d], true, nOS_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(4, datasets[d], true, nZmass_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(5, datasets[d], true, nJet_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(6, datasets[d], true, nBJet_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(7, datasets[d], true, nMWT_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(8, datasets[d], true, nMET_uue_total*Luminosity*globalSF/EquilumiSF);
      
      
      
      detailcutFlowPlot3->Fill(0, datasets[d], true, nZmass_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(1, datasets[d], true, nJet_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(2, datasets[d], true, nBJet_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(3, datasets[d], true, nMWT_uue_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(4, datasets[d], true, nMET_uue_total*Luminosity*globalSF/EquilumiSF);
      
      selecTable.Fill(d, 0, nTrigg_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 1, n3lep_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 2, nVetoMu_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 3, nVetoEl_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 4, nOS_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 5, nZmass_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 6, nJet_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 7, nBJet_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 8, nMWT_uue_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 9, nMET_uue_total *Luminosity*globalSF/EquilumiSF);
      
    }
    else if(dChan == 0) /// uuu
    {
      cutFlowPlot->Fill(0, datasets[d], true, nTrigg_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(1, datasets[d], true, n3lep_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(2, datasets[d], true, nVetoMu_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(3, datasets[d], true, nVetoEl_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(4, datasets[d], true, nOS_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(5, datasets[d], true, nZmass_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(6, datasets[d], true, nJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(7, datasets[d], true, nBJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(8, datasets[d], true, nMWT_uuu_total*Luminosity*globalSF/EquilumiSF);
      cutFlowPlot->Fill(9, datasets[d], true, nMET_uuu_total*Luminosity*globalSF/EquilumiSF);
      
      
      detailcutFlowPlot->Fill(0, datasets[d], true, n3lep_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(1, datasets[d], true, nVetoMu_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(2, datasets[d], true, nVetoEl_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(3, datasets[d], true, nOS_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(4, datasets[d], true, nZmass_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(5, datasets[d], true, nJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(6, datasets[d], true, nBJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(7, datasets[d], true, nMWT_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot->Fill(8, datasets[d], true, nMET_uuu_total*Luminosity*globalSF/EquilumiSF);
      
      
      
      detailcutFlowPlot3->Fill(0, datasets[d], true, nZmass_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(1, datasets[d], true, nJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(2, datasets[d], true, nBJet_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(3, datasets[d], true, nMWT_uuu_total*Luminosity*globalSF/EquilumiSF);
      detailcutFlowPlot3->Fill(4, datasets[d], true, nMET_uuu_total*Luminosity*globalSF/EquilumiSF);
      
      selecTable.Fill(d, 0, nTrigg_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 1, n3lep_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 2, nVetoMu_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 3, nVetoEl_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 4, nOS_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 5, nZmass_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 6, nJet_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 7, nBJet_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 8, nMWT_uuu_total *Luminosity*globalSF/EquilumiSF);
      selecTable.Fill(d, 9, nMET_uuu_total *Luminosity*globalSF/EquilumiSF);
      
    }
    
    
    
    
    
    
    
    
    
  }//for-loop datasets
  
  
  
  cout << "unloading datasets" << endl;
  cout << "Writing Cut Flow Plot" << endl;
  string pathPNG = "MSPlotsCutflow/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += date +"/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += decaychannel;
  pathPNG += "/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  TFile *outfile = new TFile((pathPNG+"Baseline_CutFlow.root").c_str(),"recreate");
  outfile->cd();
  
  cutFlowPlot->setChannel(true, decaychannel);
  cutFlowPlot->setBins(CutsselecTable);
  cutFlowPlot->doCutFlowPlot();
  cutFlowPlot->Draw("Cut Flow", 1, false, false, false, 1);
  cutFlowPlot->Write(outfile, "CutFlow_"+decaychannel, true, pathPNG, "png");
  
  detailcutFlowPlot->setChannel(true, decaychannel);
  detailcutFlowPlot->setBins(CutsselecTabledetail);
  detailcutFlowPlot->doCutFlowPlot();
  detailcutFlowPlot->Draw("Cut Flow Detail", 1, false, false, false, 1);
  detailcutFlowPlot->Write(outfile, "CutFlowDetail_"+decaychannel, true, pathPNG, "png");
  
  detailcutFlowPlot3->setChannel(true, decaychannel);
  detailcutFlowPlot3->setBins(CutsselecTabledetail3);
  detailcutFlowPlot3->doCutFlowPlot();
  detailcutFlowPlot3->Draw("Cut Flow Detail Last 3", 1, false, false, false, 1);
  detailcutFlowPlot3->Write(outfile, "CutFlowDetail3_"+decaychannel, true, pathPNG, "png");
  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
  selecTable.TableCalculator(true, true, true, true, true);
  
  // Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies
  // (false), addTotalEfficiencies (false), writeLandscape (false)
  selecTable.Write(pathPNG + "CutFlow_Table_" + decaychannel + ".tex", false, true, true, true, false, false, true);
  
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  
}

// function that converts an int into a string
string intToStr (int number)
{
  ostringstream buff;
  buff<<number;
  return buff.str();
}

