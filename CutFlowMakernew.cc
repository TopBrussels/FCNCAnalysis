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
  
  
  MultiSamplePlot* cutFlowPlot = new MultiSamplePlot(datasets, "Cut Flow", 12, -0.5, 11.5, "Cut");
  MultiSamplePlot* detailcutFlowPlot = new MultiSamplePlot(datasets, "Cut Flow Detail", 12, -0.5, 11.5, "Cut");
  MultiSamplePlot* detailcutFlowPlot3 = new MultiSamplePlot(datasets, "Cut Flow Last 3 cuts", 12, -0.5, 11.5, "");
  
  CutsselecTable.clear();
  CutsselecTabledetail.clear();
 
  
  
  CutsselecTabledetail3.clear();
  
 /* CutsselecTabledetail3.push_back(string(">0 bjet"));
  CutsselecTabledetail.push_back(string("mWt"));
  CutsselecTabledetail3.push_back(string("SMtop"));
  CutsselecTabledetail3.push_back(string("MET"));
  */
  
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
    cout<<"Dataset:  :"<<dataSetName<<endl;
    if(dataSetName.find("tZq")==string::npos) continue;
    filepath = TreePath+"/"+dataSetName + ".root";
    filepath = "NtupleMakerOutput/Ntuples/Ntuples_170115/FCNC_3L_tZq_80X_1.root";

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
    
    

    
    
    vector<Int_t> *cutstep = 0; //
    vector<Int_t> *cutstep_eee =0;
    vector<Int_t> *cutstep_eeu = 0;
    vector<Int_t> *cutstep_uuu = 0;
    vector<Int_t> *cutstep_uue = 0 ;
    
    TBranch *bcutstep = 0; //
    TBranch *bcutstep_eee =0;
    TBranch *bcutstep_eeu = 0;
    TBranch *bcutstep_uuu = 0;
    TBranch *bcutstep_uue = 0 ;
    
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
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep",&cutstep, &bcutstep);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_eee",&cutstep_eee, &bcutstep_eee);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_eeu",&cutstep_eeu, &bcutstep_eeu);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_uuu",&cutstep_uuu, &bcutstep_uuu);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("cutstep_uue",&cutstep_uue, &bcutstep_uue);
   

/*
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    int n5 = 0;
    int n6 = 0;
    int n7 = 0;
    int n8 = 0;
    int n9 = 0;
    int n10 = 0;
    int n11 = 0;
    
    
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nTrigg", &n1);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("n3lep", &n2);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoMu", &n3);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nVetoEl", &n4);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nOS", &n5);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nZmass",&n6);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nJet", &n7);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nBJet",&n8);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMWT", &n9);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nSMtop",&n10);
    globalttree[(dataSetName).c_str()]->SetBranchAddress("nMET",&n11);
*/
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
/*    int n1b = 0;
    int n2b = 0;
    int n3b = 0;
    int n4b = 0;
    int n5b = 0;
    int n6b = 0;
    int n7b = 0;
    int n8b = 0;
    int n9b = 0;
    int n10b = 0;
    int n11b = 0;*/
    for (int j = 0; j<globalnEntries; j++)
    {
      
      if(debug)
      {
        cin.get();
        cout << " " << endl;
        cout << "------------NEW File in Merged sample: " << j << " --------------" << endl;
      }
      Long64_t tentry = globalttree[dataSetName.c_str()]->GetEntry(j);
      
      bcutstep->GetEntry(tentry);
      bcutstep_eee->GetEntry(tentry);
      bcutstep_eeu->GetEntry(tentry);
      bcutstep_uuu->GetEntry(tentry);
      bcutstep_uue->GetEntry(tentry);
      
      if( decaychannel.find("eee")!=std::string::npos ){
        for(int i = 0; i< CutsselecTable.size(); i++)
        {
          
          cout << "cut num " << i << " nb of events " << cutstep_eee->at(i) << endl;
          int cutstep_eee_input = cutstep_eee->at(i);
          cutFlowPlot->Fill(i, datasets[d], true, cutstep_eee_input*Luminosity*globalSF/EquilumiSF);
          selecTable.Fill(d, i, cutstep_eee_input*Luminosity*globalSF/EquilumiSF );
          
        }
        for(int i = 1; i< CutsselecTabledetail.size()+1; i++)
        {
          cout << "cut num " << i << " nb of events " << cutstep_eee->at(i) << endl;
          int cutstep_eee_input = cutstep_eee->at(i);
          detailcutFlowPlot->Fill(i-1, datasets[d], true, cutstep_eee_input*Luminosity*globalSF/EquilumiSF);
        }
        for(int i = 6; i< CutsselecTabledetail3.size()+6; i++)
        {
          cout << "cut num " << i << " nb of events " << cutstep_eee->at(i) << endl;
          int cutstep_eee_input = cutstep_eee->at(i);
          detailcutFlowPlot->Fill(i-6, datasets[d], true, cutstep_eee_input*Luminosity*globalSF/EquilumiSF);
        }
        
      }
      else if( decaychannel.find("eeu")!=std::string::npos ){
        for(int i = 0; i< CutsselecTable.size(); i++)
        {
          
          //cout << "cut num " << i << " nb of events " << cutstep_eeu->at(i) << endl;
          int cutstep_eeu_input = cutstep_eeu->at(i);
          cutFlowPlot->Fill(i, datasets[d], true, cutstep_eeu_input*Luminosity*globalSF/EquilumiSF);
          
          selecTable.Fill(d, i, cutstep_eeu_input*Luminosity*globalSF/EquilumiSF );
          
        }
        for(int i = 1; i< CutsselecTabledetail.size()+1; i++)
        {
          //cout << "cut num " << i << " nb of events " << cutstep_eeu->at(i) << endl;
          int cutstep_eeu_input = cutstep_eeu->at(i);
          detailcutFlowPlot->Fill(i-1, datasets[d], true, cutstep_eeu_input*Luminosity*globalSF/EquilumiSF);
        }
        for(int i = 6; i< CutsselecTabledetail3.size()+6; i++)
        {
          //cout << "cut num " << i << " nb of events " << cutstep_eeu->at(i) << endl;
          int cutstep_eeu_input = cutstep_eeu->at(i);
          detailcutFlowPlot3->Fill(i-6, datasets[d], true, cutstep_eeu_input*Luminosity*globalSF/EquilumiSF);
        }
        
      }
      else if( decaychannel.find("uuu")!=std::string::npos ){
        for(int i = 0; i< CutsselecTable.size(); i++)
        {
          
          cout << "cut num " << i << " nb of events " << cutstep_uuu->at(i) << endl;
          int cutstep_uuu_input = cutstep_uuu->at(i);
          cutFlowPlot->Fill(i, datasets[d], true, cutstep_uuu_input*Luminosity*globalSF/EquilumiSF);
          selecTable.Fill(d, i, cutstep_uuu_input*Luminosity*globalSF/EquilumiSF );
          
        }
        for(int i = 1; i< CutsselecTabledetail.size()+1; i++)
        {
          cout << "cut num " << i << " nb of events " << cutstep_uuu->at(i) << endl;
          int cutstep_uuu_input = cutstep_uuu->at(i);
          detailcutFlowPlot->Fill(i-1, datasets[d], true, cutstep_uuu_input*Luminosity*globalSF/EquilumiSF);
        }
        for(int i = 6; i< CutsselecTabledetail3.size()+6; i++)
        {
          cout << "cut num " << i << " nb of events " << cutstep_uuu->at(i) << endl;
          int cutstep_uuu_input = cutstep_uuu->at(i);
          detailcutFlowPlot3->Fill(i-6, datasets[d], true, cutstep_uuu_input*Luminosity*globalSF/EquilumiSF);
        }
        
      }
      else if( decaychannel.find("uue")!=std::string::npos ){
        for(int i = 0; i< CutsselecTable.size(); i++)
        {
          
          cout << "cut num " << i << " nb of events " << cutstep_uue->at(i) << endl;
          int cutstep_uue_input = cutstep_uue->at(i);
          cutFlowPlot->Fill(i, datasets[d], true, cutstep_uue_input*Luminosity*globalSF/EquilumiSF);
          selecTable.Fill(d, i, cutstep_uue_input*Luminosity*globalSF/EquilumiSF );
          
        }
        for(int i = 1; i< CutsselecTabledetail.size()+1; i++)
        {
          //cout << "cut num " << i << " nb of events " << cutstep_uue->at(i) << endl;
          int cutstep_uue_input = cutstep_uue->at(i);
          detailcutFlowPlot->Fill(i-1, datasets[d], true, cutstep_uue_input*Luminosity*globalSF/EquilumiSF);
        }
        for(int i = 6; i< CutsselecTabledetail3.size()+6; i++)
        {
          //cout << "cut num " << i << " nb of events " << cutstep_uue->at(i) << endl;
          int cutstep_uue_input = cutstep_uue->at(i);
          detailcutFlowPlot3->Fill(i-6, datasets[d], true, cutstep_uue_input*Luminosity*globalSF/EquilumiSF);
        }
        
      }
      
      /*n1b = n1b + n1 ;
      n2b = n2b + n2 ;
      n3b = n3b + n3 ;
      n4b = n4b + n4 ;
      n5b = n5 + n5b ;
      n6b = n6 + n6b ;
      n7b = n7b + n7;
      n8b = n8b + n8;
      n9b = n9b + n9;
      n10b = n10b + n10;
      n11b = n11b + n11;*/
      
    }//for-loop events
    globalttree[(dataSetName).c_str()]->ResetSetBranchAddress();

    /*
    cutFlowPlot->Fill(0, datasets[d], true, n1b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(1, datasets[d], true, n2b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(2, datasets[d], true, n3b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(3, datasets[d], true, n4b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(4, datasets[d], true, n6b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(5, datasets[d], true, n7b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(6, datasets[d], true, n8b*Luminosity*globalSF/EquilumiSF);
   // cutFlowPlot->Fill(7, datasets[d], true, n9b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(7, datasets[d], true, n10b*Luminosity*globalSF/EquilumiSF);
    cutFlowPlot->Fill(8, datasets[d], true, n11b*Luminosity*globalSF/EquilumiSF);
    */
    
    
 /*   detailcutFlowPlot->Fill(0, datasets[d], true, n2b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(1, datasets[d], true, n3b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(2, datasets[d], true, n4b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(3, datasets[d], true, n6b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(4, datasets[d], true, n7b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(5, datasets[d], true, n8b*Luminosity*globalSF/EquilumiSF);
   // detailcutFlowPlot->Fill(7, datasets[d], true, n9b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(6, datasets[d], true, n10b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot->Fill(7, datasets[d], true, n11b*Luminosity*globalSF/EquilumiSF);
    
    detailcutFlowPlot3->Fill(0, datasets[d], true, n8b*Luminosity*globalSF/EquilumiSF);
    // detailcutFlowPlot->Fill(7, datasets[d], true, n9b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot3->Fill(1, datasets[d], true, n10b*Luminosity*globalSF/EquilumiSF);
    detailcutFlowPlot3->Fill(2, datasets[d], true, n11b*Luminosity*globalSF/EquilumiSF);
    */
    /*
    selecTable.Fill(d, 0, n1b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 1, n2b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 2, n3b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 3, n4b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 4, n6b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 5, n7b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 6, n8b *Luminosity*globalSF/EquilumiSF);
   // selecTable.Fill(d, 7, n9b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 7, n10b *Luminosity*globalSF/EquilumiSF);
    selecTable.Fill(d, 8, n11b *Luminosity*globalSF/EquilumiSF);
    */
    
  }//for-loop datasets
  
  
  
  cout << "unloading datasets" << endl;
  cout << "Writing Cut Flow Plot" << endl;
  string pathPNG = "MSPlotsCutflow/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MSPlots";
  pathPNG += decaychannel;
  pathPNG += "/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  TFile *outfile = new TFile((pathPNG+"Baseline_CutFlow.root").c_str(),"recreate");
  outfile->cd();
  
  
  cutFlowPlot->setBins(CutsselecTable);
  cutFlowPlot->doCutFlowPlot();
  cutFlowPlot->Draw("Cut Flow", 1, false, false, false, 1);
  cutFlowPlot->Write(outfile, "CutFlow_"+decaychannel, true, pathPNG, "png");
  
  detailcutFlowPlot->setBins(CutsselecTabledetail);
  detailcutFlowPlot->doCutFlowPlot();
  detailcutFlowPlot->Draw("Cut Flow Detail", 1, false, false, false, 1);
  detailcutFlowPlot->Write(outfile, "CutFlowDetail_"+decaychannel, true, pathPNG, "png");
  
  detailcutFlowPlot3->setBins(CutsselecTabledetail3);
  detailcutFlowPlot3->doCutFlowPlot();
  detailcutFlowPlot3->Draw("Cut Flow Detail Last 3", 1, false, false, false, 1);
  detailcutFlowPlot3->Write(outfile, "CutFlowDetail3_"+decaychannel, true, pathPNG, "png");
  //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
  selecTable.TableCalculator(true, true, true, true, true);
  
  // Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies
  // (false), addTotalEfficiencies (false), writeLandscape (false)
  selecTable.Write(pathPNG + "CutFlow_Table" + decaychannel + ".tex", false, true, true, true, false, false, true);
  
  
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

