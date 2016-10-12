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
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"


//inlcludes for TMVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"


using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TFile*> globalFileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,TNtuple*> ntree;
map<string,TNtuple*> otree;
map<string,TTree*> globalttree;
map<string,MultiSamplePlot*> MSPlot;


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG);
void MSPCreator (string pathPNG);
void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );
void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string BkgName ,std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel);
string ConvertIntToString(int nb, bool pad);
string MakeTimeStamp();





// CONFIGURATION
Bool_t debug = true;
bool mumumu  = false;
bool eee = false;
string channelpostfix = "";
double DataLumi = -1;
bool elecPlot = false;
bool muPlot = false;
bool doScaling = false;
//applying all appropriate scale factors for individual objects if the bool is set to true

Bool_t train_mva = false;
string dateString;

int main(int argc, char* argv[])
{
  if (debug){
    cout << "argc = " << argc << endl;
    for(int i = 0; i < argc; i++)
    {
      cout << "argv[" << i << "] = " << argv[i] << endl;
    }
  }
  
  
  
  //Placing arguments in properly typed variables
  const string channel = argv[1];
  debug = false;
  train_mva = false;
  debug = strtol(argv[2],NULL,10);
  train_mva = strtol(argv[3],NULL,10);
  const string signalName = argv[4];
  const string bkgName = argv[5];
  doScaling =  strtol(argv[6],NULL,10);
  
  
  string xmlFileName;
  string CraneenPath;
  CraneenPath = "NtupleMakerOutput/MergedTuples/";
  if(channel=="MuMuMu")
  {
    cout << " --> Using the TriMuon channel..." << endl;
    channelpostfix = "_mumumu";
    xmlFileName = "config/Run2TriLepton_samples_analyzer_mumumu.xml";
    mumumu = true;
    eee = false;
    DataLumi  = 2100 ;//2612.180735004;//  pb-1
    CraneenPath += "mumumu/";
  }
  else if(channel=="ElElEl")
  {
    cout << " --> Using the TriElectron channel..." << endl;
    channelpostfix = "_eee";
    xmlFileName = "config/Run2TriLepton_samples_analyzer_eee.xml";
    mumumu = false;
    eee = true;
    DataLumi  = 2612.180735004;//  pb-1
    CraneenPath += "eee/";
  }
  else if(channel=="All")
  {
    cout << " --> Using the all channel..." << endl;
    channelpostfix = "_all";
    xmlFileName = "config/Run2TriLepton_samples_analy.xml";
    mumumu = false;
    eee = false;
    DataLumi  = 2612.180735004;//  pb-1
    CraneenPath += "all/";
  }
  else
  {
    cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
    exit(1);
  }
  dateString = MakeTimeStamp();
  //    CraneenPath += dateString + "/";
  CraneenPath += "160718/";
  string pathPNG = "myOutput";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "/" + dateString + "/";
  mkdir(pathPNG.c_str(),0777);
  pathPNG += "MVAPlots"+channelpostfix+"/";
  mkdir(pathPNG.c_str(),0777);
  cout <<"Making directory :"<< pathPNG  <<endl;
  
  cout << "xmlFileName for training /evaluating is " << xmlFileName << endl;
  cout << "NtupleFolder is " << CraneenPath << endl;
  if(debug) cout << "debugging on" << endl;
  
  
  // calling the function that writtes all the MSPlots in a root file
  //void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel)
  MVAanalysis(train_mva, "BDT", 2, signalName, bkgName,xmlFileName, xmlFileName, CraneenPath, channel); // divide sample in 2
  if(!train_mva) MSPCreator (pathPNG);
  
}



// function that writes all the MSPlots created in a root file
void MSPCreator (string pathPNG)
{
  // cout << pathPNG.c_str() << endl;
  
  TFile *outfile = new TFile((pathPNG+"/Output_MVA.root").c_str(),"recreate");
  outfile->cd();
  cout << "created " << (pathPNG+"/Output_MVA.root").c_str() << endl;
  
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
  {
    
    string name = "MyMSP_" + it->first;
    cout << " name " << name << endl;
    MultiSamplePlot *temp = it->second;
    if (debug){
      cout << "Saving the MSP" << endl;
      cout << " and it->first is " << it->first << endl;
    }
    temp->Draw("MyMSP"+it->first, 1, false, false, false, 10);// 0 = no ratio
    //    name += "_3L";
    
    cout << "name " << name << endl;
    temp->Write(outfile, name, true,pathPNG.c_str() , "png");
    //      vector<string> temp_histo = it->GetTH1FNames();
    //      for (int i_hist=0; i_hist < temp_histo.size();i_hist++  ){
    //	cout << "hist is" << temp_histo[i_hist] << endl;
    //	cout << "integral is " << it->GetTH1F(temp_histo[i_hist].GetSum()) << endl;
    //      }
  }
  
  outfile->Write("kOverwrite");
}


void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName,std::string BkgName,std::string xmlNom_train,std::string xmlNom_evaluate , TString CraneenPath, std::string channel)
{

  // Set cuts
  //
  TString set_MET_cut = ">30";
  TString set_mTW_cut = "";
  TString set_NJets_cut = ">1"; //ONLY STRICT SIGN (> / < / ==)
  TString set_NBJets_cut = "=1"; //ONLY STRICT SIGN (> / < / ==)
  
  std::vector<TString > CUTvars;
  CUTvars.push_back("MET"); 
  
  
  MVAComputer* Eventcomputer_ =0;
  MVATrainer* Eventtrainer_ = 0;
  if(doTraining) Eventtrainer_ = new MVATrainer(MVAmethod,"TrainedEventMVA"+channel, "TrainedEventMVA"+channel+".root");
  vector<std::string> MVAvars;
  if(debug) cout << "event trainier initialised " << endl;
  
  // the name of the variables to be used, order is important!
  //  MVAvars.push_back("NumberOfElectrons"); // doesn't contribute
  //  MVAvars.push_back("NumberOfMuons"); // doesn't contribute
  //  MVAvars.push_back("Zmass"); // doesn't contribute
  //  MVAvars.push_back("TrMassW"); lowcontribution ? 
//  MVAvars.push_back("MET");
  //  MVAvars.push_back("CvsL_1");
  // MVAvars.push_back("CvsB_1");
  //  MVAvars.push_back("CvsL_2");
  //  MVAvars.push_back("CvsB_2");
  MVAvars.push_back("pt_electron_1");
  //  MVAvars.push_back("pt_electron_2");
  //  MVAvars.push_back("pt_electron_3");
  MVAvars.push_back("pt_muon_1");
  //  MVAvars.push_back("pt_muon_2");
  // MVAvars.push_back("pt_muon_3");
//  MVAvars.push_back("pt_jet_1"); // related wit FNC top  mass
  //  MVAvars.push_back("pt_jet_2");
  MVAvars.push_back("topMass");
  //  MVAvars.push_back("nCSVL");
  MVAvars.push_back("nCSVM");
  //  MVAvars.push_back("nCSVT");
//  MVAvars.push_back("bdis_1");
  MVAvars.push_back("bdis_2");
  MVAvars.push_back("FCNCtopmass");
 // MVAvars.push_back("Pt_cjet");
//  MVAvars.push_back("deltaPhiSMFCNCtop"); // related with FCNC tiop mass
//  MVAvars.push_back("deltaPhiWlepb");
//  MVAvars.push_back("deltaPhiWlepc");
//  MVAvars.push_back("deltaPhiZc"); // related with FCNC top mass 
  MVAvars.push_back("deltaPhiZb");
  //MVAvars.push_back("deltaRSMFCNCtop"); // correlated with deltaPhiWlepC deltaPhiZc
  MVAvars.push_back("deltaRWlepb");
//  MVAvars.push_back("deltaRWlepc");
 // MVAvars.push_back("deltaRZc"); // related with fcnc top mass
  MVAvars.push_back("deltaRZb");
  MVAvars.push_back("MassWlepB");
  
  
  if(doTraining){
    for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
    {
      Eventtrainer_->bookInputVar(MVAvars[N_var]);
    }
    Eventtrainer_->bookWeight("Weight");
    if(debug) cout << "input variables booked" << endl;
  }
  if(!doTraining) Eventcomputer_ = new MVAComputer(MVAmethod,"TrainedEventMVA"+channel+".root", "TrainedEventMVA"+channel,MVAvars, "test");
  
  
  
  cout<<""<<endl;
  cout<<"RUNNING MVA: ";
  if(doTraining) cout << " training " <<endl;
  else cout << "evaluating " << endl;
  cout<<""<<endl;
  
  const char *xmlfile_train = xmlNom_train.c_str();
  const char *xmlfile_eval = xmlNom_evaluate.c_str();
  if(doTraining) cout << "used config file: " << xmlfile_train << endl;
  else cout << "used config file: " << xmlfile_eval << endl;
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  if(doTraining) treeLoader.LoadDatasets (datasets, xmlfile_train);	//cout<<"datasets loaded"<<endl;
  else treeLoader.LoadDatasets(datasets, xmlfile_eval);
  
  //***************************************************CREATING PLOT****************************************************
  //  MSPlot[MVAmethod.c_str()] = new MultiSamplePlot(datasets, MVAmethod.c_str() , 50, (float) -1, (float) 1, MVAmethod.c_str());
  //  MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
  
  
  
  
  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;
  float ScaleFactor, NormFactor;
  
  
  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
  {
    dataSetName = datasets[d]->Name();
    cout<<"Dataset:  :"<<dataSetName<<endl;
    filepath = CraneenPath+ dataSetName + ".root";
    if (debug) cout<<"filepath: "<<filepath<<endl;
    bool isData= false;
    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
    {
      if(debug) cout << "Data found" << endl;
      isData =true;
    }
    bool isAMC= false;
    if(dataSetName.find("amc")!=string::npos || dataSetName.find("AMC")!=string::npos )
    {
      if(debug) cout << "AMC found" << endl;
      isAMC =true;
    }
    
    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"update"); //create TFile for each dataset
    //   globalFileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"update"); //create TFile for each dataset
    
    string TTreename = "tree"; // tree containing the variablesi
    //   string sglobaltree = "globaltree";
    
    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
    nEntries = ttree[dataSetName.c_str()]->GetEntries();
    cout<<"                 nEntries: "<<nEntries<<endl;
    
    
    
    //    globalttree[dataSetName.c_str()] = (TTree*)globalFileObj[dataSetName.c_str()]->Get(sglobaltree.c_str());
    //    int globalnEntries = (int)globalttree[dataSetName.c_str()]->GetEntries();
    //    cout<<"		    nEntries gt: "<<globalnEntries<<endl;
    if(debug) cout << " trees initialised" << endl;
    
    
    // get the SF from the corresponding branch
    Int_t PassedMET = 0;
    ttree[dataSetName.c_str()]->SetBranchAddress("PassedMETFilter",&PassedMET);
    
    Double_t	puSF	=    1.    ;
    ttree[dataSetName.c_str()]->SetBranchAddress("puSF",&puSF);
    
    Double_t	nloW;
    ttree[dataSetName.c_str()]->SetBranchAddress("nloWeight",&nloW);
    
    Double_t	electronSF[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("ElectronSF",&electronSF);
    
    Double_t	muonID[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonIDSF",	&muonID);
    
    Double_t	muonIso[10];
    ttree[dataSetName.c_str()]->SetBranchAddress("MuonIsoSF",	 &muonIso);
    
    //Int_t    nPosW;
    //globalttree[dataSetName.c_str()]->SetBranchAddress("nofPosWeights",&nPosW);
    
    // Int_t    nNegW;
    // globalttree[dataSetName.c_str()]->SetBranchAddress("nofNegWeights",&nNegW);
    
    Double_t	BSF;
    ttree[dataSetName.c_str()]->SetBranchAddress("btagSF",&BSF);
    
    
    double nloSF = 1.;
    //    int nPos = 0;
    //    int nNeg = 0;
    
    /*   if(isAMC && !isData)
     {
     for (int k = 0; k<globalnEntries; k++)
     {
     globalttree[(dataSetName).c_str()]->GetEntry(k);
     nPos += nPosW;
     nNeg += nNegW;
     }
     //          if(!isData) nloSF *= (double) Weights/(double) Ev; //
     nloSF *= ((double) (nPos + nNeg))/((double) (nPos - nNeg));
     cout << "                nloSF: " << nloSF << endl;
     }
     */
    if(dataSetName.find("tZq_amc")!=string::npos) nloSF = 3.76316;
    if(dataSetName.find("WZJets_amc")!=string::npos) nloSF = 1.52061;
    if(dataSetName.find("ST_tamc")!=string::npos) nloSF = 4.64013;
    if(dataSetName.find("TTWJetsToLNu_amc")!=string::npos) nloSF = 1.94035;
    if(dataSetName.find("TTZToQQ_amc")!=string::npos) nloSF = 2.13364;
    if(dataSetName.find("TTZToLLNuNu_amc")!=string::npos) nloSF = 2.15175;
    if(dataSetName.find("Zjets(0amc")!=string::npos) nloSF = 1.49209;
    /////////////////////////////////////////
    // Define variables relevant for MVA
    ////////////////////////////////////////
    int NumberOfElectrons, NumberOfMuons, nJets,nCSVL, nCSVM, nCSVT;
    double TrMassW, MET, CvsL_1, CvsL[20], CvsB_1, CvsB[20], pt_electron_1, pt_electron_2, pt_electron_3, pt_muon_1, pt_muon_2, pt_muon_3, pt_jet_1,pt_jet_2,topMass, Zmass, bDisc[20], FCNCtopmass, Pt_cjet, deltaPhiWlepb, deltaPhiWlepc, deltaPhiZc, deltaPhiZb, deltaPhiSMFCNCtop,deltaRWlepb, deltaRWlepc,    deltaRZc, deltaRZb, deltaRSMFCNCtop, MassWlepB ;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("nElectrons",&NumberOfElectrons);
    ttree[(dataSetName).c_str()]->SetBranchAddress("nMuons",&NumberOfMuons);
    ttree[(dataSetName).c_str()]->SetBranchAddress("Zboson_M", &Zmass);
    ttree[(dataSetName).c_str()]->SetBranchAddress("mWt", &TrMassW);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Pt", &MET);
    ttree[(dataSetName).c_str()]->SetBranchAddress("nJets", &nJets);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet_1", &CvsL_1);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet_1", &CvsB_1);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet", &CvsL);
    //    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet[nJets]", &CvsB);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron_1", &pt_electron_1);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron_2", &pt_electron_2);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron_3", &pt_electron_3);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon_1", &pt_muon_1);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon_2", &pt_muon_2);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon_3", &pt_muon_3);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet_1", &pt_jet_1);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet_2", &pt_jet_2);
    ttree[(dataSetName).c_str()]->SetBranchAddress("SMtop_M", &topMass);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet", &bDisc);
    //   ttree[(dataSetName).c_str()]->SetBranchAddress("nJets_CSVL", &nCSVL);
    ttree[(dataSetName).c_str()]->SetBranchAddress("nJets_CSVM", &nCSVM);
    ttree[(dataSetName).c_str()]->SetBranchAddress("nJets_CSVT", &nCSVT);
    ttree[(dataSetName).c_str()]->SetBranchAddress("FCNCtop_M", &FCNCtopmass);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cjet_Pt", &Pt_cjet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dPhiSMFCNCtop",&deltaPhiSMFCNCtop);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dPhiWlepb",&deltaPhiWlepb);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dPhiWlepc",&deltaPhiWlepc);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dPhiZc",&deltaPhiZc);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dPhiZb",&deltaPhiZb);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dRSMFCNCtop",&deltaRSMFCNCtop);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dRWlepb",&deltaRWlepb);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dRWlepc",&deltaRWlepc);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dRZc",&deltaRZc);
    ttree[(dataSetName).c_str()]->SetBranchAddress("dRZb",&deltaRZb);
    ttree[(dataSetName).c_str()]->SetBranchAddress("mlb",&MassWlepB);
    
    
    
    
    
    
    //    int Train_nEntries = int(nEntries/skipEvents);
    //    if(isData && !doTraining) Train_nEntries = int(nEntries/skipEvents);
    
    cout << "Number of entries: " << nEntries << endl; //", number of train Entries: " << Train_nEntries << endl;
    
    //////////////////////////////////////////////////////////
    // Running on events
    //////////////////////////////////////////////////////////
    
    if(doTraining)
    {
      cout << "signal is " << SignalName << endl;
      cout << "bkg is " << BkgName << endl;
      
      if(isData && SignalName != "Data") continue;
      
      for (int j = 0; j<nEntries; j++)
      {
        ScaleFactor = 1.; // event scale factor
        ttree[(dataSetName).c_str()]->GetEntry(j);
        
        if(doScaling){
          if(!PassedMET){continue;}
          if(debug) cout << "MET SF " << endl;
          //electron
          for(unsigned int iEl = 0; iEl < NumberOfElectrons ; iEl ++)
          {
            //cout << "entry " << j << " electron " << iEl << " SF " <<  electronSF[iEl]<< endl;
            ScaleFactor *= electronSF[iEl];
          }
          if(debug) cout << "electron SF " << endl;
          //muon
          for(unsigned int iMu = 0; iMu < NumberOfMuons ; iMu ++)
          {
            ScaleFactor *= muonID[iMu]*muonIso[iMu];
            
          }
          if(debug) cout << "muon SF " << endl;
          //PU
          ScaleFactor *= puSF;
          if(debug) cout << "PU SF " << endl;
          //btag
          ScaleFactor *= BSF;
          if(debug) cout << "btag SF " << endl;
          //AMC
          ScaleFactor *= nloSF * nloW;
          if(debug) cout << "amc SF " << endl;
          
          if(dataSetName.find("NP")!=string::npos){ ScaleFactor = 1.;}
        }
        NormFactor = DataLumi / datasets[d]->EquivalentLumi(); //data->NormFactor*Lumi
        ScaleFactor *= NormFactor;
        if(debug) cout << "normalisation SF " << endl;
        if(dataSetName.find("NP")!=string::npos){ ScaleFactor = 1.;} 
        if(ScaleFactor < 0 ) ScaleFactor *= -1;  
        if(dataSetName.find("NP_overlay")!=string::npos )
        {
          
          Eventtrainer_->FillWeight("S","Weight", ScaleFactor);
          //         Eventtrainer_->Fill("S","NumberOfElectrons", NumberOfElectrons);
          //          Eventtrainer_->Fill("S","NumberOfMuons", NumberOfMuons);
          //	  Eventtrainer_->Fill("S","Zmass", Zmass);
//          Eventtrainer_->Fill("S","TrMassW", TrMassW);
//          Eventtrainer_->Fill("S","MET", MET );
          //	  Eventtrainer_->Fill("S","CvsL_1", CvsL_1 );
          //  Eventtrainer_->Fill("S","CvsB_1", CvsB_1 );
          //	  Eventtrainer_->Fill("S","CvsL_2", CvsL[1] );
          //	  Eventtrainer_->Fill("S","CvsB_2", CvsB[1] );
          Eventtrainer_->Fill("S","pt_electron_1", pt_electron_1 );
          //	  Eventtrainer_->Fill("S","pt_electron_2", pt_electron_2 );
          //	  Eventtrainer_->Fill("S","pt_electron_3", pt_electron_3 );
          Eventtrainer_->Fill("S","pt_muon_1", pt_muon_1 );
          //	  Eventtrainer_->Fill("S","pt_muon_2", pt_muon_2 );
          //	  Eventtrainer_->Fill("S","pt_muon_3", pt_muon_3 );
//          Eventtrainer_->Fill("S","pt_jet_1", pt_jet_1 );
          //	  Eventtrainer_->Fill("S","pt_jet_2", pt_jet_2 );
          Eventtrainer_->Fill("S","topMass", topMass );
          //	  Eventtrainer_->Fill("S","nCSVL", nCSVL);
          Eventtrainer_->Fill("S","nCSVM", nCSVM);
          //	  Eventtrainer_->Fill("S","nCSVT", nCSVT);
//          Eventtrainer_->Fill("S","bdis_1", bDisc[0]);
          Eventtrainer_->Fill("S","bdis_2", bDisc[1]);
          Eventtrainer_->Fill("S","FCNCtopmass",FCNCtopmass);
//          Eventtrainer_->Fill("S","Pt_cjet",Pt_cjet);
//          Eventtrainer_->Fill("S","deltaPhiSMFCNCtop",deltaPhiSMFCNCtop);
//          Eventtrainer_->Fill("S","deltaPhiWlepb",deltaPhiWlepb);
//          Eventtrainer_->Fill("S","deltaPhiWlepc",deltaPhiWlepc);
//          Eventtrainer_->Fill("S","deltaPhiZc",deltaPhiZc);
          Eventtrainer_->Fill("S","deltaPhiZb",deltaPhiZb);
//          Eventtrainer_->Fill("S","deltaRSMFCNCtop",deltaRSMFCNCtop);
          Eventtrainer_->Fill("S","deltaRWlepb",deltaRWlepb);
//          Eventtrainer_->Fill("S","deltaRWlepc",deltaRWlepc);
//          Eventtrainer_->Fill("S","deltaRZc",deltaRZc);
          Eventtrainer_->Fill("S","deltaRZb",deltaRZb);
         Eventtrainer_->Fill("S","MassWlepB",MassWlepB);
        }
        else if((dataSetName.find("TTZ")!=string::npos || dataSetName.find("WZJets")!=string::npos || dataSetName.find("tZq")!=string::npos) && BkgName.find("all")==string::npos )
        {
          //cout << "train against 1 bkg" << endl;
          Eventtrainer_->FillWeight("B","Weight", ScaleFactor);
          //          Eventtrainer_->Fill("B","NumberOfElectrons", NumberOfElectrons);
          //          Eventtrainer_->Fill("B","NumberOfMuons", NumberOfMuons);
          //	  Eventtrainer_->Fill("B","Zmass", Zmass);
//          Eventtrainer_->Fill("B","TrMassW", TrMassW);
//          Eventtrainer_->Fill("B","MET", MET );
          //	  Eventtrainer_->Fill("B","CvsL_1", CvsL_1 );
          // Eventtrainer_->Fill("B","CvsB_1", CvsB_1 );
          //	  Eventtrainer_->Fill("B","CvsL_2", CvsL[1] );
          //	  Eventtrainer_->Fill("B","CvsB_2", CvsB[1] );
          Eventtrainer_->Fill("B","pt_electron_1", pt_electron_1 );
          //	  Eventtrainer_->Fill("B","pt_electron_2", pt_electron_2 );
          //	  Eventtrainer_->Fill("B","pt_electron_3", pt_electron_3 );
          Eventtrainer_->Fill("B","pt_muon_1", pt_muon_1 );
          //	  Eventtrainer_->Fill("B","pt_muon_2", pt_muon_2 );
          //	  Eventtrainer_->Fill("B","pt_muon_3", pt_muon_3 );
//          Eventtrainer_->Fill("B","pt_jet_1", pt_jet_1 );
          //	  Eventtrainer_->Fill("B","pt_jet_2", pt_jet_2 );
          Eventtrainer_->Fill("B","topMass", topMass );
          //	  Eventtrainer_->Fill("B","nCSVL", nCSVL);
          Eventtrainer_->Fill("B","nCSVM", nCSVM);
          //	  Eventtrainer_->Fill("B","nCSVT", nCSVT);
//          Eventtrainer_->Fill("B","bdis_1", bDisc[0]);
          Eventtrainer_->Fill("B","bdis_2", bDisc[1]);
          Eventtrainer_->Fill("B","FCNCtopmass",FCNCtopmass);
//          Eventtrainer_->Fill("B","Pt_cjet",Pt_cjet);
//          Eventtrainer_->Fill("B","deltaPhiSMFCNCtop",deltaPhiSMFCNCtop);
//         Eventtrainer_->Fill("B","deltaPhiWlepb",deltaPhiWlepb);
//          Eventtrainer_->Fill("B","deltaPhiWlepc",deltaPhiWlepc);
//          Eventtrainer_->Fill("B","deltaPhiZc",deltaPhiZc);
          Eventtrainer_->Fill("B","deltaPhiZb",deltaPhiZb);
//          Eventtrainer_->Fill("B","deltaRSMFCNCtop",deltaRSMFCNCtop);
          Eventtrainer_->Fill("B","deltaRWlepb",deltaRWlepb);
//          Eventtrainer_->Fill("B","deltaRWlepc",deltaRWlepc);
//          Eventtrainer_->Fill("B","deltaRZc",deltaRZc);
          Eventtrainer_->Fill("B","deltaRZb",deltaRZb);
          Eventtrainer_->Fill("B","MassWlepB",MassWlepB);
          
        }
        else if( BkgName.find("all")!=string::npos)
        {
          //cout << "train against all bkg" << endl;
          Eventtrainer_->FillWeight("B","Weight", ScaleFactor);
          //          Eventtrainer_->Fill("B","NumberOfElectrons", NumberOfElectrons);
          //          Eventtrainer_->Fill("B","NumberOfMuons", NumberOfMuons);
          //	  Eventtrainer_->Fill("B","Zmass", Zmass);
          Eventtrainer_->Fill("B","TrMassW", TrMassW);
          Eventtrainer_->Fill("B","MET", MET );
          ///	  Eventtrainer_->Fill("B","CvsL_1", CvsL_1 );
          Eventtrainer_->Fill("B","CvsB_1", CvsB_1 );
          //	  Eventtrainer_->Fill("B","CvsL_2", CvsL[1] );
          //	  Eventtrainer_->Fill("B","CvsB_2", CvsB[1] );
          Eventtrainer_->Fill("B","pt_electron_1", pt_electron_1 );
          //	  Eventtrainer_->Fill("B","pt_electron_2", pt_electron_2 );
          //	  Eventtrainer_->Fill("B","pt_electron_3", pt_electron_3 );
          Eventtrainer_->Fill("B","pt_muon_1", pt_muon_1 );
          //	  Eventtrainer_->Fill("B","pt_muon_2", pt_muon_2 );
          //	  Eventtrainer_->Fill("B","pt_muon_3", pt_muon_3 );
          Eventtrainer_->Fill("B","pt_jet_1", pt_jet_1 );
          //	  Eventtrainer_->Fill("B","pt_jet_2", pt_jet_2 );
          Eventtrainer_->Fill("B","topMass", topMass );
          //	  Eventtrainer_->Fill("B","nCSVL", nCSVL);
          Eventtrainer_->Fill("B","nCSVM", nCSVM);
          //	  Eventtrainer_->Fill("B","nCSVT", nCSVT);
          Eventtrainer_->Fill("B","bdis_1", bDisc[0]);
          Eventtrainer_->Fill("B","bdis_2", bDisc[1]);
          Eventtrainer_->Fill("B","FCNCtopmass",FCNCtopmass);
          Eventtrainer_->Fill("B","Pt_cjet",Pt_cjet);
          Eventtrainer_->Fill("B","deltaPhiSMFCNCtop",deltaPhiSMFCNCtop);
          Eventtrainer_->Fill("B","deltaPhiWlepb",deltaPhiWlepb);
//          Eventtrainer_->Fill("B","deltaPhiWlepc",deltaPhiWlepc);
          Eventtrainer_->Fill("B","deltaPhiZc",deltaPhiZc);
          Eventtrainer_->Fill("B","deltaPhiZb",deltaPhiZb);
          Eventtrainer_->Fill("B","deltaRSMFCNCtop",deltaRSMFCNCtop);
          Eventtrainer_->Fill("B","deltaRWlepb",deltaRWlepb);
          Eventtrainer_->Fill("B","deltaRWlepc",deltaRWlepc);
          Eventtrainer_->Fill("B","deltaRZc",deltaRZc);
          Eventtrainer_->Fill("B","deltaRZb",deltaRZb);
          Eventtrainer_->Fill("B","MassWlepB",MassWlepB);
          
        }
      }//for-loop events
      
    }//If-statement doTraining
    else if(!doTraining) //not training but computing
    {
      double BDTscore;
      TBranch *newb =  ttree[(dataSetName).c_str()]->Branch("BDTscore",&BDTscore,"BDTscore/D");;
      for (int j = 0; j<nEntries; j++)
      {
        ScaleFactor = 1.; // event scale factor
        ttree[(dataSetName).c_str()]->GetEntry(j);
        if(doScaling){
          if(!PassedMET){continue;}
          if(debug) cout << "MET SF " << endl;
          //electron
          for(unsigned int iEl = 0; iEl < NumberOfElectrons ; iEl ++)
          {
            //cout << "entry " << j << " electron " << iEl << " SF " <<  electronSF[iEl]<< endl;
            ScaleFactor *= electronSF[iEl];
          }
          if(debug) cout << "electron SF " << endl;
          //muon
          for(unsigned int iMu = 0; iMu < NumberOfMuons ; iMu ++)
          {
            ScaleFactor *= muonID[iMu]*muonIso[iMu];
            
          }
          if(debug) cout << "muon SF " << endl;
          //PU
          ScaleFactor *= puSF;
          if(debug) cout << "PU SF " << endl;
          //btag
          ScaleFactor *= BSF;
          if(debug) cout << "btag SF " << endl;
          //AMC
          ScaleFactor *= nloSF * nloW;
          if(debug) cout << "amc SF " << endl;
          
          if(dataSetName.find("NP")!=string::npos){ ScaleFactor = 1.;}
        }
        NormFactor = DataLumi / datasets[d]->EquivalentLumi(); //data->NormFactor*Lumi
        ScaleFactor *= NormFactor;
        if(debug) cout << "normalisation SF " << endl;
        
        if (Eventcomputer_ == 0) cout <<"null computer...." <<endl;

        
        
        Eventcomputer_->FillVar("Weight", ScaleFactor);
        
        Eventcomputer_->FillVar("pt_electron_1", pt_electron_1 );

        Eventcomputer_->FillVar("pt_muon_1", pt_muon_1 );
        Eventcomputer_->FillVar("topMass", topMass );
        //	  Eventcomputer_->FillVar("nCSVL", nCSVL);
        Eventcomputer_->FillVar("nCSVM", nCSVM);
        //	  Eventcomputer_->FillVar("nCSVT", nCSVT);
        //          Eventcomputer_->FillVar("bdis_1", bDisc[0]);
        Eventcomputer_->FillVar("bdis_2", bDisc[1]);
        Eventcomputer_->FillVar("FCNCtopmass",FCNCtopmass);
        //          Eventcomputer_->FillVar("Pt_cjet",Pt_cjet);
        //          Eventcomputer_->FillVar("deltaPhiSMFCNCtop",deltaPhiSMFCNCtop);
        //         Eventcomputer_->FillVar("deltaPhiWlepb",deltaPhiWlepb);
        //          Eventcomputer_->FillVar("deltaPhiWlepc",deltaPhiWlepc);
        //          Eventcomputer_->FillVar("deltaPhiZc",deltaPhiZc);
        Eventcomputer_->FillVar("deltaPhiZb",deltaPhiZb);
        //          Eventcomputer_->FillVar("deltaRSMFCNCtop",deltaRSMFCNCtop);
        Eventcomputer_->FillVar("deltaRWlepb",deltaRWlepb);
        //          Eventcomputer_->FillVar("deltaRWlepc",deltaRWlepc);
        //          Eventcomputer_->FillVar("deltaRZc",deltaRZc);
        Eventcomputer_->FillVar("deltaRZb",deltaRZb);
        Eventcomputer_->FillVar("MassWlepB",MassWlepB);

        
        
        std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
        
        for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
        {
          if(fabs(it->second > 1)) cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
          BDTscore = it->second;
        }
        //   cout << "BDT score " << BDTscore << endl;
        newb->Fill();
        //        if(isData) MSPlot[MVAmethod.c_str()]->Fill(BDTscore, datasets[d], true, 1.);
        //        else  MSPlot[MVAmethod.c_str()]->Fill(BDTscore, datasets[d], true, ScaleFactor);
      }
      
    }
    ttree[(dataSetName).c_str()]->Write(); 
    
  }//for-loop datasets
  std::string  mycutS = ""; 
  std::string mycutB = ""; 
  if(doTraining) Eventtrainer_->TrainMVA("Random",mycutS,0,0,mycutB,0,0,"test",false);
  
  delete Eventtrainer_;
  delete Eventcomputer_;
}









string ConvertIntToString(int Number, bool pad)
{
  ostringstream convert;
  convert.clear();
  if ( pad && Number < 10 ) { convert << std::setw(2) << std::setfill('0');}
  convert << Number;
  return convert.str();
};


string MakeTimeStamp()
{
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, true);
  string month_str = ConvertIntToString(month, true);
  string day_str = ConvertIntToString(day, true);
  string hour_str = ConvertIntToString(hour, true);
  string min_str = ConvertIntToString(min, true);
  //string sec_str = ConvertIntToString(sec, true);
  
  string date_str = year_str + month_str + day_str; //+ "_" + hour_str + min_str;
  return date_str;
};
