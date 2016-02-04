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



using namespace std;
using namespace TopTree;

// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,TTree*> globalttree;
map<string,MultiSamplePlot*> MSPlot;


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string pathPNG);
void MSPCreator (string pathPNG);

void TH2FPlotter (int nBinsX,float lowX, float highX, string sVarofinterestX );


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



// CONFIGURATION
Bool_t debug = false;
bool mumumu  = false;
string channelpostfix = "";
double DataLumi = -1; 

//applying all appropriate scale factors for individual objects if the bool is set to true
Bool_t applyElectronSF, applyMuonSF , applyPUSF, applyGlobalSF;




int main(int argc, char* argv[])
{
  if (debug){
    cout << "argc = " << argc << endl; 
    for(int i = 0; i < argc; i++)
      {
	cout << "argv[" << i << "] = " << argv[i] << endl; 
      }
  }
  if(argc < 5) cout << " ERROR: 5 arguments expected" << endl; 


  //Placing arguments in properly typed variables
  const string channel = argv[1];
  debug = false; 
  applyElectronSF = false;
  applyMuonSF = false;
  applyPUSF = false;
  applyGlobalSF = false;
  debug = strtol(argv[2],NULL,10); 
  applyElectronSF = strtol(argv[3],NULL,10);
  applyMuonSF = strtol(argv[4],NULL,10);
  applyPUSF = strtol(argv[5],NULL,10);
  applyGlobalSF = strtol(argv[6],NULL,10);
  string xmlFileName;
  string CraneenPath; 
  CraneenPath = "NtupleMakerOutput/MergedTuples/";
  if(channel=="MuMuMu")
    {
      cout << " --> Using the TriMuon channel..." << endl;
      channelpostfix = "_mumumu";
      xmlFileName = "config/Run2TriLepton_samples_analyzer_mumumu.xml";
      mumumu = true; 
      DataLumi  = 10000;
      CraneenPath += "mumumu/";  
    }
  else
    {
      cerr << "The channel '" << channel << "' is not in the list of authorised channels !!" << endl;
      exit(1);
    }
    string dateString = MakeTimeStamp(); 
    CraneenPath += dateString + "/"; 
    string pathPNG = "myOutput";
    mkdir(pathPNG.c_str(),0777); 
    pathPNG += "/" + dateString + "/"; 
    mkdir(pathPNG.c_str(),0777);
    pathPNG += "MSPlots"+channelpostfix+"/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;         

    cout << "xmlFileName is " << xmlFileName << endl;
    cout << "NtupleFolder is " << CraneenPath << endl; 
    if(debug) cout << "debugging on" << endl; 
    if(applyElectronSF) cout << "Electron SF on " << endl;
    if(applyMuonSF) cout << "Muon SF on " << endl; 
    if(applyPUSF) cout << "PU SF on" <<endl; 
    if(applyGlobalSF) cout << "applying SF " << endl; 


    // calling datasetPlotter to create MSPplots
          
      // event plots
    DatasetPlotter(70, -0.5, 69.5, "npu", xmlFileName,CraneenPath,pathPNG);
      //    DatasetPlotter(70, -0.5, 69.5, "nvtx", xmlFileName,CraneenPath);
          
     
      /*
      // electron plots
      //    DatasetPlotter(11, -0.5, 10.5, "nElectrons_mumu", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0, 1000, "pt_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(50, -3.15, 3.15, "eta_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(30, -3.15, 3.15, "phi_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      DatasetPlotter(100, 0.0, 0.2, "pfIso_electron_mumu[nElectrons_mumu]", xmlFileName,CraneenPath);
      */
      
 
      // muon plots
      //    DatasetPlotter(11, -0.5, 10.5, "nMuons", xmlFileName,CraneenPath);
//      DataddsetPlotter(100, 0, 1000, "pt_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      //	DatasetPlotter(20, 10, 50, "pt_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
  //    DatasetPlotter(50, -3.15, 3.15, "eta_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
   //   DatasetPlotter(30, -3.15, 3.15, "phi_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
      //    DatasetPlotter(100, -0.1, 0.1, "d0_muon[nMuons]", xmlFileName,CraneenPath);
     // DatasetPlotter(100, -10, 10, "vz_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
     // DatasetPlotter(100, -0.015, 0.015, "d0BeamSpot_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);
     // DatasetPlotter(100, 0.0, 0.2, "pfIso_muon_mumu[nMuons_mumu]", xmlFileName,CraneenPath);



  // eo selecting the right plots and/or variables depending on the final state  


  // calling the function that writtes all the MSPlots in a root file
  MSPCreator (pathPNG);

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath, string PathPNG)
{
  cout<<""<<endl;
  cout<<"RUNNING NOMINAL DATASETS"<<endl;
  cout<<""<<endl;

  const char *xmlfile = xmlNom.c_str();
  if(debug) cout << "used config file: " << xmlfile << endl;
  
  
  
  ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					
  if (debug) cout << "will start loading from xml file ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  if (debug) cout << "finished loading from xml file ..." << endl;
  
  
   //***************************************************CREATING PLOTS****************************************************
  //  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  //  outfile->cd();
  string plotname = sVarofinterest;   ///// Non Jet Split plot
  // make for loop here!!!
  MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str()); 

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath ,fileDatapath, slumi;
  
  int nEntries;
  float ScaleFactor, NormFactor;
  int varofInterest;
  double varofInterest_double [20];


  
  vector<string> v;
  // to avoid modifying original string
  // first duplicate the original string and return a char pointer then free the memory
  
  char delim[] = " []";
  char * dup = strdup(sVarofinterest.c_str());
  char * token = strtok(dup, delim);
  while(token != NULL){
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);

//  if (debug) cout << v[0] << "  " << v[1] << endl;
  
   for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
   {
     dataSetName = datasets[d]->Name();
     cout<<"Dataset:  :"<<dataSetName<<endl;
     
     // get the tree corresponding to the final state of interest
     string stree = "tree";
     string sglobaltree = "globaltree"; 
     string sbaselinetree = "baselinetree"; 
     filepath = TreePath + dataSetName + ".root"; 
     string dataN = "DataRunD";
     fileDatapath = TreePath + dataN + ".root"; 
		  
     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
     string TTreename = sbaselinetree;	

     
     ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
     nEntries = (int)ttree[dataSetName.c_str()]->GetEntries();
     cout<<"                 nEntries: "<<nEntries<<endl;     
     


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
      // eo logic to set the right branch address depending on the string given as argument of the datasetplotter

      bool isData= false;
      if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
      if(debug) cout << "isData? " << isData << endl; 
      if(isData) applyGlobalSF = false;
      ///////////////////////////////////
      // determine event scalefactor ///
      //////////////////////////////////
      if(applyGlobalSF) cout << "                 Applying scale factors " << endl; 
      
      // get the SF from the corresponding branch
      Double_t puSF = 1. ;
      ttree[dataSetName.c_str()]->SetBranchAddress("puSF",&puSF); 
      
      Double_t electronSF[10]; 
      ttree[dataSetName.c_str()]->SetBranchAddress("ElectronSF",&electronSF); 
      
      Double_t muonID[10]; 
      ttree[dataSetName.c_str()]->SetBranchAddress("MuonIDSF", &muonID); 
      
      Double_t muonIso[10]; 
      ttree[dataSetName.c_str()]->SetBranchAddress("MuonIsoSF", &muonIso);
      
      Double_t muonTrigv2[10]; 
      ttree[dataSetName.c_str()]->SetBranchAddress("MuonTrigSFv2", &muonTrigv2);  
      
      Double_t muonTrigv3[10]; 
      ttree[dataSetName.c_str()]->SetBranchAddress("MuonTrigSFv3", &muonTrigv3); 
      
      Int_t nEl; 
      ttree[dataSetName.c_str()]->SetBranchAddress("nElectrons",&nEl);
      Int_t nMu;
      ttree[dataSetName.c_str()]->SetBranchAddress("nMuons",&nMu);
      

      
      if(debug) cout << "done setting SF addresses " << endl; 
      
      // -----------
      // eo of event SF
      
      double globalScaleFactor; 
      for (int j = 0; j<nEntries; j++)
      {
	  ttree[(dataSetName).c_str()]->GetEntry(j);
          globalScaleFactor = 1.; 
	  
	  if(applyGlobalSF) // sf on and not data
	  {
	      // Electron scale factors
	      if(applyElectronSF)
	      {
	  	for(unsigned int i = 0; i < nEl ; i ++)
		{
		   if(debug) cout << "lepton sf at index " << i << " is " << electronSF[i] << endl;
		   globalScaleFactor *= electronSF[i];
		   if(debug) cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	  	}
	      }
	      if(applyMuonSF)
	      {
	  	for(unsigned int i = 0; i < nMu ; i ++)
		{
		   if(debug) cout << "Muon ID sf at index " << i << " is " << muonID[i] << endl;
		   if(debug) cout << "Muon Iso sf at index " << i << " is " << muonIso[i] << endl;
		   if(debug) cout << "Muon trig v2 sf at index " << i << " is " << muonTrigv2[i] << endl;
		   if(debug) cout << "Muon trig v3 sf at index " << i << " is " << muonTrigv3[i] << endl;
		   
		   globalScaleFactor *= muonID[i] *  muonIso[i]  ;
		   if(debug) cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	  	}
	      }
	      if(applyPUSF)
	      {
	           globalScaleFactor *= puSF;
	           if (debug){
	        	cout << "puSF is " << puSF << endl;
		        cout << "the globalScaleFactor is " << globalScaleFactor << endl;
	           }
	      
	      }
	  
	  }
	  // ----------------
	  // eo event SF
	  
	  // make MS plot for single value
	  if (v.size() == 1){
	    if (isData) 
	    {
		// for data, fill once per event, weighted with the event scale factor only ???? what??
		MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], false, 1); 	
	    }
	    else
	    {
		// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
		MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, globalScaleFactor*DataLumi);
	    }
	  }
	  // make MS plot for vector
	  if (v.size() == 2){

	    // bo of loop over the number of object per entry
	    for (int i_object =0 ; i_object < varofInterest ;i_object ++ )
	      {
		if (debug) cout << "varofInterest is " << varofInterest << endl;
		if (isData) 
		  {
		    // for data, fill once per event, weighted with the event scale factor 
		    MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], false,1); 	
		  }
		else
		  {
		    // for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
		    MSPlot[plotname.c_str()]->Fill(varofInterest_double[i_object], datasets[d], true, globalScaleFactor*DataLumi);

		  }
	      
	      }

	  }
	  
      } // nentries
      
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
   } // datasets
   
   
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
}; // datasetplotter


// function that writes all the MSPlots created in a root file
void MSPCreator (string pathPNG)
{  
 // cout << pathPNG.c_str() << endl; 
  
  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  outfile->cd();
  
  
  // Loop over all the MSPlots
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      string name = "MyMSP_" + it->first;
      MultiSamplePlot *temp = it->second;
      if (debug){
	cout << "Saving the MSP" << endl;
	cout << " and it->first is " << it->first << endl;
      }
      temp->Draw("MyMSP", 1, false, false, false, 10);
      
      if(!applyGlobalSF) name += "_noSF";
      if(!applyPUSF) name += "_noPUSF";
      if(!applyElectronSF) name += "_noElSF"; 
      if(!applyMuonSF) name+= "_noMuSF"; 

	    temp->Write(outfile, name, true,pathPNG.c_str() , "png");
     //      vector<string> temp_histo = it->GetTH1FNames();
      //      for (int i_hist=0; i_hist < temp_histo.size();i_hist++  ){
      //	cout << "hist is" << temp_histo[i_hist] << endl;
      //	cout << "integral is " << it->GetTH1F(temp_histo[i_hist].GetSum()) << endl;
      //      }
    }
  
  outfile->Write("kOverwrite");
}

