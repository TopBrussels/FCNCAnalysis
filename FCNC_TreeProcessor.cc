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

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;


float Luminosity = 569.017464; // pb-1
TString slumi = "569.017464";


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath);
void MSPCreator ();



int main()
{
    int NumberOfBins = 3;	//fixed width nBins


    string xmlFileName;
    string CraneenPath;


	xmlFileName = "config/FullMcBkgdSamplesV7TreeProc_TreeProcessor.xml";
    cout << "xmlFileName is " << xmlFileName << endl;



    // calling datasetPlotter to create MSPplots

    // electron plots
    DatasetPlotter(11, -0.5, 10.5, "nb_jets", xmlFileName,CraneenPath);
    DatasetPlotter(11, -0.5, 10.5, "nb_bjets", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0, 400, "leptonpt", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0, 400, "jet1_Pt", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0, 400, "jet2_Pt", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0, 400, "jet3_Pt", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0, 400, "MissingEt", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0., 1, "bdisc1", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0., 1, "bdisc2", xmlFileName,CraneenPath);
    DatasetPlotter(40, 0., 1, "bdisc3", xmlFileName,CraneenPath);

    // muon plots

    // calling the function that writtes all the MSPlots in a root file
	MSPCreator ();

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string sVarofinterest, string xmlNom, string TreePath)
{
  	Bool_t debug = false;
  	cout<<""<<endl;
  	cout<<"RUNNING NOMINAL DATASETS"<<endl;
  	cout<<""<<endl;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
  	string pathPNG = "myOutput";
  	pathPNG += "_MSPlots/";
  	mkdir(pathPNG.c_str(),0777);
  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  
    //***************************************************CREATING PLOTS****************************************************
  	//  TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  	//  outfile->cd();
  	string plotname = sVarofinterest;   ///// Non Jet Split plot
  	// make for loop here!!!
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str()); 

  
  	//***********************************************OPEN FILES & GET NTUPLES**********************************************
  	string dataSetName, filepath;
  	int nEntries;
  	float ScaleFactor, NormFactor;
  	float varofInterest;
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


  	TString CraneenPath = "/user/kderoove/FCNC/TopTreeFramework_Run2/CMSSW_7_4_14/src/TopBrussels/FCNCAnalysis/Trees_SelectionOutput_Mu/Trees_SelectionOutput_21_10_2015";

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
		dataSetName = datasets[d]->Name();
		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName +"_Mu"+ ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
	

		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "tree";	
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
		  
		// bo logic to set the right branch address depending on the string given as argument of the datasetplotter
		if (v.size() == 2){
			ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&varofInterest); 
			ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),varofInterest_double);
			if(debug) cout << "setting right branch adress" << endl;
		}

		else if (v.size() == 1){
			if (debug)	cout << "v.size is to 1" << " and v[0] is " << v[0] << endl ;
			ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),&varofInterest);//&varofInterest
		}
		else {
			cout << "Vector of string does not have the good size!!!" << endl;
		}
		// eo logic to set the right branch address depending on the string given as argument of the datasetplotter
		  
		  

		bool isData= false;
		if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) isData =true;
		  
		ScaleFactor = 1.; // event scale factor
		  
		//      histo1D[dataSetName.c_str()] = new TH1F((dataSetName+"_"+v[0]).c_str(),(dataSetName+"_"+v[0]).c_str(), nBins, plotLow, plotHigh);

		// bo of loop through entries and fill plots
		for (int j = 0; j<nEntries; j++)
		{
			ttree[dataSetName.c_str()]->GetEntry(j);

			if (debug) cout << "varofInterest is " << varofInterest << endl;
			if (isData) 
			{// for data, fill once per event, weighted with the event scale factor
				MSPlot[plotname.c_str()]->Fill(varofInterest/*_double[i_object]*/, datasets[d], true, ScaleFactor);
			}
			else
			{// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				MSPlot[plotname.c_str()]->Fill(varofInterest/*_double[i_object]*/, datasets[d], true, ScaleFactor*Luminosity);
			}
		}
		  
	 
		TCanvas *canv = new TCanvas(("canv_"+v[0]+dataSetName).c_str(),("canv_"+v[0]+dataSetName).c_str());
		  
		  
		//      histo1D[dataSetName.c_str()]->Draw();
		string writename = "";
		if(isData)
		{
			writename = "data_nominal";
		}
		else
		{
			writename = dataSetName +"__nominal";
		}
	}// end loop datasets


	treeLoader.UnLoadDataset();
  
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
  


};


// function that writes all the MSPlots created in a root file
void MSPCreator ()
{
  	Bool_t debug = false;

  	string pathPNG = "myOutput";
  	pathPNG += "_MSPlots/";
  	mkdir(pathPNG.c_str(),0777);
  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
  	TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  	outfile->cd();
  
  
  	// Loop over all the MSPlots
  	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      	string name = it->first;
      	MultiSamplePlot *temp = it->second;
      	if (debug){
			cout << "Saving the MSP" << endl;
			cout << " and it->first is " << it->first << endl;
      	}
      	temp->setDataLumi(Luminosity);
		temp->Draw("MyMSP_"+it->first, 0, false, false, false, 1);
      	temp->Write(outfile, it->first, false,"myOutput_MSPlots" , "png");
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


