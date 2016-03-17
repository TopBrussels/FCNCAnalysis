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


float Luminosity = 2628.727204156; // pb-1 Muon  = 2628.727204156, Electron = 2094.087
TString slumi = "2628.727204156";
std::string channel = "_El";
std::string date = "_17_3_2016";
int maxNumbObjToPlot = 5;
Bool_t debug = false;
bool applyAMC;

// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string xmlNom, string TreePath, string NTupleName);
void MSPCreator ();



int main()
{
    int NumberOfBins = 3;	//fixed width nBins


    string xmlFileName;
    string CraneenPath;
    


    xmlFileName = "config/FullMcBkgdSamplesV8_TreeProcessor.xml";
    cout << "xmlFileName is " << xmlFileName << endl;



    // calling datasetPlotter to create MSPplots

    // electron plots
    /*DatasetPlotter(11, -0.5, 10.5, "I_nJets", xmlFileName,CraneenPath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVL", xmlFileName,CraneenPath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVM", xmlFileName,CraneenPath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVT", xmlFileName,CraneenPath,"ObjectVarsTree");
    //if(channel == "_Mu") DatasetPlotter(40, 0, 400, "pt_muon", xmlFileName,CraneenPath,"ObjectVarsTree");
    //else if(channel == "_El") DatasetPlotter(40, 0, 400, "pt_electron", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(40, 0., 500., "MTlepmet", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MLepTop_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW_GenMatch", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 5., "dR_lepJet_min", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "Mbb", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW", xmlFileName,CraneenPath,"AdvancedVarsTree");
    DatasetPlotter(35, -0.5, 34.5, "I_nvtx", xmlFileName,CraneenPath,"ObjectVarsTree");
    */DatasetPlotter(70, 0, 700, "pt_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(50, -3.15, 3.15, "eta_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(30, -3.15, 3.15, "phi_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(25, 0, 1, "bdisc_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsL_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[nJets]", xmlFileName,CraneenPath,"ObjectVarsTree");

    // calling the function that writtes all the MSPlots in a root file
	MSPCreator ();

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string xmlNom, string TreePath, string NTupleName)
{
  	cout<<""<<endl;
  	cout<<"RUNNING NOMINAL DATASETS"<<endl;
  	cout<<""<<endl;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
  	string pathPNG = "";
  	pathPNG += "MSPlots";
  	pathPNG += channel;
  	mkdir(pathPNG.c_str(),0777);
  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory
  
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


  	TString CraneenPath = "Merged/Ntuples" + channel + "/Ntuples" + date;

     if (v.size() == 2)//Meaning we have a variable of the form "var[n_obj]", which is an array of values for the variable 'var'
     {

                //If plotting a variable which consists of several values (e.g.: jet_pt contains the pt of all jets), make also plots for the individual values (e.g.: plot the pt for every jet separately). For now, only done for 5 first objects
                for(int iToPlot = 1; iToPlot <= maxNumbObjToPlot; iToPlot++)
                {
                          string conv_str;
                          ostringstream conv;   // stream used for the conversion
                          conv << (iToPlot);      // insert the textual representation of 'Number' in the characters in the stream
                          conv_str = "_"+conv.str(); // set 'Result' to the contents of the stream

                          MSPlot[(plotname+conv_str).c_str()] = new MultiSamplePlot(datasets, (plotname+conv_str).c_str(), nBins, plotLow, plotHigh, s_varofInterest.c_str()); 
                }                

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = NTupleName;	
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
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
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttre for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  //ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);

                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(applyAMC && isAMC && !isData)
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
				                    if(i_obj< maxNumbObjToPlot) MSPlot[(plotname+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, 1.);//Fill MSPlot for first 5 variables
			                    }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, ScaleFactor*Luminosity);
				                    if(i_obj<maxNumbObjToPlot) MSPlot[(plotname+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true,  ScaleFactor*Luminosity);//Fill MSPlot for first 5 variables
			                    }
			                }
			                
			                
		              }
		              
		          }//for-loop datasets
               

      }//end statement on variable-plotting consisting of array		   (v.size()==2)
     else if (v.size() == 1)
     {
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
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
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttre for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  //ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);

                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(applyAMC && isAMC && !isData)
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
  


};


// function that writes all the MSPlots created in a root file
void MSPCreator ()
{
  	Bool_t debug = false;

  	string pathPNG = "myOutput";
  	pathPNG += "_MSPlots";
  	pathPNG += channel;
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


