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

//includes for MVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"


using namespace std;
using namespace TopTree;
//using namespace KINFIT;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;


float Luminosity = 12900; // pb-1
bool ManualLuminosity = false;  //Put this to true if you have produced the data n-tuples with the wrong luminosity.
std::string channel = "_El";
std::string date = "_11_10_2016";
int maxNumbObjToPlot = 5;
Bool_t debug = false;
float workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Medium = 0.800;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Tight = 0.935;//working points updated to 2016 BTV-POG recommendations.

std::string xmlNom = "config/FullMcBkgdSamples_TreeProcessor.xml";
TString TreePath = "Merged/Ntuples" + channel + "/Ntuples" + date;


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string NTupleName, bool Inclusive, int Nbjets, int Njets, string units);
void MSPCreator ();



int main(int argc, char *argv[])
{

    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);

    bool doInclusive = false;
    if(baseline_bjets == 0 && baseline_jets == 0) doInclusive = true;
    

    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    if(!doInclusive) cout << "MACRO command line arguments, category: " << baseline_bjets << "b" << baseline_jets << "j" << endl;
    else cout << "MACRO command line arguments, Inclusive " << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();

    // calling datasetPlotter to create MSPplots

//    DatasetPlotter(11, -0.5, 10.5, "I_nJets", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets,"");
//   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVL", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
//   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVM", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets,"");
//   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVT", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
//    DatasetPlotter(40, 0, 400, "pt_muon", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
//    DatasetPlotter(40, 0, 400, "pt_electron", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
//    DatasetPlotter(35, -0.5, 34.5, "I_nvtx", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets, "");
    DatasetPlotter(70, 0, 700, "pt_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets, "GeV");
/*    DatasetPlotter(50, -3.15, 3.15, "eta_Jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(30, -3.15, 3.15, "phi_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(21, -10.5, 10.5, "charge_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(25, 0, 1, "bdisc_Jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_pdgID[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_motherpdgID[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_grannypdgID[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(25,-1, 1, "cdiscCvsL_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
    DatasetPlotter(50,-1, 1, "incl_charge_jet[I_nJets]", "ObjectVarsTree", doInclusive, baseline_bjets, baseline_jets);
*/


	MSPCreator ();

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

}







void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string NTupleName, bool Inclusive, int Nbjets, int Njets, string units)
{
  	cout<<""<<endl;
  	cout<<"RUNNING NOMINAL DATASETS"<<endl;
  	cout<<""<<endl;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
  	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  
    //***************************************************CREATING PLOT****************************************************
  	string plotname = s_varofInterest;
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, s_varofInterest.c_str(),"Events", "", units); 
    MSPlot["Njets"] = new MultiSamplePlot(datasets, "Njets", 11, -0.5, 10.5, "Number of jets","Events", "", units); 
    MSPlot["NBjets"] = new MultiSamplePlot(datasets, "NBjets", 10, -0.5, 10.5, "Number of CSVv2 M","Events", "", units); 

  
  	//***********************************************OPEN FILES & GET NTUPLES**********************************************
  	string dataSetName, filepath;
  	int nEntries;
  	float ScaleFactor, NormFactor;
  	int varofInterest;
  	Double_t d_varofInterest;
 	  int n_object = 0;
  	double v_d_varofInterest_double [20];
    int NumberOfJets, NumberOfBjets;
 
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


     if (v.size() == 2)//Meaning we have a variable of the form "var[n_obj]", which is an array of values for the variable 'var'
     {

                //If plotting a variable which consists of several values (e.g.: jet_pt contains the pt of all jets), make also plots for the individual values (e.g.: plot the pt for every jet separately). For now, only done for 5 first objects
                for(int iToPlot = 1; iToPlot <= maxNumbObjToPlot; iToPlot++)
                {
                          string conv_str;
                          ostringstream conv;   // stream used for the conversion
                          conv << (iToPlot);      // insert the textual representation of 'Number' in the characters in the stream
                          conv_str = "_"+conv.str(); // set 'Result' to the contents of the stream

                          MSPlot[(v[0]+conv_str).c_str()] = new MultiSamplePlot(datasets, (v[0]+conv_str).c_str(), nBins, plotLow, plotHigh, (v[0]+conv_str).c_str(),"Events", "", units); 
                }                

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = NTupleName;	
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;
		                
		                
                   // bo logic to set the right branch address depending on the string given as argument of the datasetplotter
	                 ttree[dataSetName.c_str()]->SetBranchAddress(v[0].c_str(),v_d_varofInterest_double); //v[0] is the string of the variable you want to plot. This variable should be an array of values, according to the number of objects
	                 ttree[dataSetName.c_str()]->SetBranchAddress(v[1].c_str(),&n_object); // v[1] is the string of the variable between [] in the string. This should correspond to the number of objects


		              bool isData= false;
		              bool isAMC = false;
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
		              string TTreename_NtupleInfo = "NtupleInfoTree";	
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttree for each dataset
		              ttree[(dataSetName + "NtupleInfoTree").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_NtupleInfo.c_str()); //get ttree for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->GetEntry(0);
                  Luminosity = luminosity_;
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_CSVv2M_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);
                  
                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&NumberOfBjets);
                  

                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(isAMC && !isData)
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
		              
			                ttree[dataSetName.c_str()]->GetEntry(j);
		                  if(!Inclusive)
		                  {
		                      if(Njets == 3 && NumberOfJets != Njets && NumberOfBjets != Nbjets) continue;
		                      if(Njets == 4 && NumberOfJets < Njets && NumberOfBjets != Nbjets) continue;
		                  }
		              
                  		ScaleFactor = 1.; // event scale factor
//			                ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                if(ScaleFactor < 0) ScaleFactor = 0;
//			                ttree[dataSetName.c_str()]->GetEntry(j);
			                
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
				                    MSPlot["Njets"]->Fill(NumberOfJets, datasets[d], true, 1.);
				                    MSPlot["NBjets"]->Fill(NumberOfBjets, datasets[d], true, 1.);
				                    if(i_obj< maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, 1.);//Fill MSPlot for first 5 variables
			                    }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot["Njets"]->Fill(NumberOfJets, datasets[d], true, ScaleFactor*Luminosity);
				                    MSPlot["NBjets"]->Fill(NumberOfBjets, datasets[d], true, ScaleFactor*Luminosity);
				                    MSPlot[plotname.c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, ScaleFactor*Luminosity);
				                    if(i_obj<maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true,  ScaleFactor*Luminosity);//Fill MSPlot for first 5 variables
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
		              filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = NTupleName;	
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttre for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;


                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&NumberOfBjets);

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
		              bool isAMC = false;
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
		              string TTreename_NtupleInfo = "NtupleInfoTree";	
		              ttree[(dataSetName + "weights").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_Weights.c_str()); //get ttree for each dataset
		              ttree[(dataSetName + "NtupleInfoTree").c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_NtupleInfo.c_str()); //get ttree for each dataset
		
                  Double_t lumiweight, LeptonSF, bTagSF, luminosity_;
                  Double_t  nloweight;
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->SetBranchAddress("Luminosity_",&luminosity_);
                  ttree[(dataSetName + "NtupleInfoTree").c_str()]->GetEntry(0);
                  Luminosity = luminosity_;
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("puSF",&lumiweight);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("fleptonSF",&LeptonSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_CSVv2M_mujets_central",&bTagSF);
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("nloWeight",& nloweight);


                    double nloSF = 1.;
                    int nPos = 0; 
                    int nNeg = 0;
                    int Ev = 0; 
                    if(isAMC && !isData)
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
			                ttree[dataSetName.c_str()]->GetEntry(j);
		                  if(!Inclusive)
		                  {
		                      if(Njets == 3 && NumberOfJets != Njets && NumberOfBjets != Nbjets) continue;
		                      else if(Njets == 4 && NumberOfJets < Njets && NumberOfBjets != Nbjets) continue;
		                  }
		              
                  		ScaleFactor = 1.; // event scale factor
//			                ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                if(ScaleFactor < 0) ScaleFactor = 0;
//			                ttree[dataSetName.c_str()]->GetEntry(j);

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
				                    MSPlot["Njets"]->Fill(NumberOfJets, datasets[d], true, 1.);
				                    MSPlot["NBjets"]->Fill(NumberOfBjets, datasets[d], true, 1.);
				                  }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(d_varofInterest, datasets[d], true, ScaleFactor*Luminosity);
				                    MSPlot["Njets"]->Fill(NumberOfJets, datasets[d], true, ScaleFactor*Luminosity);
				                    MSPlot["NBjets"]->Fill(NumberOfBjets, datasets[d], true, ScaleFactor*Luminosity);
			                    }
			                }
			                		if(debug) cout << "Event " << j << endl;

			            }
			                
			                
		           }//for-loop datasets
     }//end of statement if variable is not an array of values
      
      




//	treeLoader.UnLoadDataset();
  	// clearing vector
  	v.clear();
  	if (debug){
    	cout << "after cleaning" << endl ;
    	cout << "v.size() is " << v.size() << endl;
  	}
  
cout << "MSPlot size: " << MSPlot.size() << endl;      


};






// function that writes all the MSPlots created in a root file
void MSPCreator ()
{

  	string pathPNG = "MSPlots/";
  	mkdir(pathPNG.c_str(),0777);
  	pathPNG += "MSPlots";
  	pathPNG += channel;
  	mkdir(pathPNG.c_str(),0777);
  	pathPNG += "/";
  	pathPNG += date;
  	mkdir(pathPNG.c_str(),0777);
  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  	TFile *outfile = new TFile((pathPNG+"/Output.root").c_str(),"recreate");
  	outfile->cd();


  	// Loop over all the MSPlots
  	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
      	string name = it->first;
      	MultiSamplePlot *temp = it->second;
      	//if (debug){
			cout << "Saving the MSP" << endl;
			cout << " and it->first is " << it->first << endl;
			cout << " Luminosity is " << Luminosity << endl;
      	//}
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
        cout << "Drawing MSP: " << it->second << endl;
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

