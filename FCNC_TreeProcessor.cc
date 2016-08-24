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

//includes for MVA
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit_final/kinfit.h"


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


float Luminosity = 0.; // pb-1 Muon  = 2628.727204156, Electron = 2094.087
std::string channel = "_El";
std::string date = "_26_7_2016";
int maxNumbObjToPlot = 5;
Bool_t debug = false;
bool train_mva = true;
float workingpointvalue_Loose = 0.460;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Medium = 0.800;//working points updated to 2016 BTV-POG recommendations.
float workingpointvalue_Tight = 0.935;//working points updated to 2016 BTV-POG recommendations.


// functions prototype
std::string intToStr (int number);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string xmlNom, TString CraneenPath, string NTupleName);
void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath);
void MVA_JetCombTraining(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod);
void MVA_JetCombComputer(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod);
void MCAnalysis(std::string xmlNom, TString CraneenPath, string KinFitMethod);
void MSPCreator ();



int main()
{


    string xmlFileName;
  	TString TreePath = "Merged/Ntuples" + channel + "/Ntuples" + date;
    


    xmlFileName = "config/FullMcBkgdSamplesV8_TreeProcessor.xml";
    cout << "xmlFileName is " << xmlFileName << endl;



    // calling datasetPlotter to create MSPplots

/*    DatasetPlotter(11, -0.5, 10.5, "I_nJets", xmlFileName,TreePath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVL", xmlFileName,TreePath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVM", xmlFileName,TreePath,"ObjectVarsTree");
   	DatasetPlotter(11, -0.5, 10.5, "I_nJets_CSVT", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(40, 0, 400, "pt_muon", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(40, 0, 400, "pt_electron", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(35, -0.5, 34.5, "I_nvtx", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(70, 0, 700, "pt_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(50, -3.15, 3.15, "eta_Jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(30, -3.15, 3.15, "phi_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(21, -10.5, 10.5, "charge_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(25, 0, 1, "bdisc_Jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_pdgID[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_motherpdgID[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(59,-29.5, 29.5, "jet_matchedMC_grannypdgID[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(25,-1, 1, "cdiscCvsB_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(50,-1, 1, "incl_charge_jet[I_nJets]", xmlFileName,TreePath,"ObjectVarsTree");
    DatasetPlotter(40, 0., 500., "MTlepmet", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MLepTop_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW_GenMatch", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 5., "dR_lepJet_min", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MHadTop", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaLepTop", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaHadTop", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "MassW", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, 0., 500., "Mbb", xmlFileName,TreePath,"AdvancedVarsTree");
    DatasetPlotter(40, -5., 5., "EtaW", xmlFileName,TreePath,"AdvancedVarsTree");
*/


  //MVAanalysis(train_mva, "BDT", 3, "NP_overlay_ST_tHToBB_1L_Kappa_hct", xmlFileName, TreePath);
  // calling the function that writtes all the MSPlots in a root file
//    MVA_JetCombTraining( "BDT", 3, "NP_overlay_ST_tHToBB_1L_Kappa_hct", xmlFileName, TreePath, "STSignal");
//    MVA_JetCombComputer( "BDT", 3, "NP_overlay_TTtoTHToBB-1L-Kappa-hct", xmlFileName, TreePath, "TTSignal");
    MCAnalysis(xmlFileName, TreePath, "ALL");
  //  MVA_JetCombTraining( "BDT", 10, "NP_overlay_ST_tHToBB_1L_Kappa_hct", xmlFileName, TreePath);
	MSPCreator ();

}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string s_varofInterest, string xmlNom, TString CraneenPath, string NTupleName)
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


     if (v.size() == 2)//Meaning we have a variable of the form "var[n_obj]", which is an array of values for the variable 'var'
     {

                //If plotting a variable which consists of several values (e.g.: jet_pt contains the pt of all jets), make also plots for the individual values (e.g.: plot the pt for every jet separately). For now, only done for 5 first objects
                for(int iToPlot = 1; iToPlot <= maxNumbObjToPlot; iToPlot++)
                {
                          string conv_str;
                          ostringstream conv;   // stream used for the conversion
                          conv << (iToPlot);      // insert the textual representation of 'Number' in the characters in the stream
                          conv_str = "_"+conv.str(); // set 'Result' to the contents of the stream

                          MSPlot[(v[0]+conv_str).c_str()] = new MultiSamplePlot(datasets, (v[0]+conv_str).c_str(), nBins, plotLow, plotHigh, (v[0]+conv_str).c_str()); 
                }                

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
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
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
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
		
		
		              Int_t ThreeJets = 0;
		              Int_t MoreThanThreeJets = 0;
		              //////////////////////////////////////////////////////////
		              // Making MS plots
		              //////////////////////////////////////////////////////////
		              for (int j = 0; j<nEntries; j++)
		              {
                  		ScaleFactor = 1.; // event scale factor
			                ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                if(ScaleFactor < 0) ScaleFactor = 0;
			                ttree[dataSetName.c_str()]->GetEntry(j);
			                
			                if(n_object == 3) ThreeJets++;
			                else MoreThanThreeJets++; 

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
				                    if(i_obj< maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, 1.);//Fill MSPlot for first 5 variables
			                    }
			                    else
			                    {// for MC, fill once per event and multiply by the event scale factor. Then reweigt by Lumi/Eqlumi where Eqlumi is gotten from the xml file
				                    MSPlot[plotname.c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true, ScaleFactor*Luminosity);
				                    if(i_obj<maxNumbObjToPlot) MSPlot[(v[0]+conversion_str).c_str()]->Fill(v_d_varofInterest_double[i_obj], datasets[d], true,  ScaleFactor*Luminosity);//Fill MSPlot for first 5 variables
			                    }
			                }
			                
			                
		              }
		              
		              
		              cout << dataSetName << ": Number of Events with exactly 3 jets: " << ThreeJets << endl;
		              cout << dataSetName << ": Number of Events with more than 3 jets: " << MoreThanThreeJets << endl;
		              
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
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
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
			                		if(debug) cout << "Event " << j << endl;

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
  
cout << "MSPlot size: " << MSPlot.size() << endl;      


};



// function that converts an int into a string
std::string intToStr (int number)
{
  	std::ostringstream buff;
  	buff<<number;
  	return buff.str();
}



void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath)
{

  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);
  string pathMVA_ = pathMVA+"TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);


  MVAComputer* Eventcomputer_ =0;   
  MVATrainer* Eventtrainer_ = 0;
  if(doTraining) Eventtrainer_ = new MVATrainer(MVAmethod,"TrainedEventMVA"+channel, pathMVA_+"TrainedEventMVA"+channel+"_"+SignalName+".root");
  vector<std::string> MVAvars;
  
  MVAvars.push_back("InclCharge_bjet1");
  MVAvars.push_back("InclCharge_bjet2");
  MVAvars.push_back("InclCharge_bjet3");
  MVAvars.push_back("SummedCharge_bjet1");
  MVAvars.push_back("SummedCharge_bjet2");
  MVAvars.push_back("SummedCharge_bjet3");
  MVAvars.push_back("CvsL_bjet1");
  MVAvars.push_back("CvsL_bjet2");
  MVAvars.push_back("CvsL_bjet3");
  MVAvars.push_back("CvsB_bjet1");
  MVAvars.push_back("CvsB_bjet2");
  MVAvars.push_back("CvsB_bjet3");
  
  for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
  {
      if(doTraining) Eventtrainer_->bookInputVar(MVAvars[N_var]);
  }
  
  if(!doTraining) Eventcomputer_ = new MVAComputer(MVAmethod,pathMVA_+"TrainedEventMVA"+channel+"_"+SignalName+".root", "TrainedEventMVA"+channel,MVAvars, "test");



  	cout<<""<<endl;
  	cout<<"RUNNING MVA"<<endl;
  	cout<<""<<endl;

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
  	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  
    //***************************************************CREATING PLOT****************************************************
    MultiSamplePlot* MSP_MVA = new MultiSamplePlot(datasets, MVAmethod.c_str() , 40, -1., 1., MVAmethod.c_str()); 
    //MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, s_varofInterest.c_str()); 

  
  	//***********************************************OPEN FILES & GET NTUPLES**********************************************
  	string dataSetName, filepath;
  	int nEntries;
  	float ScaleFactor, NormFactor;

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
		              dataSetName = datasets[d]->Name();
		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = "ObjectVarsTree";
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;
		                
                  /////////////////////////////////////////
                  // Define variables relevant for MVA
                  ////////////////////////////////////////
                  int NumberOfJets;
                  
                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
                  
                  double InclJetCharge[20];
                  double SummedJetCharge[20];
                  double bDiscJet[20];
                  double CvsBJet[20];
                  double CvsLJet[20];
                  double pdgID[20];
                  double MotherpdgID[20];
                  double GrandMotherpdgID[20];

                  ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("charge_jet",&SummedJetCharge);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_Jet",&bDiscJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
		                
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
                  ttree[(dataSetName + "weights").c_str()]->SetBranchAddress("btagWeight_mujets_central",&bTagSF);
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

                  int Train_nEntries = int(nEntries/skipEvents);
                  if(isData && !doTraining) Train_nEntries = int(nEntries/skipEvents);
                  
                  cout << "Number of entries: " << nEntries << ", number of train Entries: " << Train_nEntries << endl;

		              //////////////////////////////////////////////////////////
		              // Running on events
		              //////////////////////////////////////////////////////////

                  if(doTraining)
                  {

                      if(isData && SignalName != "Data") continue;

		                  for (int j = 0; j<Train_nEntries; j++)
		                  {
			                    ttree[(dataSetName).c_str()]->GetEntry(j);

			                            if(dataSetName == SignalName)
			                            {
                                      Eventtrainer_->Fill("S","InclCharge_bjet1", InclJetCharge[0]);
                                      Eventtrainer_->Fill("S","InclCharge_bjet2", InclJetCharge[1]);
                                      Eventtrainer_->Fill("S","InclCharge_bjet3", InclJetCharge[2]);
                                      Eventtrainer_->Fill("S","SummedCharge_bjet1", SummedJetCharge[0]);
                                      Eventtrainer_->Fill("S","SummedCharge_bjet2", SummedJetCharge[1]);
                                      Eventtrainer_->Fill("S","SummedCharge_bjet3", SummedJetCharge[2]);
                                      Eventtrainer_->Fill("S","CvsL_bjet1", CvsBJet[0]);
                                      Eventtrainer_->Fill("S","CvsL_bjet2", CvsBJet[1]);
                                      Eventtrainer_->Fill("S","CvsL_bjet3", CvsBJet[2]);
                                      Eventtrainer_->Fill("S","CvsB_bjet1", CvsLJet[0]);
                                      Eventtrainer_->Fill("S","CvsB_bjet2", CvsLJet[1]);
                                      Eventtrainer_->Fill("S","CvsB_bjet3", CvsLJet[2]);
			                            }
			                            else
			                            {
                                      Eventtrainer_->Fill("B","InclCharge_bjet1", InclJetCharge[0]);
                                      Eventtrainer_->Fill("B","InclCharge_bjet2", InclJetCharge[1]);
                                      Eventtrainer_->Fill("B","InclCharge_bjet3", InclJetCharge[2]);
                                      Eventtrainer_->Fill("B","SummedCharge_bjet1", SummedJetCharge[0]);
                                      Eventtrainer_->Fill("B","SummedCharge_bjet2", SummedJetCharge[1]);
                                      Eventtrainer_->Fill("B","SummedCharge_bjet3", SummedJetCharge[2]);
                                      Eventtrainer_->Fill("B","CvsL_bjet1", CvsBJet[0]);
                                      Eventtrainer_->Fill("B","CvsL_bjet2", CvsBJet[1]);
                                      Eventtrainer_->Fill("B","CvsL_bjet3", CvsBJet[2]);
                                      Eventtrainer_->Fill("B","CvsB_bjet1", CvsLJet[0]);
                                      Eventtrainer_->Fill("B","CvsB_bjet2", CvsLJet[1]);
                                      Eventtrainer_->Fill("B","CvsB_bjet3", CvsLJet[2]);
			                            }
		                  }//for-loop events
                  }//If-statement doTraining
                  else //not training but computing
			            {
			            
			                if(isData) Train_nEntries = -1;
			            		 for (int j = Train_nEntries+1; j<nEntries; j++)
		                  {
                      		ScaleFactor = 1.; // event scale factor
			                    ttree[(dataSetName + "weights").c_str()]->GetEntry(j);
			                    ScaleFactor = ScaleFactor * lumiweight * LeptonSF * bTagSF * nloSF;
			                    if(ScaleFactor < 0) ScaleFactor = 0;
			                    ttree[(dataSetName).c_str()]->GetEntry(j);

                                      if (Eventcomputer_ == 0) cout <<"null computer...." <<endl;
                                      Eventcomputer_->FillVar("InclCharge_bjet1", InclJetCharge[0]);
                                      Eventcomputer_->FillVar("InclCharge_bjet2", InclJetCharge[1]);
                                      Eventcomputer_->FillVar("InclCharge_bjet3", InclJetCharge[2]);
                                      Eventcomputer_->FillVar("SummedCharge_bjet1", SummedJetCharge[0]);
                                      Eventcomputer_->FillVar("SummedCharge_bjet2", SummedJetCharge[1]);
                                      Eventcomputer_->FillVar("SummedCharge_bjet3", SummedJetCharge[2]);
                                      Eventcomputer_->FillVar("CvsL_bjet1", CvsBJet[0]);
                                      Eventcomputer_->FillVar("CvsL_bjet2", CvsBJet[1]);
                                      Eventcomputer_->FillVar("CvsL_bjet3", CvsBJet[2]);
                                      Eventcomputer_->FillVar("CvsB_bjet1", CvsLJet[0]);
                                      Eventcomputer_->FillVar("CvsB_bjet2", CvsLJet[1]);
                                      Eventcomputer_->FillVar("CvsB_bjet3", CvsLJet[2]);
                                      
                                      
                          double BDTscore;

                              std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                          
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                  BDTscore = it->second;
                              }                       
                          if(isData)
                          {
                             MSP_MVA->Fill(BDTscore, datasets[d], true, 1.);
//                             cout << "Filling MVA for data: " << BDTscore << endl;
                          }
                          else
                          {
                             MSP_MVA->Fill(BDTscore, datasets[d], true, ScaleFactor*Luminosity*skipEvents);
//                              cout << "Filling MVA for MC: " << BDTscore << endl;
                          }
                      }
			            }		                
		              
		          }//for-loop datasets
               
      
  if(!doTraining)
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
  
  	TFile *outfile = new TFile((pathPNG+"/OutputMVA.root").c_str(),"recreate");
  	outfile->cd();

  
			cout << "Saving the MSP_MVA" << endl;
      MSP_MVA->setDataLumi(Luminosity);
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
        cout << "Drawing MSP_MVA: " << MSP_MVA << endl;
		    MSP_MVA->Draw("MyMSP_MVA", 1, false, false, false, 1);
      	MSP_MVA->Write(outfile, "MVA" , true,pathPNG, "png");
  
  	    outfile->Write("kOverwrite");
  }
      if(doTraining) Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test",false);
      

      delete Eventtrainer_;
      delete Eventcomputer_;

}



void MVA_JetCombTraining(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod) 
{


  ///////////////////////////////////////////////////////////////
  /// Initializing TopKinematicFit stuff
  //////////////////////////////////////////////////////////////
  cout << "Initializing TopKinFit for MVA training of correct Jet combination" << endl;
  int nToys = 500;

  std::string pdfFileName_SMtt = "TopKinFit_final/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TTSignal = "TopKinFit_final/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_STSignal = "TopKinFit_final/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMtt = new KINFIT::kfit();
  KINFIT::kfit *kf_STSignal = new KINFIT::kfit();
  KINFIT::kfit *kf_TTSignal = new KINFIT::kfit();

  
   if(KinFitMethod == "SMtt")
   {
      kf_SMtt->Init(TOPTOPLEPHAD);
      kf_SMtt->SetPDF("TopWMass",pdfFileName_SMtt.c_str(),"TopWM_Fit");
      kf_SMtt->SetPDF("TopMass",pdfFileName_SMtt.c_str(),"TopM_Fit");
      kf_SMtt->SetPDF("TopWHadMass",pdfFileName_SMtt.c_str(),"TopWHadRecM_Fit");
      kf_SMtt->SetPDF("TopHadMass",pdfFileName_SMtt.c_str(),"TopHadRecM_Fit");
      kf_SMtt->SetPDF("MetPx",pdfFileName_SMtt.c_str(),"dMetPx_Gaus");
      kf_SMtt->SetPDF("MetPy",pdfFileName_SMtt.c_str(),"dMetPy_Gaus");
      kf_SMtt->SetPDF("BJetPx",pdfFileName_SMtt.c_str(),"dBJetPx_Fit");
      kf_SMtt->SetPDF("BJetPy",pdfFileName_SMtt.c_str(),"dBJetPy_Fit");
      kf_SMtt->SetPDF("BJetPz",pdfFileName_SMtt.c_str(),"dBJetPz_Fit");
      kf_SMtt->SetPDF("BJetE",pdfFileName_SMtt.c_str(),"dBJetE_Fit");
      kf_SMtt->SetPDF("ElecPx",pdfFileName_SMtt.c_str(),"dElecPx_Fit");
      kf_SMtt->SetPDF("ElecPy",pdfFileName_SMtt.c_str(),"dElecPy_Fit");
      kf_SMtt->SetPDF("ElecPz",pdfFileName_SMtt.c_str(),"dElecPz_Fit");
      kf_SMtt->SetPDF("ElecE",pdfFileName_SMtt.c_str(),"dElecE_Fit");
      kf_SMtt->SetPDF("MuonPx",pdfFileName_SMtt.c_str(),"dMuonPx_Fit");
      kf_SMtt->SetPDF("MuonPy",pdfFileName_SMtt.c_str(),"dMuonPy_Fit");
      kf_SMtt->SetPDF("MuonPz",pdfFileName_SMtt.c_str(),"dMuonPz_Fit");
      kf_SMtt->SetPDF("MuonE",pdfFileName_SMtt.c_str(),"dMuonE_Fit");
      kf_SMtt->SetPDF("NonBJetPx",pdfFileName_SMtt.c_str(),"dNonBJetPx_Fit");
      kf_SMtt->SetPDF("NonBJetPy",pdfFileName_SMtt.c_str(),"dNonBJetPy_Fit");
      kf_SMtt->SetPDF("NonBJetPz",pdfFileName_SMtt.c_str(),"dNonBJetPz_Fit");
      kf_SMtt->SetPDF("NonBJetE",pdfFileName_SMtt.c_str(),"dNonBJetE_Fit");
      kf_SMtt->SetNToy(nToys);
  }
  if(KinFitMethod == "TTSignal")
  {
      kf_TTSignal->Init(TOPTOPLEPHBB);
      kf_TTSignal->SetPDF("TopWMass",pdfFileName_TTSignal.c_str(),"TopWM_Fit");
      kf_TTSignal->SetPDF("TopMass",pdfFileName_TTSignal.c_str(),"TopM_Fit");
      kf_TTSignal->SetPDF("TopWHadMass",pdfFileName_TTSignal.c_str(),"TopWHadRecM_Fit");
      kf_TTSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_TTSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_TTSignal->SetPDF("MetPx",pdfFileName_TTSignal.c_str(),"dMetPx_Gaus");
      kf_TTSignal->SetPDF("MetPy",pdfFileName_TTSignal.c_str(),"dMetPy_Gaus");
      kf_TTSignal->SetPDF("BJetPx",pdfFileName_TTSignal.c_str(),"dBJetPx_Fit");
      kf_TTSignal->SetPDF("BJetPy",pdfFileName_TTSignal.c_str(),"dBJetPy_Fit");
      kf_TTSignal->SetPDF("BJetPz",pdfFileName_TTSignal.c_str(),"dBJetPz_Fit");
      kf_TTSignal->SetPDF("BJetE",pdfFileName_TTSignal.c_str(),"dBJetE_Fit");
      kf_TTSignal->SetPDF("ElecPx",pdfFileName_TTSignal.c_str(),"dElecPx_Fit");
      kf_TTSignal->SetPDF("ElecPy",pdfFileName_TTSignal.c_str(),"dElecPy_Fit");
      kf_TTSignal->SetPDF("ElecPz",pdfFileName_TTSignal.c_str(),"dElecPz_Fit");
      kf_TTSignal->SetPDF("ElecE",pdfFileName_TTSignal.c_str(),"dElecE_Fit");
      kf_TTSignal->SetPDF("MuonPx",pdfFileName_TTSignal.c_str(),"dMuonPx_Fit");
      kf_TTSignal->SetPDF("MuonPy",pdfFileName_TTSignal.c_str(),"dMuonPy_Fit");
      kf_TTSignal->SetPDF("MuonPz",pdfFileName_TTSignal.c_str(),"dMuonPz_Fit");
      kf_TTSignal->SetPDF("MuonE",pdfFileName_TTSignal.c_str(),"dMuonE_Fit");
      kf_TTSignal->SetPDF("NonBJetPx",pdfFileName_TTSignal.c_str(),"dNonBJetPx_Fit");
      kf_TTSignal->SetPDF("NonBJetPy",pdfFileName_TTSignal.c_str(),"dNonBJetPy_Fit");
      kf_TTSignal->SetPDF("NonBJetPz",pdfFileName_TTSignal.c_str(),"dNonBJetPz_Fit");
      kf_TTSignal->SetPDF("NonBJetE",pdfFileName_TTSignal.c_str(),"dNonBJetE_Fit");
      kf_TTSignal->SetNToy(nToys);
  }
  if (KinFitMethod ==   "STSignal")
  {
      kf_STSignal->Init(TOPHLEPBB);
      kf_STSignal->SetPDF("TopWMass",pdfFileName_STSignal.c_str(),"TopWM_Fit");
      kf_STSignal->SetPDF("TopMass",pdfFileName_STSignal.c_str(),"TopM_Fit");
      kf_STSignal->SetPDF("TopWHadMass",pdfFileName_STSignal.c_str(),"TopWHadRecM_Fit");
      kf_STSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_STSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_STSignal->SetPDF("MetPx",pdfFileName_STSignal.c_str(),"dMetPx_Gaus");
      kf_STSignal->SetPDF("MetPy",pdfFileName_STSignal.c_str(),"dMetPy_Gaus");
      kf_STSignal->SetPDF("BJetPx",pdfFileName_STSignal.c_str(),"dBJetPx_Fit");
      kf_STSignal->SetPDF("BJetPy",pdfFileName_STSignal.c_str(),"dBJetPy_Fit");
      kf_STSignal->SetPDF("BJetPz",pdfFileName_STSignal.c_str(),"dBJetPz_Fit");
      kf_STSignal->SetPDF("BJetE",pdfFileName_STSignal.c_str(),"dBJetE_Fit");
      kf_STSignal->SetPDF("ElecPx",pdfFileName_STSignal.c_str(),"dElecPx_Fit");
      kf_STSignal->SetPDF("ElecPy",pdfFileName_STSignal.c_str(),"dElecPy_Fit");
      kf_STSignal->SetPDF("ElecPz",pdfFileName_STSignal.c_str(),"dElecPz_Fit");
      kf_STSignal->SetPDF("ElecE",pdfFileName_STSignal.c_str(),"dElecE_Fit");
      kf_STSignal->SetPDF("MuonPx",pdfFileName_STSignal.c_str(),"dMuonPx_Fit");
      kf_STSignal->SetPDF("MuonPy",pdfFileName_STSignal.c_str(),"dMuonPy_Fit");
      kf_STSignal->SetPDF("MuonPz",pdfFileName_STSignal.c_str(),"dMuonPz_Fit");
      kf_STSignal->SetPDF("MuonE",pdfFileName_STSignal.c_str(),"dMuonE_Fit");
      kf_STSignal->SetPDF("NonBJetPx",pdfFileName_STSignal.c_str(),"dNonBJetPx_Fit");
      kf_STSignal->SetPDF("NonBJetPy",pdfFileName_STSignal.c_str(),"dNonBJetPy_Fit");
      kf_STSignal->SetPDF("NonBJetPz",pdfFileName_STSignal.c_str(),"dNonBJetPz_Fit");
      kf_STSignal->SetPDF("NonBJetE",pdfFileName_STSignal.c_str(),"dNonBJetE_Fit");
      kf_STSignal->SetNToy(nToys);
  }
   
  ///////////////////////////////////////////////////////////////////
  // Initializing MVA
  ///////////////////////////////////////////////////////////////////
  cout << "Initializing MVA for correct Jet combination" << endl;
  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);

  string pathMVA_ = "MVA/TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);

  MVATrainer* Eventtrainer_ = 0;
  Eventtrainer_ = new MVATrainer(MVAmethod,"TrainedJetCombMVA_"+KinFitMethod+channel, pathMVA_+"TrainedJetCombMVA_"+KinFitMethod+channel+"_"+SignalName+".root");

  vector<std::string> MVAvars;

  if(KinFitMethod == "TTSignal" || KinFitMethod == "SMtt")
  {
      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
      MVAvars.push_back("SumCharge_TopJets");
      MVAvars.push_back("SumCharge_FCNHJetLep");
      MVAvars.push_back("CvsL_Hjet1");
      MVAvars.push_back("CvsL_Hjet2");
      MVAvars.push_back("CvsL_SMb");
      MVAvars.push_back("CvsL_FCNHjet");
      MVAvars.push_back("CvsB_Hjet1");
      MVAvars.push_back("CvsB_Hjet2");
      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("CvsB_FCNHjet");
      MVAvars.push_back("Hmass");
//      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("HadTopmass");
//      MVAvars.push_back("DR_H_HadTop");
//      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("DR_H_SMb");
      MVAvars.push_back("DR_Hb1_Hb2");
//      MVAvars.push_back("DR_Lep_SMb");
      MVAvars.push_back("DR_Lep_H");
      MVAvars.push_back("DR_Lep_HadTop");
      MVAvars.push_back("Chi2");
//      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HPt");
      MVAvars.push_back("HadTopPt");
  }
  else if(KinFitMethod == "STSignal")
  {
      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
      MVAvars.push_back("CvsL_Hjet1");
      MVAvars.push_back("CvsL_Hjet2");
      MVAvars.push_back("CvsL_SMb");
      MVAvars.push_back("CvsB_Hjet1");
      MVAvars.push_back("CvsB_Hjet2");
      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("DR_H_SMb");
      MVAvars.push_back("DR_Hb1_Hb2");
      MVAvars.push_back("DR_Lep_SMb");
      MVAvars.push_back("DR_Lep_H");
      MVAvars.push_back("Chi2");
      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HPt");
  }
  
  if(KinFitMethod == "TTSignal" || KinFitMethod == "SMtt")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_->bookInputVar(MVAvars[N_var]);
      }
  }
  if(KinFitMethod == "STSignal")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_->bookInputVar(MVAvars[N_var]);
      }
  }

  ////////////////////////////////////////////////////////////
  // Load Datasets
  //////////////////////////////////////////////////////////////////////
 	const char *xmlfile = xmlNom.c_str();
 	cout << "Using config file: " << xmlfile << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
	
		dataSetName = datasets[d]->Name();
    if(SignalName != dataSetName) continue;

    bool SingleTop = false;
    if(dataSetName.find("ST")!=string::npos || dataSetName.find("SingleTop")!=string::npos) SingleTop = true;

		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "ObjectVarsTree";
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
    /////////////////////////////////////////
    // Get object variables
    ////////////////////////////////////////
    int NumberOfJets;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
    
    double InclJetCharge[20];
    double bDiscJet[20];
    double CvsBJet[20];
    double CvsLJet[20];
    double pdgID[20];
    double MotherpdgID[20];
    double GrandMotherpdgID[20];
    double lepCharge;
    double lepPt;
    double lepEta;
    double lepPhi;
    double lepE;

    double pt_jet[20];
    double phi_jet[20];
    double eta_jet[20];
    double E_jet[20];

    double met_px;
    double met_py;


    ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
    if(channel == "_Mu")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
    }
    if(channel == "_El")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
    }
		  


		//////////////////////////////////////////////////////////
		// Running on events
		//////////////////////////////////////////////////////////
		for (int j = 0; j<int(nEntries/skipEvents); j++)
		{
	        std::vector<float> BJetPt;
	        std::vector<float> BJetEta;
	        std::vector<float> BJetPhi;
	        std::vector<float> BJetE;

	        std::vector<float> NonBJetPt;
	        std::vector<float> NonBJetEta;
	        std::vector<float> NonBJetPhi;
	        std::vector<float> NonBJetE;

	        std::vector<float> LeptonPt;
	        std::vector<float> LeptonEta;
	        std::vector<float> LeptonPhi;
	        std::vector<float> LeptonE;
	        
	        vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	        vector <int> MapIndex_NonBindex;
	        
			    ttree[(dataSetName).c_str()]->GetEntry(j);
          for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
          {
          
 
                if(bDiscJet[i_Jet]  > workingpointvalue_Medium)
                {
                    BJetPt.push_back(pt_jet[i_Jet]);
                    BJetEta.push_back(eta_jet[i_Jet]);
                    BJetPhi.push_back(phi_jet[i_Jet]);
                    BJetE.push_back(E_jet[i_Jet]);
                    
                    MapIndex_Bindex.push_back(i_Jet);
                }
                else
                {
                    NonBJetPt.push_back(pt_jet[i_Jet]);
                    NonBJetEta.push_back(eta_jet[i_Jet]);
                    NonBJetPhi.push_back(phi_jet[i_Jet]);
                    NonBJetE.push_back(E_jet[i_Jet]);

                    MapIndex_NonBindex.push_back(i_Jet);
                }
          }
          
          LeptonPt.push_back(lepPt);
          LeptonEta.push_back(lepEta);
          LeptonPhi.push_back(lepPhi);
          LeptonE.push_back(lepE);

          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "SMtt")
          {
              if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
              {
              

                  kf_SMtt->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SMtt->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SMtt->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SMtt->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SMtt->SetMet(met_px,met_py);
                            
                            
                  kf_SMtt->Run();
                  int nPerm_SMtt = kf_SMtt->GetNPerm();
                            
                            
                  for(int ip=0;ip<nPerm_SMtt;ip++)
                  {
                      double chi2 = kf_SMtt->GetDisc(ip);

                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SMtt->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kf_SMtt->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kf_SMtt->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kf_SMtt->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];
		                  

                      double nuPx = kf_SMtt->GetNuPx(ip,0);
                      double nuPy = kf_SMtt->GetNuPy(ip,0);
                      double nuPz = kf_SMtt->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJETHAD_, NONBJET1_, NONBJET2_;
                      TLorentzVector Higgs_, HadTop_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETHAD],eta_jet[IndexAllJetColl_BJETHAD],phi_jet[IndexAllJetColl_BJETHAD],E_jet[IndexAllJetColl_BJETHAD]);
                      NONBJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET1],eta_jet[IndexAllJetColl_NONBJET1],phi_jet[IndexAllJetColl_NONBJET1],E_jet[IndexAllJetColl_NONBJET1]);
                      NONBJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET2],eta_jet[IndexAllJetColl_NONBJET2],phi_jet[IndexAllJetColl_NONBJET2],E_jet[IndexAllJetColl_NONBJET2]);
                      
                      Higgs_ = NONBJET1_+NONBJET2_;
                      HadTop_ = Higgs_+BJETHAD_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;

                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
//                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
//                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
//                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = NONBJET1_.DeltaR(NONBJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
                      double Chi2 = chi2;
//                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
//                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("S","Chi2", Chi2);
//                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
//                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("B","Chi2", Chi2);
//                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                          {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
//                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("S","Chi2", Chi2);
//                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
//                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("B","Chi2", Chi2);
//                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHAD selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "TTSignal")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
              

                  kf_TTSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_TTSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_TTSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_TTSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_TTSignal->SetMet(met_px,met_py);
                            
                            
                  kf_TTSignal->Run();
                  int nPerm_TTSignal = kf_TTSignal->GetNPerm();
                            
                            
                  for(int ip=0;ip<nPerm_TTSignal;ip++)
                  {
                      double chi2 = kf_TTSignal->GetDisc(ip);
                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_TTSignal->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kf_TTSignal->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_TTSignal->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_TTSignal->GetIndex(BJET2_TOPTOPLEPHBB,ip)];
                      
                      double nuPx = kf_TTSignal->GetNuPx(ip,0);
                      double nuPy = kf_TTSignal->GetNuPy(ip,0);
                      double nuPz = kf_TTSignal->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, NONBJETHAD_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, HadTop_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      HadTop_ = Higgs_+NONBJETHAD_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
//                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
//                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
//                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
                      double Chi2 = chi2;
//                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
//                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("S","Chi2", Chi2);
//                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
//                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("B","Chi2", Chi2);
//                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                          {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
//                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("S","Chi2", Chi2);
//                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
//                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
//                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
//                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
                              Eventtrainer_->Fill("B","Chi2", Chi2);
//                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHBB selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "STSignal")
          {
              if(BJetPt.size()>=3)
              {
              

                  kf_STSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_STSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_STSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_STSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_STSignal->SetMet(met_px,met_py);
                            
                            
                  kf_STSignal->Run();
                  int nPerm_STSignal = kf_STSignal->GetNPerm();
                            
                            
                  for(int ip=0;ip<nPerm_STSignal;ip++)
                  {
                      double chi2 = kf_STSignal->GetDisc(ip);
                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_STSignal->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_STSignal->GetIndex(BJET1_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_STSignal->GetIndex(BJET2_TOPHLEPBB,ip)];
                      
                      double nuPx = kf_STSignal->GetNuPx(ip,0);
                      double nuPy = kf_STSignal->GetNuPy(ip,0);
                      double nuPz = kf_STSignal->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();

                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HPt", HPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HPt", HPt);
                          }
                  }
              }//TOPHLEPBB selection
          }
		}//for-loop events
  }//for-loop datasets


  Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test",false);
      

  delete Eventtrainer_;
}





void MVA_JetCombComputer(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod) 
{


  ///////////////////////////////////////////////////////////////
  /// Initializing TopKinematicFit stuff
  //////////////////////////////////////////////////////////////
  cout << "Initializing TopKinFit for MVA training of correct Jet combination" << endl;
  int nToys = 500;

  std::string pdfFileName_SMtt = "TopKinFit_final/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TTSignal = "TopKinFit_final/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_STSignal = "TopKinFit_final/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMtt = new KINFIT::kfit();
  KINFIT::kfit *kf_STSignal = new KINFIT::kfit();
  KINFIT::kfit *kf_TTSignal = new KINFIT::kfit();

  
   if(KinFitMethod == "SMtt")
   {
      kf_SMtt->Init(TOPTOPLEPHAD);
      kf_SMtt->SetPDF("TopWMass",pdfFileName_SMtt.c_str(),"TopWM_Fit");
      kf_SMtt->SetPDF("TopMass",pdfFileName_SMtt.c_str(),"TopM_Fit");
      kf_SMtt->SetPDF("TopWHadMass",pdfFileName_SMtt.c_str(),"TopWHadRecM_Fit");
      kf_SMtt->SetPDF("TopHadMass",pdfFileName_SMtt.c_str(),"TopHadRecM_Fit");
      kf_SMtt->SetPDF("MetPx",pdfFileName_SMtt.c_str(),"dMetPx_Gaus");
      kf_SMtt->SetPDF("MetPy",pdfFileName_SMtt.c_str(),"dMetPy_Gaus");
      kf_SMtt->SetPDF("BJetPx",pdfFileName_SMtt.c_str(),"dBJetPx_Fit");
      kf_SMtt->SetPDF("BJetPy",pdfFileName_SMtt.c_str(),"dBJetPy_Fit");
      kf_SMtt->SetPDF("BJetPz",pdfFileName_SMtt.c_str(),"dBJetPz_Fit");
      kf_SMtt->SetPDF("BJetE",pdfFileName_SMtt.c_str(),"dBJetE_Fit");
      kf_SMtt->SetPDF("ElecPx",pdfFileName_SMtt.c_str(),"dElecPx_Fit");
      kf_SMtt->SetPDF("ElecPy",pdfFileName_SMtt.c_str(),"dElecPy_Fit");
      kf_SMtt->SetPDF("ElecPz",pdfFileName_SMtt.c_str(),"dElecPz_Fit");
      kf_SMtt->SetPDF("ElecE",pdfFileName_SMtt.c_str(),"dElecE_Fit");
      kf_SMtt->SetPDF("MuonPx",pdfFileName_SMtt.c_str(),"dMuonPx_Fit");
      kf_SMtt->SetPDF("MuonPy",pdfFileName_SMtt.c_str(),"dMuonPy_Fit");
      kf_SMtt->SetPDF("MuonPz",pdfFileName_SMtt.c_str(),"dMuonPz_Fit");
      kf_SMtt->SetPDF("MuonE",pdfFileName_SMtt.c_str(),"dMuonE_Fit");
      kf_SMtt->SetPDF("NonBJetPx",pdfFileName_SMtt.c_str(),"dNonBJetPx_Fit");
      kf_SMtt->SetPDF("NonBJetPy",pdfFileName_SMtt.c_str(),"dNonBJetPy_Fit");
      kf_SMtt->SetPDF("NonBJetPz",pdfFileName_SMtt.c_str(),"dNonBJetPz_Fit");
      kf_SMtt->SetPDF("NonBJetE",pdfFileName_SMtt.c_str(),"dNonBJetE_Fit");
      kf_SMtt->SetNToy(nToys);
  }
  if(KinFitMethod == "TTSignal")
  {
      kf_TTSignal->Init(TOPTOPLEPHBB);
      kf_TTSignal->SetPDF("TopWMass",pdfFileName_TTSignal.c_str(),"TopWM_Fit");
      kf_TTSignal->SetPDF("TopMass",pdfFileName_TTSignal.c_str(),"TopM_Fit");
      kf_TTSignal->SetPDF("TopWHadMass",pdfFileName_TTSignal.c_str(),"TopWHadRecM_Fit");
      kf_TTSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_TTSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_TTSignal->SetPDF("MetPx",pdfFileName_TTSignal.c_str(),"dMetPx_Gaus");
      kf_TTSignal->SetPDF("MetPy",pdfFileName_TTSignal.c_str(),"dMetPy_Gaus");
      kf_TTSignal->SetPDF("BJetPx",pdfFileName_TTSignal.c_str(),"dBJetPx_Fit");
      kf_TTSignal->SetPDF("BJetPy",pdfFileName_TTSignal.c_str(),"dBJetPy_Fit");
      kf_TTSignal->SetPDF("BJetPz",pdfFileName_TTSignal.c_str(),"dBJetPz_Fit");
      kf_TTSignal->SetPDF("BJetE",pdfFileName_TTSignal.c_str(),"dBJetE_Fit");
      kf_TTSignal->SetPDF("ElecPx",pdfFileName_TTSignal.c_str(),"dElecPx_Fit");
      kf_TTSignal->SetPDF("ElecPy",pdfFileName_TTSignal.c_str(),"dElecPy_Fit");
      kf_TTSignal->SetPDF("ElecPz",pdfFileName_TTSignal.c_str(),"dElecPz_Fit");
      kf_TTSignal->SetPDF("ElecE",pdfFileName_TTSignal.c_str(),"dElecE_Fit");
      kf_TTSignal->SetPDF("MuonPx",pdfFileName_TTSignal.c_str(),"dMuonPx_Fit");
      kf_TTSignal->SetPDF("MuonPy",pdfFileName_TTSignal.c_str(),"dMuonPy_Fit");
      kf_TTSignal->SetPDF("MuonPz",pdfFileName_TTSignal.c_str(),"dMuonPz_Fit");
      kf_TTSignal->SetPDF("MuonE",pdfFileName_TTSignal.c_str(),"dMuonE_Fit");
      kf_TTSignal->SetPDF("NonBJetPx",pdfFileName_TTSignal.c_str(),"dNonBJetPx_Fit");
      kf_TTSignal->SetPDF("NonBJetPy",pdfFileName_TTSignal.c_str(),"dNonBJetPy_Fit");
      kf_TTSignal->SetPDF("NonBJetPz",pdfFileName_TTSignal.c_str(),"dNonBJetPz_Fit");
      kf_TTSignal->SetPDF("NonBJetE",pdfFileName_TTSignal.c_str(),"dNonBJetE_Fit");
      kf_TTSignal->SetNToy(nToys);
  }
  if (KinFitMethod ==   "STSignal")
  {
      kf_STSignal->Init(TOPHLEPBB);
      kf_STSignal->SetPDF("TopWMass",pdfFileName_STSignal.c_str(),"TopWM_Fit");
      kf_STSignal->SetPDF("TopMass",pdfFileName_STSignal.c_str(),"TopM_Fit");
      kf_STSignal->SetPDF("TopWHadMass",pdfFileName_STSignal.c_str(),"TopWHadRecM_Fit");
      kf_STSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_STSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_STSignal->SetPDF("MetPx",pdfFileName_STSignal.c_str(),"dMetPx_Gaus");
      kf_STSignal->SetPDF("MetPy",pdfFileName_STSignal.c_str(),"dMetPy_Gaus");
      kf_STSignal->SetPDF("BJetPx",pdfFileName_STSignal.c_str(),"dBJetPx_Fit");
      kf_STSignal->SetPDF("BJetPy",pdfFileName_STSignal.c_str(),"dBJetPy_Fit");
      kf_STSignal->SetPDF("BJetPz",pdfFileName_STSignal.c_str(),"dBJetPz_Fit");
      kf_STSignal->SetPDF("BJetE",pdfFileName_STSignal.c_str(),"dBJetE_Fit");
      kf_STSignal->SetPDF("ElecPx",pdfFileName_STSignal.c_str(),"dElecPx_Fit");
      kf_STSignal->SetPDF("ElecPy",pdfFileName_STSignal.c_str(),"dElecPy_Fit");
      kf_STSignal->SetPDF("ElecPz",pdfFileName_STSignal.c_str(),"dElecPz_Fit");
      kf_STSignal->SetPDF("ElecE",pdfFileName_STSignal.c_str(),"dElecE_Fit");
      kf_STSignal->SetPDF("MuonPx",pdfFileName_STSignal.c_str(),"dMuonPx_Fit");
      kf_STSignal->SetPDF("MuonPy",pdfFileName_STSignal.c_str(),"dMuonPy_Fit");
      kf_STSignal->SetPDF("MuonPz",pdfFileName_STSignal.c_str(),"dMuonPz_Fit");
      kf_STSignal->SetPDF("MuonE",pdfFileName_STSignal.c_str(),"dMuonE_Fit");
      kf_STSignal->SetPDF("NonBJetPx",pdfFileName_STSignal.c_str(),"dNonBJetPx_Fit");
      kf_STSignal->SetPDF("NonBJetPy",pdfFileName_STSignal.c_str(),"dNonBJetPy_Fit");
      kf_STSignal->SetPDF("NonBJetPz",pdfFileName_STSignal.c_str(),"dNonBJetPz_Fit");
      kf_STSignal->SetPDF("NonBJetE",pdfFileName_STSignal.c_str(),"dNonBJetE_Fit");
      kf_STSignal->SetNToy(nToys);
  }
   
  ///////////////////////////////////////////////////////////////////
  // Initializing MVA
  ///////////////////////////////////////////////////////////////////
  cout << "Initializing MVA for correct Jet combination" << endl;
  string pathMVA = "MVA/";
  mkdir(pathMVA.c_str(),0777);
  pathMVA += "weightstest";
  mkdir(pathMVA.c_str(),0777);

  string pathMVA_ = "MVA/TrainFiles/";
  mkdir(pathMVA_.c_str(),0777);




  MVAComputer* Eventcomputer_ =0;   
  vector<std::string> MVAvars;

  if(KinFitMethod == "TTSignal" || KinFitMethod == "SMtt")
  {
      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
      MVAvars.push_back("SumCharge_TopJets");
      MVAvars.push_back("SumCharge_FCNHJetLep");
      MVAvars.push_back("CvsL_Hjet1");
      MVAvars.push_back("CvsL_Hjet2");
      MVAvars.push_back("CvsL_SMb");
      MVAvars.push_back("CvsL_FCNHjet");
      MVAvars.push_back("CvsB_Hjet1");
      MVAvars.push_back("CvsB_Hjet2");
      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("CvsB_FCNHjet");
      MVAvars.push_back("Hmass");
//      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("HadTopmass");
//      MVAvars.push_back("DR_H_HadTop");
//      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("DR_H_SMb");
      MVAvars.push_back("DR_Hb1_Hb2");
//      MVAvars.push_back("DR_Lep_SMb");
      MVAvars.push_back("DR_Lep_H");
      MVAvars.push_back("DR_Lep_HadTop");
      MVAvars.push_back("Chi2");
//      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HPt");
      MVAvars.push_back("HadTopPt");
  }
  else if(KinFitMethod == "STSignal")
  {
      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
      MVAvars.push_back("CvsL_Hjet1");
      MVAvars.push_back("CvsL_Hjet2");
      MVAvars.push_back("CvsL_SMb");
      MVAvars.push_back("CvsB_Hjet1");
      MVAvars.push_back("CvsB_Hjet2");
      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("DR_H_SMb");
      MVAvars.push_back("DR_Hb1_Hb2");
      MVAvars.push_back("DR_Lep_SMb");
      MVAvars.push_back("DR_Lep_H");
      MVAvars.push_back("Chi2");
      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HPt");
  }
  
  Eventcomputer_ = new MVAComputer(MVAmethod,pathMVA_+"TrainedJetCombMVA_"+KinFitMethod+channel+"_"+SignalName+".root","TrainedJetCombMVA_"+KinFitMethod+channel,MVAvars, "test");

  ////////////////////////////////////////////////////////////
  // Load Datasets
  //////////////////////////////////////////////////////////////////////
 	const char *xmlfile = xmlNom.c_str();
 	cout << "Using config file: " << xmlfile << endl;
  TTreeLoader treeLoader;
  vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  //***********************************************OPEN FILES & GET NTUPLES**********************************************
  string dataSetName, filepath;
  int nEntries;

  
	for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	{
	                  int NSelectionPassedEvents_SMtt = 0;
	                  int NMCIdentifiedEvents_SMtt_STsignal = 0;
	                  int NMCIdentifiedEvents_SMtt_TTsignal = 0;
	                  int NMCIdentifiedEvents_SMtt_TTbackground = 0;
	                  int NSelectionPassedEvents_TTSignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_STsignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_TTsignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_TTbackground = 0;
	                  int NSelectionPassedEvents_STSignal = 0;
	                  int NMCIdentifiedEvents_STSignal_STsignal = 0;
	                  int NMCIdentifiedEvents_STSignal_TTsignal = 0;
	                  int NMCIdentifiedEvents_STSignal_TTbackground = 0;
	              
                    int nMCMatchedPassedEvents_SMtt_TTbackground = 0;
                    int nMCMatchedPassedEvents_SMtt_STsignal = 0;
                    int nMCMatchedPassedEvents_SMtt_TTsignal = 0;
                    int nMCMatchedPassedEvents_STSignal_TTbackground = 0;
                    int nMCMatchedPassedEvents_STSignal_STsignal = 0;
                    int nMCMatchedPassedEvents_STSignal_TTsignal = 0;
                    int nMCMatchedPassedEvents_TTSignal_TTbackground = 0;
                    int nMCMatchedPassedEvents_TTSignal_STsignal = 0;
                    int nMCMatchedPassedEvents_TTSignal_TTsignal = 0;
	              
	
		dataSetName = datasets[d]->Name();

    bool SingleTop = false;
    if(dataSetName.find("ST")!=string::npos || dataSetName.find("SingleTop")!=string::npos) SingleTop = true;

		cout<<"Dataset:  :"<<dataSetName<<endl;
		filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		if (debug) cout<<"filepath: "<<filepath<<endl;
		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		string TTreename = "ObjectVarsTree";
		ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		nEntries = ttree[dataSetName.c_str()]->GetEntries();
		cout<<"                 nEntries: "<<nEntries<<endl;
		  
    /////////////////////////////////////////
    // Get object variables
    ////////////////////////////////////////
    int NumberOfJets;
    
    ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
    
    double InclJetCharge[20];
    double bDiscJet[20];
    double CvsBJet[20];
    double CvsLJet[20];
    double pdgID[20];
    double MotherpdgID[20];
    double GrandMotherpdgID[20];
    double lepCharge;
    double lepPt;
    double lepEta;
    double lepPhi;
    double lepE;

    double pt_jet[20];
    double phi_jet[20];
    double eta_jet[20];
    double E_jet[20];

    double met_px;
    double met_py;


    ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
    ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
    ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
    ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
    if(channel == "_Mu")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
    }
    if(channel == "_El")
    {
        ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
        ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
        ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
        ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
        ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
    }
		  


		//////////////////////////////////////////////////////////
		// Running on events
		//////////////////////////////////////////////////////////
		for (int j = int(nEntries/skipEvents)+1; j<2*int(nEntries/skipEvents); j++)
		{
          int SMTT_Matched = 0;
          int Signal_Matched = 0;
	        
	        std::vector<float> BJetPt;
	        std::vector<float> BJetEta;
	        std::vector<float> BJetPhi;
	        std::vector<float> BJetE;

	        std::vector<float> NonBJetPt;
	        std::vector<float> NonBJetEta;
	        std::vector<float> NonBJetPhi;
	        std::vector<float> NonBJetE;

	        std::vector<float> LeptonPt;
	        std::vector<float> LeptonEta;
	        std::vector<float> LeptonPhi;
	        std::vector<float> LeptonE;
	        
	        vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	        vector <int> MapIndex_NonBindex;

	        
			    ttree[(dataSetName).c_str()]->GetEntry(j);
          for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
          {
          
                if( (MotherpdgID[i_Jet] == 25) || (fabs(MotherpdgID[i_Jet]) == 6) ) Signal_Matched++;
                if( (MotherpdgID[i_Jet] == 24) || (fabs(MotherpdgID[i_Jet]) == 6) ) SMTT_Matched++;
 
                if(bDiscJet[i_Jet]  > workingpointvalue_Medium)
                {
                    BJetPt.push_back(pt_jet[i_Jet]);
                    BJetEta.push_back(eta_jet[i_Jet]);
                    BJetPhi.push_back(phi_jet[i_Jet]);
                    BJetE.push_back(E_jet[i_Jet]);
                    
                    MapIndex_Bindex.push_back(i_Jet);
                }
                else
                {
                    NonBJetPt.push_back(pt_jet[i_Jet]);
                    NonBJetEta.push_back(eta_jet[i_Jet]);
                    NonBJetPhi.push_back(phi_jet[i_Jet]);
                    NonBJetE.push_back(E_jet[i_Jet]);

                    MapIndex_NonBindex.push_back(i_Jet);
                }
          }
          
          LeptonPt.push_back(lepPt);
          LeptonEta.push_back(lepEta);
          LeptonPhi.push_back(lepPhi);
          LeptonE.push_back(lepE);

          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "SMtt")
          {
              if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
              {
              
	                NSelectionPassedEvents_SMtt++;
	                if(SMTT_Matched == 4) nMCMatchedPassedEvents_SMtt_TTbackground++;
	                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_SMtt_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SMtt_STsignal++;

                  kf_SMtt->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SMtt->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SMtt->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SMtt->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SMtt->SetMet(met_px,met_py);
                            
                            
                  kf_SMtt->Run();
                  int nPerm_SMtt = kf_SMtt->GetNPerm();
                            
                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_BJETHAD = -1;
                  int HighestBDT_IndexAllJetColl_NONBJET1 = -1;
                  int HighestBDT_IndexAllJetColl_NONBJET2 = -1;
                            
                  for(int ip=0;ip<nPerm_SMtt;ip++)
                  {
                      double chi2 = kf_SMtt->GetDisc(ip);

                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SMtt->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kf_SMtt->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kf_SMtt->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kf_SMtt->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];
		                  

                      double nuPx = kf_SMtt->GetNuPx(ip,0);
                      double nuPy = kf_SMtt->GetNuPy(ip,0);
                      double nuPz = kf_SMtt->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJETHAD_, NONBJET1_, NONBJET2_;
                      TLorentzVector Higgs_, HadTop_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETHAD],eta_jet[IndexAllJetColl_BJETHAD],phi_jet[IndexAllJetColl_BJETHAD],E_jet[IndexAllJetColl_BJETHAD]);
                      NONBJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET1],eta_jet[IndexAllJetColl_NONBJET1],phi_jet[IndexAllJetColl_NONBJET1],E_jet[IndexAllJetColl_NONBJET1]);
                      NONBJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJET2],eta_jet[IndexAllJetColl_NONBJET2],phi_jet[IndexAllJetColl_NONBJET2],E_jet[IndexAllJetColl_NONBJET2]);
                      
                      Higgs_ = NONBJET1_+NONBJET2_;
                      HadTop_ = Higgs_+BJETHAD_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;

                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
//                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
//                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
//                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = NONBJET1_.DeltaR(NONBJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
                      double Chi2 = chi2;
//                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
//                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
                      Eventcomputer_->FillVar("SumCharge_TopJets", SumCharge_TopJets);
                      Eventcomputer_->FillVar("SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
                      Eventcomputer_->FillVar("CvsL_FCNHjet", CvsL_FCNHjet);
                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
                      Eventcomputer_->FillVar("CvsB_FCNHjet", CvsB_FCNHjet);
                      Eventcomputer_->FillVar("Hmass", Hmass);
//                      Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                      Eventcomputer_->FillVar("HadTopmass", HadTopmass);
//                      Eventcomputer_->FillVar("DR_H_HadTop", DR_H_HadTop);
//                      Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
                      Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
//                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
                      Eventcomputer_->FillVar("DR_Lep_HadTop", DR_Lep_HadTop);
                      Eventcomputer_->FillVar("Chi2", Chi2);
//                      Eventcomputer_->FillVar("LepTopPt", LepTopPt);
                      Eventcomputer_->FillVar("HPt", HPt);
                      Eventcomputer_->FillVar("HadTopPt", HadTopPt);


                      double BDTscore_tmp;
                              std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                          
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                  BDTscore_tmp = it->second;
                              }
                              
                              
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                              HighestBDT_IndexAllJetColl_BJETHAD = IndexAllJetColl_BJETHAD;
                              HighestBDT_IndexAllJetColl_NONBJET1 = IndexAllJetColl_NONBJET1;
                              HighestBDT_IndexAllJetColl_NONBJET2 = IndexAllJetColl_NONBJET2;
                          }


                  }

                  //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_SMtt_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_SMtt_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_SMtt_TTbackground++;
	                                        }
                                      }

              }//TOPTOPLEPHAD selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          else if(KinFitMethod == "TTSignal")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
              

                  kf_TTSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_TTSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_TTSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_TTSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_TTSignal->SetMet(met_px,met_py);
                            
                            
                  kf_TTSignal->Run();
                  int nPerm_TTSignal = kf_TTSignal->GetNPerm();
                            
	                NSelectionPassedEvents_TTSignal++;
	                if(SMTT_Matched == 4) nMCMatchedPassedEvents_TTSignal_TTbackground++;
	                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_TTSignal_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_TTSignal_STsignal++;

                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_NONBJETHAD = -1;
                  int HighestBDT_IndexAllJetColl_BJET1 = -1;
                  int HighestBDT_IndexAllJetColl_BJET2 = -1;

                            
                  for(int ip=0;ip<nPerm_TTSignal;ip++)
                  {
                      double chi2 = kf_TTSignal->GetDisc(ip);
                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_TTSignal->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kf_TTSignal->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_TTSignal->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_TTSignal->GetIndex(BJET2_TOPTOPLEPHBB,ip)];
                      
                      double nuPx = kf_TTSignal->GetNuPx(ip,0);
                      double nuPy = kf_TTSignal->GetNuPy(ip,0);
                      double nuPz = kf_TTSignal->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, NONBJETHAD_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, HadTop_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      NONBJETHAD_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_NONBJETHAD],eta_jet[IndexAllJetColl_NONBJETHAD],phi_jet[IndexAllJetColl_NONBJETHAD],E_jet[IndexAllJetColl_NONBJETHAD]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      HadTop_ = Higgs_+NONBJETHAD_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
//                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
//                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
//                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
                      double Chi2 = chi2;
//                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
//                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
                      Eventcomputer_->FillVar("SumCharge_TopJets", SumCharge_TopJets);
                      Eventcomputer_->FillVar("SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
                      Eventcomputer_->FillVar("CvsL_FCNHjet", CvsL_FCNHjet);
                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
                      Eventcomputer_->FillVar("CvsB_FCNHjet", CvsB_FCNHjet);
                      Eventcomputer_->FillVar("Hmass", Hmass);
//                      Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                      Eventcomputer_->FillVar("HadTopmass", HadTopmass);
//                      Eventcomputer_->FillVar("DR_H_HadTop", DR_H_HadTop);
//                      Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
                      Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
//                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
                      Eventcomputer_->FillVar("DR_Lep_HadTop", DR_Lep_HadTop);
                      Eventcomputer_->FillVar("Chi2", Chi2);
//                      Eventcomputer_->FillVar("LepTopPt", LepTopPt);
                      Eventcomputer_->FillVar("HPt", HPt);
                      Eventcomputer_->FillVar("HadTopPt", HadTopPt);
                          
                      double BDTscore_tmp;
                              std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                          
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                  BDTscore_tmp = it->second;
                              }
                              
                              
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                              HighestBDT_IndexAllJetColl_NONBJETHAD = IndexAllJetColl_NONBJETHAD;
                              HighestBDT_IndexAllJetColl_BJET1 = IndexAllJetColl_BJET1;
                              HighestBDT_IndexAllJetColl_BJET2 = IndexAllJetColl_BJET2;
                          }


                  }

                  //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_TTSignal_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_TTSignal_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_TTSignal_TTbackground++;
	                                        }
                                      }

              }//TOPTOPLEPHBB selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
          /////////////////////////////////////////////////////////////////////////////
          else if(KinFitMethod == "STSignal")
          {
              if(BJetPt.size()>=3)
              {
	                NSelectionPassedEvents_STSignal++;
	                if(SMTT_Matched >= 3) nMCMatchedPassedEvents_STSignal_TTbackground++;
	                if(Signal_Matched >= 3 && !SingleTop) nMCMatchedPassedEvents_STSignal_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_STSignal_STsignal++;
              

                  kf_STSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_STSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_STSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_STSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_STSignal->SetMet(met_px,met_py);
                            
                            
                  kf_STSignal->Run();
                  int nPerm_STSignal = kf_STSignal->GetNPerm();
                            
                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_BJET1 = -1;
                  int HighestBDT_IndexAllJetColl_BJET2 = -1;
                            
                  for(int ip=0;ip<nPerm_STSignal;ip++)
                  {
                      double chi2 = kf_STSignal->GetDisc(ip);
                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_STSignal->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_STSignal->GetIndex(BJET1_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_STSignal->GetIndex(BJET2_TOPHLEPBB,ip)];
                      
                      double nuPx = kf_STSignal->GetNuPx(ip,0);
                      double nuPy = kf_STSignal->GetNuPy(ip,0);
                      double nuPz = kf_STSignal->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
                      double HPt = Higgs_.Pt();

                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
//                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
                      Eventcomputer_->FillVar("Hmass", Hmass);
                      Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                      Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
                      Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
                      Eventcomputer_->FillVar("Chi2", Chi2);
                      Eventcomputer_->FillVar("LepTopPt", LepTopPt);
                      Eventcomputer_->FillVar("HPt", HPt);
                          
                      double BDTscore_tmp;
                              std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                          
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                  BDTscore_tmp = it->second;
                              }
                              
                              
                          if(BDTscore_tmp > BDTscore)
                          {
                              BDTscore = BDTscore_tmp;
                              HighestBDT_IndexAllJetColl_BJETLEP = IndexAllJetColl_BJETLEP;
                              HighestBDT_IndexAllJetColl_BJET1 = IndexAllJetColl_BJET1;
                              HighestBDT_IndexAllJetColl_BJET2 = IndexAllJetColl_BJET2;
                          }


                  }

                  //Counting whether the jet combination with the highest BDT score is matched to the GenLevel correct jet combination
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5))
                                          {
	                                          if(SingleTop) NMCIdentifiedEvents_STSignal_STsignal++;
	                                          else NMCIdentifiedEvents_STSignal_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5))
                                          {
	                                          NMCIdentifiedEvents_STSignal_TTbackground++;
	                                        }


              }//TOPHLEPBB selection
          }
		}//for-loop events
		              

    if(KinFitMethod == "SMtt")
    {
		    cout << "************ TOPTOPLEPHAD ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHAD & " << 100*double(NSelectionPassedEvents_SMtt)/double((2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_STsignal)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_STsignal)/double(nMCMatchedPassedEvents_SMtt_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_STsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_TTsignal)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_TTsignal)/double(nMCMatchedPassedEvents_SMtt_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_TTsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_TTbackground)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_TTbackground)/double(nMCMatchedPassedEvents_SMtt_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_TTbackground)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "TTSignal")
    {
		    cout << "************ TOPTOPLEPHBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHBB & " << 100*double(NSelectionPassedEvents_TTSignal)/double((2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_STsignal)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_STsignal)/double(nMCMatchedPassedEvents_TTSignal_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_STsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_TTsignal)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_TTsignal)/double(nMCMatchedPassedEvents_TTSignal_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_TTsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_TTbackground)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_TTbackground)/double(nMCMatchedPassedEvents_TTSignal_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_TTbackground)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "STSignal")
    {
		    cout << "************ TOPHLEPBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPHLEPBB & " << 100*double(NSelectionPassedEvents_STSignal)/double((2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_STsignal)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_STsignal)/double(nMCMatchedPassedEvents_STSignal_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_STsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_TTsignal)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_TTsignal)/double(nMCMatchedPassedEvents_STSignal_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_TTsignal)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_TTbackground)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_TTbackground)/double(nMCMatchedPassedEvents_STSignal_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_TTbackground)/double(2*int(nEntries/skipEvents)-int(nEntries/skipEvents)+1) << endl;//Total efficiency
        }
    }
    cout << " " << endl;		              
/*		              cout << "************ TOPTOPLEPHAD ************" << endl;
		              cout << "Number of events passed at least 2 b and 2 non b: " <<  NSelectionPassedEvents_SMtt << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_SMtt_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_SMtt_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_SMtt_STsignal << endl;
                  cout << " " << endl;		              
		              cout << "************ TOPTOPLEPHBB ************" << endl;
		              cout << "Number of events passed at least 3 b and 1 non b: " <<  NSelectionPassedEvents_TTSignal << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_TTSignal_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_TTSignal_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_TTSignal_STsignal << endl;
                  cout << " " << endl;		              
		              cout << "************ TOPHLEPBB ************" << endl;
		              cout << "Number of events passed at least 3 b: " <<  NSelectionPassedEvents_STSignal << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_STSignal_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_STSignal_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_STSignal_STsignal << endl;
                  cout << "---------------------------------------------------------" << endl;		              
                  cout << " " << endl;		              
*/
  }//for-loop datasets



  delete Eventcomputer_;
}










void MCAnalysis(std::string xmlNom, TString CraneenPath, string KinFitMethod)
{

   int nToys = 500;

 	const char *xmlfile = xmlNom.c_str();
 	cout << "used config file: " << xmlfile << endl;

  std::string pdfFileName_SMtt = "TopKinFit_final/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TTSignal = "TopKinFit_final/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_STSignal = "TopKinFit_final/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMtt = new KINFIT::kfit();
  KINFIT::kfit *kf_STSignal = new KINFIT::kfit();
  KINFIT::kfit *kf_TTSignal = new KINFIT::kfit();

  
   if(KinFitMethod == "SMtt" || KinFitMethod == "ALL")
   {
      kf_SMtt->Init(TOPTOPLEPHAD);
      kf_SMtt->SetPDF("TopWMass",pdfFileName_SMtt.c_str(),"TopWM_Fit");
      kf_SMtt->SetPDF("TopMass",pdfFileName_SMtt.c_str(),"TopM_Fit");
      kf_SMtt->SetPDF("TopWHadMass",pdfFileName_SMtt.c_str(),"TopWHadRecM_Fit");
      kf_SMtt->SetPDF("TopHadMass",pdfFileName_SMtt.c_str(),"TopHadRecM_Fit");
      kf_SMtt->SetPDF("MetPx",pdfFileName_SMtt.c_str(),"dMetPx_Gaus");
      kf_SMtt->SetPDF("MetPy",pdfFileName_SMtt.c_str(),"dMetPy_Gaus");
      kf_SMtt->SetPDF("BJetPx",pdfFileName_SMtt.c_str(),"dBJetPx_Fit");
      kf_SMtt->SetPDF("BJetPy",pdfFileName_SMtt.c_str(),"dBJetPy_Fit");
      kf_SMtt->SetPDF("BJetPz",pdfFileName_SMtt.c_str(),"dBJetPz_Fit");
      kf_SMtt->SetPDF("BJetE",pdfFileName_SMtt.c_str(),"dBJetE_Fit");
      kf_SMtt->SetPDF("ElecPx",pdfFileName_SMtt.c_str(),"dElecPx_Fit");
      kf_SMtt->SetPDF("ElecPy",pdfFileName_SMtt.c_str(),"dElecPy_Fit");
      kf_SMtt->SetPDF("ElecPz",pdfFileName_SMtt.c_str(),"dElecPz_Fit");
      kf_SMtt->SetPDF("ElecE",pdfFileName_SMtt.c_str(),"dElecE_Fit");
      kf_SMtt->SetPDF("MuonPx",pdfFileName_SMtt.c_str(),"dMuonPx_Fit");
      kf_SMtt->SetPDF("MuonPy",pdfFileName_SMtt.c_str(),"dMuonPy_Fit");
      kf_SMtt->SetPDF("MuonPz",pdfFileName_SMtt.c_str(),"dMuonPz_Fit");
      kf_SMtt->SetPDF("MuonE",pdfFileName_SMtt.c_str(),"dMuonE_Fit");
      kf_SMtt->SetPDF("NonBJetPx",pdfFileName_SMtt.c_str(),"dNonBJetPx_Fit");
      kf_SMtt->SetPDF("NonBJetPy",pdfFileName_SMtt.c_str(),"dNonBJetPy_Fit");
      kf_SMtt->SetPDF("NonBJetPz",pdfFileName_SMtt.c_str(),"dNonBJetPz_Fit");
      kf_SMtt->SetPDF("NonBJetE",pdfFileName_SMtt.c_str(),"dNonBJetE_Fit");
      kf_SMtt->SetNToy(nToys);
  }
  if(KinFitMethod == "TTSignal" || KinFitMethod == "ALL")
  {
      kf_TTSignal->Init(TOPTOPLEPHBB);
      kf_TTSignal->SetPDF("TopWMass",pdfFileName_TTSignal.c_str(),"TopWM_Fit");
      kf_TTSignal->SetPDF("TopMass",pdfFileName_TTSignal.c_str(),"TopM_Fit");
      kf_TTSignal->SetPDF("TopWHadMass",pdfFileName_TTSignal.c_str(),"TopWHadRecM_Fit");
      kf_TTSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_TTSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_TTSignal->SetPDF("MetPx",pdfFileName_TTSignal.c_str(),"dMetPx_Gaus");
      kf_TTSignal->SetPDF("MetPy",pdfFileName_TTSignal.c_str(),"dMetPy_Gaus");
      kf_TTSignal->SetPDF("BJetPx",pdfFileName_TTSignal.c_str(),"dBJetPx_Fit");
      kf_TTSignal->SetPDF("BJetPy",pdfFileName_TTSignal.c_str(),"dBJetPy_Fit");
      kf_TTSignal->SetPDF("BJetPz",pdfFileName_TTSignal.c_str(),"dBJetPz_Fit");
      kf_TTSignal->SetPDF("BJetE",pdfFileName_TTSignal.c_str(),"dBJetE_Fit");
      kf_TTSignal->SetPDF("ElecPx",pdfFileName_TTSignal.c_str(),"dElecPx_Fit");
      kf_TTSignal->SetPDF("ElecPy",pdfFileName_TTSignal.c_str(),"dElecPy_Fit");
      kf_TTSignal->SetPDF("ElecPz",pdfFileName_TTSignal.c_str(),"dElecPz_Fit");
      kf_TTSignal->SetPDF("ElecE",pdfFileName_TTSignal.c_str(),"dElecE_Fit");
      kf_TTSignal->SetPDF("MuonPx",pdfFileName_TTSignal.c_str(),"dMuonPx_Fit");
      kf_TTSignal->SetPDF("MuonPy",pdfFileName_TTSignal.c_str(),"dMuonPy_Fit");
      kf_TTSignal->SetPDF("MuonPz",pdfFileName_TTSignal.c_str(),"dMuonPz_Fit");
      kf_TTSignal->SetPDF("MuonE",pdfFileName_TTSignal.c_str(),"dMuonE_Fit");
      kf_TTSignal->SetPDF("NonBJetPx",pdfFileName_TTSignal.c_str(),"dNonBJetPx_Fit");
      kf_TTSignal->SetPDF("NonBJetPy",pdfFileName_TTSignal.c_str(),"dNonBJetPy_Fit");
      kf_TTSignal->SetPDF("NonBJetPz",pdfFileName_TTSignal.c_str(),"dNonBJetPz_Fit");
      kf_TTSignal->SetPDF("NonBJetE",pdfFileName_TTSignal.c_str(),"dNonBJetE_Fit");
      kf_TTSignal->SetNToy(nToys);
  }
  if (KinFitMethod ==   "STSignal" || KinFitMethod == "ALL")
  {
      kf_STSignal->Init(TOPHLEPBB);
      kf_STSignal->SetPDF("TopWMass",pdfFileName_STSignal.c_str(),"TopWM_Fit");
      kf_STSignal->SetPDF("TopMass",pdfFileName_STSignal.c_str(),"TopM_Fit");
      kf_STSignal->SetPDF("TopWHadMass",pdfFileName_STSignal.c_str(),"TopWHadRecM_Fit");
      kf_STSignal->SetPDF("HiggsMass",pdfFileName_TTSignal.c_str(),"MHiggs_GenMatch_Fit");
      kf_STSignal->SetPDF("TopHadMass",pdfFileName_TTSignal.c_str(),"MHadTop_GenMatch_Fit");
      kf_STSignal->SetPDF("MetPx",pdfFileName_STSignal.c_str(),"dMetPx_Gaus");
      kf_STSignal->SetPDF("MetPy",pdfFileName_STSignal.c_str(),"dMetPy_Gaus");
      kf_STSignal->SetPDF("BJetPx",pdfFileName_STSignal.c_str(),"dBJetPx_Fit");
      kf_STSignal->SetPDF("BJetPy",pdfFileName_STSignal.c_str(),"dBJetPy_Fit");
      kf_STSignal->SetPDF("BJetPz",pdfFileName_STSignal.c_str(),"dBJetPz_Fit");
      kf_STSignal->SetPDF("BJetE",pdfFileName_STSignal.c_str(),"dBJetE_Fit");
      kf_STSignal->SetPDF("ElecPx",pdfFileName_STSignal.c_str(),"dElecPx_Fit");
      kf_STSignal->SetPDF("ElecPy",pdfFileName_STSignal.c_str(),"dElecPy_Fit");
      kf_STSignal->SetPDF("ElecPz",pdfFileName_STSignal.c_str(),"dElecPz_Fit");
      kf_STSignal->SetPDF("ElecE",pdfFileName_STSignal.c_str(),"dElecE_Fit");
      kf_STSignal->SetPDF("MuonPx",pdfFileName_STSignal.c_str(),"dMuonPx_Fit");
      kf_STSignal->SetPDF("MuonPy",pdfFileName_STSignal.c_str(),"dMuonPy_Fit");
      kf_STSignal->SetPDF("MuonPz",pdfFileName_STSignal.c_str(),"dMuonPz_Fit");
      kf_STSignal->SetPDF("MuonE",pdfFileName_STSignal.c_str(),"dMuonE_Fit");
      kf_STSignal->SetPDF("NonBJetPx",pdfFileName_STSignal.c_str(),"dNonBJetPx_Fit");
      kf_STSignal->SetPDF("NonBJetPy",pdfFileName_STSignal.c_str(),"dNonBJetPy_Fit");
      kf_STSignal->SetPDF("NonBJetPz",pdfFileName_STSignal.c_str(),"dNonBJetPz_Fit");
      kf_STSignal->SetPDF("NonBJetE",pdfFileName_STSignal.c_str(),"dNonBJetE_Fit");
      kf_STSignal->SetNToy(nToys);
  }
   
   
  	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

  
  	//***********************************************OPEN FILES & GET NTUPLES**********************************************
  	string dataSetName, filepath;
  	int nEntries;

  
	              for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	              {
	              
	                  int NSelectionPassedEvents_SMtt = 0;
	                  int NMCIdentifiedEvents_SMtt_STsignal = 0;
	                  int NMCIdentifiedEvents_SMtt_TTsignal = 0;
	                  int NMCIdentifiedEvents_SMtt_TTbackground = 0;
	                  int NSelectionPassedEvents_TTSignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_STsignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_TTsignal = 0;
	                  int NMCIdentifiedEvents_TTSignal_TTbackground = 0;
	                  int NSelectionPassedEvents_STSignal = 0;
	                  int NMCIdentifiedEvents_STSignal_STsignal = 0;
	                  int NMCIdentifiedEvents_STSignal_TTsignal = 0;
	                  int NMCIdentifiedEvents_STSignal_TTbackground = 0;
	              
                    int nMCMatchedPassedEvents_SMtt_TTbackground = 0;
                    int nMCMatchedPassedEvents_SMtt_STsignal = 0;
                    int nMCMatchedPassedEvents_SMtt_TTsignal = 0;
                    int nMCMatchedPassedEvents_STSignal_TTbackground = 0;
                    int nMCMatchedPassedEvents_STSignal_STsignal = 0;
                    int nMCMatchedPassedEvents_STSignal_TTsignal = 0;
                    int nMCMatchedPassedEvents_TTSignal_TTbackground = 0;
                    int nMCMatchedPassedEvents_TTSignal_STsignal = 0;
                    int nMCMatchedPassedEvents_TTSignal_TTsignal = 0;
	              
		              dataSetName = datasets[d]->Name();

                  bool SingleTop = false;

                  if(dataSetName.find("ST")!=string::npos || dataSetName.find("SingleTop")!=string::npos) SingleTop = true;

		              cout<<"Dataset:  :"<<dataSetName<<endl;
		              filepath = CraneenPath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		              if (debug) cout<<"filepath: "<<filepath<<endl;
	

		              FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		              string TTreename = "ObjectVarsTree";
		              ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		              nEntries = ttree[dataSetName.c_str()]->GetEntries();
		              cout<<"                 nEntries: "<<nEntries<<endl;
		                
                  /////////////////////////////////////////
                  // Define variables relevant for MVA
                  ////////////////////////////////////////
                  int NumberOfJets;
                  
                  ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&NumberOfJets);
                  
                  double InclJetCharge[20];
                  double SummedJetCharge[20];
                  double bDiscJet[20];
                  double CvsBJet[20];
                  double CvsLJet[20];
                  double pdgID[20];
                  double MotherpdgID[20];
                  double GrandMotherpdgID[20];
                  double lepCharge;
                  double lepPt;
                  double lepEta;
                  double lepPhi;
                  double lepE;

                  double pt_jet[20];
                  double phi_jet[20];
                  double eta_jet[20];
                  double E_jet[20];

                  double met_px;
                  double met_py;


                  ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&InclJetCharge);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("charge_jet",&SummedJetCharge);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("bdisc_jet",&bDiscJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsB_jet",&CvsBJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("cdiscCvsL_jet",&CvsLJet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_pdgID",&pdgID);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_motherpdgID",&MotherpdgID);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("jet_matchedMC_grannypdgID",&GrandMotherpdgID);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("pt_jet",&pt_jet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("phi_jet",&phi_jet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("eta_jet",&eta_jet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("E_jet",&E_jet);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("met_Px",&met_px);
                  ttree[(dataSetName).c_str()]->SetBranchAddress("met_Py",&met_py);
                  if(channel == "_Mu")
                  {
                      ttree[(dataSetName).c_str()]->SetBranchAddress("charge_muon",&lepCharge);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("pt_muon",&lepPt);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("eta_muon",&lepEta);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("phi_muon",&lepPhi);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("E_muon",&lepE);
                  }
                  if(channel == "_El")
                  {
                      ttree[(dataSetName).c_str()]->SetBranchAddress("charge_electron",&lepCharge);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("pt_electron",&lepPt);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("eta_electron",&lepEta);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("phi_electron",&lepPhi);
                      ttree[(dataSetName).c_str()]->SetBranchAddress("E_electron",&lepE);
                  }
		                


                  cout << "Number of entries: " << nEntries << endl;

		              //////////////////////////////////////////////////////////
		              // Running on events
		              //////////////////////////////////////////////////////////
		              for (int j = 0; j<100000; j++)
		              {
	                      std::vector<float> BJetPt;
	                      std::vector<float> BJetEta;
	                      std::vector<float> BJetPhi;
	                      std::vector<float> BJetE;

	                      std::vector<float> NonBJetPt;
	                      std::vector<float> NonBJetEta;
	                      std::vector<float> NonBJetPhi;
	                      std::vector<float> NonBJetE;

	                      std::vector<float> LeptonPt;
	                      std::vector<float> LeptonEta;
	                      std::vector<float> LeptonPhi;
	                      std::vector<float> LeptonE;
	                      
	                      vector <int> MapIndex_Bindex; //first element is the b-jet index.   The second one the index in the jet-collection
	                      vector <int> MapIndex_NonBindex;
	                      
			                  ttree[(dataSetName).c_str()]->GetEntry(j);
                        int SMTT_Matched = 0;
                        int Signal_Matched = 0;
                        for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
                        {
                        
 
                            if( (MotherpdgID[i_Jet] == 25) || (fabs(MotherpdgID[i_Jet]) == 6) ) Signal_Matched++;
                            if( (MotherpdgID[i_Jet] == 24) || (fabs(MotherpdgID[i_Jet]) == 6) ) SMTT_Matched++;

                              if(bDiscJet[i_Jet]  > workingpointvalue_Medium)
                              {
                                  BJetPt.push_back(pt_jet[i_Jet]);
                                  BJetEta.push_back(eta_jet[i_Jet]);
                                  BJetPhi.push_back(phi_jet[i_Jet]);
                                  BJetE.push_back(E_jet[i_Jet]);
                                  
                                  MapIndex_Bindex.push_back(i_Jet);
                              }
                              else
                              {
                                  NonBJetPt.push_back(pt_jet[i_Jet]);
                                  NonBJetEta.push_back(eta_jet[i_Jet]);
                                  NonBJetPhi.push_back(phi_jet[i_Jet]);
                                  NonBJetE.push_back(E_jet[i_Jet]);

                                  MapIndex_NonBindex.push_back(i_Jet);
                              }
                        }
                        
                        LeptonPt.push_back(lepPt);
                        LeptonEta.push_back(lepEta);
                        LeptonPhi.push_back(lepPhi);
                        LeptonE.push_back(lepE);


                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "SMtt" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
                            {

                                  //cout <<  "Passed Event: " << j << ", Number of Bjets: " << BJetPt.size() << ", Number of NonBjets: " << NonBJetPt.size() << "Number of leptons" << LeptonPt.size() << endl;
                                  
	                                NSelectionPassedEvents_SMtt++;
	                                if(SMTT_Matched == 4) nMCMatchedPassedEvents_SMtt_TTbackground++;
	                                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_SMtt_TTsignal++;
	                                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SMtt_STsignal++;

                                      kf_SMtt->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_SMtt->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_SMtt->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_SMtt->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_SMtt->SetMet(met_px,met_py);


                                      kf_SMtt->Run();
                                      int nPerm_SMtt = kf_SMtt->GetNPerm();

                                      double LowestDisc_SMtt = kf_SMtt->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_SMtt = -1;
                                      int IndexBJETHAD_SMtt = -1;
                                      int IndexNonBJET1HAD_SMtt = -1;
                                      int IndexNonBJET2HAD_SMtt = -1;

                                      for(int ip=0;ip<nPerm_SMtt;ip++)
                                      {
                                           double d = kf_SMtt->GetDisc(ip);
                                           if(d>10E+9)  continue;
                                           if(LowestDisc_SMtt != d) continue;

		                                        IndexBJETLEP_SMtt = kf_SMtt->GetIndex(BJETLEP_TOPTOPLEPHAD,ip);
		                                        IndexBJETHAD_SMtt = kf_SMtt->GetIndex(BJETHAD_TOPTOPLEPHAD,ip);
		                                        IndexNonBJET1HAD_SMtt = kf_SMtt->GetIndex(NONBJET1_TOPTOPLEPHAD,ip);
		                                        IndexNonBJET2HAD_SMtt = kf_SMtt->GetIndex(NONBJET2_TOPTOPLEPHAD,ip);
                                      }

                                      //Seeing if the kinematic fit can be matched to the ST-signal
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMtt]] == 25) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMtt]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMtt]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SMtt]]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_SMtt_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMtt]] == 25) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMtt]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMtt]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SMtt]]) == 5) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETHAD_SMtt]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETHAD_SMtt]]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_SMtt_TTsignal++;
	                                        }
                                          if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMtt]] == 24) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMtt]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMtt]]) == 6) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETHAD_SMtt]]) == 6) )
                                          {
	                                          NMCIdentifiedEvents_SMtt_TTbackground++;
	                                        }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHAD
                        }
                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPTOPLEPHbb (selection: at least 3 b-jets and at least 1 non-b jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "TTSignal" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
                            {
	                                NSelectionPassedEvents_TTSignal++;
	                                if(SMTT_Matched == 4) nMCMatchedPassedEvents_TTSignal_TTbackground++;
	                                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_TTSignal_TTsignal++;
	                                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_TTSignal_STsignal++;

                                      kf_TTSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_TTSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_TTSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_TTSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_TTSignal->SetMet(met_px,met_py);


                                      kf_TTSignal->Run();
                                      int nPerm_TTSignal = kf_TTSignal->GetNPerm();

                                      double LowestDisc_TTSignal = kf_TTSignal->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_TTSignal = -1;
                                      int IndexBJETHAD_TTSignal = -1;
                                      int IndexNonBJET1HAD_TTSignal = -1;
                                      int IndexNonBJET2HAD_TTSignal = -1;

                                      for(int ip=0;ip<nPerm_TTSignal;ip++)
                                      {
                                           double d = kf_TTSignal->GetDisc(ip);
                                           if(d>10E+9)  continue;
                                           if(LowestDisc_TTSignal != d) continue;

		                                        IndexBJETLEP_TTSignal = kf_TTSignal->GetIndex(BJETLEP_TOPTOPLEPHBB,ip);
		                                        IndexBJETHAD_TTSignal = kf_TTSignal->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip);
		                                        IndexNonBJET1HAD_TTSignal = kf_TTSignal->GetIndex(BJET1_TOPTOPLEPHBB,ip);
		                                        IndexNonBJET2HAD_TTSignal = kf_TTSignal->GetIndex(BJET2_TOPTOPLEPHBB,ip);
                                      }

                                      //Seeing if the kinematic fit can be matched to the ST-signal
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TTSignal]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TTSignal]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TTSignal]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_TTSignal]]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_TTSignal_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TTSignal]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TTSignal]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TTSignal]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_TTSignal]]) == 5) && (fabs(MotherpdgID[MapIndex_NonBindex[IndexBJETHAD_TTSignal]]) == 6) && (fabs(pdgID[MapIndex_NonBindex[IndexBJETHAD_TTSignal]]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_TTSignal_TTsignal++;
	                                        }
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TTSignal]] == 24) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TTSignal]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TTSignal]]) == 6) && (fabs(MotherpdgID[MapIndex_NonBindex[IndexBJETHAD_TTSignal]]) == 6) )
                                          {
	                                          NMCIdentifiedEvents_TTSignal_TTbackground++;
	                                        }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHbb
                        }
                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPHLEPbb (selection: at least 3 b-jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "STSignal" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=3)
                            {
	                                NSelectionPassedEvents_STSignal++;
	                                if(SMTT_Matched >= 3) nMCMatchedPassedEvents_STSignal_TTbackground++;
	                                if(Signal_Matched >= 3 && !SingleTop) nMCMatchedPassedEvents_STSignal_TTsignal++;
	                                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_STSignal_STsignal++;

                                      kf_STSignal->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_STSignal->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_STSignal->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_STSignal->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_STSignal->SetMet(met_px,met_py);


                                      kf_STSignal->Run();
                                      int nPerm_STSignal = kf_STSignal->GetNPerm();

                                      double LowestDisc_STSignal = kf_STSignal->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_STSignal = -1;
                                      int IndexNonBJET1HAD_STSignal = -1;
                                      int IndexNonBJET2HAD_STSignal = -1;

                                      for(int ip=0;ip<nPerm_STSignal;ip++)
                                      {
                                           double d = kf_STSignal->GetDisc(ip);
                                           if(d>10E+9)  continue;
                                           if(LowestDisc_STSignal != d) continue;

		                                        IndexBJETLEP_STSignal = kf_STSignal->GetIndex(BJETLEP_TOPHLEPBB,ip);
		                                        IndexNonBJET1HAD_STSignal = kf_STSignal->GetIndex(BJET1_TOPHLEPBB,ip);
		                                        IndexNonBJET2HAD_STSignal = kf_STSignal->GetIndex(BJET2_TOPHLEPBB,ip);
                                      }
                                      //Seeing if the kinematic fit can be matched to the ST-signal
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_STSignal]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_STSignal]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_STSignal]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_STSignal]]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_STSignal_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_STSignal]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_STSignal]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_STSignal]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_STSignal]]) == 5))
                                          {
	                                          NMCIdentifiedEvents_STSignal_TTsignal++;
	                                        }
                                          if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_STSignal]] == 24) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_STSignal]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_STSignal]]) == 6) )
                                          {
	                                          NMCIdentifiedEvents_STSignal_TTbackground++;
	                                        }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHAD
                        }
		              }//for-loop events

    if(KinFitMethod == "SMtt" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPTOPLEPHAD ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHAD & " << 100*double(NSelectionPassedEvents_SMtt)/double(100000 ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_STsignal)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_STsignal)/double(nMCMatchedPassedEvents_SMtt_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_STsignal)/double(100000) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_TTsignal)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_TTsignal)/double(nMCMatchedPassedEvents_SMtt_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_TTsignal)/double(100000) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMtt != 0) cout << 100*double(nMCMatchedPassedEvents_SMtt_TTbackground)/double(NSelectionPassedEvents_SMtt) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMtt_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMtt_TTbackground)/double(nMCMatchedPassedEvents_SMtt_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMtt_TTbackground)/double(100000) << endl;//Total efficiency
        }    
    }
    if(KinFitMethod == "TTSignal" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPTOPLEPHBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHBB & " << 100*double(NSelectionPassedEvents_TTSignal)/double(100000 ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_STsignal)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_STsignal)/double(nMCMatchedPassedEvents_TTSignal_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_STsignal)/double(100000) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_TTsignal)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_TTsignal)/double(nMCMatchedPassedEvents_TTSignal_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_TTsignal)/double(100000) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TTSignal != 0) cout << 100*double(nMCMatchedPassedEvents_TTSignal_TTbackground)/double(NSelectionPassedEvents_TTSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TTSignal_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TTSignal_TTbackground)/double(nMCMatchedPassedEvents_TTSignal_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TTSignal_TTbackground)/double(100000) << endl;//Total efficiency
        }    
    }
    if(KinFitMethod == "STSignal" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPHLEPBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPHLEPBB & " << 100*double(NSelectionPassedEvents_STSignal)/double(100000) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_STsignal)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_STsignal)/double(nMCMatchedPassedEvents_STSignal_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_STsignal)/double(100000) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_TTsignal)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_TTsignal)/double(nMCMatchedPassedEvents_STSignal_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_TTsignal)/double(100000) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_STSignal != 0) cout << 100*double(nMCMatchedPassedEvents_STSignal_TTbackground)/double(NSelectionPassedEvents_STSignal) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_STSignal_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_STSignal_TTbackground)/double(nMCMatchedPassedEvents_STSignal_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_STSignal_TTbackground)/double(100000) << endl;//Total efficiency
        }
    }
    cout << " " << endl;		              
		              
/*
		              cout << "************ TOPTOPLEPHAD ************" << endl;
		              cout << "Number of events passed at least 2 b and 2 non b: " <<  NSelectionPassedEvents_SMtt << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_SMtt_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_SMtt_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_SMtt_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_SMtt_STsignal << endl;
                  cout << " " << endl;		              
		              cout << "************ TOPTOPLEPHBB ************" << endl;
		              cout << "Number of events passed at least 3 b and 1 non b: " <<  NSelectionPassedEvents_TTSignal << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_TTSignal_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_TTSignal_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_TTSignal_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_TTSignal_STsignal << endl;
                  cout << " " << endl;		              
		              cout << "************ TOPHLEPBB ************" << endl;
		              cout << "Number of events passed at least 3 b: " <<  NSelectionPassedEvents_STSignal << endl;
		              cout << "Number of SM tt with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_TTbackground << endl;
		              cout << "Number of TTsignal with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_TTsignal << endl;
		              cout << "Number of STsignal with correct jet-combination: " <<  NMCIdentifiedEvents_STSignal_STsignal << endl;
                  cout << " MATCHED BEFORE EVENT SELECTION " << endl;
	                cout << "  - SM tt: " <<  nMCMatchedPassedEvents_STSignal_TTbackground << endl;
	                cout << "  - TTsignal: " <<  nMCMatchedPassedEvents_STSignal_TTsignal << endl;
	                cout << "  - STsignal: " <<  nMCMatchedPassedEvents_STSignal_STsignal << endl;
                  cout << "---------------------------------------------------------" << endl;		              
                  cout << " " << endl;		              
*/
		         }//for-loop datasets

}

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

