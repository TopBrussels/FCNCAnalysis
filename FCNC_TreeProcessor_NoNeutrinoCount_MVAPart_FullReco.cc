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
#include "FCNCAnalysis_76X/TopKinFit/kinfit.h"


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
void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string JetCombSelection);
void MVA_JetCombTraining_FullReco(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod);
void MVA_JetCombTraining_PartReco(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod);
void MVA_JetCombComputer(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod);
void MCAnalysis(std::string xmlNom, TString CraneenPath, string KinFitMethod);
void MSPCreator ();



int main()
{

    clock_t start = clock();

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


    MVA_JetCombTraining_FullReco( "BDT", 4, "NP_overlay_TTtoTHToBB-1L-Kappa-hct", xmlFileName, TreePath, "TThypo");
      MVA_JetCombTraining_PartReco( "BDT", 2, "NP_overlay_TTtoTHToBB-1L-Kappa-hct", xmlFileName, TreePath, "TThypo");
    MVA_JetCombComputer( "BDT", 2, "NP_overlay_TTtoTHToBB-1L-Kappa-hct", xmlFileName, TreePath, "TThypo");
//    MCAnalysis(xmlFileName, TreePath, "SMttHypo");
	MSPCreator ();

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

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





void MVA_JetCombTraining_FullReco(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod) 
{


  ///////////////////////////////////////////////////////////////
  /// Initializing TopKinematicFit stuff
  //////////////////////////////////////////////////////////////
  cout << "Initializing TopKinFit for MVA training of correct Jet combination" << endl;
  int nToys = 500;

  std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMttHypo = new KINFIT::kfit();
  KINFIT::kfit *kf_SThypo = new KINFIT::kfit();
  KINFIT::kfit *kf_TThypo = new KINFIT::kfit();

  

   if(KinFitMethod == "SMttHypo")
   {
      kf_SMttHypo->Init(TOPTOPLEPHAD);
      kf_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
      kf_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
      kf_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
      kf_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
      kf_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
      kf_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
      kf_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
      kf_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
      kf_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
      kf_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
      kf_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
      kf_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
      kf_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
      kf_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
      kf_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
      kf_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
      kf_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
      kf_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
      kf_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
      kf_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
      kf_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
      kf_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
      kf_SMttHypo->SetNToy(nToys);
  }
  if(KinFitMethod == "TThypo")
  {
      kf_TThypo->Init(TOPTOPLEPHBB);
      kf_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
      kf_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
      kf_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
      kf_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
      kf_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
      kf_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
      kf_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
      kf_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
      kf_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
      kf_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
      kf_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
      kf_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
      kf_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
      kf_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
      kf_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
      kf_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
      kf_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
      kf_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
      kf_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
      kf_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
      kf_TThypo->SetNToy(nToys);
  }
  if (KinFitMethod ==   "SThypo")
  {
      kf_SThypo->Init(TOPHLEPBB);
      kf_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
      kf_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_SThypo.c_str(),"HiggsRecM_Fit");
      kf_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
      kf_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
      kf_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
      kf_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
      kf_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
      kf_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
      kf_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
      kf_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
      kf_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
      kf_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
      kf_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
      kf_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
      kf_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
      kf_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
      kf_SThypo->SetNToy(nToys);
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

  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
//      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
//      MVAvars.push_back("SumCharge_TopJets");
//      MVAvars.push_back("SumCharge_FCNHJetLep");
//      MVAvars.push_back("CvsL_Hjet1");
//      MVAvars.push_back("CvsL_Hjet2");
//      MVAvars.push_back("CvsL_SMb");
//      MVAvars.push_back("CvsL_FCNHjet");
//      MVAvars.push_back("CvsB_Hjet1");
//      MVAvars.push_back("CvsB_Hjet2");
//      MVAvars.push_back("CvsB_SMb");
//      MVAvars.push_back("CvsB_FCNHjet");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("HadTopmass");
      MVAvars.push_back("DR_H_HadTop");
      MVAvars.push_back("DR_H_LepTop");
//      MVAvars.push_back("DR_H_SMb");
//      MVAvars.push_back("DR_Hb1_Hb2");
//      MVAvars.push_back("DR_Lep_SMb");
//      MVAvars.push_back("DR_Lep_H");
//      MVAvars.push_back("DR_Lep_HadTop");
//      MVAvars.push_back("Chi2");
      MVAvars.push_back("LepTopPt");
//      MVAvars.push_back("HPt");
      MVAvars.push_back("HadTopPt");
  }
  else if(KinFitMethod == "SThypo")
  {
//      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep"); //Not part of any training
//      MVAvars.push_back("CvsL_Hjet1"); //Not part of 7 variable training
//      MVAvars.push_back("CvsL_Hjet2"); //Not part of 7 variable training
//      MVAvars.push_back("CvsL_SMb"); //Not part of 7 variable training
//      MVAvars.push_back("CvsB_Hjet1");
//      MVAvars.push_back("CvsB_Hjet2");
//      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("Hmass"); //KIRILL 4 VAR
      MVAvars.push_back("LepTopmass"); //KIRILL 4 VAR
      MVAvars.push_back("DR_H_LepTop"); //Not part of 12 variable training or lower //KIRILL 4 VAR
//      MVAvars.push_back("DR_H_SMb");  //Not part of 7 variable training
//      MVAvars.push_back("DR_Hb1_Hb2"); //Not part of 12 variable training or lower
//      MVAvars.push_back("DR_Lep_SMb"); //Not part of 12 variable training or lower
//      MVAvars.push_back("DR_Lep_H");
//      MVAvars.push_back("Chi2"); //Not part of 12 variable training or lower
      MVAvars.push_back("LepTopPt"); //Not part of 12 variable training or lower //KIRILL 4 VAR
//      MVAvars.push_back("HPt"); //Not part of 7 variable training
  }
  
  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_->bookInputVar(MVAvars[N_var]);
      }
  }
  if(KinFitMethod == "SThypo")
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
          if(KinFitMethod == "SMttHypo")
          {
              if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
              {


                  kf_SMttHypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SMttHypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SMttHypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SMttHypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SMttHypo->SetMet(met_px,met_py);


                  kf_SMttHypo->Run();
                  int nPerm_SMttHypo = kf_SMttHypo->GetNPerm();


                  for(int ip=0;ip<nPerm_SMttHypo;ip++)
                  {
                      double chi2 = kf_SMttHypo->GetDisc(ip);

                      if(chi2>10E+9)  continue;

		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];


                      double nuPx = kf_SMttHypo->GetNuPx(ip,0);
                      double nuPy = kf_SMttHypo->GetNuPy(ip,0);
                      double nuPz = kf_SMttHypo->GetNuPz(ip,0);
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

//                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
//                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
//                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
//                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
//                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
//                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
//                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
//                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
//                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
//                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
//                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
//                      double DR_Hb1_Hb2 = NONBJET1_.DeltaR(NONBJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
//                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
//                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
//                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
//                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
//                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
//                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                          {
//                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
//                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHAD selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "TThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {


                  kf_TThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_TThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_TThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_TThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_TThypo->SetMet(met_px,met_py);


                  kf_TThypo->Run();
                  int nPerm_TThypo = kf_TThypo->GetNPerm();


                  for(int ip=0;ip<nPerm_TThypo;ip++)
                  {
                      double chi2 = kf_TThypo->GetDisc(ip);
                      if(chi2>10E+9)  continue;

		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_TThypo->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kf_TThypo->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                      double nuPx = kf_TThypo->GetNuPx(ip,0);
                      double nuPy = kf_TThypo->GetNuPy(ip,0);
                      double nuPz = kf_TThypo->GetNuPz(ip,0);
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
                      
//                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
//                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
//                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
//                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
//                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
//                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
//                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
//                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
//                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
//                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
//                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
//                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
//                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
//                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
//                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
//                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
//                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
//                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                          {
//                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("S","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("S","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("S","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("S","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("S","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("S","HPt", HPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
//                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("B","SumCharge_TopJets", SumCharge_TopJets);
//                              Eventtrainer_->Fill("B","SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("B","CvsL_FCNHjet", CvsL_FCNHjet);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
//                              Eventtrainer_->Fill("B","CvsB_FCNHjet", CvsB_FCNHjet);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("B","DR_Lep_HadTop", DR_Lep_HadTop);
//                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("B","HPt", HPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHBB selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "SThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
              

                  kf_SThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SThypo->SetMet(met_px,met_py);
                            
                            
                  kf_SThypo->Run();
                  int nPerm_SThypo = kf_SThypo->GetNPerm();
                            
                            
                  for(int ip=0;ip<nPerm_SThypo;ip++)
                  {
                      double chi2 = kf_SThypo->GetDisc(ip);
                      if(chi2>10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SThypo->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET1_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET2_TOPHLEPBB,ip)];
                      
                      double nuPx = kf_SThypo->GetNuPx(ip,0);
                      double nuPy = kf_SThypo->GetNuPy(ip,0);
                      double nuPz = kf_SThypo->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
//                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
//                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
//                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
//                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
//                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
//                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
//                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
//                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
//                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
//                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
//                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
//                      double HPt = Higgs_.Pt();

                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
//                              Eventtrainer_->Fill("S","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("S","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("S","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("S","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("S","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("S","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("S","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("S","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("S","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("S","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("S","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("S","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("S","Chi2", Chi2);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("S","HPt", HPt);
                          }
                          else
                          {
//                              Eventtrainer_->Fill("B","SumCharge_Hjets", SumCharge_Hjets);
//                              Eventtrainer_->Fill("B","SumCharge_SMbLep", SumCharge_SMbLep);
//                              Eventtrainer_->Fill("B","CvsL_Hjet1", CvsL_Hjet1);
//                              Eventtrainer_->Fill("B","CvsL_Hjet2", CvsL_Hjet2);
//                              Eventtrainer_->Fill("B","CvsL_SMb", CvsL_SMb);
//                              Eventtrainer_->Fill("B","CvsB_Hjet1", CvsB_Hjet1);
//                              Eventtrainer_->Fill("B","CvsB_Hjet2", CvsB_Hjet2);
//                              Eventtrainer_->Fill("B","CvsB_SMb", CvsB_SMb);
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
//                              Eventtrainer_->Fill("B","DR_H_SMb", DR_H_SMb);
//                              Eventtrainer_->Fill("B","DR_Hb1_Hb2", DR_Hb1_Hb2);
//                              Eventtrainer_->Fill("B","DR_Lep_SMb", DR_Lep_SMb);
//                              Eventtrainer_->Fill("B","DR_Lep_H", DR_Lep_H);
//                              Eventtrainer_->Fill("B","Chi2", Chi2);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
//                              Eventtrainer_->Fill("B","HPt", HPt);
                          }
                  }
              }//TOPHLEPBB selection
          }
		}//for-loop events
  }//for-loop datasets


  Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test",false);
      

  delete Eventtrainer_;
}




void MVA_JetCombTraining_PartReco(std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string KinFitMethod) 
{


  ///////////////////////////////////////////////////////////////
  /// Initializing TopKinematicFit stuff
  //////////////////////////////////////////////////////////////
  cout << "Initializing TopKinFit for MVA training of correct Jet combination" << endl;
  int nToys = 500;

  std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMttHypo = new KINFIT::kfit();
  KINFIT::kfit *kf_SThypo = new KINFIT::kfit();
  KINFIT::kfit *kf_TThypo = new KINFIT::kfit();

  

   if(KinFitMethod == "SMttHypo")
   {
      kf_SMttHypo->Init(TOPTOPLEPHAD);
      kf_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
      kf_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
      kf_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
      kf_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
      kf_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
      kf_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
      kf_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
      kf_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
      kf_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
      kf_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
      kf_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
      kf_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
      kf_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
      kf_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
      kf_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
      kf_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
      kf_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
      kf_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
      kf_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
      kf_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
      kf_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
      kf_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
      kf_SMttHypo->SetNToy(nToys);
  }
  if(KinFitMethod == "TThypo")
  {
      kf_TThypo->Init(TOPTOPLEPHBB);
      kf_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
      kf_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
      kf_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
      kf_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
      kf_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
      kf_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
      kf_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
      kf_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
      kf_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
      kf_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
      kf_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
      kf_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
      kf_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
      kf_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
      kf_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
      kf_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
      kf_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
      kf_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
      kf_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
      kf_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
      kf_TThypo->SetNToy(nToys);
  }
  if (KinFitMethod ==   "SThypo")
  {
      kf_SThypo->Init(TOPHLEPBB);
      kf_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
      kf_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_SThypo.c_str(),"HiggsRecM_Fit");
      kf_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
      kf_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
      kf_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
      kf_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
      kf_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
      kf_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
      kf_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
      kf_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
      kf_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
      kf_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
      kf_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
      kf_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
      kf_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
      kf_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
      kf_SThypo->SetNToy(nToys);
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
  Eventtrainer_ = new MVATrainer(MVAmethod,"TrainedJetCombMVA_PartReco_"+KinFitMethod+channel, pathMVA_+"TrainedJetCombMVA_PartReco_"+KinFitMethod+channel+"_"+SignalName+".root");

  vector<std::string> MVAvars;

  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("HadTopmass");
      MVAvars.push_back("DR_H_HadTop");
      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("LepTopPt");
      MVAvars.push_back("HadTopPt");
  }
  else if(KinFitMethod == "SThypo")
  {
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("DR_H_LepTop");
      MVAvars.push_back("LepTopPt");
  }
  
  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
      for(unsigned int N_var = 0; N_var < MVAvars.size(); N_var++)
      {
          	Eventtrainer_->bookInputVar(MVAvars[N_var]);
      }
  }
  if(KinFitMethod == "SThypo")
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
          if(KinFitMethod == "SMttHypo")
          {
              if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
              {


                  kf_SMttHypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SMttHypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SMttHypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SMttHypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SMttHypo->SetMet(met_px,met_py);


                  kf_SMttHypo->Run();
                  int nPerm_SMttHypo = kf_SMttHypo->GetNPerm();


                  for(int ip=0;ip<nPerm_SMttHypo;ip++)
                  {
                      double chi2 = kf_SMttHypo->GetDisc(ip);

                      if(chi2<10E+9)  continue;

		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];


                      double nuPx = kf_SMttHypo->GetNuPx(ip,0);
                      double nuPy = kf_SMttHypo->GetNuPy(ip,0);
                      double nuPz = kf_SMttHypo->GetNuPz(ip,0);
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

                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                      double HadTopmass = HadTop_.M();
                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                      double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                      double LepTopPt = LepTop_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                          {
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHAD selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "TThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {


                  kf_TThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_TThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_TThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_TThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_TThypo->SetMet(met_px,met_py);


                  kf_TThypo->Run();
                  int nPerm_TThypo = kf_TThypo->GetNPerm();


                  for(int ip=0;ip<nPerm_TThypo;ip++)
                  {
                      double chi2 = kf_TThypo->GetDisc(ip);
                      if(chi2<10E+9)  continue;

		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_TThypo->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kf_TThypo->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                      double nuPx = kf_TThypo->GetNuPx(ip,0);
                      double nuPy = kf_TThypo->GetNuPy(ip,0);
                      double nuPz = kf_TThypo->GetNuPz(ip,0);
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
                      
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                      double HadTopmass = HadTop_.M();
                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                      double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                      double LepTopPt = LepTop_.Pt();
                      double HadTopPt = HadTop_.Pt();

                      if(SingleTop)
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                      else
                      {
                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                          {
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("S","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("S","HadTopPt", HadTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","HadTopmass", HadTopmass);
                              Eventtrainer_->Fill("B","DR_H_HadTop", DR_H_HadTop);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
                              Eventtrainer_->Fill("B","HadTopPt", HadTopPt);
                          }
                      }
                  }
              }//TOPTOPLEPHBB selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
          /////////////////////////////////////////////////////////////////////////////
          if(KinFitMethod == "SThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
              

                  kf_SThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SThypo->SetMet(met_px,met_py);
                            
                            
                  kf_SThypo->Run();
                  int nPerm_SThypo = kf_SThypo->GetNPerm();
                            
                            
                  for(int ip=0;ip<nPerm_SThypo;ip++)
                  {
                      double chi2 = kf_SThypo->GetDisc(ip);
                      if(chi2<10E+9)  continue;
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SThypo->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET1_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET2_TOPHLEPBB,ip)];
                      
                      double nuPx = kf_SThypo->GetNuPx(ip,0);
                      double nuPy = kf_SThypo->GetNuPy(ip,0);
                      double nuPz = kf_SThypo->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                      double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                      double LepTopPt = LepTop_.Pt();

                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                          {
                              Eventtrainer_->Fill("S","Hmass", Hmass);
                              Eventtrainer_->Fill("S","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("S","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("S","LepTopPt", LepTopPt);
                          }
                          else
                          {
                              Eventtrainer_->Fill("B","Hmass", Hmass);
                              Eventtrainer_->Fill("B","LepTopmass", LepTopmass);
                              Eventtrainer_->Fill("B","DR_H_LepTop", DR_H_LepTop);
                              Eventtrainer_->Fill("B","LepTopPt", LepTopPt);
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

  std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMttHypo = new KINFIT::kfit();
  KINFIT::kfit *kf_SThypo = new KINFIT::kfit();
  KINFIT::kfit *kf_TThypo = new KINFIT::kfit();

  
   if(KinFitMethod == "SMttHypo")
   {
      kf_SMttHypo->Init(TOPTOPLEPHAD);
      kf_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
      kf_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
      kf_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
      kf_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
      kf_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
      kf_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
      kf_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
      kf_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
      kf_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
      kf_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
      kf_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
      kf_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
      kf_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
      kf_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
      kf_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
      kf_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
      kf_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
      kf_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
      kf_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
      kf_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
      kf_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
      kf_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
      kf_SMttHypo->SetNToy(nToys);
  }
  if(KinFitMethod == "TThypo")
  {
      kf_TThypo->Init(TOPTOPLEPHBB);
      kf_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
      kf_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
      kf_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
      kf_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
      kf_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
      kf_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
      kf_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
      kf_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
      kf_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
      kf_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
      kf_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
      kf_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
      kf_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
      kf_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
      kf_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
      kf_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
      kf_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
      kf_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
      kf_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
      kf_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
      kf_TThypo->SetNToy(nToys);
  }
  if (KinFitMethod ==   "SThypo")
  {
      kf_SThypo->Init(TOPHLEPBB);
      kf_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
      kf_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
      kf_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
      kf_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
      kf_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
      kf_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
      kf_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
      kf_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
      kf_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
      kf_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
      kf_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
      kf_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
      kf_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
      kf_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
      kf_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
      kf_SThypo->SetNToy(nToys);
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
  MVAComputer* Eventcomputer_PartReco_ =0;   
  vector<std::string> MVAvars;

  if(KinFitMethod == "TThypo" || KinFitMethod == "SMttHypo")
  {
//      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
//      MVAvars.push_back("SumCharge_TopJets");
//      MVAvars.push_back("SumCharge_FCNHJetLep");
//      MVAvars.push_back("CvsL_Hjet1");
//      MVAvars.push_back("CvsL_Hjet2");
//      MVAvars.push_back("CvsL_SMb");
//      MVAvars.push_back("CvsL_FCNHjet");
//      MVAvars.push_back("CvsB_Hjet1");
//      MVAvars.push_back("CvsB_Hjet2");
//      MVAvars.push_back("CvsB_SMb");
//      MVAvars.push_back("CvsB_FCNHjet");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("HadTopmass");
      MVAvars.push_back("DR_H_HadTop");
      MVAvars.push_back("DR_H_LepTop");
//      MVAvars.push_back("DR_H_SMb");
//      MVAvars.push_back("DR_Hb1_Hb2");
//      MVAvars.push_back("DR_Lep_SMb");
//      MVAvars.push_back("DR_Lep_H");
//      MVAvars.push_back("DR_Lep_HadTop");
//      MVAvars.push_back("Chi2");
      MVAvars.push_back("LepTopPt");
//      MVAvars.push_back("HPt");
      MVAvars.push_back("HadTopPt");
  }
  else if(KinFitMethod == "SThypo")
  {
//      MVAvars.push_back("SumCharge_Hjets");
//      MVAvars.push_back("SumCharge_SMbLep");
//      MVAvars.push_back("CvsL_Hjet1");
//      MVAvars.push_back("CvsL_Hjet2");
//      MVAvars.push_back("CvsL_SMb");
//      MVAvars.push_back("CvsB_Hjet1");
//      MVAvars.push_back("CvsB_Hjet2");
//      MVAvars.push_back("CvsB_SMb");
      MVAvars.push_back("Hmass");
      MVAvars.push_back("LepTopmass");
      MVAvars.push_back("DR_H_LepTop");
//      MVAvars.push_back("DR_H_SMb");
//      MVAvars.push_back("DR_Hb1_Hb2");
//      MVAvars.push_back("DR_Lep_SMb");
//      MVAvars.push_back("DR_Lep_H");
//      MVAvars.push_back("Chi2");
      MVAvars.push_back("LepTopPt");
//      MVAvars.push_back("HPt");
  }
  
  Eventcomputer_ = new MVAComputer(MVAmethod,pathMVA_+"TrainedJetCombMVA_"+KinFitMethod+channel+"_"+SignalName+".root","TrainedJetCombMVA_"+KinFitMethod+channel,MVAvars, "test");
  Eventcomputer_PartReco_ = new MVAComputer(MVAmethod,pathMVA_+"TrainedJetCombMVA_PartReco_"+KinFitMethod+channel+"_"+SignalName+".root","TrainedJetCombMVA_PartReco_"+KinFitMethod+channel,MVAvars, "test");

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
	                  int NumberOfEvents = 0;
	                  
	                  //TopTopLepHad initializations
	                  int NSelectionPassedEvents_SMttHypo = 0;
	                  int NMCIdentifiedEvents_SMttHypo_STsignal = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_SMttHypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SMttHypo_TTbackground_TKF = 0;

	                  //TopTopLepHbb initializations
	                  int NSelectionPassedEvents_TThypo = 0;
	                  int NMCIdentifiedEvents_TThypo_STsignal = 0;
	                  int NMCIdentifiedEvents_TThypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_TThypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_TThypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_TThypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_TThypo_TTbackground_TKF = 0;

	                  //TopHLepbb initializations
	                  int NSelectionPassedEvents_SThypo = 0;
	                  int NMCIdentifiedEvents_SThypo_STsignal = 0;
	                  int NMCIdentifiedEvents_SThypo_TTsignal = 0;
	                  int NMCIdentifiedEvents_SThypo_TTbackground = 0;
	                  int NMCIdentifiedEvents_SThypo_STsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SThypo_TTsignal_TKF = 0;
	                  int NMCIdentifiedEvents_SThypo_TTbackground_TKF = 0;


	              
                    int nMCMatchedPassedEvents_SMttHypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_SMttHypo_STsignal = 0;
                    int nMCMatchedPassedEvents_SMttHypo_TTsignal = 0;
                    int nMCMatchedPassedEvents_SThypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_SThypo_STsignal = 0;
                    int nMCMatchedPassedEvents_SThypo_TTsignal = 0;
                    int nMCMatchedPassedEvents_TThypo_TTbackground = 0;
                    int nMCMatchedPassedEvents_TThypo_STsignal = 0;
                    int nMCMatchedPassedEvents_TThypo_TTsignal = 0;

              
	
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
		
		      NumberOfEvents++;
		
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
          if(KinFitMethod == "SMttHypo")
          {
              if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
              {
              
	                NSelectionPassedEvents_SMttHypo++;
	                if(SMTT_Matched == 4) nMCMatchedPassedEvents_SMttHypo_TTbackground++;
	                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_SMttHypo_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SMttHypo_STsignal++;

                  kf_SMttHypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SMttHypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SMttHypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SMttHypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SMttHypo->SetMet(met_px,met_py);
                            
                            
                  kf_SMttHypo->Run();
                  int nPerm_SMttHypo = kf_SMttHypo->GetNPerm();
                            
                  double LowestDisc_SMttHypo = kf_SMttHypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_BJETHAD = -1;
                  int HighestBDT_IndexAllJetColl_NONBJET1 = -1;
                  int HighestBDT_IndexAllJetColl_NONBJET2 = -1;
                            
                  for(int ip=0;ip<nPerm_SMttHypo;ip++)
                  {
                      double chi2 = kf_SMttHypo->GetDisc(ip);

		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETLEP_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_BJETHAD = MapIndex_Bindex[kf_SMttHypo->GetIndex(BJETHAD_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET1 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET1_TOPTOPLEPHAD,ip)];
		                  int IndexAllJetColl_NONBJET2 = MapIndex_NonBindex[kf_SMttHypo->GetIndex(NONBJET2_TOPTOPLEPHAD,ip)];
		                  
		                  //Counting whether the TopKinFit method can match to the correct signal
                      if(chi2 == LowestDisc_SMttHypo)
                      {
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_SMttHypo_STsignal_TKF++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_SMttHypo_TTsignal_TKF++;
	                                        }
                                          if( (MotherpdgID[IndexAllJetColl_NONBJET1] == 24) && (MotherpdgID[IndexAllJetColl_NONBJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_SMttHypo_TTbackground_TKF++;
	                                        }
                                      }                      
                      } 

                      double nuPx = kf_SMttHypo->GetNuPx(ip,0);
                      double nuPy = kf_SMttHypo->GetNuPy(ip,0);
                      double nuPz = kf_SMttHypo->GetNuPz(ip,0);
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


                      double BDTscore_tmp = -9999.;
                      if(chi2<10E+9)
                      {

    //                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_NONBJET1]-InclJetCharge[IndexAllJetColl_NONBJET2]);
    //                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
    //                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_BJETHAD]);
    //                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_BJETHAD]-lepCharge);
    //                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_NONBJET1];
    //                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_NONBJET2];
    //                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
    //                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_BJETHAD];
    //                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_NONBJET1];
    //                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_NONBJET2];
    //                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
    //                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_BJETHAD];
                          double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                          double LepTopmass = LepTop_.M();
                          double HadTopmass = HadTop_.M();
                          double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                          double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
    //                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
    //                      double DR_Hb1_Hb2 = NONBJET1_.DeltaR(NONBJET2_);
    //                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
    //                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
    //                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
    //                      double Chi2 = chi2;
                          double LepTopPt = LepTop_.Pt();
    //                      double HPt = Higgs_.Pt();
                          double HadTopPt = HadTop_.Pt();

    //                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
    //                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
    //                      Eventcomputer_->FillVar("SumCharge_TopJets", SumCharge_TopJets);
    //                      Eventcomputer_->FillVar("SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
    //                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
    //                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
    //                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
    //                      Eventcomputer_->FillVar("CvsL_FCNHjet", CvsL_FCNHjet);
    //                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
    //                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
    //                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
    //                      Eventcomputer_->FillVar("CvsB_FCNHjet", CvsB_FCNHjet);
                          Eventcomputer_->FillVar("Hmass", Hmass);
                          Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                          Eventcomputer_->FillVar("HadTopmass", HadTopmass);
                          Eventcomputer_->FillVar("DR_H_HadTop", DR_H_HadTop);
                          Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
    //                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
     //                     Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
    //                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
    //                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
    //                      Eventcomputer_->FillVar("DR_Lep_HadTop", DR_Lep_HadTop);
    //                      Eventcomputer_->FillVar("Chi2", Chi2);
                          Eventcomputer_->FillVar("LepTopPt", LepTopPt);
    //                      Eventcomputer_->FillVar("HPt", HPt);
                          Eventcomputer_->FillVar("HadTopPt", HadTopPt);


                                  std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                              
                                  for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                                  {                          
                                      //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                      BDTscore_tmp = it->second;
                                  }
                      }                              
                      else//MVA_PartReco
                      {

                          double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                          double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                          double HadTopmass = HadTop_.M();
                          double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                          double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                          double LepTopPt = LepTop_.Pt();
                          double HadTopPt = HadTop_.Pt();

                          Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                          Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                          Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                          Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                          Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                          Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                          Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                                  std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              
                                  for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                                  {                          
                                      //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                      BDTscore_tmp = it->second;
                                  }
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
	                                            NMCIdentifiedEvents_SMttHypo_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_SMttHypo_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_NONBJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_SMttHypo_TTbackground++;
	                                        }
                                      }

              }//TOPTOPLEPHAD selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPTOPLEPHBB (selection: at least 3 b-jets and at least 1 non-b jets)
          /////////////////////////////////////////////////////////////////////////////
          else if(KinFitMethod == "TThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
	                NSelectionPassedEvents_TThypo++;
	                if(SMTT_Matched == 4) nMCMatchedPassedEvents_TThypo_TTbackground++;
	                if(Signal_Matched == 4 && !SingleTop) nMCMatchedPassedEvents_TThypo_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_TThypo_STsignal++;
              

                  kf_TThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_TThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_TThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_TThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_TThypo->SetMet(met_px,met_py);
                            
                            
                  kf_TThypo->Run();
                  int nPerm_TThypo = kf_TThypo->GetNPerm();
                            

                  double LowestDisc_TThypo = kf_TThypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_NONBJETHAD = -1;
                  int HighestBDT_IndexAllJetColl_BJET1 = -1;
                  int HighestBDT_IndexAllJetColl_BJET2 = -1;

                            
                  for(int ip=0;ip<nPerm_TThypo;ip++)
                  {
                      double chi2 = kf_TThypo->GetDisc(ip);
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_TThypo->GetIndex(BJETLEP_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_NONBJETHAD = MapIndex_NonBindex[kf_TThypo->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET1_TOPTOPLEPHBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_TThypo->GetIndex(BJET2_TOPTOPLEPHBB,ip)];

                      if(chi2 == LowestDisc_TThypo)
                      {
                                      if(SingleTop)
                                      {
                                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) )
                                          {
	                                            NMCIdentifiedEvents_TThypo_STsignal_TKF++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_TThypo_TTsignal_TKF++;
	                                        }
                                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 24) && (MotherpdgID[IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[IndexAllJetColl_NONBJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_TThypo_TTbackground_TKF++;
	                                        }
                                      }
                      }

                      
                      double nuPx = kf_TThypo->GetNuPx(ip,0);
                      double nuPy = kf_TThypo->GetNuPy(ip,0);
                      double nuPz = kf_TThypo->GetNuPz(ip,0);
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
                      
                      double BDTscore_tmp = -9999.;
                      if(chi2<10E+9)
                      {
//                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
//                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
//                      double SumCharge_TopJets = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-InclJetCharge[IndexAllJetColl_NONBJETHAD]);
//                      double SumCharge_FCNHJetLep = fabs(InclJetCharge[IndexAllJetColl_NONBJETHAD]-lepCharge);
//                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
//                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
//                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
//                      double CvsL_FCNHjet = CvsLJet[IndexAllJetColl_NONBJETHAD];
//                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
//                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
//                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
//                      double CvsB_FCNHjet = CvsBJet[IndexAllJetColl_NONBJETHAD];
                      double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                      double LepTopmass = LepTop_.M();
                      double HadTopmass = HadTop_.M();
                      double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                      double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
//                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
//                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
//                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
//                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
//                      double DR_Lep_HadTop = LEP_.DeltaR(HadTop_);
//                      double Chi2 = chi2;
                      double LepTopPt = LepTop_.Pt();
//                      double HPt = Higgs_.Pt();
                      double HadTopPt = HadTop_.Pt();

//                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
//                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
//                      Eventcomputer_->FillVar("SumCharge_TopJets", SumCharge_TopJets);
//                      Eventcomputer_->FillVar("SumCharge_FCNHJetLep", SumCharge_FCNHJetLep);
//                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
//                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
//                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
//                      Eventcomputer_->FillVar("CvsL_FCNHjet", CvsL_FCNHjet);
//                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
//                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
//                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
//                      Eventcomputer_->FillVar("CvsB_FCNHjet", CvsB_FCNHjet);
                      Eventcomputer_->FillVar("Hmass", Hmass);
                      Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                      Eventcomputer_->FillVar("HadTopmass", HadTopmass);
                      Eventcomputer_->FillVar("DR_H_HadTop", DR_H_HadTop);
                      Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
//                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
//                      Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
//                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
//                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
//                      Eventcomputer_->FillVar("DR_Lep_HadTop", DR_Lep_HadTop);
//                      Eventcomputer_->FillVar("Chi2", Chi2);
                      Eventcomputer_->FillVar("LepTopPt", LepTopPt);
//                      Eventcomputer_->FillVar("HPt", HPt);
                      Eventcomputer_->FillVar("HadTopPt", HadTopPt);
                          
                              std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                          
                              for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                              {                          
                                  //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                  BDTscore_tmp = it->second;
                              }
                      }                              
                      else//MVA_PartReco
                      {

                          double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                          double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                          double HadTopmass = HadTop_.M();
                          double DR_H_HadTop = Higgs_.DeltaR(HadTop_);
                          double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                          double LepTopPt = LepTop_.Pt();
                          double HadTopPt = HadTop_.Pt();

                          Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                          Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                          Eventcomputer_PartReco_->FillVar("HadTopmass", HadTopmass);
                          Eventcomputer_PartReco_->FillVar("DR_H_HadTop", DR_H_HadTop);
                          Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                          Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                          Eventcomputer_PartReco_->FillVar("HadTopPt", HadTopPt);


                                  std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              
                                  for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                                  {                          
                                      //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                      BDTscore_tmp = it->second;
                                  }
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
	                                            NMCIdentifiedEvents_TThypo_STsignal++;
                                          }
                                      }
                                      else
                                      {
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 25) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) != 5) )
                                          {
	                                          NMCIdentifiedEvents_TThypo_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_NONBJETHAD]) == 5) )
                                          {
	                                          NMCIdentifiedEvents_TThypo_TTbackground++;
	                                        }
                                      }

              }//TOPTOPLEPHBB selection
          }
          /////////////////////////////////////////////////////////////////////////////
          // Section for the kinematic fit using the TOPHLEPBB (selection: at least 3 b-jets)
          /////////////////////////////////////////////////////////////////////////////
          else if(KinFitMethod == "SThypo")
          {
              if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
              {
	                NSelectionPassedEvents_SThypo++;
	                if(SMTT_Matched >= 3) nMCMatchedPassedEvents_SThypo_TTbackground++;
	                if(Signal_Matched >= 3 && !SingleTop) nMCMatchedPassedEvents_SThypo_TTsignal++;
	                if(Signal_Matched == 3 && SingleTop) nMCMatchedPassedEvents_SThypo_STsignal++;
              

                  kf_SThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                  kf_SThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                  if(channel == "_El") kf_SThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  else if(channel == "_Mu") kf_SThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                  kf_SThypo->SetMet(met_px,met_py);
                            
                            
                  kf_SThypo->Run();
                  int nPerm_SThypo = kf_SThypo->GetNPerm();
                            
                  double LowestDisc_SThypo = kf_SThypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                  double BDTscore = -9999.;
                  int HighestBDT_IndexAllJetColl_BJETLEP = -1;
                  int HighestBDT_IndexAllJetColl_BJET1 = -1;
                  int HighestBDT_IndexAllJetColl_BJET2 = -1;
                            
                  for(int ip=0;ip<nPerm_SThypo;ip++)
                  {
                      double chi2 = kf_SThypo->GetDisc(ip);
		                   
		                  int IndexAllJetColl_BJETLEP = MapIndex_Bindex[kf_SThypo->GetIndex(BJETLEP_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET1 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET1_TOPHLEPBB,ip)];
		                  int IndexAllJetColl_BJET2 = MapIndex_Bindex[kf_SThypo->GetIndex(BJET2_TOPHLEPBB,ip)];

                      if(chi2 == LowestDisc_SThypo)
                      {
                                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 25) && (MotherpdgID[IndexAllJetColl_BJET2]==25) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5))
                                          {
	                                          if(SingleTop) NMCIdentifiedEvents_SThypo_STsignal_TKF++;
	                                          else NMCIdentifiedEvents_SThypo_TTsignal_TKF++;
	                                        }
                                          if( (MotherpdgID[IndexAllJetColl_BJET1] == 24) && (MotherpdgID[IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[IndexAllJetColl_BJETLEP]) == 5))
                                          {
	                                          NMCIdentifiedEvents_SThypo_TTbackground_TKF++;
	                                        }
                      }
                      
                      double nuPx = kf_SThypo->GetNuPx(ip,0);
                      double nuPy = kf_SThypo->GetNuPy(ip,0);
                      double nuPz = kf_SThypo->GetNuPz(ip,0);
                      TLorentzVector Nu_, LEP_, BJETLEP_, BJET1_, BJET2_;
                      TLorentzVector Higgs_, LepTop_;
                      
                      Nu_.SetPxPyPzE(nuPx,nuPy,nuPz,sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz));
                      LEP_.SetPtEtaPhiE(lepPt,lepEta,lepPhi,lepE);
                      BJETLEP_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJETLEP],eta_jet[IndexAllJetColl_BJETLEP],phi_jet[IndexAllJetColl_BJETLEP],E_jet[IndexAllJetColl_BJETLEP]);
                      BJET1_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET1],eta_jet[IndexAllJetColl_BJET1],phi_jet[IndexAllJetColl_BJET1],E_jet[IndexAllJetColl_BJET1]);
                      BJET2_.SetPtEtaPhiE(pt_jet[IndexAllJetColl_BJET2],eta_jet[IndexAllJetColl_BJET2],phi_jet[IndexAllJetColl_BJET2],E_jet[IndexAllJetColl_BJET2]);
                      
                      Higgs_ = BJET1_+BJET2_;
                      LepTop_ = Nu_ + LEP_ + BJETLEP_;
                      
                      double BDTscore_tmp = -9999.;
                      if(chi2<10E+9)
                      {
    //                      double SumCharge_Hjets = fabs(InclJetCharge[IndexAllJetColl_BJET1]-InclJetCharge[IndexAllJetColl_BJET2]);
    //                      double SumCharge_SMbLep = fabs(InclJetCharge[IndexAllJetColl_BJETLEP]-lepCharge);
    //                      double CvsL_Hjet1 = CvsLJet[IndexAllJetColl_BJET1];
    //                      double CvsL_Hjet2 = CvsLJet[IndexAllJetColl_BJET2];
    //                      double CvsL_SMb = CvsLJet[IndexAllJetColl_BJETLEP];
    //                      double CvsB_Hjet1 = CvsBJet[IndexAllJetColl_BJET1];
    //                      double CvsB_Hjet2 = CvsBJet[IndexAllJetColl_BJET2];
    //                      double CvsB_SMb = CvsBJet[IndexAllJetColl_BJETLEP];
                          double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                          double LepTopmass = LepTop_.M();
                          double DR_H_LepTop = Higgs_.DeltaR(LepTop_);
    //                      double DR_H_SMb = Higgs_.DeltaR(BJETLEP_);
    //                      double DR_Hb1_Hb2 = BJET1_.DeltaR(BJET2_);
    //                      double DR_Lep_SMb = LEP_.DeltaR(BJETLEP_);
    //                      double DR_Lep_H = LEP_.DeltaR(Higgs_);
    //                      double Chi2 = chi2;
                          double LepTopPt = LepTop_.Pt();
    //                      double HPt = Higgs_.Pt();

    //                      Eventcomputer_->FillVar("SumCharge_Hjets", SumCharge_Hjets);
    //                      Eventcomputer_->FillVar("SumCharge_SMbLep", SumCharge_SMbLep);
    //                      Eventcomputer_->FillVar("CvsL_Hjet1", CvsL_Hjet1);
    //                      Eventcomputer_->FillVar("CvsL_Hjet2", CvsL_Hjet2);
    //                      Eventcomputer_->FillVar("CvsL_SMb", CvsL_SMb);
    //                      Eventcomputer_->FillVar("CvsB_Hjet1", CvsB_Hjet1);
    //                      Eventcomputer_->FillVar("CvsB_Hjet2", CvsB_Hjet2);
    //                      Eventcomputer_->FillVar("CvsB_SMb", CvsB_SMb);
                          Eventcomputer_->FillVar("Hmass", Hmass);
                          Eventcomputer_->FillVar("LepTopmass", LepTopmass);
                          Eventcomputer_->FillVar("DR_H_LepTop", DR_H_LepTop);
    //                      Eventcomputer_->FillVar("DR_H_SMb", DR_H_SMb);
    //                      Eventcomputer_->FillVar("DR_Hb1_Hb2", DR_Hb1_Hb2);
    //                      Eventcomputer_->FillVar("DR_Lep_SMb", DR_Lep_SMb);
    //                      Eventcomputer_->FillVar("DR_Lep_H", DR_Lep_H);
    //                      Eventcomputer_->FillVar("Chi2", Chi2);
                          Eventcomputer_->FillVar("LepTopPt", LepTopPt);
    //                      Eventcomputer_->FillVar("HPt", HPt);
                              
                                  std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
                              
                                  for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                                  {                          
                                      //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                      BDTscore_tmp = it->second;
                                  }
                      }
                      else        
                      {
                          double Hmass = Higgs_.M(); //The second index indicates for idx=0 leptonic part and idx = 1 hadronic part
                          double LepTopmass = sqrt(2*Higgs_.Pt() * LepTop_.Pt() * (1-cos( Higgs_.DeltaPhi( LepTop_ )) ) );
                          double DR_H_LepTop = Higgs_.DeltaPhi(LepTop_);
                          double LepTopPt = LepTop_.Pt();

                          Eventcomputer_PartReco_->FillVar("Hmass", Hmass);
                          Eventcomputer_PartReco_->FillVar("LepTopmass", LepTopmass);
                          Eventcomputer_PartReco_->FillVar("DR_H_LepTop", DR_H_LepTop);
                          Eventcomputer_PartReco_->FillVar("LepTopPt", LepTopPt);
                              
                          double BDTscore_tmp;
                                  std::map<std::string,Float_t> MVAVals = Eventcomputer_PartReco_->GetMVAValues();
                              
                                  for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
                                  {                          
                                      //cout <<"MVA Method : "<< it->first    <<" Score: "<< it->second <<endl;
                                      BDTscore_tmp = it->second;
                                  }
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
	                                          if(SingleTop) NMCIdentifiedEvents_SThypo_STsignal++;
	                                          else NMCIdentifiedEvents_SThypo_TTsignal++;
	                                        }
                                          if( (MotherpdgID[HighestBDT_IndexAllJetColl_BJET1] == 24) && (MotherpdgID[HighestBDT_IndexAllJetColl_BJET2]==24) && (fabs(MotherpdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 6) && (fabs(pdgID[HighestBDT_IndexAllJetColl_BJETLEP]) == 5))
                                          {
	                                          NMCIdentifiedEvents_SThypo_TTbackground++;
	                                        }


              }//TOPHLEPBB selection
          }
		}//for-loop events
		              

    if(KinFitMethod == "SMttHypo")
    {
		    cout << "************ TOPTOPLEPHAD ************" << endl;
        
        //cout efficiencies for TopKinFit method
        cout << " TOPTOPLEPHAD & " << 100*double(NSelectionPassedEvents_SMttHypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_STsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal_TKF)/double(nMCMatchedPassedEvents_SMttHypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_SMttHypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTbackground)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_SMttHypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }    


        //cout efficiencies for MVA method
        cout << " MVA & " << 100*double(NSelectionPassedEvents_SMttHypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_STsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(nMCMatchedPassedEvents_SMttHypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(nMCMatchedPassedEvents_SMttHypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTbackground)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(nMCMatchedPassedEvents_SMttHypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "TThypo")
    {
		    cout << "************ TOPTOPLEPHBB ************" << endl;

        //cout TopKinFit method
        cout << " TOPTOPLEPHBB & " << 100*double(NSelectionPassedEvents_TThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_STsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal_TKF)/double(nMCMatchedPassedEvents_TThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_TThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTbackground)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_TThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }    


        //cout MVA method
        cout << " MVA & " << 100*double(NSelectionPassedEvents_TThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_STsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(nMCMatchedPassedEvents_TThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(nMCMatchedPassedEvents_TThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTbackground)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(nMCMatchedPassedEvents_TThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }    
    }
    else if(KinFitMethod == "SThypo")
    {
		    cout << "************ TOPHLEPBB ************" << endl;

        //cout TopKinFit method efficiencies
        cout << " TOPHLEPBB & " << 100*double(NSelectionPassedEvents_SThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_STsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal_TKF)/double(nMCMatchedPassedEvents_SThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal_TKF)/double(nMCMatchedPassedEvents_SThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal_TKF)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTbackground)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground_TKF)/double(nMCMatchedPassedEvents_SThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground_TKF)/double(NumberOfEvents) << endl;//Total efficiency
        }


        //cout MVA method efficiencies
        cout << " MVA & " << 100*double(NSelectionPassedEvents_SThypo)/double(NumberOfEvents) << " & ";//selection efficiency
        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_STsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(nMCMatchedPassedEvents_SThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(nMCMatchedPassedEvents_SThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTbackground)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(nMCMatchedPassedEvents_SThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }
    }
    cout << " " << endl;		              
  }//for-loop datasets



  delete Eventcomputer_;
}










void MCAnalysis(std::string xmlNom, TString CraneenPath, string KinFitMethod)
{

   int nToys = 500;

 	const char *xmlfile = xmlNom.c_str();
 	cout << "used config file: " << xmlfile << endl;

  std::string pdfFileName_SMttHypo = "TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SMttHypo = new KINFIT::kfit();
  KINFIT::kfit *kf_SThypo = new KINFIT::kfit();
  KINFIT::kfit *kf_TThypo = new KINFIT::kfit();

  
   if(KinFitMethod == "SMttHypo" || KinFitMethod == "ALL")
   {
      kf_SMttHypo->Init(TOPTOPLEPHAD);
      kf_SMttHypo->SetPDF("TopWMass",pdfFileName_SMttHypo.c_str(),"TopLepWM_Fit");
      kf_SMttHypo->SetPDF("TopMass",pdfFileName_SMttHypo.c_str(),"TopLepRecM_Fit");
      kf_SMttHypo->SetPDF("TopWHadMass",pdfFileName_SMttHypo.c_str(),"TopHadWRecM_Fit");
      kf_SMttHypo->SetPDF("TopHadMass",pdfFileName_SMttHypo.c_str(),"TopHadRecM_Fit");
      kf_SMttHypo->SetPDF("MetPx",pdfFileName_SMttHypo.c_str(),"dMetPx_Gaus");
      kf_SMttHypo->SetPDF("MetPy",pdfFileName_SMttHypo.c_str(),"dMetPy_Gaus");
      kf_SMttHypo->SetPDF("BJetPx",pdfFileName_SMttHypo.c_str(),"dBJetPx_Fit");
      kf_SMttHypo->SetPDF("BJetPy",pdfFileName_SMttHypo.c_str(),"dBJetPy_Fit");
      kf_SMttHypo->SetPDF("BJetPz",pdfFileName_SMttHypo.c_str(),"dBJetPz_Fit");
      kf_SMttHypo->SetPDF("BJetE",pdfFileName_SMttHypo.c_str(),"dBJetE_Fit");
      kf_SMttHypo->SetPDF("ElecPx",pdfFileName_SMttHypo.c_str(),"dElecPx_Fit");
      kf_SMttHypo->SetPDF("ElecPy",pdfFileName_SMttHypo.c_str(),"dElecPy_Fit");
      kf_SMttHypo->SetPDF("ElecPz",pdfFileName_SMttHypo.c_str(),"dElecPz_Fit");
      kf_SMttHypo->SetPDF("ElecE",pdfFileName_SMttHypo.c_str(),"dElecE_Fit");
      kf_SMttHypo->SetPDF("MuonPx",pdfFileName_SMttHypo.c_str(),"dMuonPx_Fit");
      kf_SMttHypo->SetPDF("MuonPy",pdfFileName_SMttHypo.c_str(),"dMuonPy_Fit");
      kf_SMttHypo->SetPDF("MuonPz",pdfFileName_SMttHypo.c_str(),"dMuonPz_Fit");
      kf_SMttHypo->SetPDF("MuonE",pdfFileName_SMttHypo.c_str(),"dMuonE_Fit");
      kf_SMttHypo->SetPDF("NonBJetPx",pdfFileName_SMttHypo.c_str(),"dNonBJetPx_Fit");
      kf_SMttHypo->SetPDF("NonBJetPy",pdfFileName_SMttHypo.c_str(),"dNonBJetPy_Fit");
      kf_SMttHypo->SetPDF("NonBJetPz",pdfFileName_SMttHypo.c_str(),"dNonBJetPz_Fit");
      kf_SMttHypo->SetPDF("NonBJetE",pdfFileName_SMttHypo.c_str(),"dNonBJetE_Fit");
      kf_SMttHypo->SetNToy(nToys);
  }
  if(KinFitMethod == "TThypo" || KinFitMethod == "ALL")
  {
      kf_TThypo->Init(TOPTOPLEPHBB);
      kf_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
      kf_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
      kf_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
      kf_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
      kf_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
      kf_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
      kf_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
      kf_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
      kf_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
      kf_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
      kf_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
      kf_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
      kf_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
      kf_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
      kf_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
      kf_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
      kf_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
      kf_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
      kf_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
      kf_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
      kf_TThypo->SetNToy(nToys);
  }
  if (KinFitMethod ==   "SThypo" || KinFitMethod == "ALL")
  {
      kf_SThypo->Init(TOPHLEPBB);
      kf_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
      kf_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
      kf_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
      kf_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
      kf_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
      kf_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
      kf_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
      kf_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
      kf_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
      kf_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
      kf_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
      kf_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
      kf_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
      kf_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
      kf_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
      kf_SThypo->SetNToy(nToys);
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
	                int NumberOfEvents = 0;
	              
	              
	                int NSelectionPassedEvents_SMttHypo = 0;
	                int NMCIdentifiedEvents_SMttHypo_STsignal = 0;
	                int NMCIdentifiedEvents_SMttHypo_TTsignal = 0;
	                int NMCIdentifiedEvents_SMttHypo_TTbackground = 0;
	                int NSelectionPassedEvents_TThypo = 0;
	                int NMCIdentifiedEvents_TThypo_STsignal = 0;
	                int NMCIdentifiedEvents_TThypo_TTsignal = 0;
	                int NMCIdentifiedEvents_TThypo_TTbackground = 0;
	                int NSelectionPassedEvents_SThypo = 0;
	                int NMCIdentifiedEvents_SThypo_STsignal = 0;
	                int NMCIdentifiedEvents_SThypo_TTsignal = 0;
	                int NMCIdentifiedEvents_SThypo_TTbackground = 0;
	              
                  int nMCMatchedPassedEvents_SMttHypo_TTbackground = 0;
                  int nMCMatchedPassedEvents_SMttHypo_STsignal = 0;
                  int nMCMatchedPassedEvents_SMttHypo_TTsignal = 0;
                  int nMCMatchedPassedEvents_SThypo_TTbackground = 0;
                  int nMCMatchedPassedEvents_SThypo_STsignal = 0;
                  int nMCMatchedPassedEvents_SThypo_TTsignal = 0;
                  int nMCMatchedPassedEvents_TThypo_TTbackground = 0;
                  int nMCMatchedPassedEvents_TThypo_STsignal = 0;
                  int nMCMatchedPassedEvents_TThypo_TTsignal = 0;
                    
                  int fracNoNeutrino_SMttHypo = 0;
                  int fracNoNeutrino_TThypo = 0;
                  int fracNoNeutrino_SThypo = 0;
                  int fracNoNeutrino_SMttHypo_signalMatched = 0;
                  int fracNoNeutrino_TThypo_signalMatched = 0;
                  int fracNoNeutrino_SThypo_signalMatched = 0;
	              
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
		              
		                    NumberOfEvents++;
		              
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

                        bool SMTT_Matched = false;
                        bool SignalTT_Matched = false;
                        bool SignalST_Matched = false;
                        int SMTT_TopbJets = 0;
                        int SMTT_WJets = 0;
                        int Signal_TopbJet = 0;
                        int Signal_TopFCNCJet = 0;
                        int Signal_HJets = 0;


                        for(int i_Jet = 0; i_Jet < NumberOfJets; i_Jet++)
                        {
                        
                            if(dataSetName.find("NP")!=string::npos && dataSetName.find("hct")!=string::npos)
                            {
                                if( (MotherpdgID[i_Jet] == 25) && (fabs(pdgID[i_Jet]) == 5)) Signal_HJets++;
                                else if( (MotherpdgID[i_Jet] == 6) && (fabs(pdgID[i_Jet]) == 5)) Signal_TopbJet++;
                                else if( (MotherpdgID[i_Jet] == 6) && (fabs(pdgID[i_Jet]) == 4)) Signal_TopFCNCJet++;
                            }
                            else if(dataSetName.find("NP")!=string::npos && dataSetName.find("hut")!=string::npos)
                            {
                                if( (MotherpdgID[i_Jet] == 25) && (fabs(pdgID[i_Jet]) == 5)) Signal_HJets++;
                                else if( (MotherpdgID[i_Jet] == 6) && (fabs(pdgID[i_Jet]) == 5)) Signal_TopbJet++;
                                else if( (MotherpdgID[i_Jet] == 6) && (fabs(pdgID[i_Jet]) == 2)) Signal_TopFCNCJet++;
                            }
                            else if(dataSetName.find("TTJets")!=string::npos)
                            {
                                if( (MotherpdgID[i_Jet] == 24) ) SMTT_WJets++;
                                else if( (MotherpdgID[i_Jet] == 6) && (fabs(pdgID[i_Jet]) == 5)) SMTT_TopbJets++;
                            }
                            

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
                        
                        if(SMTT_TopbJets == 2 && SMTT_WJets == 2) SMTT_Matched = true;
                        if(Signal_TopbJet == 1 && Signal_TopFCNCJet == 1 && Signal_HJets == 2) SignalTT_Matched = true;           
                        if(Signal_TopbJet == 1 && Signal_TopFCNCJet == 0 && Signal_HJets == 2) SignalST_Matched = true;           
                        
                        LeptonPt.push_back(lepPt);
                        LeptonEta.push_back(lepEta);
                        LeptonPhi.push_back(lepPhi);
                        LeptonE.push_back(lepE);


                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPTOPLEPHAD (selection: at least 2 -jets and at least 2 non-b jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "SMttHypo" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=2 && NonBJetPt.size()>=2)
                            {

                                  //cout <<  "Passed Event: " << j << ", Number of Bjets: " << BJetPt.size() << ", Number of NonBjets: " << NonBJetPt.size() << "Number of leptons" << LeptonPt.size() << endl;
                                  
	                                NSelectionPassedEvents_SMttHypo++;
	                                if(SMTT_Matched) nMCMatchedPassedEvents_SMttHypo_TTbackground++;
	                                if(SignalTT_Matched) nMCMatchedPassedEvents_SMttHypo_TTsignal++;
	                                if(SignalST_Matched) nMCMatchedPassedEvents_SMttHypo_STsignal++;

                                      kf_SMttHypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_SMttHypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_SMttHypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_SMttHypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_SMttHypo->SetMet(met_px,met_py);


                                      kf_SMttHypo->Run();
                                      int nPerm_SMttHypo = kf_SMttHypo->GetNPerm();

                                      double LowestDisc_SMttHypo = kf_SMttHypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_SMttHypo = -1;
                                      int IndexBJETHAD_SMttHypo = -1;
                                      int IndexNonBJET1HAD_SMttHypo = -1;
                                      int IndexNonBJET2HAD_SMttHypo = -1;

                                      for(int ip=0;ip<nPerm_SMttHypo;ip++)
                                      {
                                           double d = kf_SMttHypo->GetDisc(ip);
//                                           if(d>10E+9)  continue;
//                                           if(LowestDisc_SMttHypo != d) continue;

                                            if(d>10E+9) fracNoNeutrino_SMttHypo++;

		                                        IndexBJETLEP_SMttHypo = kf_SMttHypo->GetIndex(BJETLEP_TOPTOPLEPHAD,ip);
		                                        IndexBJETHAD_SMttHypo = kf_SMttHypo->GetIndex(BJETHAD_TOPTOPLEPHAD,ip);
		                                        IndexNonBJET1HAD_SMttHypo = kf_SMttHypo->GetIndex(NONBJET1_TOPTOPLEPHAD,ip);
		                                        IndexNonBJET2HAD_SMttHypo = kf_SMttHypo->GetIndex(NONBJET2_TOPTOPLEPHAD,ip);


                                          //Seeing if the kinematic fit can be matched to the ST-signal
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMttHypo]] == 25) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMttHypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMttHypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SMttHypo]]) == 5) )
                                              {
	                                                if(LowestDisc_SMttHypo == d) NMCIdentifiedEvents_SMttHypo_STsignal++;
                                                  if(d>10E+9) fracNoNeutrino_SMttHypo_signalMatched++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMttHypo]] == 25) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMttHypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMttHypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SMttHypo]]) == 5) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETHAD_SMttHypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETHAD_SMttHypo]]) != 5) )
                                              {
	                                                if(LowestDisc_SMttHypo == d) NMCIdentifiedEvents_SMttHypo_TTsignal++;
                                                  if(d>10E+9) fracNoNeutrino_SMttHypo_signalMatched++;
	                                            }
                                              if( (MotherpdgID[MapIndex_NonBindex[IndexNonBJET1HAD_SMttHypo]] == 24) && (MotherpdgID[MapIndex_NonBindex[IndexNonBJET2HAD_SMttHypo]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SMttHypo]]) == 6) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETHAD_SMttHypo]]) == 6) )
                                              {
	                                                if(LowestDisc_SMttHypo == d) NMCIdentifiedEvents_SMttHypo_TTbackground++;
                                                  if(d>10E+9) fracNoNeutrino_SMttHypo_signalMatched++;
	                                            }
                                          }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHAD
                        }
                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPTOPLEPHbb (selection: at least 3 b-jets and at least 1 non-b jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "TThypo" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=3 && NonBJetPt.size()>=1)
                            {
	                                NSelectionPassedEvents_TThypo++;
	                                if(SMTT_Matched) nMCMatchedPassedEvents_SMttHypo_TTbackground++;
	                                if(SignalTT_Matched) nMCMatchedPassedEvents_SMttHypo_TTsignal++;
	                                if(SignalST_Matched) nMCMatchedPassedEvents_SMttHypo_STsignal++;

                                      kf_TThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_TThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_TThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_TThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_TThypo->SetMet(met_px,met_py);


                                      kf_TThypo->Run();
                                      int nPerm_TThypo = kf_TThypo->GetNPerm();

                                      double LowestDisc_TThypo = kf_TThypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_TThypo = -1;
                                      int IndexBJETHAD_TThypo = -1;
                                      int IndexNonBJET1HAD_TThypo = -1;
                                      int IndexNonBJET2HAD_TThypo = -1;

                                      for(int ip=0;ip<nPerm_TThypo;ip++)
                                      {
                                           double d = kf_TThypo->GetDisc(ip);
//                                           if(d>10E+9)  continue;
//                                           if(LowestDisc_TThypo != d) continue;

                                            if(d>10E+9) fracNoNeutrino_TThypo++;

		                                        IndexBJETLEP_TThypo = kf_TThypo->GetIndex(BJETLEP_TOPTOPLEPHBB,ip);
		                                        IndexBJETHAD_TThypo = kf_TThypo->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip);
		                                        IndexNonBJET1HAD_TThypo = kf_TThypo->GetIndex(BJET1_TOPTOPLEPHBB,ip);
		                                        IndexNonBJET2HAD_TThypo = kf_TThypo->GetIndex(BJET2_TOPTOPLEPHBB,ip);

                                          //Seeing if the kinematic fit can be matched to the ST-signal
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TThypo]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TThypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TThypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_TThypo]]) == 5) )
                                              {
	                                                if(LowestDisc_TThypo == d) NMCIdentifiedEvents_TThypo_STsignal++;
                                                  if(d>10E+9) fracNoNeutrino_TThypo_signalMatched++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TThypo]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TThypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TThypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_TThypo]]) == 5) && (fabs(MotherpdgID[MapIndex_NonBindex[IndexBJETHAD_TThypo]]) == 6) && (fabs(pdgID[MapIndex_NonBindex[IndexBJETHAD_TThypo]]) != 5) )
                                              {
	                                                if(LowestDisc_TThypo == d) NMCIdentifiedEvents_TThypo_TTsignal++;
                                                  if(d>10E+9) fracNoNeutrino_TThypo_signalMatched++;
	                                            }
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_TThypo]] == 24) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_TThypo]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_TThypo]]) == 6) && (fabs(MotherpdgID[MapIndex_NonBindex[IndexBJETHAD_TThypo]]) == 6) )
                                              {
	                                                if(LowestDisc_TThypo == d) NMCIdentifiedEvents_TThypo_TTbackground++;
                                                  if(d>10E+9) fracNoNeutrino_TThypo_signalMatched++;
	                                            }
                                          }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHbb
                        }
                        /////////////////////////////////////////////////////////////////////////////
                        // Section for the kinematic fit using the TOPHLEPbb (selection: at least 3 b-jets)
                        /////////////////////////////////////////////////////////////////////////////
                        if(KinFitMethod ==   "SThypo" || KinFitMethod == "ALL")
                        {
                            if(BJetPt.size()>=3)
                            {
	                                NSelectionPassedEvents_SThypo++;
	                                if(SMTT_Matched) nMCMatchedPassedEvents_SMttHypo_TTbackground++;
	                                if(SignalST_Matched) nMCMatchedPassedEvents_SMttHypo_TTsignal++;
	                                if(SignalST_Matched) nMCMatchedPassedEvents_SMttHypo_STsignal++;

                                      kf_SThypo->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
                                      kf_SThypo->SetNonBJet(NonBJetPt,NonBJetEta,NonBJetPhi,NonBJetE);
                                      if(channel == "_El") kf_SThypo->SetElectron(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      else if(channel == "_Mu") kf_SThypo->SetMuon(LeptonPt,LeptonEta,LeptonPhi,LeptonE);
                                      kf_SThypo->SetMet(met_px,met_py);


                                      kf_SThypo->Run();
                                      int nPerm_SThypo = kf_SThypo->GetNPerm();

                                      double LowestDisc_SThypo = kf_SThypo->GetDisc(); //The minimum of likelihood == the best jet-combination
                                      int IndexBJETLEP_SThypo = -1;
                                      int IndexNonBJET1HAD_SThypo = -1;
                                      int IndexNonBJET2HAD_SThypo = -1;

                                      for(int ip=0;ip<nPerm_SThypo;ip++)
                                      {
                                           double d = kf_SThypo->GetDisc(ip);
//                                           if(d>10E+9)  continue;
//                                           if(LowestDisc_SThypo != d) continue;

                                            if(d>10E+9) fracNoNeutrino_SThypo++;

		                                        IndexBJETLEP_SThypo = kf_SThypo->GetIndex(BJETLEP_TOPHLEPBB,ip);
		                                        IndexNonBJET1HAD_SThypo = kf_SThypo->GetIndex(BJET1_TOPHLEPBB,ip);
		                                        IndexNonBJET2HAD_SThypo = kf_SThypo->GetIndex(BJET2_TOPHLEPBB,ip);
                                          //Seeing if the kinematic fit can be matched to the ST-signal
                                          if(SingleTop)
                                          {
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_SThypo]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_SThypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SThypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SThypo]]) == 5) )
                                              {
	                                                if(LowestDisc_SThypo == d) NMCIdentifiedEvents_SThypo_STsignal++;
                                                  if(d>10E+9) fracNoNeutrino_SThypo_signalMatched++;
                                              }
                                          }
                                          else
                                          {
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_SThypo]] == 25) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_SThypo]]==25) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SThypo]]) == 6) && (fabs(pdgID[MapIndex_Bindex[IndexBJETLEP_SThypo]]) == 5))
                                              {
      	                                          if(LowestDisc_SThypo == d) NMCIdentifiedEvents_SThypo_TTsignal++;
                                                  if(d>10E+9) fracNoNeutrino_SThypo_signalMatched++;
	                                            }
                                              if( (MotherpdgID[MapIndex_Bindex[IndexNonBJET1HAD_SThypo]] == 24) && (MotherpdgID[MapIndex_Bindex[IndexNonBJET2HAD_SThypo]]==24) && (fabs(MotherpdgID[MapIndex_Bindex[IndexBJETLEP_SThypo]]) == 6) )
                                              {
      	                                          if(LowestDisc_SThypo == d) NMCIdentifiedEvents_SThypo_TTbackground++;
                                                  if(d>10E+9) fracNoNeutrino_SThypo_signalMatched++;
	                                            }
                                          }
                                      }
                            }//End kinematic fit reconstruction for TOPTOPLEPHAD
                        }
		              }//for-loop events

    if(KinFitMethod == "SMttHypo" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPTOPLEPHAD ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHAD & " << 100*double(NSelectionPassedEvents_SMttHypo)/double(NumberOfEvents ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_STsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(nMCMatchedPassedEvents_SMttHypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTsignal)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(nMCMatchedPassedEvents_SMttHypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SMttHypo != 0) cout << 100*double(nMCMatchedPassedEvents_SMttHypo_TTbackground)/double(NSelectionPassedEvents_SMttHypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SMttHypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(nMCMatchedPassedEvents_SMttHypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SMttHypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }
        
        cout << "Number of combinations where no neutrino is reconstructed: " << fracNoNeutrino_SMttHypo << endl;
        cout << "Number of combinations where no neutrino is reconstructed and matched: " << fracNoNeutrino_SMttHypo_signalMatched << endl;
    }
    if(KinFitMethod == "TThypo" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPTOPLEPHBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPTOPLEPHBB & " << 100*double(NSelectionPassedEvents_TThypo)/double(NumberOfEvents ) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_STsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(nMCMatchedPassedEvents_TThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTsignal)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(nMCMatchedPassedEvents_TThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_TThypo != 0) cout << 100*double(nMCMatchedPassedEvents_TThypo_TTbackground)/double(NSelectionPassedEvents_TThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_TThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(nMCMatchedPassedEvents_TThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_TThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }    
        cout << "Number of combinations where no neutrino is reconstructed: " << fracNoNeutrino_TThypo << endl;
        cout << "Number of combinations where no neutrino is reconstructed and matched: " << fracNoNeutrino_TThypo_signalMatched << endl;
    }
    if(KinFitMethod == "SThypo" || KinFitMethod == "ALL")
    {
		    cout << "************ TOPHLEPBB ************" << endl;
		    cout << " TopKinFit-method$_{bTagWP}$ & $\epsilon_{selection}$ (\%) & $\epsilon_{MCmatched}$ (\%) & $\epsilon_{Method}$ (\%) & $\epsilon_{Total}$ (\%)" << endl;
        cout << " TOPHLEPBB & " << 100*double(NSelectionPassedEvents_SThypo)/double(NumberOfEvents) << " & ";//selection efficiency

        if(dataSetName.find("NP")!=string::npos)
        {
            if(SingleTop)
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_STsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_STsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(nMCMatchedPassedEvents_SThypo_STsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_STsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
            else
            {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTsignal)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTsignal != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(nMCMatchedPassedEvents_SThypo_TTsignal) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTsignal)/double(NumberOfEvents) << endl;//Total efficiency
            }
        }
        else
        {
                if(NSelectionPassedEvents_SThypo != 0) cout << 100*double(nMCMatchedPassedEvents_SThypo_TTbackground)/double(NSelectionPassedEvents_SThypo) << " & ";//reconstruction efficiency
                if(nMCMatchedPassedEvents_SThypo_TTbackground != 0) cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(nMCMatchedPassedEvents_SThypo_TTbackground) << " & ";//method efficiency
                cout << 100*double(NMCIdentifiedEvents_SThypo_TTbackground)/double(NumberOfEvents) << endl;//Total efficiency
        }
        cout << "Number of combinations where no neutrino is reconstructed: " << fracNoNeutrino_SThypo << endl;
        cout << "Number of combinations where no neutrino is reconstructed and matched: " << fracNoNeutrino_SThypo_signalMatched << endl;
    }
    cout << " " << endl;		              
		              
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







void MVAanalysis(bool doTraining, std::string MVAmethod, int skipEvents, std::string SignalName, std::string xmlNom, TString CraneenPath, std::string JetCombSelection)
{

  cout << "------------------------------------------------" << endl;
  cout << "-- Start of analysis macro ------------" << endl;
  cout << "------------------------------------------------" << endl;

  //////////////////////////////////////////////////////////////////////////////////
  /// Initializing all MVA instances for training and computing
  //////////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "Initializing MVA instances" << endl;
  
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
  if(debug) cout << "MVA initialized" << endl;

  //////////////////////////////////////////////////////////////////////////////////
  /// Initializing all KinFit instances for jet-combination selection
  //////////////////////////////////////////////////////////////////////////////////
  if(debug) cout << "Initializing KinFit instances" << endl;
   int nToys = 500;

  std::string pdfFileName_TThypo = "TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
  std::string pdfFileName_SThypo = "TopKinFit/test/GenAnalysis/TopHLepbb/pdf.root";

  KINFIT::kfit *kf_SThypo = new KINFIT::kfit();
  KINFIT::kfit *kf_TThypo = new KINFIT::kfit();

  
  if(JetCombSelection == "TThypo")
  {
      kf_TThypo->Init(TOPTOPLEPHBB);
      kf_TThypo->SetPDF("TopWMass",pdfFileName_TThypo.c_str(),"TopLepWM_Fit");
      kf_TThypo->SetPDF("TopMass",pdfFileName_TThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_TThypo->SetPDF("TopHadMass",pdfFileName_TThypo.c_str(),"TopHadRecM_Fit");
      kf_TThypo->SetPDF("MetPx",pdfFileName_TThypo.c_str(),"dMetPx_Gaus");
      kf_TThypo->SetPDF("MetPy",pdfFileName_TThypo.c_str(),"dMetPy_Gaus");
      kf_TThypo->SetPDF("BJetPx",pdfFileName_TThypo.c_str(),"dBJetPx_Fit");
      kf_TThypo->SetPDF("BJetPy",pdfFileName_TThypo.c_str(),"dBJetPy_Fit");
      kf_TThypo->SetPDF("BJetPz",pdfFileName_TThypo.c_str(),"dBJetPz_Fit");
      kf_TThypo->SetPDF("BJetE",pdfFileName_TThypo.c_str(),"dBJetE_Fit");
      kf_TThypo->SetPDF("ElecPx",pdfFileName_TThypo.c_str(),"dElecPx_Fit");
      kf_TThypo->SetPDF("ElecPy",pdfFileName_TThypo.c_str(),"dElecPy_Fit");
      kf_TThypo->SetPDF("ElecPz",pdfFileName_TThypo.c_str(),"dElecPz_Fit");
      kf_TThypo->SetPDF("ElecE",pdfFileName_TThypo.c_str(),"dElecE_Fit");
      kf_TThypo->SetPDF("MuonPx",pdfFileName_TThypo.c_str(),"dMuonPx_Fit");
      kf_TThypo->SetPDF("MuonPy",pdfFileName_TThypo.c_str(),"dMuonPy_Fit");
      kf_TThypo->SetPDF("MuonPz",pdfFileName_TThypo.c_str(),"dMuonPz_Fit");
      kf_TThypo->SetPDF("MuonE",pdfFileName_TThypo.c_str(),"dMuonE_Fit");
      kf_TThypo->SetPDF("NonBJetPx",pdfFileName_TThypo.c_str(),"dNonBJetPx_Fit");
      kf_TThypo->SetPDF("NonBJetPy",pdfFileName_TThypo.c_str(),"dNonBJetPy_Fit");
      kf_TThypo->SetPDF("NonBJetPz",pdfFileName_TThypo.c_str(),"dNonBJetPz_Fit");
      kf_TThypo->SetPDF("NonBJetE",pdfFileName_TThypo.c_str(),"dNonBJetE_Fit");
      kf_TThypo->SetNToy(nToys);
  }
  else if (JetCombSelection ==   "SThypo")
  {
      kf_SThypo->Init(TOPHLEPBB);
      kf_SThypo->SetPDF("TopWMass",pdfFileName_SThypo.c_str(),"TopLepWM_Fit");
      kf_SThypo->SetPDF("TopMass",pdfFileName_SThypo.c_str(),"TopLepRecM_Fit");
      kf_SThypo->SetPDF("HiggsMass",pdfFileName_TThypo.c_str(),"HiggsRecM_Fit");
      kf_SThypo->SetPDF("MetPx",pdfFileName_SThypo.c_str(),"dMetPx_Gaus");
      kf_SThypo->SetPDF("MetPy",pdfFileName_SThypo.c_str(),"dMetPy_Gaus");
      kf_SThypo->SetPDF("BJetPx",pdfFileName_SThypo.c_str(),"dBJetPx_Fit");
      kf_SThypo->SetPDF("BJetPy",pdfFileName_SThypo.c_str(),"dBJetPy_Fit");
      kf_SThypo->SetPDF("BJetPz",pdfFileName_SThypo.c_str(),"dBJetPz_Fit");
      kf_SThypo->SetPDF("BJetE",pdfFileName_SThypo.c_str(),"dBJetE_Fit");
      kf_SThypo->SetPDF("ElecPx",pdfFileName_SThypo.c_str(),"dElecPx_Fit");
      kf_SThypo->SetPDF("ElecPy",pdfFileName_SThypo.c_str(),"dElecPy_Fit");
      kf_SThypo->SetPDF("ElecPz",pdfFileName_SThypo.c_str(),"dElecPz_Fit");
      kf_SThypo->SetPDF("ElecE",pdfFileName_SThypo.c_str(),"dElecE_Fit");
      kf_SThypo->SetPDF("MuonPx",pdfFileName_SThypo.c_str(),"dMuonPx_Fit");
      kf_SThypo->SetPDF("MuonPy",pdfFileName_SThypo.c_str(),"dMuonPy_Fit");
      kf_SThypo->SetPDF("MuonPz",pdfFileName_SThypo.c_str(),"dMuonPz_Fit");
      kf_SThypo->SetPDF("MuonE",pdfFileName_SThypo.c_str(),"dMuonE_Fit");
      kf_SThypo->SetNToy(nToys);
  }
  if(debug) cout << "KinFit initialized" << endl;



  
  	/////////////////////////////////////////////////////////////
  	// Loading datasets
  	/////////////////////////////////////////////////////////////
    if(debug) cout << "Loading datasets" << endl;

  	const char *xmlfile = xmlNom.c_str();
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets;
  	treeLoader.LoadDatasets (datasets, xmlfile);
    MultiSamplePlot* MSP_MVA = new MultiSamplePlot(datasets, MVAmethod.c_str() , 40, -1., 1., MVAmethod.c_str()); 
  
  	string dataSetName, filepath;
  	int nEntries;
  	float ScaleFactor, NormFactor;
    if(debug) cout << "Datasets loaded" << endl;

  
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


