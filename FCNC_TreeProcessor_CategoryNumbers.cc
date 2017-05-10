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
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"



using namespace std;


/// Normal Plots (TH1F* and TH2F*)
map<string,double> SystScaleFactor;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;
map<string,MultiSamplePlot*> MSPlot_nPV;



// functions prototype
string intToStr (int number);
void MakeNPV_Distributions(string channel, string date, bool debug);
void MakeTotalSystErrorBand_Distributions(string  outfilename, vector< string > systematics, vector <string> datasetNames, vector<string> NominalVariableNames, string outputFile);

inline bool FileExists (const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{

    if(argc < 6)
    {
        cerr << "INVALID number of arguments. The necessary arguments are: " << endl;
        cout << "    string channel            = argv[1];" << endl;
        cout << "    string date            = argv[2];" << endl;
        cout << "    bool PVreweighing = strtol(argv[3], NULL,10);" << endl;
        cout << "    bool doJESSys  = strtol(argv4], NULL,10);" << endl;
        cout << "    bool doJERSys  = strtol(argv5], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[6], NULL,10);" << endl;

        return 1;
    }


    string channel            = argv[1];
    string date            = argv[2];
    bool PVreweighing = strtol(argv[3], NULL,10);
    bool doJESSys  = strtol(argv[4], NULL,10);
    bool doJERSys  = strtol(argv[5], NULL,10);
    bool debug         =strtol(argv[6], NULL,10);

    bool split_ttbar = true;   
    

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
    if(doJESSys) WhatSysts.push_back("JESPlus");
    if(doJESSys) WhatSysts.push_back("JESMinus");
    if(doJERSys) WhatSysts.push_back("JERPlus");
    if(doJERSys) WhatSysts.push_back("JERMinus");
    WhatSysts.push_back("");   

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


    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();

    cout << " ... Making the TreeProcessor .xml files " << endl;
    system("python scripts/MakeXMLforTreeProcessor.py");

    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    else if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    else if(channel == "_All")
    {
        xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
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
    Dataset* ttbar_bb = new Dataset();
    Dataset* ttbar_cc = new Dataset();
    Dataset* ttbar_ll = new Dataset();

    
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
            
            ttbar_bb->SetName("TTJets_bb");
            ttbar_cc->SetName("TTJets_cc");
            ttbar_ll->SetName("TTJets_ll");

            ttbar_bb->SetTitle("tt+bb");
            ttbar_cc->SetTitle("tt+cc");
            ttbar_ll->SetTitle("tt+lf");

            ttbar_bb->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            ttbar_cc->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            ttbar_ll->SetEquivalentLuminosity(datasets[d]->EquivalentLumi());
            
            ttbar_bb->SetColor(kAzure+3);
            ttbar_cc->SetColor(kAzure+5);
            ttbar_ll->SetColor(kAzure+7);
            
            datasets_splittedTTbar.pop_back();
            datasets_splittedTTbar.push_back(ttbar_bb);
            datasets_splittedTTbar.push_back(ttbar_cc);
            datasets_splittedTTbar.push_back(ttbar_ll);

            cout << " - split TTBar dataset into ..."  << ttbar_ll->Name() << ", " << ttbar_cc->Name() << " and " << ttbar_ll->Name()  << endl;

        }     
    }

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
    for(int iSyst = 0; iSyst<WhatSysts.size();iSyst++)
    {
        MSPlot[("CategoryRates"+WhatSysts[iSyst]).c_str()] = new MultiSamplePlot(datasets_splittedTTbar, ("CategoryRates"+WhatSysts[iSyst]).c_str(), 5, -0.5, 4.5, "BDT output","Events", "");
    }
  
 


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

        if(dataSetName.find("Private") != string::npos) continue;//Do not read out on private signal samples
        
        for(int JecCounter = WhatSysts_noJECs.size(); JecCounter < WhatSysts.size(); JecCounter++)
        {
            string postfix = "";
            if(isData && WhatSysts[JecCounter] != "") continue;
            if(!isData) postfix = WhatSysts[JecCounter];
	      

		        cout<<"Dataset:  :"<<(dataSetName+WhatSysts[JecCounter]).c_str()<<endl;
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
                pathPlot += "Output_NPV.root";            

                if(!FileExists(pathPlot))
                {
                    MakeNPV_Distributions(channel, date, debug);
                }
                
                W_nPV = reweight::LumiReWeighting( pathPlot.c_str(), pathPlot.c_str(), ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_Data");
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
            Double_t W_TopPtReweighing;

            Int_t nvtx;
            Int_t npu;
            Int_t genTTX;


            //variable for jets 
            Int_t nJets;
	          Int_t nJets_CSVM; 
	          
            // Weights
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Plus",&W_fleptonSF_Plus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF_Minus",&W_fleptonSF_Minus); //Contains, if muon, the  isoSF, idSF & trigSF
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Minus",&W_puSF_Minus);
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF_Plus",&W_puSF_Plus);
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
            ttree[(dataSetName).c_str()]->SetBranchAddress("W_TopPtReweighing",&W_TopPtReweighing);  

            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_npu",&npu);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_genTTX",&genTTX);

            // jets
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&nJets);
            ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&nJets_CSVM);
           
            double nloSF = 1.;
            int nPos = 0; 
            int nNeg = 0;
            if(isAMC && !isData)
            {
                for (int k = 0; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);
                    if( W_nloWeight > 0) nPos++;
                    else if( W_nloWeight < 0) nNeg ++;
                }
                if(nPos + nNeg != 0) nloSF *=((double) (nPos + nNeg))/ ((double) (nPos - nNeg));
            }		

            Double_t average_TopPtWeight = 0.;
            if(dataSetName.find("TTJets") != string::npos)
            {
                int nEventsPassed = 0;
                for (int k = 0; k<nEntries; k++)
                {
                    ttree[dataSetName.c_str()]->GetEntry(k);
		                double TopPtReweighing_Up = 1+ 2*(1-W_TopPtReweighing);

                    average_TopPtWeight = average_TopPtWeight + W_TopPtReweighing;
                    nEventsPassed++;
                }
                if(nEventsPassed != 0)
                {
                    average_TopPtWeight = average_TopPtWeight/nEventsPassed;
                }
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
                double ScaleFactor = 1.;
                double W_puSF_applied = 1.;
			          if(!isData)
			          {
			              if(!PVreweighing) W_puSF_applied = W_puSF;
			              else
			              {
			                  W_puSF_applied = W_nPV.ITweight( (int)nvtx );
			              }

                    //Nominal scale factor -- scale factors for systematic shifts are calculated below
                    ScaleFactor *= W_puSF_applied;
                    ScaleFactor *= W_fleptonSF;
                    ScaleFactor *= W_btagWeight_shape;
                    ScaleFactor *= nloSF * W_nloWeight;
                }
                else ScaleFactor = 1.;    
		          
      	        //***********************************************FILLING PLOTS**********************************************
                bool ScalePlots = true;
                if(isData) ScalePlots = false;

                if(filepath.find("JESMinus") == string::npos && filepath.find("JESPlus") == string::npos  && filepath.find("JERMinus") == string::npos && filepath.find("JERPlus") == string::npos)
                {
                    for(int iSyst_ = 0; iSyst_ < WhatSysts_noJECs.size(); iSyst_++)
                    {
                    
                        //-----------------------------------------------------------------------------------------------------------
                        // Calculate Scale factors
                        //-----------------------------------------------------------------------------------------------------------
                        SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] = 1.;
                    
			                  if(!isData)
			                  {
                            if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hf;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_lfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_lfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_hfstats2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_hfstats2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr1Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr1;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Plus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_up_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "iterativefit_cferr2Minus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape_down_cferr2;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "pileupMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Plus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "leptonMinus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF_Minus;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtPlus")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                                if(dataSetName.find("TTJets") != string::npos) SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_TopPtReweighing/average_TopPtWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "TopPtMinus")//Apply no TopPt reweighing
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
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
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoBTagSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoPUSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoLepSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= nloSF * W_nloWeight;
                            }
                            else if(WhatSysts_noJECs[iSyst_] == "NoNLOSF")
                            {
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_puSF_applied;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_fleptonSF;
                                SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] *= W_btagWeight_shape;
                            }
                        }//if(!isData)
                        else SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()] = 1.;



                        //-----------------------------------------------------------------------------------------------------------
                        // Fill Plots
                        //-----------------------------------------------------------------------------------------------------------

                        if(nJets_CSVM == 2 && nJets == 3) MSPlot[("CategoryRates"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(0, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        if(nJets_CSVM == 2 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(1, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        if(nJets_CSVM == 3 && nJets == 3) MSPlot[("CategoryRates"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(2, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        if(nJets_CSVM == 3 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(3, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                        if(nJets_CSVM == 4 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts_noJECs[iSyst_]).c_str()]->Fill(4, Sample, ScalePlots, Luminosity * SystScaleFactor[WhatSysts_noJECs[iSyst_].c_str()]);
                    }
               }
               if(filepath.find("JESMinus") != string::npos || filepath.find("JESPlus") != string::npos  || filepath.find("JERMinus") != string::npos || filepath.find("JERPlus") != string::npos || isData || WhatSysts[JecCounter] == "")
               {
                        if(nJets_CSVM == 2 && nJets == 3) MSPlot[("CategoryRates"+WhatSysts[JecCounter]).c_str()]->Fill(0, Sample, ScalePlots, Luminosity * ScaleFactor);
                        if(nJets_CSVM == 2 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts[JecCounter]).c_str()]->Fill(1, Sample, ScalePlots, Luminosity * ScaleFactor);
                        if(nJets_CSVM == 3 && nJets == 3) MSPlot[("CategoryRates"+WhatSysts[JecCounter]).c_str()]->Fill(2, Sample, ScalePlots, Luminosity * ScaleFactor);
                        if(nJets_CSVM == 3 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts[JecCounter]).c_str()]->Fill(3, Sample, ScalePlots, Luminosity * ScaleFactor);
                        if(nJets_CSVM == 4 && nJets >= 4) MSPlot[("CategoryRates"+WhatSysts[JecCounter]).c_str()]->Fill(4, Sample, ScalePlots, Luminosity * ScaleFactor);
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
  cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

  string outfilename = pathPNG+"/OutputCategory.root";

  TFile *outfile = new TFile(outfilename.c_str(),"recreate");
  outfile->cd();

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
      cout << "Drawing MSP: " << name << endl;
      temp->showNumberEntries(false);

      if(name.find("Minus") == string::npos && name.find("Plus")== string::npos 
        && name.find("noSF") == string::npos && name.find("OnlyTopPtSF") == string::npos && name.find("OnlyBTagSF") == string::npos && 
        name.find("OnlyPUSF") == string::npos && name.find("OnlyLepSF") == string::npos && name.find("OnlyNLOSF") == string::npos && 
        name.find("NoTopPtSF") == string::npos && name.find("NoBTagSF") == string::npos && name.find("NoPUSF") == string::npos && 
        name.find("NoLepSF") == string::npos && name.find("NoNLOSF") == string::npos)//Do not save the pictures of the systematics
      {
          NominalVariableNames.push_back(name);
      }
      temp->Write(outfile, name, false,pathPNG, "png");
	}

  outfile->Write("kOverwrite");
  outfile->Close();
  
  cout << "  - Making total systematic bands " << endl;
  string errorbandfile = (pathPNG+"/Systematics_BareHistos.root");
  WhatSysts.pop_back();//Delete the last entry (which should be "") for the systematics plotting
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


      if(name.find("CategoryRates") != string::npos)
      {
          vector<string> label;
          label.push_back("b2j3");
          label.push_back("b2j4");
          label.push_back("b3j3");
          label.push_back("b3j4");
          label.push_back("b4j4");
          temp->setBins(label);

          temp->showNumberEntries(false);
          temp->Draw("MyMSP_"+name, 1, true, true, true, 1);
          temp->Write(outfile_errorbands, name, true,pathPNG, "eps");

      }
      else
      {
           	
          if (debug)
          {
              cout << "Saving the MSP" << endl;
              cout << " and it->first is " << name << endl;
              cout << " Luminosity is " << Luminosity << endl;
          }
          cout << "Drawing MSP: " << name << endl;
          temp->showNumberEntries(false);
          temp->Draw("MyMSP_"+name, 1, true, true, true, 1);
          temp->Write(outfile_errorbands, name, true,pathPNG, "eps");//You can only call 1 format for saving the plots. The second time you want to draw, the THStacks are empty, because the object has been written to the root-file.
    //      temp->Write(outfile_errorbands, name, true,pathPNG, "png");
    //      temp->Write(outfile_errorbands, name, true,pathPNG, "pdf");//.pdf files are corrputed due to the #backslash symbol defined in MultiSamplePlot.cc
     }
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

void MakeNPV_Distributions( string channel, string date, bool debug)
{

    cout << ".. ..Making nPV_unw distributions for all samples.. .." << endl;


    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
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

    MSPlot_nPV["NPV_unw"] = new MultiSamplePlot(datasets, "NPV_unw", 51, -0.5, 50.5, "Number of PV","Events", ""); 

  
 

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

    //Define rate uncertainties
    Double_t LumiUncPlus = 0.026;
    Double_t LumiUncMinus = 0.026;
    Double_t XSecTTJetPlus = 0.055;
    Double_t XSecTTJetMinus = 0.055;
    Double_t XSecOtherPlus = 0.1;
    Double_t XSecOtherMinus = 0.1;

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
//cout << "h_tmp->GetBinContent(h_tmp->GetNbinsX()): " << h_tmp->GetBinContent(h_tmp->GetNbinsX()) << endl;
/*
            //making sure that the overflow is transferred to the last 'visible' bin; analogously for underflow...
            TH1F* h_tmp = (TH1F*) h_tmp_->Clone();
            int Nbins_ = h_tmp->GetNbinsX();
            h_tmp->SetBinContent(Nbins_,h_tmp->GetBinContent(Nbins_)+h_tmp->GetBinContent(Nbins_+1));
            h_tmp->SetBinContent(Nbins_+1,0);
            h_tmp->SetBinContent(1,h_tmp->GetBinContent(0)+h_tmp->GetBinContent(1));
            h_tmp->SetBinContent(0,0);
*/
            TH1F* h_tmp__scaleup = (TH1F*) h_tmp->Clone();//Make a new tmp which will be scaled according to the cross section uncertainty 
            TH1F* h_tmp__scaledown = (TH1F*) h_tmp->Clone();//Make a new tmp which will be scaled according to the cross section uncertainty



            if(datasetNames[iDataName].find("TTJets")!= string::npos)
            {
                h_tmp__scaleup->Scale(1+XSecTTJetPlus);
                h_tmp__scaledown->Scale(1-XSecTTJetMinus);
            }
            else
            {
                h_tmp__scaleup->Scale(1+XSecOtherPlus);
                h_tmp__scaledown->Scale(1-XSecOtherMinus);
            }

            if(iDataName == 0)
            {
                histo1D_nominal[NominalVariableNames[iVar].c_str()] = (TH1F*) h_tmp->Clone();
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()] = (TH1F*) h_tmp->Clone();//PlusPlus in the object name due to convention down below
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()] = (TH1F*) h_tmp->Clone();

                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncPlusPlus").c_str()] = (TH1F*) h_tmp__scaleup->Clone();//PlusPlus in the object name due to convention down below
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncMinusMinus").c_str()] = (TH1F*) h_tmp__scaledown->Clone();
            }
            else
            {
                histo1D_nominal[NominalVariableNames[iVar].c_str()]->Add(h_tmp);
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()]->Add(h_tmp);
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()]->Add(h_tmp);
                histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncPlusPlus").c_str()]->Add(h_tmp__scaleup);
                histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"XSecUncMinusMinus").c_str()]->Add(h_tmp__scaledown);
            }
        }

        //Scale the  Lumi Uncertainties
        histo1D_Up_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncPlusPlus").c_str()]->Scale(1+LumiUncPlus);
        histo1D_Down_SamplesAdded[(NominalVariableNames[iVar]+"LumiUncMinusMinus").c_str()]->Scale(1-LumiUncMinus);
        
            
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
//cout << "h_tmp->GetBinContent(h_tmp->GetNbinsX()): " << h_tmp->GetBinContent(h_tmp->GetNbinsX()) << endl;
                
                if(systematics[iSys].find("Plus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()] = (TH1F*) h_tmp->Clone();
                    else histo1D_Up_SamplesAdded[(varNameSys+"Plus").c_str()]->Add(h_tmp);
                }
                else if(systematics[iSys].find("Minus")!= string::npos)
                {
                    if(iDataName == 0) histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()] = (TH1F*) h_tmp->Clone();
                    else histo1D_Down_SamplesAdded[(varNameSys+"Minus").c_str()]->Add(h_tmp);
                }
            }
        }
            

        // Add the rate systematic undertainties to the systematics list, as well as the MC statistical uncertainty
        systematics.push_back("LumiUncPlus");
        systematics.push_back("LumiUncMinus");
        systematics.push_back("XSecUncPlus");
        systematics.push_back("XSecUncMinus");
       
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
            
//if(iBin == nBins) cout << "Content of last bin for Plus: " << bincontent_nominal + sqrt(bincontent_up_squared + bincontent_nominal) << endl;
//if(iBin == nBins) cout << "Content of last bin for Minus: " << bincontent_nominal  - sqrt(bincontent_down_squared + bincontent_nominal) << endl;
//if(iBin == nBins) cout << "Content of last bin for Nominal: " << bincontent_nominal << endl;

            histo1D_TotalUp[(NominalVariableNames[iVar]+"Plus").c_str()]->SetBinContent(iBin,bincontent_nominal + sqrt(bincontent_up_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
            histo1D_TotalDown[(NominalVariableNames[iVar]+"Minus").c_str()]->SetBinContent(iBin,bincontent_nominal  - sqrt(bincontent_down_squared + bincontent_nominal));//Also add once the statistical uncertainty on the MC
        }
        
        //Delete the last entries for the systematics that were added in MakeTotalSystErrorBand_Distributions
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
        systematics.pop_back();    
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


