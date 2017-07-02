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

//includes for MVA
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//includes for Kinematic fitting
#include "FCNCAnalysis/TopKinFit/kinfit.h"


using namespace std;
using namespace TopTree;
//using namespace KINFIT;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TTree*> ttree;
map<string,MultiSamplePlot*> MSPlot;






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


    string channel            = argv[1];
    string date            = argv[2];
    bool debug         =strtol(argv[3], NULL,10);
    
   
    
    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();




    ///////////////////////////////////////////////////////////////////////////////////
    // **************** Preparing samples ********************//
    ///////////////////////////////////////////////////////////////////////////////////
//    cout << " ... Making the TreeProcessor .xml files " << endl;
//    system("python scripts/MakeXMLforTreeProcessor.py");

    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    else if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
    else if(channel == "_All") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";//The xml file for the combined lepton channel doesn't really matter. Just make sure the correct data-lumi is in there
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
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  	string dataSetName, filepath;
    float Luminosity = 0;

    vector<string> CutsselecTable;
    MultiSamplePlot* cutFlowPlot = new MultiSamplePlot(datasets, "Cut Flow", 8, -0.5, 7.5, "Cut Number");

    if(channel == "_Mu")
    {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("At least 1 PV"));
        CutsselecTable.push_back(string("Event cleaning"));
        CutsselecTable.push_back(string("Trigger"));
        CutsselecTable.push_back(string("Exactly 1 Tight Isolated Muon"));
        CutsselecTable.push_back(string("Veto on Loose Electrons"));
        CutsselecTable.push_back(string("Veto on extra Loose Muons"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("At least 2 CSVv2M Jet"));
    }
    else if(channel == "_El")
    {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("At least 1 PV"));
        CutsselecTable.push_back(string("Event cleaning"));
        CutsselecTable.push_back(string("Trigger"));
        CutsselecTable.push_back(string("Exactly 1 Tight Electron"));
        CutsselecTable.push_back(string("Veto on Loose Muons"));
        CutsselecTable.push_back(string("Veto on extra Loose Electrons"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("At least 2 CSVv2M Jet"));
    }
    else if(channel == "_All")
    {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("At least 1 PV"));
        CutsselecTable.push_back(string("Event cleaning"));
        CutsselecTable.push_back(string("Trigger"));
        CutsselecTable.push_back(string("Exactly 1 Lepton"));
        CutsselecTable.push_back(string("Veto on Loose opposite flavour leptons"));
        CutsselecTable.push_back(string("Veto on extra (loose) same flavour leptons"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("At least 2 CSVv2M Jet"));
    }
    
    
    //***************************************************GETTING LUMI FROM DATA IN XML****************************************************
	  for (int d = 0; d < datasets.size(); d++)   //Loop through datasets  
	  {
		    dataSetName = datasets[d]->Name();
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        Luminosity = datasets[d]->EquivalentLumi();
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
		    filepath = TreePath+"/FCNC_1L3B__Run2_TopTree_Study_"+dataSetName + ".root";
		    if (debug)
		    {
		        cout<<"filepath: "<<filepath<<endl;
            cout <<"Equivalent luminosity of the dataset is: " << datasets[d]->EquivalentLumi() << endl;
		    }
	
		    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
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
		    string TTreename_info = "NtupleInfoTree";	
		    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_info.c_str()); //get ntuple creation information

        int nEntries;

		    nEntries = ttree[dataSetName.c_str()]->GetEntries();
		    cout<<"                 nEntries: "<<nEntries<<endl;//Should be the number of files that were merged, not the number of events
		
        Double_t cutstep[10]; //0: no cut, 1: PV cleaning, 2:event cleaning, 3: trigger, 4: lepton selection, 5: loose other-flavoured lepton veto, 6: veto extra loose leptons, 7: nb jets, 8: nb b-jets
        Int_t nCuts = 8; //REDEFINE if ncuts change
        Int_t nofPosWeights = 0;
        Int_t nofNegWeights = 0;

        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nofPosWeights",&nofPosWeights);  
	      ttree[(dataSetName).c_str()]->SetBranchAddress("I_nofNegWeights",&nofNegWeights);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nCuts",&nCuts); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("cutstep",&cutstep);


        double nloSF = 1.;
        int nPos = 0; 
        int nNeg = 0;
        if(isAMC)
        {
            for (int k = 0; k<nEntries; k++)
            {
                ttree[dataSetName.c_str()]->GetEntry(k);
                nPos = nPos + nofPosWeights;
                nNeg = nNeg + nofNegWeights;
            }
            nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
        }		

  	    //***********************************************RUNNING OVER EVENTS**********************************************
		    for (int j = 0; j<nEntries; j++)
		    {
		              
            if(debug)
            {
                cin.get();
                cout << " " << endl;
                cout << "------------NEW File in Merged sample: " << j << " --------------" << endl;
            }
			      ttree[dataSetName.c_str()]->GetEntry(j);
			      cuts.resize(nCuts);
			      
			      for(int cutnum = 0; cutnum < nCuts; cutnum++)
			      {
			          cuts[cutnum] = cuts[cutnum] + cutstep[cutnum];
			      }
		  }//for-loop events

		  for(int i = 0; i< cuts.size(); i++)
		  {
		      if(isData)
		      {
              cutFlowPlot->Fill(i, datasets[d], true, cuts[i]/datasets[d]->NormFactor());
              selecTable.Fill(d, i, cuts[i]/datasets[d]->NormFactor()/Luminosity);
          }
          else
          {
              cutFlowPlot->Fill(i, datasets[d], true, cuts[i]*Luminosity);
              selecTable.Fill(d, i, cuts[i]);
          }
      }

    }//for-loop datasets
               
    cout << "unloading datasets" << endl;
    cout << "Writing Cut Flow Plot" << endl;
    string pathPNG = "MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    pathPNG += "MSPlots";
    pathPNG += channel;
    pathPNG += "/";
    mkdir(pathPNG.c_str(),0777);
    pathPNG += date;
    pathPNG += "/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    TFile *outfile = new TFile((pathPNG+"Baseline_CutFlow.root").c_str(),"recreate");
    outfile->cd();



    cutFlowPlot->Draw("Cut Flow", 0, false, true, false, 100);
    cutFlowPlot->Write(outfile, "CutFlow", true, pathPNG, "png");
    //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    selecTable.TableCalculator(true, true, true, true, true,true,true);

    // Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies
    // (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTable.Write(pathPNG + "CutFlow_Table" + channel + ".tex", false, true, true, true, false, false, true);
  

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


