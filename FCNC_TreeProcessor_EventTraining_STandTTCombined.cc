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
map<string,MultiSamplePlot*> MSPlot_nPV;

int nTrainingVars = 0;
bool PrivateSampleTraining = false;





// functions prototype
string intToStr (int number);
void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug);
TMVA::Factory *TrainingFACTORY(string trName, TFile *outfile);

inline bool FileExists (const string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


int main(int argc, char *argv[])
{

    if(argc < 9)
    {
        cerr << "INVALID number of arguments. The necessary arguments are: " << endl;
        cout << "    int baseline_bjets             = strtol(argv[1], NULL,10);" << endl;
        cout << "    int baseline_jets                 = strtol(argv[2], NULL,10);" << endl;
        cout << "    string SignalSample            = argv[3];" << endl;
        cout << "    string channel            = argv[4];" << endl;
        cout << "    string date            = argv[5];" << endl;
        cout << "    bool PVreweighing = strtol(argv[6], NULL,10);" << endl;
        cout << "    bool debug         =strtol(argv[7], NULL,10);" << endl;
        cout << "    int khut         =strtod(argv[8], NULL,10);" << endl;
        cout << "    int khct         =strtod(argv[9], NULL,10);" << endl;

        return 1;
    }


    int baseline_bjets             = strtol(argv[1], NULL,10);
    int baseline_jets                 = strtol(argv[2], NULL,10);
    string SignalSample  = argv[3];//Valid arguments are: hut & hct, 2D
    string channel            = argv[4];
    string date            = argv[5];
    bool PVreweighing = strtol(argv[6], NULL,10);
    bool debug         =strtol(argv[7], NULL,10);
    int khut         =strtol(argv[8], NULL,10); //Divide this number by 100
    int khct         =strtol(argv[9], NULL,10);
    
   string coupling_hut = "khut0p" + intToStr(khut);
   if(khut < 10) coupling_hut = "khut0p0" + intToStr(khut);
   string coupling_hct = "khct0p" + intToStr(khct);
   if(khct < 10) coupling_hct = "khct0p0" + intToStr(khct);
    
    bool doInclusive = false;
    string category;
    if(baseline_bjets == 0 && baseline_jets == 0)
    {
        doInclusive = true;
        category = "Inclusive";
    }
    else
    {
        category = "b"+intToStr(baseline_bjets)+"j"+intToStr(baseline_jets);
    }
    string TrainingName = "";
    if(khut == 0 && khct == 0) TrainingName = "CombTraining_" + SignalSample + channel + "_" +  category;
    else TrainingName = "CombTraining_" + SignalSample + "_" + coupling_hut + "_" + coupling_hct + "_" + channel + "_" +  category;

    cout << "------------------------------------------------------------------------------------------------" << endl;
    cout << "Begin program" << endl;
    cout << " - Category: " << category << endl;
    cout << "------------------------------------------------------------------------------------------------" << endl;


    clock_t start = clock();

    ///////////////////////////////////////////////////////////////////////////////////////
    // ******* Defining the training factory **********************//
    ///////////////////////////////////////////////////////////////////////////////////////
	  TFile* outfile = TFile::Open(("weights/"+TrainingName+".root").c_str(),"RECREATE");
	  TMVA::Factory *factory = TrainingFACTORY(TrainingName,outfile);
    if(factory == 0)
    {
        cerr << "Training factory is empty" << endl;
        return 1;
    }

	  TRandom3 *r = new TRandom3(666);
	  std::vector<double> vars(nTrainingVars);






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

  	const char *xmlfile = xmlNom.c_str();
  	cout << "used config file: " << xmlfile << endl;
  
    //***************************************************LOADING DATASETS****************************************************
  	TTreeLoader treeLoader;
  	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
  	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;
  	string dataSetName, filepath;
    float Luminosity = 0;

  
 
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
            pathPlot += category;
            pathPlot += "/";
            mkdir(pathPlot.c_str(),0777);            
            pathPlot += "Output_NPV.root";            

            if(!FileExists(pathPlot))
            {
                MakeNPV_Distributions(baseline_jets, baseline_bjets, channel, date, debug);
            }
            
            W_nPV = reweight::LumiReWeighting( pathPlot.c_str(), pathPlot.c_str(), ("MultiSamplePlot_NPV_unw/NPV_unw_"+dataSetName).c_str(), "MultiSamplePlot_NPV_unw/NPV_unw_Data");    
        }

		    FileObj[dataSetName.c_str()] = new TFile((filepath).c_str(),"READ"); //create TFile for each dataset      
		                
		                
		    bool isAMC = false;
		    if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
		    {
		        if(debug) cout << "Data found" << endl;
		        continue;
	      }
        else if(dataSetName.find("NLO") != string::npos || dataSetName.find("nlo") !=string::npos || dataSetName.find("amc") !=string::npos) isAMC = true;

		    bool isSignal = false;
		    double SignalWeight_2D = 1.;
		    if(dataSetName.find("NP_overlay")!=string::npos)
		    {
		        if(SignalSample == "hut" && dataSetName.find("hut")!=string::npos) isSignal = true;
		        else if(SignalSample == "hct" && dataSetName.find("hct")!=string::npos) isSignal = true;
		        else if(SignalSample   == "2D")
		        {
		            isSignal = true;
		            if(dataSetName.find("hct")!=string::npos) SignalWeight_2D = pow(double(khct)/100,2);
		            else if(dataSetName.find("hut")!=string::npos) SignalWeight_2D = pow(double(khut)/100,2);
		        }
		        
		        if(!isSignal) continue;
		        
		        if(PrivateSampleTraining && dataSetName.find("Private")==string::npos) continue;//Don't train over non-private signal samples
	      }


  	    //***********************************************IMPORTING VARIABLES**********************************************
		    string TTreename = "ObjectVarsTree";	
		    string TTreename_info = "NtupleInfoTree";	
		    ttree[dataSetName.c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename.c_str()); //get ttree for each dataset
		    ttree[(dataSetName+TTreename_info).c_str()] = (TTree*)FileObj[dataSetName.c_str()]->Get(TTreename_info.c_str()); //get ntuple creation information

        int nEntries;

		    nEntries = ttree[dataSetName.c_str()]->GetEntries();
		    cout<<"                 nEntries: "<<nEntries<<endl;
		
        //Weights
        Double_t W_puSF;
        Double_t W_fleptonSF;
        Double_t W_btagWeight_shape;
        Double_t W_nloWeight;// for amc@nlo samples
      
        Int_t run_num;
        Int_t evt_num;
        Int_t lumi_num;
        Int_t nvtx;
        Int_t npu;

        //variable for  leptons
        Int_t LepCharge;
  
        //variable for jets 
	      Int_t nJets; 
	      Int_t nJets_CSVM; 
            Double_t incl_charge_jet[20];
	      
	      //JetIndices_correctJetComb
        Double_t MVA_TOPTOPLEPHAD = -999.;
        Double_t MVA_TOPTOPLEPHBB = -999.;
        Double_t MVA_TOPHLEPBB_hut = -999.;
        Double_t MVA_TOPHLEPBB_hct = -999.;
	          Int_t TOPTOPLEPHAD_JetIdx_LepTop;
	          Int_t TOPTOPLEPHAD_JetIdx_HadTop;
	          Int_t TOPTOPLEPHAD_JetIdx_W1;
	          Int_t TOPTOPLEPHAD_JetIdx_W2;
	          Int_t TOPTOPLEPHBB_JetIdx_LepTop;
	          Int_t TOPTOPLEPHBB_JetIdx_HadTop;
	          Int_t TOPTOPLEPHBB_JetIdx_H1;
	          Int_t TOPTOPLEPHBB_JetIdx_H2;
        //Variables for signal/background training
	      Double_t HiggsMass_TOPHLEPBB_hut;
	      Double_t HiggsMass_TOPHLEPBB_hct;
	      Double_t HiggsEta_TOPHLEPBB_hut;
	      Double_t HiggsEta_TOPHLEPBB_hct;
	      Double_t TopLepMass_TOPHLEPBB_hut;
	      Double_t TopLepMass_TOPHLEPBB_hct;
        Double_t TopLepPt_TOPHLEPBB_hut;
        Double_t TopLepPt_TOPHLEPBB_hct;
        Double_t TopLepEta_TOPHLEPBB_hut;
        Double_t TopLepEta_TOPHLEPBB_hct;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
        Double_t TopLepHiggsDr_TOPHLEPBB_hut;
        Double_t TopLepHiggsDr_TOPHLEPBB_hct;
        Double_t HiggsBJet1CSVv2_TOPHLEPBB_hut;
        Double_t HiggsBJet1CSVv2_TOPHLEPBB_hct;
        Double_t HiggsBJet2CSVv2_TOPHLEPBB_hut;
        Double_t HiggsBJet2CSVv2_TOPHLEPBB_hct;
        Double_t TopLepBJetCSVv2_TOPHLEPBB_hut;
        Double_t TopLepBJetCSVv2_TOPHLEPBB_hct;
        Double_t TopHadMass_TOPTOPLEPHAD;
        Double_t TopLepMass_TOPTOPLEPHAD;
        Double_t TopLepTopHadDr_TOPTOPLEPHAD;
        Double_t TopLepBJetCSVv2_TOPTOPLEPHAD;
        Double_t TopHadBJetCSVv2_TOPTOPLEPHAD;
        Double_t TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
        Double_t TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;
        Double_t HiggsMass_TOPTOPLEPHBB;
        Double_t TopLepMass_TOPTOPLEPHBB;
        Double_t HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
        Double_t TopLepHiggsDr_TOPTOPLEPHBB;
        Double_t HiggsBJet1CSVv2_TOPTOPLEPHBB;
        Double_t HiggsBJet2CSVv2_TOPTOPLEPHBB;
        Double_t TopLepBJetCSVv2_TOPTOPLEPHBB;
        Double_t TopHadNonBJetCSVv2_TOPTOPLEPHBB;

       
        // Weights
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_fleptonSF",&W_fleptonSF); //Contains, if muon, the  isoSF, idSF & trigSF
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_puSF",&W_puSF);
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_btagWeight_shape",&W_btagWeight_shape); 
        ttree[(dataSetName).c_str()]->SetBranchAddress("W_nloWeight",&W_nloWeight); 

        ttree[(dataSetName).c_str()]->SetBranchAddress("I_run_num",&run_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_evt_num",&evt_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_lumi_num",&lumi_num);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nvtx",&nvtx);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_npu",&npu);


        //SelectedLepton
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_LepCharge",&LepCharge);
        
        // jets
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets",&nJets);
        ttree[(dataSetName).c_str()]->SetBranchAddress("I_nJets_CSVM",&nJets_CSVM);
            ttree[(dataSetName).c_str()]->SetBranchAddress("incl_charge_jet",&incl_charge_jet);	    
       
        // Jet-indices associated to the jet-assignment in the bMVA method
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHAD",&MVA_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPTOPLEPHBB",&MVA_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hut",&MVA_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("MVA_TOPHLEPBB_hct",&MVA_TOPHLEPBB_hct);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_LepTop",&TOPTOPLEPHAD_JetIdx_LepTop);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_HadTop",&TOPTOPLEPHAD_JetIdx_HadTop);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_W1",&TOPTOPLEPHAD_JetIdx_W1);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHAD_JetIdx_W2",&TOPTOPLEPHAD_JetIdx_W2);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_LepTop",&TOPTOPLEPHBB_JetIdx_LepTop);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_HadTop",&TOPTOPLEPHBB_JetIdx_HadTop);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_H1",&TOPTOPLEPHBB_JetIdx_H1);
           ttree[(dataSetName).c_str()]->SetBranchAddress("I_TOPTOPLEPHBB_JetIdx_H2",&TOPTOPLEPHBB_JetIdx_H2);
        //Variables for signal/background training
	      ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPHLEPBB_hut",&HiggsMass_TOPHLEPBB_hut);
	      ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPHLEPBB_hct",&HiggsMass_TOPHLEPBB_hct);
	      ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsEta_TOPHLEPBB_hut",&HiggsEta_TOPHLEPBB_hut);
	      ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsEta_TOPHLEPBB_hct",&HiggsEta_TOPHLEPBB_hct);
	      ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPHLEPBB_hut",&TopLepMass_TOPHLEPBB_hut);
	      ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPHLEPBB_hct",&TopLepMass_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepPt_TOPHLEPBB_hut",&TopLepPt_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepPt_TOPHLEPBB_hct",&TopLepPt_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepEta_TOPHLEPBB_hut",&TopLepEta_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepEta_TOPHLEPBB_hct",&TopLepEta_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct",&HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPHLEPBB_hut",&TopLepHiggsDr_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPHLEPBB_hct",&TopLepHiggsDr_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB_hut",&HiggsBJet1CSVv2_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPHLEPBB_hct",&HiggsBJet1CSVv2_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB_hut",&HiggsBJet2CSVv2_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPHLEPBB_hct",&HiggsBJet2CSVv2_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB_hut",&TopLepBJetCSVv2_TOPHLEPBB_hut);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPHLEPBB_hct",&TopLepBJetCSVv2_TOPHLEPBB_hct);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadMass_TOPTOPLEPHAD",&TopHadMass_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPTOPLEPHAD",&TopLepMass_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepTopHadDr_TOPTOPLEPHAD",&TopLepTopHadDr_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHAD",&TopLepBJetCSVv2_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadBJetCSVv2_TOPTOPLEPHAD",&TopHadBJetCSVv2_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet1CSVv2_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",&TopHadWNonBJet2CSVv2_TOPTOPLEPHAD);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsMass_TOPTOPLEPHBB",&HiggsMass_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepMass_TOPTOPLEPHBB",&TopLepMass_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",&HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepHiggsDr_TOPTOPLEPHBB",&TopLepHiggsDr_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet1CSVv2_TOPTOPLEPHBB",&HiggsBJet1CSVv2_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("HiggsBJet2CSVv2_TOPTOPLEPHBB",&HiggsBJet2CSVv2_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopLepBJetCSVv2_TOPTOPLEPHBB",&TopLepBJetCSVv2_TOPTOPLEPHBB);
        ttree[(dataSetName).c_str()]->SetBranchAddress("TopHadNonBJetCSVv2_TOPTOPLEPHBB",&TopHadNonBJetCSVv2_TOPTOPLEPHBB);

        int nTrainingEntries = nEntries;
        nTrainingEntries = int(nEntries/2);

                  
        double nloSF = 1.;
        int nPos = 0; 
        int nNeg = 0;
        if(isAMC)
        {
            for (int k = 0; k<nTrainingEntries; k++)
            {
                ttree[dataSetName.c_str()]->GetEntry(k);
                if( W_nloWeight > 0) nPos++;
                else if( W_nloWeight < 0) nNeg ++;
            }
            nloSF *= ((double) (nPos - nNeg))/((double) (nPos + nNeg));
        }		


        
        
  	    //***********************************************RUNNING OVER EVENTS**********************************************
		    for (int j = 0; j<nTrainingEntries; j++)
		    {
		              
            if(debug)
            {
                cin.get();
                cout << " " << endl;
                cout << "------------NEW EVENT: " << j << " --------------" << endl;
            }
			      ttree[dataSetName.c_str()]->GetEntry(j);
		        if(!doInclusive)
		        {
		            if(nJets_CSVM != baseline_bjets)  continue;

		            if(baseline_jets == 3 && nJets != baseline_jets) continue;
		            else if(baseline_jets == 4 && nJets < baseline_jets) continue;
		        }

                //////////////////////////////////////
                //Applying the scale factors
                ///////////////////////////////////////
                float ScaleFactor = 1.; // event scale factor
                double W_puSF_applied = 1.;
			          if(!PVreweighing) W_puSF_applied = W_puSF;
			          else
			          {
			              W_puSF_applied = W_nPV.ITweight( (int)nvtx );
			          }

                //Safety triggers in case there are strange things happening in the event weights
                if(W_fleptonSF < 0 || W_btagWeight_shape < 0 || nloSF < 0 || Luminosity < 0 || W_puSF_applied < 0)
                {
                      cout << "----- Event " << j << " has a negative weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << "; Luminosity=" << Luminosity << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_fleptonSF != W_fleptonSF || W_btagWeight_shape != W_btagWeight_shape || nloSF != nloSF || W_puSF_applied != W_puSF_applied)
                {
                      cout << "----- Event " << j << " has a Nan weight. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
                      cout << "----- The event will be skipped....." << endl;
                      continue;
                }
                else if(W_fleptonSF >= 40 || W_btagWeight_shape >= 40 || nloSF >= 40 || W_puSF_applied >= 40)
                {
                      cout << "----- Event " << j << " has a weight larger than 20. Weights are: W_puSF=" << W_puSF_applied << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                      cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
//                      cout << "----- The event will be skipped....." << endl;
//                      continue;
                }


		           ScaleFactor = ScaleFactor * W_puSF_applied;
               ScaleFactor = ScaleFactor * W_fleptonSF;
  			       ScaleFactor = ScaleFactor * W_btagWeight_shape;
               ScaleFactor = ScaleFactor * nloSF;
               
               double weight = ScaleFactor * Luminosity * datasets[d]->NormFactor();

			      if(debug)
			      {
                cout << "----- Event " << j << "Weights are: W_puSF=" << W_puSF << "; W_fleptonSF=" << W_fleptonSF << "; W_btagWeight_shape=" << W_btagWeight_shape << "; nloSF=" << nloSF << endl;
                cout << "----- event number: " << evt_num << ", lumi_num: " << lumi_num << endl;
			          cout << "   SCALE FACTOR is: " << ScaleFactor << endl;
			      }
			      
            Double_t TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB = fabs(incl_charge_jet[TOPTOPLEPHBB_JetIdx_HadTop]+incl_charge_jet[TOPTOPLEPHBB_JetIdx_LepTop]);
            Double_t TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB = fabs(incl_charge_jet[TOPTOPLEPHBB_JetIdx_HadTop]+LepCharge);
            Double_t TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD = fabs(incl_charge_jet[TOPTOPLEPHAD_JetIdx_HadTop]+incl_charge_jet[TOPTOPLEPHAD_JetIdx_LepTop]);
            Double_t TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD = fabs(incl_charge_jet[TOPTOPLEPHAD_JetIdx_HadTop]+LepCharge);


	          if( HiggsMass_TOPHLEPBB_hut > 500. ) HiggsMass_TOPHLEPBB_hut = 500.;
	          if( TopLepMass_TOPHLEPBB_hut > 500. ) TopLepMass_TOPHLEPBB_hut = 500.;
	          if( TopLepPt_TOPHLEPBB_hut > 1000. ) TopLepPt_TOPHLEPBB_hut = 1000.;
	          if( HiggsMass_TOPHLEPBB_hct > 500. ) HiggsMass_TOPHLEPBB_hct = 500.;
	          if( TopLepMass_TOPHLEPBB_hct > 500. ) TopLepMass_TOPHLEPBB_hct = 500.;
	          if( TopLepPt_TOPHLEPBB_hct > 1000. ) TopLepPt_TOPHLEPBB_hct = 1000.;
	          if( TopLepMass_TOPTOPLEPHAD > 500.) TopLepMass_TOPTOPLEPHAD = 500.;
	          if( HiggsMass_TOPTOPLEPHBB > 500. ) HiggsMass_TOPTOPLEPHBB = 500.;
	          if( TopLepMass_TOPTOPLEPHBB > 500. ) TopLepMass_TOPTOPLEPHBB = 500.;
	          if( TopHadMass_TOPTOPLEPHAD > 1000. ) TopHadMass_TOPTOPLEPHAD = 1000.;
            if( HiggsBJet1CSVv2_TOPHLEPBB_hut < 0.) HiggsBJet1CSVv2_TOPHLEPBB_hut= 0.;
            if( HiggsBJet1CSVv2_TOPHLEPBB_hct < 0.) HiggsBJet1CSVv2_TOPHLEPBB_hct= 0.;
            if( HiggsBJet2CSVv2_TOPHLEPBB_hut < 0.) HiggsBJet2CSVv2_TOPHLEPBB_hut= 0.;
            if( HiggsBJet2CSVv2_TOPHLEPBB_hct < 0.) HiggsBJet2CSVv2_TOPHLEPBB_hct= 0.;
            if( TopLepBJetCSVv2_TOPHLEPBB_hut < 0.) TopLepBJetCSVv2_TOPHLEPBB_hut= 0.;
            if( TopLepBJetCSVv2_TOPHLEPBB_hct < 0.) TopLepBJetCSVv2_TOPHLEPBB_hct= 0.;
            if( TopLepBJetCSVv2_TOPTOPLEPHAD < 0.) TopLepBJetCSVv2_TOPTOPLEPHAD= 0.;
            if( TopHadBJetCSVv2_TOPTOPLEPHAD < 0.) TopHadBJetCSVv2_TOPTOPLEPHAD= 0.;
            if( TopHadWNonBJet1CSVv2_TOPTOPLEPHAD < 0.) TopHadWNonBJet1CSVv2_TOPTOPLEPHAD = 0.;
            if( TopHadWNonBJet2CSVv2_TOPTOPLEPHAD < 0.) TopHadWNonBJet2CSVv2_TOPTOPLEPHAD= 0.;
            if( HiggsBJet1CSVv2_TOPTOPLEPHBB < 0.) HiggsBJet1CSVv2_TOPTOPLEPHBB= 0.;
            if( HiggsBJet2CSVv2_TOPTOPLEPHBB < 0.) HiggsBJet2CSVv2_TOPTOPLEPHBB= 0.;
            if( TopLepBJetCSVv2_TOPTOPLEPHBB < 0.) TopLepBJetCSVv2_TOPTOPLEPHBB= 0.;
            if( TopHadNonBJetCSVv2_TOPTOPLEPHBB < 0.) TopHadNonBJetCSVv2_TOPTOPLEPHBB = 0.;



	         
            if(TrainingName.find("hut")!=string::npos && TrainingName.find("j3")!=string::npos)
            {
	               vars[0] = HiggsMass_TOPHLEPBB_hut;
	               vars[1] = MVA_TOPHLEPBB_hut;
	               vars[2] = LepCharge;
	               vars[3] = HiggsEta_TOPHLEPBB_hut;
	               vars[4] = TopLepMass_TOPHLEPBB_hut;
	               vars[5] = TopLepPt_TOPHLEPBB_hut;
	               vars[6] = TopLepEta_TOPHLEPBB_hut;
	               vars[7] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
	               vars[8] = TopLepHiggsDr_TOPHLEPBB_hut;
	               vars[9] = HiggsBJet1CSVv2_TOPHLEPBB_hut;
	               vars[10] = HiggsBJet2CSVv2_TOPHLEPBB_hut;
	               vars[11] = TopLepBJetCSVv2_TOPHLEPBB_hut;
	          }
            else if(TrainingName.find("hct")!=string::npos && TrainingName.find("j3")!=string::npos)
            {
	               vars[0] = HiggsMass_TOPHLEPBB_hct;
	               vars[1] = MVA_TOPHLEPBB_hct;
	               vars[2] = HiggsEta_TOPHLEPBB_hct;
	               vars[3] = TopLepMass_TOPHLEPBB_hct;
	               vars[4] = TopLepPt_TOPHLEPBB_hct;
	               vars[5] = TopLepEta_TOPHLEPBB_hct;
	               vars[6] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
	               vars[7] = TopLepHiggsDr_TOPHLEPBB_hct;
	               vars[8] = HiggsBJet1CSVv2_TOPHLEPBB_hct;
	               vars[9] = HiggsBJet2CSVv2_TOPHLEPBB_hct;
	               vars[10] = TopLepBJetCSVv2_TOPHLEPBB_hct;
	          }
            else if(TrainingName.find("hut")!=string::npos && TrainingName.find("j4")!=string::npos)
            {
	               vars[0] = HiggsMass_TOPHLEPBB_hut;
	               vars[1] = TopHadMass_TOPTOPLEPHAD;
	               vars[2] = MVA_TOPHLEPBB_hut;
	               vars[3] = MVA_TOPTOPLEPHAD;
	               vars[4] = LepCharge;
	               vars[5] = HiggsEta_TOPHLEPBB_hut;
//	               vars[6] = TopLepMass_TOPHLEPBB_hut;
//	               vars[7] = TopLepMass_TOPTOPLEPHAD;
//	               vars[8] = TopLepPt_TOPHLEPBB_hut;
//	               vars[9] = TopLepEta_TOPHLEPBB_hut;
	               vars[6] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hut;
//	               vars[11] = TopLepHiggsDr_TOPHLEPBB_hut;
//	               vars[12] = TopLepTopHadDr_TOPTOPLEPHAD;
	               vars[7] = HiggsBJet1CSVv2_TOPHLEPBB_hut;
	               vars[8] = HiggsBJet2CSVv2_TOPHLEPBB_hut;
	               vars[9] = TopLepBJetCSVv2_TOPHLEPBB_hut;
//	               vars[16] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	               vars[10] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	               vars[11] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	               vars[12] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

//	               vars[20] = HiggsMass_TOPTOPLEPHBB;
//	               vars[21] = MVA_TOPTOPLEPHBB;
//	               vars[22] = TopLepMass_TOPTOPLEPHBB;
//	               vars[23] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
//	               vars[24] = TopLepHiggsDr_TOPTOPLEPHBB;
//	               vars[25] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
//	               vars[26] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
//	               vars[27] = TopLepBJetCSVv2_TOPTOPLEPHBB;
//	               vars[28] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
/*
	               vars[29] = TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB;
	               vars[30] = TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB;
	               vars[31] = TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD;
	               vars[32] = TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD;

*/	          }
            else if(TrainingName.find("hct")!=string::npos && TrainingName.find("j4")!=string::npos)
            {
	               vars[0] = HiggsMass_TOPHLEPBB_hct;
	               vars[1] = TopHadMass_TOPTOPLEPHAD;
	               vars[2] = MVA_TOPHLEPBB_hct;
	               vars[3] = MVA_TOPTOPLEPHAD;
	               vars[4] = HiggsEta_TOPHLEPBB_hct;
//	               vars[5] = TopLepMass_TOPHLEPBB_hct;
//	               vars[6] = TopLepMass_TOPTOPLEPHAD;
//	               vars[7] = TopLepPt_TOPHLEPBB_hct;
//	               vars[8] = TopLepEta_TOPHLEPBB_hct;
	               vars[5] = HiggsBJet1HiggsBJet2Dr_TOPHLEPBB_hct;
//	               vars[10] = TopLepHiggsDr_TOPHLEPBB_hct;
//	               vars[11] = TopLepTopHadDr_TOPTOPLEPHAD;
	               vars[6] = HiggsBJet1CSVv2_TOPHLEPBB_hct;
	               vars[7] = HiggsBJet2CSVv2_TOPHLEPBB_hct;
	               vars[8] = TopLepBJetCSVv2_TOPHLEPBB_hct;
//	               vars[15] = TopLepBJetCSVv2_TOPTOPLEPHAD;
	               vars[9] = TopHadBJetCSVv2_TOPTOPLEPHAD;
	               vars[10] = TopHadWNonBJet1CSVv2_TOPTOPLEPHAD;
	               vars[11] = TopHadWNonBJet2CSVv2_TOPTOPLEPHAD;

//	               vars[19] = HiggsMass_TOPTOPLEPHBB;
//	               vars[20] = MVA_TOPTOPLEPHBB;
//	               vars[21] = TopLepMass_TOPTOPLEPHBB;
//	               vars[22] = HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB;
//	               vars[23] = TopLepHiggsDr_TOPTOPLEPHBB;
//	               vars[24] = HiggsBJet1CSVv2_TOPTOPLEPHBB;
//	               vars[25] = HiggsBJet2CSVv2_TOPTOPLEPHBB;
//	               vars[26] = TopLepBJetCSVv2_TOPTOPLEPHBB;
//	               vars[27] = TopHadNonBJetCSVv2_TOPTOPLEPHBB;
/*
	               vars[28] = TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB;
	               vars[29] = TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB;
	               vars[30] = TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD;
	               vars[31] = TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD;

*/	          }



            float rnd = r->Rndm();
	          if(isSignal)
            {
	              if( rnd < 0.5 )
	                factory->AddSignalTrainingEvent(vars,weight*SignalWeight_2D);
	              else
	                factory->AddSignalTestEvent(vars,weight*SignalWeight_2D);
			      }
			      else
			      {
	              if( rnd < 0.5 )
	                factory->AddBackgroundTrainingEvent(vars,weight);
	              else
	              {
	                factory->AddBackgroundTestEvent(vars,weight);
                  }
			      } 
			                
		  }//for-loop events
		              
    }//for-loop datasets
               

	  factory->PrepareTrainingAndTestTree("","",
					      "SplitMode=Random:NormMode=NumEvents:!V" );
	
	  factory->BookMethod(TMVA::Types::kBDT,"BDT",
			      "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:IgnoreNegWeightsInTraining" );
	
	  factory->TrainAllMethods();
	
	  factory->TestAllMethods();
	
	  factory->EvaluateAllMethods();

	  outfile->Write();
	  outfile->Close();


  

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

void MakeNPV_Distributions(int baseline_jets, int baseline_bjets, string channel, string date, bool debug)
{
    bool doInclusive = false;
    string category;
    if(baseline_bjets == 0 && baseline_jets == 0)
    {
        doInclusive = true;
        category = "Inclusive";
    }
    else
    {
        category = "b"+intToStr(baseline_bjets)+"j"+intToStr(baseline_jets);
    }    

    cout << ".. ..Making nPV_unw distributions for all samples.. .." << endl;


    string xmlNom;
    if(channel == "_El") xmlNom = "config/FullMcBkgdSamples_El_TreeProcessor.xml";
    if(channel == "_Mu") xmlNom = "config/FullMcBkgdSamples_Mu_TreeProcessor.xml";
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

    MSPlot_nPV["NPV_unw"] = new MultiSamplePlot(datasets, "NPV_unw", 51, -0.5, 50.5, "Number of PV","Events", category); 

  
 

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
  pathPNG += "/";
  pathPNG += category;
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


TMVA::Factory *TrainingFACTORY(string trName, TFile *outfile)
{
	TMVA::Factory *factory = new TMVA::Factory((trName).c_str(),outfile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );

  if(trName.find("hut")!=string::npos && trName.find("j3")!=string::npos)
  {
	    factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	    factory->AddVariable("MVA_TOPHLEPBB",'D');
	    factory->AddVariable("LepCharge",'I');
	    factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
	    nTrainingVars = 12;
  }
  else if(trName.find("hct")!=string::npos && trName.find("j3")!=string::npos)
  {
	    factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	    factory->AddVariable("MVA_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
	    nTrainingVars = 11;
  }  
  else if(trName.find("hut")!=string::npos && trName.find("j4")!=string::npos)
  {
	    factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	    factory->AddVariable("TopHadMass_TOPTOPLEPHAD",'D');
	    factory->AddVariable("MVA_TOPHLEPBB",'D');
	    factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	    factory->AddVariable("LepCharge",'I');
	    factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
//	    factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	    factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');

//	    factory->AddVariable("HiggsMass_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("MVA_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepMass_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",'D');

/*	    factory->AddVariable("TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB",'D');
	    factory->AddVariable("TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB",'D');
	    factory->AddVariable("TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD",'D');
*/
	    nTrainingVars = 13;
  }
  else if(trName.find("hct")!=string::npos && trName.find("j4")!=string::npos)
  {
	    factory->AddVariable("HiggsMass_TOPHLEPBB",'D');
	    factory->AddVariable("TopHadMass_TOPTOPLEPHAD",'D');
	    factory->AddVariable("MVA_TOPHLEPBB",'D');
	    factory->AddVariable("MVA_TOPTOPLEPHAD",'D');
	    factory->AddVariable("HiggsEta_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepMass_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepMass_TOPTOPLEPHAD",'D');
//	    factory->AddVariable("TopLepPt_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepEta_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepHiggsDr_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepTopHadDr_TOPTOPLEPHAD",'D');
	    factory->AddVariable("HiggsBJet1CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("HiggsBJet2CSVv2_TOPHLEPBB",'D');
	    factory->AddVariable("TopLepBJetCSVv2_TOPHLEPBB",'D');
//	    factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadBJetCSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadWNonBJet1CSVv2_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadWNonBJet2CSVv2_TOPTOPLEPHAD",'D');

//	    factory->AddVariable("HiggsMass_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("MVA_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepMass_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet1HiggsBJet2Dr_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepHiggsDr_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet1CSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("HiggsBJet2CSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopLepBJetCSVv2_TOPTOPLEPHBB",'D');
//	    factory->AddVariable("TopHadNonBJetCSVv2_TOPTOPLEPHBB",'D');

/*	    factory->AddVariable("TopHadNonBJetTopLepBJet_SumInclCharge_TOPTOPLEPHBB",'D');
	    factory->AddVariable("TopHadNonBJetLep_SumInclCharge_TOPTOPLEPHBB",'D');
	    factory->AddVariable("TopHadBJetTopLepBJet_SumInclCharge_TOPTOPLEPHAD",'D');
	    factory->AddVariable("TopHadBJetLep_SumInclCharge_TOPTOPLEPHAD",'D');
*/
	    nTrainingVars = 12;
  }
  else
  {
      cerr << "No correct signal selected" << endl;
      return 0;
  }

  return factory;

}



