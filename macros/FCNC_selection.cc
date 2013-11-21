// isis.marina.van.parijs@cern.ch 
// kevin.deroover@cern.ch
// 2013
// This is a program that runs over the toptrees and calculates the 
// efficiencies of certain cuts in the datasamples. 

#include "TStyle.h"
#include "TH3F.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

#include "../../TopTreeProducer/interface/TRootRun.h"
#include "../../TopTreeProducer/interface/TRootEvent.h"

#include "../../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"

#include "../../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "../../TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

#include "../../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../../TopTreeAnalysisBase/Content/interface/Dataset.h"

#include "../../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"



#include "../macros/Style.C"

using namespace std;	//needed for cout and stuff
using namespace TopTree;	//needed for TT
using namespace reweight;  //needed for PUreweighting


/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;





int main(int argc, char *argv[]){
	//Make plots nicer: color, style, ... 
	setMyStyle();

	//see how long the program takes to run with a clock 
	clock_t start = clock(); 
	std::cout << "******************************************" << std::endl; 
	std::cout << " Starting clock" << endl; 
	std::cout << "******************************************"<<std::endl; 
	std::cout << " Beginning of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl; 

	// bool for debugging
	bool debug = false; 
	bool warnings = true; 
	bool information = true; 

        //set the xml file
	string xmlfile = "config/FCNC_config.xml";     //place of the xml file 
	
	//set the channel 
	string channel = "undefined";
	
	//set a default luminosity in pb^-1
	float Luminosity = 20000;


	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// Controll flags for scale factor shifts etc.. ///////////
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
  	cout << "dobTagEffShift: " << dobTagEffShift << endl;

  	int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  	cout << "domisTagEffShift: " << domisTagEffShift << endl;
	

	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   different options for executing this macro          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	std::string tempxml; 
	bool foundxml = false; 
    
    	for(int iarg = 0; iarg < argc && argc>1 ; iarg++)
	{
        	std::string argval=argv[iarg];
		
        	if(argval=="--help" || argval =="--h")
		{
			cout << "--xml myxml.xml: change Xml file" << endl;
			cout << "--1gamma: use the 1 gamma channel" << endl;
			cout << "--2gamma: use the 2 gamma channel" << endl; 
			cout << "--1L3B: use the 1 lepton + 3 b-tags channel" << endl; 
			cout << "--SSdilepton: use the same sign dilepton channel" << endl;
			cout << "--OSdilepton: use the opposite sign dilepton channel" << endl;
			cout << "--3L: use the 3 lepton channel (exactly 3)" << endl;
			cout << "--4L: use the 4 lepton channel (at least 4)" << endl;
                	return 0;
        	}
		if (argval=="--xml") {
			iarg++;
			tempxml = argv[iarg];
			foundxml = true; 
		}
		if (argval=="--1gamma") {
                	channel = "1gamma";
			xmlfile = "config/FCNC_1gamma_config.xml";
        	}
		if (argval=="--2gamma") {
                	channel = "2gamma";
			xmlfile = "config/FCNC_2gamma_config.xml";
        	}
		if (argval=="--1L3B") {
                	channel = "1L3B";
			xmlfile = "config/FCNC_1L3B_config.xml";
        	}
		if (argval=="--SSdilepton") {
                	channel = "SSdilepton";
			xmlfile = "config/FCNC_SSdilepton_config.xml";
        	}
		if (argval=="--OSdilepton") {
                	channel = "OSdilepton";
			xmlfile = "config/FCNC_OSdilepton_config.xml";
        	}
		if (argval=="--3L") {
                	channel = "3L";
			xmlfile = "config/FCNC_3L_config.xml";
        	}
		if (argval=="--4L") {
                	channel = "4L";
			xmlfile = "config/FCNC_4L_config.xml";
        	}


    	} 
	
    	if (foundxml)
	{
		xmlfile = tempxml; 
	}
	if(information)	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl;
	if(information)	std::cout << "[INFO]	Used channel: " << channel << endl;
	if(channel.find("undefined")!=string::npos && warnings) std:cout << "[WARNING]	No channel was defined" << endl; 
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   end different options for executing this macro      //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// Options for different b-tagging algorithms	      /////
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//Choose which b-tag algorithm will be used.
	string btagger = "CSVM";
	
	// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10, TCHEM: data/MC scalefactor = 0.94 +- 0.09
	// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12, TCHEM: data/MC scalefactor = 1.21 +- 0.17
 	 float scalefactorbtageff, mistagfactor;
  	if(btagger == "TCHPM" || btagger == "TCHET" || btagger == "SSV" ){
    	cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    	exit(1);
	}
  	else if(btagger == "TCHEM") //redundant for now (these values need updating), but will use as skeleton for CSVM
  	{
           if(dobTagEffShift == 0)
                scalefactorbtageff = 0.94;
         if(dobTagEffShift == 1)
                scalefactorbtageff = 0.85;
         if(dobTagEffShift == 2)
                scalefactorbtageff = 1.03;
                
         if(domisTagEffShift == 0)
                mistagfactor = 1.21;
         if(domisTagEffShift == 1)
                mistagfactor = 1.04;
         if(domisTagEffShift == 2)
                mistagfactor = 1.38;
  	}
    
  	float workingpointvalue = 9999; //working points updated to 2012 BTV-POG recommendations.
 
  	if(btagger == "TCHPM" || btagger == "TCHET" || btagger == "SSV" ){
    	cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    	exit(1);
  	}
   	else if(btagger == "CSVL")
     	workingpointvalue = .244;        
  	else if(btagger == "CSVM")
    	workingpointvalue = .679;
  	else if(btagger == "CSVT")
    	workingpointvalue = .898;
	
	/////////////////////////////////////////////////////////
	// End b-tagging working points      ////////////////////
	/////////////////////////////////////////////////////////

	//Load the analysisenvironment
	AnalysisEnvironment anaEnv; 
	if(debug) std::cout << "[PROCES]	Loading the analysisenvironment" << endl; 
	AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile.c_str());    //load via the xml file the environment
	
	
	//Load the datasets
	TTreeLoader treeLoader; 
	vector <Dataset*> datasets; //vector that will contain all datasets
	if(debug) std::cout << "[PROCES]	Loading the datasets " << endl; 
	treeLoader.LoadDatasets(datasets, xmlfile.c_str()); //put datasets via xmlfile in the dataset vector

	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                       output stuff                    //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	// Set an output rootfile
	char rootFileName[900];
	char channelchar[900];
	if(channel.find("1gamma")!=string::npos)	sprintf(channelchar, "1gamma");
	if(channel.find("2gamma")!=string::npos)	sprintf(channelchar, "2gamma");
	if(channel.find("3L")!=string::npos)		sprintf(channelchar, "3L");	
	if(channel.find("4L")!=string::npos)		sprintf(channelchar, "4L");
	if(channel.find("1L3B")!=string::npos)		sprintf(channelchar, "1L3B");
	if(channel.find("SSdilepton")!=string::npos)	sprintf(channelchar, "SSdilepton");
	if(channel.find("OSdilepton")!=string::npos)	sprintf(channelchar, "OSdilepton");
	
	sprintf(rootFileName,"Output/FCNC_selection_%s.root",channelchar);
	TFile *fout = new TFile (rootFileName, "RECREATE");
	if(debug) cout << "[PROCES]	Declared output rootfiles  "<< endl;
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                        end output stuff               //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	
	// Declare variables: 
	if(debug) cout << "[PROCES]	Variable declaration  "<< endl;
  	vector < TRootVertex* >   vertex;
  	vector < TRootMuon* >     init_muons;
  	vector < TRootElectron* > init_electrons;
  	vector < TRootJet* >      init_jets;
  	vector < TRootMET* >      mets;

	
	//Define an event (global variable)
	TRootEvent* event = 0;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////// MultiSample plots: convenient class which combines multiple MC and DATA histograms into single plots. //////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	MSPlot["NbOfSelectedJets"] = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0., 15., "Nb. of jets");
    	MSPlot["NbOfSelectedLightJets"] = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 15, 0., 15, "Nb. of jets");
    	MSPlot["NbOfSelectedBJets"] = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0., 8., "Nb. of jets");
    	MSPlot["JetEta"] = new MultiSamplePlot(datasets, "JetEta", 30,-3., 3., "Jet #eta");
    	MSPlot["JetPhi"] = new MultiSamplePlot(datasets, "JetPhi", 50, -4., 4., "Jet #phi");
	MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 40, 0., 700., "MET");
	
	

	//////////////////  Cut flow histograms	/////////////////////////////
	char plotTitle_total_B[900];
	sprintf(plotTitle_total_B,"The total cutflow for %s channel (B)",channelchar); 
	histo1D["cutflow_total_B"] = new TH1F("cutflow_total_B", plotTitle_total_B, 6, -0.5,5.5);
	//histo1D["cutflow_total_B"]->Sumw2();
	histo1D["cutflow_total_B"]->GetYaxis()->SetTitle("Eff.");

	char plotTitle_total_S[900];
	sprintf(plotTitle_total_S,"The total cutflow for %s channel (S)",channelchar); 
	histo1D["cutflow_total_S"] = new TH1F("cutflow_total_S", plotTitle_total_S, 6, -0.5,5.5);
	//histo1D["cutflow_total_S"]->Sumw2();
	histo1D["cutflow_total_S"]->GetYaxis()->SetTitle("Eff.");


	// Define different cutflow plots for each channel and dataset	
	for(unsigned int d = 0; d < datasets.size();d++){ //loop over datasets in order to pre-define cutflow histograms for every process
		
		//Load datasets
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		char datasetNamechar[900];
		if(datasetName.find("ttbar")!=string::npos) {sprintf(datasetNamechar,"ttbar");}
		if(datasetName.find("Wjets")!=string::npos) {sprintf(datasetNamechar,"wjets");}
		if(datasetName.find("ttt")!=string::npos) {sprintf(datasetNamechar,"ttt");}
		if(datasetName.find("ttW")!=string::npos) {sprintf(datasetNamechar,"ttw");}
		if(datasetName.find("WZ")!=string::npos) {sprintf(datasetNamechar,"wz");}
		if(datasetName.find("ZZ")!=string::npos) {sprintf(datasetNamechar,"zz");}
		if(datasetName.find("ttZ")!=string::npos) {sprintf(datasetNamechar,"ttz");}
		if(datasetName.find("Zjets")!=string::npos) {sprintf(datasetNamechar,"Zjets");}
		
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToBB_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToBB_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToBB_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToBB_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL_HctL")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToLL_HctL");}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL_HctR")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToLL_HctR");}
		if(datasetName.find("TTJetsTocZbW")!=string::npos) {sprintf(datasetNamechar,"TTJetsTocZbW");}


		// Define different plots for each channel and dataset
		char plotTitle[900];
		char NamePlot[900];
		sprintf(plotTitle,"The cutflow for %s channel: %s dataset",channelchar,datasetNamechar); 
		sprintf(NamePlot,"cutflow_%s",datasetNamechar);
				
		string Process_cutflow = "cutflow_";
		Process_cutflow +=datasetNamechar;
		
		histo1D[Process_cutflow] = new TH1F(NamePlot, plotTitle, 6, -0.5,5.5);
		//histo1D[Process_cutflow]->Sumw2();
		histo1D[Process_cutflow]->GetYaxis()->SetTitle("Eff.");
	}
	
	
	if(debug) cout << "[PROCES]	Declared cutflow histograms  "<< endl;
  	
	//Defining a directory in which .png files of all the plots created will be stored.
	char pathPNG[900];
	sprintf(pathPNG,"FCNC_%s_MSPlots_MCStudy/",channelchar);
  	//string pathPNG = "FCNC_%s";
  	//pathPNG += "_MSPlots_MCStudy/";
  	//mkdir(pathPNG.c_str(),0777);	
	mkdir(pathPNG,0777);
	if(debug) cout << "[PROCES]	Declared PNG directory  "<< endl;
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                START LOOPING OVER THE DATASETS        //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	if(information)	cout << "[PROCES]	Looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	for(unsigned int d = 0; d < datasets.size();d++)
	{
	bool is_signal = false;
		//Load datasets
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		char datasetNamechar[900];
		if(datasetName.find("ttbar")!=string::npos) {sprintf(datasetNamechar,"ttbar");}
		if(datasetName.find("Wjets")!=string::npos) {sprintf(datasetNamechar,"wjets");}
		if(datasetName.find("ttt")!=string::npos) {sprintf(datasetNamechar,"ttt");}
		if(datasetName.find("ttW")!=string::npos) {sprintf(datasetNamechar,"ttw");}
		if(datasetName.find("WZ")!=string::npos) {sprintf(datasetNamechar,"wz");}
		if(datasetName.find("ZZ")!=string::npos) {sprintf(datasetNamechar,"zz");}
		if(datasetName.find("ttZ")!=string::npos) {sprintf(datasetNamechar,"ttz");}
		if(datasetName.find("Zjets")!=string::npos) {sprintf(datasetNamechar,"Zjets");}
				
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToWW_WToLNuL_HctR");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToBB_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToBB_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToBB_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToBB_HctR");
			is_signal = true;
		}		
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR");
			is_signal = true;
		}		
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR");
			is_signal = true;
		}		
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR");
			is_signal = true;
		}		
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL_HctL")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToLL_HctL");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL_HctR")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocHbW_HToZZ_ZToLL_HctR");
			is_signal = true;
		}
		if(datasetName.find("TTJetsTocZbW")!=string::npos) {
			sprintf(datasetNamechar,"TTJetsTocZbW");
			is_signal = true;
		}
		
		string Process_cutflow = "cutflow_";
		Process_cutflow += datasetNamechar;
		
		if(information) cout << "[INFO]	Dataset " << d << " name : " << datasetName << " / title : " << datasets[d]->Title() << endl;
		
		
		///////////////////////////////////////////////////////////
		//                START LOOPING OVER THE EVENTS          //
		///////////////////////////////////////////////////////////

		int NofEvts = 100000;

		int NofRuns = 0; 
		if( NofEvts > datasets[d]->NofEvtsToRunOver()) 
		{
			NofRuns = datasets[d]->NofEvtsToRunOver(); 
		}
		else
		{
			NofRuns = NofEvts; 
		} 
		
		if(information) cout << "[PROCES]	looping over " << NofRuns <<" events "<< endl;
		
		for(int ievent = 0; ievent <NofRuns; ievent++)
		{
			if(ievent%1000 == 0 && information)
			{
				// << flush << "\r" means this line will be overwritten next time 
				std::cout << "[PROCES]	Processing the " << ievent << "th event" << flush << "\r";  
			}    
			
			
			//Load the event 
			event = treeLoader.LoadEvent(ievent, vertex, init_muons, init_electrons, init_jets, mets);

			
			histo1D[Process_cutflow]->Fill(1);
			histo1D["cutflow_total_B"]->Fill(1);
			if(is_signal) histo1D["cutflow_total_S"]->Fill(1);
			
			histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(2, "initial");
			if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(2, "initial");
			histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(2, "initial");
			
			
			//Make a selection 
			Selection selection(init_jets, init_muons,init_electrons,mets);
			//define selection cuts --> have to be validated!!!
			// From the class Selection the following functions are used: 
				// 	void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon)
				//	void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits) 
				//	void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid) 
				//	void Selection::setDiMuonCuts(float Pt, float Eta, float RelIso, float d0) 
				// 	void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
			selection.setJetCuts(20.,5.,0.01,1.,0.98,0.3,0.1);
			selection.setDiMuonCuts(20.,2.4,0.20,999.);
			selection.setDiElectronCuts(20.,2.5,0.15,0.04,0.5,1,0.3,1); 
			selection.setLooseMuonCuts(10.,2.5,0.2);
			selection.setLooseDiElectronCuts(15.0,2.5,0.2,0.5); 
			
			//select the right objects and put them in a vector
			vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
			vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
			vector<TRootMuon*> looseMuons = selection.GetSelectedLooseMuons();
			vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
			vector<TRootElectron*> looseElectrons = selection.GetSelectedLooseDiElectrons();
			vector<TRootJet*> selectedBJets; // B-Jets, to be filled after b-tagging
    			vector<TRootJet*> selectedLightJets; // light-Jets, to be filled afer b-tagging
			// vector<TRootPhoton*> selectedPhotons = selection.GetSelecetedPhotons(); Photons not yet included in the selection class!!!!
			
			
			//order the jets according to the Pt 
			sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order jets wrt Pt.  
			
			
			
			int bjets = 0;
			int nTags = 0;
			float Passed_selection = false;
			
			
			// scale factor for the event
        		float scaleFactor = 1.;

        /*	ONLY HAS TO BE APPLIED WHEN INCLUDING DATA AND/OR PILE-UP REWEIGHTING
	
	
	
        		//A simple filter to remove events which arise from the interaction LHC beam with
        		// the beampipe known as "beam scraping events".
        		if(datasetName == "Data" || datasetName == "data" || datasetName == "DATA")
        		{
                	// Apply the scraping veto. (Is it still needed?)
                	bool isBeamBG = true;
                	if(event->nTracks() > 10)
                	{
                        	if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
                        	isBeamBG = false;
                	}
                      		if(isBeamBG) continue;
        		}
        		else{
        		double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
        		double lumiWeightOLD=lumiWeight;
        		if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
         		lumiWeight=1;
        		scaleFactor = scaleFactor*lumiWeight;
        
        		}
	*/
			
					
      			for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++){
				int pdgID = selectedJets[iJet]->partonFlavour();
				
				if(fabs(pdgID) == 5 && datasetName != "Data" && datasetName != "data" && datasetName != "DATA")
				{
	      				bjets++;
								
	 			}
				
				//filling vector of b-jets
				if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue){
					nTags++;
					selectedBJets.push_back(selectedJets[iJet]);
				}
				else selectedLightJets.push_back(selectedJets[iJet]);
			} 

			
			if(debug) cout << "looseElectrons.size() = " << looseElectrons.size() << endl; 
			if(debug) cout << "looseMuons.size() = " << looseMuons.size() << endl; 
			
			//exactly 3 leptons
			if(channel.find("3L")!=string::npos)
			{
				if(debug) cout << "in 3L channel" << endl;
				if(looseElectrons.size() + looseMuons.size() ==3)
				{ 
					if(debug) cout << "fill 3L" << endl;
					histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "3L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "3L");		
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "3L");
					if(debug) cout << "filled 3L" << endl;
					
					Passed_selection = true;
				}
				if(debug)	cout << "out fill 3L loop" << endl; 
			}
			//more than 4 leptons
			if(channel.find("4L")!=string::npos)
			{
				if(debug) cout << "in 3L channel" << endl;
				if(looseElectrons.size() + looseMuons.size() > 3)
				{ 
					if(debug) cout << "fill 4L" << endl;
					histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "4L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "4L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "4L");
					if(debug) cout << "filled 4L" << endl;
					
					Passed_selection = true;
				}
				if(debug)	cout << "out fill 4L loop" << endl; 
			}
			//1 lepton + 3 b-jets
			if(channel.find("1L3B")!=string::npos)
			{
				if(debug) cout << "in 1L3B channel" << endl;
				if(looseElectrons.size() +  looseMuons.size() == 1)
				{
					if(debug) cout << "in fill 1l3b loop" << endl;
					histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "1L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "1L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "1L");
					if(debug) cout << "selectedJets.size() = " << selectedJets.size() << endl;
					
					if(selectedJets.size() >= 3)
					{
						if(debug) cout << "in fill 1l3b loop: 3jets" << endl;
						histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, ">= 3jets");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, ">= 3jets");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, ">=3jets");
						if(nTags == 3)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							histo1D["cutflow_total_B"]->Fill(4);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(4);
							histo1D[Process_cutflow]->Fill(4);
							histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							
							Passed_selection = true;
						}
					}
				
					if(debug) cout << "out fill 1l3b loop" << endl;
				}
			}
			if(channel.find("SSdilepton")!=string::npos)
			{
				if(debug) cout << "in SSdilepton channel" << endl;
				if(looseElectrons.size() + looseMuons.size() == 2)
				{
					if(debug) cout << "in fill SS dilepton " << endl; 
					histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "2L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "2L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "2L");
					
					bool electron = false; 
					bool muon = false; 
					bool EMu = false;

					if(looseElectrons.size() == 2)
					{
						if(looseElectrons[0]->charge() == looseElectrons[1]->charge()) electron = true; 
						if(debug) cout << "Electron boolean defined" << endl;
					}
				
					if(looseMuons.size() == 2)
					{
						if(looseMuons[0]->charge() == looseMuons[1]->charge()) muon = true; 
						if(debug) cout << "Muon boolean defined" << endl;
					}
					if(looseMuons.size() == 1 && looseElectrons.size() == 1)
					{
						if(looseMuons[0]->charge() == looseElectrons[0]->charge()) EMu = true; 
						if(debug) cout << "EMu boolean defined" << endl;
					}
				
					if(muon || electron || EMu)
					{
						if(debug) cout << "in fill SS dilepton: same sign " << endl;
						histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, "2 SS L");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, "2 SS L");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, "2 SS L");
						
						Passed_selection = true;
					}
					if(debug) cout << "out fill SS dilepton " << endl;
				}
			}
			if(channel.find("OSdilepton")!=string::npos)
			{
				if(debug) cout << "in OSdilepton channel" << endl;
				if(looseElectrons.size() + looseMuons.size() == 2)
				{
					if(debug) cout << "in fill OS dilepton " << endl; 
					histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "2L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "2L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "2L");
					
					bool electron = false; 
					bool muon = false; 
					bool EMu = false;

					if(looseElectrons.size() == 2)
					{
						if(looseElectrons[0]->charge() != looseElectrons[1]->charge()) electron = true; 
						if(debug) cout << "Electron boolean defined" << endl;
					}
				
					if(looseMuons.size() == 2)
					{
						if(looseMuons[0]->charge() != looseMuons[1]->charge()) muon = true; 
						if(debug) cout << "Muon boolean defined" << endl; 
					}
					if(looseMuons.size() == 1 && looseElectrons.size() == 1)
					{
						if(looseMuons[0]->charge() != looseElectrons[0]->charge()) EMu = true; 
						if(debug) cout << "EMu boolean defined" << endl; 
					}
				
					if(muon || electron || EMu)
					{
						if(debug) cout << "in fill OS dilepton: same sign " << endl;
						histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, "2 OS L");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, "2 OS L");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, "2 OS L");
						
						Passed_selection = true;
					}
					if(debug) cout << "out fill OS dilepton " << endl;
				}
			}

			if(channel.find("1gamma")!=string::npos)
			{
				cout << "Photons are not yet included in selection class !!! " << endl; 
			/*	if(debug) cout << "in fill 1gamma " << endl;
				
				if(selectedPhotons.size() == 1)
				{
					if(debug) cout << "in fill 1 gamma: selected photons loop " << endl;
					histo1D[Process_cutflow]->Fill(100/NofRuns);
					histo1D["cutflow_total_B"]->Fill(100/NofRuns)
					if(is_signal) histo1D["cutflow_total_S"]->Fill(100/NofRuns)
					
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "1 photons");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "1 photons");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "1 photons");
				}	
				if(debug) cout << "out fill 1gamma " << endl;
			*/
			}
			if(channel.find("2gamma")!=string::npos)
			{
				cout << "Photons are not yet included in selection class !!! " << endl; 
			/*	if(debug) cout << "in fill 2gamma " << endl;
				
				if(selectedPhotons.size() == 2)
				{
					if(debug) cout << "in fill 2 gamma: selected photons loop " << endl;
					histo1D[Process_cutflow]->Fill(2);
					histo1D["cutflow_total_B"]->Fill(2)
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2)
					
					histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "2 photons");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "2 photons");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "2 photons");
				}	
				if(debug) cout << "out fill 2gamma " << endl;
			*/
			}                                                                
    			
	

			//////////////////////////////////////////////////////////////////////////////////
			// Filling histograms 							//////////
			//////////////////////////////////////////////////////////////////////////////////
			if(Passed_selection){
			MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
		        MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
			MSPlot["MET"]->Fill(mets[0]->E(), datasets[d], true, Luminosity*scaleFactor);
		

			for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){

				MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                 		MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);
			}
			}
	
		}
		
		///////////////////////////////////////////////////////////
		//                END LOOPING OVER THE EVENTS            //
		///////////////////////////////////////////////////////////
		
		
	
	}
	if(information)	cout << "[PROCES]	End of looping over the datasets:  " << datasets.size()<< " datasets" << endl;
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//                END LOOPING OVER THE DATASETS          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////

	
	fout->cd();
	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    	{
        
      
        	MultiSamplePlot *temp = it->second;
        	TH1F *tempHisto_data;
        	TH1F *tempHisto_TTTT;
        	//        temp->addText("CMS preliminary");
        	string name = it->first;
		temp->Draw( name, 0, false, false, false, 1);
      
      	cout <<" looping plots..., name ... "<< name<<endl;
        
        	temp->Write(fout, name, false/*, pathPNG, "pdf"*/);
        	cout <<" written...."<<endl;

  	}
	
	TDirectory* th1dir = fout->mkdir("Histos1D_cutflows");
  	th1dir->cd();
  	for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  	{

    
        	TH1F *temp = it->second;
        	temp->Write();
        	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
        	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  	}
	
	std::cout << "******************************************"<<std::endl; 
	std::cout << " End of the program for the FCNC selection " << std::endl; 
	std::cout << "******************************************"<<std::endl;

	// display how long it took for the program to run
	std::cout << "******************************************" << std::endl; 
	std::cout << " It took the program " << ((double)clock() - start) /CLOCKS_PER_SEC << " to run. " << endl; 
	std::cout << "******************************************" << std::endl; 
	return 0; 
	
}
