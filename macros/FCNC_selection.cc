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
	bool debug = false;  // can be set to true using the options for executing
	bool warnings = true; // can be set to false using the options for executing
	bool information = true; // can be set to false using the options for executing

        //set the xml file: default is bigxml
	string xmlfile = "../config/FCNC_config.xml";     //place of the xml file 
	
	//set the channel 
	string channel = "undefined";
	
	//set a default luminosity in pb^-1
	float Luminosity = 20000;

	//set btag algorithm:CSVM
	string btagger = "CSVM";
	
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	//   different options for executing this macro          //
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	std::string tempxml;
	string tempbtagger;
	bool foundxml = false;
	bool signal_xml = false; 
	bool foundbtag = false;
	bool Big_xml = false;  //Can be set to true with the options
    
    	for(int iarg = 0; iarg < argc && argc>1 ; iarg++)
	{
        	std::string argval=argv[iarg];
		
        	if(argval=="--help" || argval =="--h")
		{
			cout << "--NoWarnings: put warnings off " << endl; 
			cout << "--NoInfo: put information off " << endl; 
			cout << "--debug: put debug output on" << endl; 
			cout << "--xml myxml.xml: change Xml file" << endl;
			cout << "--btag CSVM: change btag algorithm" << endl; 
			cout << "--1L3B: use the 1 lepton + 3 b-tags channel" << endl; 
			cout << "--SSdilepton: use the same sign dilepton channel" << endl;
			cout << "--OSdilepton: use the opposite sign dilepton channel" << endl;
			cout << "--3L: use the 3 lepton channel (exactly 3)" << endl;
			cout << "--4L: use the 4 lepton channel (exactly 4)" << endl;
			cout << "--5L: use the 5 lepton channel (at least 5)" << endl;
			cout << "--45: at least 4 leptons " << endl; 
			cout << "--3L4L: use the 3 and 4 lepton channel (at least 3)" << endl;
			cout << "--Bigxml: use the xml file containing all samples (not channel dependent)" << endl; 
                	return 0;
        	}
		if (argval=="--NoInfo") {
			iarg++;
			information = false;  
		}
		if (argval=="--NoWarnings") {
			iarg++;
			warnings = false;  
		}
		if (argval=="--debug") {
			iarg++;
			debug = true;  
		}		
		if (argval=="--btag") {
			iarg++;
			tempbtagger = argv[iarg];
			foundbtag = true; 
		}
		if (argval=="--xml") {
			iarg++;
			tempxml = argv[iarg];
			foundxml = true; 
		}
		if (argval=="--1gamma") {
                	channel = "1gamma";
			xmlfile = "../config/FCNC_1gamma_config.xml";
        	}
		if (argval=="--2gamma") {
                	channel = "2gamma";
			xmlfile = "../config/FCNC_2gamma_config.xml";
        	}
		if (argval=="--3L4L") {
                	channel = "5";
			xmlfile = "../config/FCNC_3L4L_config.xml";
        	}
		if (argval=="--1L3B") {
                	channel = "1L3B";
			xmlfile = "../config/FCNC_1L3B_config.xml";
        	}
		if (argval=="--SSdilepton") {
                	channel = "SSdilepton";
			xmlfile = "../config/FCNC_SSdilepton_config.xml";
        	}
		if (argval=="--OSdilepton") {
                	channel = "OSdilepton";
			xmlfile = "../config/FCNC_OSdilepton_config.xml";
        	}
		if (argval=="--3L") {
                	channel = "3L";
			xmlfile = "../config/FCNC_3L_config.xml";
        	}
		if (argval=="--4L") {
                	channel = "4L";
			xmlfile = "../config/FCNC_4L_config.xml";
        	}
		if (argval=="--45") {
                	channel = "45";
			xmlfile = "../config/FCNC_45_config.xml";
        	}
		if (argval=="--5L") {
                	channel = "5L";
			xmlfile = "../config/FCNC_5L_config.xml";
        	}
		if (argval=="--Bigxml"){
			Big_xml= true; 
			xmlfile = "../config/FCNC_config.xml";
		
		}
		if (argval=="--signal"){
			signal_xml= true; 
			
		
		}
		


    	} 
	
	double mH = 125; //GeV - higgs mass
	double mW = 80; //GeV - W mass
	double mZ = 90; //GeV - Z mass

    	if (Big_xml)	xmlfile = "../config/FCNC_config.xml"; 
	if (foundbtag) btagger = tempbtagger;
	if (foundxml)  xmlfile = tempxml;
	if(information)	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl;
	if(information)	std::cout << "[INFO]	Used channel: " << channel << endl;
	if(information) std::cout << "[INFO]	Used btag algorithm: " << btagger << endl; 
	if(information) std::cout << "[INFO]	Assumed masses: mH = " << mH << " GeV - mW = " << mW << " GeV - mZ = " << mZ <<" GeV"<< endl;
	if(channel.find("undefined")!=string::npos && warnings) std:cout << "[WARNING]	No channel was defined" << endl; 
	if(Big_xml && warnings) std::cout << "[WARNING]   Using the big xml file" << endl;
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
	float workingpointvalue = 9999; //working points updated to 2012 BTV-POG recommendations.
	float Tightworkingpoint = 0.898;
 
  	if(btagger == "TCHPM" || btagger == "TCHET" || btagger == "SSV" ){
    		cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    		exit(1);
  	}
   	else if(btagger == "CSVL")	workingpointvalue = .244;        
  	else if(btagger == "CSVM")   	workingpointvalue = .679;
  	else if(btagger == "CSVT")    	workingpointvalue = .898;
	
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
	if(channel.find("5L")!=string::npos)		sprintf(channelchar, "5L");
	if(channel.find("45")!=string::npos)		sprintf(channelchar, "4L5L");
	if(channel.find("1L3B")!=string::npos)		sprintf(channelchar, "1L3B");
	if(channel.find("SSdilepton")!=string::npos)	sprintf(channelchar, "SSdilepton");
	if(channel.find("OSdilepton")!=string::npos)	sprintf(channelchar, "OSdilepton");
	
	if(signal_xml)
	{
	   sprintf(rootFileName,"../data/FCNC_selection_%s_signal.root",channelchar);
	}
	else{
	   sprintf(rootFileName,"../data/FCNC_selection_%s.root",channelchar);
	}
	TFile *fout = new TFile (rootFileName, "RECREATE");
	if(information) std::cout << "[INFO]	Made output rootfile: " << rootFileName << endl; 
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
	MSPlot["NbOfSelectedJets"] = new MultiSamplePlot(datasets, "NbOfSelectedJets", 6,-0.5, 5.5, "Nb. of jets");
	MSPlot["NbOfSelectedLightJets"] = new MultiSamplePlot(datasets,"NbOfSelectedLightJets", 8, -0.5, 7.5, "Nb. of light jets");
	MSPlot["NbOfSelectedBJets_CSVM"] = new MultiSamplePlot(datasets,"NbOfSelectedBJets_CSVM", 5, -0.5, 4.5, "Nb. of medium Bjets");
    	MSPlot["NbOfSelectedBJets_CSVT"] = new MultiSamplePlot(datasets,"NbOfSelectedBJets_CSVT", 4, -0.5, 3.5, "Nb. of tight Bjets");
	MSPlot["NbOfSelectedLeptons"] = new MultiSamplePlot(datasets,"NbOfSelectedLeptons", 7, 0.5, 7.5, "Nb. of leptons");
    	MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 40, 0., 160., "MET");
	MSPlot["JetEta"] = new MultiSamplePlot(datasets, "JetEta", 15,-3., 3., "Jet #eta");
    	MSPlot["JetPhi"] = new MultiSamplePlot(datasets, "JetPhi", 25, -4., 4., "Jet #phi");
	MSPlot["JetHt"] = new MultiSamplePlot(datasets, "JetHt", 125, 0., 500, "Sum of Jet Pt");
	MSPlot["LepHt"] = new MultiSamplePlot(datasets, "LepHt", 125, 0., 500, "Sum of Lep Pt");
	MSPlot["JetPx"] = new MultiSamplePlot(datasets, "JetPx", 100, -200, 200, "Sum of Jet Px");
	MSPlot["LepPx"] = new MultiSamplePlot(datasets, "LepPx", 100, -200, 200, "Sum of Lep Px");
	MSPlot["JetPy"] = new MultiSamplePlot(datasets, "JetPy", 100, -200, 200, "Sum of Jet Py");
	MSPlot["LepPy"] = new MultiSamplePlot(datasets, "LepPy", 100, -200, 200, "Sum of Lep Py");
	MSPlot["DeltaPhi_LJ"]= new MultiSamplePlot(datasets, "DeltaPhi_LJ", 50,0, 6, "DeltaPhi Jets-Leptons");
	
	MSPlot["Pt_2nd_leading_Bjet_CSVM"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_Bjet_CSVM",10,30,120,"Pt 2nd leading Bjet M"); 
	MSPlot["Pt_2nd_leading_Bjet_CSVT"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_Bjet_CSVT",10,30,120,"Pt 2nd leading Bjet T"); 
	MSPlot["Pt_2nd_leading_jet"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_jet",10,20,100,"Pt 2nd leading jet");
	MSPlot["Pt_2nd_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_lepton",50,20,150,"Pt 2nd leading lepton");
	MSPlot["Pt_leading_Bjet_CSVM"] = new MultiSamplePlot(datasets,"Pt_leading_Bjet_CSVM",30,20,140,"Pt leading Bjet M"); 
	MSPlot["Pt_leading_Bjet_CSVT"] = new MultiSamplePlot(datasets,"Pt_leading_Bjet_CSVT",30,20,140,"Pt leading Bjet T"); 
	MSPlot["Pt_leading_jet"] = new MultiSamplePlot(datasets,"Pt_leading_jet",20,20,160,"Pt leading jet"); 
	MSPlot["Pt_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_leading_lepton",25,30,140,"Pt leading lepton");

	MSPlot["Pt_3d_leading_jet"] = new MultiSamplePlot(datasets,"Pt_3d_leading_jet",100,0,200,"Pt third leading jet"); 
	MSPlot["Pt_4th_leading_jet"] = new MultiSamplePlot(datasets,"Pt_4th_leading_jet",100,0,200,"Pt fourth leading jet"); 
	
	MSPlot["Pt_3d_leading_Bjet_CSVM"] = new MultiSamplePlot(datasets,"Pt_3d_leading_Bjet_CSVM",100,0,200,"Pt third leading Bjet_CSVM"); 
	MSPlot["Pt_4th_leading_Bjet_CSVM"] = new MultiSamplePlot(datasets,"Pt_4th_leading_Bjet_CSVM",100,0,200,"Pt fourth leading Bjet_CSVM"); 
	MSPlot["Pt_3d_leading_Bjet_CSVT"] = new MultiSamplePlot(datasets,"Pt_3d_leading_Bjet_CSVT",100,0,200,"Pt third leading Bjet_CSVT"); 
	MSPlot["Pt_4th_leading_Bjet_CSVT"] = new MultiSamplePlot(datasets,"Pt_4th_leading_Bjet_CSVT",100,0,200,"Pt fourth leading Bjet_CSVT"); 

	
	if(channel.find("4L")!=string::npos || channel.find("45")!=string::npos)
	{
		MSPlot["MllW"] = new MultiSamplePlot(datasets,"M_{l#nu}",100,80,180,"M_{l#nu} ~ m_{W}");
		 MSPlot["Mllll"] = new MultiSamplePlot(datasets,"M_{llll}",100,80,180,"Mllll ~ m_{H}");
		 MSPlot["NofElectrons_bf"]=new MultiSamplePlot(datasets,"NofElectrons_bf",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_bf"]=new MultiSamplePlot(datasets,"NofMuons_bf",7, -0.5, 6.5,"Nb of muons");
		MSPlot["NofElectrons_af"]=new MultiSamplePlot(datasets,"NofElectrons_af",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_af"]=new MultiSamplePlot(datasets,"NofMuons_af",7, -0.5, 6.5,"Nb of muons");
		MSPlot["NofElectrons_bf_4L"]=new MultiSamplePlot(datasets,"NofElectrons_bf_4L",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_bf_4L"]=new MultiSamplePlot(datasets,"NofMuons_bf_4L",7, -0.5, 6.5,"Nb of muons");
		MSPlot["NofElectrons_af_4L"]=new MultiSamplePlot(datasets,"NofElectrons_af_4L",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_af_4L"]=new MultiSamplePlot(datasets,"NofMuons_af_4L",7, -0.5, 6.5,"Nb of muons");
		MSPlot["NofElectrons_bf_5L"]=new MultiSamplePlot(datasets,"NofElectrons_bf_5L",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_bf_5L"]=new MultiSamplePlot(datasets,"NofMuons_bf_5L",7, -0.5, 6.5,"Nb of muons");
		MSPlot["NofElectrons_af_5L"]=new MultiSamplePlot(datasets,"NofElectrons_af_5L",7, -0.5, 6.5,"Nb of electrons");
		MSPlot["NofMuons_af_5L"]=new MultiSamplePlot(datasets,"NofMuons_af_5L",7, -0.5, 6.5,"Nb of muons");
	}
	
	
	//////////////////  Cut flow histograms	/////////////////////////////
	MSPlot["MScutflow"] = new MultiSamplePlot(datasets,"MScutflow",20,-0.5,19.5, "cutflow"); 
	
	
	if(debug) cout << "[PROCES]	Declared MS histograms  "<< endl;

	char plotTitle_total_B[900];
	sprintf(plotTitle_total_B,"The total cutflow for %s channel (B)",channelchar); 
	histo1D["cutflow_total_B"] = new TH1F("cutflow_total_B", plotTitle_total_B, 11, -0.5,10.5);
	histo1D["cutflow_total_B"]->Sumw2();
	histo1D["cutflow_total_B"]->GetYaxis()->SetTitle("#evts.");


	char plotTitle_total_S[900];
	sprintf(plotTitle_total_S,"The total cutflow for %s channel (S)",channelchar); 
	histo1D["cutflow_total_S"] = new TH1F("cutflow_total_S", plotTitle_total_S, 11, -0.5,10.5);
	histo1D["cutflow_total_S"]->Sumw2();
	histo1D["cutflow_total_S"]->GetYaxis()->SetTitle("Eff.");

	

	
	if(debug) cout << "[PROCES]	Declared cutflow histograms  "<< endl;
  	
	//Defining a directory in which .png files of all the plots created will be stored.
	char pathPNG[900];
	sprintf(pathPNG,"../data/FCNC_%s_MSPlots_MCStudy/",channelchar);
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
		if(datasetName.find("GluGluHiggs4lep")!=string::npos) {sprintf(datasetNamechar,"GluGluHiggs4lep");}
		if(datasetName.find("Wjets")!=string::npos || datasetName.find("wjets")!=string::npos) {sprintf(datasetNamechar,"Wjets");}
		if(datasetName.find("W_1Jets")!=string::npos) {sprintf(datasetNamechar,"W_1Jets");}
		if(datasetName.find("W_2Jets")!=string::npos) {sprintf(datasetNamechar,"W_2Jets");}
		if(datasetName.find("W_3Jets")!=string::npos) {sprintf(datasetNamechar,"W_3Jets");}
		if(datasetName.find("W_4Jets")!=string::npos) {sprintf(datasetNamechar,"W_4Jets");}
		
		if(datasetName.find("WW")!=string::npos || datasetName.find("ww")!=string::npos) {sprintf(datasetNamechar,"WW");}
                if(datasetName.find("WZ")!=string::npos || datasetName.find("wz")!=string::npos) {sprintf(datasetNamechar,"WZ");}
                if(datasetName.find("ZZ")!=string::npos || datasetName.find("zz")!=string::npos) {sprintf(datasetNamechar,"ZZ");}
		if(datasetName.find("WW_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"WW_To2L2Nu");}
		if(datasetName.find("WZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"WZ_To2L2Q");}
		if(datasetName.find("WZ_To3LNu")!=string::npos) {sprintf(datasetNamechar,"WZ_To3LNu");}
		if(datasetName.find("ZZ_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Nu");}
		if(datasetName.find("ZZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Q");}
		if(datasetName.find("ZZ_To4L")!=string::npos) {sprintf(datasetNamechar,"ZZ_To4L");}
		
		if(datasetName.find("ST_T_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_s-ch");}
                if(datasetName.find("ST_TBar_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_s-ch");}
		if(datasetName.find("ST_T_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_t-ch");}
		if(datasetName.find("ST_Tbar_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_Tbar_t-ch");}
		if(datasetName.find("ST_T_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_tW-ch");}
                if(datasetName.find("ST_TBar_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_tW-ch");}
		if(datasetName.find("ST_TToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToDilepton_tW-ch");}
		if(datasetName.find("ST_TToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToThadWlep_tW-ch");}
		if(datasetName.find("ST_TBarToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToDilepton_tW-ch");}
		if(datasetName.find("ST_TBarToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TBarToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToThadWlep_tW-ch");}
		
		if(datasetName.find("TT_SemiLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_SemiLeptMGDecays");}
		if(datasetName.find("TT_FullLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_FullLeptMGDecays");}
		if(datasetName.find("TT_HadronicMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_HadronicMGDecays");}
		if(datasetName.find("TTJets")!=string::npos) {sprintf(datasetNamechar,"TTJets");}
		
		if(datasetName.find("TTZ")!=string::npos) {sprintf(datasetNamechar,"TTZ");}
		if(datasetName.find("TTW")!=string::npos) {sprintf(datasetNamechar,"TTW");}
		
		if(datasetName.find("Z_M-10To50")!=string::npos) {sprintf(datasetNamechar,"Z_M-10To50");}
		if(datasetName.find("Z_M-50")!=string::npos) {sprintf(datasetNamechar,"Z_M-50");}
		if(datasetName.find("Z_1Jets")!=string::npos) {sprintf(datasetNamechar,"Z_1Jets");}
		if(datasetName.find("Z_2Jets")!=string::npos) {sprintf(datasetNamechar,"Z_2Jets");}
		if(datasetName.find("Z_3Jets")!=string::npos) {sprintf(datasetNamechar,"Z_3Jets");}
		if(datasetName.find("Z_4Jets")!=string::npos) {sprintf(datasetNamechar,"Z_4Jets");}
			
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
		
		
		
		if(information) cout << "[INFO]	Dataset " << d << " name : " << datasetName << " / title : " << datasets[d]->Title() << endl;

		///////////////////////////////////////////////////////////
		//                START LOOPING OVER THE EVENTS          //
		///////////////////////////////////////////////////////////

		
		if(information) cout << "[PROCES]	looping over " << datasets[d]->NofEvtsToRunOver() <<" events "<< endl;
		
		for(int ievent = 0; ievent <datasets[d]->NofEvtsToRunOver(); ievent++)
		{
			if(ievent%1000 == 0 && information)
			{
				// << flush << "\r" means this line will be overwritten next time 
				std::cout << "[PROCES]	Processing the " << ievent << "th event" << flush << "\r";  
			}    
			// scale factor for the event
        		float scaleFactor = 1.;
			
			//Load the event 
			event = treeLoader.LoadEvent(ievent, vertex, init_muons, init_electrons, init_jets, mets);

			MSPlot["MScutflow"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
			//histo1D[Process_cutflow]->Fill(1);
			if(!is_signal) histo1D["cutflow_total_B"]->Fill(1);
			if(is_signal) histo1D["cutflow_total_S"]->Fill(1);
			
			if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(2, "initial");
			if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(2, "initial");
			//histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(2, "initial");
			
			
			//Make a preliminary selection 
			Selection selection(init_jets, init_muons,init_electrons,mets);
			//define selection cuts --> have to be validated!!!
			// From the class Selection the following functions are used: 
				// 	void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon) 
				//	void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid)  
				// 	void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
			selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1); 
			selection.setDiMuonCuts(10.,2.5,0.2,0.04);
			selection.setDiElectronCuts(15.0,2.4,0.15,0.04,0.5,1,0.3,1); 
			//void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits)
			//select the right objects and put them in a vector
			vector<TRootJet*> selectedJets = selection.GetSelectedJets(true);
			vector<TRootMuon*> looseMuons = selection.GetSelectedDiMuons();
			vector<TRootElectron*> looseElectrons = selection.GetSelectedDiElectrons();
			
			vector<TRootJet*> selectedBJets_CSVM; // B-Jets, to be filled after b-tagging
    			vector<TRootJet*> selectedBJets_CSVT; // B-jets at the Tight working point
			vector<TRootJet*> selectedLightJets; // light-Jets, to be filled afer b-tagging
			// vector<TRootPhoton*> selectedPhotons = selection.GetSelecetedPhotons(); Photons not yet included in the selection class!!!!
			
			
			//order the jets according to the Pt 
			sort(selectedJets.begin(),selectedJets.end(),HighestPt());   
			sort(looseElectrons.begin(),looseElectrons.end(),HighestPt());
			sort(looseMuons.begin(),looseMuons.end(),HighestPt());
			
			// Selection boolean
			bool Passed_selection = false;
			
			//check missing Et 
			double met_px = 0; 
			double met_py = 0; 
			double met_pt = 0; 
			met_px = mets[0]->Px(); 
			met_py = mets[0]->Py();
			met_pt = sqrt(met_px*met_px + met_py*met_py); 
			if(debug) cout << "[INFO]	met_px = " << met_px << endl; 
			if(debug) cout << "[INFO]	met_py = " << met_py << endl; 
			if(debug) cout << "[INFO]	met_pt = " << met_pt << endl; 
			
			
			
			
	
			//implement btagging	
      			for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++){
				if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue)
				{
					selectedBJets_CSVM.push_back(selectedJets[iJet]);
				}
				else selectedLightJets.push_back(selectedJets[iJet]);
				
				if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > Tightworkingpoint)
				{
					selectedBJets_CSVT.push_back(selectedJets[iJet]);
				}
				
				
			} 

			
			if(debug) cout << "[INFO]	looseElectrons.size() = " << looseElectrons.size() << endl; 
			if(debug) cout << "[INFO]	looseMuons.size() = " << looseMuons.size() << endl; 
			
			

    			//exactly 3 leptons
			if(channel.find("3L")!=string::npos)
			{
				
				if(debug) cout << "[PROCES]	in 3L channel" << endl;
				if(looseElectrons.size() + looseMuons.size() ==3)
				{ 
					if(debug) cout << "[PROCES]	fill 3L" << endl;
					
					//fill histograms
					MSPlot["MScutflow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
					if(!is_signal)histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					
					
					//set labels
					if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "3L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "3L");		
					
					
					Passed_selection = true;
					
					
					if(debug) cout << "[PROCES]	filled 3L" << endl;
					
					
				}
				if(debug)	cout << "[PROCES]	out fill 3L loop" << endl; 
			}
			
			bool additional_3L = false; 
			
			if(channel.find("3L")!=string::npos && Passed_selection)
			{
				
				if(met_pt>45)
				{
					MSPlot["MScutflow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
							
				}
				else
				{
					
				    
				      MSPlot["MScutflow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
				
				}
			     
			        if(selectedJets.size()>3)
				{
					// in t(bW)t(qH) --> b vl qZ(ll)Z(qq)  loop
				  	MSPlot["MScutflow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
				
				
				}	
				else 
				{
					 MSPlot["MScutflow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
				
				}
				
				
				additional_3L = true; 
			
			
			
			
			}
			
			
			
			
			
			
			
			
			
			
			
			if(channel.find("45")!=string::npos)
			{
				
				if(debug) cout << "[PROCES]	in 45 channel" << endl;
				
				if(looseElectrons.size() + looseMuons.size() > 3)
				{ 
					if(debug) cout << "[PROCES]	fill 45" << endl;
					
					//fill histograms
					if(!is_signal)	histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, ">3L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, ">3L");
					
					MSPlot["MScutflow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
					//label histograms
					
					
					Passed_selection = true;
									
					if(debug) cout << "[PROCES]	filled 45" << endl;
					
					
				}
				if(debug)	cout << "[PROCES]	out fill 45 loop" << endl; 
			}
			
			bool additional_45 = false;
			
			bool leptons_45 = false; 
			
			double mllll_45 =0; 
			double mllW_45=0;
			
			TLorentzVector leptonpair_45_mll_1;
			TLorentzVector leptonpair_45_mll_2;
			TLorentzVector leptonpair_45_mllll;
			TLorentzVector leptonpair_45_mllW;
			TLorentzVector lepton_45_0;
			TLorentzVector lepton_45_1;
			TLorentzVector lepton_45_2;
			TLorentzVector lepton_45_3;
			TLorentzVector lepton_45_4;
			TLorentzVector Met_W; 
			Met_W.Clear();
			lepton_45_0.Clear(); 
			lepton_45_1.Clear();
			lepton_45_2.Clear();
			lepton_45_3.Clear();
			lepton_45_4.Clear();
			leptonpair_45_mll_1.Clear();
			leptonpair_45_mll_2.Clear();
			leptonpair_45_mllll.Clear();

			if(channel.find("45")!=string::npos && Passed_selection)
			{
			    if(debug) cout << "in leptoncounting" << endl; 
			    MSPlot["NofElectrons_bf"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			    MSPlot["NofMuons_bf"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
			   if(looseMuons.size()+looseElectrons.size() == 4){
				 	MSPlot["NofElectrons_bf_4L"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			   	 	MSPlot["NofMuons_bf_4L"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
			   }
			   if(looseMuons.size()+looseElectrons.size() == 5){
				 	MSPlot["NofElectrons_bf_5L"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			   	 	MSPlot["NofMuons_bf_5L"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
			   }
			   
			   //select ZZ events, muon flavours come in pairs
			   if(looseMuons.size()+ looseElectrons.size() == 5)
			   {
			   	Met_W.SetPxPyPzE(met_px,met_py,mets[0]->Pz(),mets[0]->Energy());
				if(debug) cout << "defined met" << endl;
				
			   	if( looseMuons.size() == 2 & looseElectrons.size() == 3){ 
			   		leptons_45=true; 
				
					lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			     		lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseElectrons[2]->Px(),looseElectrons[2]->Py(),looseElectrons[2]->Pz(),looseElectrons[2]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			       		lepton_45_3.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			   
			       		leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
					leptonpair_45_mllW = lepton_45_4 + Met_W;
				
			   	} 
			   	if( looseMuons.size() == 3 & looseElectrons.size() == 2)
			   	{
			   		leptons_45=true; 
				
					lepton_45_3.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseMuons[2]->Px(),looseMuons[2]->Py(),looseMuons[2]->Pz(),looseMuons[2]->Energy());
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			        	lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			   
			       		leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_4 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
					leptonpair_45_mllW = lepton_45_2 + Met_W;
				
			   	}
			   	if( looseMuons.size() == 5 & looseElectrons.size() == 0)
			   	{
			   		leptons_45=true; 
				
					lepton_45_3.SetPxPyPzE(looseMuons[3]->Px(),looseMuons[3]->Py(),looseMuons[3]->Pz(),looseMuons[3]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseMuons[4]->Px(),looseMuons[4]->Py(),looseMuons[4]->Pz(),looseMuons[4]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseMuons[2]->Px(),looseMuons[2]->Py(),looseMuons[2]->Pz(),looseMuons[2]->Energy());
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			        	lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			   
			                leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
					leptonpair_45_mllW = lepton_45_4 + Met_W;
			       
				
			   	}
			   	if( looseMuons.size() == 0 & looseElectrons.size() == 5)
			   	{
			   		leptons_45=true; 
				
					lepton_45_3.SetPxPyPzE(looseElectrons[3]->Px(),looseElectrons[3]->Py(),looseElectrons[3]->Pz(),looseElectrons[3]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseElectrons[4]->Px(),looseElectrons[4]->Py(),looseElectrons[4]->Pz(),looseElectrons[4]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseElectrons[2]->Px(),looseElectrons[2]->Py(),looseElectrons[2]->Pz(),looseElectrons[2]->Energy());
			     		lepton_45_0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			        	lepton_45_1.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			   
			                leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
					leptonpair_45_mllW = lepton_45_4 + Met_W;
			       
				
			   	}
			   	if( looseMuons.size() == 4 & looseElectrons.size() == 1)
			   	{
			   		leptons_45=true; 
				
					lepton_45_3.SetPxPyPzE(looseMuons[3]->Px(),looseMuons[3]->Py(),looseMuons[3]->Pz(),looseMuons[3]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseMuons[2]->Px(),looseMuons[2]->Py(),looseMuons[2]->Pz(),looseMuons[2]->Energy());
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			        	lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			   
			   		leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
			        	leptonpair_45_mllW = lepton_45_4 + Met_W;
			       
				
			   	}
			   	if( looseMuons.size() == 1 & looseElectrons.size() == 4)
			   	{
			   		leptons_45=true; 
				
					lepton_45_3.SetPxPyPzE(looseElectrons[2]->Px(),looseElectrons[2]->Py(),looseElectrons[2]->Pz(),looseElectrons[2]->Energy());
			     		lepton_45_4.SetPxPyPzE(looseElectrons[3]->Px(),looseElectrons[3]->Py(),looseElectrons[3]->Pz(),looseElectrons[3]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			        	lepton_45_1.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			   
			        	leptonpair_45_mll_1 = lepton_45_2 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_4 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
			       		leptonpair_45_mllW = lepton_45_0 + Met_W;
				
			   	}
			  
			   
			   
			}
			else if(looseMuons.size()+ looseElectrons.size() == 4)
			{
		            if(debug) cout << "in 45 4L loop" << endl; 
			   if( looseMuons.size() == 2 || looseMuons.size() == 0 || looseMuons.size() == 4)
			   {
			        
				leptons_45= true;
				if(looseMuons.size() == 0)
			   	{
			     		lepton_45_0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			     		lepton_45_1.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseElectrons[2]->Px(),looseElectrons[2]->Py(),looseElectrons[2]->Pz(),looseElectrons[2]->Energy());
			     		lepton_45_3.SetPxPyPzE(looseElectrons[3]->Px(),looseElectrons[3]->Py(),looseElectrons[3]->Pz(),looseElectrons[3]->Energy());
			   
			    	 	leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
			   	}
			   	if(looseElectrons.size() == 0)
			   	{
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			     		lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseMuons[2]->Px(),looseMuons[2]->Py(),looseMuons[2]->Pz(),looseMuons[2]->Energy());
			     		lepton_45_3.SetPxPyPzE(looseMuons[3]->Px(),looseMuons[3]->Py(),looseMuons[3]->Pz(),looseMuons[3]->Energy());
			   
			   		leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
			   	}
			   	if(looseElectrons.size() == 2)
			   	{
			     		lepton_45_0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
			     		lepton_45_1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
			     		lepton_45_2.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
			     		lepton_45_3.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
			      	
					leptonpair_45_mll_1 = lepton_45_0 + lepton_45_1; 
			     		leptonpair_45_mll_2 = lepton_45_2 + lepton_45_3;
			     		leptonpair_45_mllll = leptonpair_45_mll_1 + leptonpair_45_mll_2;
				}
					
				  	
			   }
			
			
			 }
			}
			
			if(leptons_45)
			{	
				if(debug) cout << "in cutflow" << endl; 
				 MSPlot["MScutflow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
				 if(!is_signal)	histo1D["cutflow_total_B"]->Fill(3);
				  if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
				  if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, "Zdecay");
				  if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, "Zdecay");
				 
				 MSPlot["NofElectrons_af"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			   	 MSPlot["NofMuons_af"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
				 if(looseMuons.size()+looseElectrons.size() == 4){
				 	MSPlot["NofElectrons_af_4L"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			   	 	MSPlot["NofMuons_af_4L"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
				}
				if(looseMuons.size()+looseElectrons.size() == 5){
				 	MSPlot["NofElectrons_af_5L"]->Fill(looseElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
			   	 	MSPlot["NofMuons_af_5L"]->Fill(looseMuons.size(), datasets[d], true, Luminosity*scaleFactor);
				}
				
				mllll_45 = leptonpair_45_mllll.M();
				mllW_45 = leptonpair_45_mllW.M();
			  	MSPlot["Mllll"]->Fill(mllll_45, datasets[d], true, Luminosity*scaleFactor);
				MSPlot["MllW"]->Fill(mllW_45, datasets[d], true, Luminosity*scaleFactor);
				 	
			}
			if(fabs(mllll_45-mH) < 10 & leptons_45)
			{
				MSPlot["MScutflow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
				if(!is_signal)	histo1D["cutflow_total_B"]->Fill(4);
				if(is_signal) histo1D["cutflow_total_S"]->Fill(4);
				if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(5, "M_{H}^{reco}");
				if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(5, "M_{H}^{reco}");
				 
			}
			
			if(channel.find("45")!=string::npos && Passed_selection)
			{
				if(debug) cout << "In pt px py loop" << endl; 
				TLorentzVector Leptons; 
				TLorentzVector Jets;
				Leptons.Clear();
				Jets.Clear(); 
				double Ht_jet = 0; 
				double Ht_lep = 0; 
				double Px_jet = 0;
				double Px_lep = 0;
				double Py_jet = 0;
				double Py_lep = 0; 
				double DeltaPhi = 0; 
				for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
				     if(debug) cout << "in jets" << endl; 
			             Ht_jet += selectedJets[seljet1]->Pt();
				     TLorentzVector tempJet; 
				     tempJet.Clear();
				     tempJet.SetPxPyPzE(selectedJets[seljet1]->Px(),selectedJets[seljet1]->Py(),selectedJets[seljet1]->Pz(),selectedJets[seljet1]->Energy());
				     Jets += tempJet; 
			        }
				for (Int_t selLep =0; selLep < looseMuons.size(); selLep++ ){
				     if(debug) cout << "in muons" << endl; 
			             Ht_lep += looseMuons[selLep]->Pt();
				     TLorentzVector tempLep; 
				     tempLep.Clear();
				     tempLep.SetPxPyPzE(looseMuons[selLep]->Px(),looseMuons[selLep]->Py(),looseMuons[selLep]->Pz(),looseMuons[selLep]->Energy());
				    
				     Leptons += tempLep; 
			        }
				if(debug) cout << "looseElectrons.size() " << looseElectrons.size() << endl; 
				for (Int_t selLep =0; selLep < looseElectrons.size(); selLep++ ){
				     if(debug) cout << "in electrons" << endl; 
			             Ht_lep += looseElectrons[selLep]->Pt();
				     TLorentzVector tempLep1; 
				     tempLep1.Clear();
				     tempLep1.SetPxPyPzE(looseElectrons[selLep]->Px(),looseElectrons[selLep]->Py(),looseElectrons[selLep]->Pz(),looseElectrons[selLep]->Energy());
				    
				     Leptons += tempLep1; 
			        }
				
				Px_jet = Jets.Px(); 
				Py_jet = Jets.Py(); 
				Px_lep = Leptons.Px(); 
				Py_lep = Leptons.Py(); 
				 
				DeltaPhi = fabs(Jets.Phi() - Leptons.Phi());
				if(debug){
					cout << "ht lep " << Ht_lep << endl; 
					cout << "ht jet " << Ht_jet << endl; 
					cout << "Px lep " << Px_lep << endl; 
					cout << "Px jet " << Px_jet << endl;
					cout << "Py lep " << Py_lep << endl; 
					cout << "Py jet " << Py_jet << endl;
					cout << "delta phi " << DeltaPhi << endl; 
				}
				
				MSPlot["JetHt"]->Fill(Ht_jet , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "jetht filled " << endl; 
				MSPlot["LepHt"]->Fill(Ht_lep , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "lepht filled " << endl; 
				MSPlot["JetPx"]->Fill(Px_jet , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "jetpx filled " << endl; 
				MSPlot["LepPx"]->Fill(Px_lep , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "leppx filled " << endl; 
				MSPlot["JetPy"]->Fill(Py_jet , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "jetpy filled " << endl; 
				MSPlot["LepPy"]->Fill(Py_lep , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "leppy filled " << endl; 
                 		MSPlot["DeltaPhi_LJ"]->Fill(DeltaPhi , datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "dphi filled " << endl; 
				
			}
			
			if(channel.find("45")!=string::npos && Passed_selection)
			{
			    additional_45 = true;
			    if(selectedJets.size()>1)
			    {
			        MSPlot["MScutflow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
				if(selectedBJets_CSVM.size() == 1)
				{
					MSPlot["MScutflow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
					
				}
			    }
			    
			}   
			   
			
			
			
			
			
			

			//////////////////////////////////////////////////////////////////////////////////
			// Filling histograms 							//////////
			//////////////////////////////////////////////////////////////////////////////////





			
			if(additional_3L ||additional_45)
			{
				if(debug) cout << "[PROCES]	In passed_selection loop" << endl; 
				
				MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "[PROCES]	Filled NbOfSelectedJets" << endl; 
				MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
		        	if(debug) cout << "[PROCES]	Filled NbOfSelectedLightJets" << endl; 
				MSPlot["NbOfSelectedBJets_CSVM"]->Fill(selectedBJets_CSVM.size(), datasets[d], true, Luminosity*scaleFactor);
		        	if(debug) cout << "[PROCES]	Filled NbOfSelectedBJets_CSVM" << endl; 
				MSPlot["NbOfSelectedBJets_CSVT"]->Fill(selectedBJets_CSVT.size(), datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "[PROCES]	Filled NbOfSelectedJets_CSVT" << endl; 
				MSPlot["NbOfSelectedLeptons"]->Fill(looseMuons.size()+looseElectrons.size(),datasets[d],true,Luminosity*scaleFactor);
				if(debug) cout << "[PROCES]	Filled NbOfSelectedLeptons" << endl; 
				if(debug) cout << "[PROCES]	Filling MET with " << (float) met_pt << endl;
				MSPlot["MET"]->Fill((float) met_pt, datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "[PROCES]	Filled MET" << endl; 
				
				for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
					
					MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                 			MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);
				}
                                if( selectedJets.size() > 0) {
					if(debug) cout << "[PROCES]	In selectedJets.size() > 0" << endl; 
					MSPlot["Pt_leading_jet"]->Fill(selectedJets[0]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
                                	if(debug) cout << "[PROCES]	Out selectedJets.size() > 0" << endl; 
				}
				if( selectedJets.size() > 1)
				{
					MSPlot["Pt_2nd_leading_jet"]->Fill(selectedJets[1]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedJets.size() > 2) {
					MSPlot["Pt_3d_leading_jet"]->Fill(selectedJets[2]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
                                	
				}
				if( selectedJets.size() > 3)
				{
					MSPlot["Pt_4th_leading_jet"]->Fill(selectedJets[3]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				
				if( selectedBJets_CSVM.size() > 0) {
					MSPlot["Pt_leading_Bjet_CSVM"]->Fill(selectedBJets_CSVM[0]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
                                }
				if( selectedBJets_CSVT.size() > 0)
				{
					MSPlot["Pt_leading_Bjet_CSVT"]->Fill(selectedBJets_CSVM[0]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 1)
				{
					MSPlot["Pt_2nd_leading_Bjet_CSVM"]->Fill(selectedBJets_CSVM[1]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVT.size() > 1)
				{
					MSPlot["Pt_2nd_leading_Bjet_CSVT"]->Fill(selectedBJets_CSVT[1]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 2)
				{
				  	MSPlot["Pt_3d_leading_Bjet_CSVM"]->Fill(selectedBJets_CSVM[2]->Pt(), datasets[d],true,Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVT.size() > 2)
				{	
				  	MSPlot["Pt_3d_leading_Bjet_CSVT"]->Fill(selectedBJets_CSVT[2]->Pt(), datasets[d],true,Luminosity*scaleFactor);
				}
				
				
				TLorentzVector lepton0;
				TLorentzVector lepton1;
				
				lepton0.Clear(); 
				lepton1.Clear(); 
				if(looseElectrons.size()==0)
				{
					if(looseMuons.size()!=0)
					{
						lepton0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						
						
						if(looseMuons.size()>1 )
						{
							lepton1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
							
							MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
							MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
						}
					}
					 
				}
				else if (looseMuons.size()==0)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size()==0" << endl; 
					if(looseElectrons.size()!=0)
					{
						lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						
						if(looseElectrons.size()>1 )
						{
							lepton1.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
							
							MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
							MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
						}
					}
				}
				else if(looseMuons.size()>1 && looseElectrons.size()>1)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size()>1 && looseElectrons.size()>1" << endl; 
					double px12 = (looseElectrons[1]->Px())*(looseElectrons[1]->Px());
					double py12 = (looseElectrons[1]->Py())*(looseElectrons[1]->Py());
					double pt_electron1=TMath::Sqrt(px12+py12);
					if(debug) cout << "[INFO]	pt_electron1 = " << pt_electron1 << endl; 
					double pt_electron0=TMath::Sqrt((looseElectrons[0]->Px())*(looseElectrons[0]->Px())+(looseElectrons[0]->Py())*(looseElectrons[0]->Py()));
					if(debug) cout << "[INFO]	pt_electron0 = " << pt_electron0 << endl;
					double pt_muon1=TMath::Sqrt((looseMuons[1]->Px())*(looseMuons[1]->Px())+(looseMuons[1]->Py())*(looseMuons[1]->Py()));
					if(debug) cout << "[INFO]	pt_muon1 = " << pt_muon1 << endl;
					double pt_muon0=TMath::Sqrt((looseMuons[0]->Px())*(looseMuons[0]->Px())+(looseMuons[0]->Py())*(looseMuons[0]->Py()));
					if(debug) cout << "[INFO]	pt_muon0 = " << pt_muon0 << endl;
					if(pt_electron0 > pt_muon0)
					{
						if(debug) cout << "[PROCES]	in pt_electron0 > pt_muon0" << endl; 
						lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						if(pt_muon0 >= pt_electron1)
						{
							lepton1.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						}
						else if(pt_muon0 < pt_electron1)
						{
							lepton1.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
						}
						
					}
					else if (pt_muon0 > pt_electron0)
					{
						lepton0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						if(pt_electron0 > pt_muon1)
						{
							lepton1.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						}
						else if (pt_electron0 < pt_muon1)
						{
							lepton1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
						
						}
					
					}
					
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				else if(looseMuons.size() == 1 && looseElectrons.size()>1)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size() == 1 && looseElectrons.size()>1" << endl; 
					double px12 = (looseElectrons[1]->Px())*(looseElectrons[1]->Px());
					double py12 = (looseElectrons[1]->Py())*(looseElectrons[1]->Py());
					double pt_electron1=TMath::Sqrt(px12+py12);
					if(debug) cout << "[INFO]	pt_electron1 = " << pt_electron1 << endl; 
					double pt_electron0=TMath::Sqrt((looseElectrons[0]->Px())*(looseElectrons[0]->Px())+(looseElectrons[0]->Py())*(looseElectrons[0]->Py()));
					if(debug) cout << "[INFO]	pt_electron0 = " << pt_electron0 << endl;
					double pt_muon0=TMath::Sqrt((looseMuons[0]->Px())*(looseMuons[0]->Px())+(looseMuons[0]->Py())*(looseMuons[0]->Py()));
					if(debug) cout << "[INFO]	pt_muon0 = " << pt_muon0 << endl;
					if(pt_electron0 > pt_muon0)
					{
						if(debug) cout << "[PROCES]	in pt_electron0 > pt_muon0" << endl; 
						lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						if(pt_muon0 >= pt_electron1)
						{
							lepton1.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						}
						else if(pt_muon0 < pt_electron1)
						{
							lepton1.SetPxPyPzE(looseElectrons[1]->Px(),looseElectrons[1]->Py(),looseElectrons[1]->Pz(),looseElectrons[1]->Energy());
						}
						
					}
					else if (pt_muon0 > pt_electron0)
					{
						lepton0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						lepton1.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
					}
					
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				else if(looseMuons.size() > 1 && looseElectrons.size() == 1)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size() > 1 && looseElectrons.size() == 1" << endl; 
					double pt_electron0=TMath::Sqrt((looseElectrons[0]->Px())*(looseElectrons[0]->Px())+(looseElectrons[0]->Py())*(looseElectrons[0]->Py()));
					if(debug) cout << "[INFO]	pt_electron0 = " << pt_electron0 << endl;
					double pt_muon0=TMath::Sqrt((looseMuons[0]->Px())*(looseMuons[0]->Px())+(looseMuons[0]->Py())*(looseMuons[0]->Py()));
					if(debug) cout << "[INFO]	pt_muon0 = " << pt_muon0 << endl;
					double pt_muon1=TMath::Sqrt((looseMuons[1]->Px())*(looseMuons[1]->Px())+(looseMuons[1]->Py())*(looseMuons[1]->Py()));
					if(debug) cout << "[INFO]	pt_muon1 = " << pt_muon1 << endl;
					if(pt_electron0 > pt_muon0)
					{
						if(debug) cout << "[PROCES]	in pt_electron0 > pt_muon0" << endl; 
						lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						lepton1.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
					}
					else if (pt_muon0 > pt_electron0)
					{
						lepton0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						if(pt_electron0 > pt_muon1)
						{
							lepton1.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
						}
						else if (pt_electron0 < pt_muon1)
						{
							lepton1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
						
						}
					
					}
					
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				else if(looseMuons.size() == 1 && looseElectrons.size() == 1)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size() == 1 && looseElectrons.size() == 1" << endl; 
					lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
					lepton1.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
					
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				
				
				
				
			}
			
			//end passed selection loop
			
			
			
			
			

			
	
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
	if(information)	cout << "[PROCES]	Looping over MSplots" << endl;  
	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    	{
  
  
   
               
        	MultiSamplePlot *temp = it->second;
        	TH1F *tempHisto_data;
        	TH1F *tempHisto_back;
        	temp->addText("CMS simulation");
        	string name = it->first;
		temp->Draw( name, 0, false, false, false, 1);
      		
		
      		if(debug) cout <<" looping MS plots..., name ... "<< name<<endl;
        
        	temp->Write(fout, name, false);//, pathPNG, "pdf");
		
        	if(debug) cout <<" written MSplot " << name <<endl;

  	}
	
	TDirectory* th1dir = fout->mkdir("Histos1D");
  	th1dir->cd();
	if(information) cout << "[PROCES]	Looping over 1D plots" << endl; 
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
