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

struct HighestCVSBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
            return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};


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
			cout << "--4L: use the 4 lepton channel (at least 4)" << endl;
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
			xmlfile = argv[iarg];
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
		if (argval=="--Bigxml"){
			Big_xml= true; 
			xmlfile = "../config/FCNC_config.xml";
		
		}
		


    	} 
	

    	if (Big_xml)	xmlfile = "../config/FCNC_config.xml"; 
	if (foundbtag) btagger = tempbtagger;
	
	if(information)	std::cout << "[INFO]	Used configuration file: " << xmlfile << endl;
	if(information)	std::cout << "[INFO]	Used channel: " << channel << endl;
	if(information) std::cout << "[INFO]	Used btag algorithm: " << btagger << endl; 
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
	float Tightworkingpoint = .898;
 
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
	if(channel.find("1L3B")!=string::npos)		sprintf(channelchar, "1L3B");
	if(channel.find("SSdilepton")!=string::npos)	sprintf(channelchar, "SSdilepton");
	if(channel.find("OSdilepton")!=string::npos)	sprintf(channelchar, "OSdilepton");
	
	sprintf(rootFileName,"../data/FCNC_selection_%s.root",channelchar);
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
	MSPlot["NbOfSelectedJets"] = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, -0.5, 14.5, "Nb. of jets");
    	MSPlot["NbOfSelectedLightJets"] = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 15, -0.5, 14.5, "Nb. of jets");
	MSPlot["NbOfSelectedLeptons"] = new MultiSamplePlot(datasets, "NbOfSelectedLeptons", 10, -0.5, 9.5, "Nb. of leptons");
    	MSPlot["NbOfSelectedBJets_CSVM"] = new MultiSamplePlot(datasets, "NbOfSelectedBJets_CSVM", 15, -0.5, 14.5, "Nb. of jets");
    	MSPlot["NbOfSelectedBJets_CSVT"] = new MultiSamplePlot(datasets, "NbOfSelectedBJets_CSVT", 15, -0.5, 14.5, "Nb. of jets");
    	MSPlot["JetEta"] = new MultiSamplePlot(datasets, "JetEta", 30,-3., 3., "Jet #eta");
    	MSPlot["JetPhi"] = new MultiSamplePlot(datasets, "JetPhi", 50, -4., 4., "Jet #phi");
	MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 40, 0., 700., "MET");
	if(channel.find("3L")!=string::npos) MSPlot["MET_3LcH"]= new MultiSamplePlot(datasets, "MET_3LcH", 40, 0., 700., "MET");
	//MSPlot["mll_z"] = new MultiSamplePlot(datasets,"mll_z",50,0,100,"Invariant mass of the leptons that make the Z boson");
	MSPlot["Pt_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_leading_lepton",50,0,100,"Pt leading lepton");
	MSPlot["Pt_2nd_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_lepton",50,0,100,"Pt 2nd leading lepton");
	MSPlot["Pt_3d_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_3d_leading_lepton",50,0,100,"Pt third leading lepton");
	MSPlot["Pt_4th_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_4th_leading_lepton",50,0,100,"Pt fourth leading lepton");
	MSPlot["Pt_5th_leading_lepton"] = new MultiSamplePlot(datasets,"Pt_5th_leading_lepton",50,0,100,"Pt fifth leading lepton");
	MSPlot["Pt_leading_jet"] = new MultiSamplePlot(datasets,"Pt_leading_jet",100,0,200,"Pt leading jet"); 
	MSPlot["Pt_2nd_leading_jet"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_jet",100,0,200,"Pt 2nd leading jet"); 
	MSPlot["Pt_3d_leading_jet"] = new MultiSamplePlot(datasets,"Pt_3d_leading_jet",100,0,200,"Pt third leading jet"); 
	MSPlot["Pt_4th_leading_jet"] = new MultiSamplePlot(datasets,"Pt_4th_leading_jet",100,0,200,"Pt fourth leading jet"); 
	MSPlot["Pt_5th_leading_jet"] = new MultiSamplePlot(datasets,"Pt_5th_leading_jet",100,0,200,"Pt fifth leading jet"); 
	MSPlot["Pt_6th_leading_jet"] = new MultiSamplePlot(datasets,"Pt_6th_leading_jet",100,0,200,"Pt sixth leading jet");
	MSPlot["Pt_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_leading_Bjet",100,0,200,"Pt leading Bjet"); 
	MSPlot["Pt_2nd_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_2nd_leading_Bjet",100,0,200,"Pt 2nd leading Bjet"); 
	MSPlot["Pt_3d_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_3d_leading_Bjet",100,0,200,"Pt third leading Bjet"); 
	MSPlot["Pt_4th_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_4th_leading_Bjet",100,0,200,"Pt fourth leading Bjet"); 
	MSPlot["Pt_5th_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_5th_leading_Bjet",100,0,200,"Pt fifth leading Bjet"); 
	MSPlot["Pt_6th_leading_Bjet"] = new MultiSamplePlot(datasets,"Pt_6th_leading_Bjet",100,0,200,"Pt sixth leading Bjet"); 
	MSPlot["Mll"] = new MultiSamplePlot(datasets,"Mll",50,0,200,"Mll of leading and second leading lepton");
	if(channel.find("3L")!=string::npos) MSPlot["Mll_3LcH"]= new MultiSamplePlot(datasets,"Mll_3LcH",50,0,200,"Mll of leading and second leading lepton 3LcH");
	//MSPlot["Mllq"] = new MultiSamplePlot(datasets,"Mllq",50,0,100,"Invariant mass of llq ~ mtop");
	//if(channel.find("4L")!=string::npos) MSPlot["Mllll"] = new MultiSamplePlot(datasets,"Mllll",50,0,250,"Invariant mass of llll ~ 2mZ");
	MSPlot["DR_toplepton_MET"] = new MultiSamplePlot(datasets,"DR_toplepton_MET",50,0,100,"DR between toplepton and neutrino");
	MSPlot["DR_toplepton_bjet"] = new MultiSamplePlot(datasets,"DR_toplepton_bjet",50,0,100,"DR between toplepton and bjet");
	MSPlot["Mt_toplepton_MET"] = new MultiSamplePlot(datasets,"Mt_toplepton_MET",50,0,100,"Transverse mass of toplepton and neutrino");
	MSPlot["Mt_toplepton_MET_bjet"] = new MultiSamplePlot(datasets,"Mt_toplepton_MET_bjet",50,0,100,"Transverse mass of toplepton, bjet and neutrino");
	MSPlot["Mbqq"] = new MultiSamplePlot(datasets,"Mbqq",50,0,100,"Invariant mass of bqq ~ mtop");
	MSPlot["Mllqq"] = new MultiSamplePlot(datasets,"Mllqq",50,0,100,"Invariant mass of llqq ~ mH");
	MSPlot["Mllqqq"] = new MultiSamplePlot(datasets,"Mllqqq",50,0,100,"Invariant mass of llqqq ~ mtop");
	MSPlot["Mbb"]= new MultiSamplePlot(datasets,"Mbb",50,0,200,"Invariant mass of bb ~ Higgs");
	MSPlot["DeltaPhi_bb"]= new MultiSamplePlot(datasets,"DeltaPhi_bb",30,0,5,"DeltaPhi_bb");
	MSPlot["DR_bb"]= new MultiSamplePlot(datasets,"DR_bb",30,0,5,"DR_bb");
	//////////////////  Cut flow histograms	/////////////////////////////
	MSPlot["MScutflow"] = new MultiSamplePlot(datasets,"MScutflow",20,-0.5,19.5, "cutflow"); 
	if(channel.find("4L")!=string::npos) MSPlot["NbOfJets_4L4"]= new MultiSamplePlot(datasets,"NbOfJets_4L4",50,0,250,"#jets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["NbOfBJets_4L4_CSVM"]= new MultiSamplePlot(datasets,"NbOfBJets_4L4_CSVM",50,0,250,"#Bjets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["NbOfBJets_4L4_CSVT"]= new MultiSamplePlot(datasets,"NbOfBJets_4L4_CSVT",50,0,250,"#Bjets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["MET_4L4"]= new MultiSamplePlot(datasets, "MET_4L4", 40, 0., 700., "MET");
	if(channel.find("4L")!=string::npos) MSPlot["NbOfJets_4L5"]= new MultiSamplePlot(datasets,"NbOfJets_4L5",50,0,250,"#jets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["NbOfBJets_4L5_CSVM"]= new MultiSamplePlot(datasets,"NbOfBJets_4L5_CSVM",50,0,250,"#Bjets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["NbOfBJets_4L5_CSVT"]= new MultiSamplePlot(datasets,"NbOfBJets_4L5_CSVT",50,0,250,"#Bjets for exactly 4 leptons");
	if(channel.find("4L")!=string::npos) MSPlot["MET_4L5"] = new MultiSamplePlot(datasets, "MET_4L5", 40, 0., 700., "MET");
	
	
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

	

	// Define different cutflow plots for each channel and dataset	
	for(unsigned int d = 0; d < datasets.size();d++){ 
		//Load datasets
		treeLoader.LoadDataset(datasets[d], anaEnv); 
		string datasetName = datasets[d]->Name(); 
		
		char datasetNamechar[900];
		if(datasetName.find("W_1Jets")!=string::npos) {sprintf(datasetNamechar,"W_1Jets");}
		if(datasetName.find("W_2Jets")!=string::npos) {sprintf(datasetNamechar,"W_2Jets");}
		if(datasetName.find("W_3Jets")!=string::npos) {sprintf(datasetNamechar,"W_3Jets");}
		if(datasetName.find("W_4Jets")!=string::npos) {sprintf(datasetNamechar,"W_4Jets");}
		if(datasetName.find("WW_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"WW_To2L2Nu");}
		if(datasetName.find("WZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"WZ_To2L2Q");}
		if(datasetName.find("WZ_To3LNu")!=string::npos) {sprintf(datasetNamechar,"WZ_To3LNu");}
		if(datasetName.find("ZZ_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Nu");}
		if(datasetName.find("ZZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Q");}
		if(datasetName.find("ZZ_To4L")!=string::npos) {sprintf(datasetNamechar,"ZZ_To4L");}
		if(datasetName.find("ST_T_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_t-ch");}
		if(datasetName.find("ST_Tbar_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_Tbar_t-ch");}
		if(datasetName.find("ST_TToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToDilepton_tW-ch");}
		if(datasetName.find("ST_TToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToThadWlep_tW-ch");}
		if(datasetName.find("ST_TBarToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToDilepton_tW-ch");}
		if(datasetName.find("ST_TBarToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TBarToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToThadWlep_tW-ch");}
		if(datasetName.find("TT_SemiLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_SemiLeptMGDecays");}
		if(datasetName.find("TT_FullLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_FullLeptMGDecays");}
		if(datasetName.find("TT_HadronicMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_HadronicMGDecays");}
		if(datasetName.find("Z_M-10To50")!=string::npos) {sprintf(datasetNamechar,"Z_M-10To50");}
		if(datasetName.find("Z_M-50")!=string::npos) {sprintf(datasetNamechar,"Z_M-50");}
		if(datasetName.find("Z_1Jets")!=string::npos) {sprintf(datasetNamechar,"Z_1Jets");}
		if(datasetName.find("Z_2Jets")!=string::npos) {sprintf(datasetNamechar,"Z_2Jets");}
		if(datasetName.find("Z_3Jets")!=string::npos) {sprintf(datasetNamechar,"Z_3Jets");}
		if(datasetName.find("Z_4Jets")!=string::npos) {sprintf(datasetNamechar,"Z_4Jets");}
		if(datasetName.find("TTZ")!=string::npos) {sprintf(datasetNamechar,"TTZ");}
		if(datasetName.find("TTW")!=string::npos) {sprintf(datasetNamechar,"TTW");}
		if(datasetName.find("ttbar")!=string::npos) {sprintf(datasetNamechar,"ttbar");}
                if(datasetName.find("ttbar_fullLept")!=string::npos) {sprintf(datasetNamechar,"ttbar_fullLept");}
                if(datasetName.find("ttbar_semiLept")!=string::npos) {sprintf(datasetNamechar,"ttbar_semiLept");}
                if(datasetName.find("Wjets")!=string::npos || datasetName.find("wjets")!=string::npos) {sprintf(datasetNamechar,"wjets");}
                if(datasetName.find("ttt")!=string::npos) {sprintf(datasetNamechar,"ttt");}
                if(datasetName.find("ttW")!=string::npos) {sprintf(datasetNamechar,"ttw");}
                if(datasetName.find("WW")!=string::npos || datasetName.find("ww")!=string::npos) {sprintf(datasetNamechar,"ww");}
                if(datasetName.find("WZ")!=string::npos || datasetName.find("wz")!=string::npos) {sprintf(datasetNamechar,"wz");}
                if(datasetName.find("ZZ")!=string::npos || datasetName.find("zz")!=string::npos) {sprintf(datasetNamechar,"zz");}
                if(datasetName.find("ttZ")!=string::npos || datasetName.find("ttz")!=string::npos) {sprintf(datasetNamechar,"ttz");}
                if(datasetName.find("Zjets")!=string::npos || datasetName.find("zjets")!=string::npos) {sprintf(datasetNamechar,"Zjets");}
                if(datasetName.find("ST_T_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_tW-ch");}
                if(datasetName.find("ST_TBar_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_tW-ch");}
		if(datasetName.find("ST_T_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_s-ch");}
                if(datasetName.find("ST_TBar_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_s-ch");}

		
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
		
		histo1D[Process_cutflow] = new TH1F(NamePlot, plotTitle, 15, -0.5,14.5);
		//histo1D[Process_cutflow]->Sumw2();
		histo1D[Process_cutflow]->GetYaxis()->SetTitle("#evts.");

               
	
	}
	
	
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
		if(datasetName.find("W_1Jets")!=string::npos) {sprintf(datasetNamechar,"W_1Jets");}
		if(datasetName.find("W_2Jets")!=string::npos) {sprintf(datasetNamechar,"W_2Jets");}
		if(datasetName.find("W_3Jets")!=string::npos) {sprintf(datasetNamechar,"W_3Jets");}
		if(datasetName.find("W_4Jets")!=string::npos) {sprintf(datasetNamechar,"W_4Jets");}
		if(datasetName.find("WW_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"WW_To2L2Nu");}
		if(datasetName.find("WZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"WZ_To2L2Q");}
		if(datasetName.find("WZ_To3LNu")!=string::npos) {sprintf(datasetNamechar,"WZ_To3LNu");}
		if(datasetName.find("ZZ_To2L2Nu")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Nu");}
		if(datasetName.find("ZZ_To2L2Q")!=string::npos) {sprintf(datasetNamechar,"ZZ_To2L2Q");}
		if(datasetName.find("ZZ_To4L")!=string::npos) {sprintf(datasetNamechar,"ZZ_To4L");}
		if(datasetName.find("ST_T_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_t-ch");}
		if(datasetName.find("ST_Tbar_t-ch")!=string::npos) {sprintf(datasetNamechar,"ST_Tbar_t-ch");}
		if(datasetName.find("ST_TToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToDilepton_tW-ch");}
		if(datasetName.find("ST_TToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TToThadWlep_tW-ch");}
		if(datasetName.find("ST_TBarToDilepton_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToDilepton_tW-ch");}
		if(datasetName.find("ST_TBarToTlepWhad_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToTlepWhad_tW-ch");}
		if(datasetName.find("ST_TBarToThadWlep_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBarToThadWlep_tW-ch");}
		if(datasetName.find("TT_SemiLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_SemiLeptMGDecays");}
		if(datasetName.find("TT_FullLeptMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_FullLeptMGDecays");}
		if(datasetName.find("TT_HadronicMGDecays")!=string::npos) {sprintf(datasetNamechar,"TT_HadronicMGDecays");}
		if(datasetName.find("Z_M-10To50")!=string::npos) {sprintf(datasetNamechar,"Z_M-10To50");}
		if(datasetName.find("Z_M-50")!=string::npos) {sprintf(datasetNamechar,"Z_M-50");}
		if(datasetName.find("Z_1Jets")!=string::npos) {sprintf(datasetNamechar,"Z_1Jets");}
		if(datasetName.find("Z_2Jets")!=string::npos) {sprintf(datasetNamechar,"Z_2Jets");}
		if(datasetName.find("Z_3Jets")!=string::npos) {sprintf(datasetNamechar,"Z_3Jets");}
		if(datasetName.find("Z_4Jets")!=string::npos) {sprintf(datasetNamechar,"Z_4Jets");}
		if(datasetName.find("TTZ")!=string::npos) {sprintf(datasetNamechar,"TTZ");}
		if(datasetName.find("TTW")!=string::npos) {sprintf(datasetNamechar,"TTW");}
		if(datasetName.find("ttbar")!=string::npos) {sprintf(datasetNamechar,"ttbar");}
                if(datasetName.find("ttbar_fullLept")!=string::npos) {sprintf(datasetNamechar,"ttbar_fullLept");}
                if(datasetName.find("ttbar_semiLept")!=string::npos) {sprintf(datasetNamechar,"ttbar_semiLept");}
                if(datasetName.find("Wjets")!=string::npos || datasetName.find("wjets")!=string::npos) {sprintf(datasetNamechar,"wjets");}
                if(datasetName.find("ttt")!=string::npos) {sprintf(datasetNamechar,"ttt");}
                if(datasetName.find("ttW")!=string::npos) {sprintf(datasetNamechar,"ttw");}
                if(datasetName.find("WW")!=string::npos || datasetName.find("ww")!=string::npos) {sprintf(datasetNamechar,"ww");}
                if(datasetName.find("WZ")!=string::npos || datasetName.find("wz")!=string::npos) {sprintf(datasetNamechar,"wz");}
                if(datasetName.find("ZZ")!=string::npos || datasetName.find("zz")!=string::npos) {sprintf(datasetNamechar,"zz");}
                if(datasetName.find("ttZ")!=string::npos || datasetName.find("ttz")!=string::npos) {sprintf(datasetNamechar,"ttz");}
                if(datasetName.find("Zjets")!=string::npos || datasetName.find("zjets")!=string::npos) {sprintf(datasetNamechar,"Zjets");}
                if(datasetName.find("ST_T_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_tW-ch");}
                if(datasetName.find("ST_TBar_tW-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_tW-ch");}
		if(datasetName.find("ST_T_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_T_s-ch");}
                if(datasetName.find("ST_TBar_s-ch")!=string::npos) {sprintf(datasetNamechar,"ST_TBar_s-ch");}
		
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
		
		//def
		string Process_cutflow = "cutflow_";
		Process_cutflow += datasetNamechar;
		
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
			histo1D[Process_cutflow]->Fill(1);
			if(!is_signal) histo1D["cutflow_total_B"]->Fill(1);
			if(is_signal) histo1D["cutflow_total_S"]->Fill(1);
			
			if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(2, "initial");
			if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(2, "initial");
			histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(2, "initial");
			
			
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
			
			//Start btagging 
			int nTags = 0;
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
					nTags++;
					selectedBJets_CSVM.push_back(selectedJets[iJet]);
				}
				else selectedLightJets.push_back(selectedJets[iJet]);
				
				if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > Tightworkingpoint)
				{
					selectedBJets_CSVT.push_back(selectedJets[iJet]);
				}
				
			}

		sort(selectedBJets_CSVM.begin(), selectedBJets_CSVM.end(), HighestCVSBtag());

			
			if(debug) cout << "[INFO]	looseElectrons.size() = " << looseElectrons.size() << endl; 
			if(debug) cout << "[INFO]	looseMuons.size() = " << looseMuons.size() << endl; 
			
			
			bool OneLepton_4Jets = false;
			bool chan3L4L = false; 
			
			
			//exactly 3 leptons
			if(channel.find("3L")!=string::npos)
			{
				chan3L4L = true;
				if(debug) cout << "[PROCES]	in 3L channel" << endl;
				if(looseElectrons.size() + looseMuons.size() ==3)
				{ 
					if(debug) cout << "[PROCES]	fill 3L" << endl;
					
					//fill histograms
					MSPlot["MScutflow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
					if(!is_signal)histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					
					//set labels
					if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "3L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "3L");		
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "3L");
					
					Passed_selection = true;
					
					
					if(debug) cout << "[PROCES]	filled 3L" << endl;
					
					
				}
				if(debug)	cout << "[PROCES]	out fill 3L loop" << endl; 
			}
			//more than 4 leptons
			if(channel.find("4L")!=string::npos)
			{
				chan3L4L = true;
				if(debug) cout << "[PROCES]	in 4L channel" << endl;
				
				if(looseElectrons.size() + looseMuons.size() > 3)
				{ 
					if(debug) cout << "[PROCES]	fill 4L" << endl;
					
					//fill histograms
					if(!is_signal)	histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					MSPlot["MScutflow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
					//label histograms
					if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "4L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "4L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "4L");
					Passed_selection = true;
									
					if(debug) cout << "[PROCES]	filled 4L" << endl;
					
					
				}
				if(debug)	cout << "[PROCES]	out fill 4L loop" << endl; 
			}
			//1 lepton + 3 b-jets
			if(channel.find("1L3B")!=string::npos)
			{
				if(debug) cout << "in 1L3B channel" << endl;
				
				if(looseElectrons.size() +  looseMuons.size() == 1)
				{
					if(debug) cout << "in fill 1l3b loop" << endl;
					
					//fill histograms
					MSPlot["MScutflow"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
					if(!is_signal) histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					//label histograms
					if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "1L");
					if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(3, "1L");
					histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(3, "1L");
					
					if(debug) cout << "selectedJets.size() = " << selectedJets.size() << endl;
					
					if(selectedJets.size() >= 4)
					{
						OneLepton_4Jets = true;
						if(debug) cout << "in fill 1l3b loop: 3jets" << endl;
						//fill histograms
						MSPlot["MScutflow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
						if(!is_signal) histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						//label histograms
						if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, ">= 4jets");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, ">= 4jets");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, ">=4jets");
						
						if(nTags == 3)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							//fill histograms
							MSPlot["MScutflow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
							if(!is_signal) histo1D["cutflow_total_B"]->Fill(4);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(4);
							histo1D[Process_cutflow]->Fill(4);
							//label histograms
							if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(5, "== 3 bjets");
							
							
						}
						if(nTags > 0)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							//fill histograms
							MSPlot["MScutflow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
							if(!is_signal) histo1D["cutflow_total_B"]->Fill(5);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(5);
							histo1D[Process_cutflow]->Fill(5);
							//label histograms
							if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(6, ">= 1 bjets");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(6, ">= 1 bjets");
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(6, ">= 1 bjets");
							
						}
						if(nTags > 1)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							//fill histograms
							MSPlot["MScutflow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
							if(!is_signal) histo1D["cutflow_total_B"]->Fill(6);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(6);
							histo1D[Process_cutflow]->Fill(6);
							//label histograms
							if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(7, ">= 2 bjets");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(7, ">= 2 bjets");
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(7, ">= 2 bjets");
							
						}
						if(nTags > 2)
						{
							if(debug) cout << "in fill 1l3b loop: 3bjets" << endl;
							//fill histograms
							MSPlot["MScutflow"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
							if(!is_signal) histo1D["cutflow_total_B"]->Fill(7);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(7);
							histo1D[Process_cutflow]->Fill(7);
							//label histograms
							if(!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(8, ">= 3 bjets");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(8, ">= 3 bjets");
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(8, ">= 3 bjets");
							
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
					
					if(!is_signal)histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					
					if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "2L");
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
						
						if(!is_signal)histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						
						if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, "2 SS L");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, "2 SS L");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, "2 SS L");
						
						if(selectedJets.size()>=1)
						{
							//fill histograms
							if(!is_signal)histo1D["cutflow_total_B"]->Fill(4);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(4);
							histo1D[Process_cutflow]->Fill(4);
					
							//set labels
							if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(5,">=1j");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(5,">=1j");		
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(5, ">=1j");
							Passed_selection = true;
						}						
					}
					if(debug) cout << "out fill SS dilepton " << endl;
				}

			}
			if(channel.find("OSdilepton")!=string::npos)
			{
			float Zmass = 91.1876;  // ref-> pdg
		                float massDiff = 9999;
		                TLorentzVector leptonPairMass;
	                	bool mll = false;
				if(debug) cout << "in OSdilepton channel" << endl;
				if(looseElectrons.size() + looseMuons.size() == 2)
				{
					if(debug) cout << "in fill OS dilepton " << endl; 
					if(!is_signal)histo1D["cutflow_total_B"]->Fill(2);
					if(is_signal) histo1D["cutflow_total_S"]->Fill(2);
					histo1D[Process_cutflow]->Fill(2);
					
					if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(3, "2L");
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
						if(!is_signal) histo1D["cutflow_total_B"]->Fill(3);
						if(is_signal) histo1D["cutflow_total_S"]->Fill(3);
						histo1D[Process_cutflow]->Fill(3);
						if (!is_signal) histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(4, "2 OS L");
						if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(4, "2 OS L");
						histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(4, "2 OS L");
						
						Passed_selection = true;
						// Additional cuts
						if (selectedJets.size()>=4){
						
						if (electron) {
						leptonPairMass = looseElectrons[0]->M() + looseElectrons[1]->M();
						mll = true;}
						if (muon){
						leptonPairMass = looseMuons[0]->M() + looseMuons[1]->M();
						mll =true;}
						if (EMu){      
						leptonPairMass = looseMuons[0]->M() + looseElectrons[0]->M();
						mll = true;
						}
						if (mll){
						massDiff = leptonPairMass.M() - Zmass ;
						if (massDiff > 15)
						{
						MSPlot["massDiff"]->Fill(massDiff, datasets[d],true,Luminosity*scaleFactor);
						}
						}}
						 	
					}
					
						if(selectedJets.size()>=1  && Passed_selection)
						{
							//fill histograms
							if(!is_signal)histo1D["cutflow_total_B"]->Fill(4);
							if(is_signal) histo1D["cutflow_total_S"]->Fill(4);
							histo1D[Process_cutflow]->Fill(4);
					
							//set labels
							if(!is_signal)histo1D["cutflow_total_B"]->GetXaxis()->SetBinLabel(5, ">1j");
							if(is_signal) histo1D["cutflow_total_S"]->GetXaxis()->SetBinLabel(5, ">1j");		
							histo1D[Process_cutflow]->GetXaxis()->SetBinLabel(5, ">=1j");
							Passed_selection = true;
						}		
					if(debug) cout << "out fill OS dilepton " << endl;
				}

			}

			if(channel.find("1gamma")!=string::npos)
			{
			}
			if(channel.find("2gamma")!=string::npos)
			{
			}                                                                
    			
	

			//////////////////////////////////////////////////////////////////////////////////
			// Filling histograms 							//////////
			//////////////////////////////////////////////////////////////////////////////////
			int nLeptons = looseElectrons.size() +  looseMuons.size();
			if(nLeptons == 1){
				MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbOfSelectedBJets_CSVM"]->Fill(selectedBJets_CSVM.size(), datasets[d], true, Luminosity*scaleFactor);
				MSPlot["NbOfSelectedBJets_CSVT"]->Fill(selectedBJets_CSVT.size(), datasets[d], true, Luminosity*scaleFactor);
				
			}
			if(!Passed_selection) continue;
			if(Passed_selection){
				if(debug) cout << "[PROCES]	In passed_selection loop" << endl; 
				
				if(nLeptons != 1) MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
				if(debug) cout << "[PROCES]	Filled NbOfSelectedJets" << endl; 
				MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
		        	if(debug) cout << "[PROCES]	Filled NbOfSelectedLightJets" << endl; 
				if(nLeptons != 1) MSPlot["NbOfSelectedBJets_CSVM"]->Fill(selectedBJets_CSVM.size(), datasets[d], true, Luminosity*scaleFactor);
		        	if(debug) cout << "[PROCES]	Filled NbOfSelectedBJets_CSVM" << endl; 
				if(nLeptons != 1) MSPlot["NbOfSelectedBJets_CSVT"]->Fill(selectedBJets_CSVT.size(), datasets[d], true, Luminosity*scaleFactor);
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
                                	if(chan3L4L) MSPlot["MScutflow"]->Fill(3, datasets[d], true, Luminosity*scaleFactor);
					if(debug) cout << "[PROCES]	Out selectedJets.size() > 0" << endl; 
				}
				if( selectedJets.size() > 1)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(4, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_2nd_leading_jet"]->Fill(selectedJets[1]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedJets.size() > 2)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(5, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_3d_leading_jet"]->Fill(selectedJets[2]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedJets.size() > 3)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(6, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_4th_leading_jet"]->Fill(selectedJets[3]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedJets.size() > 4)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(7, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_5th_leading_jet"]->Fill(selectedJets[4]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedJets.size() > 5)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(8, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_6th_leading_jet"]->Fill(selectedJets[5]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 0) {
					if(chan3L4L) MSPlot["MScutflow"]->Fill(9, datasets[d], true, Luminosity*scaleFactor);
                                	MSPlot["Pt_leading_Bjet"]->Fill(selectedBJets_CSVM[0]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
                                }
				if( selectedBJets_CSVM.size() > 1)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(10, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_2nd_leading_Bjet"]->Fill(selectedBJets_CSVM[1]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 2)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(11, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_3d_leading_Bjet"]->Fill(selectedBJets_CSVM[2]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 3)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(12, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_4th_leading_Bjet"]->Fill(selectedBJets_CSVM[3]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 4)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(13, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_5th_leading_Bjet"]->Fill(selectedBJets_CSVM[4]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				if( selectedBJets_CSVM.size() > 5)
				{
					if(chan3L4L) MSPlot["MScutflow"]->Fill(14, datasets[d], true, Luminosity*scaleFactor);
				  	MSPlot["Pt_6th_leading_Bjet"]->Fill(selectedBJets_CSVM[5]->Pt(), datasets[d],true,	Luminosity*scaleFactor);
				}
				
				
				double mll = 0; 
				TLorentzVector leptonpair_mll;
				TLorentzVector lepton0;
				TLorentzVector lepton1;
				leptonpair_mll.Clear(); 
				lepton0.Clear(); 
				lepton1.Clear(); 
				if(looseElectrons.size()==0)
				{
					if(debug) cout << "[PROCES]	in looseElectrons.size()==0" << endl; 
					if(looseMuons.size()!=0)
					{
						lepton0.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
						
						
						if(looseMuons.size()>1 )
						{
							lepton1.SetPxPyPzE(looseMuons[1]->Px(),looseMuons[1]->Py(),looseMuons[1]->Pz(),looseMuons[1]->Energy());
							leptonpair_mll = lepton0+lepton1;
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
							leptonpair_mll = lepton0+lepton1;
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
					leptonpair_mll = lepton0+lepton1;
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
					leptonpair_mll = lepton0+lepton1;
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
					leptonpair_mll = lepton0+lepton1;
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				else if(looseMuons.size() == 1 && looseElectrons.size() == 1)
				{
					if(debug) cout << "[PROCES]	in looseMuons.size() == 1 && looseElectrons.size() == 1" << endl; 
					lepton0.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->Energy());
					lepton1.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->Energy());
					leptonpair_mll = lepton0+lepton1;
					MSPlot["Pt_leading_lepton"]->Fill(lepton0.Pt(),datasets[d], true, Luminosity*scaleFactor);
					MSPlot["Pt_2nd_leading_lepton"]->Fill(lepton1.Pt(),datasets[d], true, Luminosity*scaleFactor);
				}
				bool empty = false; 
				if(looseMuons.size() == 0 && looseElectrons.size() == 0) empty = true; 
				if(looseMuons.size() == 0 && looseElectrons.size() == 1) empty = true; 
				if(looseMuons.size() == 1 && looseElectrons.size() == 0) empty = true; 
				if(!empty && (leptonpair_mll!=(0,0,0,0)))
				{
					mll = leptonpair_mll.M();
					if(debug) cout << "[INFO]	mll = 	" << mll << endl; 
					MSPlot["Mll"]->Fill(mll,datasets[d],true,Luminosity*scaleFactor);
					if(channel.find("3L")!=string::npos)
			 		{
			 			if(selectedBJets_CSVM.size()>1)
						{
							if(debug) cout << "[PROCES]	In 3L (#bjets > 1)" <<  endl; 
							MSPlot["Mll_3LcH"]->Fill(mll,datasets[d],true,Luminosity*scaleFactor);
							MSPlot["MET_3LcH"]->Fill((float) met_pt, datasets[d], true, Luminosity*scaleFactor);
						}
			 		}
					
				}
				
				
				
				/*
				double mllll = 0; 
				TLorentzVector leptonpair_mllll;
				TLorentzVector lepton2;
				TLorentzVector lepton3;
				leptonpair_mllll.Clear(); 
				lepton2.Clear(); 
				lepton3.Clear(); 
				if(looseMuons.size()>3 && looseElectrons.size() == 0)
				{
					lepton2.SetPxPyPzE(looseMuons[2]->Px(),looseMuons[2]->Py(),looseMuons[2]->Pz(),looseMuons[2]->Energy());
					lepton3.SetPxPyPzE(looseMuons[3]->Px(),looseMuons[3]->Py(),looseMuons[3]->Pz(),looseMuons[3]->Energy());
					leptonpair_mllll = leptonpair_mll + lepton2 + lepton3;
					mllll = leptonpair_mllll.M(); 
				}
				if(looseMuons.size()==0&& looseElectrons.size() > 3)
				{
					lepton2.SetPxPyPzE(looseElectrons[2]->Px(),looseElectrons[2]->Py(),looseElectrons[2]->Pz(),looseElectrons[2]->Energy());
					lepton3.SetPxPyPzE(looseElectrons[3]->Px(),looseElectrons[3]->Py(),looseElectrons[3]->Pz(),looseElectrons[3]->Energy());
					leptonpair_mllll = leptonpair_mll + lepton2 + lepton3;
					mllll = leptonpair_mllll.M(); 
				}
				
				
				
				
				
				
				if(mllll != 0 && channel.find("4L")!=string::npos ) MSPlot["Mllll"]->Fill(mllll,datasets[d],true,Luminosity*scaleFactor);
				*/
				
			}
			
			
			
			

			
			// adding additional cuts
			 if(channel.find("4L")!=string::npos && Passed_selection)
			 {
			 	
			 	if(looseMuons.size()+looseElectrons.size() < 6)
				{
					// 3 < #leptons < 6
					MSPlot["MScutflow"]->Fill(15, datasets[d], true, Luminosity*scaleFactor);
					
				
				}
				if(looseMuons.size()+looseElectrons.size() < 5)
				{
					// #lepons = 4
					MSPlot["MScutflow"]->Fill(16, datasets[d], true, Luminosity*scaleFactor);
					MSPlot["NbOfBJets_4L4_CSVM"]->Fill(selectedBJets_CSVM.size(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["NbOfBJets_4L4_CSVT"]->Fill(selectedBJets_CSVT.size(),datasets[d],true,Luminosity*scaleFactor);
					// #bjets should be 1 + Z(bb) = 3
					MSPlot["NbOfJets_4L4"]->Fill(selectedJets.size(),datasets[d],true,Luminosity*scaleFactor);  
					// #jets should be 4
					MSPlot["MET_4L4"]->Fill((float) met_pt, datasets[d], true, Luminosity*scaleFactor);
					// MET should be zero
				}
				if(looseMuons.size()+looseElectrons.size() == 5)
				{
					// #leptons = 5
					MSPlot["MScutflow"]->Fill(17, datasets[d], true, Luminosity*scaleFactor);
					MSPlot["NbOfBJets_4L5_CSVM"]->Fill(selectedBJets_CSVM.size(),datasets[d],true,Luminosity*scaleFactor);
					MSPlot["NbOfBJets_4L5_CSVT"]->Fill(selectedBJets_CSVT.size(),datasets[d],true,Luminosity*scaleFactor);
					// #bjets should be 1 
					MSPlot["NbOfJets_4L5"]->Fill(selectedJets.size(),datasets[d],true,Luminosity*scaleFactor);  
					// #jets should be 2
					MSPlot["MET_4L5"]->Fill((float) met_pt, datasets[d], true, Luminosity*scaleFactor);
					// MET should be half the mass of a W boson
				
				}
			 	
			 
			 }
			 if(channel.find("3L")!=string::npos && Passed_selection)
			 {
			 		
			 }
			
				
			if(channel.find("1L3B")!=string::npos && Passed_selection){
					TLorentzVector Lepton;
					TLorentzVector bb_cand;

					Lepton.Clear();
					bb_cand.Clear();
					
					vector <float> frac_b_topcandidate;
					if (looseElectrons.size() == 1) Lepton.SetPxPyPzE(looseElectrons[0]->Px(),looseElectrons[0]->Py(),looseElectrons[0]->Pz(),looseElectrons[0]->E());
					if (looseMuons.size() == 1) Lepton.SetPxPyPzE(looseMuons[0]->Px(),looseMuons[0]->Py(),looseMuons[0]->Pz(),looseMuons[0]->E());

				
					float TransvM_lept_MET = 0;
					float TransvM_lept_MET_b = 0;
					float InvM_bb = 0;
					float DeltaR_lepton_b_min = 9999; //The b which is closest to the lepton is more probable to come from the SM top decay
					float DeltaR_bb_min = 9999; 
					float DeltaR_bb_max = -9999; //The bigger DeltaR_bb, the more probable those 2 b-jets do not come from H->bb (and therefore the max is a good handle to determine which of the 3 b-jets comes from SM top decay
					float dR_lb_temp = 9999;
					float dR_bb_temp = 9999;
					int index_BTopcandidate = 9999; // based on the minimal DeltaR_lepton_b + DeltaR_bb
					
					//Determine the minimal Delta R
					for(int iBjet = 0; iBjet < selectedBJets_CSVM.size(); iBjet++){
						
						dR_lb_temp = sqrt(pow(Lepton.Eta() - selectedBJets_CSVM[iBjet]->Eta(),2)+pow(Lepton.Phi() - selectedBJets_CSVM[iBjet]->Phi(),2));
						if(DeltaR_lepton_b_min>dR_lb_temp){
							DeltaR_lepton_b_min = dR_lb_temp;
							//index_BTopcandidate = iBjet;
						}
						
						for(int iBjet2 = 0; iBjet2 < selectedBJets_CSVM.size(); iBjet2++){
							if(iBjet != iBjet2) dR_bb_temp = sqrt(pow(selectedBJets_CSVM[iBjet]->Eta() - selectedBJets_CSVM[iBjet2]->Eta(),2)+pow(selectedBJets_CSVM[iBjet]->Phi() - selectedBJets_CSVM[iBjet2]->Phi(),2));
							if(DeltaR_bb_min>dR_bb_temp) DeltaR_bb_min = dR_bb_temp;
							if(DeltaR_bb_max<dR_bb_temp) DeltaR_bb_max = dR_bb_temp;
						}
						
						frac_b_topcandidate.push_back(DeltaR_lepton_b_min/DeltaR_bb_max);
					}
					
					//Get the index of the b-jet which most probably comes from the top quark through the fraction DeltaR_lepton_b_min/DeltaR_bb_max
					float min_frac = 999999;
					for(int iBjet = 0; iBjet < frac_b_topcandidate.size(); iBjet++){
						if(min_frac > frac_b_topcandidate[iBjet]){
							min_frac = frac_b_topcandidate[iBjet];
							index_BTopcandidate = iBjet;
						}
					}
					
					
					
					
					float DeltaPhi_bb = 0;
					float InvMass_bb = 0;
					for(int iBjet = 0; iBjet < selectedBJets_CSVM.size(); iBjet++){
						for(int iBjet2 = 0; iBjet2 < selectedBJets_CSVM.size(); iBjet2++){
							if((iBjet != index_BTopcandidate) && (iBjet2 != index_BTopcandidate) && (iBjet != iBjet2)){
								
								bb_cand.SetPxPyPzE(selectedBJets_CSVM[iBjet]->Px()+selectedBJets_CSVM[iBjet2]->Px(),selectedBJets_CSVM[iBjet]->Py()+selectedBJets_CSVM[iBjet2]->Py(),selectedBJets_CSVM[iBjet]->Pz()+selectedBJets_CSVM[iBjet2]->Pz(),selectedBJets_CSVM[iBjet]->E()+selectedBJets_CSVM[iBjet]->E());
								
								DeltaPhi_bb = sqrt( pow(( selectedBJets_CSVM[iBjet]->Phi() - selectedBJets_CSVM[iBjet2]->Phi() ), 2));
								InvMass_bb = bb_cand.M();
							}
						}
					}
					
					
				MSPlot["DR_toplepton_bjet"]->Fill(DeltaR_lepton_b_min, datasets[d], true, Luminosity*scaleFactor);
				MSPlot["DR_bb"]->Fill(DeltaR_bb_min, datasets[d],true,Luminosity*scaleFactor);
				MSPlot["Mbb"]->Fill(InvMass_bb, datasets[d],true,Luminosity*scaleFactor);
				MSPlot["DeltaPhi_bb"]->Fill(DeltaPhi_bb, datasets[d],true,Luminosity*scaleFactor);
				MSPlot["Pt_leading_lepton"]->Fill(Lepton.Pt(), datasets[d], true, Luminosity*scaleFactor);
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
