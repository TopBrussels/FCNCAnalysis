#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TFile.h"
#include <TH1.h>
#include <iostream>
#include <TCanvas.h>
#include <vector>


using namespace std;

void FCNC_CutEfficiencies(string channel = "3L"){


	string rootFileName = "../Output/FCNC_selection_";
	rootFileName += channel;
	rootFileName += ".root";
	
	
	TFile *file = new TFile(rootFileName.c_str(),"read");


	vector <string> Vector_SampleName;
	vector <TH1F*> Efficiency_cutflows;
	
	//Put in the totals for signal and background
	Vector_SampleName.push_back("total_B");
	Vector_SampleName.push_back("total_S");
	
	//Put in the individual samples
	Vector_SampleName.push_back("ttbar");
	Vector_SampleName.push_back("wjets");
	Vector_SampleName.push_back("ttt");
	Vector_SampleName.push_back("ttw");
	Vector_SampleName.push_back("wz");
	Vector_SampleName.push_back("zz");
	Vector_SampleName.push_back("ttz");
	Vector_SampleName.push_back("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToWW_WToLNuL_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToWW_WToLNuL_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToBB_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToBB_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToLL_HctL");
	Vector_SampleName.push_back("TTJetsTocHbW_HToZZ_ZToLL_HctR");


	//loop over all TH1F in the input rootfile and check if the corresponding samples are present in them
	for(int j = 0; j < Vector_SampleName.size(); j++){

		string Path_To_Histo = "Histos1D_cutflows/cutflow_";
		Path_To_Histo += Vector_SampleName[j];
		
		TH1F *Histo( (TH1F*) file->Get(Path_To_Histo.c_str()) );
		
		if(Histo) Efficiency_cutflows.push_back(Histo);
	}
	
	string outputName = "Output/Efficiency_cutflows_";
	outputName += channel;
	outputName += ".root";
	TFile *outputFile = new TFile(outputName.c_str(),"RECREATE"); 
	
	for( int i = 0; i < Efficiency_cutflows.size(); i++){
		int NbOfEvents = Efficiency_cutflows[i]->GetBinContent(2);
		
		Efficiency_cutflows[i]->Scale(1./NbOfEvents);
		Efficiency_cutflows[i]->Write();
	}


}
