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

void FCNC_CutEfficiencies(string channel = "45"){
	cout << "------------------------------------------------------------------" << endl;
	cout << "         [LOOPING OVER " << channel << " CHANNEL]" << endl; 
	
	string	Vector_cutflow_3L[19] = {"","initial","3L","3L & MET>45 GeV","3L & MET < 45 GeV","3L & >3 jets","3L & <4 jets","","","","","","","","","","", "",""};
	string Vector_cutflow_4L[19] = {"","initial","4L","4L & Z(e/mu)","4L & Z(e/mu) & |mllll - mh|<10GeV","4L & Z(e/mu) &|mllll - mh|<10GeV & >3 jets","4L & Z(e/mu) & |mllll - mh|<5 GeV","4L & Z(e/mu) & |mllll - mh|<5 GeV & >3jets","","","","","","","","","","",""};
	string Vector_cutflow_5L[19] = {"","initial","5L","5L & >1 jets","5L & >1 jets & >0 bjets","","","","","","","","","","","","","",""};
	string Vector_cutflow_45[19] = {"","initial",">3L",">3L & >0 jets",">3L & >0 jets & >0 bjets","","","","","","","","","","","","","",""};
	
	
	

	string rootFileName = "../data/FCNC_selection_";
        if(channel.find("45")!=string::npos)
        { rootFileName += "4L5L";
        }
        else
        {	
           rootFileName += channel;
	}
        rootFileName += ".root";
	
	
	TFile *file = new TFile(rootFileName.c_str(),"read");


	vector <string> Vector_SampleName;
	vector <TH1F*> Efficiency_cutflows;
	
	//Put in the totals for signal and background
	Vector_SampleName.push_back("total_B");
	Vector_SampleName.push_back("total_S");
	
	//Put in the individual samples
	Vector_SampleName.push_back("GluGluHiggs4lep");
	Vector_SampleName.push_back("Wjets");
	Vector_SampleName.push_back("W_1Jets");
	Vector_SampleName.push_back("W_2Jets");
	Vector_SampleName.push_back("W_3Jets");
	Vector_SampleName.push_back("W_4Jets");
	
	Vector_SampleName.push_back("WW");
        Vector_SampleName.push_back("WZ");
        Vector_SampleName.push_back("ZZ");
	Vector_SampleName.push_back("WW_To2L2Nu");
	Vector_SampleName.push_back("WZ_To2L2Q");
	Vector_SampleName.push_back("WZ_To3LNu");
	Vector_SampleName.push_back("ZZ_To2L2Nu");
	Vector_SampleName.push_back("ZZ_To2L2Q");
	Vector_SampleName.push_back("ZZ_To4L");
	
	Vector_SampleName.push_back("ST_T_s-ch");
        Vector_SampleName.push_back("ST_TBar_s-ch");
	Vector_SampleName.push_back("ST_T_t-ch");
	Vector_SampleName.push_back("ST_Tbar_t-ch");
	Vector_SampleName.push_back("ST_T_tW-ch");
        Vector_SampleName.push_back("ST_TBar_tW-ch");
	Vector_SampleName.push_back("ST_TToDilepton_tW-ch");
	Vector_SampleName.push_back("ST_TToTlepWhad_tW-ch");
	Vector_SampleName.push_back("ST_TToThadWlep_tW-ch");
	Vector_SampleName.push_back("ST_TBarToDilepton_tW-ch");
	Vector_SampleName.push_back("ST_TBarToTlepWhad_tW-ch");
	Vector_SampleName.push_back("ST_TBarToThadWlep_tW-ch");
	
	Vector_SampleName.push_back("TT_SemiLeptMGDecays");
	Vector_SampleName.push_back("TT_FullLeptMGDecays");
	Vector_SampleName.push_back("TT_HadronicMGDecays");
	Vector_SampleName.push_back("TTJets");
	
	Vector_SampleName.push_back("TTZ");
	Vector_SampleName.push_back("TTW");
	
	Vector_SampleName.push_back("Z_M-10To50");
	Vector_SampleName.push_back("Z_M-50");
	Vector_SampleName.push_back("Z_1Jets");
	Vector_SampleName.push_back("Z_2Jets");
	Vector_SampleName.push_back("Z_3Jets");
	Vector_SampleName.push_back("Z_4Jets");
	
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
	Vector_SampleName.push_back("TTJetsTocZbW");



	//loop over all TH1F in the input rootfile and check if the corresponding samples are present in them
	for(unsigned int j = 0; j < Vector_SampleName.size(); j++){

		string Path_To_Histo = "MultiSamplePlot_MScutflow/MScutflow_";
		Path_To_Histo += Vector_SampleName[j];
		
		TH1F *Histo( (TH1F*) file->Get(Path_To_Histo.c_str()) );
		
		if(Histo) Efficiency_cutflows.push_back(Histo);
	}
	
	string outputName = "../data/Efficiency_cutflows_";
	outputName += channel;
	outputName += ".root";
	TFile *outputFile = new TFile(outputName.c_str(),"RECREATE"); 
	
	for( unsigned int i = 0; i < Efficiency_cutflows.size(); i++){
		double NbOfEvents = Efficiency_cutflows[i]->GetBinContent(2);
		
		//cout << "NbOfEvents " << NbOfEvents << endl; 
		Efficiency_cutflows[i]->Scale(1./NbOfEvents);
		Efficiency_cutflows[i]->Write();
		 
	}
	
	for( unsigned int k = 0; k < Efficiency_cutflows.size(); k++)
	{
		TH1F *histogram = (TH1F*) Efficiency_cutflows[k]; 
		//cout << "Histogram: " << k << " with name " << histogram->GetName() << " and title " << histogram->GetTitle() << endl;
		double NbOfbins = histogram->GetNbinsX(); 
		string HistoTitle = histogram->GetName();
		if(histogram->GetBinContent(2) != 0)
		{
		cout << "------------------------------------------------------------------" << endl;
		cout << "*** "<< HistoTitle  << " : " << "***"<< endl;
		cout << "" << endl; 
		for(unsigned int i = 1; i < NbOfbins; i ++)
		{
			int iBin = i; 
			char *BinName;
			
			if(channel.find("4L")!=string::npos)
			{
				BinName = (char*)Vector_cutflow_4L[i-1].c_str();	
			}
			if(channel.find("3L")!=string::npos)
			{
				BinName = (char*)Vector_cutflow_3L[i-1].c_str();
			}
			if(channel.find("5L")!=string::npos)
			{
				BinName = (char*)Vector_cutflow_5L[i-1].c_str();
			}
			if(channel.find("45")!=string::npos)
			{
				BinName = (char*)Vector_cutflow_45[i-1].c_str();
			}
			histogram->GetXaxis()->SetBinLabel(iBin, BinName);
			
			string binlabel = histogram->GetXaxis()->GetBinLabel(i);
			double icontent = histogram->GetBinContent(i); 
			if(!binlabel.empty())
			{
				cout << binlabel << " : " << icontent << " efficiency -- " << icontent*100 <<" % efficiency"<< endl;
			}
		}
		 
		}		
	}
	outputFile->Write();


}
