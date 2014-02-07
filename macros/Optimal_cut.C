#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include <TCanvas.h>
#include <vector>
#include "TGraph.h"
#include "TVectorT.h"

using namespace std;

string infile = "../data/FCNC_selection_4L5L.root";
string eclusiveinfile = "";
bool debug = false;
int jetcut = 5 + 1; //+1 for getting the right binnumber
int bcut = 3 +1; //+1 for the binnumber

void Optimal_cut(){

	TFile *file = new TFile(infile.c_str(),"read");
	//TFile *exclusiveFile = new TFile(eclusiveinfile.c_str(),"read");
	TFile *outputfile = new TFile("../data/OptCut/OptimalCuts_4L5L.root","RECREATE");
	
	vector <string> Variables;
		
	Variables.push_back("NbOfSelectedJets");
	Variables.push_back("NbOfSelectedBJets_CSVM");
	Variables.push_back("Pt_leading_jet");
	Variables.push_back("Pt_leading_Bjet_CSVM");
	Variables.push_back("Pt_leading_lepton");
	Variables.push_back("Pt_2nd_leading_jet");
	Variables.push_back("Pt_2nd_leading_Bjet_CSVM");
	Variables.push_back("Pt_2nd_leading_lepton");
	Variables.push_back("MET");
	Variables.push_back("Mllll");
	Variables.push_back("MllW");
	Variables.push_back("JetHt");
	Variables.push_back("LepHt");
	Variables.push_back("JetPx");
	Variables.push_back("LepPx");
	Variables.push_back("JetPy");
	Variables.push_back("LepPy");
	Variables.push_back("DeltaPhi");
	Variables.push_back("NofElectrons_af");
	Variables.push_back("NofMuons_af");
	
	


	
	vector <string> signalname;
	vector <string> backgroundnames;
		
	backgroundnames.push_back("GluGluHiggs4lep");	
/*	backgroundnames.push_back("W_1Jets");
	backgroundnames.push_back("W_2Jets");
	backgroundnames.push_back("W_3Jets");
	backgroundnames.push_back("W_4Jets");
	backgroundnames.push_back("WW_To2L2Nu");
	backgroundnames.push_back("WZ_To2L2Q");
*/	backgroundnames.push_back("WZ_To3LNu");
/*	backgroundnames.push_back("ZZ_To2L2Nu");
	backgroundnames.push_back("ZZ_To2L2Q");
*/	backgroundnames.push_back("ZZ_To4L");
/*	backgroundnames.push_back("ST_TToDilepton_tW-ch");
	backgroundnames.push_back("ST_TToTlepWhad_tW-ch");
	backgroundnames.push_back("ST_TToThadWlep_tW-ch");
	backgroundnames.push_back("ST_TBarToDilepton_tW-ch");
	backgroundnames.push_back("ST_TBarToTlepWhad_tW-ch");
	backgroundnames.push_back("ST_TBarToThadWlep_tW-ch");
	backgroundnames.push_back("TT_SemiLeptMGDecays");
	backgroundnames.push_back("TT_FullLeptMGDecays");
	backgroundnames.push_back("TT_HadronicMGDecays");
	backgroundnames.push_back("Z_M-10To50");
	backgroundnames.push_back("Z_M-50");
	backgroundnames.push_back("Z_1Jets");
	backgroundnames.push_back("Z_2Jets");
	backgroundnames.push_back("Z_3Jets");
	backgroundnames.push_back("Z_4Jets");
*/	backgroundnames.push_back("TTZ");
	backgroundnames.push_back("TTW");
/*	backgroundnames.push_back("ttbar");
	backgroundnames.push_back("ttbar_fullLept");
      	backgroundnames.push_back("ttbar_semiLept");
       	backgroundnames.push_back("wjets");
       	backgroundnames.push_back("ttt");
       	backgroundnames.push_back("ttw");
       	backgroundnames.push_back("WW");
       	backgroundnames.push_back("WZ");
       	backgroundnames.push_back("ZZ");
       	backgroundnames.push_back("ttz");
        backgroundnames.push_back("Zjets");
       	backgroundnames.push_back("ST_T_tW-ch");
       	backgroundnames.push_back("ST_TBar_tW-ch");
	backgroundnames.push_back("ST_T_s-ch");
	backgroundnames.push_back("ST_Tbar_s-ch");
	backgroundnames.push_back("ST_T_t-ch");
	backgroundnames.push_back("ST_Tbar_t-ch");
*/
//	signalname.push_back("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctR");
	signalname.push_back("TTJetsTocHbW_HToWW_WToLNuL_WToJets_HctL");
	signalname.push_back("TTJetsTocHbW_HToWW_WToLNuL_HctL");
//	signalname.push_back("TTJetsTocHbW_HToWW_WToLNuL_HctR");
	signalname.push_back("TTJetsTocHbW_HToBB_HctL");
//	signalname.push_back("TTJetsTocHbW_HToBB_HctR");
	signalname.push_back("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctL");
//	signalname.push_back("TTJetsTocHbW_HToZZ_ZToBB_ZToLL_HctR");
	signalname.push_back("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctL");
//	signalname.push_back("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL_HctR");
	signalname.push_back("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctL");
//	signalname.push_back("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL_HctR");
	signalname.push_back("TTJetsTocHbW_HToZZ_ZToLL_HctL");
//	signalname.push_back("TTJetsTocHbW_HToZZ_ZToLL_HctR");
	signalname.push_back("TTJetsTocZbW");
	
	for(unsigned int iVar = 0; iVar<Variables.size(); iVar++){
	
		TDirectory* th1dir = outputfile->mkdir(Variables[iVar].c_str());
		th1dir->cd();


		vector <TH1F*> Histo_samples;
		vector <string> Name_samples;
		
		string Path_To_Histo = "MultiSamplePlot_";
		Path_To_Histo += Variables[iVar];
		Path_To_Histo += "/";
		Path_To_Histo += Variables[iVar];
		Path_To_Histo += "_";
	
		int first_index_S = -1; // Define the first indices for S and B, which indicate which histograms are filled and which not from all samples
		int first_index_B = -1;
		
		TH1F* HistoSignal=0;
		
		for(unsigned int iSignal = 0; iSignal < signalname.size(); iSignal++){
			if(HistoSignal) continue; //We only want to go further if the HistoSignal isn't filled yet
			
			string histoName = Path_To_Histo;
			cout << "Path_To_Histo: " << Path_To_Histo << endl; 
			histoName += signalname[iSignal];
			cout << "histoName; " << histoName << endl; 
			TH1F *histo( (TH1F*) file->Get(histoName.c_str()) );
			//TH1F *histo_exclusive( (TH1F*) exclusiveFile->Get(histoName.c_str()) );
			
			TH1F* Signal_histogram = 0; 
			//TH1F* Histo = 0; 
			

			if(histo){
				histo->Clone("Signal_histogram");
				HistoSignal = (TH1F*) histo->Clone("HistoSignal");//Signal_histogram;
				cout << "Made HistoSignal? " << HistoSignal << endl; 
				first_index_S = iSignal;
				Histo_samples.push_back(histo);
				Name_samples.push_back(signalname[iSignal]);
				
				
				cout << "debug" << endl; 
				cout << histo->Integral(6,15) << endl;			
				cout << Histo_samples.back()->Integral(6,15) << endl;
				//cout << HistoSignal->Integral(6,15) << endl; 
				cout << " " << endl;
				
				continue;
			}		
			/*if(histo_exclusive){
				histo_exclusive->Clone("Histo");
				HistoSignal = Histo;
				first_index_S = iSignal;
				Histo_samples.push_back(histo_exclusive);
				Name_samples.push_back(signalname[iSignal]);			
				
				if(debug){
				cout << histo_exclusive->Integral(6,15) << endl;			
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				}
			}*/
			
		}		
		TH1F *HistoBackground = 0;
		for(unsigned int iBackgr = 0; iBackgr < backgroundnames.size(); iBackgr++){
			if(HistoBackground) continue; //We only want to go further if the HistoBackground isn't filled yet
			string histoName = Path_To_Histo;
			cout << "Path_To_Histo" << Path_To_Histo << endl; 
			histoName += backgroundnames[iBackgr];
			cout << "histoName" << histoName << endl; 
			TH1F *histo( (TH1F*) file->Get(histoName.c_str()) );
			//TH1F *histo_exclusive( (TH1F*) exclusiveFile->Get(histoName.c_str()) );

			TH1F* Background_histogram =0;
			//TH1F* Histogr=0;
			
			if(histo){
				histo->Clone("Background_histogram");
				HistoBackground = (TH1F*) histo->Clone("HistoBackground");
				first_index_B = iBackgr;
				Histo_samples.push_back(histo);
				Name_samples.push_back(backgroundnames[iBackgr]);
				
				cout << "debug" << endl; 
				cout << histo->Integral(6,15) << endl;
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				
				continue;
			}		
		/*	if(histo_exclusive){
				histo_exclusive->Clone("Histogr");
				HistoBackground = Histogr;
				first_index_B = iBackgr;
				Histo_samples.push_back(histo_exclusive);
				Name_samples.push_back(backgroundnames[iBackgr]);
				
				if(debug){
				cout << histo_exclusive->Integral(6,15) << endl;
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				}
			}
		*/
		}
	
                 cout << "before loop" << endl; 
		 
		cout << "first_index_S " << first_index_S << endl; 
		cout << "signalname.size() " << signalname.size() << endl;
		//Add the histograms into 1 for S and B seperately
		if(signalname.size()>1)
		{
		for(unsigned int iSignal = first_index_S+1; iSignal < signalname.size(); iSignal++){
			
			cout << "in loop" << endl; 
			string histoName = Path_To_Histo;
			cout << "Path_To_Histo " << Path_To_Histo << endl; 
			histoName += signalname[iSignal];
			cout << "histoName " << histoName << endl;
			TH1F *histo_0( (TH1F*) file->Get(histoName.c_str()) );
			//TH1F *histo_exclusive( (TH1F*) exclusiveFile->Get(histoName.c_str()) );
			
			cout << "before if statement" << endl; 
			if(histo_0){
				cout << "add" << endl; 
				HistoSignal->Add(histo_0);
				cout << "added" << endl; 
				Histo_samples.push_back(histo_0);
				Name_samples.push_back(signalname[iSignal]);

				cout << "debug" << endl; 
				cout << histo_0->Integral(6,15) << endl;		
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				
				continue;
			}
			/*if(histo_exclusive){
				HistoSignal->Add( histo_exclusive );
				Histo_samples.push_back(histo_exclusive);
				Name_samples.push_back(signalname[iSignal]);			
				
				if(debug){
				cout << histo_exclusive->Integral(6,15) << endl;
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				}
			}*/
		}
		}
		
		
		if(backgroundnames.size()>1){
		for(unsigned int iBackgr = first_index_B+1; iBackgr < backgroundnames.size(); iBackgr++){
			

			string histoName = Path_To_Histo;
			histoName += backgroundnames[iBackgr];
			TH1F *histo( (TH1F*) file->Get(histoName.c_str()) );
			//OTH1F *histo_exclusive((TH1F*) exclusiveFile->Get(histoName.c_str()));

			if(histo){
				HistoBackground->Add( histo );
				Histo_samples.push_back(histo);
				Name_samples.push_back(backgroundnames[iBackgr]);

				if(debug){
				cout << histo->Integral(6,15) << endl;
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				}
				continue;			
			}
			/*if(histo_exclusive){
				HistoBackground->Add( histo_exclusive );
				Histo_samples.push_back(histo_exclusive);
				Name_samples.push_back(backgroundnames[iBackgr]);			

				if(debug){
				cout << histo_exclusive->Integral(6,15) << endl;
				cout << Histo_samples.back()->Integral(6,15) << endl;
				cout << " " << endl;
				}
			}*/
		}
		}
	
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		//// Now we get into the real part where we find the optimal cut, where the significance is the highest
		////////////////////////////////////////////////////////////////////////////////////////////////////
		cout << "end" << endl; 
		cout << HistoSignal->GetNbinsX() << endl; 
		int end = HistoSignal->GetNbinsX();
		cout << "end" << end << endl; 
	
		//Efficiencies calculating as #events_passing_cut/#Total_events
		double Total_signal = 0;
		double Total_background = 0;
		Total_signal = HistoSignal->Integral();
		Total_background = HistoBackground->Integral();
	
		double * Signal_Integral_PerBin = new double [end];
		double * Background_Integral_PerBin = new double [end];
		double * Eff_Signal = new double [end];
		double * RejectionEff_Background = new double [end];
		
		for(unsigned int i = 0; i< end; i++){
			double s = 0;
			double b = 0;
		
			s = (HistoSignal->Integral(i , end));
			b = (HistoBackground->Integral(i,end));
		
			Signal_Integral_PerBin[i] = s;
			Background_Integral_PerBin[i] = b;
			Eff_Signal[i] = s/Total_signal;
			RejectionEff_Background[i] = (1- b/Total_background);
		}
	

	
		//Determine the optimal cut-value for a cut-and-count experiment
		string optcutName = "Opt_cut";
		TH1F *Opt_cut = new TH1F(optcutName.c_str(),optcutName.c_str(), end, HistoSignal->GetXaxis()->GetXmin(), HistoSignal->GetXaxis()->GetXmax());
		for(unsigned int i = 0; i<end; i++){
			double signal_significance = Signal_Integral_PerBin[i]/sqrt(Background_Integral_PerBin[i]);
			if(Background_Integral_PerBin[i] == 0) signal_significance = 1;
	
			Opt_cut->SetBinContent(i, signal_significance);
		}
		Opt_cut->GetXaxis()->SetTitle(Variables[iVar].c_str());
		Opt_cut->GetYaxis()->SetTitle("Signif.");
		
		string signaleeffName = "Signal_eff";
		TH1F *Signal_eff = new TH1F(signaleeffName.c_str(),signaleeffName.c_str(), end, HistoSignal->GetXaxis()->GetXmin(), HistoSignal->GetXaxis()->GetXmax());
		for(unsigned int i = 0; i<end; i++){
			Signal_eff->SetBinContent(i, Eff_Signal[i]);
		}
		Signal_eff->GetXaxis()->SetTitle(Variables[iVar].c_str());
		Signal_eff->GetYaxis()->SetTitle("Eff.");
		
		string brejName = "B_rej";
		TH1F *B_rej = new TH1F(brejName.c_str(),brejName.c_str(), end, HistoSignal->GetXaxis()->GetXmin(), HistoSignal->GetXaxis()->GetXmax());
		for(unsigned int i = 0; i<end; i++){
			B_rej->SetBinContent(i, RejectionEff_Background[i]);
		}
		B_rej->GetXaxis()->SetTitle(Variables[iVar].c_str());
		B_rej->GetYaxis()->SetTitle("Rejection eff.");

		
		cout << "******************************************************************" << endl;
		cout << "Optimal cut efficiencies for each sample " << Variables[iVar] << endl;
		if(Variables[iVar] == "NbOfSelectedJets")cout << ". Cutvalue: " << jetcut - 1 << endl;
		if(Variables[iVar] == "NbOfSelectedBJets_CSVM")cout << ". Cutvalue: " << bcut -1 << endl;
		cout << "******************************************************************" << endl;

		for(unsigned int i =0; i<Histo_samples.size(); i++){
			double Nevents = -999;
			if(Variables[iVar] == "NbOfSelectedJets") Nevents = Histo_samples[i]->Integral(jetcut,end);
			if(Variables[iVar] == "NbOfSelectedBJets_CSVM")Nevents = Histo_samples[i]->Integral(bcut,end);
			
			Histo_samples[i]->Write();
			
			cout << Variables[iVar] << "... Efficiency for optimal cut " << Name_samples[i] << ": " << Nevents << endl;
		}

		
		Histo_samples.clear();
		Name_samples.clear();
	}
	
		
	outputfile->Write();
	

}





