#include <iostream>
#include <stdio.h>
#include <cmath>

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TList.h"

#include "TSystem.h"
#include "TKey.h"
#include "TString.h"

 
using namespace std;

int main(int argc, char *argv[])
{
  cout<<"DISCLAIMER: this code to create a root file with error bands for MultiSamplePlots was used for a vector-like quark search at 8 TeV"<<endl;
	cout<<"The systematics sources and the SM processes to be treated are analysis dependent!"<<endl;
  cout<<"This macro will calculate the systematic error bands as the quadratic sum of all different sources."<<endl;
	cout<<"Note that it makes use of the fixed structure of MultiSamplePlots + some hardcoded strings to recognize systematics and physics-process names!"<<endl;

  unsigned int verbose = 0; //0 = minimal, 1 = normal, 2 = loud

  //string basedir = "/user/gvonsem/VectorLikeQuarkSearch/Git/CMSSW_5_3_6_patch1/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_VLQAnalyzer_withoutQCDestimation_7Jan13_signal500overlay_FittedWjetsTogetherMu_ZjetsFlavorSplitting/"; //OutputFiles_VLQAnalyzer_withoutQCDestimation_15Dec13_signal700overlay_allsignal_nomodelrescaling/";
	string basedir = "/user/gvonsem/VectorLikeQuarkSearch/Git/CMSSW_5_3_6_patch1/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_VLQAnalyzer_5Jul15_withQCDestimation_noModelScaling_COPYFORERRORBAND/";//"/user/gvonsem/VectorLikeQuarkSearch/Git/CMSSW_5_3_6_patch1/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_VLQAnalyzer_29Oct14_newTTbarXS_smallerbins/"   //  OutputFiles_VLQAnalyzer_withQCDestimation_15Mar14_FlatBinning_NoModelScaling_v4_VJetsSystFix/ -> probably the one used for thesis                        //"/user/gvonsem/VectorLikeQuarkSearch/Git/CMSSW_5_3_6_patch1/src/TopBrussels/TopTreeAnalysis/macros/OutputFiles_VLQAnalyzer_withQCDestimation_10Mar14_FlatBinning_NoModelScaling_v3/";

	string InputfilenameNominal = basedir+"VLQTreeAnalyzer_Nominal.root";
	string InputfilenameJESPlus = basedir+"VLQTreeAnalyzer_JESPlus.root";
	string InputfilenameJESMinus = basedir+"VLQTreeAnalyzer_JESMinus.root";
	string InputfilenameMETUnClusteredEnergyPlus = basedir+"VLQTreeAnalyzer_METUnClusteredEnergyPlus.root";
	string InputfilenameMETUnClusteredEnergyMinus = basedir+"VLQTreeAnalyzer_METUnClusteredEnergyMinus.root";
	string InputfilenameJERPlus = basedir+"VLQTreeAnalyzer_JERPlus.root";
	string InputfilenameJERMinus = basedir+"VLQTreeAnalyzer_JERMinus.root";
	string InputfilenamebTagPlus = basedir+"VLQTreeAnalyzer_bTagPlus.root";
	string InputfilenamebTagMinus = basedir+"VLQTreeAnalyzer_bTagMinus.root";
	string InputfilenamemisTagPlus = basedir+"VLQTreeAnalyzer_misTagPlus.root";
	string InputfilenamemisTagMinus = basedir+"VLQTreeAnalyzer_misTagMinus.root";	
	string InputfilenamePUPlus = basedir+"VLQTreeAnalyzer_PUPlus.root";
	string InputfilenamePUMinus = basedir+"VLQTreeAnalyzer_PUMinus.root";
	string InputfilenameMuonSFPlus = basedir+"VLQTreeAnalyzer_MuonSFPlus.root";
	string InputfilenameMuonSFMinus = basedir+"VLQTreeAnalyzer_MuonSFMinus.root";
	string InputfilenameElectronSFPlus = basedir+"VLQTreeAnalyzer_ElectronSFPlus.root";
	string InputfilenameElectronSFMinus = basedir+"VLQTreeAnalyzer_ElectronSFMinus.root";
	
	string Outputpath = "VLQSearch_ErrorBand/";
	string Outputfilename = "ErrorBandFile_15Jul15.root"; //"ErrorBandFile_29Oct14.root" //"ErrorBandFile_15Mar14.root"
//	mkdir(Outputpath.c_str(),0777); //compilation problem suddenly???? why???
	//bool mergeSignal = false;
	
	//systematics you want to consider in the errorband. The effects will be quadratically summed bin by bin (i.e. assumption = uncorrelated systematics, which is in most cases rather unrealistic).
  vector<string> syst;
  syst.push_back("JES");
//	syst.push_back("METUnClusteredEnergy");
  syst.push_back("JER");
	syst.push_back("bTag");
	syst.push_back("misTag");
//  syst.push_back("PU");
	syst.push_back("MuonSF");
	syst.push_back("ElectronSF");
	syst.push_back("top_XS");
  syst.push_back("w_XS");
  syst.push_back("z_XS");
  syst.push_back("st_XS");
	syst.push_back("ww_XS");
  syst.push_back("wz_XS");
	syst.push_back("zz_XS");
	syst.push_back("ttw_XS");
	syst.push_back("ttz_XS");
	syst.push_back("ssww_XS");
	syst.push_back("vvv_XS");
	syst.push_back("qcd_XS");
	//syst.push_back("lep_eff"); //now no overall factor anymore...
	syst.push_back("lumi");
	
	
	//numbers are to be revisited
	float rel_unc_top = 0.15;//0.05 //0.14 -> used in thesis //0.043 = sqrt(max(plus,minus variations)^2 + PDF^2)... 8% is the difference of NNLO with other calculations (according to the T2/3 search paper I think). 14% unc CMS measurement
	float rel_unc_w_LF = 0.20; //0.14 //5% ~SMP-12-011 https://cdsweb.cern.ch/record/1460098, via sqrt(stat + syst + lumi quadratic sum) (including lumi... ok or not?) //14% is difference between NLO (?) and fitted value (LF+HF)
	float rel_unc_z_LF = 0.20; //0.14 //5% ~SMP-12-011 but I should rerun everything with theory XS, the measurement is not entirely applicable
	float rel_unc_w_HF = 0.30; //5% ~SMP-12-011 https://cdsweb.cern.ch/record/1460098, via sqrt(stat + syst + lumi quadratic sum) (including lumi... ok or not?) //14% is difference between NLO (?) and fitted value (LF+HF)
	float rel_unc_z_HF = 0.30; //5% ~SMP-12-011 but I should rerun everything with theory XS, the measurement is not entirely applicable
	float rel_unc_ww = 0.10; //CMS measurement, http://arxiv.org/abs/1301.4698 (again taken lumi uncert too...)
	float rel_unc_wz = 0.04; //twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV doesn't give uncertainties for this... for the moment just taking the same as AN2013-006-v7 (T' in multilepton) i.e. 17%... Update, on the twiki now they give uncertainties... I get 4%...
	float rel_unc_zz = 0.04; //twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV doesn't give uncertainties for this... for the moment just taking the same as AN2013-006-v7 (T' in multilepton)... Is it really only 5%, not 50%?? Update, on the twiki now they give uncertainties... I get 4%...
  float rel_unc_ttw = 0.32; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
	float rel_unc_ttz = 0.12; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
	float rel_unc_ssww = 0.50; //for the moment just taking the same as AN2013-006-v7 (T' in multilepton)... conservative
	float rel_unc_st = 0.034; //~calculated from from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
	float rel_unc_vvv = 0.50; //for the moment just taking the same as AN2013-006-v7 (T' in multilepton)... conservative
	float rel_unc_qcd = 1.; //100%
	//now some overall scale systematics
	//float rel_unc_mu_eff 	= 0.03; //now no overall factor anymore...
	//float rel_unc_el_eff 	= 0.05; //now no overall factor anymore...
	float rel_unc_lumi = 0.026; //2012 rereco


  bool doIgnore = true; //true: ignore MSPlots with names including one of the strings in MSPlotsIgnore, false: strings in MSPlotsIgnore will not be used to ignore certain MSPlots.
  vector<string> MSPlotsIgnore;
	MSPlotsIgnore.push_back("BtagWeight");
	MSPlotsIgnore.push_back("prebox");
  bool doOnlyConsider = true; //true: only consider MSPlots with names including one of the strings in MSPlotsIgnore, false: strings in MSPlotsIgnore will not be used to 'only consider' certain MSPlots
  vector<string> MSPlotsOnlyConsider;
	MSPlotsOnlyConsider.push_back("box");
  //note: if plot name should both be 'ignored' and 'only considered', it will be ignored.




  TFile* outFile = new TFile((Outputpath+Outputfilename).c_str(),"RECREATE");

	TFile* inFileNominal = new TFile(InputfilenameNominal.c_str(),"READ");
	if(verbose >= 0) cout<<"Nominal file read = "<<inFileNominal->GetName()<<endl;
	TIter nextkey(inFileNominal->GetListOfKeys());
	TKey *key;
  if(verbose >= 0)
	{
	   cout<<"Loop over de directories (~the multisampleplots) in the file"<<endl;
	   cout<<"------------------------------------------------------------"<<endl;
	}
  while (key = (TKey*)nextkey()) 
	{			
	  TH1F* sumNominal = 0;
	  TH1F* errorMinus = 0; //lower values of error band
    TH1F* errorPlus = 0; //upper values of error band
				
	  TDirectoryFile* subdir = (TDirectoryFile*) key->ReadObj();
		string subdirname = subdir->GetName();
		if(verbose >= 0) cout<<"*** "<<subdirname<<" ***"<<endl;

		if(subdirname.find("MultiSamplePlot")==0)
		{
		 bool doPlot = true;
		 if(doOnlyConsider)
		 {
		   doPlot = false;
		   for(unsigned int j=0; j<MSPlotsOnlyConsider.size(); j++)
		   {
		     if(subdirname.find(MSPlotsOnlyConsider[j])!=string::npos)
				 {
				   doPlot = true;
					 break;
				 }
		   }
		 }
		 if(doIgnore)
		 {
		   for(unsigned int j=0; j<MSPlotsIgnore.size(); j++)
		   {
		     if(subdirname.find(MSPlotsIgnore[j])!=string::npos)
				 {
				   doPlot = false;
					 break;
				 }
		   }
		 }		
		
		 if(!doPlot)
		 {
		    if(verbose >= 0) cout<<" -> IGNORING "<<subdirname<<endl;
		 }
		 else
		 {
		  TIter nextplotkey(subdir->GetListOfKeys());
	    TKey *plotkey;
			if(verbose >= 0) cout<<"loop over histograms (in multisampleplot directory) and obtain total nominal histogram"<<endl;
			bool firstdataset = true;
      while (plotkey = (TKey*)nextplotkey()) 
	    {
				    TH1F* histo = 0;
						if(((string) (plotkey->ReadObj())->ClassName()) == "TH1F")
						{
						   histo = (TH1F*) plotkey->ReadObj();
							 if(!( ((string) histo->GetName()).find("Data")!=string::npos || ((string) histo->GetName()).find("data")!=string::npos || ((string) histo->GetName()).find("DATA")!=string::npos || ((string) histo->GetName()).find("_NP_")!=string::npos ) )
							 {
							    //WARNING: slightly dangerous to search for a string _NP_ in the plot 'name', since it is in principle possible to write this string in the title of a plot without realizing this is not allowed..
							    if(verbose >= 1) cout<<"      histo->GetName() = "<<histo->GetName()<<", integral = "<<histo->Integral(0,histo->GetNbinsX()+1)<<endl;
									if(firstdataset)
						      {
									  if(verbose >= 2) cout<<"         First dataset..."<<endl;
						        sumNominal = histo;
							      firstdataset = false;
						      }
				          else sumNominal->Add( histo );
							 }
						}
			}
			
			//maybe not needed, but to be sure; add the overflow bin to the last bin, (and underflow to first bin) because MSPlot does this too
		  sumNominal->SetBinContent(sumNominal->GetNbinsX(),sumNominal->GetBinContent(sumNominal->GetNbinsX()) + sumNominal->GetBinContent(sumNominal->GetNbinsX()+1));
		  sumNominal->SetBinContent(sumNominal->GetNbinsX()+1, 0);
			sumNominal->SetBinContent(1,sumNominal->GetBinContent(0) + sumNominal->GetBinContent(1));
		  sumNominal->SetBinContent(0, 0);
			if(verbose >= 1) cout<<"  => sumNominal integral = "<<sumNominal->Integral(1,sumNominal->GetNbinsX())<<endl;
					
			if(verbose >= 0) cout<<"loop over systematics"<<endl;
			//TFile* inFilePlus = 0; //testing
			//TFile* inFileMinus = 0; //testing
			for(unsigned int iSyst=0; iSyst<syst.size(); iSyst++)
      {			
			   if(verbose >= 0) cout<<" -> systematic "<<syst[iSyst]<<endl;				 
				 if(verbose >= 1) cout<<"    up variations"<<endl;
		     TFile* inFilePlus = 0;
				 if(syst[iSyst] == "JES") inFilePlus = new TFile(InputfilenameJESPlus.c_str(),"READ");
				 else if(syst[iSyst] == "METUnClusteredEnergy") inFilePlus = new TFile(InputfilenameMETUnClusteredEnergyPlus.c_str(),"READ");
		     else if(syst[iSyst] == "JER") inFilePlus = new TFile(InputfilenameJERPlus.c_str(),"READ");
         else if(syst[iSyst] == "bTag") inFilePlus = new TFile(InputfilenamebTagPlus.c_str(),"READ");
		     else if(syst[iSyst] == "misTag") inFilePlus = new TFile(InputfilenamemisTagPlus.c_str(),"READ");
		     else if(syst[iSyst] == "PU") inFilePlus = new TFile(InputfilenamePUPlus.c_str(),"READ");
				 else if(syst[iSyst] == "MuonSF") inFilePlus = new TFile(InputfilenameMuonSFPlus.c_str(),"READ");
				 else if(syst[iSyst] == "ElectronSF") inFilePlus = new TFile(InputfilenameElectronSFPlus.c_str(),"READ");
		     else inFilePlus = new TFile(InputfilenameNominal.c_str(),"READ");
		     if(verbose >= 1) cout<<"    Plusfile read = "<<inFilePlus->GetName()<<endl;
				 TDirectory* subdirPlus = (TDirectory*) inFilePlus->Get(subdirname.c_str());
         TH1F* sumPlus = 0;
				 
				 TIter nextplotkeyPlus(subdirPlus->GetListOfKeys());
	       TKey *plotkeyPlus;
				 bool firstdatasetPlus = true;
				 				 				 
				 if(syst[iSyst] == "JES" || syst[iSyst] == "METUnClusteredEnergy" || syst[iSyst] == "JER" || syst[iSyst] == "bTag" || syst[iSyst] == "misTag" || syst[iSyst] == "PU" || syst[iSyst] == "MuonSF" || syst[iSyst] == "ElectronSF") //WARNING: this is some necessary hardcoding, unfortunately...
         {
						if(verbose >= 1) cout<<"    loop over histograms (in multisampleplot directory) and obtain total up-varied histogram"<<endl;
            while (plotkeyPlus = (TKey*)nextplotkeyPlus()) 
	          {
				       TH1F* histoPlus = 0;
						   if(((string) (plotkeyPlus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoPlus = (TH1F*) plotkeyPlus->ReadObj();
							    if(!( ((string) histoPlus->GetName()).find("Data")!=string::npos || ((string) histoPlus->GetName()).find("data")!=string::npos || ((string) histoPlus->GetName()).find("DATA")!=string::npos || ((string) histoPlus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       if(verbose >= 1) cout<<"       histoPlus->GetName() = "<<histoPlus->GetName()<<", integral = "<<histoPlus->Integral(0,histoPlus->GetNbinsX()+1)<<endl;										 
								     if(firstdatasetPlus)
						         {
										   if(verbose >= 2) cout<<"          First dataset..."<<endl;
						           sumPlus = histoPlus;
							         firstdatasetPlus = false;
						         }
				             else sumPlus->Add( histoPlus );
							    }
						   }
				    }	 //end loop over histograms	
				 }
				 else
				 {
						if(verbose >= 1) cout<<"   loop over histograms (in multisampleplot directories) and obtain total up-varied histogram"<<endl;
            while (plotkeyPlus = (TKey*)nextplotkeyPlus()) 
	          {
				       TH1F* histoPlus = 0;
						   if(((string) (plotkeyPlus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoPlus = (TH1F*) plotkeyPlus->ReadObj();
							    if(!( ((string) histoPlus->GetName()).find("Data")!=string::npos || ((string) histoPlus->GetName()).find("data")!=string::npos || ((string) histoPlus->GetName()).find("DATA")!=string::npos || ((string) histoPlus->GetName()).find("_NP_")!=string::npos ) )
							    {
								     if(syst[iSyst]=="top_XS" && ((string) histoPlus->GetName()).find("TTbarJets")!=string::npos) //WARNING: this is some necessary hardcoding, unfortunately...
						         {
										    histoPlus->Scale(1+rel_unc_top);
						         }
										 else if(syst[iSyst] == "w_XS" && ((string) histoPlus->GetName()).find("WJets_LF")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_w_LF);
										 }
										 else if(syst[iSyst] == "w_XS" && ((string) histoPlus->GetName()).find("WJets_HF")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_w_HF);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoPlus->GetName()).find("ZJets_LF")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_z_LF);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoPlus->GetName()).find("ZJets_HF")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_z_HF);
										 }
										 else if(syst[iSyst] == "st_XS" && ((string) histoPlus->GetName()).find("ST_")!=string::npos)
										 {
										    histoPlus->Scale(1+rel_unc_st);
										 }
										 else if(syst[iSyst] == "ww_XS" && ((string) histoPlus->GetName()).find("WW")!=string::npos && ((string) histoPlus->GetName()).find("WWW")==string::npos && ((string) histoPlus->GetName()).find("WWZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ww);
										 }
										 else if(syst[iSyst] == "wz_XS" && ((string) histoPlus->GetName()).find("WZ")!=string::npos && ((string) histoPlus->GetName()).find("WWZ")==string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_wz);
										 }
										 else if(syst[iSyst] == "zz_XS" && ((string) histoPlus->GetName()).find("ZZ")!=string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos && ((string) histoPlus->GetName()).find("ZZZ")==string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_zz);
										 }
										 else if(syst[iSyst] == "ttw_XS" && ((string) histoPlus->GetName()).find("ttW")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ttw);
										 }
										 else if(syst[iSyst] == "ttz_XS" && ((string) histoPlus->GetName()).find("ttZ")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_ttz);
										 }
										 else if(syst[iSyst] == "ssww_XS" && ( ((string) histoPlus->GetName()).find("WpWp")!=string::npos || ((string) histoPlus->GetName()).find("WmWm")!=string::npos ) )
			               {
										    histoPlus->Scale(1+rel_unc_ssww);
										 }
										 else if(syst[iSyst] == "vvv_XS" && ( ((string) histoPlus->GetName()).find("WWW")!=string::npos || ((string) histoPlus->GetName()).find("WWZ")!=string::npos || ((string) histoPlus->GetName()).find("WZZ")!=string::npos || ((string) histoPlus->GetName()).find("ZZZ")!=string::npos ) )
			               {
										    histoPlus->Scale(1+rel_unc_vvv);
										 }
										 else if(syst[iSyst] == "qcd_XS" && ((string) histoPlus->GetName()).find("InvIso")!=string::npos)
			               {
										    histoPlus->Scale(1+rel_unc_qcd);
										 }
										 else if(syst[iSyst] == "lumi")
			               {
										    histoPlus->Scale(1+rel_unc_lumi);
										 }
										 
										 if(verbose >= 1) cout<<"       histoPlus->GetName() = "<<histoPlus->GetName()<<", integral = "<<histoPlus->Integral(0,histoPlus->GetNbinsX()+1)<<endl;
										 
										 if(firstdatasetPlus)
						         {
										   if(verbose >= 2) cout<<"          First dataset..."<<endl;
						           sumPlus = histoPlus;
							         firstdatasetPlus = false;
						         }
				             else sumPlus->Add( histoPlus );										 
							    }
						   }
				    } //end loop over histograms				 
				 }
				 //maybe not needed, but to be sure; add the overflow bin to the last bin, (and underflow to first bin) because MSPlot does this too
		     sumPlus->SetBinContent(sumPlus->GetNbinsX(),sumPlus->GetBinContent(sumPlus->GetNbinsX()) + sumPlus->GetBinContent(sumPlus->GetNbinsX()+1));
		     sumPlus->SetBinContent(sumPlus->GetNbinsX()+1, 0);
			   sumPlus->SetBinContent(1,sumPlus->GetBinContent(0) + sumPlus->GetBinContent(1));
		     sumPlus->SetBinContent(0, 0);
				 if(verbose >= 1) cout<<"   => sumPlus integral = "<<sumPlus->Integral(1,sumPlus->GetNbinsX())<<endl;
				 sumPlus->SetDirectory(0);
         inFilePlus->Close();
         delete inFilePlus;
				 
		     if(verbose >= 1) cout<<"    minus variations"<<endl;
		     TFile* inFileMinus = 0;
				 if(syst[iSyst] == "JES") inFileMinus = new TFile(InputfilenameJESMinus.c_str(),"READ");
				 else if(syst[iSyst] == "METUnClusteredEnergy") inFileMinus = new TFile(InputfilenameMETUnClusteredEnergyMinus.c_str(),"READ");
		     else if(syst[iSyst] == "JER") inFileMinus = new TFile(InputfilenameJERMinus.c_str(),"READ");
         else if(syst[iSyst] == "bTag") inFileMinus = new TFile(InputfilenamebTagMinus.c_str(),"READ");
		     else if(syst[iSyst] == "misTag") inFileMinus = new TFile(InputfilenamemisTagMinus.c_str(),"READ");
		     else if(syst[iSyst] == "PU") inFileMinus = new TFile(InputfilenamePUMinus.c_str(),"READ");
				 else if(syst[iSyst] == "MuonSF") inFileMinus = new TFile(InputfilenameMuonSFMinus.c_str(),"READ");
				 else if(syst[iSyst] == "ElectronSF") inFileMinus = new TFile(InputfilenameElectronSFMinus.c_str(),"READ");
		     else inFileMinus = new TFile(InputfilenameNominal.c_str(),"READ");
				 if(verbose >= 1) cout<<"   Minusfile read = "<<inFileMinus->GetName()<<endl;
		     TDirectory* subdirMinus = (TDirectory*) inFileMinus->Get(subdirname.c_str());
         TH1F* sumMinus = 0;
				 
				 TIter nextplotkeyMinus(subdirMinus->GetListOfKeys());
	       TKey *plotkeyMinus;
				 bool firstdatasetMinus = true;
				 
				 if(syst[iSyst] == "JES" || syst[iSyst] == "METUnClusteredEnergy" || syst[iSyst] == "JER" || syst[iSyst] == "bTag" || syst[iSyst] == "misTag" || syst[iSyst] == "PU" || syst[iSyst] == "MuonSF" || syst[iSyst] == "ElectronSF")
         {
						if(verbose >= 1) cout<<"   loop over histograms (in multisampleplot directories) and obtain total down-varied histogram"<<endl;
            while (plotkeyMinus = (TKey*)nextplotkeyMinus()) 
	          {
				       TH1F* histoMinus = 0;
						   if(((string) (plotkeyMinus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoMinus = (TH1F*) plotkeyMinus->ReadObj();
							    if(!( ((string) histoMinus->GetName()).find("Data")!=string::npos || ((string) histoMinus->GetName()).find("data")!=string::npos || ((string) histoMinus->GetName()).find("DATA")!=string::npos || ((string) histoMinus->GetName()).find("_NP_")!=string::npos ) )
							    {
							       if(verbose >= 1) cout<<"       histoMinus->GetName() = "<<histoMinus->GetName()<<", integral = "<<histoMinus->Integral(0,histoMinus->GetNbinsX()+1)<<endl;
								     if(firstdatasetMinus)
						         {
										   if(verbose >= 2) cout<<"         First dataset..."<<endl;
						           sumMinus = histoMinus;
							         firstdatasetMinus = false;
						         }
				             else sumMinus->Add( histoMinus );
							    }
						   }
				    }	
				}
				else
				{
						if(verbose >= 1) cout<<"   loop over histograms (in multisampleplot directories) and obtain total down-varied histogram"<<endl;
            while (plotkeyMinus = (TKey*)nextplotkeyMinus()) 
	          {
				       TH1F* histoMinus = 0;
						   if(((string) (plotkeyMinus->ReadObj())->ClassName()) == "TH1F")
						   {
						      histoMinus = (TH1F*) plotkeyMinus->ReadObj();
							    if(!( ((string) histoMinus->GetName()).find("Data")!=string::npos || ((string) histoMinus->GetName()).find("data")!=string::npos || ((string) histoMinus->GetName()).find("DATA")!=string::npos || ((string) histoMinus->GetName()).find("_NP_")!=string::npos ) )
							    {
								     if(syst[iSyst]=="top_XS" && ((string) histoMinus->GetName()).find("TTbarJets")!=string::npos) //WARNING: this is some necessary hardcoding, unfortunately...
						         {
										    histoMinus->Scale(1-rel_unc_top);
						         }
										 else if(syst[iSyst] == "w_XS" && ((string) histoMinus->GetName()).find("WJets_LF")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_w_LF);
										 }
										 else if(syst[iSyst] == "w_XS" && ((string) histoMinus->GetName()).find("WJets_HF")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_w_HF);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoMinus->GetName()).find("ZJets_LF")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_z_LF);
										 }
										 else if(syst[iSyst] == "z_XS" && ((string) histoMinus->GetName()).find("ZJets_HF")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_z_HF);
										 }
										 else if(syst[iSyst] == "st_XS" && ((string) histoMinus->GetName()).find("ST_")!=string::npos)
										 {
										    histoMinus->Scale(1-rel_unc_st);
										 }
										 else if(syst[iSyst] == "ww_XS" && ((string) histoMinus->GetName()).find("WW")!=string::npos && ((string) histoMinus->GetName()).find("WWW")==string::npos && ((string) histoMinus->GetName()).find("WWZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ww);
										 }
										 else if(syst[iSyst] == "wz_XS" && ((string) histoMinus->GetName()).find("WZ")!=string::npos && ((string) histoMinus->GetName()).find("WWZ")==string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_wz);
										 }
										 else if(syst[iSyst] == "zz_XS" && ((string) histoMinus->GetName()).find("ZZ")!=string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos && ((string) histoMinus->GetName()).find("ZZZ")==string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_zz);
										 }
										 else if(syst[iSyst] == "ttw_XS" && ((string) histoMinus->GetName()).find("ttW")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ttw);
										 }
										 else if(syst[iSyst] == "ttz_XS" && ((string) histoMinus->GetName()).find("ttZ")!=string::npos)
			               {
										    histoMinus->Scale(1-rel_unc_ttz);
										 }
										 else if(syst[iSyst] == "ssww_XS" && ( ((string) histoMinus->GetName()).find("WpWp")!=string::npos || ((string) histoMinus->GetName()).find("WmWm")!=string::npos ) )
			               {
										    histoMinus->Scale(1-rel_unc_ssww);
										 }
										 else if(syst[iSyst] == "vvv_XS" && ( ((string) histoMinus->GetName()).find("WWW")!=string::npos || ((string) histoMinus->GetName()).find("WWZ")!=string::npos || ((string) histoMinus->GetName()).find("WZZ")!=string::npos || ((string) histoMinus->GetName()).find("ZZZ")!=string::npos ) )
			               {
										    histoMinus->Scale(1-rel_unc_vvv);
										 }
										 else if(syst[iSyst] == "qcd_XS" && ((string) histoMinus->GetName()).find("InvIso")!=string::npos)
			               {
										    histoMinus->Scale(1+rel_unc_qcd);
										 }
										 else if(syst[iSyst] == "lumi")
			               {
										    histoMinus->Scale(1-rel_unc_lumi);
										 }
										 
										 if(verbose >= 1) cout<<"       histoMinus->GetName() = "<<histoMinus->GetName()<<", integral = "<<histoMinus->Integral(0,histoMinus->GetNbinsX()+1)<<endl;
										 										 
										 if(firstdatasetMinus)
						         {
										   if(verbose >= 2) cout<<"          First dataset..."<<endl;
						           sumMinus = histoMinus;
							         firstdatasetMinus = false;
						         }
				             else sumMinus->Add( histoMinus );
										 
							    }
						   }
				    }				 //end loop over histograms		
				}	 
				//maybe not needed, but to be sure; add the overflow bin to the last bin, (and underflow to first bin) because MSPlot does this too
		    sumMinus->SetBinContent(sumMinus->GetNbinsX(),sumMinus->GetBinContent(sumMinus->GetNbinsX()) + sumMinus->GetBinContent(sumMinus->GetNbinsX()+1));
		    sumMinus->SetBinContent(sumMinus->GetNbinsX()+1, 0);
			  sumMinus->SetBinContent(1,sumMinus->GetBinContent(0) + sumMinus->GetBinContent(1));
		    sumMinus->SetBinContent(0, 0);
				if(verbose >= 1) cout<<"   => sumMinus integral = "<<sumMinus->Integral(1,sumMinus->GetNbinsX())<<endl;
				sumMinus->SetDirectory(0);
        inFileMinus->Close();
        delete inFileMinus;
		
		
		    if(verbose >= 1) cout << "   Current integrals for " << subdirname << ", syst " << syst[iSyst] << " -> Nominal: " << sumNominal->Integral(1,sumNominal->GetNbinsX()) << "  Plus: " << sumPlus->Integral(1,sumPlus->GetNbinsX()) << "  Minus: " << sumMinus->Integral(1,sumMinus->GetNbinsX()) << endl;
		 		
		    //maybe not needed, but to be sure; add the overflow bin to the last bin, because MSPlot does this as well
		    sumMinus->SetBinContent(sumMinus->GetNbinsX(),sumMinus->GetBinContent(sumMinus->GetNbinsX()) + sumMinus->GetBinContent(sumMinus->GetNbinsX()+1));
		    sumMinus->SetBinContent(sumMinus->GetNbinsX()+1, 0);
		    sumPlus->SetBinContent(sumPlus->GetNbinsX(),sumPlus->GetBinContent(sumPlus->GetNbinsX()) + sumPlus->GetBinContent(sumPlus->GetNbinsX()+1));
		    sumPlus->SetBinContent(sumPlus->GetNbinsX()+1, 0);	 
		  
			  if(verbose >= 1) cout<<"   looping over bins to set bin contents of error histograms"<<endl;
		    if(iSyst == 0) //only the first time in the loop over systematics...
        {
           errorMinus = (TH1F*) sumMinus->Clone();
           errorPlus = (TH1F*) sumPlus->Clone();
					 for(unsigned int iBin=0; iBin<=errorMinus->GetNbinsX()+1; iBin++)
           { 
             if(sumMinus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					   {
               errorMinus->SetBinContent(iBin, sumMinus->GetBinContent(iBin));
               if(verbose >= 2) cout<<"          (Normal hierarchy) errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					   }
					   else
					   {
               errorPlus->SetBinContent(iBin, sumMinus->GetBinContent(iBin));
               if(verbose >= 2) cout<<"          (Inverted hierarchy) errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
				     }
			
             if(sumPlus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					   {
               errorMinus->SetBinContent(iBin, sumPlus->GetBinContent(iBin));
               if(verbose >= 2) cout<<"          (Inverted hierarchy) errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					   }
					   else
					   {
               errorPlus->SetBinContent(iBin, sumPlus->GetBinContent(iBin));
               if(verbose >= 2) cout<<"          (Normal hierarchy) errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
					   }
				   }
        }
				else
				{
         for(unsigned int iBin=0; iBin<=errorMinus->GetNbinsX()+1; iBin++)
         { 
          if(sumMinus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					{
            errorMinus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) - sqrt( pow(errorMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) ); //basically 'adding' the effect (sumMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin)) in quadrature to the previous effects (errorMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin))
            if(verbose >= 2) cout<<"          (Normal hierarchy) errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					}
					else
					{
            errorPlus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) + sqrt( pow(errorPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            if(verbose >= 2) cout<<"          (Inverted hierarchy) errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
				  }
			
          if(sumPlus->GetBinContent(iBin) < sumNominal->GetBinContent(iBin))
					{
            errorMinus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) - sqrt( pow(errorMinus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            if(verbose >= 2) cout<<"          (Inverted hierarchy) errorMinus->GetBinContent("<<iBin<<") = "<<errorMinus->GetBinContent(iBin)<<endl;
					}
					else
					{
            errorPlus->SetBinContent(iBin, sumNominal->GetBinContent(iBin) + sqrt( pow(errorPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) + pow(sumPlus->GetBinContent(iBin) - sumNominal->GetBinContent(iBin),2) ) );
            if(verbose >= 2) cout<<"          (Normal hierarchy) errorPlus->GetBinContent("<<iBin<<") = "<<errorPlus->GetBinContent(iBin)<<endl;
					}
				 }
				}				
		    if(verbose >= 1) cout <<"       >>>> Updated summary for " << subdirname << " after syst "<<syst[iSyst]<<" -> errorPlus: " << errorPlus->Integral(1,errorPlus->GetNbinsX()) << "  errorMinus: " << errorMinus->Integral(1,errorMinus->GetNbinsX()) << " <<<<" << endl;
		 
		  } 

			if(verbose >= 1)
			{
			   cout << "Loop over systematics done." << endl;
			   cout << "INTEGRALS for plot " << subdirname << " ->   errorPlus: " << errorPlus->Integral(1,errorPlus->GetNbinsX()) << ", nominal: " << sumNominal->Integral(1,sumNominal->GetNbinsX()) << ",  errorMinus: " << errorMinus->Integral(1,errorMinus->GetNbinsX()) << endl << endl;
      }
			
      outFile->cd();
		  if(outFile->Get(subdirname.c_str())==0)
		    outFile->mkdir(subdirname.c_str());
		  outFile->cd(subdirname.c_str());
		  errorMinus->SetNameTitle("Minus","Minus");
			errorMinus->SetLineColor(1);
      errorMinus->Write();
      errorPlus->SetNameTitle("Plus","Plus");
			errorPlus->SetLineColor(1);
      errorPlus->Write();
	    sumNominal->SetNameTitle("Nominal","Nominal");
			sumNominal->SetLineColor(1);
	    sumNominal->Write();

			
     }
	  }					
  }
	outFile->Close();
  delete outFile;
	inFileNominal->Close();
	delete inFileNominal;
	
	cout<<"Error band histograms written in output file: "<<Outputpath+Outputfilename<<endl;

}
