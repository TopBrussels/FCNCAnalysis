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
#include <cstdlib>

using namespace std;

//int main(unsigned int argc, char *argv[])
int CreateErrorBands( unsigned int verbose = 0, string basedir = "ErrorBandInput/", bool doOnlyConsider = true) //0 = silent, 1 = normal, 2 = loud
{
  cout<<"This macro will calculate the systematic error bands."<<endl;
  cout<<"Note that it makes use of the fixed structure of MultiSamplePlots + some hardcoded strings to recognize systematics and physics-process names!"<<endl;
  
  // unsigned int verbose = 0; //0 = silent, 1 = normal, 2 = loud
  
  //string basedir = "OutputPlots/170823_1854/";
  string InputfilenameNominal = basedir+"NtuplePlots.root"; // nominal
  
  // systematics
  string InputfilenameJESPlus = basedir+"NtuplePlots_JESPlus.root";
  string InputfilenameJESMinus = basedir+"NtuplePlots_JESMinus.root";
  
  string InputfilenameJERPlus = basedir+"NtuplePlots_JERPlus.root";
  string InputfilenameJERMinus = basedir+"NtuplePlots_JERMinus.root";
  
  string InputfilenamebTagPlus = basedir+"NtuplePlots_bTagPlus.root";
  string InputfilenamebTagMinus = basedir+"NtuplePlots_bTagMinus.root";
  
  string InputfilenamePUPlus = basedir+"NtuplePlots_PUPlus.root";
  string InputfilenamePUMinus = basedir+"NtuplePlots_PUMinus.root";
  
  string InputfilenameMuonSFPlus = basedir+"NtuplePlots_MuonSFPlus.root";
  string InputfilenameMuonSFMinus = basedir+"NtuplePlots_MuonSFMinus.root";
  
  string InputfilenameElectronSFPlus = basedir+"NtuplePlots_ElectronSFPlus.root";
  string InputfilenameElectronSFMinus = basedir+"NtuplePlots_ElectronSFMinus.root";
  
  //btag
  string Inputfilenamecferr1Minus = basedir + "NtuplePlots_cferr1Minus.root";
  string Inputfilenamecferr1Plus = basedir + "NtuplePlots_cferr1Plus.root";
  string Inputfilenamecferr2Minus = basedir + "NtuplePlots_cferr2Minus.root";
  string Inputfilenamecferr2Plus = basedir + "NtuplePlots_cferr2Plus.root";
  string InputfilenamehfMinus = basedir + "NtuplePlots_hfMinus.root";
  string InputfilenamehfPlus = basedir + "NtuplePlots_hfPlus.root";
  string Inputfilenamehfstats1Minus = basedir + "NtuplePlots_hfstats1Minus.root";
  string Inputfilenamehfstats1Plus = basedir + "NtuplePlots_hfstats1Plus.root";
  string Inputfilenamehfstats2Minus = basedir + "NtuplePlots_hfstats2Minus.root";
  string Inputfilenamehfstats2Plus = basedir + "NtuplePlots_hfstats2Plus.root";
  string InputfilenamelfMinus = basedir + "NtuplePlots_lfMinus.root";
  string InputfilenamelfPlus = basedir + "NtuplePlots_lfPlus.root";
  string Inputfilenamelfstats1Minus = basedir + "NtuplePlots_lfstats1Minus.root";
  string Inputfilenamelfstats1Plus = basedir + "NtuplePlots_lfstats1Plus.root";
  string Inputfilenamelfstats2Minus = basedir + "NtuplePlots_lfstats2Minus.root";
  string Inputfilenamelfstats2Plus = basedir + "NtuplePlots_lfstats2Plus.root";
  
  
  string Outputpath = "ErrorBand/";
  string Outputfilename = "ErrorBandFile_new.root";
  // mkdir(Outputpath.c_str(),0777);
  //bool mergeSignal = false;
  
  //systematics you want to consider in the errorband. The effects will be quadratically summed bin by bin (i.e. assumption = uncorrelated systematics, which is in most cases rather unrealistic).
  vector<string> syst;
  syst.push_back("JES");
  syst.push_back("JER");
  syst.push_back("cferr1");
  syst.push_back("cferr2");
  syst.push_back("hf");
  syst.push_back("hfstats1");
  syst.push_back("hfstats2");
  syst.push_back("lf");
  syst.push_back("lfstats1");
  syst.push_back("lfstats2");
  syst.push_back("PU");
  syst.push_back("MuonSF");
  syst.push_back("ElectronSF");
  // rates
  syst.push_back("zz_XS");
  syst.push_back("tzq_XS");
  syst.push_back("ttz_XS");
  syst.push_back("other_XS");
  syst.push_back("fake_XS");
  syst.push_back("wz_XS");
  syst.push_back("lumi");
  
  
  //numbers are to be revisited
  float rel_unc_wz = 0.30;
  float rel_unc_zz = 0.30;
  float rel_unc_tzq = 0.30;
  float rel_unc_other = 0.30;
  float rel_unc_ttz = 0.30;
  float rel_unc_fake = 0.50; //
  //now some overall scale systematics
  float rel_unc_lumi = 0.025; //2016
  
  
  bool doIgnore = true; //true: ignore MSPlots with names including one of the strings in MSPlotsIgnore, false: strings in MSPlotsIgnore will not be used to ignore certain MSPlots.
  vector<string> MSPlotsIgnore;
  //MSPlotsIgnore.push_back("BtagWeight");
  //bool doOnlyConsider = true; //true: only consider MSPlots with names including one of the strings in MSPlotsIgnore, false: strings in MSPlotsIgnore will not be used to 'only consider' certain MSPlots
  vector<string> MSPlotsOnlyConsider;
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_mWt");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_SMtop_eta");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_mlb");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dPhiWlepb");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_deltaRWlepJet_min");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_Zboson_M");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dPhiZWlep");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dRWlepb");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_NJets_CSVv2M");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_FCNCtop_M");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dRZc");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dRZb");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dRSMjetLightjet");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_charge_asym");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_bdiscCSVv2_jet_0");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_TotalHt_lep");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_ptWQ");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_TotalInvMass");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_TotalInvMass_lep");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_dRZWlep");
  MSPlotsOnlyConsider.push_back("wzcontrol_MVA_Bdis_LightJet");
  
  //note: if plot name should both be 'ignored' and 'only considered', it will be ignored.
  
  
  TFile* outFile = new TFile((Outputpath+Outputfilename).c_str(),"RECREATE");
  
  TFile* inFileNominal = new TFile(InputfilenameNominal.c_str(),"READ");
  if(verbose >= 1) cout<<"Nominal file read = "<<inFileNominal->GetName()<<endl;
  TIter nextkey(inFileNominal->GetListOfKeys());
  TKey *key;
  if(verbose >= 1)
  {
	   cout<<"Loop over de directories (~the multisampleplots) in the file"<<endl;
	   cout<<"------------------------------------------------------------"<<endl;
  }
  while  ((key = (TKey*)nextkey()))
  {
    TH1F* sumNominal = 0;
    TH1F* errorMinus = 0; //lower values of error band
    TH1F* errorPlus = 0; //upper values of error band
				
    TDirectoryFile* subdir = (TDirectoryFile*) key->ReadObj();
    string subdirname = subdir->GetName();
    if(verbose >= 1) cout<<"*** "<<subdirname<<" ***"<<endl;
    
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
        if(verbose >= 1) cout<<" -> IGNORING "<<subdirname<<endl;
      }
      else
      {
        TIter nextplotkey(subdir->GetListOfKeys());
        TKey *plotkey;
        if(verbose >= 1) cout<<"loop over histograms (in multisampleplot directory) and obtain total nominal histogram"<<endl;
        bool firstdataset = true;
        while ((plotkey = (TKey*)nextplotkey()))
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
        
        if(verbose >= 1) cout<<"loop over systematics"<<endl;
        for(unsigned int iSyst=0; iSyst<syst.size(); iSyst++)
        {
          if(verbose >= 1) cout<<" -> systematic "<<syst[iSyst]<<endl;
          if(verbose >= 1) cout<<"    up variations"<<endl;
          TFile* inFilePlus = 0;
          if(syst[iSyst] == "JES")            inFilePlus = new TFile(InputfilenameJESPlus.c_str(),"READ");
          else if(syst[iSyst] == "JER")       inFilePlus = new TFile(InputfilenameJERPlus.c_str(),"READ");
          else if(syst[iSyst] == "bTag")      inFilePlus = new TFile(InputfilenamebTagPlus.c_str(),"READ");
          else if(syst[iSyst] == "PU")        inFilePlus = new TFile(InputfilenamePUPlus.c_str(),"READ");
          else if(syst[iSyst] == "MuonSF")    inFilePlus = new TFile(InputfilenameMuonSFPlus.c_str(),"READ");
          else if(syst[iSyst] == "ElectronSF")inFilePlus = new TFile(InputfilenameElectronSFPlus.c_str(),"READ");
          else if(syst[iSyst] ==  "cferr1"  ) inFilePlus = new TFile(Inputfilenamecferr1Plus.c_str(),"READ");
          else if(syst[iSyst] == "cferr2"   ) inFilePlus = new TFile(Inputfilenamecferr2Plus.c_str(),"READ");
          else if(syst[iSyst] == "hf"     )   inFilePlus = new TFile(InputfilenamehfPlus.c_str(),"READ");
          else if(syst[iSyst] ==  "hfstats1") inFilePlus = new TFile(Inputfilenamehfstats1Plus.c_str(),"READ");
          else if(syst[iSyst] ==  "hfstats2") inFilePlus = new TFile(Inputfilenamehfstats2Plus.c_str(),"READ");
          else if(syst[iSyst] ==  "lf"    )   inFilePlus = new TFile(InputfilenamelfPlus.c_str(),"READ");
          else if(syst[iSyst] == "lfstats1" ) inFilePlus = new TFile(Inputfilenamelfstats1Plus.c_str(),"READ");
          else if(syst[iSyst] ==  "lfstats2") inFilePlus = new TFile(Inputfilenamelfstats2Plus.c_str(),"READ");
          else                                inFilePlus = new TFile(InputfilenameNominal.c_str(),"READ");
          
          if(verbose >= 1) cout<<"    Plusfile read = "<<inFilePlus->GetName()<<endl;
          TDirectory* subdirPlus = (TDirectory*) inFilePlus->Get(subdirname.c_str());
          TH1F* sumPlus = 0;
          
          TIter nextplotkeyPlus(subdirPlus->GetListOfKeys());
          TKey *plotkeyPlus;
          bool firstdatasetPlus = true;
          
          // shapes needed
          if(syst[iSyst] == "JES"  || syst[iSyst] == "JER" || syst[iSyst] == "cferr1" || syst[iSyst] == "cferr2" || syst[iSyst] == "hf" || syst[iSyst] == "hfstats1" || syst[iSyst] == "hfstats2" || syst[iSyst] == "lf" || syst[iSyst] == "lfstats1" || syst[iSyst] == "lfstats2" ||  syst[iSyst] == "PU" || syst[iSyst] == "MuonSF" || syst[iSyst] == "ElectronSF") //WARNING: HARDCODED
          {
            if(verbose >= 1) cout<<"    loop over histograms (in multisampleplot directory) and obtain total up-varied histogram"<<endl;
            while ((plotkeyPlus = (TKey*)nextplotkeyPlus()))
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
          else  // rate unc
          {
            if(verbose >= 1) cout<<"   loop over histograms (in multisampleplot directories) and obtain total up-varied histogram"<<endl;
            while ((plotkeyPlus = (TKey*)nextplotkeyPlus()))
            {
              TH1F* histoPlus = 0;
              if(((string) (plotkeyPlus->ReadObj())->ClassName()) == "TH1F")
              {
                histoPlus = (TH1F*) plotkeyPlus->ReadObj();
                if(!( ((string) histoPlus->GetName()).find("Data")!=string::npos || ((string) histoPlus->GetName()).find("data")!=string::npos || ((string) histoPlus->GetName()).find("DATA")!=string::npos || ((string) histoPlus->GetName()).find("_NP_")!=string::npos ) )
                {
                  if(syst[iSyst] == "wz_XS" && ((string) histoPlus->GetName()).find("WZ")!=string::npos && ((string) histoPlus->GetName()).find("WWZ")==string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos)
			               {
                       histoPlus->Scale(1+rel_unc_wz);
                     }
                  else if(syst[iSyst] == "ttz_XS" && ((string) histoPlus->GetName()).find("TTZ")!=string::npos)
			               {
                       histoPlus->Scale(1+rel_unc_ttz);
                     }
                  else if( syst[iSyst] == "zz_XS" && ((string) histoPlus->GetName()).find("ZZ")!=string::npos && ((string) histoPlus->GetName()).find("ZZZ")==string::npos && ((string) histoPlus->GetName()).find("WZZ")==string::npos )
			               {
                       histoPlus->Scale(1+rel_unc_zz);
                     }
                  else if(syst[iSyst] == "tzq_XS" && ((string) histoPlus->GetName()).find("tZq")!=string::npos)
			               {
                       histoPlus->Scale(1+rel_unc_ttz);
                     }
                  else if(syst[iSyst] == "other_XS" && ((string) histoPlus->GetName()).find("TTZ")==string::npos &&  ((string) histoPlus->GetName()).find("tZq")==string::npos   && ((string) histoPlus->GetName()).find("TTZ")!=string::npos  && ((string) histoPlus->GetName()).find("WZT")!=string::npos  && ((string) histoPlus->GetName()).find("ZZT")!=string::npos)
			               {
                       histoPlus->Scale(1+rel_unc_other);
                     }
                  else if(syst[iSyst] == "fake_XS" &&  ((string) histoPlus->GetName()).find("fake")!=string::npos  )
			               {
                       histoPlus->Scale(1+rel_unc_fake);
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
          if(syst[iSyst] == "JES")              inFileMinus = new TFile(InputfilenameJESMinus.c_str(),"READ");
          else if(syst[iSyst] == "JER")         inFileMinus = new TFile(InputfilenameJERMinus.c_str(),"READ");
          else if(syst[iSyst] == "bTag")        inFileMinus = new TFile(InputfilenamebTagMinus.c_str(),"READ");
          else if(syst[iSyst] == "PU")          inFileMinus = new TFile(InputfilenamePUMinus.c_str(),"READ");
          else if(syst[iSyst] == "MuonSF")      inFileMinus = new TFile(InputfilenameMuonSFMinus.c_str(),"READ");
          else if(syst[iSyst] == "ElectronSF")  inFileMinus = new TFile(InputfilenameElectronSFMinus.c_str(),"READ");
          else if(syst[iSyst] ==  "cferr1"    ) inFileMinus = new TFile(Inputfilenamecferr1Minus.c_str(),"READ");
          else if(syst[iSyst] == "cferr2"     ) inFileMinus = new TFile(Inputfilenamecferr2Minus.c_str(),"READ");
          else if(syst[iSyst] == "hf"     )     inFileMinus = new TFile(InputfilenamehfMinus.c_str(),"READ");
          else if(syst[iSyst] ==  "hfstats1"  ) inFileMinus = new TFile(Inputfilenamehfstats1Minus.c_str(),"READ");
          else if(syst[iSyst] ==  "hfstats2"  ) inFileMinus = new TFile(Inputfilenamehfstats2Minus.c_str(),"READ");
          else if(syst[iSyst] ==  "lf"    )     inFileMinus = new TFile(InputfilenamelfMinus.c_str(),"READ");
          else if(syst[iSyst] == "lfstats1"   ) inFileMinus = new TFile(Inputfilenamelfstats1Minus.c_str(),"READ");
          else if(syst[iSyst] ==  "lfstats2"  ) inFileMinus = new TFile(Inputfilenamelfstats2Minus.c_str(),"READ");
          else                                  inFileMinus = new TFile(InputfilenameNominal.c_str(),"READ");
          if(verbose >= 1) cout<<"   Minusfile read = "<<inFileMinus->GetName()<<endl;
          TDirectory* subdirMinus = (TDirectory*) inFileMinus->Get(subdirname.c_str());
          TH1F* sumMinus = 0;
          
          TIter nextplotkeyMinus(subdirMinus->GetListOfKeys());
          TKey *plotkeyMinus;
          bool firstdatasetMinus = true;
          
          if(syst[iSyst] == "JES"  || syst[iSyst] == "JER" || syst[iSyst] == "cferr1" || syst[iSyst] == "cferr2" || syst[iSyst] == "hf" || syst[iSyst] == "hfstats1" || syst[iSyst] == "hfstats2" || syst[iSyst] == "lf" || syst[iSyst] == "lfstats1" || syst[iSyst] == "lfstats2" || syst[iSyst] == "PU" || syst[iSyst] == "MuonSF" || syst[iSyst] == "ElectronSF") // WARNING HARDCODED
          {
            if(verbose >= 1) cout<<"   loop over histograms (in multisampleplot directories) and obtain total down-varied histogram"<<endl;
            while ( (plotkeyMinus = (TKey*)nextplotkeyMinus()) )
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
            while ((plotkeyMinus = (TKey*)nextplotkeyMinus()))
            {
              TH1F* histoMinus = 0;
              if(((string) (plotkeyMinus->ReadObj())->ClassName()) == "TH1F")
              {
                histoMinus = (TH1F*) plotkeyMinus->ReadObj();
                if(!( ((string) histoMinus->GetName()).find("Data")!=string::npos || ((string) histoMinus->GetName()).find("data")!=string::npos || ((string) histoMinus->GetName()).find("DATA")!=string::npos || ((string) histoMinus->GetName()).find("_NP_")!=string::npos ) )
                {
                  if(syst[iSyst] == "wz_XS" && ((string) histoMinus->GetName()).find("WZ")!=string::npos && ((string) histoMinus->GetName()).find("WWZ")==string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos)
			               {
                       histoMinus->Scale(1-rel_unc_wz);
                     }
                  else if(syst[iSyst] == "ttz_XS" && ((string) histoMinus->GetName()).find("TTZ")!=string::npos)
			               {
                       histoMinus->Scale(1-rel_unc_ttz);
                     }
                  else if( syst[iSyst] == "zz_XS" && ((string) histoMinus->GetName()).find("ZZ")!=string::npos && ((string) histoMinus->GetName()).find("ZZZ")==string::npos && ((string) histoMinus->GetName()).find("WZZ")==string::npos )
			               {
                       histoMinus->Scale(1-rel_unc_zz);
                     }
                  else if(syst[iSyst] == "tzq_XS" && ((string) histoMinus->GetName()).find("tZq")!=string::npos)
			               {
                       histoMinus->Scale(1-rel_unc_ttz);
                     }
                  else if(syst[iSyst] == "other_XS" && ((string) histoMinus->GetName()).find("TTZ")==string::npos && ((string) histoMinus->GetName()).find("tZq")==string::npos   && ((string) histoMinus->GetName()).find("TTZ")!=string::npos  && ((string) histoMinus->GetName()).find("WZT")!=string::npos  && ((string) histoMinus->GetName()).find("ZZT")!=string::npos)
			               {
                       histoMinus->Scale(1-rel_unc_other);
                     }
                  else if(syst[iSyst] == "fake_XS" &&  ((string) histoMinus->GetName()).find("fake")!=string::npos  )
			               {
                       histoMinus->Scale(1-rel_unc_fake);
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
  return 0;
}
