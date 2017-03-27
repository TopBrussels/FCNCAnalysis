#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include "TStyle.h"
#include "TPaveText.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"
#include <TFile.h>
#include <TLeaf.h>
#include <utility>
#include "Style.C"
#include <map>

// used TopTreeAnalysis classes
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"

#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/LHAPDF.h"
#include "/user/kderoove/Software/LHAPDF/LHAPDF-6.1.6/include/LHAPDF/Reweighting.h"

using namespace std;
using namespace TopTree;


///////////////////////////////////// PLOT MAPPING /////////////////////////////////////////
// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,MultiSamplePlot*> MSPlot;

////////////////////////////////// mapping ///////////////////////////////
map<string,TFile*> tFileMap;
map<string,TTree*> tTree;
map<string,TTree*> tStatsTree;
vector < Dataset* > datasets;

////////////////////////////////// functions ////////////////////////////////////////////
// bookkeeping
string ConvertIntToString(int Number, int pad);
string MakeTimeStamp();
string intToStr (int number);
double maximumValue(vector<double> array);
double minimumValue(vector<double> array);
// initialisations
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair);
void InitCalculatePDFWeightHisto(string dataSetName);
void InitSystematicHisto(string dataSetName, string systematic);
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut);
// functions
vector<double> BDTCUT(string region, string coupling);
void CalculatePDFWeight(string dataSetName, double BDT, double MVA_weight_nom, int MVA_channel);
void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData, bool toppair, double weight_, int MVA_channel);
void GetPDFEnvelope(string dataSetName);

//////////////////////////////// settings ////////////////////////////////
bool makePlots = false;
bool PlotSystematics = false;
bool doPseudoData = false;
bool doSystematics = false;
int channel = -999;
bool datafound = false;
bool testing = false;
bool doZut = false;
bool toppair = false;
bool addData = false;
bool DetermineCut = false;
bool doPDFunc  = false;
string placeNtup = "";
string tempstring = "";
string systematic = "";
string decaystring = "";
string coupling = "Zct";
string region = "singletop";
string placeOutputReading = "";
string combinetemplate_filename = "";
Double_t Luminosity = 36000.;
string dataSetName = "";
bool isData = false;
string tTreeName = "";
string postfix = "";
string output_histo_name = "";
string ntupleFileName ="";
int nbin = 100;
int nEntries = -1;


//////////////////////////// branches //////////////////////////////
// Declaration of leaf types
Float_t         MVA_Zboson_pt;
Float_t         MVA_dPhiWlepb;
Float_t         MVA_charge_asym;
Float_t         MVA_bdiscCSVv2_jet_0;

Float_t         MVA_cdiscCvsB_jet_1;
Float_t         MVA_cdiscCvsB_jet_0;
Float_t         MVA_dRZc;
Float_t         MVA_dRWlepb;
Float_t         MVA_dRZWlep;
Float_t         MVA_mlb;
Float_t         MVA_FCNCtop_M;
Float_t           MVA_nJets_CharmL;
Float_t           MVA_NJets_CSVv2M;
Int_t         MVA_region;
Double_t       MVA_x1;
Double_t       MVA_x2;
Int_t       MVA_id1;
Int_t       MVA_id2;
Double_t       MVA_q;
//Float_t         MVA_weight;
Double_t       MVA_weight_nom;
Int_t         MVA_channel;
Double_t       MVA_BDT;
Double_t        MVA_EqLumi;
Double_t       MVA_weight_puSF_up;
Double_t       MVA_weight_puSF_down;
Double_t       MVA_weight_electronSF_up;
Double_t       MVA_weight_electronSF_down;
Double_t       MVA_weight_muonSF_up;
Double_t       MVA_weight_muonSF_down;
Double_t       MVA_weight_btagSF_cferr1_up;
Double_t       MVA_weight_btagSF_cferr1_down;
Double_t       MVA_weight_btagSF_cferr2_up;
Double_t       MVA_weight_btagSF_cferr2_down;
Double_t       MVA_weight_btagSF_hf_up;
Double_t       MVA_weight_btagSF_hf_down;
Double_t       MVA_weight_btagSF_hfstats1_up;
Double_t       MVA_weight_btagSF_hfstats1_down;
Double_t       MVA_weight_btagSF_hfstats2_up;
Double_t       MVA_weight_btagSF_hfstats2_down;
Double_t       MVA_weight_btagSF_lf_up;
Double_t       MVA_weight_btagSF_lf_down;
Double_t       MVA_weight_btagSF_lfstats1_up;
Double_t       MVA_weight_btagSF_lfstats1_down;
Double_t       MVA_weight_btagSF_lfstats2_up;
Double_t       MVA_weight_btagSF_lfstats2_down;

// List of branches
TBranch        *b_MVA_Zboson_pt;   //!
TBranch        *b_MVA_dPhiWlepb;   //!
TBranch        *b_MVA_charge_asym;   //!
TBranch        *b_MVA_bdiscCSVv2_jet_0;   //!

TBranch        *b_MVA_cdiscCvsB_jet_1;   //!
TBranch        *b_MVA_cdiscCvsB_jet_0;   //!
TBranch        *b_MVA_dRZc;   //!
TBranch        *b_MVA_dRWlepb;   //!
TBranch        *b_MVA_dRZWlep;   //!
TBranch        *b_MVA_mlb;   //!
TBranch        *b_MVA_FCNCtop_M;   //!
TBranch        *b_MVA_nJets_CharmL;   //!
TBranch        *b_MVA_NJets_CSVv2M;   //!
TBranch        *b_MVA_region;   //!
TBranch        *b_MVA_x1;   //!
TBranch        *b_MVA_x2;   //!
TBranch        *b_MVA_id1;   //!
TBranch        *b_MVA_id2;   //!
TBranch        *b_MVA_q;   //!
TBranch        *b_MVA_weight;   //!
TBranch        *b_MVA_channel;   //!
TBranch        *b_MVA_BDT;   //!
TBranch        *b_MVA_EqLumi;   //!
TBranch        *b_MVA_weight_puSF_up;   //!
TBranch        *b_MVA_weight_puSF_down;   //!
TBranch        *b_MVA_weight_electronSF_up;   //!
TBranch        *b_MVA_weight_electronSF_down;   //!
TBranch        *b_MVA_weight_muonSF_up;   //!
TBranch        *b_MVA_weight_muonSF_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr1_down;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_up;   //!
TBranch        *b_MVA_weight_btagSF_cferr2_down;   //!
TBranch        *b_MVA_weight_btagSF_hf_up;   //!
TBranch        *b_MVA_weight_btagSF_hf_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_hfstats2_down;   //!
TBranch        *b_MVA_weight_btagSF_lf_up;   //!
TBranch        *b_MVA_weight_btagSF_lf_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats1_down;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_up;   //!
TBranch        *b_MVA_weight_btagSF_lfstats2_down;   //!

///////////////////////////////////// MAIN CODE /////////////////////////////////////////
int main(int argc, char* argv[]){
  string dateString = MakeTimeStamp();
  cout << "*********************************************" << endl;
  cout << "***         Beginning of program          ***" << endl;
  cout << "*********************************************" << endl;
  cout << "*********************************************" << endl;
  cout << "* Current time: " << dateString << "                 *" << endl;
  cout << "*********************************************" << endl;
  
  clock_t start = clock();
  //setTDRStyle(); // TO FIX
  //setMyStyle(); // TO FIX stat box title
  
  
  int testnr = -1;
  //////////// Settings of the analysis //////////////////
  for(int i = 0; i <argc; i++){
    if(string(argv[i]).find("help")!=string::npos) {
      std::cout << "****** help ******" << endl;
      std::cout << " run code with ./AnalyzerBDT [options]" << endl;
      std::cout << "Options: " << endl;
      std::cout << "   Zct / Zut: make plots for this one" << endl;
      std::cout << "   toppair: make plots for this one" << endl;
      std::cout << "   MakePlots: make plots" << endl;
      std::cout << "   Ntup placeNtup: set where ntuples are stored. MVAoutput/placeNtup " << endl;
      std::cout << "   Test nb: loop over nb events" << endl;
      std::cout << "   Data: add CRcut" << endl;
      std::cout << "   PSdata: generate pseudo data" << endl;
      std::cout << "   Systematics: loop over systematics" << endl;
      std::cout << "   doPDFunc: calculate PDF unc" << endl;
      std::cout << "   PlotSystematics: calculate PDF unc" << endl;
      return 0;
    }
    if(string(argv[i]).find("doPDFunc")!=std::string::npos){
      doPDFunc = true;
      cout << "******* calculating pdf uncertainties *********" << endl;
    }
    if(string(argv[i]).find("Systematics")!=std::string::npos) {
      doSystematics= true;
    }
    if(string(argv[i]).find("PlotSystematics")!=std::string::npos) {
      PlotSystematics= true;
    }
    if(string(argv[i]).find("PSdata")!=std::string::npos) {
      doPseudoData = true;
    }
    if(string(argv[i]).find("Data")!=std::string::npos) {
      addData = true;
    }
    if(string(argv[i]).find("toppair")!=std::string::npos) {
      toppair = true;
      region = "toppair";
    }
    if(string(argv[i]).find("Zut")!=std::string::npos) {
      doZut = true;
      coupling = "Zut";
    }
    if(string(argv[i]).find("Test")!=std::string::npos) {
      testing = true;
      testnr = strtol(argv[i+1], NULL, 10);
      i++;
    }
    if(string(argv[i]).find("Ntup")!=std::string::npos) {
      placeNtup = argv[i+1];
      i++;
    }
    if(string(argv[i]).find("MakePlots")!=string::npos) {
      makePlots = true;
    }
    if(string(argv[i]).find("DetermineCut")!=string::npos) {
      DetermineCut= true;
    }
    
  }
  string xmlFileName = "";
  xmlFileName = "config/Run2TriLepton_samples_analyBDT.xml" ;
  const char* xmlFile = xmlFileName.c_str();
  cout << " - Using config file " << xmlFile << endl;
  cout << " - Using mvatrees of " << placeNtup << endl;
  placeOutputReading = "MVAoutput/";
  mkdir(placeOutputReading.c_str(), 0777);
  placeOutputReading += "outputtemplates";
  mkdir(placeOutputReading.c_str(), 0777);
  placeOutputReading += "/" + coupling + "_" + region;
  mkdir(placeOutputReading.c_str(), 0777);
  combinetemplate_filename = placeOutputReading+"/Reader_" + coupling +"_" + region + ".root";
  cout <<" - Combine templates stored at " << combinetemplate_filename.c_str() << endl;
  
  std::vector < int>  decayChannels = {0,1,2,3,-9}; // uuu uue eeu eee all
  vector <string> thesystlist;
  thesystlist.clear();
  thesystlist.push_back(""); // nominal
  
  ///////////////// Initialisation ////////////////////
  if(makePlots){
    
       tempstring = region + "_"+coupling;
       InitMSPlots(tempstring, decayChannels, toppair);
    
  }
  
  
  
  
  ///////////////:  load datasets
  datasets.clear();
  TTreeLoader treeLoader;
  cout << "loading " << endl;
  treeLoader.LoadDatasets(datasets, xmlFile);
  cout << "datasets loaded" <<endl;
  for (int d = 0; d < datasets.size(); d++){   //Loop through datasets to get lumi setting
    
    dataSetName = datasets[d]->Name();
    if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  )
    {
      Luminosity = datasets[d]->EquivalentLumi();
      datafound = true;
      cout << "data found" <<endl;
    }
  }
  
  
  //////////// START LOOPING ON SYS - DATASETS - EVENTS //////////
  ///*****************///
  ///   MAIN CODE  ///
  ///*****************///
  
    string systematic = "";
    MSPlot["test"] = new MultiSamplePlot(datasets, "test", 10, -1.,1., "test");
    for (int d = 0; d < datasets.size(); d++)   //Loop through datasets
    {
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " / title : " << datasets[d]->Title() << endl;
      // settings
      isData = false;
      dataSetName = datasets[d]->Name();
      if (dataSetName.find("Data")!=std::string::npos || dataSetName.find("data")!=std::string::npos|| dataSetName.find("DATA")!=std::string::npos  ){
        isData = true;
      }
      ntupleFileName = placeNtup;
      tFileMap[dataSetName.c_str()] = new TFile((ntupleFileName).c_str(),"READ"); //create TFile for each dataset
      
      postfix = "";
      if(systematic.find("JES_down")!=std::string::npos) postfix = "_JESdown";
      else if(systematic.find("JES_up")!=std::string::npos) postfix = "_JESup";
      else if(systematic.find("JER_down")!=std::string::npos) postfix = "_JERdown";
      else if(systematic.find("JER_up")!=std::string::npos) postfix = "_JERup";
      else postfix = "";
      tTreeName = "Control_"+dataSetName + postfix;
      /// Get data
      cout << "   treename " << tTreeName << " from " << ntupleFileName <<  endl;
      tTree[dataSetName.c_str()] = (TTree*)tFileMap[dataSetName.c_str()]->Get(tTreeName.c_str()); //get ttree for each dataset
      nEntries = -1;
      nEntries = (int)tTree[dataSetName.c_str()]->GetEntries();
      cout << "                nEntries: " << nEntries << endl;
      
      // Initialise tree
      InitTree(tTree[dataSetName.c_str()], isData, toppair, doZut);
      
      
      

      
      int endEvent =nEntries;
      if(testing){
        if(endEvent > testnr) endEvent = testnr;
      }
      /// loop on events
      double weight = 1.;
      for (int ievt = 0; ievt < endEvent; ievt++)
      {
        if (ievt%100 == 0)
          std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)nEntries)*100  << "%)" << flush << "\r";
        
        
        /// Load event
        tTree[(dataSetName).c_str()]->GetEntry(ievt);
        
        // cout << "mva channel " << MVA_channel << endl;
        
        weight = MVA_weight_nom;
        
        // if(addData && BDT < cut) continue;
        if (makePlots)
        {
          //cout << "ievt " << ievt << endl;
          tempstring = region + "_"+coupling;
          
          FillGeneralPlots(d, tempstring, decayChannels, isData, toppair, weight, MVA_channel);
        }
        MSPlot["test"]->Fill(0., datasets[d], true, 1.);
        
      } // events
      cout << endl;
      
        }// datasets
    
  
  
  ///*****************///
  ///   Write plots   ///
  ///*****************///
  if(makePlots || doPDFunc){
    string pathOutput = "OutputPlots/";
    mkdir(pathOutput.c_str(),0777);
    string pathOutputdate = pathOutput + dateString + "/"  ;
    mkdir(pathOutputdate.c_str(),0777);
    string place =pathOutputdate+"MSPlot/";
    string placeTH1F = pathOutputdate+"TH1F/";
    string placeTH2F = pathOutputdate+"TH2F/";
    vector <string> vlabel_chan = {"3#mu", "1e2#mu", "2e1#mu", "3e"};
    mkdir(place.c_str(),0777);
    mkdir(placeTH1F.c_str(),0777);
    mkdir(placeTH2F.c_str(),0777);
    string rootFileName ="NtuplePlotsMVA.root";
    
    cout << " - Recreate output file ..." << endl;
    TFile *fout = new TFile ((pathOutputdate+rootFileName).c_str(), "RECREATE");
    cout << "   Output file is " << pathOutputdate+rootFileName << endl;
    
    ///Write histograms
    fout->cd();
    if(makePlots){
      for (map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
      {
        cout << "MSPlot: " << it->first << endl;
        MultiSamplePlot *temp = it->second;
        string name = it->first;
        if(!datafound) temp->setDataLumi(Luminosity);
        /*if(name.find("all")!=std::string::npos) temp->setChannel(true, "all");
         if(name.find("eee")!=std::string::npos) temp->setChannel(true, "3e");
         if(name.find("eeu")!=std::string::npos) temp->setChannel(true, "2e1#mu");
         if(name.find("uue")!=std::string::npos) temp->setChannel(true, "1e2#mu");
         if(name.find("uuu")!=std::string::npos) temp->setChannel(true, "3#mu");
         if(name.find("Decay")!=std::string::npos) temp->setBins(vlabel_chan);*/
        temp->Draw(name, 1, false, false, false, 10);  // string label, unsigned int RatioType, bool addRatioErrorBand, bool addErrorBand, bool ErrorBandAroundTotalInput, int scaleNPSignal
        cout << "writing to " << pathOutputdate+"MSPlot" << endl;
        cout << "plot " << name << endl;
        cout << "temp " << temp << endl;
        temp->Write(fout, name, true, (pathOutputdate+"MSPlot").c_str(), "png");  // TFile* fout, string label, bool savePNG, string pathPNG, string ext
      }
      
    }
    TDirectory* th1dir = fout->mkdir("1D_histograms");
    th1dir->cd();
    gStyle->SetOptStat(1110);
    for (std::map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
    {
      TH1F *temp = it->second;
      int N = temp->GetNbinsX();
      temp->SetBinContent(N,temp->GetBinContent(N)+temp->GetBinContent(N+1));
      temp->SetBinContent(N+1,0);
      temp->SetEntries(temp->GetEntries()-2); // necessary since each SetBinContent adds +1 to the number of entries...
      temp->Write();
      if(it->first.find("nominal")!=std::string::npos || it->first.find("Up")!=std::string::npos || it->first.find("Down")!=std::string::npos) {
        TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
        tempCanvas->SaveAs( (placeTH1F+it->first+".png").c_str() );
      }
    }
    

    fout->Write();
    fout->Close();
    
    delete fout;
    
  }
  
  
  
  
  
  double time = ((double)clock() - start) / CLOCKS_PER_SEC;
  cout << "It took us " << time << " s to run the program" << endl;
  if ( time >= 60 )
  {
    int mins = time/60;
    float secs = time - mins*60;
    
    if (mins >= 60 )
    {
      int hours = mins/60;
      mins = mins - hours*60;
      cout << "(This corresponds to " << hours << " hours, " << mins << " min and " << secs << " s)" << endl;
    }
    else
      cout << "(This corresponds to " << mins << " min and " << secs << " s)" << endl;
  }
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  cout << " - Goodbye" << endl;
  return 0;
}



///// BOOK KEEPING FUNCTIONS
string ConvertIntToString(int Number, int pad){
  ostringstream convert;
  convert.clear();  // clear bits
  convert.str(std::string());  // clear content
  if ( pad > 1 ) { convert << std::setw(pad) << std::setfill('0');}
  convert << Number;
  return convert.str();
}
string MakeTimeStamp(){
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  int year = now->tm_year - 100;  /// + 1900 to get current year
  int month = now->tm_mon + 1;
  int day = now->tm_mday;
  int hour = now->tm_hour;
  int min = now->tm_min;
  //int sec = now->tm_sec;
  
  string year_str = ConvertIntToString(year, 2);
  string month_str = ConvertIntToString(month, 2);
  string day_str = ConvertIntToString(day, 2);
  string hour_str = ConvertIntToString(hour, 2);
  string min_str = ConvertIntToString(min, 2);
  //string sec_str = ConvertIntToString(sec, 2);
  
  string date_str = year_str + month_str + day_str + "_" + hour_str + min_str;
  return date_str;
}
string intToStr (int number){
  ostringstream buff;
  buff<<number;
  return buff.str();
}
double maximumValue(vector<double> array){
  int length = array.size();  // establish size of array
  double max = array[0];       // start with max = first element
  
  for(int i = 1; i<length; i++)
  {
    if(array[i] > max)
      max = array[i];
  }
  return max;                // return highest value in array
}
double minimumValue(vector<double> array){
  int length = array.size();  // establish size of array
  double max = array[0];       // start with max = first element
  
  for(int i = 1; i<length; i++)
  {
    if(array[i] < max)
      max = array[i];
  }
  return max;                // return highest value in array
}

//// INITIALISATIONS
void InitMSPlots(string prefix, vector <int> decayChannels , bool istoppair){
  clock_t start_sub = clock();
  
  
  
  // control plots
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    
    
    //cout << "init " << (prefix+"_BDT_"+decaystring).c_str() << endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()] = new MultiSamplePlot(datasets, (prefix+"_BDT_"+decaystring).c_str(), nbin, -1.,1., "BDT");
    
  }
  
  
  
  
  
}
void InitCalculatePDFWeightHisto(string dataSetName){
  TH1::SetDefaultSumw2();
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  for(int iChan = 0; iChan < channel_list.size(); iChan++){
    channel = channel_list[iChan];
    for ( int i=0; i<101; i++)
    {
      output_histo_name = dataSetName+"_BDT_"+channel+"_"+intToStr(i);
      histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    }
    output_histo_name = dataSetName+"_BDT_"+channel+"_nominal";
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    output_histo_name = dataSetName+"_BDT_"+channel+"_PDFEnvelopeUp";
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    output_histo_name = dataSetName+"_BDT_"+channel+"_PDFEnvelopeDown";
    histo1D[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
    
    output_histo_name = "";
  }
}
void InitSystematicHisto(string dataSetName, string systematic){
  TH1::SetDefaultSumw2();
  
  /*output_histo_name = dataSetName+"_BDT_nominal";
   histo1DSys[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
   output_histo_name = dataSetName+"_BDT_"+systematic;
   histo1DSys[output_histo_name] = new TH1F(output_histo_name.c_str(), output_histo_name.c_str(), nbin,-1.,1.);
   output_histo_name = "";*/
  
}
void InitTree(TTree* tree, bool isData, bool istoppair, bool doZut){
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetMakeClass(1);
  if(istoppair && doZut){
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_1", &MVA_cdiscCvsB_jet_1, &b_MVA_cdiscCvsB_jet_1);
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
    tree->SetBranchAddress("MVA_dRZc", &MVA_dRZc, &b_MVA_dRZc);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    tree->SetBranchAddress("MVA_FCNCtop_M", &MVA_FCNCtop_M, &b_MVA_FCNCtop_M);
    tree->SetBranchAddress("MVA_nJets_CharmL", &MVA_nJets_CharmL, &b_MVA_nJets_CharmL);
    tree->SetBranchAddress("MVA_NJets_CSVv2M", &MVA_NJets_CSVv2M, &b_MVA_NJets_CSVv2M);
  }
  else if(!istoppair && doZut){
    
    tree->SetBranchAddress("MVA_Zboson_pt", &MVA_Zboson_pt, &b_MVA_Zboson_pt);
    tree->SetBranchAddress("MVA_dRWlepb", &MVA_dRWlepb, &b_MVA_dRWlepb);
    tree->SetBranchAddress("MVA_dPhiWlepb", &MVA_dPhiWlepb, &b_MVA_dPhiWlepb);
    tree->SetBranchAddress("MVA_charge_asym", &MVA_charge_asym, &b_MVA_charge_asym);
    tree->SetBranchAddress("MVA_dRZWlep", &MVA_dRZWlep, &b_MVA_dRZWlep);
    tree->SetBranchAddress("MVA_bdiscCSVv2_jet_0", &MVA_bdiscCSVv2_jet_0, &b_MVA_bdiscCSVv2_jet_0);
    tree->SetBranchAddress("MVA_cdiscCvsB_jet_0", &MVA_cdiscCvsB_jet_0, &b_MVA_cdiscCvsB_jet_0);
    tree->SetBranchAddress("MVA_mlb", &MVA_mlb, &b_MVA_mlb);
    
  }
  
  tree->SetBranchAddress("MVA_region", &MVA_region, &b_MVA_region);
  //tree->SetBranchAddress("MVA_weight", &MVA_weight, &b_MVA_weight);
  tree->SetBranchAddress("MVA_channel", &MVA_channel, &b_MVA_channel);
  tree->SetBranchAddress("MVA_BDT", &MVA_BDT, &b_MVA_BDT);
  tree->SetBranchAddress("MVA_EqLumi", &MVA_EqLumi, &b_MVA_EqLumi);
  
  
  tree->SetBranchAddress("MVA_x1", &MVA_x1, &b_MVA_x1);
  tree->SetBranchAddress("MVA_x2", &MVA_x2, &b_MVA_x2);
  tree->SetBranchAddress("MVA_id1", &MVA_id1, &b_MVA_id1);
  tree->SetBranchAddress("MVA_id2", &MVA_id2, &b_MVA_id2);
  tree->SetBranchAddress("MVA_q", &MVA_q, &b_MVA_q);
  
  tree->SetBranchAddress("MVA_weight_puSF_up", &MVA_weight_puSF_up, &b_MVA_weight_puSF_up);
  tree->SetBranchAddress("MVA_weight_puSF_down", &MVA_weight_puSF_down, &b_MVA_weight_puSF_down);
  tree->SetBranchAddress("MVA_weight_electronSF_up", &MVA_weight_electronSF_up, &b_MVA_weight_electronSF_up);
  tree->SetBranchAddress("MVA_weight_electronSF_down", &MVA_weight_electronSF_down, &b_MVA_weight_electronSF_down);
  tree->SetBranchAddress("MVA_weight_muonSF_up", &MVA_weight_muonSF_up, &b_MVA_weight_muonSF_up);
  tree->SetBranchAddress("MVA_weight_muonSF_down", &MVA_weight_muonSF_down, &b_MVA_weight_muonSF_down);
  tree->SetBranchAddress("MVA_weight_btagSF_cferr1_up", &MVA_weight_btagSF_cferr1_up, &b_MVA_weight_btagSF_cferr1_up);
  tree->SetBranchAddress("MVA_weight_btagSF_cferr1_down", &MVA_weight_btagSF_cferr1_down, &b_MVA_weight_btagSF_cferr1_down);
  tree->SetBranchAddress("MVA_weight_btagSF_cferr2_up", &MVA_weight_btagSF_cferr2_up, &b_MVA_weight_btagSF_cferr2_up);
  tree->SetBranchAddress("MVA_weight_btagSF_cferr2_down", &MVA_weight_btagSF_cferr2_down, &b_MVA_weight_btagSF_cferr2_down);
  tree->SetBranchAddress("MVA_weight_btagSF_hf_up", &MVA_weight_btagSF_hf_up, &b_MVA_weight_btagSF_hf_up);
  tree->SetBranchAddress("MVA_weight_btagSF_hf_down", &MVA_weight_btagSF_hf_down, &b_MVA_weight_btagSF_hf_down);
  tree->SetBranchAddress("MVA_weight_btagSF_hfstats1_up", &MVA_weight_btagSF_hfstats1_up, &b_MVA_weight_btagSF_hfstats1_up);
  tree->SetBranchAddress("MVA_weight_btagSF_hfstats1_down", &MVA_weight_btagSF_hfstats1_down, &b_MVA_weight_btagSF_hfstats1_down);
  tree->SetBranchAddress("MVA_weight_btagSF_hfstats2_up", &MVA_weight_btagSF_hfstats2_up, &b_MVA_weight_btagSF_hfstats2_up);
  tree->SetBranchAddress("MVA_weight_btagSF_hfstats2_down", &MVA_weight_btagSF_hfstats2_down, &b_MVA_weight_btagSF_hfstats2_down);
  tree->SetBranchAddress("MVA_weight_btagSF_lf_up", &MVA_weight_btagSF_lf_up, &b_MVA_weight_btagSF_lf_up);
  tree->SetBranchAddress("MVA_weight_btagSF_lf_down", &MVA_weight_btagSF_lf_down, &b_MVA_weight_btagSF_lf_down);
  tree->SetBranchAddress("MVA_weight_btagSF_lfstats1_up", &MVA_weight_btagSF_lfstats1_up, &b_MVA_weight_btagSF_lfstats1_up);
  tree->SetBranchAddress("MVA_weight_btagSF_lfstats1_down", &MVA_weight_btagSF_lfstats1_down, &b_MVA_weight_btagSF_lfstats1_down);
  tree->SetBranchAddress("MVA_weight_btagSF_lfstats2_up", &MVA_weight_btagSF_lfstats2_up, &b_MVA_weight_btagSF_lfstats2_up);
  tree->SetBranchAddress("MVA_weight_btagSF_lfstats2_down", &MVA_weight_btagSF_lfstats2_down, &b_MVA_weight_btagSF_lfstats2_down);
  
  
}
////////// FUNCTIONS
vector<double> BDTCUT(string region, string coupling){
  
  cout << "Determine BDT cut " << endl;
  
  
  string bdtinput_name = placeOutputReading+"/Reader_" + coupling + "_" + region + ".root";
  TFile* bdt_file = TFile::Open( bdtinput_name.c_str(), "READ" );
  
  cout << bdtinput_name << endl;
  
  
  
  
  string histo_name = "";
  string template_fake_name = "";
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  
  double cut_eee, cut_eeu, cut_uuu, cut_uue;
  cut_eee = cut_uue = cut_eeu = cut_uuu = -1.;
  for(int ichan=0; ichan<channel_list.size(); ichan++)
  {
    string channel = channel_list[ichan];
    TH1F *h_sum_bkg(0), *h_sum_sig(0), *h_tmp(0);
    histo_name = "";
    if(channel.find("uuu") != std::string::npos || channel.find("eeu") != std::string::npos) {template_fake_name = "FakeMu";}
    else {template_fake_name = "FakeEl";}
    
    
    for(int isample = 0; isample < datasets.size(); isample++)
    {
      string dataSetName = datasets[isample]->Name();
      //cout << dataSetName << endl;
      if(datasets[isample]->Name().find("data")!=std::string::npos) {continue; } // no signal in data
      if(datasets[isample]->Name().find("fake")==std::string::npos && datasets[isample]->Name().find("FCNC")==std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << endl;
        h_tmp = 0;
        histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
        if(h_tmp->GetEntries() != 0){
          if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
          else {h_sum_bkg->Add(h_tmp);}}
      }
      else if(datasets[isample]->Name().find("FCNC")!=std::string::npos) {
        //cout << "  -- sample " << datasets[isample]->Name() << " is signal "<< endl;
        h_tmp = 0;
        histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name();
        //histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + datasets[isample]->Name() + "_"  + systematic ;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; break;}
        h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
        if(h_tmp->GetEntries() != 0){
          if(h_sum_sig == 0) {h_sum_sig = (TH1F*) h_tmp->Clone();}
          else {h_sum_sig->Add(h_tmp);}
        }
      }
      else{
        h_tmp = 0;
        histo_name = coupling + "_" + region + "_" + channel_list[ichan] + "_" + template_fake_name;
        //cout << "  --- histo " << histo_name << endl;
        if(!bdt_file->GetListOfKeys()->Contains(histo_name.c_str())) {cout<<histo_name<<" : not found"<<endl;}
        else
        {
          h_tmp = (TH1F*) bdt_file->Get(histo_name.c_str());
          //cout << "h_tmp->GetEntries() " << h_tmp->GetEntries()  << endl;
          if(h_tmp->GetEntries() != 0){
            if(h_sum_bkg == 0) {h_sum_bkg = (TH1F*) h_tmp->Clone();}
            else {h_sum_bkg->Add(h_tmp);}
          }
        }
        
      }
    }
    //cout << "h_sum_bkg " << h_sum_bkg << " h_sum_sig "<< h_sum_sig << endl;
    if(h_sum_bkg == 0 || h_sum_sig == 0) {cout<<endl<<"--- Empty histogram (Reader empty ?) ! Exit !"<<endl<<endl; }
    //S+B histogram
    TH1F* h_total = (TH1F*) h_sum_bkg->Clone();
    h_total->Add(h_sum_sig);
    
    //Normalization
    h_total->Scale(1/h_total->Integral()); //Integral() : Return integral of bin contents in range [binx1,binx2] (inclusive !)
    h_sum_bkg->Scale(1/h_sum_bkg->Integral());
    h_sum_sig->Scale(1/h_sum_sig->Integral());
    
    double sig_over_total = 100; //initialize to unreasonable value
    int bin_cut = -1; //initialize to false value
    int nofbins = h_total->GetNbinsX();
    
    for(int ibin=nofbins; ibin>0; ibin--)
    {
      //Search the bin w/ lowest sig/total, while keeping enough bkg events (criterion needs to be optimized/tuned)
      if( (h_sum_sig->Integral(1, ibin) / h_total->Integral(1, ibin)) < sig_over_total && (h_sum_bkg->Integral(1, ibin) / h_sum_bkg->Integral()) >= 0.6 )
      {
        bin_cut = ibin;
        sig_over_total = h_sum_sig->Integral(1, bin_cut) / h_total->Integral(1,bin_cut);
      }
    }
    
    double cut = h_total->GetBinLowEdge(bin_cut+1); //Get the BDT cut value to apply to create a BDT CR control tree
    double min = TMath::Min(h_sum_bkg->GetMinimum(), h_sum_sig->GetMinimum());
    double max = TMath::Max(h_sum_bkg->GetMaximum(), h_sum_sig->GetMaximum());
    //Create plot to represent the cut on BDT
    TCanvas* c = new TCanvas("c", "Signal vs Background");
    gStyle->SetOptStat(0);
    h_sum_bkg->GetYaxis()->SetRange(min,1.049*max);
    h_sum_bkg->GetXaxis()->SetTitle("BDT Discriminant");
    h_sum_bkg->SetTitle("Signal vs Background");
    h_sum_bkg->SetLineColor(kBlue);
    h_sum_sig->SetLineColor(kRed-3);
    h_sum_bkg->SetLineWidth(3);
    h_sum_sig->SetLineWidth(3);
    //h_sum_sig->Scale(100);
    h_sum_bkg->Draw("HIST");
    h_sum_sig->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.7,0.75,0.88,0.85);
    leg->SetHeader("");
    leg->AddEntry(h_sum_sig,"Signal","L");
    leg->AddEntry(h_sum_bkg,"Background","L");
    leg->Draw();
    //Draw vertical line at cut value
    
    TLine* l = new TLine(cut,min, cut, 1.049*max);
    l->SetLineWidth(3);
    l->SetLineColor(921);
    l->Draw("");
    TString outputsaving = "CutPlots";
    mkdir(outputsaving,0777);
    outputsaving += "/"+region + "_" + coupling + "_";
    outputsaving += "Signal_Background_BDT_"+channel+".png";
    c->SaveAs(outputsaving.Data());
    
    //Cout some results
    cout<<"---------------------------------------"<<endl;
    cout<<"* Cut Value = "<<cut<<endl;
    cout<<"-> BDT_CR defined w/ all events inside bins [1 ; "<<bin_cut<<"] of the BDT distribution!"<<endl<<endl;
    cout<<"* Signal integral = "<<h_sum_sig->Integral(1, bin_cut)<<" / Total integral "<<h_total->Integral(1, bin_cut)<<endl;
    cout<<"Signal contamination in CR --> Sig/Total = "<<sig_over_total<<endl;
    cout<<"Bkg(CR) / Bkg(Total) = "<<h_sum_bkg->Integral(1,bin_cut) / h_sum_bkg->Integral()<<endl;
    cout<<"---------------------------------------"<<endl<<endl;
    
    //for(int i=0; i<h_sig->GetNbinsX(); i++) {cout<<"bin content "<<i+1<<" = "<<h_sig->GetBinContent(i+1)<<endl;} //If want to verify that the signal is computed correctly
    delete c; delete leg; delete l;
    if(channel.find("uuu") != std::string::npos) cut_uuu = cut;
    else if(channel.find("eeu") != std::string::npos) cut_eeu = cut;
    else if(channel.find("eee") != std::string::npos) cut_eee = cut;
    else if(channel.find("uue") != std::string::npos) cut_uue = cut;
    
    delete h_sum_sig;
    delete h_sum_bkg;
    delete h_tmp;
  }
  bdt_file->Close();
  vector<double> v_return;
  v_return.push_back(cut_uuu);
  v_return.push_back(cut_uue);
  v_return.push_back(cut_eeu);
  v_return.push_back(cut_eee);
  return v_return;
  
  
}
void CalculatePDFWeight(string dataSetName, double BDT, double MVA_weight_nom, int MVA_channel){
  // cout << "calculate pdf" << endl;
  //std::vector<double> pdfweights;
  //cout << "MVA channel " << MVA_channel << endl;
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  
  //PDF weights calculation
  LHAPDF::setVerbosity(0);
  string sbase = "";
  if(dataSetName.find("WZTo3LNu")!=std::string::npos) sbase = "NNPDF30_nlo_as_0118"; //kevin (ttbar)
  //cout << "base set " << sbase << endl;
  LHAPDF::PDFSet basepdfSet(sbase.c_str());
  //LHAPDF::PDFSet basepdfSet("NNPDF30_lo_as_0130"); // base from main MC // WZ
  //LHAPDF::PDFSet basepdfSet("NNPDF23_lo_as_0130_qed");
  LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_100"); // give the correct name see https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#PDF_uncertainties
  //LHAPDF::PDFSet newpdfSet("PDF4LHC15_nlo_30");
  const LHAPDF::PDF* basepdf = basepdfSet.mkPDF(0);
  channel = channel_list[MVA_channel];
  // cout << "MVA_x1 "<< MVA_x1 <<" MVA_x2 "<< MVA_x2 <<" MVA_id1 "<< MVA_id1 <<" MVA_id2 "<< MVA_id2 <<" MVA_q "<< MVA_q << endl;
  
  
  for ( size_t i=0; i<newpdfSet.size(); i++)
  {
    const LHAPDF::PDF* newpdf = newpdfSet.mkPDF(i);
    
    double weightpdf = LHAPDF::weightxxQ(MVA_id1, MVA_id2, MVA_x1, MVA_x2, MVA_q, *basepdf, *newpdf);
    output_histo_name = dataSetName+"_BDT_"+channel+"_"+intToStr(i);
    histo1D[output_histo_name]->Fill(MVA_BDT, MVA_weight_nom*weightpdf);
    // cout << "fill " << (dataSetName+"_BDT"+"_"+intToStr(i)).c_str() << " with " << BDT << " " <<  weightpdf << endl;
    //cout << "pdf weight " << weight << endl;
    //pdfweights.push_back(weight);
    
    delete newpdf;
    
  }
  output_histo_name = dataSetName+"_BDT_"+channel+"_nominal";
  histo1D[output_histo_name]->Fill(MVA_BDT, MVA_weight_nom);
  output_histo_name = "";
  delete basepdf;
  
}
void FillGeneralPlots(int d, string prefix, vector <int> decayChannels, bool isData, bool toppair, double weight_, int MVA_channel){
  //cout << "fill plots" << endl;
  decaystring = "";
  Double_t eventW = 1.;
  eventW = weight_;
  
  
  // if(datasets[d]->Name().find("fake")!=std::string::npos) eventW *= 0.0001;
  
  for(int iChan =0; iChan < decayChannels.size() ; iChan++){
    decaystring = "";
    
    if((decayChannels[iChan] != -9) && (decayChannels[iChan] != MVA_channel)) continue;
    if(decayChannels[iChan] == 0) decaystring = "uuu";
    if(decayChannels[iChan] == 1) decaystring = "uue";
    if(decayChannels[iChan] == 2) decaystring = "eeu";
    if(decayChannels[iChan] == 3) decaystring = "eee";
    if(decayChannels[iChan] == -9) decaystring = "all";
    //cout << "filling " << datasets[d]->Name() << endl;
    
    //cout << "bdt " << MVA_BDT << " in " << (prefix+"_BDT_"+decaystring).c_str()<< endl;
    MSPlot[(prefix+"_BDT_"+decaystring).c_str()]->Fill(0. , datasets[d], true, 1.);
  }
}
void GetPDFEnvelope(string dataSetName){
  double binContentMax = -1.;
  double binContentMin = 1000000000000000.;
  vector<double> bincontents;
  
  std::vector<string> channel_list;
  channel_list.push_back("eee");
  channel_list.push_back("uue");
  channel_list.push_back("eeu");
  channel_list.push_back("uuu");
  string channel = "";
  // loop over channels
  for(int iChan = 0; iChan < channel_list.size(); iChan++){
    channel = channel_list[iChan];
    output_histo_name = dataSetName+"_BDT_" + channel + "_nominal";
    // get nominal th1F
    TH1F* histo_nom = (TH1F*) histo1D[output_histo_name]->Clone();
    
    // loop over bins
    for( int ibin = 1; ibin <histo_nom->GetNbinsX(); ibin++)
    {
      binContentMax = -1.;
      binContentMin = 1000000000000000.;
      bincontents.clear();
      // get bincontents of each histo for this bin
      for(int iCount = 0; iCount < 101; iCount++)
      {
        output_histo_name = dataSetName+"_BDT_"+channel + "_" +intToStr(iCount);
        bincontents.push_back(histo1D[output_histo_name]->GetBinContent(ibin));
      }
      if(binContentMin > minimumValue(bincontents)) binContentMin = minimumValue(bincontents);
      else binContentMin = histo_nom->GetBinContent(ibin);
      if(binContentMax < maximumValue(bincontents)) binContentMax = maximumValue(bincontents);
      else binContentMax = histo_nom->GetBinContent(ibin);
      
      output_histo_name = dataSetName+"_BDT_"+channel + "_PDFEnvelopeUp";
      histo1D[output_histo_name]->SetBinContent(ibin, binContentMax);
      output_histo_name = dataSetName+"_BDT_"+channel + "_PDFEnvelopeDown";
      histo1D[output_histo_name]->SetBinContent(ibin, binContentMin);;
      output_histo_name = "";
    }// bins
  }//channels
}
