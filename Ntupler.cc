////////////////////////////////////////////////////////////
/////                                                  /////
/////  Preliminary macro to try and make some plots    /////
/////                                                  /////
////////////////////////////////////////////////////////////



#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <TLorentzVector.h>

//used TopTreeAnalysis classes
#include "../TopTreeProducer/interface/TRootRun.h"
#include "../TopTreeProducer/interface/TRootEvent.h"
#include "../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
#include "../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"


using namespace std;
using namespace reweight;
using namespace TopTree;


int main (int argc, char *argv[])
{
  
  string rootFileName = "testAnalyser_output.root";
  
  clock_t start = clock();
  
  
  /////////////////////
  ///  Configuration
  /////////////////////
  
  bool eventSelected = false;
  int nofSelectedEvents = 0;
  
  
  /// xml file
  string xmlFileName ="config/Run2TriLepton_samples.xml";
  
  if (argc > 1)
    xmlFileName = (string)argv[1];
  
  const char *xmlfile = xmlFileName.c_str();
  
  cout << " - Using config file " << xmlfile << endl;
  
  //Configuration output format
  TTree *configTree = new TTree("configTree","configuration Tree");
  TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
  configTree->Branch("Datasets","TClonesArray",&tcdatasets);
  TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
  configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
  
  
  
  ////////////////////////////////////
  ///  AnalysisEnvironment
  ////////////////////////////////////
  
  AnalysisEnvironment anaEnv;
  cout << " - Loading environment ..." << endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  
  cout << anaEnv.JetCollection << " " <<  anaEnv.METCollection << " "
    << anaEnv.ElectronCollection << " " << anaEnv.MuonCollection << " "
    << anaEnv.PrimaryVertexCollection << " " << anaEnv.GenJetCollection << " "
    << endl;
  
  new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
  int verbose = anaEnv.Verbose;
  float oldLuminosity = anaEnv.Luminosity;  // in 1/pb
  cout << "Analysis environment luminosity for rescaling " << oldLuminosity << endl;
  
  
  
  
  
  
  
  /////////////////////
  ///  Load Datasets
  /////////////////////
  
  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;

  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); 
  cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;


  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    

  if( dataSetName.find("WZ") == 0 )
  {
      datasets[d]->SetColor(kBlue); 
  }
  }
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  cout << " - Recreate output file ..." << endl;
 // if(dataSetName == "WZ") { sprintf(rootFileName, "wz.root"); 
  
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  
  //Global variable
  TRootEvent* event = 0;
  
  //nof selected events
  double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  
  
  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTable;
  CutsSelecTable.push_back(string("preselected"));
  CutsSelecTable.push_back(string("trigged"));
  CutsSelecTable.push_back(string("Good PV"));
  CutsSelecTable.push_back(string("3 selected leptons"));
   
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
    
  SelectionTable selecTable(CutsSelecTable, datasets);
  selecTable.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;

  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    
    int isdata; 
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << " - title : " << datasets[d]->Title() << endl;
    if (verbose > 1)
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    
    // make root tree file name
    string roottreename = datasets[d]->Name();
    roottreename+="_tree.root";
    
    // create the output file that is used for further analysis. This file can contain histograms and/or trees (in this case only a tree)
    TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
    fileout->cd();
    
     //////////////////////////////
     // My tree - variables //
     //////////////////////////////
     // various weights
     Double_t pu_weight;
     Int_t run_num;
     Int_t evt_num;
     Int_t lumi_num;
     Int_t nvtx;
     Int_t npu;
     
     //other
     Int_t nElectrons;
    
    
     ///////////////////////////////
     // My trees                  //
     ///////////////////////////////
     TTree *bookkeeping = new TTree("startevents","startevents");
     bookkeeping->Branch("run_num",&run_num,"run_num/I");
     bookkeeping->Branch("evt_num",&evt_num,"evt_num/I");
     bookkeeping->Branch("lumi_num",&lumi_num,"lumi_num/I");
     bookkeeping->Branch("nvtx",&nvtx,"nvtx/I");
     bookkeeping->Branch("npu",&npu,"npu/I");


     // define the output tree
     TTree* myTree = new TTree("tree","tree");
     myTree->Branch("isdata",&isdata,"isdata/I");
     myTree->Branch("run_num",&run_num,"run_num/I");
     myTree->Branch("evt_num",&evt_num,"evt_num/I");
     myTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
     myTree->Branch("nvtx",&nvtx,"nvtx/I");
     myTree->Branch("npu",&npu,"npu/I");
     
     myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
    
    
     TTree *noselTree = new TTree("noselTree","noselTree");
     noselTree->Branch("isdata",&isdata,"isdata/I");
     noselTree->Branch("run_num",&run_num,"run_num/I");
     noselTree->Branch("evt_num",&evt_num,"evt_num/I");
     noselTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
     noselTree->Branch("nvtx",&nvtx,"nvtx/I");
     noselTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
    
    
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    
     nofSelectedEvents = 0; 
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    //some bookkeeping variables
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
   // for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
    for (unsigned int ievt = 0; ievt < 2000; ievt++)
    {
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      //vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%1000 == 0)
        std::cout << "Processing the " << ievt << "th event (" << ((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100  << "%)" << flush << "\r";
      
      
      
      ////////////////////
      ///  LOAD EVENT  ///
      ////////////////////
      
      TRootEvent* event = treeLoader.LoadEvent(ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);


      ////////////////////////////
      ///  Include trigger set up here when using data
      ////////////////////////////
      
      string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
      if (previousFilename != currentFilename)
      {
        previousFilename = currentFilename;
        iFile++;
        cout << "File changed!!! => iFile = " << iFile << endl;
      }
      
      int currentRun = event->runId();
      
      if (previousRun != currentRun)
        previousRun = currentRun;
      
      
      run_num=event->runId();
      evt_num=event->eventId();
      lumi_num=event->lumiBlockId();
      nvtx=vertex.size();
      npu=(int)event->nTruePU();
      if( run_num > 10000){//data
         isdata=1;
      }
      bookkeeping->Fill(); 
      
      
      ////////////////////////////////////
      ///  DETERMINE EVENT SCALEFACTOR  ///
      /////////////////////////////////////
      
      // scale factor for the event
      float scaleFactor = 1.;
      
      
      // PU reweighting
      
      // old method
      //cout << "scalefactor " << scaleFactor << endl;
      double lumiWeight = 1; //LumiWeights.ITweight( (int)event->nTruePU() ); // currently no pile-up reweighting applied
      
      if (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
        lumiWeight=1;
      
      // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
      // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
      
      scaleFactor = scaleFactor*lumiWeight;
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.); //isPVSelected(const std::vector<TRootVertex*>& vertex, int NdofCut, float Zcut, float RhoCut)
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(20, 2.5, true, "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(20, 2.5, 0.12);  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(20, 2.5, "Tight", "Spring15_50ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
   
      
      
      eventSelected = false;
      
      selecTable.Fill(d,0,scaleFactor);
      /// At the moment do not use trigger
      selecTable.Fill(d,1,scaleFactor);
      if (isGoodPV)
      {
     
        selecTable.Fill(d,2,scaleFactor);
	nElectrons = selectedElectrons.size(); 
	
	noselTree->Fill();
        if (selectedMuons.size() + selectedElectrons.size()== 3)
        {
          selecTable.Fill(d,3,scaleFactor);
	  
         
        }  // 3 leptons
      }  // good PV
      
      
      if (! eventSelected )
      {
        //cout << "Event no. " << ievt << " was not selected. " << endl;
        continue;
      }
      
      nofSelectedEvents++;
      myTree->Fill();
      
      if (verbose > 2)
        cout << "  Event " << ievt << " is selected" << endl;
      

      
      
      
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  //loop on events
   
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    
    
    myTree->Write();
    fileout->Write();
    fileout->Close();
    // delete myTree;
    delete fileout;
    
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // Loop on datasets
  
  

  
  fout->Close();
  
  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
  
}
