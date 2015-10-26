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
  string xmlFileName ="config/test.xml";
  
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
  //vector < Dataset* > datasetsMu;
  //vector < Dataset* > datasetsEl;
  
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets(datasets, xmlfile); cout << "Number of datasets: " << datasets.size() << endl;
  for (unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
  
  float Luminosity = oldLuminosity;
  //float LuminosityMu = oldLuminosity;
  //float LuminosityEl = oldLuminosity;
  
  //bool foundMu = false;
  //bool foundEl = false;
  
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
    if (Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
    string dataSetName = datasets[d]->Name();
    
    //if (dataSetName.find("Data_Mu") == 0 || dataSetName.find("data_Mu") == 0 || dataSetName.find("DATA_Mu") == 0) {
    //  LuminosityMu = datasets[d]->EquivalentLumi();
    //  foundMu=true;
    //}
    //if (dataSetName.find("Data_El") == 0 || dataSetName.find("data_El") == 0 || dataSetName.find("DATA_El") == 0) {
    //  LuminosityEl = datasets[d]->EquivalentLumi();
    //  foundEl=true;
    //}
    
   // if ( dataSetName.find("QCD") == 0 ) datasets[d]->SetColor(kYellow);
   // if ( dataSetName.find("TT") == 0 ) datasets[d]->SetColor(kRed+1);
   //if ( dataSetName.find("TTbarJets_Other") == 0 ) datasets[d]->SetColor(kRed-7);
   // if ( dataSetName.find("WJets") == 0 )
   /* {
      datasets[d]->SetTitle("W#rightarrowl#nu");
      datasets[d]->SetColor(kGreen-3);
    }
    if ( dataSetName.find("ZJets") == 0 )
    {
      datasets[d]->SetTitle("Z/#gamma*#rightarrowl^{+}l^{-}");
      datasets[d]->SetColor(kMagenta);
    }
    if ( dataSetName.find("ST") == 0 || dataSetName.find("SingleTop") ==0 )
      datasets[d]->SetColor(kBlue-2);
    //if (dataSetName.find("NP") == 0 )
    //{
    //	datasets[d]->SetTitle("Signal");
    //	datasets[d]->SetColor(kGreen+4);
    //}
  }
  */
  if( dataSetName.find("WZ") == 0 )
  {
      datasets[d]->SetColor(kBlue); 
  }
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  
  cout << " - Recreate output file ..." << endl;
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  
  //Global variable
  //TRootEvent* event = 0;
  
  //nof selected events
  //double NEvtsData = 0;
  Double_t *nEvents = new Double_t[datasets.size()];
  
  
  
  ////////////////////////////////////
  ///  Selection table
  ////////////////////////////////////
  
  vector<string> CutsSelecTableSemiMu;
  //CutsSelecTableSemiMu.push_back(string("preselected"));
  //CutsSelecTableSemiMu.push_back(string("trigged"));
  //CutsSelecTableSemiMu.push_back(string("Good PV"));
  //CutsSelecTableSemiMu.push_back(string("1 selected muon"));
  //CutsSelecTableSemiMu.push_back(string("Veto 2nd muon"));
  //CutsSelecTableSemiMu.push_back(string("Veto electron"));
  
  //char LabelNJets[100];
  //sprintf(LabelNJets,"$\\geq$ %d jets", anaEnv.NofJets);
  //CutsSelecTableSemiMu.push_back(string(LabelNJets));
  
  //CutsSelecTableSemiMu.push_back("Missing $E_T$");
  //CutsSelecTableSemiMu.push_back("$H_T$ cut");
  //CutsSelecTableSemiMu.push_back("$\\geq$ 1 b-jet (CSVM)");
  //CutsSelecTableSemiMu.push_back("$\\geq$ 2 b-jets (CSVM)");
  
  if (verbose > 0)
    cout << " - CutsSelectionTable instantiated ..." << endl;
  SelectionTable selecTableSemiMu(CutsSelecTableSemiMu, datasets);
  //selecTableSemiMu.SetLuminosity(LuminosityMu);
  selecTableSemiMu.SetLuminosity(Luminosity);
  if (verbose > 0)
    cout << " - SelectionTable instantiated ..." << endl;
  
  
  
  ////////////////////////////////////
  ///  Loop on datasets
  ////////////////////////////////////
  
  if (verbose > 0)
    cout << " - Loop over datasets ... " << datasets.size() << " datasets !" << endl;
/*  
  for (unsigned int d = 0; d < datasets.size(); d++)
  { 
    cout << "equivalent luminosity of dataset " << datasets[d]->EquivalentLumi() << endl;
    string previousFilename = "";
    int iFile = -1;
    string dataSetName = datasets[d]->Name();
    nofSelectedEvents = 0;
    
    if (verbose > 1)
      cout << "   Dataset " << d << ": " << datasets[d]->Name() << "/ title : " << datasets[d]->Title() << endl;
    if (verbose > 1)
      cout << "      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
    
    //open files and load
    cout << "LoadEvent" << endl;
    treeLoader.LoadDataset(datasets[d], anaEnv);
    cout << "LoadEvent" << endl;
    
    
    
    ////////////////////////////////////
    ///  Loop on events
    ////////////////////////////////////
    
    nEvents[d] = 0;
    int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
    if (verbose > 1)
      cout << "	Loop over events " << endl;
    
    //for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
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
      
      //if (! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
      //  genjets = treeLoader.LoadGenJet(ievt,false);
      //  //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
      //}
      
      
      // BE CAREFUL: TRootGenEvent is now obsolete!
      
      
      
      /////////////////////////////////////
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
      
      
      
      /////////////////////////
      ///  EVENT SELECTION  ///
      /////////////////////////
      
      //Declare selection instance
      Run2Selection selection(init_jets_corrected, init_muons, init_electrons, mets);
      
      bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2.);
      
      vector<TRootPFJet*> selectedJets = selection.GetSelectedJets(20, 2.5, true, "Tight");  // GetSelectedJets(float PtThr, float EtaThr, bool applyJetID, std::string TightLoose)
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(20, 2.5, 0.12);  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
      //vector<TRootElectron*> selectedElectrons = selection.GetSelectedTightElectronsCutsBasedSpring15_50ns(20, 2.5);  // GetSelectedTightElectronsCutsBasedSpring15_50ns(float PtThr, float EtaThr)
      vector<TRootMuon*> vetoMuons = selection.GetSelectedLooseMuons(10, 2.5, 0.2);  // GetSelectedLooseMuons(float PtThr, float EtaThr,float MuonRelIso)
      //vector<TRootElectron*> vetoElectronsSemiMu = selection.GetSelectedLooseElectronsCutsBasedSpring15_50ns(20, 2.5);  // GetSelectedLooseElectronsCutsBasedSpring15_50ns(float PtThr, float EtaThr)
      
      
      //if (selectedJets.size() >= 4)
      //  if (selectedJets[3]->Pt() < 30) selectedJets.clear();
      
      vector<TRootMCParticle*> mcParticles;
      
      if ( dataSetName.find("TT") == 0 )
      {
        treeLoader.LoadMCEvent(ievt, 0, 0, mcParticles,false);
        sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
      }
      
      eventSelected = false;
      
      selecTableSemiMu.Fill(d,0,scaleFactor);
      /// At the moment do not use trigger
      selecTableSemiMu.Fill(d,1,scaleFactor);
      if (isGoodPV)
      {
        selecTableSemiMu.Fill(d,2,scaleFactor);
        if (selectedMuons.size() == 1)
        {
          selecTableSemiMu.Fill(d,3,scaleFactor);
          //if (vetoMuons.size() == 1) {
ent at 14:14 on Wednesday
 Isabelle:  hij kan niet en Joyce werkt
 Sent at 14:17 on Wednesday
 lieselotte moreels has left.
 me:  ok
 
          //  if (vetoElectronsSemiMu.size() == 0) {
          //    selecTableSemiMu.Fill(d,5,scaleFactor);
              if ( selectedJets.size() >= 4 )
              {
                //selecTableSemiMu.Fill(d,6,scaleFactor);
                selecTableSemiMu.Fill(d,4,scaleFactor);
                eventSelected = true;
              }  // at least 4 jets
          //  }  // no loose electrons
          //}  // no loose muons
        }  // 1 good muon
      }  // good PV
      
      
      if (! eventSelected )
      {
        //cout << "Event no. " << ievt << " was not selected. " << endl;
        continue;
      }
      
      nofSelectedEvents++;
      
      if (verbose > 2)
        cout << "  Event " << ievt << " is selected" << endl;
      
      
      
      /// Save selected events in new root file
      //......
      //......
      
      
      
      
      
      
      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  //loop on events
  */  
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // Loop on datasets
  
  
  
  ///Selection tables
  //selecTableSemiMu.TableCalculator(false, true, true, true, true);
  //string selectiontableMu = "SelectionTable_test_SemiMu.tex";
  //selecTableSemiMu.Write(selectiontableMu.c_str());
  //selecTableSemiMu.Write(selectiontableMu.c_str(), true, true, true, true, true, true, false);
  
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
