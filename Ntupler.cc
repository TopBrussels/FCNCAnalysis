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

std::vector <int> OSSFLeptonPairCalculator(std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 
TLorentzVector CreateZboson(std::vector<int> Lep, std::vector<TRootElectron*> Elec, std::vector<TRootMuon*> Mu, int verb); 

int main (int argc, char *argv[])
{
  
//  string rootFileName = "testAnalyser_output.root";
  
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
    

    if( dataSetName.find("WZ") == 0 ){      datasets[d]->SetColor(kBlue);   }
  }
  if ( Luminosity != oldLuminosity ) cout << "Changed analysis environment luminosity to "<< Luminosity << endl;
  
  


  
  
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
  CutsSelecTable.push_back(string("At least 2 jets")); 
  CutsSelecTable.push_back(string("At least 1 CSV loose jet")); 
  CutsSelecTable.push_back(string("At least 1 OSSF lepton pair"));
  CutsSelecTable.push_back(string("Z mass window of 15 GeV")); 
   
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
    string roottreename = "Ntuples/";
    roottreename+= datasets[d]->Name();
    roottreename+="_tree.root";

    cout << "  - Recreate outputfile ... " << roottreename.c_str() << endl; 
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
    
     //doubles
     double ptZboson;
     double pxZboson;
     double pyZboson;
     double pzZboson;
     double etaZboson;
     double eZboson; 
     double mZboson; 

     double ptWboson_lep;
     double pxWboson_lep;
     double pyWboson_lep;
     double pzWboson_lep;
     double etaWboson_lep;
     double eWboson_lep; 
 
     //vectors
     std::vector<double> *ptMuon; 
     std::vector<double> *pxMuon; 
     std::vector<double> *pyMuon; 
     std::vector<double> *pzMuon; 
     std::vector<double> *etaMuon;
     std::vector<double> *eMuon; 
     std::vector<double> *qMuon;
     
     std::vector<double> *ptElectron; 
     std::vector<double> *pxElectron; 
     std::vector<double> *pyElectron; 
     std::vector<double> *pzElectron; 
     std::vector<double> *etaElectron; 
     std::vector<double> *eElectron; 
     std::vector<double> *qElectron;
     
     std::vector<double> *ptJet; 
     std::vector<double> *pxJet; 
     std::vector<double> *pyJet; 
     std::vector<double> *pzJet; 
     std::vector<double> *eJet; 
     std::vector<double> *etaJet;
     std::vector<double> *qJet;  
     std::vector<double> *BtagCSVjet;
     std::vector<bool> *BtagCSVL; 

    
    
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
     
     //Set branches for doubles 
 ///   myTree -> Branch("metPt", &metPt, "metPt/D");
 //    myTree -> Branch("metPx", &metPx, "metPx/D");
 //    myTree -> Branch("metPy", &metPy, "metPy/D");
     
     myTree->Branch("ptZboson", &ptZboson,"ptZboson/D");
     myTree->Branch("pxZboson", &pxZboson,"pxZboson/D");
     myTree->Branch("pyZboson", &pyZboson,"pyZboson/D");
     myTree->Branch("pzZboson", &pzZboson,"pzZboson/D");  
     myTree->Branch("etaZboson", &etaZboson,"etaZboson/D");
     myTree->Branch("eZboson", &eZboson,"eZboson/D");
     myTree->Branch("mZboson", &mZboson,"mZboson/D");
     
     myTree->Branch("ptWboson_lep", &ptWboson_lep,"ptWboson_lep/D");
     myTree->Branch("pxWboson_lep", &pxWboson_lep,"pxWboson_lep/D");
     myTree->Branch("pyWboson_lep", &pyWboson_lep,"pyWboson_lep/D");
     myTree->Branch("pzWboson_lep", &pzWboson_lep,"pzWboson_lep/D");
     myTree->Branch("etaWboson_lep", &etaWboson_lep,"etaWboson_lep/D");
     myTree->Branch("eWboson_lep", &eWboson_lep,"eWboson_lep/D");
    
     // Set branches for vectors
     myTree->Branch("ptMuon","std::vector<double>",&ptMuon);   // Make a branch with name ptMuon, from type vector(double), on loaction ptMuon
     myTree->Branch("pxMuon","std::vector<double>",&pxMuon);
     myTree->Branch("pyMuon","std::vector<double>",&pyMuon);
     myTree->Branch("pzMuon","std::vector<double>",&pzMuon);
     myTree->Branch("eMuon","std::vector<double>",&eMuon);
     myTree->Branch("etaMuon","std::vector<double>",&etaMuon);
     myTree->Branch("qMuon","std::vector<double>",&qMuon);
     
     myTree->Branch("ptElectron","std::vector<double>",&ptElectron);   // Make a branch with name ptElectron, from type vector(double), on loaction ptElectron
     myTree->Branch("pxElectron","std::vector<double>",&pxElectron);
     myTree->Branch("pyElectron","std::vector<double>",&pyElectron);
     myTree->Branch("pzElectron","std::vector<double>",&pzElectron);
     myTree->Branch("eElectron","std::vector<double>",&eElectron);
     myTree->Branch("etaElectron","std::vector<double>",&etaElectron);
     myTree->Branch("qElectron","std::vector<double>",&qElectron);     
     
     myTree->Branch("ptJet","std::vector<double>",&ptJet);
     myTree->Branch("pxJet","std::vector<double>",&pxJet);
     myTree->Branch("pyJet","std::vector<double>",&pyJet);
     myTree->Branch("pzJet","std::vector<double>",&pzJet);
     myTree->Branch("eJet","std::vector<double>",&eJet);
     myTree->Branch("etaJet","std::vector<double>",&etaJet);
     myTree->Branch("qJet","std::vector<double>",&qJet);
     myTree->Branch("BtagCSVjet", "std::vector<double>",&BtagCSVjet);
     myTree->Branch("BtagCSVL","std::vector<bool>",&BtagCSVL);
    
    
    
//     TTree *noselTree = new TTree("noselTree","noselTree");
/*      noselTree->Branch("isdata",&isdata,"isdata/I");
     noselTree->Branch("run_num",&run_num,"run_num/I");
     noselTree->Branch("evt_num",&evt_num,"evt_num/I");
     noselTree->Branch("lumi_num",&lumi_num,"lumi_num/I");
     noselTree->Branch("nvtx",&nvtx,"nvtx/I");
 */    // noselTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
    
    
    
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
    
    for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++)
//    for (unsigned int ievt = 0; ievt < 2000; ievt++)
    {
      vector < TRootVertex* > vertex;
      vector < TRootMuon* > init_muons;
      vector < TRootElectron* > init_electrons;
      vector < TRootJet* > init_jets_corrected;
      vector < TRootJet* > init_jets;
      vector < TRootMET* > mets;
      //vector < TRootGenJet* > genjets;
      
      nEvents[d]++;
      
      if (ievt%500 == 0)
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
      vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(20, 2.5, 0.12,"Tight","Spring15");  // GetSelectedMuons(float PtThr, float EtaThr,float MuonRelIso)
      vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(20, 2.5, "Tight", "Spring15_50ns", true);  // GetSelectedElectrons(float PtThr, float etaThr, string WorkingPoint, string ProductionCampaign, bool CutsBased)
      
      sort(selectedJets.begin(), selectedJets.end(),HighestPt()); 
      sort(selectedMuons.begin(), selectedMuons.end(), HighestPt()); 
      sort(selectedElectrons.begin(), selectedElectrons.end(), HighestPt()); 
      
      vector<bool> BtagBooleans; 
      BtagBooleans.clear(); 
      
      // Start analysis selection
      eventSelected = false;
      TLorentzVector Zboson;
      TLorentzVector Wlep; 
       selecTable.Fill(d,0,scaleFactor);
      /// At the moment do not use trigger
      selecTable.Fill(d,1,scaleFactor);
      if (isGoodPV)
      {
        if(verbose>3) cout << "GoodPV" << endl; 
        selecTable.Fill(d,2,scaleFactor);

        if (selectedMuons.size() + selectedElectrons.size()== 3)
        {
            if(verbose>3) cout << "3 electrons "<< endl; 
	    selecTable.Fill(d,3,scaleFactor);
	    
	    
	    if(selectedJets.size() > 1)
	    {
	       if(verbose>3) cout << " at least 2 jets " << endl; 
	       selecTable.Fill(d,4,scaleFactor); 
	       int nBtagged = 0; 
	       for(unsigned int i = 0; i < selectedJets.size() ; i++)
	       {
	         
	         bool Btagged = false; 
		 TRootJet* tempJet = (TRootJet*) selectedJets[i];
		 if(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.605)//loose WP
		 {
		    Btagged = true; 
		    nBtagged++; 
	         }
		 BtagBooleans.push_back(Btagged);
	       } 
	       if(verbose > 3) cout << "btagging done" << endl; 
	       if(nBtagged>1)
	       {
	         if(verbose>3) cout << " at least 1 bjet " << endl; 
	         selecTable.Fill(d,5,scaleFactor); 
		 
		 std::vector<int> Leptons; 
		 Leptons.clear(); 
		 Leptons = OSSFLeptonPairCalculator(selectedElectrons, selectedMuons, verbose);
		 if(verbose>3) cout <<   Leptons[0]<< " , " <<  Leptons[1]<< " , " <<  Leptons[2]<< " , " <<  Leptons[3]<< " , " <<  Leptons[4]<< " , " <<  Leptons[5]   << endl; 
		 
		 bool OSSFpair = false; 
		 if( (Leptons[0] != -5 && Leptons[1] != -5) | (Leptons[3] != -5 && Leptons[4] != -5) ) OSSFpair = true; 
		 if(OSSFpair)
	         {
		    if(verbose>3) cout << " OSSF "<< endl; 
		    selecTable.Fill(d,6,scaleFactor); 
		    
		    Zboson.Clear(); 
		    Zboson = CreateZboson(Leptons, selectedElectrons, selectedMuons, verbose); 
		    
                    Wlep.Clear(); 
                    if(fabs(Leptons[2]) != 5)
		    {
		      if(verbose>3) cout << " the W lepton is an electron " << endl; 
		      Wlep.SetPxPyPzE(selectedElectrons[Leptons[2]]->Px(), selectedElectrons[Leptons[2]]->Py(), selectedElectrons[Leptons[2]]->Pz(), selectedElectrons[Leptons[2]]->Energy()); 
                    }
		    else if(fabs(Leptons[5]) != 5)
		    {
		      if(verbose>3) cout << " the W lepton is a muon " << endl; 
		      Wlep.SetPxPyPzE(selectedMuons[Leptons[5]]->Px(), selectedMuons[Leptons[5]]->Py(), selectedMuons[Leptons[5]]->Pz(), selectedMuons[Leptons[5]]->Energy()); 		    
		    }
		    bool ZmassWindow = false; 
		    if(fabs(Zboson.M()-90.0) < 15.0) ZmassWindow = true; 
		    if(ZmassWindow)
		    {
		      if(verbose>3) cout << " Zmass window " << endl; 
		      selecTable.Fill(d,7,scaleFactor);
		      eventSelected = true;
	            }
		 }
	       } 
	    }
         
        }  // 3 leptons
      }  // good PV
      
      
      if (! eventSelected )
      {
        continue;
      }
      
      if(verbose>3) cout << "filling the tree" << endl; 
      ptMuon = new std::vector<double>; 
      pxMuon = new std::vector<double>; 
      pyMuon = new std::vector<double>; 
      pzMuon = new std::vector<double>; 
      etaMuon = new std::vector<double>;
      eMuon = new std::vector<double>; 
      qMuon = new std::vector<double>;
      
      ptElectron = new std::vector<double>; 
      pxElectron = new std::vector<double>; 
      pyElectron = new std::vector<double>; 
      pzElectron = new std::vector<double>; 
      etaElectron = new std::vector<double>; 
      eElectron = new std::vector<double>; 
      qElectron = new std::vector<double>;
      
 
      ptJet = new std::vector<double>; 
      pxJet = new std::vector<double>; 
      pyJet = new std::vector<double>; 
      pzJet = new std::vector<double>; 
      eJet = new std::vector<double>; 
      etaJet = new std::vector<double>; 
      qJet = new std::vector<double>; 
      BtagCSVjet = new std::vector<double>;
      BtagCSVL = new std::vector<bool>; 
      
      
      ptZboson = Zboson.Pt();
      pxZboson = Zboson.Px();
      pyZboson = Zboson.Py();
      pzZboson = Zboson.Pt();
      etaZboson = Zboson.Eta();
      eZboson = Zboson.Energy();
      mZboson = Zboson.M();  

      ptWboson_lep = Wlep.Pt();
      pxWboson_lep = Wlep.Px();
      pyWboson_lep = Wlep.Py();
      pzWboson_lep = Wlep.Pz();
      etaWboson_lep = Wlep.Eta();
      eWboson_lep = Wlep.Energy();
      
      for (unsigned int i = 0; i < selectedElectrons.size(); i++) 
      {
      	TRootElectron* tempElectron = (TRootElectron*) selectedElectrons[i];
        ptElectron->push_back(tempElectron->Pt()); 
	pxElectron->push_back(tempElectron->Px()); 
	pyElectron->push_back(tempElectron->Py()); 
	pzElectron->push_back(tempElectron->Pz()); 
	eElectron->push_back(tempElectron->Energy());
	etaElectron->push_back(tempElectron->Eta()); 
	qElectron->push_back(tempElectron->charge());
      }
      
      for (unsigned int i = 0; i < selectedMuons.size(); i++) 
      {
      	TRootMuon* tempMuon = (TRootMuon*) selectedMuons[i];
        ptMuon->push_back(tempMuon->Pt()); 
	pxMuon->push_back(tempMuon->Px()); 
	pyMuon->push_back(tempMuon->Py()); 
	pzMuon->push_back(tempMuon->Pz()); 
	eMuon->push_back(tempMuon->Energy());
	etaMuon->push_back(tempMuon->Eta()); 
	qMuon->push_back(tempMuon->charge());
      }
            
      for (unsigned int i =0; i < selectedJets.size(); i ++)
      {
        TRootJet* tempJet = (TRootJet*) selectedJets[i];
        ptJet->push_back(tempJet->Pt());
        pxJet->push_back(tempJet->Px());
        pyJet->push_back(tempJet->Py());
        pzJet->push_back(tempJet->Pz());
        eJet->push_back(tempJet->Energy());
        etaJet->push_back(tempJet->Eta());
        qJet->push_back(tempJet->charge());
        BtagCSVjet->push_back(tempJet->btag_combinedInclusiveSecondaryVertexV2BJetTags()); 
        BtagCSVL->push_back(BtagBooleans[i]); 
      }
 
      nofSelectedEvents++; 
      myTree->Fill();
      
      delete pxElectron;
      delete pyElectron;
      delete pzElectron;
      delete etaElectron; 
      delete eElectron;
      delete qElectron;
      
      delete ptMuon;
      delete pxMuon;
      delete pyMuon;
      delete pzMuon;
      delete etaMuon; 
      delete eMuon;
      delete qMuon;
            
      delete ptJet;
      delete pxJet;
      delete pyJet;
      delete pzJet;
      delete eJet;
      delete etaJet;
      delete qJet;
      delete BtagCSVjet;
      delete BtagCSVL; 
      
      
      if (verbose > 2)
        cout << "  Event " << ievt << " is selected" << endl;
      

      //////////////////////
      ///  END OF EVENT  ///
      //////////////////////
      
    }  //loop on events
   
    cout << endl;
    cout << "Data set " << datasets[d]->Title() << " has " << nofSelectedEvents << " selected events." << endl;
    
    selecTable.TableCalculator(false, true, true, true, true);
    string selectionTable = "SelectionTables/SelectionTable_"+datasets[d]->Name() +".tex";
  //selecTableSemiMu.Write(selectiontableMu.c_str());
    selecTable.Write(selectionTable.c_str(), true, true, true, true, true, true, false);
    
    myTree->Write();
    fileout->Write();
    fileout->Close();
    delete fileout; 
 //   delete myTree;
   
    
    
    
    //important: free memory
    treeLoader.UnLoadDataset();
    
  }  // Loop on datasets
  
  
   selecTable.TableCalculator(false, true, true, true, true); 
   string selectionTableAll = "SelectionTables/SelectionTable_allSamples.tex";
   selecTable.Write(selectionTableAll.c_str(), true, true, true, true, true, true, false);
  
//  fout->Close();
  
//  delete fout;
  delete tcdatasets;
  delete tcAnaEnv;
  delete configTree;
 
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " s to run the program" << endl;
    
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
  
}


std::vector<int> OSSFLeptonPairCalculator(vector<TRootElectron*> electrons, vector<TRootMuon*> muons, int verbose)
{
  if(verbose > 3) cout << "In OSSFLeptonPairCalculator " << endl; 
  std::vector<int> leptons;
  leptons.clear();  
  int ElecZ0 = -5; 
  int ElecZ1 = -5; 
  int ElecW = -5; 
  int MuZ0 = -5; 
  int MuZ1 = -5; 
  int MuW = -5;
  
  if(electrons.size() == 3)
  {
    if(electrons[0]->charge() != electrons[1]->charge())
    {
      ElecZ0 = 0; 
      ElecZ1 = 1;
      ElecW = 2; 
      if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
    }    
    else if (electrons[0]->charge() != electrons[2]->charge())
    {
      ElecZ0 = 0; 
      ElecZ1 = 2; 
      ElecW = 1; 
      if(verbose>3) cout << " the Zboson consistes of electrons " << endl; 
    }
    else if (electrons[1]->charge() != electrons[2]->charge())
    {
      ElecZ0 = 1; 
      ElecZ1 = 2; 
      ElecW = 0;
      if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
    }

  }
  else if ( (electrons.size() == 2) && (electrons[0]->charge() != electrons[1]->charge()) ) 
  {
    ElecZ0 = 0; 
    ElecZ1 = 1; 
    MuW = 0; 
    if(verbose>3) cout << " the Zboson consists of electrons " << endl; 
  }
  else if ( (muons.size() == 2) && (muons[0]->charge() != muons[1]->charge()) )
  {
    ElecW = 0; 
    MuZ0 = 0; 
    MuZ1 = 1;
    if(verbose>3) cout << " the Zboson consists of muons " << endl; 
  }
  else if (muons.size() == 3)
  {
    if(muons[0]->charge() != muons[1]->charge())
    {
      MuZ0 = 0; 
      MuZ1 = 1;
      MuW = 2; 
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }    
    else if (muons[0]->charge() != muons[2]->charge())
    {
      MuZ0 = 0; 
      MuZ1 = 2; 
      MuW = 1; 
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }
    else if (muons[1]->charge() != muons[2]->charge())
    {
      MuZ0 = 1; 
      MuZ1 = 2; 
      MuW = 0;
      if(verbose>3) cout << " the Zboson consists of muons " << endl; 
    }
     
  }
  leptons.push_back(ElecZ0); 
  leptons.push_back(ElecZ1); 
  leptons.push_back(ElecW); 
  leptons.push_back(MuZ0); 
  leptons.push_back(MuZ1); 
  leptons.push_back(MuW); 
  if(verbose>3) cout << " out OSSF.. " << endl; 
  if(verbose>3) cout <<   leptons[0]<< " , " <<  leptons[1]<< " , " <<  leptons[2]<< " , " <<  leptons[3]<< " , " <<  leptons[4]<< " , " <<  leptons[5]   << endl; 
  
  return leptons; 
}


TLorentzVector CreateZboson(std::vector<int> leptons, std::vector<TRootElectron*> electrons, std::vector<TRootMuon*> muons, int verbose)
{
  if(verbose>3) cout << " in Zboson creator " << endl; 
  TLorentzVector Zbos;
  Zbos.Clear(); 
  TLorentzVector Zbos_lep0; 
  Zbos_lep0.Clear(); 
  TLorentzVector Zbos_lep1; 
  Zbos_lep1.Clear(); 
  
  if(verbose>3) cout <<   leptons[0]<< " , " <<  leptons[1]<< " , " <<  leptons[2]<< " , " <<  leptons[3]<< " , " <<  leptons[4]<< " , " <<  leptons[5]   << endl; 
  
  if(fabs(leptons[0]) < 3 && fabs(leptons[1]) < 3) 
  {
      Zbos_lep0.SetPxPyPzE(electrons[leptons[0]]->Px(), electrons[leptons[0]]->Py(), electrons[leptons[0]]->Pz(), electrons[leptons[0]]->Energy()); 
      Zbos_lep1.SetPxPyPzE(electrons[leptons[1]]->Px(), electrons[leptons[1]]->Py(),electrons[leptons[1]]->Pz(), electrons[leptons[1]]->Energy()); 
  }
  else if(fabs(leptons[3]) < 3 && fabs(leptons[4]) < 3)
  {
     Zbos_lep0.SetPxPyPzE(muons[leptons[3]]->Px(),muons[leptons[3]]->Py(), muons[leptons[3]]->Pz(),muons[leptons[3]]->Energy()); 
     Zbos_lep1.SetPxPyPzE(muons[leptons[4]]->Px(),muons[leptons[4]]->Py(), muons[leptons[4]]->Pz(),muons[leptons[4]]->Energy());     
  }

  Zbos = Zbos_lep0 + Zbos_lep1; 
  if(verbose>3) cout << " out Zboson creator " <<endl; 
  return Zbos; 
}
