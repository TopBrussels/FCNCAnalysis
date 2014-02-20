#include <cmath>
#include <fstream>
#include <sstream>
#include "TDCacheFile.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TopBrussels/TopTreeProducer/interface/TRootRun.h"
#include "TopBrussels/TopTreeProducer/interface/TRootEvent.h"
#include "TopBrussels/TopTreeProducer/interface/TRootLepton.h"
#include "TopBrussels/TopTreeProducer/interface/TRootMuon.h"
#include "TopBrussels/TopTreeProducer/interface/TRootElectron.h"
#include "TopBrussels/TopTreeProducer/interface/TRootPhoton.h"
#include "TopBrussels/TopTreeProducer/interface/TRootJet.h"
#include "TopBrussels/TopTreeProducer/interface/TRootPFJet.h"
#include "TopBrussels/TopTreeProducer/interface/TRootMET.h"

using namespace std;
using namespace TopTree;

struct HighestCSVBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
        return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};

int main (int argc, char *argv[])
{
  gSystem->Load("pluginTopBrusselsTopTreeProducer.so");
  AutoLibraryLoader::enable();

  double nevents = 0;

  TChain * t = new TChain("eventTree");
  TChain * r = new TChain("runTree");

  t->Add(argv[1]);
  r->Add(argv[1]);

  //t->Add("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/tjkim/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/V5_0_2/28102013_230503/TOPTREE/TOPTREE_10_*.root");
  //r->Add("dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/tjkim/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/V5_0_2/28102013_230503/TOPTREE/TOPTREE_10_*.root");

  TBranch* run_br = 0;
  TRootRun* runInfos = 0;
  run_br = (TBranch *) r->GetBranch("runInfos");
  run_br->SetAddress(&runInfos);
  int ntotal = runInfos->nHLTEvents();

  double nevt = t->GetEntries();
 
  TFile *fout = new TFile (argv[2], "RECREATE");

  /////////////////////
  // My tree         //
  ////////////////////
  int prePathCounter;
  double diel1_pt;
  double diel1_eta;
  double diel1_phi;
  double diel2_pt;
  double diel2_eta;
  double diel2_phi;
  double diel_mass;
  double dimu1_pt;
  double dimu1_eta;
  double dimu1_phi;
  double dimu2_pt;
  double dimu2_eta;
  double dimu2_phi;
  double dimu_mass;
  double emu1_pt;
  double emu1_eta;
  double emu1_phi;
  double emu2_pt;
  double emu2_eta;
  double emu2_phi;
  double emu_mass;

  double photon1_pt;
  double photon1_eta;
  double photon1_phi;
  double photon1_chIso;
  double photon1_phIso;
  double photon1_nhIso;
  double photon1_relIso;
  double photon2_pt;
  double photon2_eta;
  double photon2_phi;
  double photon2_chIso;
  double photon2_phIso;
  double photon2_nhIso;
  double photon2_relIso;
  double diphoton_mass;
  int nphotons;

  double jet1_pt;
  double jet1_eta;
  double jet1_phi;
  double jet1_csv;
  double jet2_pt;
  double jet2_eta;
  double jet2_phi;
  double jet2_csv;
  double jet3_pt;
  double jet3_eta;
  double jet3_phi;
  double jet3_csv;
  double jet4_pt;
  double jet4_eta;
  double jet4_phi;
  double jet4_csv;

  int njets;
  int nbjets_CSVL;
  int nbjets_CSVM;
  int nbjets_CSVT;
  double met;
  
  TTree* myTree = new TTree("tree","tree");
  
  myTree->Branch("prePathCounter", &prePathCounter, "prePathCounter/I");

  //dileptons
  myTree->Branch("diel1_pt", &diel1_pt, "diel1_pt/D");
  myTree->Branch("diel1_eta", &diel1_eta, "diel1_eta/D");
  myTree->Branch("diel1_phi", &diel1_phi, "diel1_phi/D");
  myTree->Branch("diel2_pt", &diel2_pt, "diel2_pt/D");
  myTree->Branch("diel2_eta", &diel2_eta, "diel2_eta/D");
  myTree->Branch("diel2_phi", &diel2_phi, "diel2_phi/D");
  myTree->Branch("diel_mass", &diel_mass, "diel_mass/D");
  myTree->Branch("dimu1_pt", &dimu1_pt, "dimu1_pt/D");
  myTree->Branch("dimu1_eta", &dimu1_eta, "dimu1_eta/D");
  myTree->Branch("dimu1_phi", &dimu1_phi, "dimu1_phi/D");
  myTree->Branch("dimu2_pt", &dimu2_pt, "dimu2_pt/D");
  myTree->Branch("dimu2_eta", &dimu2_eta, "dimu2_eta/D");
  myTree->Branch("dimu2_phi", &dimu2_phi, "dimu2_phi/D");
  myTree->Branch("dimu_mass", &dimu_mass, "dimu_mass/D");
  myTree->Branch("emu1_pt", &emu1_pt, "emu1_pt/D");
  myTree->Branch("emu1_eta", &emu1_eta, "emu1_eta/D");
  myTree->Branch("emu1_phi", &emu1_phi, "emu1_phi/D");
  myTree->Branch("emu2_pt", &emu2_pt, "emu2_pt/D");
  myTree->Branch("emu2_eta", &emu2_eta, "emu2_eta/D");
  myTree->Branch("emu2_phi", &emu2_phi, "emu2_phi/D");
  myTree->Branch("emu_mass", &emu_mass, "emu_mass/D");

  //photon
  myTree->Branch("nphotons", &nphotons, "nphotons/I");
  myTree->Branch("photon1_pt", &photon1_pt, "photon1_pt/D");
  myTree->Branch("photon1_eta", &photon1_eta, "photon1_eta/D");
  myTree->Branch("photon1_phi", &photon1_phi, "photon1_phi/D");
  myTree->Branch("photon1_chIso", &photon1_chIso, "photon1_chIso/D");
  myTree->Branch("photon1_phIso", &photon1_phIso, "photon1_phIso/D");
  myTree->Branch("photon1_nhIso", &photon1_nhIso, "photon1_nhIso/D");
  myTree->Branch("photon1_relIso", &photon1_relIso, "photon1_relIso/D");   
  myTree->Branch("photon2_pt", &photon2_pt, "photon2_pt/D");
  myTree->Branch("photon2_eta", &photon2_eta, "photon2_eta/D");
  myTree->Branch("photon2_phi", &photon2_phi, "photon2_phi/D");
  myTree->Branch("photon2_chIso", &photon2_chIso, "photon2_chIso/D");
  myTree->Branch("photon2_phIso", &photon2_phIso, "photon2_phIso/D");
  myTree->Branch("photon2_nhIso", &photon2_nhIso, "photon2_nhIso/D");
  myTree->Branch("photon2_relIso", &photon2_relIso, "photon2_relIso/D");
  myTree->Branch("diphoton_mass", &diphoton_mass, "diphoton_mass/D");

  //jets
  myTree->Branch("jet1_pt", &jet1_pt, "jet1_pt/D");
  myTree->Branch("jet1_eta", &jet1_eta, "jet1_eta/D");
  myTree->Branch("jet1_phi", &jet1_phi, "jet1_phi/D");
  myTree->Branch("jet1_csv", &jet1_csv, "jet1_csv/D");
  myTree->Branch("jet2_pt", &jet2_pt, "jet2_pt/D");
  myTree->Branch("jet2_eta", &jet2_eta, "jet2_eta/D");
  myTree->Branch("jet2_phi", &jet2_phi, "jet2_phi/D");
  myTree->Branch("jet2_csv", &jet2_csv, "jet2_csv/D");
  myTree->Branch("jet3_pt", &jet3_pt, "jet3_pt/D");
  myTree->Branch("jet3_eta", &jet3_eta, "jet3_eta/D");
  myTree->Branch("jet3_phi", &jet3_phi, "jet3_phi/D");
  myTree->Branch("jet3_csv", &jet3_csv, "jet3_csv/D");
  myTree->Branch("jet4_pt", &jet4_pt, "jet4_pt/D");
  myTree->Branch("jet4_eta", &jet4_eta, "jet4_eta/D");
  myTree->Branch("jet4_phi", &jet4_phi, "jet4_phi/D");
  myTree->Branch("jet4_csv", &jet4_csv, "jet4_csv/D");
  myTree->Branch("njets", &njets, "njets/I");
  myTree->Branch("nbjets_CSVL", &nbjets_CSVL, "nbjets_CSVL/I");
  myTree->Branch("nbjets_CSVM", &nbjets_CSVM, "nbjets_CSVM/I");
  myTree->Branch("nbjets_CSVT", &nbjets_CSVT, "nbjets_CSVT/I");
  
  //met
  myTree->Branch("met", &met, "met/D");

  TClonesArray *tcmuons = new TClonesArray ("TopTree::TRootMuon", 0);
  t->SetBranchAddress("Muons_selectedPatMuonsPF2PAT",&tcmuons);
  TClonesArray *tcelectrons = new TClonesArray ("TopTree::TRootElectron", 0);
  t->SetBranchAddress("Electrons_selectedPatElectronsPF2PAT",&tcelectrons);
  TClonesArray *tcphotons = new TClonesArray ("TopTree::TRootPhoton", 0);
  t->SetBranchAddress("Photons_selectedPatPhotons",&tcphotons);
  TClonesArray *tcjets = new TClonesArray ("TopTree::TRootPFJet", 0);
  t->SetBranchAddress("PFJets_selectedPatJetsPF2PAT",&tcjets);
  TClonesArray *tcmets = new TClonesArray ("TopTree::TRootPFMET", 0);
  t->SetBranchAddress("PFMET_patType1CorrectedPFMetPF2PAT",&tcmets);

  TH1F * EventSummary = new TH1F("EventSummary","EventSummary",2,0,2);
 
  prePathCounter = ntotal; 
  EventSummary->SetBinContent(1, nevt);
  ////////////////////////////////////
  //	Loop on events
  ////////////////////////////////////
 
  for (unsigned int ievt = 0; ievt < nevt; ievt++)
  {
    if( ievt % 10000 == 0 ) cout << ievt << "/" << nevt << endl;

    vector < TRootMuon* > muons;
    vector < TRootElectron* > electrons;
    vector < TRootPhoton* > photons;
    vector < TRootPFJet* > jets;
    vector < TRootMET* > mets;
      
    nevents++;
    ////////////////
    // LOAD EVENT //
    ////////////////

    t->GetEntry(ievt);

    //init
    //dileptons
    diel1_pt = -999;
    diel1_eta = -999;
    diel1_phi = -999;
    diel2_pt = -999;
    diel2_eta = -999;
    diel2_phi = -999;
    diel_mass = -999;
    dimu1_pt = -999;
    dimu1_eta = -999;
    dimu1_phi = -999;
    dimu2_pt = -999;
    dimu2_eta = -999;
    dimu2_phi = -999;
    dimu_mass = -999;
    emu1_pt = -999;
    emu1_eta = -999;
    emu1_phi = -999;
    emu2_pt = -999;
    emu2_eta = -999;
    emu2_phi = -999;
    emu_mass = -999;
    //photons
    photon1_pt = -999;
    photon1_eta = -999;
    photon1_phi = -999;
    photon1_chIso = -999;
    photon1_phIso = -999;
    photon1_nhIso = -999;
    photon1_relIso = -999;
    photon2_pt = -999;
    photon2_eta = -999;
    photon2_phi = -999;
    photon2_chIso = -999;
    photon2_phIso = -999;
    photon2_nhIso = -999;
    photon2_relIso = -999;
    diphoton_mass = -999;
    nphotons = -999;
    //jets
    jet1_pt = -999;
    jet1_eta = -999;
    jet1_phi = -999;
    jet1_csv = -999;
    jet2_pt = -999;
    jet2_eta = -999;
    jet2_phi = -999;
    jet2_csv = -999;
    jet3_pt = -999;
    jet3_eta = -999;
    jet3_phi = -999;
    jet3_csv = -999;
    jet4_pt = -999;
    jet4_eta = -999;
    jet4_phi = -999;
    jet4_csv = -999;
    njets = -999;
    //b-jets
    nbjets_CSVL = -999;
    nbjets_CSVM = -999;
    nbjets_CSVT = -999;
    //met
    met = -999;

    //clear vectors   
    muons.clear();
    electrons.clear();
    photons.clear();
    jets.clear();

    for (int i = 0; i < tcmuons->GetEntriesFast(); i++)
      muons.push_back ((TRootMuon *) tcmuons->At(i));
    for (int i = 0; i < tcelectrons->GetEntriesFast(); i++){

      TRootElectron * tcelectron = (TRootElectron *) tcelectrons->At(i);

      bool pass = tcelectron->Pt() > 20 && tcelectron->mvaTrigId() > 0.5 && tcelectron->relPfIso(3,0.5) < 0.15;
      if( pass ) {
        electrons.push_back ((TRootElectron *) tcelectrons->At(i));
      }
    }

    for (int i = 0; i < tcphotons->GetEntriesFast(); i++){

     TRootPhoton * tcphoton = (TRootPhoton *) tcphotons->At(i);
 
      bool pass = false;

      if( abs(tcphoton->Eta()) < 1.479 ){
        pass = tcphoton->Pt() > 20 && tcphoton->sigmaIetaIeta() < 0.011 && tcphoton->hadronicOverEm() < 0.05 && tcphoton->passelectronveto();
      }else{
        pass = tcphoton->Pt() > 20 && tcphoton->sigmaIetaIeta() < 0.033 && tcphoton->hadronicOverEm() < 0.05 && tcphoton->passelectronveto();
      }

      if( pass ){
        photons.push_back ((TRootPhoton *) tcphotons->At(i));
      }
    }

    if( photons.size() < 2 ) continue; 

    for (int i = 0; i < tcjets->GetEntriesFast(); i++){
      TRootPFJet * tcjet = (TRootPFJet *) tcjets->At(i);
      bool pass = true;
 
      if( tcjet->Pt() <= 20 ) continue;

      for(int j=0; j < (int) electrons.size()-1; j++){
        double dR = sqrt( (electrons[j]->Eta()-tcjet->Eta())*(electrons[j]->Eta()-tcjet->Eta()) + (electrons[j]->Phi()-tcjet->Phi())*(electrons[j]->Phi()-tcjet->Phi()) );
        if( dR < 0.5 ) {
          pass = false; 
          break;
        }
      }
      
      if( pass ){
        jets.push_back ((TRootPFJet *) tcjets->At(i));
      }
    }
    //------------------------------//
    // re-arrange jets in CSV order //
    //------------------------------//
    sort(jets.begin(),jets.end(),HighestCSVBtag());
   
    for (int i = 0; i < tcmets->GetEntriesFast(); i++){
      mets.push_back ((TRootMET *) tcmets->At(i));
    }

    // scale factor for the event
    //float scaleFactor = 1.;
      
    //-----------------//
    // test saving tree//
    //-----------------//
    if( electrons.size() > 1){
      for(int i=0; i < (int) electrons.size()-1; i++){
        for(int j=i+1; j < (int) electrons.size(); j++){
          bool electron1_pass = electrons[i]->Pt() > 20 && electrons[i]->mvaTrigId() > 0.5 && electrons[i]->relPfIso(3,0.5) < 0.15;
          bool electron2_pass = electrons[j]->Pt() > 20 && electrons[j]->mvaTrigId() > 0.5 && electrons[j]->relPfIso(3,0.5) < 0.15;
          if( electron1_pass && electron2_pass ){
            diel1_pt = electrons[i]->Pt();
            diel1_eta = electrons[i]->Eta();
            diel1_phi = electrons[i]->Phi();
            diel2_pt = electrons[j]->Pt();
            diel2_eta = electrons[j]->Eta();
            diel2_phi = electrons[j]->Phi();
            TLorentzVector v1(electrons[i]->Px(), electrons[i]->Py(), electrons[i]->Pz(), electrons[i]->Energy());
            TLorentzVector v2(electrons[j]->Px(), electrons[j]->Py(), electrons[j]->Pz(), electrons[j]->Energy());
            diel_mass = (v1+v2).M();
          }
        }
      } 
    }

    if( muons.size() > 1){
      for(int i=0; i < (int) muons.size()-1; i++){
        for(int j=i+1; j < (int) muons.size(); j++){
          bool muon1_pass = muons[i]->Pt() > 20 && muons[i]->isPFMuon() && ( muons[i]->isGlobalMuon() || muons[i]->isTrackerMuon() ) && muons[i]->relPfIso(3,0.5) < 0.15;
          bool muon2_pass = muons[j]->Pt() > 20 && muons[j]->isPFMuon() && ( muons[j]->isGlobalMuon() || muons[j]->isTrackerMuon() ) && muons[j]->relPfIso(3,0.5) < 0.15;
          if( muon1_pass && muon2_pass ){
            dimu1_pt = muons[i]->Pt();
            dimu1_eta = muons[i]->Eta();
            dimu1_phi = muons[i]->Phi();
            dimu2_pt = muons[j]->Pt();
            dimu2_eta = muons[j]->Eta();
            dimu2_phi = muons[j]->Phi();
            TLorentzVector v1(muons[i]->Px(), muons[i]->Py(), muons[i]->Pz(), muons[i]->Energy());
            TLorentzVector v2(muons[j]->Px(), muons[j]->Py(), muons[j]->Pz(), muons[j]->Energy());
            dimu_mass = (v1+v2).M(); 
          }
        }
      } 
    }

    if( muons.size() > 0 && electrons.size() > 0){
      for(int i=0; i < (int) electrons.size(); i++){
        for(int j=0; j < (int) muons.size(); j++){
          bool electron_pass = electrons[i]->Pt() > 20 && electrons[i]->mvaTrigId() > 0.5 && electrons[i]->relPfIso(3,0.5) < 0.15;
          bool muon_pass = muons[j]->Pt() > 20 && muons[j]->isPFMuon() && ( muons[j]->isGlobalMuon() || muons[j]->isTrackerMuon() ) && muons[j]->relPfIso(3,0.5) < 0.15;

          if( electron_pass && muon_pass ){
            emu1_pt = electrons[i]->Pt();
            emu1_eta = electrons[i]->Eta();
            emu1_phi = electrons[i]->Phi();
            emu2_pt = muons[j]->Pt();
            emu2_eta = muons[j]->Eta();
            emu2_phi = muons[j]->Phi();
            TLorentzVector v1(electrons[i]->Px(), electrons[i]->Py(), electrons[i]->Pz(), electrons[i]->Energy());
            TLorentzVector v2(muons[j]->Px(), muons[j]->Py(), muons[j]->Pz(), muons[j]->Energy());
            emu_mass = (v1+v2).M();
          }
        }
      } 
    }

    //if( dimu_mass <= 0 && diel_mass <= 0 && emu_mass <= 0) continue;
  
    if( photons.size() > 1){
      for(int i=0; i < (int) photons.size()-1; i++){
        for(int j=i+1; j < (int) photons.size(); j++){
     
          bool photon1_pass = false;
          bool photon2_pass = false;
          //not recommended but useful to keep as it is reported to be very efficient cut
          bool photon1_hasPixelSeed_pass =false; 
          bool photon2_hasPixelSeed_pass =false; 

          if( abs(photons[i]->Eta()) < 1.479 ){
            photon1_pass = photons[i]->sigmaIetaIeta() < 0.011 && photons[i]->hadronicOverEm() < 0.05 && photons[i]->passelectronveto();
          }else{
            photon1_pass = photons[i]->sigmaIetaIeta() < 0.033 && photons[i]->hadronicOverEm() < 0.05 && photons[i]->passelectronveto();
          }  

          if( abs(photons[j]->Eta()) < 1.479 ){
            photon2_pass = photons[j]->sigmaIetaIeta() < 0.011 && photons[j]->hadronicOverEm() < 0.05 && photons[j]->passelectronveto();
          }else{
            photon2_pass = photons[j]->sigmaIetaIeta() < 0.033 && photons[j]->hadronicOverEm() < 0.05 && photons[j]->passelectronveto();
          }
          
          if( photon1_pass && photon2_pass){
            photon1_pt = photons[i]->Pt();
            photon1_eta = photons[i]->Eta();
            photon1_phi = photons[i]->Phi();
            photon1_chIso = photons[i]->chargedHadronIso();
            photon1_nhIso = photons[i]->neutralHadronIso();
            photon1_phIso = photons[i]->photonIso();
            photon1_relIso = photons[i]->relPfIso();
            photon2_pt = photons[j]->Pt();
            photon2_eta = photons[j]->Eta();
            photon2_phi = photons[j]->Phi();
            photon2_chIso = photons[j]->chargedHadronIso();
            photon2_nhIso = photons[j]->neutralHadronIso();
            photon2_phIso = photons[j]->photonIso();
            photon2_relIso = photons[j]->relPfIso();
            TLorentzVector v1(photons[i]->Px(), photons[i]->Py(), photons[i]->Pz(), photons[i]->Energy());
            TLorentzVector v2(photons[j]->Px(), photons[j]->Py(), photons[j]->Pz(), photons[j]->Energy());
            diphoton_mass = (v1+v2).M();
          }
        }
      }
    }

    nphotons = photons.size();

    int nb_CSVL = 0;
    int nb_CSVM = 0;
    int nb_CSVT = 0;

    for(unsigned int i=0; i < jets.size() ; i++){
      if( jets[i]->btag_combinedSecondaryVertexBJetTags() > 0.244 ){
        nb_CSVL = nb_CSVL+1;
      }
      if( jets[i]->btag_combinedSecondaryVertexBJetTags() > 0.679 ){
        nb_CSVM = nb_CSVM+1;
      }
      if( jets[i]->btag_combinedSecondaryVertexBJetTags() > 0.898 ){
        nb_CSVT = nb_CSVT+1;
      }
    }

    nbjets_CSVL = nb_CSVL;
    nbjets_CSVM = nb_CSVM;
    nbjets_CSVT = nb_CSVT;

    if( jets.size() > 0){
      jet1_pt = jets[0]->Pt();
      jet1_eta = jets[0]->Eta();
      jet1_phi = jets[0]->Phi();
      jet1_csv = jets[0]->btag_combinedSecondaryVertexBJetTags(); 
    }
    if( jets.size() > 1){  
      jet2_pt = jets[1]->Pt();
      jet2_eta = jets[1]->Eta();
      jet2_phi = jets[1]->Phi();
      jet2_csv = jets[1]->btag_combinedSecondaryVertexBJetTags();
    }
    if( jets.size() > 2){
      jet3_pt = jets[2]->Pt();
      jet3_eta = jets[2]->Eta();
      jet3_phi = jets[2]->Phi();
      jet3_csv = jets[2]->btag_combinedSecondaryVertexBJetTags();
    }
    if( jets.size() > 3){
      jet4_pt = jets[3]->Pt();
      jet4_eta = jets[3]->Eta();
      jet4_phi = jets[3]->Phi();
      jet4_csv = jets[3]->btag_combinedSecondaryVertexBJetTags();
    }

    njets = jets.size();
    met = mets[0]->Pt();

    myTree->Fill();
  }			//loop on events
  //cout << "nevents= " << nevents << endl;

  fout->Write();

  delete fout; 

  tcmuons->Delete();
  tcelectrons->Delete();

  return 0;
}
