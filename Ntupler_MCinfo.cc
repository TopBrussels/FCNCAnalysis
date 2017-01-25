#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <iostream>
#include <map>
#include <cstdlib>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"

#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"
#include "TopTreeAnalysisBase/Tools/interface/SourceDate.h"
#include "TopTreeAnalysisBase/Tools/interface/Trigger.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"

#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

map<string,TH1F*> histo1D;

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 11)
    {
        std::cerr << "INVALID INPUT FROM XMLFILE.  CHECK XML IMPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
        return 1;
    }
    //Working command: ./MACRO_DrawMC_PU TT_powheg TTjets 3 1 1 1 1 1 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_80X_v0/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_ring1680XmcRun2asymptotic2016miniAODv2v0ext4v180XmcRun2asymptoticRealisticBS25ns13TeV2016v1mccrab26/160727_222013/0000/TOPTREE_ 500 1000000

    //Placing arguments in properly typed variables for Dataset creation

    const string dName              = argv[1];
    const string dTitle             = argv[2];
    const int color                 = strtol(argv[3], NULL, 10);
    const int ls                    = strtol(argv[4], NULL, 10);
    const int lw                    = strtol(argv[5], NULL, 10);
    const float normf               = strtod(argv[6], NULL);
    const float EqLumi              = strtod(argv[7], NULL);
    const float xSect               = strtod(argv[8], NULL);
    const string location = argv[9];
    const int NumberOfFiles = strtol(argv[10], NULL, 10);
    
   
    
    


  vector<string> vecfileNames;
  for(int i = 1; i <= NumberOfFiles; i++)
  {
        std::ostringstream intToStr;
    	  intToStr<<i;
    	  
    	  std::string FileName = location+intToStr.str()+".root";

        vecfileNames.push_back(FileName);
  }





    /////////////////////////////////
    //  Set up AnalysisEnvironment 
    /////////////////////////////////
    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_calibratedPatElectrons";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;

    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_fatjets;
    vector < TRootMET* >      mets;
    TRootEvent* event = 0;

 
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset

	      // variables for ntuple
        histo1D["eta_Higgs"] = new TH1F("eta_Higgs", "eta_Higgs", 100, -3, 3);
        histo1D["eta_Top"] = new TH1F("eta_Top", "eta_Top", 100, -3, 3);
        histo1D["eta_AntiTop"] = new TH1F("eta_AntiTop", "eta_AntiTop", 100, -3, 3);

        histo1D["phi_Higgs"] = new TH1F("phi_Higgs", "phi_Higgs", 100, -3.2, 3.2);
        histo1D["phi_Top"] = new TH1F("phi_Top", "phi_Top", 100, -3.2, 3.2);
        histo1D["phi_AntiTop"] = new TH1F("phi_AntiTop", "phi_AntiTop", 100, -3.2, 3.2);

        histo1D["pt_Higgs"] = new TH1F("pt_Higgs", "pt_Higgs", 100, 0., 500.);
        histo1D["pt_Top"] = new TH1F("pt_Top", "pt_Top", 100, 0., 500.);
        histo1D["pt_AntiTop"] = new TH1F("pt_AntiTop", "pt_AntiTop", 100, 0., 500);

        histo1D["DeltaR_Higgs_Top"] = new TH1F("DeltaR_Higgs_Top", "DeltaR_Higgs_Top", 100, 2., 5.5);
        histo1D["DeltaR_Higgs_AntiTop"] = new TH1F("DeltaR_Higgs_AntiTop", "DeltaR_Higgs_AntiTop", 100, 2., 5.5);

        histo1D["eta_u_TopMother"] = new TH1F("eta_u_TopMother", "eta_u_TopMother", 100, -3, 3);
        histo1D["eta_b_TopMother"] = new TH1F("eta_b_TopMother", "eta_b_TopMother", 100, -3, 3);
        histo1D["eta_c_TopMother"] = new TH1F("eta_c_TopMother", "eta_c_TopMother", 100, -3, 3);

        histo1D["phi_u_TopMother"] = new TH1F("phi_u_TopMother", "phi_u_TopMother", 100, -3.2, 3.2);
        histo1D["phi_b_TopMother"] = new TH1F("phi_b_TopMother", "phi_b_TopMother", 100, -3.2, 3.2);
        histo1D["phi_c_TopMother"] = new TH1F("phi_c_TopMother", "phi_c_TopMother", 100, -3.2, 3.2);

        histo1D["pt_u_TopMother"] = new TH1F("pt_u_TopMother", "pt_u_TopMother", 100, 0., 500.);
        histo1D["pt_b_TopMother"] = new TH1F("pt_b_TopMother", "pt_b_TopMother", 100, 0., 500.);
        histo1D["pt_c_TopMother"] = new TH1F("pt_c_TopMother", "pt_c_TopMother", 100, 0., 500);

        histo1D["eta_u_AntiTopMother"] = new TH1F("eta_u_AntiTopMother", "eta_u_AntiTopMother", 100, -3, 3);
        histo1D["eta_b_AntiTopMother"] = new TH1F("eta_b_AntiTopMother", "eta_b_AntiTopMother", 100, -3, 3);
        histo1D["eta_c_AntiTopMother"] = new TH1F("eta_c_AntiTopMother", "eta_c_AntiTopMother", 100, -3, 3);

        histo1D["phi_u_AntiTopMother"] = new TH1F("phi_u_AntiTopMother", "phi_u_AntiTopMother", 100, -3.2, 3.2);
        histo1D["phi_b_AntiTopMother"] = new TH1F("phi_b_AntiTopMother", "phi_b_AntiTopMother", 100, -3.2, 3.2);
        histo1D["phi_c_AntiTopMother"] = new TH1F("phi_c_AntiTopMother", "phi_c_AntiTopMother", 100, -3.2, 3.2);

        histo1D["pt_u_AntiTopMother"] = new TH1F("pt_u_AntiTopMother", "pt_u_AntiTopMother", 100, 0., 500.);
        histo1D["pt_b_AntiTopMother"] = new TH1F("pt_b_AntiTopMother", "pt_b_AntiTopMother", 100, 0., 500.);
        histo1D["pt_c_AntiTopMother"] = new TH1F("pt_c_AntiTopMother", "pt_c_AntiTopMother", 100, 0., 500);

        histo1D["pt_b_HiggsMother"] = new TH1F("pt_b_HiggsMother", "pt_b_HiggsMother", 100, 0., 500.);
        histo1D["eta_b_HiggsMother"] = new TH1F("eta_b_HiggsMother", "eta_b_HiggsMother", 100, -3, 3);
        histo1D["phi_b_HiggsMother"] = new TH1F("phi_b_HiggsMother", "phi_b_HiggsMother", 100, -3.2, 3.2);

        histo1D["pt_antib_HiggsMother"] = new TH1F("pt_antib_HiggsMother", "pt_antib_HiggsMother", 100, 0., 500.);
        histo1D["eta_antib_HiggsMother"] = new TH1F("eta_antib_HiggsMother", "eta_antib_HiggsMother", 100, -3, 3);
        histo1D["phi_antib_HiggsMother"] = new TH1F("phi_antib_HiggsMother", "phi_antib_HiggsMother", 100, -3.2, 3.2);

        histo1D["pt_AllMCParticles"] = new TH1F("pt_AllMCParticles", "pt_AllMCParticles", 100, 0., 500.);

        histo1D["DeltaR_b_antib_FromHiggs"] = new TH1F("DeltaR_b_antib_FromHiggs", "DeltaR_b_antib_FromHiggs", 100, 0., 5.5);
        histo1D["DeltaR_bFromHiggs_Higgs"] = new TH1F("DeltaR_bFromHiggs_Higgs", "DeltaR_bFromHiggs_Higgs", 100, 0., 5.5);
        histo1D["DeltaR_antibFromHiggs_Higgs"] = new TH1F("DeltaR_antibFromHiggs_Higgs", "DeltaR_antibFromHiggs_Higgs", 100, 0., 5.5);

        
        int numberOfEventsToRunOver = datasets[d]->NofEvtsToRunOver();

        for (unsigned int ievt = 0; ievt < numberOfEventsToRunOver; ievt++)
        {
            vector<TLorentzVector> mcParticlesTLV, mcMuonsTLV, mcPartonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_,mcParticles;
            vector<int> mcMuonIndex, mcPartonIndex;
            JetPartonMatching muonMatching, jetMatching;


            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, false);  //load event
            datasets[d]->eventTree()->LoadTree(ievt);

            mcParticles.clear();
            treeLoader.LoadMCEvent(ievt, 0, mcParticles, false);
            sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
                  
//cout  << " mcParticles.size(): " << mcParticles.size() << endl;      
                bool foundTop = false;
                bool foundAntiTop = false;
                bool foundHiggs = false;
                bool foundbFromHiggs = false;
                bool foundantibFromHiggs = false;
                
                double eta_Top = 0;
                double eta_AntiTop = 0;
                double eta_Higgs = 0;
                double phi_Top = 0;
                double phi_AntiTop = 0;
                double phi_Higgs = 0;
                double phi_b = 0;
                double phi_antib = 0;
                double eta_b = 0;
                double eta_antib = 0;
                
                for (unsigned int i = 0; i < mcParticles.size(); i++)
                {
                    if ( (mcParticles[i]->status() > 1 && mcParticles[i]->status() <= 20) || mcParticles[i]->status() >= 30 ) continue;  /// Final state particle or particle from hardest process                   
//                    if ( mcParticles[i]->status() != 1 ) continue;  /// Final state particle or particle from hardest process                   
                    histo1D["pt_AllMCParticles"]->Fill(mcParticles[i]->Pt());
//cout << "mcParticles[i]->Pt(): " << mcParticles[i]->Pt() << endl;
                    if(mcParticles[i]->type() == 6)
                    {
                        histo1D["eta_Top"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_Top"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_Top"]->Fill(mcParticles[i]->Pt());
                        foundTop = true;
                        
                        eta_Top = mcParticles[i]->Eta();
                        phi_Top = mcParticles[i]->Phi();
                    }
                    else if(mcParticles[i]->type() == -6)
                    {
                        histo1D["eta_AntiTop"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_AntiTop"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_AntiTop"]->Fill(mcParticles[i]->Pt());
                        foundAntiTop = true;
                        
                        eta_AntiTop = mcParticles[i]->Eta();
                        phi_AntiTop = mcParticles[i]->Phi();
                    }
                    else if(mcParticles[i]->type() == 25)
                    {
                        histo1D["eta_Higgs"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_Higgs"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_Higgs"]->Fill(mcParticles[i]->Pt());
                        foundHiggs = true;
                        
                        eta_Higgs = mcParticles[i]->Eta();
                        phi_Higgs = mcParticles[i]->Phi();
                    }
                    else if(mcParticles[i]->type() == 5 && mcParticles[i]->motherType() == 6)
                    {
                        histo1D["eta_b_TopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_b_TopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_b_TopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == -5  && mcParticles[i]->motherType() == -6)
                    {
                        histo1D["eta_b_AntiTopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_b_AntiTopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_b_AntiTopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == 4 && mcParticles[i]->motherType() == 6)
                    {
                        histo1D["eta_c_TopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_c_TopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_c_TopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == -4 && mcParticles[i]->motherType() == -6)
                    {
                        histo1D["eta_c_AntiTopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_c_AntiTopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_c_AntiTopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == 2 && mcParticles[i]->motherType() == 6)
                    {
                        histo1D["eta_u_TopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_u_TopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_u_TopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == -2 && mcParticles[i]->motherType() == -6)
                    {
                        histo1D["eta_u_AntiTopMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_u_AntiTopMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_u_AntiTopMother"]->Fill(mcParticles[i]->Pt());
                    }
                    else if(mcParticles[i]->type() == 5 && mcParticles[i]->motherType() == 25)
                    {
                        histo1D["eta_b_HiggsMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_b_HiggsMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_b_HiggsMother"]->Fill(mcParticles[i]->Pt());
                        
                        foundbFromHiggs = true;
                        phi_b = mcParticles[i]->Phi();
                        eta_b = mcParticles[i]->Eta();
                    }
                    else if(mcParticles[i]->type() == -5  && mcParticles[i]->motherType() == 25)
                    {
                        histo1D["eta_antib_HiggsMother"]->Fill(mcParticles[i]->Eta());
                        histo1D["phi_antib_HiggsMother"]->Fill(mcParticles[i]->Phi());
                        histo1D["pt_antib_HiggsMother"]->Fill(mcParticles[i]->Pt());

                        foundantibFromHiggs = true;
                        phi_antib = mcParticles[i]->Phi();
                        eta_antib = mcParticles[i]->Eta();
                    }
                }
                
                if(foundTop && foundHiggs)
                {
                    double deltaR = sqrt( pow(eta_Top-eta_Higgs,2) + pow(phi_Top-phi_Higgs,2));
                    histo1D["DeltaR_Higgs_Top"]->Fill(deltaR);
                }
                if(foundAntiTop && foundHiggs)
                {
                    double deltaR = sqrt( pow(eta_AntiTop-eta_Higgs,2) + pow(phi_AntiTop-phi_Higgs,2));
                    histo1D["DeltaR_Higgs_AntiTop"]->Fill(deltaR);
                }
                if(foundbFromHiggs && foundHiggs)
                {
                    double deltaR = sqrt( pow(eta_b-eta_Higgs,2) + pow(phi_b-phi_Higgs,2));
                    histo1D["DeltaR_bFromHiggs_Higgs"]->Fill(deltaR);
                }
                if(foundantibFromHiggs && foundHiggs)
                {
                    double deltaR = sqrt( pow(eta_antib-eta_Higgs,2) + pow(phi_antib-phi_Higgs,2));
                    histo1D["DeltaR_antibFromHiggs_Higgs"]->Fill(deltaR);
                }
                if(foundbFromHiggs && foundantibFromHiggs)
                {
                    double deltaR = sqrt( pow(eta_b-eta_antib,2) + pow(phi_b-phi_antib,2));
                    histo1D["DeltaR_b_antib_FromHiggs"]->Fill(deltaR);
                }

        }
      
  }

  TFile * fout = new TFile(("MC_ComparisonOutput/MC_Comparison_"+dName+".root").c_str(),"RECREATE");

  fout->cd();
  
  for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {
     	string name = it->first;
     	TH1F *temp = it->second;

      temp->Write();

  }
  fout->Write();
}
