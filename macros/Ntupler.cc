///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// Ntupler.cc: This macro is intended to be an example analysis macro which works out of the box.           /////
/////       It should serve as the first port of call for new users of the TopBrussels framework.             /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TStyle.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

//user code
#include "../../TopTreeProducer/interface/TRootRun.h"
#include "../../TopTreeProducer/interface/TRootEvent.h"
#include "../../TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "../../TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "../../TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "../../TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "../../TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "../../TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "../../TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/MCWeighter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/ElectronPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/MuonPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/JetPlotter.h"
#include "../../TopTreeAnalysisBase/Selection/interface/VertexPlotter.h"
#include "../../TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"
#include "../../TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "../../TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "../../TopTreeAnalysis/macros/Style.C"
#include "../../TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"

using namespace std;
using namespace reweight;
using namespace TopTree;

//for ordering the leptons
bool cmp_big_first (pair<int,int> i, pair<int,int> j) { return i.second > j.second; }

//Calculating MEz
pair<double,double> MEzCalculator(TLorentzVector lepton_ , double missingEt_pX_, double missingEt_pY_, bool muon)
{
  double M_W = 80.4;
  double M_mu = 0.10566; //105 MeV 
  double M_e = 0.000511; //0.511 MeV 
  double emu = lepton_.E();
  double pxmu = lepton_.Px();
  double pymu = lepton_.Py();
  double pzmu = lepton_.Pz();
  double pxnu = missingEt_pX_; //MET_.Px();
  double pynu = missingEt_pY_; //MET_.Py();
  double pznu = 0.;
  double pznu1 = 0.;
  double pznu2 = 0.;
  double M_lep = 0.;
  bool isComplex_ = false; 
  
  if(muon)   M_lep = M_mu; 
  else M_lep = M_e;
  
  double a = M_W*M_W - M_lep*M_lep + 2.0*pxmu*pxnu + 2.0*pymu*pynu;
  double A = 4.0*(emu*emu - pzmu*pzmu);
  double B = -4.0*a*pzmu;
  double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
  
  double tmproot = B*B - 4.0*A*C;
  
  if (tmproot<0) {
    isComplex_= true;
    pznu = - B/(2*A); // take real part of complex roots
    pznu1 = pznu; 
    pznu2 = pznu; 
  }
  else {
    isComplex_ = false;
    double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
    double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
    
    pznu1 = tmpsol1; 
    pznu2 = tmpsol2; 
/*    if (type == 0 ) {
      // two real roots, pick the one closest to pz of muon
      if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
      else pznu = tmpsol1;
      // if pznu is > 300 pick the most central root
      if ( pznu > 300. ) {
		if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
		else pznu = tmpsol2;
      }
    }
   if (type == 1 ) {
     // two real roots, pick the one closest to pz of muon
      if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
      else pznu = tmpsol1;
    }
    if (type == 2 ) {
      // pick the most central root.
      if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
      else pznu = tmpsol2;
    }
*/  }
  
  //Particle neutrino;
  //neutrino.setP4( LorentzVector(pxnu, pynu, pznu, TMath::Sqrt(pxnu*pxnu + pynu*pynu + pznu*pznu ))) ;
    pair <double, double> tempPair; 
    tempPair.first = pznu1; 
    tempPair.second = pznu2;
  return tempPair;
}

//Leptons ordering according to Z decay for 45 channel
vector <TLorentzVector> Channel_45_5leptons( bool _debug, vector <TRootElectron*> _Electrons, vector <TRootMuon*> _Muons,double _missingEt_pX, double _missingEt_pY, double _missingEt_Theta, double _missingEt)
{
      bool debug = _debug; 	     
      TLorentzVector leptonpair_1; // two OSSF with highest pt
      TLorentzVector leptonpair_2; // two OSSF with 2nd highest pt
      TLorentzVector leptonfour; // decay of Z
      TLorentzVector lepton_0;
      TLorentzVector lepton_1;
      TLorentzVector lepton_2;
      TLorentzVector lepton_3;
      TLorentzVector lepton_4;
      
      lepton_0.Clear(); 
      lepton_1.Clear();
      lepton_2.Clear();
      lepton_3.Clear();
      lepton_4.Clear();
      leptonpair_1.Clear();
      leptonpair_2.Clear();
      leptonfour.Clear();
      
      
      
      bool leptons = false; 
      bool chargereq1=false; 
      bool chargereq2=false;
      bool chargereq3=false;
      bool chargereq4=false;
      bool chargereq5=false;
      bool chargereq6=false;
      bool chargereq7=false;
      bool chargereq8=false;
      bool chargereq9=false;
      bool chargereq10=false;
      
      double _nMuons = _Muons.size(); 
      double _nElectrons = _Electrons.size(); 

      pair <double, double> missingEt_pZ_pair; 


      if( _nMuons == 2 & _nElectrons == 3)
      { 
	   if(debug) cout << "In _nMuons == 2 & _nElectrons == 3 " << endl; 
	   if(_Muons[0]->charge()!=_Muons[1]->charge())	chargereq1 = true; 
	   if(_Electrons[0]->charge()!=_Electrons[1]->charge())	chargereq2 = true;
	   if(_Electrons[0]->charge()!=_Electrons[2]->charge())	chargereq3 = true;  
	   if(_Electrons[2]->charge()!=_Electrons[1]->charge())	chargereq4 = true; 
	   bool isMuon = false; 
	 
	   if(chargereq1 && chargereq2 && !leptons) // Mu0 and Mu1    El0 and El1
	   {
	 	  leptons=true; 
	 	   
	 	  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 	   
	 	  lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	   
	 	  lepton_2.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 	   
	 	  lepton_3.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 	   
	 	  lepton_4.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	  
	 
	   }
	   if(chargereq1 && chargereq3 && !leptons) // Mu0 and Mu1    El0 and El2
	   {
	 	  leptons=true;    
	 	  
	 	  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 	  lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	  lepton_2.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 	  lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	  lepton_4.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	  
	   }
	   if(chargereq1 && chargereq4 && !leptons) // Mu0 and Mu1    El2 and El1
	   {
	 	  leptons=true; 
	 	  
	 	  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 	  lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	  lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 	  lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	  lepton_4.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	  
	    }
	     
	    leptonpair_1 = lepton_0 + lepton_1; 
	    leptonpair_2 = lepton_2 + lepton_3;
	    leptonfour = leptonpair_1 + leptonpair_2;
	    missingEt_pZ_pair = MEzCalculator(lepton_4, _missingEt_pX, _missingEt_pY, isMuon); 
	 
	 
      } 
      if( _nMuons == 3 & _nElectrons == 2)
      {
	      if(debug) cout << "In _nMuons == 3 & _nElectrons == 2 " << endl; 
	      if(_Electrons[0]->charge()!=_Electrons[1]->charge())   chargereq1 = true; 
	      if(_Muons[0]->charge()!=_Muons[1]->charge())   chargereq2 = true;
	      if(_Muons[0]->charge()!=_Muons[2]->charge())   chargereq3 = true;  
	      if(_Muons[2]->charge()!=_Muons[1]->charge())   chargereq4 = true; 
	      bool isMuon = true; 
	      
	      if(chargereq1 && chargereq2 && !leptons) //El0 and El1 Mu0 and Mu1
	      {
	    	     leptons=true; 

	    	     lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	    	     lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	    	     lepton_2.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	    	     lepton_3.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	    	     lepton_4.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());

	    	     
	      
	      }
	      else if(chargereq1 && chargereq3 && !leptons) //El0 and El1 Mu0 and Mu2
	      {
	    	     leptons=true; 

	 		 lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	  
	 	  }
	 	  else if(chargereq1 && chargereq4 && !leptons) //El0 and El1 Mu2 and Mu1
	 	  {
	 		 leptons=true; 
	  
	 		 lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	  
	 	  }
	 	  leptonpair_1 = lepton_0 + lepton_1; 
	 	  leptonpair_2 = lepton_2 + lepton_3;
	 	  leptonfour = leptonpair_1 + leptonpair_2;
	 	  missingEt_pZ_pair = MEzCalculator(lepton_4, _missingEt_pX, _missingEt_pY, isMuon);
	  }
	  if( _nMuons == 5 )
	  {
	 	  if(debug) cout << "In _nMuons == 5 " << endl; 
	 	  if(_Muons[0]->charge()!=_Muons[1]->charge())   chargereq1 = true; 
	 	  if(_Muons[0]->charge()!=_Muons[2]->charge())   chargereq2 = true;
	 	  if(_Muons[0]->charge()!=_Muons[3]->charge())   chargereq3 = true;
	 	  if(_Muons[0]->charge()!=_Muons[4]->charge())   chargereq4 = true;  
	 	  if(_Muons[2]->charge()!=_Muons[1]->charge())   chargereq5 = true;
	 	  if(_Muons[3]->charge()!=_Muons[1]->charge())   chargereq6 = true;
	 	  if(_Muons[4]->charge()!=_Muons[1]->charge())   chargereq7 = true;
	 	  if(_Muons[2]->charge()!=_Muons[3]->charge())   chargereq8 = true; 
	 	  if(_Muons[2]->charge()!=_Muons[4]->charge())   chargereq9 = true;   
	 	  if(_Muons[4]->charge()!=_Muons[3]->charge())   chargereq10 = true; 
	 	  bool isMuon = true; 
	 	  
	 	  if(chargereq1 && chargereq8  && !leptons) //0 and 1	  2 and 3
	 	  {
	 		 leptons=true; 
	  
	 		 lepton_3.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 
	 	  }
	 	  if(chargereq1 && chargereq9  && !leptons) //0 and 1	  2 and 4
	 	  {
	 		 leptons=true; 
	  
	 		 lepton_3.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 	  }
	 	  if(chargereq1 && chargereq10  && !leptons) //0 and 1     3 and 4
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 	  }
	 	  if(chargereq2 && chargereq6  && !leptons) //0 and 2	  1 and 3
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 	  }
	 	  if(chargereq2 && chargereq7  && !leptons) //0 and 2	  1 and 4
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 	  }
	 	  if(chargereq2 && chargereq10  && !leptons) //0 and 2     3 and 4
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	  }
	 	  if(chargereq3 && chargereq5  && !leptons) //0 and 3	  1 and 2
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 	  }
	 	  if(chargereq3 && chargereq7  && !leptons) //0 and 3	  1 and 4
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 	  }
	 	  if(chargereq3 && chargereq9  && !leptons) //0 and 3	  4 and 2
	 	  {
	 		 leptons=true; 
	 	  
	 		 lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		 lepton_2.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		 lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		 lepton_1.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		 lepton_4.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	  }
	 	  if(chargereq4 && chargereq5  && !leptons) //0 and 4	   1 and 2
	 	  {
	 		  leptons=true; 
	 	  
	 		  lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		  lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		  lepton_1.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		  lepton_4.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 	  }
	 	  if(chargereq4 && chargereq6  && !leptons) //0 and 4	   1 and 3
	 	  {
	 		  leptons=true; 
	 	  
	 		  lepton_3.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		  lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 		  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		  lepton_1.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		  lepton_4.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 	  }
	 	  if(chargereq4 && chargereq8  && !leptons) //0 and 4	   3 and 2
	 	  {
	 		  leptons=true; 
	 	  
	 		  lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	 		  lepton_2.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	 		  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	 		  lepton_1.SetPxPyPzE(_Muons[4]->Px(),_Muons[4]->Py(),_Muons[4]->Pz(),_Muons[4]->Energy());
	 		  lepton_4.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	 	  }
	 	  
	 	  leptonpair_1 = lepton_0 + lepton_1; 
	 	  leptonpair_2 = lepton_2 + lepton_3;
	 	  leptonfour = leptonpair_1 + leptonpair_2;
	 	  missingEt_pZ_pair = MEzCalculator(lepton_4, _missingEt_pX, _missingEt_pY, isMuon);
	 	  
	  }
	   if( _nElectrons == 5)
	   {
	 	 if(debug) cout << "_nElectrons == 5 " << endl; 
	 	 if(_Electrons[0]->charge()!=_Electrons[1]->charge())   chargereq1 = true; 
	 	 if(_Electrons[0]->charge()!=_Electrons[2]->charge())   chargereq2 = true;
	 	 if(_Electrons[0]->charge()!=_Electrons[3]->charge())   chargereq3 = true;
	 	 if(_Electrons[0]->charge()!=_Electrons[4]->charge())   chargereq4 = true;  
	 	 if(_Electrons[2]->charge()!=_Electrons[1]->charge())   chargereq5 = true;
	 	 if(_Electrons[3]->charge()!=_Electrons[1]->charge())   chargereq6 = true;
	 	 if(_Electrons[4]->charge()!=_Electrons[1]->charge())   chargereq7 = true;
	 	 if(_Electrons[2]->charge()!=_Electrons[3]->charge())   chargereq8 = true; 
	 	 if(_Electrons[2]->charge()!=_Electrons[4]->charge())   chargereq9 = true;   
	 	 if(_Electrons[4]->charge()!=_Electrons[3]->charge())   chargereq10 = true; 
	 	 bool isMuon = false; 
	 	 
	 	 if(chargereq1 && chargereq8  && !leptons) //0 and 1	 2 and 3
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 	 }
	 	 if(chargereq1 && chargereq9  && !leptons) //0 and 1	 2 and 4
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 	 }
	 	 if(chargereq1 && chargereq10  && !leptons) //0 and 1	  3 and 4
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	 
	 	 }
	 	 if(chargereq2 && chargereq6  && !leptons) //0 and 2	 1 and 3
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 	 
	 	 }
	 	 if(chargereq2 && chargereq7  && !leptons) //0 and 2	 1 and 4
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 	 }
	 	 if(chargereq2 && chargereq10  && !leptons) //0 and 2	  3 and 4
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 	 }
	 	 if(chargereq3 && chargereq5  && !leptons) //0 and 3	 1 and 2
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 	 }
	 	 if(chargereq3 && chargereq7  && !leptons) //0 and 3	 1 and 4
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	 }
	 	 if(chargereq3 && chargereq9  && !leptons) //0 and 3	 4 and 2
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 	 }
	 	 if(chargereq4 && chargereq5  && !leptons) //0 and 4	 1 and 2
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 	 }
	 	 if(chargereq4 && chargereq6  && !leptons) //0 and 4	 1 and 3
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 	 }
	 	 if(chargereq4 && chargereq8  && !leptons) //0 and 4	 3 and 2
	 	 {
	 		leptons=true; 
	 	 
	 		lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	 		lepton_2.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	 		lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	 		lepton_1.SetPxPyPzE(_Electrons[4]->Px(),_Electrons[4]->Py(),_Electrons[4]->Pz(),_Electrons[4]->Energy());
	 		lepton_4.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	 	 }
	 	 
	 	 leptonpair_1 = lepton_0 + lepton_1; 
	 	 leptonpair_2 = lepton_2 + lepton_3;
	 	 leptonfour = leptonpair_1 + leptonpair_2;
	 	 missingEt_pZ_pair = MEzCalculator(lepton_4, _missingEt_pX, _missingEt_pY, isMuon);
	  }
	  
	  TLorentzVector missingEt_vector1; 
	   missingEt_vector1.SetPxPyPzE(_missingEt_pX, _missingEt_pY, missingEt_pZ_pair.first,_missingEt/TMath::Sin(_missingEt_Theta));
			
	  TLorentzVector missingEt_vector2; 
	  missingEt_vector2.SetPxPyPzE(_missingEt_pX, _missingEt_pY, missingEt_pZ_pair.second,_missingEt/TMath::Sin(_missingEt_Theta));
			
	  TLorentzVector boolVector; 
	  if(leptons) 
	  {
		boolVector.SetPxPyPzE(1,1,1,1);
		if(debug) cout << "leptons set to true" << endl; 
          }
	  else boolVector.SetPxPyPzE(0,0,0,0);
	  
	  
	  vector <TLorentzVector> TempVector; 
	  TempVector.push_back(leptonfour); 
	  TempVector.push_back(lepton_4);
	  
	  TempVector.push_back(missingEt_vector1); 
	  TempVector.push_back(missingEt_vector2); 
	  TempVector.push_back(boolVector); 
	  
	  return TempVector; 
}




//Leptons ordering according to Z decay for 45 channel
vector <TLorentzVector> Channel_45_4leptons( bool _debug, vector <TRootElectron*> _Electrons, vector <TRootMuon*> _Muons)
{
	
	double nMuons = _Muons.size(); 
	double nElectrons = _Electrons.size(); 
	         
	bool leptons = false; 
	bool chargereq1=false; 
	bool chargereq2=false;
	bool chargereq3=false;
	bool chargereq4=false;
	bool chargereq5=false;
	bool chargereq6=false;
	
	bool debug = _debug;

	TLorentzVector leptonpair_1; // two OSSF with highest pt
	TLorentzVector leptonpair_2; // two OSSF with 2nd highest pt
	TLorentzVector leptonfour; // decay of Z
	TLorentzVector lepton_0;
	TLorentzVector lepton_1;
	TLorentzVector lepton_2;
	TLorentzVector lepton_3;
	
	lepton_0.Clear(); 
	lepton_1.Clear();
	lepton_2.Clear();
	lepton_3.Clear();
	
	leptonpair_1.Clear();
	leptonpair_2.Clear();
	leptonfour.Clear();
	
 
	
	//select ZZ events, muon flavours come in pairs and OS
	if( nMuons == 2 & nElectrons == 2)
	{ 
	     if(_Muons[0]->charge()!=_Muons[1]->charge())   chargereq1 = true; 
	     if(_Electrons[0]->charge()!=_Electrons[1]->charge())   chargereq2 = true;
	
	     if(debug) cout << "In nMuons == 2 & nElectrons == 2 " << endl; 
	     if(chargereq1 && chargereq2  && !leptons)
	     {
		    leptons=true; 
	
		    lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
		    lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
		    lepton_2.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
		    lepton_3.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	
		    
	
	     }
	     
	     leptonpair_1 = lepton_0 + lepton_1; 
	     leptonpair_2 = lepton_2 + lepton_3;
	     leptonfour = leptonpair_1 + leptonpair_2;
	     
	     
	} 

	if( nMuons == 4)
	{
	   if(debug) cout << "In nMuons == 4 " << endl; 
	   if(_Muons[0]->charge()!=_Muons[1]->charge())	chargereq1 = true; 
	   if(_Muons[0]->charge()!=_Muons[2]->charge())	chargereq2 = true;
	   if(_Muons[0]->charge()!=_Muons[3]->charge())	chargereq3 = true;
	   if(_Muons[2]->charge()!=_Muons[1]->charge())	chargereq4 = true;
	   if(_Muons[3]->charge()!=_Muons[1]->charge())	chargereq5 = true;
	   if(_Muons[3]->charge()!=_Muons[2]->charge())	chargereq6 = true;
	   
	   if(chargereq1 && chargereq6  && !leptons) //0 and 1     2 and 3
	   {
	   	  leptons=true; 
	
	   	  lepton_3.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	   	  lepton_2.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	   	  lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	   	  lepton_1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	   }
	   if(chargereq2 && chargereq5  && !leptons) //0 and 2     1 and 3
	   {
	          leptons=true; 
	   
	          lepton_3.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	          lepton_1.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	          lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	          lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	   }
	   if(chargereq3 && chargereq4  && !leptons) //0 and 3     2 and 1
	   {
	          leptons=true; 
	   
	          lepton_1.SetPxPyPzE(_Muons[3]->Px(),_Muons[3]->Py(),_Muons[3]->Pz(),_Muons[3]->Energy());
	          lepton_3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(),_Muons[2]->Energy());
	          lepton_0.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(),_Muons[0]->Energy());
	          lepton_2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(),_Muons[1]->Energy());
	   }
	   
	   
	   leptonpair_1 = lepton_0 + lepton_1; 
	   leptonpair_2 = lepton_2 + lepton_3;
	   leptonfour = leptonpair_1 + leptonpair_2;
	   	    
	     
	
	   }
	   if(nElectrons == 4)
	   {
	   	   if(debug) cout << "In nElectrons == 4 " << endl; 
	   	   if(_Electrons[0]->charge()!=_Electrons[1]->charge())	chargereq1 = true; 
	   	   if(_Electrons[0]->charge()!=_Electrons[2]->charge())	chargereq2 = true;
	   	   if(_Electrons[0]->charge()!=_Electrons[3]->charge())	chargereq3 = true;
	   	   if(_Electrons[2]->charge()!=_Electrons[1]->charge())	chargereq4 = true;
	   	   if(_Electrons[3]->charge()!=_Electrons[1]->charge())	chargereq5 = true;
	   	   if(_Electrons[3]->charge()!=_Electrons[2]->charge())	chargereq6 = true;
	   	   
	   	   if(chargereq1 && chargereq6  && !leptons) //0 and 1     2 and 3
	   	   {
	   		  leptons=true; 
	
	   		  lepton_3.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	   		  lepton_2.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	   		  lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	   		  lepton_1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	   	   }
	   	   if(chargereq2 && chargereq5  && !leptons) //0 and 2     1 and 3
	   	   {
	   		  leptons=true; 
	   		  lepton_3.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	   		  lepton_1.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	   		  lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	   		  lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	   	   }
	   	   if(chargereq3 && chargereq4  && !leptons) //0 and 3     2 and 1
	   	   {
	   		  leptons=true; 
	
	   		  lepton_1.SetPxPyPzE(_Electrons[3]->Px(),_Electrons[3]->Py(),_Electrons[3]->Pz(),_Electrons[3]->Energy());
	   		  lepton_3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(),_Electrons[2]->Energy());
	   		  lepton_0.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(),_Electrons[0]->Energy());
	   		  lepton_2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(),_Electrons[1]->Energy());
	   	   }
	   	   
	   	   
	   	   leptonpair_1 = lepton_0 + lepton_1; 
	   	   leptonpair_2 = lepton_2 + lepton_3;
	   	   leptonfour = leptonpair_1 + leptonpair_2;
	
	   }
		 
	  TLorentzVector boolVector; 
	  if(leptons)
	  {
		 boolVector.SetPxPyPzE(1,1,1,1);
		if(debug) cout << "leptons set to true"<< endl; 
	  }
	  else{ boolVector.SetPxPyPzE(0,0,0,0);}
	  
     	  vector <TLorentzVector> TempVector; 
	  TempVector.clear();
	  TempVector.push_back(leptonfour); 
	  TempVector.push_back(boolVector); 
	  
	  return TempVector; 




}

//Leptons ordering according to Z decay for 3L channel
vector <TLorentzVector> Channel_3L_Zcandidate(bool _debug, vector <TRootElectron*> _Electrons, vector <TRootMuon*> _Muons)
{
	 
	bool OSSF_3L = false; 
	int Nb_OSSF_3L =0;
	bool debug = _debug; 
	
	TLorentzVector leptonpair_3L_mll;
	TLorentzVector lepton_3L1;
	TLorentzVector lepton_3L2;
	TLorentzVector lepton_3L3;
	
	leptonpair_3L_mll.Clear(); 
	lepton_3L1.Clear();
	lepton_3L2.Clear();
	
	
	 if(_Electrons.size() == 3)
	 {
	 	 if(debug) cout << "_Electrons[0]->charge() " << _Electrons[0]->charge() << endl;  
	 	 if(_Electrons[0]->charge()!=_Electrons[1]->charge())
	 	 {
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(), _Electrons[0]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(), _Electrons[1]->Energy());
			 lepton_3L3.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(), _Electrons[2]->Energy());
	 	 
	 	 }
	 	 else if(_Electrons[2]->charge()!=_Electrons[0]->charge()) 
	 	 {
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(), _Electrons[0]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(), _Electrons[2]->Energy());
			 lepton_3L3.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(), _Electrons[1]->Energy());
	 	 
	 	 }
	 	 else if(_Electrons[1]->charge()!=_Electrons[2]->charge())
	 	 { 
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(), _Electrons[1]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Electrons[2]->Px(),_Electrons[2]->Py(),_Electrons[2]->Pz(), _Electrons[2]->Energy());
			 lepton_3L3.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(), _Electrons[0]->Energy());
	 	 
	 	 }
	 	 
	 }
	 if(_Electrons.size() ==2)
	 {
	 	 if(debug) cout << "_Electrons[0]->charge() " << _Electrons[0]->charge() << endl; 
	 	 if(_Electrons[0]->charge()!=_Electrons[1]->charge())
	 	 {
	 		  OSSF_3L = true; 
	 		  Nb_OSSF_3L =1;
	 		  lepton_3L1.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(), _Electrons[0]->Energy());
	 		  lepton_3L2.SetPxPyPzE(_Electrons[1]->Px(),_Electrons[1]->Py(),_Electrons[1]->Pz(), _Electrons[1]->Energy());
	 	 
	 	 }
	 	 
	 	 lepton_3L3.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(), _Muons[0]->Energy());
	 }
	 if(_Electrons.size() == 1)
	 {
	 	 if(debug) cout << "_Electrons[0]->charge() " << _Electrons[0]->charge() << endl; 
	 	 if(_Muons[0]->charge()!=_Muons[1]->charge())
	 	 {
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(), _Muons[0]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(), _Muons[1]->Energy());
	 	 
	 	 }
	 	 lepton_3L3.SetPxPyPzE(_Electrons[0]->Px(),_Electrons[0]->Py(),_Electrons[0]->Pz(), _Electrons[0]->Energy());
	 	
	 }
	 if(_Electrons.size() ==  0)
	 {
	 	 if(debug) cout << "_Muons[0]->charge() " << _Muons[0]->charge() << endl; 
	 	 if(_Muons[0]->charge()!=_Muons[1]->charge())
	 	 {
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(), _Muons[0]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(), _Muons[1]->Energy());
	 	 	 lepton_3L3.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(), _Muons[2]->Energy());
	 	 }
	 	 else if(_Muons[0]->charge()!=_Muons[2]->charge())
	 	 {
	 		  OSSF_3L = true; 
	 		  Nb_OSSF_3L =1;
	 		  lepton_3L1.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(), _Muons[0]->Energy());
	 		  lepton_3L2.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(), _Muons[2]->Energy());
			  lepton_3L3.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(), _Muons[1]->Energy());
	 	 	 
	 	 
	 	 }
	 	 else if(_Muons[1]->charge()!=_Muons[2]->charge())
	 	 {
	 		 OSSF_3L = true; 
	 		 Nb_OSSF_3L =1;
	 		 lepton_3L1.SetPxPyPzE(_Muons[1]->Px(),_Muons[1]->Py(),_Muons[1]->Pz(), _Muons[1]->Energy());
	 		 lepton_3L2.SetPxPyPzE(_Muons[2]->Px(),_Muons[2]->Py(),_Muons[2]->Pz(), _Muons[2]->Energy());
	 	 	 lepton_3L3.SetPxPyPzE(_Muons[0]->Px(),_Muons[0]->Py(),_Muons[0]->Pz(), _Muons[0]->Energy());
	 		  
	 	 }
	 }	 
	 
	   
	 TLorentzVector boolVector; 
	  if(OSSF_3L)
	  {
		 boolVector.SetPxPyPzE(1,1,1,1);
		if(debug) cout << "OSSF_3L set to true"<< endl; 
	  }
	  else{ boolVector.SetPxPyPzE(0,0,0,0);}
	 
	 vector <TLorentzVector> TempVector; 
	 TempVector.clear();
	 TempVector.push_back(lepton_3L1); 
	 TempVector.push_back(lepton_3L2); 
	 TempVector.push_back(lepton_3L3); 
	 TempVector.push_back(boolVector);
	 
	 return TempVector; 
        

}

vector<pair<int,bool> > Channel_45_FCNC_cjet(bool _debug, vector<pair<int,bool> > _LightJets_Paired,vector<TRootJet*> _LightJets, TLorentzVector _leptonfour)
{
	bool debug = _debug; 
	if(debug) cout << "in FCNC top reconstruction" << endl; 
	
	
	Double_t DeltaR_LightJet_Higgs = 1000; 
	Double_t DeltaPhi_LightJet_Higgs = 1000; 
	Double_t Pt = 0.;
	double tag = 0; 
		
	//take the jet closest to the higgs in delta R	
	for( int i = 0; i <_LightJets.size(); i++)
	 {
	        
	 	 if(debug) cout << "in for loop" << endl;
	 	
	 	 TLorentzVector LightJet;
	 	 LightJet.SetPxPyPzE(_LightJets[i]->Px(),_LightJets[i]->Py(),_LightJets[i]->Pz(), _LightJets[i]->Energy());
	 	 
		 
		 
		 if(debug) cout << "defined lightjet" << endl; 
		 
	 	 Double_t Delta_R = LightJet.DeltaR(_leptonfour); 
		 //Double_t Delta_R = LightJet.DrEtaPhi(_leptonfour); 
	 	 if(debug) cout << "DeltaR: " << Delta_R<< endl;
		 
		 Double_t Delta_Phi = LightJet.DeltaPhi(_leptonfour); 
	 	 
		 Double_t tempPt = LightJet.Pt(); 
		 
		 double temptag = _LightJets[i]->btag_combinedSecondaryVertexBJetTags(); 
		 
		 /*
		 if(Delta_Phi < DeltaPhi_LightJet_Higgs)  // has 34% selection efficiency
	 	 {
	 		 DeltaPhi_LightJet_Higgs = Delta_Phi; 
			 _LightJets_Paired[i].second = true;
			 
			 
			 if(i != 0)
			 {
			    for(int k = 0; k<i; k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
			    for(int k = i+1; k<_LightJets.size(); k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
	 		 }
			 if(debug) cout << "set aPair.second to true" << endl; 
			 
			
	 	 }
		 */
		 /*
		 if(temptag > tag)  // 37% selection efficiency 
		 {
		 	tag = temptag;
			_LightJets_Paired[i].second = true;
			 
			 
			 if(i != 0)
			 {
			    for(int k = 0; k<i; k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
			    for(int k = i+1; k<_LightJets.size(); k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
	 		 }
			 if(debug) cout << "set aPair.second to true" << endl; 
			
		 }
		 */
		 
	 	 if(Delta_R < DeltaR_LightJet_Higgs)  // has 52% selection efficiency
	 	 {
	 		 DeltaR_LightJet_Higgs = Delta_R; 
			 _LightJets_Paired[i].second = true;
			 
			 
			 if(i != 0)
			 {
			    for(int k = 0; k<i; k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
			    for(int k = i+1; k<_LightJets.size(); k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
	 		 }
			 if(debug) cout << "set aPair.second to true" << endl; 
			 
			
	 	 }
		 
		 /*
		 if(tempPt > Pt)   // has 29% selection efficiency
		 {
		        Pt = tempPt; 
		 	_LightJets_Paired[i].second = true;
			 
			 
			 if(i != 0)
			 {
			    for(int k = 0; k<i; k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
			    for(int k = i+1; k<_LightJets.size(); k++)
			    {
			       _LightJets_Paired[k].second = false; 
			    }
	 		 }
			 if(debug) cout << "set aPair.second to true" << endl; 
		 
		 }
		 */
	 	 
	   
	 }
	
	 if(debug) cout << "out Ljets > 0" << endl;
	 
	/*
	for(int ii = 0; ii< _LightJets_Paired.size(); ii++)
	{
	   pair<int,bool> qaPair = _LightJets_Paired[ii]; 
	   cout << qaPair.first << " " << qaPair.second << endl;
			
	}

        cout << "return" << endl;
	*/
	return _LightJets_Paired; 
}

TLorentzVector Channel_45_FCNC_top(bool _debug, vector<pair<int,bool> > _Pairs, TLorentzVector _leptonfour, vector<TRootJet*> _LightJets)
{
	double InvMass = 0; 
	bool debug = _debug; 
	if(debug) cout << "in top fcnc " << endl; 
	int counter = 0; 
	TLorentzVector temp;
	
	for( int i = 0; i <_LightJets.size(); i++)
	 {
	 	 if(debug) cout << "in for loop" << endl;
	 	 pair<int,int> aPair = _Pairs[i];
	 	 
		 TLorentzVector LightJet;
	 	 LightJet.SetPxPyPzE(_LightJets[i]->Px(),_LightJets[i]->Py(),_LightJets[i]->Pz(), _LightJets[i]->Energy());
	 	 
		// cout << aPair.second << endl; 
		 if(aPair.second)
		 {
		 	temp.Clear();  
	             	temp = _leptonfour + LightJet; 
			InvMass = temp.M(); 
			counter++; 
	 	 }
		 if(counter>1) cout << "WARNING: something wrong with c-jet selection" << endl; 
	   
	 }
	return temp; 
}


// take b -jet with highest dicriminating value
vector <TLorentzVector> SM_b(bool _debug, vector<TRootJet*> _BJets)
{
	if(_debug) cout << "in SM bjet detemination " << endl; 
	vector <TLorentzVector> tempV; 
	TLorentzVector TempJet; 
	TempJet.Clear(); 
	double tag = 0; 
	int ID = -1; 
	for(int i = 0; i < _BJets.size(); i++)
	{
		
	        double temptag = _BJets[i]->btag_combinedSecondaryVertexBJetTags(); 
		if(temptag > tag)
		{
		  tag = temptag; 
		  TempJet.SetPxPyPzE(_BJets[i]->Px(),_BJets[i]->Py(),_BJets[i]->Pz(), _BJets[i]->Energy());
		  ID = i; 
		}
	
	}
	
        tempV.push_back(TempJet); 
	
	TLorentzVector tagV; 
	tagV.SetPxPyPzE(tag,0,0,0);
	tempV.push_back(tagV); 
        
	TLorentzVector tempID; 
	tempID.SetPxPyPzE(ID,0,0,0); 
	tempV.push_back(tempID); 
	if(_debug) cout << "out SM bjet detemination " << endl;
	return tempV; 
}

//find jets closest to the bjet with highest discr
vector <TLorentzVector> Channel_45_SM_Wqq(bool _debug,TLorentzVector _Bjet, vector <pair<int,bool> > _Paired, vector <TRootJet*> _LightJets)
{
    vector <TLorentzVector> tempV; 
    tempV.clear();
    double dR = 10000; 
    double Mass = 10000; 
    int JetID_1 = -1; 
    int JetID_2 = -1; 
    TLorentzVector jets;
    TLorentzVector Rjets; 
    jets.Clear(); 
    Rjets.Clear();
    	     
    	     
    
    for(int iJet = 0; iJet < _LightJets.size()-1; iJet++)
    {
    	    pair<int,int> aPair = _Paired[iJet];
	    
	    TLorentzVector LightJet;
    	    LightJet.SetPxPyPzE(_LightJets[iJet]->Px(),_LightJets[iJet]->Py(),_LightJets[iJet]->Pz(),_LightJets[iJet]->Energy());
    	    
	    
	    for(int i = 1; i < _LightJets.size(); i++)
	    {    
	        pair<int,int> bPair = _Paired[i];
		TLorentzVector LightJet2;
    	   	LightJet2.SetPxPyPzE(_LightJets[i]->Px(),_LightJets[i]->Py(),_LightJets[i]->Pz(),_LightJets[i]->Energy());
    	    
    	    	if(!aPair.second && !bPair.second )
    	    	{	    
    		    jets =  LightJet + LightJet2;
		    if(fabs(jets.M() - 80.3) < fabs(Mass - 80.3) )
		    {
		    	Mass = jets.M(); 
		    	double tempDeltaR = _Bjet.DeltaR(jets);
    		    	if(tempDeltaR < dR)
		    	{
		    		dR = tempDeltaR;
				Rjets = jets; 
				JetID_1 = aPair.first;  // place of the jet in selectedJets
				JetID_2 = bPair.first; 
			 
    		    	}
		    }

    	    	}
	   }
    }
    
    /*
    TLorentzVector firstJet; 
    int lastSet = -1; 
    for(int iJet = 0; iJet < _LightJets.size(); iJet++)
    {
    	    cout << "iJet : " << iJet << endl; 
	    pair<int,int> aPair = _Paired[iJet];
	    
	    TLorentzVector LightJet;
    	    LightJet.SetPxPyPzE(_LightJets[iJet]->Px(),_LightJets[iJet]->Py(),_LightJets[iJet]->Pz(),_LightJets[iJet]->Energy());
    	    
	   if(aPair.second) cout << "ID of true: " << iJet; 
	    
    	    if(!aPair.second)
    	    {		
    	    	
	    	double tempDeltaR = _Bjet.DeltaR(LightJet);
    	    	if(tempDeltaR < dR)
	    	{
			if(iJet!=0) 
			{
				//cout << "Pair with id " << lastSet << " set to false" << endl; 
				
				_Paired[lastSet].second = false; 
			}
			
	    	    	dR = tempDeltaR;
	    	       // Rjets = jets; 
	    	    	JetID_1 = aPair.first;  // place of the jet in selectedJets
	    	       // JetID_2 = bPair.first; 
	    	 	firstJet = LightJet; 
			_Paired[iJet].second = true; 
			lastSet = iJet; 
			//cout << "lastSet changed to " << lastSet << endl; 
    	    	}
	    	

    	    }
	   
    }
     
    int counter = 0; 
    for (int i = 0; i<_Paired.size(); i++)
    {
    	if(_Paired[i].second) counter++; 
    
    }
    //cout << "Nb set to true: " << counter << endl; 
    //cout << endl;
    //cout << endl;
    dR = 10000;
    TLorentzVector secondJet; 
    for(int iJet = 0; iJet < _LightJets.size(); iJet++)
    {
    	    pair<int,int> aPair = _Paired[iJet];
	    
	    TLorentzVector LightJet;
    	    LightJet.SetPxPyPzE(_LightJets[iJet]->Px(),_LightJets[iJet]->Py(),_LightJets[iJet]->Pz(),_LightJets[iJet]->Energy());
    	    
	    
	    
    	    if(!aPair.second)
    	    {		
    	    	
	    	double tempDeltaR = _Bjet.DeltaR(LightJet);
    	    	if(tempDeltaR < dR)
	    	{
	    	    	dR = tempDeltaR;
	    	       
	    	    	//JetID_1 = aPair.first;  // place of the jet in selectedJets
	    	        JetID_2 = aPair.first; 
	    	 	secondJet = LightJet; 
    	    	}
	    	

    	    }
	   
    }
    
    Rjets = firstJet + secondJet; 
    */
    TLorentzVector ID_1; 
    ID_1.SetPxPyPzE(JetID_1,0,0,0); 
    
    TLorentzVector ID_2; 
    ID_2.SetPxPyPzE(JetID_2,0,0,0); 
    
    
    tempV.push_back(Rjets); 
    tempV.push_back(ID_1); 
    tempV.push_back(ID_2); 
    
    return tempV; 
}

//Leptons ordering according to Z decay for 3L channel
double Channel_3L_Zcandidate_Mll(bool _debug, vector <TLorentzVector> _leptons)
{
	if(_debug)  cout << "In Channel_3L_Zcandidate_Mll " << endl; 
	double mll_3L = 0;
	TLorentzVector leptonpair_3L_mll;
	TLorentzVector lepton_3L1 = _leptons[0];
	TLorentzVector lepton_3L2 = _leptons[1];
	
	TLorentzVector Sum = lepton_3L1 + lepton_3L2; 
	mll_3L = Sum.M(); 
	 
	if(_debug)  cout << "Out Channel_3L_Zcandidate_Mll " << endl; 
	
	return mll_3L; 
        

}

TLorentzVector Channel_3L_FCNC_top_candidate(bool _debug, TLorentzVector _Z_candidate, vector<TRootJet*> _LightJets)
{
	double DeltaR = 10000; 
	TLorentzVector candidate; 
	for(int i = 0; i< _LightJets.size(); i++)
	{
		TLorentzVector tempJet; 
		tempJet.SetPxPyPzE(_LightJets[i]->Px(), _LightJets[i]->Py(), _LightJets[i]->Pz(), _LightJets[i]->Energy()); 
		if(tempJet.DeltaR(_Z_candidate) < DeltaR)
		{
			DeltaR = tempJet.DeltaR(_Z_candidate); 
			candidate = tempJet + _Z_candidate; 
			if(_debug) cout << "DeltaR = " << DeltaR << endl; 
		
		}
	
	
	}



	return candidate; 
}

//calculate Higgs for tcH(ZZ(llqq))
TLorentzVector Channel_3L_Higgs_candidate(bool _debug, TLorentzVector _Z_candidate,  vector<TRootJet*> _LightJets)
{
	double DeltaR = 10000; 
	TLorentzVector candidate; 
	for(int i = 0; i< _LightJets.size()-1; i++)
	{
		TLorentzVector tempJet; 
		tempJet.SetPxPyPzE(_LightJets[i]->Px(), _LightJets[i]->Py(), _LightJets[i]->Pz(), _LightJets[i]->Energy()); 
		
		for(int k = 1; k < _LightJets.size(); k++)
		{
			TLorentzVector tempJet2; 
			tempJet2.SetPxPyPzE(_LightJets[k]->Px(), _LightJets[k]->Py(), _LightJets[k]->Pz(), _LightJets[k]->Energy()); 
			
			TLorentzVector combi; 
			combi = tempJet + tempJet2; 
			
			if(DeltaR > combi.DeltaR(_Z_candidate))
			{
				DeltaR = combi.DeltaR(_Z_candidate); 
				candidate = combi + _Z_candidate; 
			
			}
		
		}
	}
	
	return candidate; 
}

//calculate FCNC top for tcH(ZZ(llqq))
TLorentzVector Channel_3L_FCNC_top_candidate2(bool _debug, TLorentzVector _H_candidate,  vector<TRootJet*> _LightJets)
{
	double DeltaR = 10000; 
	TLorentzVector candidate; 
	for(int i = 0; i< _LightJets.size(); i++)
	{
		TLorentzVector tempJet; 
		tempJet.SetPxPyPzE(_LightJets[i]->Px(), _LightJets[i]->Py(), _LightJets[i]->Pz(), _LightJets[i]->Energy()); 
		
		if(DeltaR > tempJet.DeltaR(_H_candidate))
		{
			DeltaR = tempJet.DeltaR(_H_candidate); 
			candidate = tempJet + _H_candidate; 
		
		}
		
	}
	
	return candidate; 
}
/*
vector <pair<TLorentzVector,bool> > Channel_3L_SM_lep(bool _debug,vector<TRootMuon*> _Muons, vector<TRootElectron*> _Electrons,TLorentzVector Bjet)
{
	
	double deltaR = 10000; 
	vector <pair<TLorentzVector,bool> > candidate; 
	
	for(int i = 0; i< _Muons.size(); i++)
	{
		TLorentzVector tempV;
		tempV.SetPxPyPzE(_Muons[i]->Px(), _Muons[i]->Py(), _Muons[i]->Pz(), _Muons[i]->Energy());
		pair <TLorentzVector,bool> Pair = make_pair<tempV,false>;
		candidate.push_back(Pair); 
		
		
	}
	
	for(int i = 0; i< _Electrons.size(); i++)
	{
		TLorentzVector tempV;
		tempV.SetPxPyPzE(_Electrons[i]->Px(), _Electrons[i]->Py(), _Electrons[i]->Pz(), _Electrons[i]->Energy());
		pair <TLorentzVector,bool> Pair = make_pair<tempV,false>;
		candidate.push_back(Pair); 
	
	}

	
	for(int i = 0; i< _Muons.size(); i++)
	{
		
		TLorentzVector tempV = candidate[i].first;
		
		if(Bjet.DeltaR(tempV) < deltaR)
		{
			deltaR = Bjet.DeltaR(tempV); 
			candidate[i].second = true;
			if(i>0)  candidate[i-1].second = false;
		}
	
	}
	
	for(int i = _Muons.size(); i< _Electrons.size()+_Muons.size(); i++)
	{
		TLorentzVector tempV = candidate[i].first;
		
		if(Bjet.DeltaR(tempV) < deltaR)
		{
			deltaR = Bjet.DeltaR(tempV); 
			candidate[i].second = true;
			if(i>0) candidate[i-1].second = false;
		}
	
	}

	return candidate;
}
*/
TLorentzVector Channel_3L_Higgs_lep(bool _debug,vector<TRootMuon*> _Muons, vector<TRootElectron*> _Electrons,vector<pair<TLorentzVector,bool> > _Wleptons)
{
	TLorentzVector candidate; 
	candidate.Clear(); 
	 
	for(int i = 0; i < _Wleptons.size();i++)
	{
		pair<TLorentzVector,bool> Pair = _Wleptons[i]; 
		if (!Pair.second)
		{
			candidate = candidate + Pair.first; 
		
		}
	
	}

	return candidate; 
}		    	    		
			

int main (int argc, char *argv[])
{
    
    
    clock_t start = clock();
    
         
    ///////////////////////////////////////////////////////////
    //   different options for executing this macro	     //
    ///////////////////////////////////////////////////////////

    std::string tempxmlName;
    std::string channelName = "undefined";
    bool information = true; 
    bool warnings = true; 
    bool debug = false; 
    int nb_Leptons = 0; 
    bool foundxml= false;
    bool validate = false; 

    for(int iarg = 0; iarg < argc && argc>1 ; iarg++)
    {
	    std::string argval=argv[iarg];
    	    
	    if(argval=="--help" || argval =="--h")
    	    {
    		    cout << "--NoWarnings: put warnings off " << endl; 
    		    cout << "--NoInfo: put information off " << endl; 
    		    cout << "--debug: put debug output on" << endl; 
    		    cout << "--xml myxml.xml: change Xml file" << endl;
    		    cout << " " << endl; 
    		    cout << "--1L3B: use the 1 lepton + 3 b-tags channel" << endl; 
    		    cout << "--SSdilepton: use the same sign dilepton channel" << endl;
    		    cout << "--OSdilepton: use the opposite sign dilepton channel" << endl;
    		    cout << "--3L: use the 3 lepton channel (exactly 3)" << endl;
    		    cout << "--45: at least 4 leptons " << endl; 
    		    
    		    return 0;
	    }
	    if (argval=="--xml") {
    		    iarg++;
    		    tempxmlName = argv[iarg];
		    foundxml = true; 
		    if(debug) cout << "foundxml = true" << endl; 
    		    
    	    }
    	    if (argval=="--NoInfo") {
    		    iarg++;
    		    information = false;  
    	    }
    	    if (argval=="--NoWarnings") {
    		    iarg++;
    		    warnings = false;  
    	    }
    	    if (argval=="--debug") {
    		    iarg++;
    		    debug = true;  
    	    }
    	    		    
    	    if (argval=="--1L3B") {
    		    nb_Leptons = 1; 
		    channelName = "1L3B";
    		    if(!foundxml) tempxmlName = "../config/FCNC_1L3B_config.xml";
	    }
    	    if (argval=="--SSdilepton") {
    		    nb_Leptons = 2; 
		    channelName = "SSdilepton";
    		    if(!foundxml) tempxmlName = "../config/FCNC_SSdilepton_config.xml";
	    }
    	    if (argval=="--OSdilepton") {
    		    nb_Leptons = 2; 
		    channelName = "OSdilepton";
    		    if(!foundxml) tempxmlName = "../config/FCNC_OSdilepton_config.xml";
	    }
    	    if (argval=="--3L") {
    		    nb_Leptons = 3;
		    channelName = "3L";
    		    if(!foundxml) tempxmlName = "../config/FCNC_3L_config.xml";
	    }
    	    if (argval=="--45") {
    		    nb_Leptons = 3;
		    channelName = "45";
    		    if(!foundxml) tempxmlName = "../config/FCNC_45_config.xml";
	    }
	    if(argval == "--validate"){
	    	validate = true; 
	    
	    }
	    
    	    
    	    

    } 
    //put in a warning 
    if(channelName.find("undefined")!=string::npos && warnings) std::cout << "[WARNING]      No channel was defined" << endl; 
    if(nb_Leptons == 0 && warnings) std::cout << "[WARNING]      No nb of leptons was defined, default setting is 0" << endl; 





    //SetStyle if needed
    //setTDRStyle();
    setMyStyle();
    
    /////////////////////
    // Configuration
    /////////////////////
    
    //xml file
    string xmlFileName = tempxmlName;
    
 
    const char *xmlfile = xmlFileName.c_str();
    const char *channel = channelName.c_str();
    if(information)
    {
    cout << "********************************************************" << endl;
    cout << "used config file: " << xmlfile << endl;
    cout << "used channel: " << channel << endl; 
    cout << "********************************************************" << endl;
    }
    //Configuration output format
    TTree *configTree = new TTree("configTree","configuration Tree");
    TClonesArray* tcdatasets = new TClonesArray("Dataset",1000);
    configTree->Branch("Datasets","TClonesArray",&tcdatasets);
    TClonesArray* tcAnaEnv = new TClonesArray("AnalysisEnvironment",1000);
    configTree->Branch("AnaEnv","TClonesArray",&tcAnaEnv);
    
    ////////////////////////////////////
    /// AnalysisEnvironment
    ////////////////////////////////////
    
    AnalysisEnvironment anaEnv;
    if(debug) std::cout << "Loading the analysisenvironment" << endl; 
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    if(debug) std::cout << "done creating AnalysisEnvironmentLoader" << endl; 

    new ((*tcAnaEnv)[0]) AnalysisEnvironment(anaEnv);
    int verbose = anaEnv.Verbose;
    float oldLuminosity = anaEnv.Luminosity;	// in 1/pb
    
    if(debug)   cout << "analysis environment luminosity for rescaling "<< oldLuminosity << endl;
    
    /////////////////////
    // Load Datasets
    /////////////////////
    
    TTreeLoader treeLoader;
    if(debug)    cout << " - Load datasets ..." << endl;
    vector < Dataset* > datasets;
    
    treeLoader.LoadDatasets (datasets, xmlfile);
    for(unsigned int i=0;i<datasets.size();i++) new ((*tcdatasets)[i]) Dataset(*datasets[i]);
    
    float Luminosity = oldLuminosity;
    
    cout << "********************************************************" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        if(Luminosity > datasets[d]->EquivalentLumi() ) Luminosity = datasets[d]->EquivalentLumi();
        
        string dataSetName = datasets[d]->Name();
        cout << "datasets: " << dataSetName << endl;
    }
    cout << "********************************************************" << endl;
    
    
    
    //Global variable
    //TRootEvent* event = 0;
    
    //nof selected events
    double NEvtsData = 0;
    Double_t *nEvents = new Double_t[datasets.size()];
    Double_t *nEvents_Selected = new Double_t[datasets.size()];
    
    ////////////////////////
    // PileUp Reweighting //

    ////////////////////////
    
    //cout << Luminosity << endl;
    
    LumiReWeighting LumiWeights, LumiWeightsUp, LumiWeightsDown;
    
    LumiWeights = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root","../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/nominal.root", "pileup", "pileup");
    LumiWeightsUp = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_up.root", "pileup", "pileup");
    LumiWeightsDown = LumiReWeighting("../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_Summer12_S10.root", "../../TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2012Data53X_UpToRun208357/sys_down.root", "pileup", "pileup");
    
   
    if(debug)    cout << " Initialized LumiReWeighting stuff" << endl;
    
    ////////////////////////////////////
    //	Loop on datasets
    ////////////////////////////////////
    
    if(debug) cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        
        string previousFilename = "";
        int iFile = -1;
        string dataSetName = datasets[d]->Name();
        
        cout << "   Dataset " << d << ": " << datasets[d]->Name () << "/ title : " << datasets[d]->Title () << endl;
        if (debug)
            std::cout<<"      -> This sample contains, " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        
//      make root tree file name
        string roottreename = "../ntuples/";
	roottreename+= channelName; 
	roottreename+= "_"; 
	roottreename+= datasets[d]->Name();
        roottreename+="_tree.root";
	//        cout << "creating tree in file " << roottreename << endl;
        
        TFile *fileout = new TFile (roottreename.c_str(), "RECREATE");
        fileout->cd();
        /////////////////////
        // My tree         //
        ////////////////////
        Int_t nElectrons;
        Double_t pX_electron[10];
        Double_t pY_electron[10];
        Double_t pZ_electron[10];
        Double_t E_electron[10];
	Double_t pfIso_electron[10];
        Int_t charge_electron[10];
        
        Int_t nMuons;
        Double_t pX_muon[10];
        Double_t pY_muon[10];
        Double_t pZ_muon[10];
	Double_t E_muon[10];
	Double_t pfIso_muon[10];
        Int_t charge_muon[10];
        
        Int_t nJets;
        Double_t pX_jet[10];
        Double_t pY_jet[10];
        Double_t pZ_jet[10];
        Double_t E_jet[10];
	
	Int_t nBJets;
        Double_t pX_Bjet[10];
        Double_t pY_Bjet[10];
        Double_t pZ_Bjet[10];
        Double_t E_Bjet[10];
	
	Int_t nLJets;
        Double_t pX_Ljet[10];
        Double_t pY_Ljet[10];
        Double_t pZ_Ljet[10];
	Double_t Phi_Ljet[10];
	Double_t E_Ljet[10];

        
        Double_t missingEt;
	Double_t missingEt_Phi;
	Double_t missingEt_Theta;
	Double_t missingEt_pX;
	Double_t missingEt_pY;
	Double_t missingEt_pZ;
	

        //45 channel variables
	
	Double_t InvMass_4lept_Zdecay;
	Double_t InvMass_FCNC_top_Zdecay; 
	Double_t InvMass_SM_lb; 
	Double_t InvMass_SM_W_lv;
	Double_t InvMass_SM_W_qq;
	Double_t InvMass_SM_W;
	Double_t InvMass_SM_top_blv;
	Double_t InvMass_SM_top_bqq;
	Double_t InvMass_SM_top;
	Double_t TrMass_W; 
	Double_t TrMass_W_qq; 
	Double_t TrMass_W_lv; 
	
	Double_t Phi_Higgs; 
	Double_t Eta_Higgs; 
	
	Double_t FCNC_top_eta; 
	Double_t FCNC_top_phi; 
        Double_t FCNC_top_Pt;
	
	Double_t Bdiscr;
	
	//3L channel variables
	/*
	Double_t InvMass_Z; 
	Double_t InvMass_H;
	Double_t Z_candidate_pT;
	Double_t Z_candidate_Eta;
	Double_t H_candidate_pT;
	Double_t H_candidate_Eta;
	
	Double_t pT_FCNC_top_tcZ;
	Double_t Eta_FCNC_top_tcZ; 
	Double_t InvMass_FCNC_top_tcZ; 
	
	
	
	Double_t pT_FCNC_top_tcH_ZZ_llqq;
	Double_t Eta_FCNC_top_tcH_ZZ_llqq; 
	Double_t InvMass_FCNC_top_tcH_ZZ_llqq; 
	
	
	Double_t Bjet_Eta; 
	Double_t Bjet_Phi; 
	Double_t Bjet_Px; 
	Double_t Bjet_Pt;
	Double_t Bjet_Py;
	Double_t Bjet_Pz;
	
	Double_t FCNC_ll_Eta; 
	Double_t FCNC_ll_Phi; 
	Double_t FCNC_ll_Px; 
	Double_t FCNC_ll_Pt;
	Double_t FCNC_ll_Py;
	Double_t FCNC_ll_Pz;
	
	Double_t InvMass_FCNC_ll; 
	Double_t DeltaR_SMlb_FCNCll; 
	Double_t DeltaPhi_SMlb_FCNCll;
	*/
	Double_t InvMass_Leptons;
	Double_t Eta_SM_Lepton;
	Double_t Pt_SM_Lepton; 
	Double_t InvMass_SM_Lepton;
	Double_t Eta_SM_top;
	Double_t Pt_SM_top; 
	//Double_t InvMass_SM_top;
	Double_t pT_FCNC_top_candidate;
	Double_t Eta_FCNC_top_candidate; 
	Double_t InvMass_FCNC_top_candidate;
	Double_t InvMass_FCNC_Leptons;
	Double_t Theta_FCNC_Leptons;
	Double_t Eta_FCNC_Leptons;
	Double_t Phi_FCNC_Leptons;
	Double_t Pt_FCNC_Leptons;
	Double_t DeltaR; 
	Double_t DeltaRmin; 
	
	Int_t nEvents_Tree; 
        Int_t isdata;
        // various weights
        Double_t pu_weight;
        
        
        TTree* myTree = new TTree("tree","tree");
        myTree->Branch("isdata",&isdata,"isdata/I");
	myTree->Branch("nEvents_Tree",&nEvents_Tree,"nEvents_Tree/I");
        
        myTree->Branch("nElectrons",&nElectrons, "nElectrons/I");
        myTree->Branch("pX_electron",pX_electron,"pX_electron[nElectrons]/D");
        myTree->Branch("pY_electron",pY_electron,"pY_electron[nElectrons]/D");
        myTree->Branch("pZ_electron",pZ_electron,"pZ_electron[nElectrons]/D");
        myTree->Branch("E_electron",E_electron,"E_electron[nElectrons]/D");
        myTree->Branch("pfIso_electron",pfIso_electron,"pfIso_electron[nElectrons]/D");
        myTree->Branch("charge_electron",charge_electron,"charge_electron[nElectrons]/I");
        
        myTree->Branch("nMuons",&nMuons, "nMuons/I");
        myTree->Branch("pX_muon",pX_muon,"pX_muon[nMuons]/D");
        myTree->Branch("pY_muon",pY_muon,"pY_muon[nMuons]/D");
        myTree->Branch("pZ_muon",pZ_muon,"pZ_muon[nMuons]/D");
        myTree->Branch("E_muon",E_muon,"E_muon[nMuons]/D");
        myTree->Branch("pfIso_muon",pfIso_muon,"pfIso_muon[nMuons]/D");
        myTree->Branch("charge_muon",charge_muon,"charge_muon[nMuons]/I");
        
        myTree->Branch("nJets",&nJets, "nJets/I");
        myTree->Branch("pX_jet",pX_jet,"pX_jet[nJets]/D");
        myTree->Branch("pY_jet",pY_jet,"pY_jet[nJets]/D");
        myTree->Branch("pZ_jet",pZ_jet,"pZ_jet[nJets]/D");
        myTree->Branch("E_jet",E_jet,"E_jet[nJets]/D");
	
	myTree->Branch("nBJets",&nJets, "nBJets/I");
        myTree->Branch("pX_Bjet",pX_Bjet,"pX_Bjet[nBJets]/D");
        myTree->Branch("pY_Bjet",pY_Bjet,"pY_Bjet[nBJets]/D");
        myTree->Branch("pZ_Bjet",pZ_Bjet,"pZ_Bjet[nBJets]/D");
        myTree->Branch("E_Bjet",E_Bjet,"E_Bjet[nBJets]/D");
	
	myTree->Branch("nLJets",&nJets, "nLJets/I");
        myTree->Branch("pX_Ljet",pX_Ljet,"pX_Ljet[nLJets]/D");
        myTree->Branch("pY_Ljet",pY_Ljet,"pY_Ljet[nLJets]/D");
        myTree->Branch("pZ_Ljet",pZ_Ljet,"pZ_Ljet[nLJets]/D");
        myTree->Branch("E_Ljet",E_Ljet,"E_Ljet[nLJets]/D");
        
        myTree->Branch("missingEt",&missingEt,"missingEt/D");
	myTree->Branch("missingEt_Phi",&missingEt_Phi,"missingEt_Phi/D");
	myTree->Branch("missingEt_Theta",&missingEt_Theta,"missingEt_Theta/D");
	myTree->Branch("missingEt_pX",&missingEt_pX,"missingEt_pX/D");
	myTree->Branch("missingEt_pY",&missingEt_pY,"missingEt_pY/D");

        myTree->Branch("pu_weight",&pu_weight,"pu_weight/D");
	
	if(channelName.find("45")!=string::npos)
	{
		//45 channel variables
		myTree->Branch("InvMass_4lept_Zdecay",&InvMass_4lept_Zdecay,"InvMass_4lept_Zdecay/D");
		myTree->Branch("InvMass_FCNC_top_Zdecay",&InvMass_FCNC_top_Zdecay,"InvMass_FCNC_top_Zdecay/D");
		myTree->Branch("InvMass_SM_lb",&InvMass_SM_lb,"InvMass_SM_lb/D");
		myTree->Branch("InvMass_SM_W_lv",&InvMass_SM_W_lv,"InvMass_SM_W_lv/D");
		myTree->Branch("InvMass_SM_W_qq",&InvMass_SM_W_qq,"InvMass_SM_W_qq/D");
		myTree->Branch("InvMass_SM_W",&InvMass_SM_W,"InvMass_SM_W/D");
		myTree->Branch("InvMass_SM_top_blv",&InvMass_SM_top_blv,"InvMass_SM_top_blv/D");
		myTree->Branch("InvMass_SM_top_bqq",&InvMass_SM_top_bqq,"InvMass_SM_top_bqq/D");
		myTree->Branch("InvMass_SM_top",&InvMass_SM_top,"InvMass_SM_top/D");
		myTree->Branch("TrMass_W",&TrMass_W,"TrMass_W/D");
		myTree->Branch("TrMass_W_qq",&TrMass_W_qq,"TrMass_W_qq/D");
		myTree->Branch("TrMass_W_lv",&TrMass_W_lv,"TrMass_W_lv/D");

		myTree->Branch("Phi_Higgs",&Phi_Higgs,"Phi_Higgs/D"); 
		myTree->Branch("Eta_Higgs",&Eta_Higgs,"Eta_Higgs/D");
		
		myTree->Branch("FCNC_top_eta", &FCNC_top_eta, "FCNC_top_eta/D"); 
		myTree->Branch("FCNC_top_phi", &FCNC_top_phi, "FCNC_top_phi/D"); 
		myTree->Branch("FCNC_top_Pt", &FCNC_top_Pt, "FCNC_top_Pt/D");
		
		myTree->Branch("Bdiscr",&Bdiscr,"Bdiscr/D");
	} 
       	if(channelName.find("3L")!=string::npos)
	{
        	//3L channel variables
		/*
        	myTree->Branch("InvMass_Z",&InvMass_Z,"InvMass_Z/D"); 
		myTree->Branch("InvMass_H",&InvMass_H,"InvMass_H/D");
		myTree->Branch("Z_candidate_pT",&Z_candidate_pT,"Z_candidate_pT/D");
		myTree->Branch("Z_candidate_Eta",&Z_candidate_Eta,"Z_candidate_Eta/D");
		myTree->Branch("H_candidate_pT",&H_candidate_pT,"H_candidate_pT/D");
		myTree->Branch("H_candidate_Eta",&H_candidate_Eta,"H_candidate_Eta/D");
		myTree->Branch("pT_FCNC_top_tcZ",&pT_FCNC_top_tcZ,"pT_FCNC_top_tcZ/D");
		myTree->Branch("Eta_FCNC_top_tcZ",&Eta_FCNC_top_tcZ,"Eta_FCNC_top_tcZ/D");
		myTree->Branch("InvMass_FCNC_top_tcZ",&InvMass_FCNC_top_tcZ,"InvMass_FCNC_top_tcZ/D");
		
		myTree->Branch("pT_FCNC_top_tcH_ZZ_llqq",&pT_FCNC_top_tcH_ZZ_llqq,"pT_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("Eta_FCNC_top_tcH_ZZ_llqq",&Eta_FCNC_top_tcH_ZZ_llqq,"Eta_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("InvMass_FCNC_top_tcH_ZZ_llqq",&InvMass_FCNC_top_tcH_ZZ_llqq,"InvMass_FCNC_top_tcH_ZZ_llqq/D");
		myTree->Branch("Bjet_Eta",&Bjet_Eta, "Bjet_Eta/D"); 
		myTree->Branch("Bjet_Phi",&Bjet_Phi, "Bjet_Phi/D"); 
		myTree->Branch("Bjet_Px",&Bjet_Px, "Bjet_Px/D"); 
		myTree->Branch("Bjet_Pt",&Bjet_Pt, "Bjet_Pt/D");
		myTree->Branch("Bjet_Py",&Bjet_Py, "Bjet_Py/D");
		myTree->Branch("Bjet_Pz",&Bjet_Pz, "Bjet_Pz/D");
		myTree->Branch("FCNC_ll_Eta",&FCNC_ll_Eta, "FCNC_ll_Eta/D"); 
		myTree->Branch("FCNC_ll_Phi",&FCNC_ll_Phi, "FCNC_ll_Phi/D"); 
		myTree->Branch("FCNC_ll_Px",&FCNC_ll_Px, "FCNC_ll_Px/D"); 
		myTree->Branch("FCNC_ll_Pt",&FCNC_ll_Pt, "FCNC_ll_Pt/D");
		myTree->Branch("FCNC_ll_Py",&FCNC_ll_Py, "FCNC_ll_Py/D");
		myTree->Branch("FCNC_ll_Pz",&FCNC_ll_Pz, "FCNC_ll_Pz/D");
		myTree->Branch("InvMass_SM_lb",&InvMass_SM_lb, "InvMass_SM_lb/D"); 
		myTree->Branch("InvMass_FCNC_ll",&InvMass_FCNC_ll, "InvMass_FCNC_ll/D"); 
		myTree->Branch("DeltaR_SMlb_FCNCll",&DeltaR_SMlb_FCNCll, "DeltaR_SMlb_FCNCll/D"); 
		myTree->Branch("DeltaPhi_SMlb_FCNCll",&DeltaPhi_SMlb_FCNCll, "DeltaPhi_SMlb_FCNCll/D");
		myTree->Branch("Bdiscr",&Bdiscr,"Bdiscr/D");
		myTree->Branch("InvMass_Leptons", &InvMass_Leptons,"InvMass_Leptons/D"); 
		*/
		myTree->Branch("DeltaR", &DeltaR, "DeltaR/D"); 
		myTree->Branch("DeltaRmin", &DeltaRmin, "DeltaRmin/D"); 
		
		myTree->Branch("Theta_FCNC_Leptons", &Theta_FCNC_Leptons,"Theta_FCNC_Leptons/D"); 
		myTree->Branch("Phi_FCNC_Leptons", &Phi_FCNC_Leptons,"Phi_FCNC_Leptons/D"); 
		myTree->Branch("Eta_FCNC_Leptons", &Eta_FCNC_Leptons,"Eta_FCNC_Leptons/D"); 
		myTree->Branch("Pt_FCNC_Leptons", &Pt_FCNC_Leptons,"Pt_FCNC_Leptons/D"); 
		myTree->Branch("InvMass_FCNC_Leptons", &InvMass_FCNC_Leptons,"InvMass_FCNC_Leptons/D"); 
		
		myTree->Branch("Eta_SM_Lepton",&Eta_SM_Lepton, "Eta_SM_Lepton/D"); 
		myTree->Branch("Pt_SM_Lepton",&Pt_SM_Lepton,"Pt_SM_Lepton/D"); 
		myTree->Branch("InvMass_SM_Lepton",&InvMass_SM_Lepton,"InvMass_SM_Lepton/D"); 
		
		myTree->Branch("Eta_SM_top",&Eta_SM_top, "Eta_SM_top/D"); 
		myTree->Branch("Pt_SM_top",&Pt_SM_top,"Pt_SM_top/D"); 
		myTree->Branch("InvMass_SM_top",&InvMass_SM_top,"InvMass_SM_top/D"); 
		
		
		myTree->Branch("pT_FCNC_top_candidate",&pT_FCNC_top_candidate,"pT_FCNC_top_candidate/D");
		myTree->Branch("Eta_FCNC_top_candidate",&Eta_FCNC_top_candidate,"Eta_FCNC_top_candidate/D");
		myTree->Branch("InvMass_FCNC_top_candidate",&InvMass_FCNC_top_candidate,"InvMass_FCNC_top_candidate/D");
	}
       
	//        myTree->Print();

        TH1F * EventSummary = new TH1F("EventSummary","EventSummary",2,0,2);
        TH1F * Xsection = new TH1F("Xsection","Xsection",2,0,2);
	TH1F * SM_b_selection_efficiency = new TH1F("SM_b_selection_efficiency","SM_b_selection_efficiency",2,0,2);
	TH1F * SM_Wjet_selection_efficiency = new TH1F("SM_Wjet_selection_efficiency","SM_Wjet_selection_efficiency",2,0,2);
	TH1F * FCNC_c_selection_efficiency = new TH1F("FCNC_c_selection_efficiency","FCNC_c_selection_efficiency",2,0,2);
	TH1F * Selected_bjet_Pt = new TH1F("Selected_bjet_Pt","Selected_bjet_Pt",200,0,200); 
	TH1F * True_bjet_Pt = new TH1F("True_bjet_Pt","True_bjet_Pt",200,0,200); 
	TH1F * Selected_cjet_Pt = new TH1F("Selected_cjet_Pt","Selected_cjet_Pt",200,0,200); 
	TH1F * True_cjet_Pt = new TH1F("True_cjet_Pt","True_cjet_Pt",200,0,200); 
	
        //open files and load
	if(debug)       cout<<"LoadEvent"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);
	if(debug)       cout<<"LoadEvent"<<endl;
        
        
        
        /////////////////////////////////////
        /// Initialize JEC factors
        /////////////////////////////////////
   	    
        vector<JetCorrectorParameters> vCorrParam;
        
        /*JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L3Absolute_AK5PFchs.txt");
         JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L2Relative_AK5PFchs.txt");
         JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_MC_L1FastJet_AK5PFchs.txt");
         
         //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
         vCorrParam.push_back(*L1JetPar);
         vCorrParam.push_back(*L2JetPar);
         vCorrParam.push_back(*L3JetPar);
         vector<TRootJet*> selectedLightJets
         if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) { // DATA!
         JetCorrectorParameters *ResJetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Summer12_V3_DATA_L2L3Residual_AK5PFchs.txt");
         vCorrParam.push_back(*ResJetCorPar);
         }*/
        
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
        
        // true means redo also the L1
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);
        
        
        ////////////////////////////////////
        //	Loop on events
        ////////////////////////////////////
        double MatchedCounter = 0.; 
	double EventsToMatch = 0.;
	double MatchedCounter_c = 0.; 
	double EventsToMatch_c = 0.; 
	double MatchedCounter_W = 0.; 
	double EventsToMatch_W = 0.;
	nEvents[d] = 0;
	nEvents_Selected[d] = 0; 
        int itriggerSemiMu = -1,itriggerSemiEl = -1, previousRun = -1;
	
	// store number of events in ntuple 
        nEvents_Tree = datasets[d]->NofEvtsToRunOver(); 
	EventSummary->SetBinContent(1, datasets[d]->NofEvtsToRunOver());
	
	
	//define all cross sections
	string datasetName = datasets[d]->Name(); 
		
	
	if(datasetName.find("GluGluHiggs4lep")!=string::npos) {Xsection->SetBinContent(1,0.005109893);}
	if(datasetName.find("VBHiggs4lep")!=string::npos) {Xsection->SetBinContent(1,0.00414341777);}
	
	if(datasetName.find("WW_To2L2Nu")!=string::npos) {Xsection->SetBinContent(1,5.757);}
	if(datasetName.find("WZ_To2L2Q")!=string::npos) {Xsection->SetBinContent(1,2.267);}
	if(datasetName.find("WZ_To3LNu")!=string::npos) {Xsection->SetBinContent(1,1.087);}
	if(datasetName.find("ZZ_To2L2Nu")!=string::npos) {Xsection->SetBinContent(1,0.713);}
	if(datasetName.find("ZZ_To2L2Q")!=string::npos) {Xsection->SetBinContent(1,2.492);}
	if(datasetName.find("ZZ_To4L")!=string::npos) {Xsection->SetBinContent(1,0.18);}
	
	if(datasetName.find("TBZ_ToLL_4F")!=string::npos) {Xsection->SetBinContent(1,0.0114);}
	if(datasetName.find("ttH")!=string::npos) {Xsection->SetBinContent(1,0.1293);}
	if(datasetName.find("TTZ")!=string::npos) {Xsection->SetBinContent(1,0.172);}
	if(datasetName.find("TTW")!=string::npos) {Xsection->SetBinContent(1,0.2148);}
		
	if(datasetName.find("ST_T_s-ch")!=string::npos) {Xsection->SetBinContent(1,3.79);}
        if(datasetName.find("ST_TBar_s-ch")!=string::npos) {Xsection->SetBinContent(1,1.76);}
	if(datasetName.find("ST_T_tW-ch")!=string::npos) {Xsection->SetBinContent(1,11.1);}
        if(datasetName.find("ST_TBar_tW-ch")!=string::npos) {Xsection->SetBinContent(1,11.1);}
	
	if(datasetName.find("TT_SemiLeptMGDecays")!=string::npos) {Xsection->SetBinContent(1,110.26);}
	if(datasetName.find("TT_FullLeptMGDecays")!=string::npos) {Xsection->SetBinContent(1,26.42);}
		
	if(datasetName.find("Z_1Jets")!=string::npos) {Xsection->SetBinContent(1,671.83);}
	if(datasetName.find("Z_2Jets")!=string::npos) {Xsection->SetBinContent(1,216.76);}
	if(datasetName.find("Z_3Jets")!=string::npos) {Xsection->SetBinContent(1,61.20);}
	if(datasetName.find("Z_4Jets")!=string::npos) {Xsection->SetBinContent(1,27.59);}
		
	
	if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL_WToJets")!=string::npos) {Xsection->SetBinContent(1,0.090636);	}
	if(datasetName.find("TTJetsTocHbW_HToWW_WToLNuL")!=string::npos) {Xsection->SetBinContent(1,0.022659);	}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToBB_ZToLL")!=string::npos) {	Xsection->SetBinContent(1, 0.0005135);	}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToJetsUDC_ZToLL")!=string::npos) {Xsection->SetBinContent(1, 0.0018609);}
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL_ZToLL")!=string::npos) {	Xsection->SetBinContent(1, 0.00067929);	}		
	if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL")!=string::npos) {Xsection->SetBinContent(1, 0.00016516);	}
	if(datasetName.find("TTJetsTocZbW")!=string::npos) {Xsection->SetBinContent(1, 0.1575);	}
	
        
        for (unsigned int ievt = 0; ievt < datasets[d]->NofEvtsToRunOver(); ievt++) // event loop
            //for (unsigned int ievt = 0; ievt < 20000; ievt++)
        {
            
            vector < TRootVertex* > vertex;
            vector < TRootMuon* > init_muons;
            vector < TRootElectron* > init_electrons;
            vector < TRootJet* > init_jets_corrected;
            vector < TRootJet* > init_jets;
            vector < TRootMET* > mets;
            vector < TRootGenJet* > genjets;
	    vector<TRootMCParticle*> mcParticles; // to check mother, daughter,..
	    
	    
	    if(datasetName.find("TTJetsTocHbW_HToZZ_ZToLL")!=string::npos && channelName.find("45")!=string::npos && validate  )
	    {
	    	if(debug) cout << "Loading mc particles" << endl; 
	    	treeLoader.LoadMCEvent(ievt,0,0,mcParticles,false); 
		
	    	sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
	    	
	    }
	    if(validate && (datasetName.find("TTJetsTocHbW_HToZZ_ZToJ")!=string::npos ||datasetName.find("TTJetsTocHbW_HToZZ_ZToB")!=string::npos ||datasetName.find("TTJetsTocHbW_HToZZ_ZToNuL")!=string::npos || datasetName.find("TTJetsTocZ")!=string::npos || datasetName.find("TTJetsTocHbW_HToWW_WToNuL")!=string::npos)  && channelName.find("3L")!=string::npos  )
	    {
	    	if(debug) cout << "Loading mc particles" << endl; 
	    	treeLoader.LoadMCEvent(ievt,0,0,mcParticles,false); 
		
	    	sort(mcParticles.begin(),mcParticles.end(),HighestPt()); // HighestPt() is included from the Selection class
	    	
	    }
	    
	   
	    
            
            nEvents[d]++;
            
            if(ievt%1000 == 0)
                std::cout<<"Processing the "<<ievt<<"th event (" <<
		((double)ievt/(double)datasets[d]->NofEvtsToRunOver())*100 <<"%)" << " +> "<< (nEvents_Selected[d]/(double)datasets[d]->NofEvtsToRunOver())*100 << "% selected of the total events" << flush<<"\r";
            
            ////////////////
            // LOAD EVENT //
            ////////////////
            
            TRootEvent* event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets_corrected, mets);
            isdata=0;
            if(! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) ) {
                genjets = treeLoader.LoadGenJet(ievt,false);
                sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
            }
            else{
                isdata=1;
            }
            
            
	    
            /////////////////////////////////
            // DETERMINE EVENT SCALEFACTOR //
            /////////////////////////////////
            
            // scale factor for the event
            float scaleFactor = 1.;
            
            // PU reweighting
            
            double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );

            if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0)
                lumiWeight=1;
            
            // up syst -> lumiWeight = LumiWeightsUp.ITweight( (int)event->nTruePU() );
            // down syst -> lumiWeight = LumiWeightsDown.ITweight( (int)event->nTruePU() );
            
            pu_weight=lumiWeight;
            
            scaleFactor = scaleFactor*lumiWeight;
            
            ///////////////////
            // TRIGGER SETUP //
            ///////////////////
            
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename){
                previousFilename = currentFilename;
                iFile++;
               // cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << dataSetName << endl;
            }
            
            int currentRun = event->runId();
            
            if(previousRun != currentRun)
                previousRun = currentRun;
            
	    //triggering only for data, I only have events that either have at least 2 muons or 2 electrons so use double, 
	    // otherwise include EMU trigger as well
            /*bool triggered = false; 
	    int trigger1; 
	    int trigger2; 
	    if(isdata == 1)
	    {
	    	cout << "isdata with currentRun " << currentRun << " iFile " << iFile << endl; 
		trigger1 = treeLoader.iTrigger("HLT_Mu17_Mu8_v17",currentRun,iFile); //double muon
		cout << "trigger1" << endl; 
		trigger2 = treeLoader.iTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",currentRun,iFile); //double electron
		cout << "trigger2" << endl; 
		
		triggered = (treeLoader.EventTrigged(trigger1) || treeLoader.EventTrigged(trigger2));
		
	 	cout << "triggered is true" << endl; 
		
	    }
            else triggered = true; 
	    
	    
	    if(!triggered) continue; //if the events isn't triggered go to the next event
	    */
	    /////////////////////////////////////////////////////////////////////////////
            // JES SYSTEMATICS && SMEAR JET RESOLUTION TO MIMIC THE RESOLUTION IN DATA //
            /////////////////////////////////////////////////////////////////////////////
            
            if( ! (dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) )
                
                jetTools->correctJetJER(init_jets_corrected, genjets, mets[0], "nominal",false);
            
            /////////////////////
            // EVENT SELECTION //
            /////////////////////
            
            //Declare selection instance
            Selection selection(init_jets_corrected, init_muons, init_electrons, mets, event->kt6PFJets_rho());
/*          selection.setJetCuts(20,2.5,0.01,1.,0.98,0.3,0.1); // standard TOP jet selection
            selection.setMuonCuts(5,2.5,0.4,0.2,0.3,1,0.5,5,0); // standard mu selection but with looser iso
            selection.setElectronCuts(10,2.5,0.4,0.02,0.5,0.3,0); // standard ele selection but with looser iso
*/            
	    //define selection cuts --> have to be validated!!!
	    // From the class Selection the following functions are used: 
	    //      void Selection::setJetCuts(float Pt, float Eta, float EMF, float n90Hits, float fHPD, float dRJetElectron, float dRJetMuon) 
	    //      void Selection::setLooseDiElectronCuts(float ptt, float Eta, float RelIso, MVAid)  
	    //      void Selection::setLooseMuonCuts(float Pt, float Eta, float RelIso) 
	    //void Selection::setDiElectronCuts(float Et, float Eta, float RelIso, float d0, float MVAId, float DistVzPVz, float DRJets, int MaxMissingHits)
		
	    selection.setJetCuts(20.,2.4,0.01,1.,0.98,0.3,0.1); 
	    selection.setDiMuonCuts(10.,2.5,0.2,0.04);
	    selection.setDiElectronCuts(15.0,2.4,0.15,0.04,0.5,1,0.3,1); 
	    	
	    
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24, 2.);
            
            if(!isGoodPV)
                continue;
            
            missingEt=mets[0]->Pt();
	    missingEt_pX=mets[0]->Px();
	    missingEt_pY=mets[0]->Py();
	    missingEt_Phi = mets[0]->Phi(); 
	    missingEt_Theta = mets[0]->Theta();
	    
            
            vector<TRootJet*> selectedJets= selection.GetSelectedJets(true);
            
//            vector<TRootMuon*> selectedMuons = selection.GetSelectedMuons(vertex[0],selectedJets);
            vector<TRootMuon*> selectedMuons = selection.GetSelectedDiMuons();
	    
//            vector<TRootElectron*> selectedElectrons = selection.GetSelectedElectrons(selectedJets);
	    vector<TRootElectron*> selectedElectrons = selection.GetSelectedDiElectrons();
	    
	    vector<TRootJet*> selectedBJets_CSVM; // B-jets at the Tight working point
	    vector<TRootJet*> selectedLightJets; // light-Jets, to be filled afer b-tagging
	    
	     
	   
	    
	     
            
            nElectrons=0;
            for(int iele=0; iele<selectedElectrons.size() && nElectrons<10; iele++){
                pX_electron[nElectrons]=selectedElectrons[iele]->Px();
                pY_electron[nElectrons]=selectedElectrons[iele]->Py();
                pZ_electron[nElectrons]=selectedElectrons[iele]->Pz();
                E_electron[nElectrons]=selectedElectrons[iele]->E();
                Double_t isocorr=0;
                
                // get isolation out, start by getting pu corrections
                if(selectedElectrons[iele]->puChargedHadronIso()>0){
                    isocorr = selectedElectrons[iele]->puChargedHadronIso();
                    
                }
                else{
                    // go through loads of pain to get rho correction, no function available. code below taken from TRootElectron selector in TopTreeAnalysisBase/*/Selector.cc
                    double EffectiveArea = 0.;
                    
                    // HCP 2012 updated for electron conesize = 0.3, taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 0.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 1.0 ) EffectiveArea = 0.130;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 1.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 1.479 ) EffectiveArea = 0.137;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 1.479 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.0 ) EffectiveArea = 0.067;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.0 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.2 ) EffectiveArea = 0.089;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.2 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.3 ) EffectiveArea = 0.107;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.3 && fabs(selectedElectrons[iele]->superClusterEta()) < 2.4 ) EffectiveArea = 0.110;
                    if (fabs(selectedElectrons[iele]->superClusterEta()) >= 2.4) EffectiveArea = 0.138;
                    isocorr = event->kt6PFJets_rho()*EffectiveArea;
                }
                
                pfIso_electron[nElectrons]=(selectedElectrons[iele]->chargedHadronIso() + max( selectedElectrons[iele]->neutralHadronIso() + selectedElectrons[iele]->photonIso()  - isocorr, 0.) )/ selectedElectrons[iele]->Pt();
                charge_electron[nElectrons]=selectedElectrons[iele]->charge();
                nElectrons++;
            }
            nMuons=0;
            for(int imuo=0; imuo<selectedMuons.size() && nMuons<10; imuo++){
                pX_muon[nMuons]=selectedMuons[imuo]->Px();
                pY_muon[nMuons]=selectedMuons[imuo]->Py();
                pZ_muon[nMuons]=selectedMuons[imuo]->Pz();
                E_muon[nMuons]=selectedMuons[imuo]->E();
                pfIso_muon[nMuons]=(selectedMuons[imuo]->chargedHadronIso() + max( 0.0, selectedMuons[imuo]->neutralHadronIso() + selectedMuons[imuo]->photonIso() - 0.5*selectedMuons[imuo]->puChargedHadronIso() ) ) / selectedMuons[imuo]->Pt(); // dBeta corrected


                charge_muon[nMuons]=selectedMuons[imuo]->charge();
                nMuons++;
            }
            nJets=0;
            for(int ijet=0; ijet<selectedJets.size() && nJets<10; ijet++){
                pX_jet[nJets]=selectedJets[ijet]->Px();
                pY_jet[nJets]=selectedJets[ijet]->Py();
                pZ_jet[nJets]=selectedJets[ijet]->Pz();
                E_jet[nJets]=selectedJets[ijet]->E();
                nJets++;
            }
	    
	    vector<int> SelectedBjets_ID; 
	    vector<int> SelectedLightjets_ID; 
	    SelectedBjets_ID.clear(); 
	    SelectedLightjets_ID.clear(); 
	    
	    for(unsigned int iJet=0; iJet<selectedJets.size(); iJet++){
		if (selectedJets[iJet]->btag_combinedSecondaryVertexBJetTags() > .679)
		{
			selectedBJets_CSVM.push_back(selectedJets[iJet]);
			SelectedBjets_ID.push_back(iJet); 
		}
		else
		{
		 	selectedLightJets.push_back(selectedJets[iJet]);
		 	SelectedLightjets_ID.push_back(iJet);
		}
			
	     } 
	     
	    nBJets=0;
            for(int ijet=0; ijet<selectedBJets_CSVM.size() && nBJets<10; ijet++){
                pX_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Px();
                pY_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Py();
                pZ_Bjet[nBJets]=selectedBJets_CSVM[ijet]->Pz();
                E_Bjet[nBJets]=selectedBJets_CSVM[ijet]->E();
                nBJets++;
            }
	    
	    nLJets=0;
            for(int ijet=0; ijet<selectedLightJets.size() && nLJets<10; ijet++){
                pX_Ljet[nLJets]=selectedLightJets[ijet]->Px();
                pY_Ljet[nLJets]=selectedLightJets[ijet]->Py();
                pZ_Ljet[nLJets]=selectedLightJets[ijet]->Pz();
                E_Ljet[nLJets]=selectedLightJets[ijet]->E();
                nLJets++;
            }
	    
	    
	    
	    vector<pair<int,int> >  HighestPtLept;
	    HighestPtLept.clear();
	
	    
	    for(int leptonIt = 0; leptonIt < selectedElectrons.size(); leptonIt++)
	    {
	    	    HighestPtLept.push_back(make_pair(leptonIt, selectedElectrons[leptonIt]->Pt()));
	    }
	    for(int k =0; k < selectedMuons.size() ; k++)
	    {
	    	    HighestPtLept.push_back(make_pair(20+k, selectedMuons[k]->Pt()));
	    	    
	    }		
	    
	    sort(HighestPtLept.begin(), HighestPtLept.end(),cmp_big_first);
	    if(debug)
	    {
	    	cout << "****************New event*****************" << endl; 
	    	for(int i = 0; i<HighestPtLept.size(); i++)
	    	{
	         	pair<int,int> aPair = HighestPtLept[i];
	         	if(aPair.first>19) 
		 	{
		 		int number = aPair.first-20;
				cout << "selectedMuons[number]->Pt: " << selectedMuons[number]->Pt() << endl; 
		 	}
		 	else
		 	{
		 		int number = aPair.first;
				cout << "selectedElectrons[number]->Pt: " << selectedElectrons[number]->Pt() << endl; 
		 
		 	}
	    	}
	    }
	    
	    /////////////////////////////
	    /// Jet parton matching ////
	    /////////////////////////////()
	    
	    std::vector<TLorentzVector> B_partons;
	    std::vector<TLorentzVector> C_partons;
	     std::vector<TLorentzVector> W_partons;
	    std::vector<TLorentzVector> tljets;
	    std::vector<int> MatchingId_b; 
	    std::vector<int> MatchingId_c; 
	    std::vector<int> MatchingId_W; 
	    MatchingId_b.clear(); 
	    MatchingId_c.clear();
	    MatchingId_W.clear();
	    B_partons.clear(); 
	    C_partons.clear(); 
	    W_partons.clear();
	    tljets.clear(); 
	    
	    
	    if(validate)
	    {
	       for(unsigned int i=0;i<selectedJets.size();i++)
	       {
	            tljets.push_back((TLorentzVector)*selectedJets[i]);
		
	       }
	    
	       TLorentzVector bquark; 
	       TLorentzVector cquark;
	        TLorentzVector Wquark;
	       for(unsigned int iMC = 0; iMC< mcParticles.size(); iMC++)
	       {
	       
	    	    //if(debug) cout << iMC << ":  status: " << mcParticles[iMC]->status() << "  pdgId: " << mcParticles[iMC]->type() << " motherPdgId: " << mcParticles[iMC]->motherType() << "  grannyPdgId: " << mcParticles[iMC]->grannyType() << endl;
            	    if( mcParticles[iMC]->status() != 3) continue;    // only those from hard interaction
	      
	    	    if(fabs(mcParticles[iMC]->type()) == 5 && fabs(mcParticles[iMC]->motherType()) == 6 && channelName.find("45")!=string::npos)  
		    {
		        bquark = *mcParticles[iMC]; 
		        
		        if(debug) cout << iMC << ":  status: " << mcParticles[iMC]->status() << "  pdgId: " << mcParticles[iMC]->type() << " motherPdgId: " << mcParticles[iMC]->motherType() << "  grannyPdgId: " << mcParticles[iMC]->grannyType() << endl;
            	        B_partons.push_back(bquark); 
		    }
		    else if(fabs(mcParticles[iMC]->type()) == 4 && fabs(mcParticles[iMC]->motherType()) == 6 )  
		    {
		        cquark = *mcParticles[iMC]; 
		        
		        if(debug) cout << iMC << ":  status: " << mcParticles[iMC]->status() << "  pdgId: " << mcParticles[iMC]->type() << " motherPdgId: " << mcParticles[iMC]->motherType() << "  grannyPdgId: " << mcParticles[iMC]->grannyType() << endl;
            	        C_partons.push_back(cquark); 
		    }
		    else if((fabs(mcParticles[iMC]->type()) == 4 || fabs(mcParticles[iMC]->type()) == 1 || fabs(mcParticles[iMC]->type()) == 2 || fabs(mcParticles[iMC]->type()) == 3 ) &&
		    fabs(mcParticles[iMC]->motherType()) == 24 )  
		    {
		        Wquark = *mcParticles[iMC]; 
		        
		        if(debug) cout << iMC << ":  status: " << mcParticles[iMC]->status() << "  pdgId: " << mcParticles[iMC]->type() << " motherPdgId: " << mcParticles[iMC]->motherType() << "  grannyPdgId: " << mcParticles[iMC]->grannyType() << endl;
            	        W_partons.push_back(Wquark); 
		    }
	        }
		
		if(B_partons.size() == 0 || C_partons.size() == 0 || W_partons.size() == 0) continue;   //when no partons are found, dismiss event
		
		const int TotalMinDist= JetPartonMatching::totalMinDist; 
		
		JetPartonMatching myJetPartonMatcher_b = JetPartonMatching(B_partons,tljets,TotalMinDist, true,true ,0.3);
		JetPartonMatching myJetPartonMatcher_c = JetPartonMatching(C_partons,tljets,TotalMinDist, true,true ,0.3);
		JetPartonMatching myJetPartonMatcher_W = JetPartonMatching(W_partons,tljets,TotalMinDist, true,true ,0.3);
		if(debug) 
		{
			cout << "B PARTONS" << endl; 
			myJetPartonMatcher_b.print();
			cout << "C PARTONS" << endl; 
			myJetPartonMatcher_c.print();
			cout << "PARTONS FROM W" << endl; 
			myJetPartonMatcher_W.print();
		        cout << endl; 
		}
			
		// get the matching jet IDs
		for(int iID = 0; iID < B_partons.size(); iID++)
		{
		    int ID = -1; 
		    
		    ID = myJetPartonMatcher_b.getMatchForParton(iID,0); 
		    MatchingId_b.push_back(ID); 
		}
		for(int iID = 0; iID < C_partons.size(); iID++)
		{
		    int ID = -1; 
		    
		    ID = myJetPartonMatcher_c.getMatchForParton(iID,0); 
		    MatchingId_c.push_back(ID); 
		}
		for(int iID = 0; iID < W_partons.size(); iID++)
		{
		    int ID = -1; 
		    
		    ID = myJetPartonMatcher_W.getMatchForParton(iID,0); 
		    MatchingId_W.push_back(ID); 
		}
		
		
		
	    }
	    
	    
	    //////////////////////////////////
	    /// variables for 45          /////
	    ///////////////////////////////////
	    
	    TLorentzVector Reco_FCNC_top_combi; 
	    TLorentzVector leptonpair_1; // two OSSF with highest pt
	    TLorentzVector leptonpair_2; // two OSSF with 2nd highest pt
	    TLorentzVector leptonfour; // decay of Z
	    TLorentzVector lepton_0;
	    TLorentzVector lepton_1;
	    TLorentzVector lepton_2;
	    TLorentzVector lepton_3;
	    TLorentzVector lepton_4;
	    TLorentzVector missingEt_vector1; 
	    TLorentzVector missingEt_vector2; 
	    

	    
	    
	    if(channelName.find("45")!=string::npos && (nElectrons+nMuons>3))
	    { 		
	    	 if(debug) cout << "In nElectrons + nMuons > 3" << endl; 
		 
		
		 
		 bool Zdecay = false ; // defines if the leptons are coming from a Z decay  
		 
		 //define all leptons
		 if(nElectrons+nMuons > 4)
		 {
		     if(debug) cout << "In nElectrons + nMuons > 4" << endl; 
		     //define all leptons
		     vector <TLorentzVector> vectors_5leptons = Channel_45_5leptons( debug, selectedElectrons,selectedMuons,missingEt_pX,missingEt_pY, missingEt_Theta, missingEt);
		     leptonfour = vectors_5leptons[0]; //vector of the four leptons voming from the Z ~higgs
		     lepton_4 = vectors_5leptons[1];  // extra lepton from SM W decay
		     missingEt_vector1 = vectors_5leptons[2];  // met candidate
		     missingEt_vector2 = vectors_5leptons[3];  // met candidate
		     
		     TLorentzVector boolV = vectors_5leptons[4]; 
		     if(boolV.Px() == 1) {
		     	Zdecay = true; 
			
		     }
		     
		     if(debug) cout << "out nElectrons + nMuons > 4" << endl; 
		 
		 }
		 else if(nElectrons+nMuons == 4)
		 {
		     if(debug) cout << "In nElectrons + nMuons = 4" << endl; 
		     vector <TLorentzVector> vectors_4leptons = Channel_45_4leptons( debug, selectedElectrons, selectedMuons);
		     leptonfour = vectors_4leptons[0]; // vector of the four leptons voming from the Z ~higgs
		     
		     TLorentzVector boolV = vectors_4leptons[1]; 
		     if(boolV.Px() == 1){
		     	 Zdecay = true; 
			 
		     }
		      
		     if(debug) cout << "Out nElectrons + nMuons = 4" << endl; 
		 }
		 
		
		 
		 bool requirements = false; 
		 if(Zdecay ) requirements = true; 
		 
		 //if the requirements aren't fullfilled , go to the next event
		 if(channelName.find("45")!=string::npos && !requirements ) continue; 
		 
		 
		 
		 
                 //define the higss invariant mass H -> ZZ* -> ll ll
	    	 if(debug) cout << "Zdecay leptonfour.M()= " << leptonfour.M() << endl; 
	    	 InvMass_4lept_Zdecay = leptonfour.M();
		 Phi_Higgs = leptonfour.Phi(); 
		 Eta_Higgs = leptonfour.Eta();
		
		 //make a pair of light jets with booleans such that I know which ones are matched
		 vector<pair<int,bool> >  LightJets_Paired;
	    	 LightJets_Paired.clear();
	    	 for(int It = 0; It < selectedLightJets.size(); It++)
	    	 {
	    	 	 LightJets_Paired.push_back(make_pair(SelectedLightjets_ID[It], false));
	    	 }
		 /*
		 // check the pairing
		 for(int It = 0; It < selectedLightJets.size(); It++)
	    	 {
	    	 	cout << "Iterator: " << It << " JetID: " <<  LightJets_Paired[It].first << " Bool: " << LightJets_Paired[It].second << endl; 
		 }
		 */
		 
		 //Search fcnc c jet
		 vector<pair<int,bool> > LightJets_Paired_c; 
		 
		 
		 if(nLJets>0)
		 {
		 	//find the right c-jet from FCNC decay
		 	LightJets_Paired_c = Channel_45_FCNC_cjet(debug, LightJets_Paired, selectedLightJets, leptonfour); 
		 	/*
			// check pairing
			for(int It = 0; It < selectedLightJets.size(); It++)
	    		{
	    	 	    cout << "Iterator: " << It << " JetID: " <<  LightJets_Paired_c[It].first << " Bool: " << LightJets_Paired_c[It].second << endl; 
		 	}
			*/ 
			
			//Calculate FCNC top from the chosen c-jet and Higgs
		 	TLorentzVector FCNC_top = Channel_45_FCNC_top(debug, LightJets_Paired_c,leptonfour,selectedLightJets);
			InvMass_FCNC_top_Zdecay = FCNC_top.M(); 
			FCNC_top_eta = FCNC_top.Eta(); 
			FCNC_top_phi = FCNC_top.Phi(); 
			FCNC_top_Pt = FCNC_top.Pt(); 
			
			if(validate && MatchingId_c[0] != -1 )
			{
			   EventsToMatch_c++;
			   
			   int ChosenID = -1; 
			   
			   for( int iJet = 0; iJet < LightJets_Paired_c.size(); iJet++)
	 		   {
	 	 		
	 	 		pair<int,int> Pair = LightJets_Paired_c[iJet];
	 	 	
				if(Pair.second)
		 		{	
		 			ChosenID = Pair.first; 
	 	 		}
		 		
	 		    }
			    if(ChosenID == MatchingId_c[0])
			    {
			    	MatchedCounter_c++; 
			    
			    
			    
			    	//extra safety 
			    	Selected_cjet_Pt -> Fill(selectedJets[ChosenID]->Pt()); 
			    	True_cjet_Pt -> Fill(selectedJets[MatchingId_c[0]]->Pt()); 
			    }
		 	 }
			 
		 }
		
		 
		 TLorentzVector tempBjet ;		 
		 if(nBJets > 0)
		 {
		 	 // tag b jet with highest discrimanting power as the SM one
		 	 
			 vector <TLorentzVector> highestDisc = SM_b(debug,selectedBJets_CSVM); 
		         tempBjet = highestDisc[0];
			 Bdiscr = highestDisc[1].Px(); 
			 
			 
			 if(validate && MatchingId_b[0] != -1)
			 {
			 	if(debug) cout << "in validation" << endl; 
				int Bjet_ID = highestDisc[2].Px();
				if(debug)cout << "Place in bjets " << Bjet_ID << endl; 
			 	int BjetPlace_in_selectedJets = SelectedBjets_ID[Bjet_ID]; 
			 	if(debug)cout << "Place in jets " << BjetPlace_in_selectedJets << endl; 
				if(debug)cout << "MatchingId size " << MatchingId_b.size() << endl; 
				if(debug)cout << "MatchingId " << MatchingId_b[0] << endl; 
				
				EventsToMatch++;
				
			 	if(BjetPlace_in_selectedJets == MatchingId_b[0] )
				{
				   if(debug) cout << "right match" << endl; 
				   MatchedCounter++;
				   //extra safety 
					Selected_bjet_Pt -> Fill(selectedJets[BjetPlace_in_selectedJets]->Pt()); 
					True_bjet_Pt -> Fill(selectedJets[MatchingId_b[0]]->Pt()); 
				   
			 	}
				else
				{
					if(debug) cout << "wrong match" << endl; 
				}
				
			 }
		 	 // W boson decays hadronically
		 	 if(nLJets>2 && (nMuons + nElectrons == 4))
		 	 {
		 	  
		 	 
		 		 //find the two jets closest to this b 
				 vector <TLorentzVector> WjetCalc = Channel_45_SM_Wqq(debug,tempBjet, LightJets_Paired_c, selectedLightJets);
		 		 
		 		 TLorentzVector Wjets = WjetCalc[0]; 
				 int Wjet_ID_1 = WjetCalc[1].Px();
				 int Wjet_ID_2 = WjetCalc[2].Px();
				 
				 
				 if(validate)
				 { 
				 	for(int i = 0; i < MatchingId_W.size(); i++)
					{
					    for(int k = 0; k < MatchingId_W.size(); k++)
					    {
					      if(k != i)
					      {
						if(MatchingId_W[i] != -1  && MatchingId_W[k] != -1)
						{
							EventsToMatch_W++;
			   
			   				
			   
			   				if(Wjet_ID_1 == MatchingId_W[i] && Wjet_ID_2 == MatchingId_W[k])
			    				{
			    					MatchedCounter_W++; 
			    				}
							else if(Wjet_ID_2 == MatchingId_W[i] && Wjet_ID_1 == MatchingId_W[k])
			    				{
			    					MatchedCounter_W++; 
			    				}
		 	 			}
					      }
					    }
					}
				}
				 
				 TrMass_W_qq = Wjets.Mt();
		 		 TrMass_W = TrMass_W_qq; 
		 		 
		 		 InvMass_SM_W_qq = Wjets.M(); 
				 InvMass_SM_W = InvMass_SM_W_qq; 
				
				 TLorentzVector TopJet = Wjets + tempBjet;  
		 		 InvMass_SM_top_bqq = TopJet.M();
		 		 InvMass_SM_top = InvMass_SM_top_bqq;
		 		 	
				 
		 		 
		 	 }
		 	 // W boson decays leptonically 
		 	 else if(nMuons + nElectrons > 4)
		 	 {
		 		 //define the missing et 
		 		 TLorentzVector missingEt_vector; 
		 	         bool NonNegMET = false; 
		
		 		 TLorentzVector SumInv1; 
		 		 SumInv1 = missingEt_vector1 + lepton_4; 
		
		 		 TLorentzVector SumInv2; 
		 		 SumInv2 = missingEt_vector2 + lepton_4;
		 	 
		 		 if(SumInv1.M() <0 && SumInv2.M()> 0)
		 		 {
				    NonNegMET = true; 
		 		    missingEt_vector = missingEt_vector2; 
		 		    InvMass_SM_W_lv = SumInv2.M();
		 		 }
		 		 else if(SumInv1.M() >0 && SumInv2.M()< 0)
				 {
				 	NonNegMET = true;
					missingEt_vector = missingEt_vector1; 
		 			InvMass_SM_W_lv = SumInv1.M();
				 }
				 else if(SumInv1.M() >0 && SumInv2.M() > 0)
		 		 {
		 			 NonNegMET = true;
					 if(fabs(SumInv1.M() - 80.4) < fabs(SumInv2.M() - 80.4))
		 			 {
		 				 InvMass_SM_W_lv = SumInv1.M();
		 				 missingEt_vector = missingEt_vector1;
		 			 }
		 			 else
		 			 {
		 				 InvMass_SM_W_lv = SumInv2.M(); 
		 				 missingEt_vector = missingEt_vector2; 
		 			 }
		 		 }
				 TLorentzVector combi; 
		 		 combi = lepton_4 + tempBjet; 
		 		 InvMass_SM_lb =combi.M();
				 
				 //when there is a good neutrino candidate found
				 if(NonNegMET)
				 {
				 	InvMass_SM_W = InvMass_SM_W_lv;

				 	TLorentzVector sum; 
		 		 	sum = lepton_4+missingEt_vector;
	    	 		 	TrMass_W_lv = sum.Mt();
		 		 	TrMass_W = TrMass_W_lv; 

		 		 	//make SM top from blv
		 		 	TLorentzVector combi2;
		 		 	combi2 = missingEt_vector + lepton_4 + tempBjet; 
		 		 	InvMass_SM_top_blv = combi2.M();
		 		 	InvMass_SM_top = InvMass_SM_top_blv;
				}
		 	 } // leptonic W
		
		 }  // nbjets > 0
		
		 myTree->Fill(); 
		 nEvents_Selected[d]++;
		 

	    } // > 3 leptons
	    
	    if(channelName.find("3L")!=string::npos && (nElectrons+nMuons == 3))
	    {
	    	if(debug)cout << " In nleptons == 3 " << endl; 
		
		TLorentzVector tempLeptons; 
		tempLeptons.Clear(); 
		
		
		
		for(int i = 0; i < selectedElectrons.size(); i++)
		{
		 	TLorentzVector electron; 
			electron.SetPxPyPzE(selectedElectrons[i]->Px(), selectedElectrons[i]->Py(), selectedElectrons[i]->Pz(), selectedElectrons[i]->Energy()); 
			
			tempLeptons = tempLeptons + electron; 
			
		
		
		}
		for(int i = 0; i< selectedMuons.size(); i++)
		{
		 	 TLorentzVector Muon; 
		 	 Muon.SetPxPyPzE(selectedMuons[i]->Px(), selectedMuons[i]->Py(), selectedMuons[i]->Pz(), selectedMuons[i]->Energy()); 
		 	 tempLeptons = tempLeptons+Muon; 
		}
		InvMass_Leptons = tempLeptons.M(); 
		
		if(debug)cout << " determined invariant mass of all 3 leptons " << endl; 
		
		
		Double_t temp_DR = 0.; 
		Double_t temp_DRmin = 100000.; 
		TLorentzVector FCNC_lepton1; 
		TLorentzVector FCNC_lepton2; 
		
		bool FCNC_lepton1_isMuon = false; 
		bool FCNC_lepton1_isElectron = false; 
		bool FCNC_lepton2_isMuon = false; 
		bool FCNC_lepton2_isElectron = false;
		bool FCNC_lepton1_isM= false; 
		bool FCNC_lepton1_isE = false; 
		bool FCNC_lepton2_isM= false; 
		bool FCNC_lepton2_isE = false; 
		int Number_1 = -1; 
		int Number_2 = -1; 
		
			
	    	for(int i = 0; i<HighestPtLept.size(); i++)
	    	{
	         	pair<int,int> aPair = HighestPtLept[i];
			
			FCNC_lepton1_isM= false; 
			FCNC_lepton1_isE = false; 
			 
			
			TLorentzVector lepton1; 
			TLorentzVector lepton2; 
			
			
			 
			
			int number1 = -1;
			int number2 = -1;
	         	if(aPair.first>19) 
		 	{
		 		number1 = aPair.first-20;
				//if(debug) cout << "selectedMuons[number1]->Pt: " << selectedMuons[number1]->Pt() << endl; 
				lepton1.SetPxPyPzE(selectedMuons[number1]->Px(), selectedMuons[number1]->Py(), selectedMuons[number1]->Pz(), selectedMuons[number1]->Energy());
		 		FCNC_lepton1_isM = true; 
			}
		 	else if(aPair.first <20)
		 	{
		 		number1 = aPair.first;
				//if(debug) cout << "selectedElectrons[number1]->Pt: " << selectedElectrons[number1]->Pt() << endl; 
		 		lepton1.SetPxPyPzE(selectedElectrons[number1]->Px(), selectedElectrons[number1]->Py(), selectedElectrons[number1]->Pz(), selectedElectrons[number1]->Energy());
		 		FCNC_lepton1_isE = true; 
			}
			if(FCNC_lepton1_isM && FCNC_lepton1_isE) cout << "WARNING: something is wrong 1" << endl; 
			
			for(int k = 0; k < HighestPtLept.size(); k++)
			{
				FCNC_lepton2_isM= false; 
			        FCNC_lepton2_isE = false;
				if(k != i)
				{
				  pair<int,int> pPair = HighestPtLept[k];
				  //cout << "k: " << k << " i: " << i << " pPair.first " << pPair.first << " aPair.first " << aPair.first <<endl; 
				  if(pPair.first>19) 
		 		  {
		 			number2 = pPair.first-20;
					//if(debug) cout << "selectedMuons[number2]->Pt: " << selectedMuons[number2]->Pt() << endl; 
					lepton2.SetPxPyPzE(selectedMuons[number2]->Px(), selectedMuons[number2]->Py(), selectedMuons[number2]->Pz(), selectedMuons[number2]->Energy());
		 			FCNC_lepton2_isM = true; 
				  }
		 		  else if(pPair.first <20)
		 		  {
		 			number2 = pPair.first;
					//if(debug) cout << "selectedElectrons[number2]->Pt: " << selectedElectrons[number2]->Pt() << endl; 
		 			lepton2.SetPxPyPzE(selectedElectrons[number2]->Px(), selectedElectrons[number2]->Py(), selectedElectrons[number2]->Pz(), selectedElectrons[number2]->Energy());
		 			FCNC_lepton2_isE = true; 
				  }
				  if(FCNC_lepton2_isM && FCNC_lepton2_isE) cout << "WARNING: something is wrong 2" << endl; 
			
			
				  double DR = lepton1.DeltaR(lepton2);
				  if(DR>temp_DR)
				  { 	temp_DR = DR;
			        	if(debug) cout << "changed DR to " << temp_DR << endl; 
				  }
				  if(DR<temp_DRmin)
				  { 
					temp_DRmin = DR; 
					if(debug) cout << "changed DRmin to " << temp_DRmin << endl; 
					FCNC_lepton1 = lepton1; 
					FCNC_lepton2 = lepton2; 
					Number_1 = number1; 
					Number_2 = number2; 
					FCNC_lepton1_isMuon = FCNC_lepton1_isM;
					FCNC_lepton1_isElectron = FCNC_lepton1_isE; 
					FCNC_lepton2_isMuon = FCNC_lepton2_isM; 
					FCNC_lepton2_isElectron = FCNC_lepton2_isE; 
				  }
				}
			}
	    	}
	        if(debug)cout << " out DR " << endl; 
		DeltaR = temp_DR; 
		DeltaRmin = temp_DRmin; 
		TLorentzVector tempFCNC = FCNC_lepton1 + FCNC_lepton2; 
		InvMass_FCNC_Leptons = tempFCNC.M();
		Eta_FCNC_Leptons = tempFCNC.Eta(); 
		Pt_FCNC_Leptons = tempFCNC.Pt(); 
		Phi_FCNC_Leptons = tempFCNC.Phi(); 
		Theta_FCNC_Leptons = tempFCNC.Theta(); 
		
		
		TLorentzVector SMLepton; 
		
		//define remaining lepton
		if(selectedMuons.size() == 3)
		{
			if(Number_1 != 0 && Number_2 != 0) SMLepton.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(),selectedMuons[0]->Pz(), selectedMuons[0]->Energy()); 
			if(Number_1 != 1 && Number_2 != 1) SMLepton.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(),selectedMuons[1]->Pz(), selectedMuons[1]->Energy()); 
			if(Number_1 != 2 && Number_2 != 2) SMLepton.SetPxPyPzE(selectedMuons[2]->Px(), selectedMuons[2]->Py(),selectedMuons[2]->Pz(), selectedMuons[2]->Energy()); 
			
		}
		else if(selectedElectrons.size() == 3)
		{
			if(Number_1 != 0 && Number_2 != 0) SMLepton.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy()); 
			if(Number_1 != 1 && Number_2 != 1) SMLepton.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy()); 
			if(Number_1 != 2 && Number_2 != 2) SMLepton.SetPxPyPzE(selectedElectrons[2]->Px(), selectedElectrons[2]->Py(),selectedElectrons[2]->Pz(), selectedElectrons[2]->Energy()); 
			
		}
		else if(selectedMuons.size()==1)
		{
			if(FCNC_lepton1_isMuon)
			{
				 
				if(Number_2 != 0) SMLepton.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy()); 
				if(Number_2 != 1) SMLepton.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy()); 
			}
			else if(FCNC_lepton2_isMuon)
			{
				if(Number_1 != 0 ) SMLepton.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy()); 
				if(Number_1 != 1 ) SMLepton.SetPxPyPzE(selectedElectrons[1]->Px(), selectedElectrons[1]->Py(),selectedElectrons[1]->Pz(), selectedElectrons[1]->Energy()); 
				
			}
			else
			{
				SMLepton.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(),selectedMuons[0]->Pz(), selectedMuons[0]->Energy()); 
			}
		}
		else if(selectedElectrons.size()==1)
		{
			if(FCNC_lepton1_isElectron)
			{
				 
				if(Number_2 != 0) SMLepton.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(),selectedMuons[0]->Pz(), selectedMuons[0]->Energy()); 
				if(Number_2 != 1) SMLepton.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(),selectedMuons[1]->Pz(), selectedMuons[1]->Energy()); 
			}
			else if(FCNC_lepton2_isElectron)
			{
				if(Number_1 != 0 ) SMLepton.SetPxPyPzE(selectedMuons[0]->Px(), selectedMuons[0]->Py(),selectedMuons[0]->Pz(), selectedMuons[0]->Energy()); 
				if(Number_1 != 1 ) SMLepton.SetPxPyPzE(selectedMuons[1]->Px(), selectedMuons[1]->Py(),selectedMuons[1]->Pz(), selectedMuons[1]->Energy()); 
				
			}
			else
			{
				SMLepton.SetPxPyPzE(selectedElectrons[0]->Px(), selectedElectrons[0]->Py(),selectedElectrons[0]->Pz(), selectedElectrons[0]->Energy()); 
			}
		}
		
		Eta_SM_Lepton = SMLepton.Eta(); 
		Pt_SM_Lepton = SMLepton.Pt(); 
		InvMass_SM_Lepton = SMLepton.M(); 
		
		bool tcZcandidate = false; 
		bool tcWW = false; 
		bool tcZZcandidate = false; 
		
		
		if(FCNC_lepton1_isMuon && FCNC_lepton2_isElectron) 
		{
				tcWW = true; 
				if(debug) cout << "OF" << endl; 
		}
		else if(FCNC_lepton2_isMuon && FCNC_lepton1_isElectron)
		{
				tcWW = true; 
				if(debug) cout << "OF" << endl;
		}
		else 
		{
			bool OS = false; 
			bool SS = false; 
			
			if(FCNC_lepton1_isMuon && FCNC_lepton2_isMuon)
			{
				if(selectedMuons[Number_1]->charge() == selectedMuons[Number_2]->charge()) SS = true; 
				else OS = true; 
				//cout << Number_1 << " selectedMuons[Number_1]->charge() " << selectedMuons[Number_1]->charge() << endl; 
				//cout << Number_2 << " selectedMuons[Number_2]->charge() " << selectedMuons[Number_2]->charge() << endl; 
			}
			if(FCNC_lepton1_isElectron && FCNC_lepton2_isElectron)
			{
				if(selectedElectrons[Number_1]->charge() == selectedElectrons[Number_2]->charge()) SS = true; 
				else OS = true; 
				//cout << Number_1 << " selectedElectrons[Number_1]->charge() " << selectedElectrons[Number_1]->charge() << endl; 
				//cout << Number_2 <<" selectedElectrons[Number_2]->charge() " << selectedElectrons[Number_2]->charge() << endl; 
			}
			
			if(OS)
			{
				tcZcandidate = true; 
				tcZZcandidate = true; 
				if(debug) cout << "OS" << endl; 
			}
			else if (SS)
			{
				tcWW = true; 
				if(debug) cout << "SS" << endl; 
			}
		}
		
		if(debug)
		{
		  if(tcZcandidate) cout << "Classified as tcZcandidate" << endl; 
		  if(tcZZcandidate) cout << "Classified as tcZZcandidate" << endl; 
		  if(tcWW) cout << "Classified as tcWW" << endl; 
		  cout << endl; 
		}
		
		
		//// add the c jet
		if(tcWW && nLJets >0)
		{
		    TLorentzVector Top_FCNC_Candidate = Channel_3L_FCNC_top_candidate(debug, tempFCNC, selectedLightJets); 
				
		    if(debug) cout << "Top_FCNC_ candidate defined" << endl; 
		    pT_FCNC_top_candidate = Top_FCNC_Candidate.Pt(); 
		    Eta_FCNC_top_candidate = Top_FCNC_Candidate.Eta();
		    InvMass_FCNC_top_candidate = Top_FCNC_Candidate.M(); 
		}	    
		
		//add b jet 
		if(tcWW && nBJets > 0)
		{
			vector <TLorentzVector> highestDisc = SM_b(debug,selectedBJets_CSVM); 
		    	TLorentzVector tempBjet; 
		    	tempBjet = highestDisc[0];
		    	Bdiscr = highestDisc[1].Px(); 
		    	
			TLorentzVector SMtop; 
			SMtop = tempBjet + SMLepton; 
			InvMass_SM_top = SMtop.M(); 
			Eta_SM_top = SMtop.Eta(); 
			Pt_SM_top = SMtop.Pt(); 
			
		
		}
		bool tcZ = false; 
		bool tcZZ = false; 
		if((tcZcandidate || tcZZcandidate) && nLJets > 0)
		{
		    TLorentzVector Top_FCNC_Candidate = Channel_3L_FCNC_top_candidate(debug, tempFCNC, selectedLightJets); 
				
		    if(debug) cout << "Top_FCNC_ candidate defined" << endl; 
		    pT_FCNC_top_candidate = Top_FCNC_Candidate.Pt(); 
		    Eta_FCNC_top_candidate = Top_FCNC_Candidate.Eta();
		    InvMass_FCNC_top_candidate = Top_FCNC_Candidate.M(); 
		    
		    if(InvMass_FCNC_top_candidate > 160)   tcZ = true; 
		    else tcZZ = true; 
		
		}
		if((tcZcandidate ||tcZZcandidate) && nBJets > 0)
		{
			vector <TLorentzVector> highestDisc = SM_b(debug,selectedBJets_CSVM); 
		    	TLorentzVector tempBjet; 
		    	tempBjet = highestDisc[0];
		    	Bdiscr = highestDisc[1].Px(); 
		    	
			TLorentzVector SMtop; 
			SMtop = tempBjet + SMLepton; 
			InvMass_SM_top = SMtop.M(); 
			Eta_SM_top = SMtop.Eta(); 
			Pt_SM_top = SMtop.Pt(); 
		
		}
		
		
		
		
		
		
		/*
		//tcZ 
		vector <TLorentzVector> Zcandidates =  Channel_3L_Zcandidate(debug, selectedElectrons, selectedMuons);
		if(debug)cout << " defined Zcandidates " << endl; 
		
		
		
		if(Zcandidates[3].Px() == 1)  //when candidates are found
		{
			InvMass_Z = Channel_3L_Zcandidate_Mll(debug, Zcandidates); 
			if(debug)cout << " InvMass_Z " << endl; 
			TLorentzVector Z_candidate = Zcandidates[0] + Zcandidates[1]; 
			Z_candidate_pT = Z_candidate.Pt(); 
			if(debug)cout << " pt Z " << endl; 
			Z_candidate_Eta = Z_candidate.Eta(); 
			if(debug)cout << " eta Z " << endl; 
			
			if(InvMass_Z < 120 && InvMass_Z > 70)  // removes H -> WW
			{
			    // find the light jet closest to the Z candidate in deltaR
			    if(nLJets > 0)
			    {
			    	TLorentzVector Top_FCNC_Candidate = Channel_3L_FCNC_top_candidate(debug, Z_candidate, selectedLightJets); 
				
			    	if(debug) cout << "Top_FCNC_ candidate defined" << endl; 
			    	pT_FCNC_top_candidate = Top_FCNC_Candidate.Pt(); 
			    	if(debug) cout << "Top_FCNC_ candidate pt" << endl; 
			    	Eta_FCNC_top_candidate = Top_FCNC_Candidate.Eta();
			    	if(debug) cout << "Top_FCNC_ candidate eta" << endl;  
			    	InvMass_FCNC_top_candidate = Top_FCNC_Candidate.M(); 
			    	if(debug) cout << "Top_FCNC_ candidate invariant mass" << endl; 
			    
			    
			    
			    
			    	if(Top_FCNC_Candidate.M() > 160 )
				{
					pT_FCNC_top_tcZ = Top_FCNC_Candidate.Pt(); 
			    		Eta_FCNC_top_tcZ = Top_FCNC_Candidate.Eta();
			    		InvMass_FCNC_top_tcZ = Top_FCNC_Candidate.M(); 
			    		
								
				
				} //inside FCNC top mass window
				else 
				{
					if(nLJets > 2) //Z -> qq
					{
						TLorentzVector HiggsCandidate = Channel_3L_Higgs_candidate(debug, Z_candidate, selectedLightJets); 
						if(debug) cout << "Higgs candidate defined" << endl; 
						H_candidate_pT = HiggsCandidate.Pt(); 
						if(debug) cout << "Higgs candidate pt" << endl; 
						H_candidate_Eta = HiggsCandidate.Eta();
						if(debug) cout << "Higgs candidate eta" << endl;  
						InvMass_H = HiggsCandidate.M(); 
						if(debug) cout << "Higgs candidate invariant mass" << endl; 
				
						TLorentzVector Top_FCNC_Candidate_tcH_ZZ_llqq =	Channel_3L_FCNC_top_candidate2(debug, HiggsCandidate, selectedLightJets);
						if(debug) cout << "Top_FCNC_ candidate defined" << endl; 
						pT_FCNC_top_tcH_ZZ_llqq = Top_FCNC_Candidate_tcH_ZZ_llqq.Pt(); 
						if(debug) cout << "Top_FCNC_ candidate pt" << endl; 
						Eta_FCNC_top_tcH_ZZ_llqq = Top_FCNC_Candidate_tcH_ZZ_llqq.Eta();
						if(debug) cout << "Top_FCNC_ candidate eta" << endl;  
						InvMass_FCNC_top_tcH_ZZ_llqq = Top_FCNC_Candidate.M(); 
						if(debug) cout << "Top_FCNC_ candidate invariant mass" << endl; 
					}
					else //Z -> vv
					{
					
					
					
					}
				} // outside fcnc topmass window
			    
			    
			    } // nLJets > 0
			
			
			
			
			} // inside Zmass window
			else   // in H --> WW
			{	
			
				//choose b-jet with highest bdisc
		    		if(nBJets > 0)
		    		{
		    	    		vector <TLorentzVector> highestDisc = SM_b(debug,selectedBJets_CSVM); 
		    	    		TLorentzVector tempBjet; 
		    	    		tempBjet = highestDisc[0];
		    	    		Bdiscr = highestDisc[1].Px(); 
		    	    		Bjet_Eta= tempBjet.Eta(); 
		    	    		Bjet_Phi = tempBjet.Phi(); 
		    	    		Bjet_Pt = tempBjet.Pt();
		    	    		Bjet_Px = tempBjet.Px(); 
		    	   		Bjet_Py = tempBjet.Py();
		    	    		Bjet_Pz = tempBjet.Pz();
		    	    
		    	    		//choose lepton closest to this bjet 
		    	   		vector <pair<TLorentzVector,bool> > Wlepton = Channel_3L_SM_lep(debug,selectedMuons, selectedElectrons,tempBjet);
		    	    		TLorentzVector SM_lepton; 
		    	    
		    	    		for(int i = 0; i< Wlepton.size();i++)
		    	    		{
		    		   		 pair<TLorentzVector,bool> Pair = Wlepton[i]; 
		    		    		if(Pair.second)  SM_lepton = Pair.first; 
		    	    
		    	    		}
		    	    
		    	    		TLorentzVector combi; 
		    	    		combi = tempBjet + SM_lepton;
		    	    		InvMass_SM_lb = combi.M(); 
		    	    
		    	    		TLorentzVector Hleptons = Channel_3L_Higgs_lep(debug,selectedMuons,selectedElectrons,Wlepton);
		    	    		InvMass_FCNC_ll = Hleptons.M(); 
		    	    		FCNC_ll_Phi=Hleptons.Phi(); 
		    	    		FCNC_ll_Eta= Hleptons.Eta(); 
		    	    		FCNC_ll_Px= Hleptons.Px(); 
		    	    		FCNC_ll_Pz = Hleptons.Pz();
		    	    		FCNC_ll_Py= Hleptons.Py(); 
		    	    		FCNC_ll_Pt= Hleptons.Pt(); 
		    	    
		    	    		DeltaR_SMlb_FCNCll = combi.DeltaR(Hleptons); 
		    	    		DeltaPhi_SMlb_FCNCll = combi.DeltaPhi(Hleptons);
		    	    
		    	    
		    		} // nBjets > 0	
				
			} // end H--> WW   (outside Z mass window

				
			
			
		} //Zcandidate found
		
		
		
		
		*/
	    	myTree->Fill();
		nEvents_Selected[d]++;
		if(debug) cout << "filled tree for 3l channel" << endl; 
		
	    }
	    
	    
	    
	    






        }			//loop on events
        
	
	
	cout<<endl;
	cout<<endl;	
	if(validate) cout << "Bjet matching efficiency is " << (double) (MatchedCounter/EventsToMatch)*100 << " % or " << MatchedCounter << " of " << EventsToMatch  << " SM b jets" <<endl; 
	if(validate) cout << "Cjet matching efficiency is " << (double) (MatchedCounter_c/EventsToMatch_c)*100 << " % or " <<MatchedCounter_c << " of " << EventsToMatch_c  << " FCNC c jets" <<endl; 
	if(validate) cout << "jet matching from W efficiency is " << (double) (MatchedCounter_W/EventsToMatch_W)*100 << " % or " << MatchedCounter_W << " of " << EventsToMatch_W  << " SM W jets" <<endl; 
	
	if(validate){
	 SM_b_selection_efficiency->SetBinContent(1,MatchedCounter/EventsToMatch); 
	 FCNC_c_selection_efficiency->SetBinContent(1,MatchedCounter_c/EventsToMatch_c);
	  SM_Wjet_selection_efficiency->SetBinContent(1,MatchedCounter_W/EventsToMatch_W); 
	}
	else
	{
	 SM_b_selection_efficiency->SetBinContent(1,-1);
	 FCNC_c_selection_efficiency->SetBinContent(1,-1);
	}
	
	cout<<endl;
        
        
        //////////////
        // CLEANING //
        //////////////
        
        if (jecUnc) delete jecUnc;
        if (jetTools) delete jetTools;
        
        //myTree->Write();
	configTree->Write("", TObject::kOverwrite);
	myTree->Write("", TObject::kOverwrite);
        fileout->Write();
        fileout->Close();
	//        delete myTree;
        delete fileout;
        
        //important: free memory
        treeLoader.UnLoadDataset();
        
    }				//loop on datasets
    
    //Once everything is filled 
    //Once everything is filled ...
        if (debug)
            cout << " We ran over all the data ;-)" << endl;
    
    // Do some special things with certain plots (normalize, BayesDivide, ... )
    //    if (debug)
    //        cout << "Treating the special plots." << endl;
    
    delete tcdatasets;
    delete tcAnaEnv;
    delete configTree;
    
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;
    
    return 0;
}
