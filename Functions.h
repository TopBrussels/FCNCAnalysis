//
//  Functions.hpp
//  
//
//  Created by Isis Van Parijs on 12/08/16.
//
//

#ifndef Functions_hpp
#define Functions_hpp

#include <stdio.h>

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
#include <iostream>
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
#include "Functions.h"
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
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
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
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


class Functions
{
  
  public :
  
  //Methods
  Functions();
  ~Functions(){};
  
  void Set_Luminosity(double); //Set the luminosity re-scaling factor to be used thoughout the code


  
  //Initializing CSVv2 b-tag WP
  float workingpointvalue_Loose = -1;
  float workingpointvalue_Medium = -1;
  float workingpointvalue_Tight = -1;
  
  
  //What you want to do
  bool synchex = false;
  bool Assigned = false;
  
  
  // home made functions
  int FCNCjetCalculator(std::vector<TRootPFJet*> Jets, TLorentzVector recoZ ,int index, int verb);
  int SMjetCalculator(std::vector<TRootPFJet*> Jets,int verb);
  double MEtz(bool mu, bool el, TLorentzVector Wlep, double MetPx, double MetPy);
  float EffectiveAreaRho(TRootElectron *el, float _rho) ;
  float EffectiveArea(TRootElectron *el) ;
  float relPfIsoEl(TRootElectron *el, float _rho);
  float IsoDBeta(TRootMuon *mu);
  vector<TLorentzVector> LeptonAssigner(std::vector<TRootElectron*> electrons,std::vector<TRootMuon*> muons);
  TLorentzVector MetzCalculator(TLorentzVector leptW, TLorentzVector v_met);
  vector< pair<unsigned int, unsigned int> > JetPartonPair;
  
  
  // administration functions
  string ConvertIntToString(int Number, bool pad);
  string MakeTimeStamp();
  
  
  
  
  
  
    // members
//   bool stop_program;
  double M_W  = 80.4;
  double M_mu =  0.10566; // 105.66 MeV/c^2
  double M_el = 0.000510999; // 0.510998910 Mev/c^2
  int nMatched = 0;
  int nNonMatched = 0;
  bool matching = false;
};

















#endif /* Functions_hpp */
