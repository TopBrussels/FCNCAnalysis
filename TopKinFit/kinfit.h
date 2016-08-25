// KINFIT class definition and includes
#ifndef KINFIT_H
#define KINFIT_H

#define PI 3.14159

#include "TROOT.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "Riostream.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1D.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <memory>

const int NTERMMAX = 20;
const int NPARMAX = 50;
const int NNUMAX = 2;

extern float getProb(TF1 *hPDF,float var);

typedef enum 
{
   TOPTOPLEPLEP,
   TOPLEP,
   TOPTOPLEPHAD,
   TOPTOPLEPHBB,
   TOPHLEPBB
} 
HYPO;

typedef enum 
{
   ELECTRON1_TOPTOPLEPLEP,
   MUON1_TOPTOPLEPLEP,
   ELECTRON2_TOPTOPLEPLEP,
   MUON2_TOPTOPLEPLEP,
   BJET1_TOPTOPLEPLEP,
   BJET2_TOPTOPLEPLEP,
   
   ELECTRON_TOPLEP,
   MUON_TOPLEP,
   BJET_TOPLEP,

   ELECTRON_TOPTOPLEPHAD,
   MUON_TOPTOPLEPHAD,
   BJETLEP_TOPTOPLEPHAD,
   BJETHAD_TOPTOPLEPHAD,
   NONBJET1_TOPTOPLEPHAD,
   NONBJET2_TOPTOPLEPHAD,

   ELECTRON_TOPTOPLEPHBB,
   MUON_TOPTOPLEPHBB,
   BJETLEP_TOPTOPLEPHBB,
   NONBJETHAD_TOPTOPLEPHBB,
   BJET1_TOPTOPLEPHBB,
   BJET2_TOPTOPLEPHBB,

   ELECTRON_TOPHLEPBB,
   MUON_TOPHLEPBB,
   BJETLEP_TOPHLEPBB,
   BJET1_TOPHLEPBB,
   BJET2_TOPHLEPBB,
} 
OBJ;

namespace KINFIT
{
   class kfit
     {
	
      public:
	
	kfit();
	virtual ~kfit();

	// user methods

	void Init(HYPO hypoMode=TOPTOPLEPLEP);
	void Run();
	void Reset();
	
	void TopTopLepLep();
	void TopLep();
	void TopTopLepHad();
	void TopTopLepHbb();
	void TopHLepbb();
	
	void SetBJet(std::vector<float> pt,
		     std::vector<float> eta,
		     std::vector<float> phi,
		     std::vector<float> E);

	void SetNonBJet(std::vector<float> pt,
			std::vector<float> eta,
			std::vector<float> phi,
			std::vector<float> E);

	void SetElectron(std::vector<float> pt,
			 std::vector<float> eta,
			 std::vector<float> phi,
			 std::vector<float> E);

	void SetMuon(std::vector<float> pt,
		     std::vector<float> eta,
		     std::vector<float> phi,		     		     
		     std::vector<float> E);

	void SetMet(float px,
		    float py);

	void SetPDF(std::string obj,std::string fileName,std::string hName);

	void SetWMassBW(float wMass,float wRMS) {*WMassBW = wMass; *WRMSBW = wRMS;};
	
	void SetNToy(int nMC) {NToy_ = nMC;};
	void SetNWRMS(int nMC) {NWRMS_ = nMC;};
	void SetNMetRMS(int nMC) {NMetRMS_ = nMC;};

	void SetNBJetPxRMS(int nMC) {NBJetPxRMS_ = nMC;};
	void SetNBJetPyRMS(int nMC) {NBJetPyRMS_ = nMC;};
	void SetNBJetPzRMS(int nMC) {NBJetPzRMS_ = nMC;};
	void SetNBJetERMS(int nMC) {NBJetERMS_ = nMC;};

	void SetNNonBJetPxRMS(int nMC) {NNonBJetPxRMS_ = nMC;};
	void SetNNonBJetPyRMS(int nMC) {NNonBJetPyRMS_ = nMC;};
	void SetNNonBJetPzRMS(int nMC) {NNonBJetPzRMS_ = nMC;};
	void SetNNonBJetERMS(int nMC) {NNonBJetERMS_ = nMC;};
	
	void SetNElecPxRMS(int nMC) {NElecPxRMS_ = nMC;};
	void SetNElecPyRMS(int nMC) {NElecPyRMS_ = nMC;};
	void SetNElecPzRMS(int nMC) {NElecPzRMS_ = nMC;};
	void SetNElecERMS(int nMC) {NElecERMS_ = nMC;};

	void SetNMuonPxRMS(int nMC) {NMuonPxRMS_ = nMC;};
	void SetNMuonPyRMS(int nMC) {NMuonPyRMS_ = nMC;};
	void SetNMuonPzRMS(int nMC) {NMuonPzRMS_ = nMC;};
	void SetNMuonERMS(int nMC) {NMuonERMS_ = nMC;};
	
	void SetMetXRMS(float met) {*MetXRMS = met;}
	void SetMetYRMS(float met) {*MetYRMS = met;}

	void SetTopMass(float top) {*TopMass = top;}
	void SetTopRMS(float top) {*TopRMS = top;}

	void NoBTag(bool btag) {NoBTag_ = btag;}
	
	float GetDrTopTop(int idxPerm) {return drTopTop_[idxMin_[idxPerm]];};
	float GetMTopTop(int idxPerm) {return mTopTop_[idxMin_[idxPerm]];};
	float GetPtTopTop(int idxPerm) {return ptTopTop_[idxMin_[idxPerm]];};
	float GetPTopTop(int idxPerm) {return pTopTop_[idxMin_[idxPerm]];};
	float GetEtaTopTop(int idxPerm) {return etaTopTop_[idxMin_[idxPerm]];};
	float GetPhiTopTop(int idxPerm) {return phiTopTop_[idxMin_[idxPerm]];};
	float GetRapTopTop(int idxPerm) {return rapTopTop_[idxMin_[idxPerm]];};
	
	float GetDisc(int idx = -1)
	  {
	     if( idx < 0 )
	       return chiPerm_[idxMin_[0]];
	     else if( idx >= NPerm_ )
	       {
		  std::cout << "Discriminator does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return chiPerm_[idxMin_[idx]];
	  };

	float GetDiscTerm(int idxPerm = -1,int idxTerm = -1)
	  {
	     if( idxPerm < 0 || idxTerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation and likelihood term" << std::endl;
		  exit(1);
	       }	     
	     else if( idxTerm >= NTerm_ )
	       {
		  std::cout << "Exceeding max number of terms in the likelihood: " << NTerm_ << std::endl;
		  exit(1);
	       }	  
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Discriminator does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return chiTerm_[idxMin_[idxPerm]][idxTerm];
	  };


	float GetPar(int idxPerm = -1,int idxPar = -1)
	  {
	     if( idxPar < 0 || idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of parameter and permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPar >= NPar_ )
	       {
		  std::cout << "Exceeding max number of parameters: " << NPar_ << std::endl;
		  exit(1);
	       }	  
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return parPerm_[idxMin_[idxPerm]][idxPar];
	  };

	float GetNuPx(int idxPerm = -1,int idxNu = -1)
	  {
	     if( idxPerm < 0 || idxNu < 0 )
	       {
		  std::cout << "Please specify the index of permutation and neutrino" << std::endl;
		  exit(1);
	       }	     
	     else if( idxNu >= NNu_ )
	       {
		  std::cout << "Exceeding number of neutrinos: " << NNu_ << std::endl;
		  exit(1);
	       }	  
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return nuPxPerm_[idxMin_[idxPerm]][idxNu];
	  };

	float GetNuPy(int idxPerm = -1,int idxNu = -1)
	  {
	     if( idxPerm < 0 || idxNu < 0 )
	       {
		  std::cout << "Please specify the index of permutation and neutrino" << std::endl;
		  exit(1);
	       }	     
	     else if( idxNu >= NNu_ )
	       {
		  std::cout << "Exceeding number of neutrinos: " << NNu_ << std::endl;
		  exit(1);
	       }	  
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return nuPyPerm_[idxMin_[idxPerm]][idxNu];
	  };

	float GetNuPz(int idxPerm = -1,int idxNu = -1)
	  {
	     if( idxPerm < 0 || idxNu < 0 )
	       {
		  std::cout << "Please specify the index of permutation and neutrino" << std::endl;
		  exit(1);
	       }	     
	     else if( idxNu >= NNu_ )
	       {
		  std::cout << "Exceeding number of neutrinos: " << NNu_ << std::endl;
		  exit(1);
	       }	  
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return nuPzPerm_[idxMin_[idxPerm]][idxNu];
	  };
	
	float GetMetX(int idxPerm) 
	  {
	     if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }
	     else
	       return MetPx_[idxMin_[idxPerm]];
	  };
	float GetMetY(int idxPerm)
	  {
	     if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }
	     else
	       return MetPy_[idxMin_[idxPerm]];
	  };
	float GetWMassGen(int index=0);
	float GetWMass(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WMass_[idxMin_[idxPerm]][index];
	  };

	float GetWP(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WP_[idxMin_[idxPerm]][index];
	  };

	float GetWPt(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WPt_[idxMin_[idxPerm]][index];
	  };

	float GetWEta(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WEta_[idxMin_[idxPerm]][index];
	  };

	float GetWPhi(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WPhi_[idxMin_[idxPerm]][index];
	  };

	float GetWRap(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return WRap_[idxMin_[idxPerm]][index];
	  };
	
	float GetFitPar(int index=0);

	float GetTopMass(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return TopMass_[idxMin_[idxPerm]][index];
	  };

	float GetTopPt(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return TopPt_[idxMin_[idxPerm]][index];
	  };

	float GetTopP(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return TopP_[idxMin_[idxPerm]][index];
	  };
	
	float GetTopEta(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return TopEta_[idxMin_[idxPerm]][index];
	  };

	float GetTopRap(int idxPerm = -1,int index = 0)
	  {
	     if( idxPerm < 0 )
	       {
		  std::cout << "Please specify the index of permutation" << std::endl;
		  exit(1);
	       }	     
	     else if( idxPerm >= NPerm_ )
	       {
		  std::cout << "Permutation does not exist, max index is " << NPerm_-1 << std::endl;
		  exit(1);
	       }	  
	     else
	       return TopRap_[idxMin_[idxPerm]][index];
	  };
	
	int GetIndex(OBJ objType=BJET1_TOPTOPLEPLEP,int idx = -1);
	
	int GetNPerm() {return NPerm_;};
	int GetNTerm() {return NTerm_;};
	
	void close();
	
      protected:
	
	TRandom3 *rnd;
	
	int hypoMode;
	
	int nBJet;
	std::vector<float> BJetPt;
	std::vector<float> BJetEta;
	std::vector<float> BJetPhi;
	std::vector<float> BJetE;
	std::vector<float> BJetPx;
	std::vector<float> BJetPy;
	std::vector<float> BJetPz;

	int nNonBJet;
	std::vector<float> NonBJetPt;
	std::vector<float> NonBJetEta;
	std::vector<float> NonBJetPhi;
	std::vector<float> NonBJetE;
	std::vector<float> NonBJetPx;
	std::vector<float> NonBJetPy;
	std::vector<float> NonBJetPz;

	int nElectron;
	std::vector<float> ElectronPt;
	std::vector<float> ElectronEta;
	std::vector<float> ElectronPhi;
	std::vector<float> ElectronE;
	std::vector<float> ElectronPx;
	std::vector<float> ElectronPy;
	std::vector<float> ElectronPz;

	int nMuon;
	std::vector<float> MuonPt;
	std::vector<float> MuonEta;
	std::vector<float> MuonPhi;
	std::vector<float> MuonE;
	std::vector<float> MuonPx;
	std::vector<float> MuonPy;
	std::vector<float> MuonPz;

	int nLepton;
	std::vector<float> LeptonPt;
	std::vector<float> LeptonEta;
	std::vector<float> LeptonPhi;
	std::vector<float> LeptonE;
	std::vector<float> LeptonPx;
	std::vector<float> LeptonPy;
	std::vector<float> LeptonPz;
	std::vector<int> LeptonLabel;
	std::vector<int> LeptonIdx;
	
	float *MetPx_;
	float *MetPy_;
	
	float disc_;
	float PxNu1_;
	float PxNu2_;
	float PyNu1_;
	float PyNu2_;
	float PzNu1_;
	float PzNu2_;
	float WMassGen1_;
	float WMassGen2_;
	float **WMass_;
	float **WP_;
	float **WPt_;
	float **WEta_;
	float **WPhi_;
	float **WRap_;
	float **TopMass_;
	float **TopPt_;
	float **TopP_;
	float **TopEta_;
	float **TopRap_;
	float TopMass1_;
	float TopMass2_;
	float MetPx;
	float MetPy;
	
	int NToy_;
	int NWRMS_;
	int NMetRMS_;

	int NBJetPxRMS_;
	int NBJetPyRMS_;
	int NBJetPzRMS_;
	int NBJetERMS_;

	int NNonBJetPxRMS_;
	int NNonBJetPyRMS_;
	int NNonBJetPzRMS_;
	int NNonBJetERMS_;
	
	int NElecPxRMS_;
	int NElecPyRMS_;
	int NElecPzRMS_;
	int NElecERMS_;

	int NMuonPxRMS_;
	int NMuonPyRMS_;
	int NMuonPzRMS_;
	int NMuonERMS_;
	
	int NPerm_;
	int NTerm_;
	int NPar_;
	int NNu_;
	float *chiPerm_;
	float **chiTerm_;
	float **parPerm_;
	float **nuPxPerm_;
	float **nuPyPerm_;
	float **nuPzPerm_;
	int *idxMin_;
	
	float FitPar[100];
	double chis[100];
	
	int *TopTopLepLep_Electron1Idx;
	int *TopTopLepLep_Muon1Idx;
	int *TopTopLepLep_Electron2Idx;
	int *TopTopLepLep_Muon2Idx;
	int *TopTopLepLep_BJet1Idx;
	int *TopTopLepLep_BJet2Idx;

	int *TopLep_ElectronIdx;
	int *TopLep_MuonIdx;
	int *TopLep_BJetIdx;

	int *TopTopLepHad_ElectronIdx;
	int *TopTopLepHad_MuonIdx;
	int *TopTopLepHad_BJetLepIdx;
	int *TopTopLepHad_BJetHadIdx;
	int *TopTopLepHad_NonBJet1Idx;
	int *TopTopLepHad_NonBJet2Idx;

	int *TopTopLepHbb_ElectronIdx;
	int *TopTopLepHbb_MuonIdx;
	int *TopTopLepHbb_BJetLepIdx;
	int *TopTopLepHbb_NonBJetHadIdx;
	int *TopTopLepHbb_BJet1Idx;
	int *TopTopLepHbb_BJet2Idx;

	int *TopHLepbb_ElectronIdx;
	int *TopHLepbb_MuonIdx;
	int *TopHLepbb_BJetLepIdx;
	int *TopHLepbb_BJet1Idx;
	int *TopHLepbb_BJet2Idx;
	
	float *drTopTop_;
	float *mTopTop_;
	float *ptTopTop_;
	float *pTopTop_;
	float *etaTopTop_;
	float *rapTopTop_;
	float *phiTopTop_;

	float getWmassBW(TRandom3 *rnd,float mWmean,float GammaW,float nSigma,double &proc);
	float getWmassFlat(TRandom3 *rnd,float mWmean,float GammaW,float nSigma);
	
	float getProbGaus(TF1 *hPDF,float max,float mean,float sigma,TRandom3 *rnd,float nSigma,double &proc);
	float getProbGausSub(TF1 *hPDF,TRandom3 *rnd,float nSigma);
	
	float BW(float mW,float mWmean,float GammaW);
	
	bool getNuMom(float Wmass,float px_l,float py_l,float pz_l,float E_l,
		      float px_nu,float py_nu,float &pz_nu1,float &pz_nu2);
	
	int getPermIndex(int idx,int *arr);
	
	void sortPermIndex();
	
	void sortPermVector(std::vector<double> vRef,std::vector<double> &vSort);
	
	float getDeltaR(float eta1,float phi1,float eta2,float phi2);
	
	float getEta(float pt,float pz);
	
	float getRap(float E,float pz);

	void checkPDF(TF1* hf,std::string tfname);
	
	float maxPDFBJetPx;
	float meanPDFBJetPx;
	float sigmaPDFBJetPx;

	float maxPDFBJetPy;
	float meanPDFBJetPy;
	float sigmaPDFBJetPy;

	float maxPDFBJetPz;
	float meanPDFBJetPz;
	float sigmaPDFBJetPz;

	float maxPDFBJetE;
	float meanPDFBJetE;
	float sigmaPDFBJetE;

	float maxPDFMetPx;
	float meanPDFMetPx;
	float sigmaPDFMetPx;
	
	float maxPDFMetPy;
	float meanPDFMetPy;
	float sigmaPDFMetPy;
	
	float maxPDFElecPx;
	float meanPDFElecPx;
	float sigmaPDFElecPx;
	
	float maxPDFElecPy;
	float meanPDFElecPy;
	float sigmaPDFElecPy;
	
	float maxPDFElecPz;
	float meanPDFElecPz;
	float sigmaPDFElecPz;
	
	float maxPDFElecE;
	float meanPDFElecE;
	float sigmaPDFElecE;
	
	float maxPDFMuonPx;
	float meanPDFMuonPx;
	float sigmaPDFMuonPx;
	
	float maxPDFMuonPy;
	float meanPDFMuonPy;
	float sigmaPDFMuonPy;
	
	float maxPDFMuonPz;
	float meanPDFMuonPz;
	float sigmaPDFMuonPz;
	
	float maxPDFMuonE;
	float meanPDFMuonE;
	float sigmaPDFMuonE;

	float maxPDFNonBJetPx;
	float meanPDFNonBJetPx;
	float sigmaPDFNonBJetPx;

	float maxPDFNonBJetPy;
	float meanPDFNonBJetPy;
	float sigmaPDFNonBJetPy;

	float maxPDFNonBJetPz;
	float meanPDFNonBJetPz;
	float sigmaPDFNonBJetPz;

	float maxPDFNonBJetE;
	float meanPDFNonBJetE;
	float sigmaPDFNonBJetE;
	
	bool NoBTag_;
	
	int NPERMEVENT;
	
      public:
	
	static std::unique_ptr<double> CHISQ;
	static std::unique_ptr<double> CHISQT1;
	static std::unique_ptr<double> CHISQT2;
	static std::unique_ptr<double> EtMissX;
	static std::unique_ptr<double> EtMissY;
	static std::unique_ptr<double> PxLepton1;
	static std::unique_ptr<double> PyLepton1;
	static std::unique_ptr<double> PzLepton1;
	static std::unique_ptr<double> ELepton1;
	static std::unique_ptr<int> LabelLepton1;
	static std::unique_ptr<double> PxLepton2;
	static std::unique_ptr<double> PyLepton2;
	static std::unique_ptr<double> PzLepton2;
	static std::unique_ptr<double> ELepton2;
	static std::unique_ptr<int> LabelLepton2;
	static std::unique_ptr<double> PxBJet1;
	static std::unique_ptr<double> PyBJet1;
	static std::unique_ptr<double> PzBJet1;
	static std::unique_ptr<double> EBJet1;
	static std::unique_ptr<double> PxBJet2;
	static std::unique_ptr<double> PyBJet2;
	static std::unique_ptr<double> PzBJet2;
	static std::unique_ptr<double> EBJet2;
	static std::unique_ptr<double> PxBJet3;
	static std::unique_ptr<double> PyBJet3;
	static std::unique_ptr<double> PzBJet3;
	static std::unique_ptr<double> EBJet3;
	static std::unique_ptr<double> PxNonBJet1;
	static std::unique_ptr<double> PyNonBJet1;
	static std::unique_ptr<double> PzNonBJet1;
	static std::unique_ptr<double> ENonBJet1;
	static std::unique_ptr<double> PxNonBJet2;
	static std::unique_ptr<double> PyNonBJet2;
	static std::unique_ptr<double> PzNonBJet2;
	static std::unique_ptr<double> ENonBJet2;
	static std::unique_ptr<double> WMassGen1;
	static std::unique_ptr<double> WMassGen2;
	static std::unique_ptr<double> PzNu1;
	static std::unique_ptr<double> PzNu2;
	static std::unique_ptr<double> TopMass1;
	static std::unique_ptr<double> TopMass2;
	static std::unique_ptr<std::vector<double> > FitParam;
	static std::unique_ptr<std::vector<double> > ChiTerm;
	static std::unique_ptr<double> WMassBW;
	static std::unique_ptr<double> WRMSBW;
	static std::unique_ptr<double> MetXRMS;
	static std::unique_ptr<double> MetYRMS;
	static std::unique_ptr<double> TopMass;
	static std::unique_ptr<double> TopRMS;

	static std::unique_ptr<TF1> hPDFTopWMass;
	static std::unique_ptr<TF1> hPDFTopMass;
	static std::unique_ptr<TF1> hPDFTopWHadMass;
	static std::unique_ptr<TF1> hPDFTopHadMass;
	static std::unique_ptr<TF1> hPDFHiggsMass;
	static std::unique_ptr<TF1> hPDFMetPx;
	static std::unique_ptr<TF1> hPDFMetPy;
	static std::unique_ptr<TF1> hPDFBJetPx;
	static std::unique_ptr<TF1> hPDFBJetPy;
	static std::unique_ptr<TF1> hPDFBJetPz;
	static std::unique_ptr<TF1> hPDFBJetE;
	static std::unique_ptr<TF1> hPDFElecPx;
	static std::unique_ptr<TF1> hPDFElecPy;
	static std::unique_ptr<TF1> hPDFElecPz;
	static std::unique_ptr<TF1> hPDFElecE;
	static std::unique_ptr<TF1> hPDFMuonPx;
	static std::unique_ptr<TF1> hPDFMuonPy;
	static std::unique_ptr<TF1> hPDFMuonPz;
	static std::unique_ptr<TF1> hPDFMuonE;
	static std::unique_ptr<TF1> hPDFNonBJetPx;
	static std::unique_ptr<TF1> hPDFNonBJetPy;
	static std::unique_ptr<TF1> hPDFNonBJetPz;
	static std::unique_ptr<TF1> hPDFNonBJetE;
	
      ClassDef(KINFIT::kfit,1)
     };
}

#endif
