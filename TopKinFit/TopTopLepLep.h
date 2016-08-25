#ifndef TOPTOPLEPLEP_H
#define TOPTOPLEPLEP_H

#include "kinfit.h"

extern void fcnTopTopLepLep(int &npar, double *gin, double &f, double *par, int iflag);

extern double funcTopTopLepLep(double PxLepton1,double PyLepton1,double PzLepton1,double ELepton1,
			       double PxLepton2,double PyLepton2,double PzLepton2,double ELepton2,
			       double PxBJet1,double PyBJet1,double PzBJet1,double EBJet1,
			       double PxBJet2,double PyBJet2,double PzBJet2,double EBJet2,
			       double &chi2W1, double &chi2W2, double &chi2Top1, double &chi2Top2,
			       double &chi2EtMissX, double &chi2EtMissY,
			       double &chi2BJet1Px, double &chi2BJet1Py, double &chi2BJet1Pz, double &chi2BJet1E,
			       double &chi2BJet2Px, double &chi2BJet2Py, double &chi2BJet2Pz, double &chi2BJet2E,
			       double &chi2Lepton1Px, double &chi2Lepton1Py, double &chi2Lepton1Pz, double &chi2Lepton1E,
			       double &chi2Lepton2Px, double &chi2Lepton2Py, double &chi2Lepton2Pz, double &chi2Lepton2E,
			       double *par);

namespace KINFIT
{         
   class TopTopLepLep : public KINFIT::kfit
     {		                        
      public:
	
	TopTopLepLep();
	virtual ~TopTopLepLep();
	
	void TopTopLepLepRun();
	
	void fit(double *chis,
		 double *par);
	
	void calcVar(int iPerm);
	
	void calcNuGrid(std::vector<double> &p1,
			std::vector<double> &p2,
			std::vector<double> &p3,
			std::vector<double> &p4,
			std::vector<double> &p5,
			std::vector<double> &p6,
			std::vector<double> &p7,
			std::vector<double> &p8,
			std::vector<double> &p9,
			std::vector<double> &p10,
			std::vector<double> &p11,
			std::vector<double> &p12,
			std::vector<double> &p13,
			std::vector<double> &p14,
			std::vector<double> &p15,
			std::vector<double> &p16,
			std::vector<double> &p17,
			std::vector<double> &p18,
			std::vector<double> &p19,
			std::vector<double> &p20,
			std::vector<double> &p21,
			std::vector<double> &p22,
			std::vector<double> &p23,
			std::vector<double> &p24,
			std::vector<double> &chi2);
	
      private:

	void CalcNPerm();
	
	ClassDef(KINFIT::TopTopLepLep,1)
     };
}

#endif
