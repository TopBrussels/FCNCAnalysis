#ifndef TOPLEP_H
#define TOPLEP_H

#include "kinfit.h"

extern void fcnTopLep(int &npar, double *gin, double &f, double *par, int iflag);

extern double funcTopLep(double PxLepton,double PyLepton,double PzLepton,double ELepton,
			 double PxBJet,double PyBJet,double PzBJet,double EBJet,
			 double &chi2W, double &chi2Top,
			 double &chi2EtMissX, double &chi2EtMissY,
			 double &chi2BJetPx, double &chi2BJetPy, double &chi2BJetPz, double &chi2BJetE,
			 double &chi2LeptonPx, double &chi2LeptonPy, double &chi2LeptonPz, double &chi2LeptonE,
			 double *par);

namespace KINFIT
{         
   class TopLep : public KINFIT::kfit
     {		                        
      public:
	
	TopLep();
	virtual ~TopLep();
	
	void TopLepRun();
	
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
			std::vector<double> &chi2);
	
      private:

	void CalcNPerm();
	
	ClassDef(KINFIT::TopLep,1)
     };
}

#endif
