#ifndef TOPHLEPBB_H
#define TOPHLEPBB_H

#include "kinfit.h"

extern void fcnTopHLepbb(int &npar, double *gin, double &f, double *par, int iflag);

extern double funcTopHLepbb(double PxBJetLep,double PyBJetLep,double PzBJetLep,double EBJetLep,
			    double PxLepton,double PyLepton,double PzLepton,double ELepton,int LabelLepton,
			    double PxBJetHig1,double PyBJetHig1,double PzBJetHig1,double EBJetHig1,
			    double PxBJetHig2,double PyBJetHig2,double PzBJetHig2,double EBJetHig2,
			    double &chi2WLep, double &chi2TopLep,
			    double &chi2Higgs,
			    double &chi2EtMissX, double &chi2EtMissY,
			    double &chi2BJetPxLep, double &chi2BJetPyLep, double &chi2BJetPzLep, double &chi2BJetELep,
			    double &chi2LeptonPx, double &chi2LeptonPy, double &chi2LeptonPz, double &chi2LeptonE,
			    double &chi2BJet1Px, double &chi2BJet1Py, double &chi2BJet1Pz, double &chi2BJet1E,
			    double &chi2BJet2Px, double &chi2BJet2Py, double &chi2BJet2Pz, double &chi2BJet2E,
			    double *par);

namespace KINFIT
{         
   class TopHLepbb : public KINFIT::kfit
     {		                        
      public:
	
	TopHLepbb();
	virtual ~TopHLepbb();
	
	void TopHLepbbRun();
	
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
			std::vector<double> &chi2);
	
      private:
	
	void EmergencyAssign();
	void CalcNPerm();

	ClassDef(KINFIT::TopHLepbb,1)
     };
}

#endif
