#ifndef TOPTOPLEPHAD_H
#define TOPTOPLEPHAD_H

#include "kinfit.h"

extern void fcnTopTopLepHad(int &npar, double *gin, double &f, double *par, int iflag);

extern double funcTopTopLepHad(double PxLepton,double PyLepton,double PzLepton,double ELepton,
			       double PxBJetLep,double PyBJetLep,double PzBJetLep,double EBJetLep,
			       double PxBJetHad,double PyBJetHad,double PzBJetHad,double EBJetHad,
			       double PxNonBJet1,double PyNonBJet1,double PzNonBJet1,double ENonBJet1,
			       double PxNonBJet2,double PyNonBJet2,double PzNonBJet2,double ENonBJet2,
			       double &chi2WLep, double &chi2TopLep,
			       double &chi2WHad, double &chi2TopHad,
			       double &chi2EtMissX, double &chi2EtMissY,
			       double &chi2BJetPxLep, double &chi2BJetPyLep, double &chi2BJetPzLep, double &chi2BJetELep,
			       double &chi2BJetPxHad, double &chi2BJetPyHad, double &chi2BJetPzHad, double &chi2BJetEHad,
			       double &chi2LeptonPx, double &chi2LeptonPy, double &chi2LeptonPz, double &chi2LeptonE,
			       double &chi2NonBJet1Px, double &chi2NonBJet1Py, double &chi2NonBJet1Pz, double &chi2NonBJet1E,
			       double &chi2NonBJet2Px, double &chi2NonBJet2Py, double &chi2NonBJet2Pz, double &chi2NonBJet2E,
			       double *par);

namespace KINFIT
{         
   class TopTopLepHad : public KINFIT::kfit
     {		                        
      public:
	
	TopTopLepHad();
	virtual ~TopTopLepHad();
	
	void TopTopLepHadRun();
	
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
	
	void EmergencyAssign();
	void CalcNPerm();

	ClassDef(KINFIT::TopTopLepHad,1)
     };
}

#endif
