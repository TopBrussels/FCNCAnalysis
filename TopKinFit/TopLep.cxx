// Top leptonic decay: t/tbar -> (W->lnu) b
#include "TopLep.h"

ClassImp(KINFIT::TopLep)

void fcnTopLep(int &npar, double *gin, double &f, double *par, int iflag);

double funcTopLep(double PxLepton,double PyLepton,double PzLepton,double ELepton,int LabelLepton,
		  double PxBJet,double PyBJet,double PzBJet,double EBJet,
		  double &chi2W, double &chi2Top,
		  double &chi2EtMissX, double &chi2EtMissY,
		  double &chi2BJetPx, double &chi2BJetPy, double &chi2BJetPz, double &chi2BJetE,
		  double &chi2LeptonPx, double &chi2LeptonPy, double &chi2LeptonPz, double &chi2LeptonE,
		  double *par, double &proc);

KINFIT::TopLep::TopLep()
{
}

KINFIT::TopLep::~TopLep()
{
}

void KINFIT::TopLep::TopLepRun()
{
   checkPDF(hPDFTopWMass.get(),"PDFTopWMass");
   checkPDF(hPDFTopMass.get(),"PDFTopMass");
   
   checkPDF(hPDFMetPx.get(),"PDFMetPx");
   checkPDF(hPDFMetPy.get(),"PDFMetPy");
   
   checkPDF(hPDFBJetPx.get(),"PDFBJetPx");
   checkPDF(hPDFBJetPy.get(),"PDFBJetPy");
   checkPDF(hPDFBJetPz.get(),"PDFBJetPz");
   checkPDF(hPDFBJetE.get(),"PDFBJetE");

   checkPDF(hPDFElecPx.get(),"PDFElecPx");
   checkPDF(hPDFElecPy.get(),"PDFElecPy");
   checkPDF(hPDFElecPz.get(),"PDFElecPz");
   checkPDF(hPDFElecE.get(),"PDFElecE");

   checkPDF(hPDFMuonPx.get(),"PDFMuonPx");
   checkPDF(hPDFMuonPy.get(),"PDFMuonPy");
   checkPDF(hPDFMuonPz.get(),"PDFMuonPz");
   checkPDF(hPDFMuonE.get(),"PDFMuonE");

   checkPDF(hPDFBJetPx.get(),"PDFBJetPx");
   checkPDF(hPDFBJetPy.get(),"PDFBJetPy");
   checkPDF(hPDFBJetPz.get(),"PDFBJetPz");
   checkPDF(hPDFBJetE.get(),"PDFBJetE");
   
   disc_ = 10E+10;
   NPerm_ = 0;
   NTerm_ = 12;
   NPar_ = 100;
   NNu_ = 1;
   
   for(int i=0;i<100;i++) FitPar[i] = 10E+10;
   
   PxNu1_ = 10E+10;
   PyNu1_ = 10E+10;   
   PzNu1_ = 10E+10;

   WMassGen1_ = 10E+10;
   
   TopMass1_ = 10E+10;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}

   if( TopLep_ElectronIdx ) {delete[] TopLep_ElectronIdx; TopLep_ElectronIdx = 0;}
   if( TopLep_MuonIdx ) {delete[] TopLep_MuonIdx; TopLep_MuonIdx = 0;}
   if( TopLep_BJetIdx ) {delete[] TopLep_BJetIdx; TopLep_BJetIdx = 0;}
   
   if( chiPerm_ ) {delete[] chiPerm_; chiPerm_ = 0;}
   
   if( MetPx_ ) {delete[] MetPx_; MetPx_ = 0;}
   if( MetPy_ ) {delete[] MetPy_; MetPy_ = 0;}

   if( parPerm_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] parPerm_[i]; delete[] parPerm_; parPerm_ = 0;}
   
   if( chiTerm_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] chiTerm_[i]; delete[] chiTerm_; chiTerm_ = 0;}

   if( WMass_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WMass_[i]; delete[] WMass_; WMass_ = 0;}
   if( WP_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WP_[i]; delete[] WP_; WP_ = 0;}
   if( WPt_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WPt_[i]; delete[] WPt_; WPt_ = 0;}
   if( WEta_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WEta_[i]; delete[] WEta_; WEta_ = 0;}
   if( WPhi_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WPhi_[i]; delete[] WPhi_; WPhi_ = 0;}
   if( WRap_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] WRap_[i]; delete[] WRap_; WRap_ = 0;}
   if( TopMass_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] TopMass_[i]; delete[] TopMass_; TopMass_ = 0;}
   if( TopPt_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] TopPt_[i]; delete[] TopPt_; TopPt_ = 0;}
   if( TopP_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] TopP_[i]; delete[] TopP_; TopP_ = 0;}
   if( TopEta_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] TopEta_[i]; delete[] TopEta_; TopEta_ = 0;}
   if( TopRap_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] TopRap_[i]; delete[] TopRap_; TopRap_ = 0;}
   if( nuPxPerm_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] nuPxPerm_[i]; delete[] nuPxPerm_; nuPxPerm_ = 0;}
   if( nuPyPerm_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] nuPyPerm_[i]; delete[] nuPyPerm_; nuPyPerm_ = 0;}
   if( nuPzPerm_ ) {for(int i=0;i<NPERMEVENT;i++) delete[] nuPzPerm_[i]; delete[] nuPzPerm_; nuPzPerm_ = 0;}

   CalcNPerm();
   
   idxMin_ = new int[NPERMEVENT];
   
   TopLep_ElectronIdx = new int[NPERMEVENT];
   TopLep_MuonIdx = new int[NPERMEVENT];
   TopLep_BJetIdx = new int[NPERMEVENT];
   chiPerm_ = new float[NPERMEVENT];
   
   MetPx_ = new float[NPERMEVENT];
   MetPy_ = new float[NPERMEVENT];

   parPerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) parPerm_[i] = new float[NPARMAX];
   
   chiTerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) chiTerm_[i] = new float[NTERMMAX];

   WMass_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WMass_[i] = new float[NNUMAX];
   WP_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WP_[i] = new float[NNUMAX];
   WPt_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPt_[i] = new float[NNUMAX];
   WEta_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WEta_[i] = new float[NNUMAX];
   WPhi_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WPhi_[i] = new float[NNUMAX];
   WRap_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) WRap_[i] = new float[NNUMAX];
   TopMass_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopMass_[i] = new float[NNUMAX];
   TopPt_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopPt_[i] = new float[NNUMAX];
   TopP_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopP_[i] = new float[NNUMAX];
   TopEta_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopEta_[i] = new float[NNUMAX];
   TopRap_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) TopRap_[i] = new float[NNUMAX];
   nuPxPerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPxPerm_[i] = new float[NNUMAX];
   nuPyPerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPyPerm_[i] = new float[NNUMAX];
   nuPzPerm_ = new float*[NPERMEVENT]; for(int i=0;i<NPERMEVENT;i++) nuPzPerm_[i] = new float[NNUMAX];
   
   for(int i=0;i<NPERMEVENT;i++)
     {	
	TopLep_ElectronIdx[i] = -1;
	TopLep_MuonIdx[i] = -1;
	TopLep_BJetIdx[i] = -1;
	chiPerm_[i] = 10E+10;

	MetPx_[i] = 10E+10;
	MetPy_[i] = 10E+10;
	
	for(int ipar=0;ipar<NPARMAX;ipar++)
	  {	     
	     parPerm_[i][ipar] = 10E+10;
	  }	

	for(int it=0;it<NTERMMAX;it++)
	  {	     
	     chiTerm_[i][it] = 10E+10;
	  }	
	
	for(int inu=0;inu<NNUMAX;inu++)
	  {	     
	     WMass_[i][inu] = 10E+10;
	     WP_[i][inu] = 10E+10;
	     WPt_[i][inu] = 10E+10;
	     WEta_[i][inu] = 10E+10;
	     WPhi_[i][inu] = 10E+10;
	     WRap_[i][inu] = 10E+10;
	     TopMass_[i][inu] = 10E+10;
	     TopPt_[i][inu] = 10E+10;
	     TopP_[i][inu] = 10E+10;
	     TopEta_[i][inu] = 10E+10;
	     TopRap_[i][inu] = 10E+10;
	     nuPxPerm_[i][inu] = 10E+10;
	     nuPyPerm_[i][inu] = 10E+10;
	     nuPzPerm_[i][inu] = 10E+10;
	  }	
     }   
   
   for(int il=0;il<nLepton;il++)
     {	
	*KINFIT::kfit::EtMissX = MetPx;
	*KINFIT::kfit::EtMissY = MetPy;
	
	*KINFIT::kfit::PxLepton1 = LeptonPx[il];
	*KINFIT::kfit::PyLepton1 = LeptonPy[il];
	*KINFIT::kfit::PzLepton1 = LeptonPz[il];
	*KINFIT::kfit::ELepton1 = LeptonE[il];
	*KINFIT::kfit::LabelLepton1 = LeptonLabel[il];
	     
	int label_l = LeptonLabel[il];
	int idx_l = LeptonIdx[il];
		       
	for(int ib=0;ib<nBJet;ib++)
	  {
	     *KINFIT::kfit::EBJet1 = BJetE[ib];
	     *KINFIT::kfit::PxBJet1 = BJetPx[ib];
	     *KINFIT::kfit::PyBJet1 = BJetPy[ib];
	     *KINFIT::kfit::PzBJet1 = BJetPz[ib];
	     
	     chiPerm_[NPerm_] = 10E+10;
			    
	     std::vector<double> vp1;
	     std::vector<double> vp2;
	     std::vector<double> vp3;
	     std::vector<double> vp4;
	     std::vector<double> vp5;
	     std::vector<double> vp6;
	     std::vector<double> vp7;
	     std::vector<double> vp8;
	     std::vector<double> vp9;
	     std::vector<double> vp10;
	     std::vector<double> vp11;
	     std::vector<double> vp12;
	     std::vector<double> vchi2;
	     
	     calcNuGrid(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,vp9,vp10,
			vp11,vp12,
			vchi2);
	     
	     int nSolutions = vchi2.size();
	     disc_ = 10E+11;

	     for(int ig=0;ig<nSolutions;ig++)
	       {		
		  double par[NPARMAX];
		  par[0] = vp1[ig];
		  par[1] = vp2[ig];
		  par[2] = vp3[ig];
		  par[3] = vp4[ig];
		  par[4] = vp5[ig];
		  par[5] = vp6[ig];
		  par[6] = vp7[ig];
		  par[7] = vp8[ig];
		  par[8] = vp9[ig];
		  par[9] = vp10[ig];
		  par[10] = vp11[ig];
		  par[11] = vp12[ig];
		  
		  *WMassGen1 = par[3];
		       
		  fit(chis,par);
		  
		  float res = (*CHISQ < vchi2[ig]) ? *CHISQ : vchi2[ig];
		  
		  if( res < disc_ )
		    {
		       disc_ = res;
		       PxNu1_ = par[0];
		       PyNu1_ = par[1];
		       PzNu1_ = *PzNu1;
		       WMassGen1_ = *WMassGen1;
		       TopMass1_ = *TopMass1;
		       MetPx_[NPerm_] = par[0];
		       MetPy_[NPerm_] = par[1];
				 
		       FitPar[0] = (*FitParam).at(0);
		       FitPar[1] = (*FitParam).at(1);
		       FitPar[2] = (*FitParam).at(2);
		       FitPar[3] = (*FitParam).at(3);
		       FitPar[4] = (*FitParam).at(4);
		       FitPar[5] = (*FitParam).at(5);
		       FitPar[6] = (*FitParam).at(6);
		       FitPar[7] = (*FitParam).at(7);
		       FitPar[8] = (*FitParam).at(8);
		       FitPar[9] = (*FitParam).at(9);
		       FitPar[10] = (*FitParam).at(10);
		       FitPar[11] = (*FitParam).at(11);
			    
		       TopLep_BJetIdx[NPerm_] = ib;
			    
		       if( label_l == 0 )
			 {
			    TopLep_ElectronIdx[NPerm_] = idx_l;
			    TopLep_MuonIdx[NPerm_] = -1;
			 }							  
		       else if( label_l == 1 )
			 {
			    TopLep_MuonIdx[NPerm_] = idx_l;
			    TopLep_ElectronIdx[NPerm_] = -1;
			 }							  

		       chiPerm_[NPerm_] = disc_;
		       
		       parPerm_[NPerm_][0] = FitPar[0];
		       parPerm_[NPerm_][1] = FitPar[1];
		       parPerm_[NPerm_][2] = FitPar[2];
		       parPerm_[NPerm_][3] = FitPar[3];
		       parPerm_[NPerm_][4] = FitPar[4];
		       parPerm_[NPerm_][5] = FitPar[5];
		       parPerm_[NPerm_][6] = FitPar[6];
		       parPerm_[NPerm_][7] = FitPar[7];
		       parPerm_[NPerm_][8] = FitPar[8];
		       parPerm_[NPerm_][9] = FitPar[9];
		       parPerm_[NPerm_][10] = FitPar[10];
		       parPerm_[NPerm_][11] = FitPar[11];
			    
		       nuPxPerm_[NPerm_][0] = PxNu1_;
		       nuPyPerm_[NPerm_][0] = PyNu1_;
		       nuPzPerm_[NPerm_][0] = PzNu1_;

		       WMass_[NPerm_][0] = FitPar[3];
		       
		       chiTerm_[NPerm_][0] = (*ChiTerm).at(0);
		       chiTerm_[NPerm_][1] = (*ChiTerm).at(1);
		       chiTerm_[NPerm_][2] = (*ChiTerm).at(2);
		       chiTerm_[NPerm_][3] = (*ChiTerm).at(3);
		       chiTerm_[NPerm_][4] = (*ChiTerm).at(4);
		       chiTerm_[NPerm_][5] = (*ChiTerm).at(5);		       
		    }		
	       }
	     
	     NPerm_++;
	  }	
     } // end permutations
   
   // choose the best permutation
   for(int ip=0;ip<NPerm_;ip++)
     {
	calcVar(ip);
     }   
}

void KINFIT::TopLep::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int il=0;il<nLepton;il++)
     {	
	for(int ib=0;ib<nBJet;ib++)
	  {
	     NPERMEVENT++;
	  }
     }   
}

double funcTopLep(double PxLepton,double PyLepton,double PzLepton,double ELepton,int LabelLepton,
		  double PxBJet,double PyBJet,double PzBJet,double EBJet,
		  double &chi2W, double &chi2Top,
		  double &chi2EtMissX, double &chi2EtMissY,
		  double &chi2BJetPx, double &chi2BJetPy, double &chi2BJetPz, double &chi2BJetE,
		  double &chi2LeptonPx, double &chi2LeptonPy, double &chi2LeptonPz, double &chi2LeptonE,
		  double *par, double &proc)
{
   clock_t tStart = clock();
   
   double val = 10E+10;

   float PxNu = par[0];
   float PyNu = par[1];
   
   float a = sqrt(par[8]*par[8]+par[9]*par[9]);
   float b = par[10];
   float d = sqrt(PxNu*PxNu+PyNu*PyNu);
   float f = par[11];

   float c = par[3]*par[3]/2+par[8]*PxNu+par[9]*PyNu;
   
   float rac = c*c*b*b-a*a*(d*d*f*f-c*c);
   
   float racAbs = fabs(rac);
   
   float PzNu = 0.;
   
   if( rac >= 0 )
     {
	PzNu = (c*b+par[2]*sqrt(rac))/a/a;
     }   
   
   *KINFIT::kfit::PzNu1 = PzNu;
   
   float ENu = sqrt(PxNu*PxNu+PyNu*PyNu+PzNu*PzNu);
   
   float totPx = par[8]+PxNu+par[4];
   float totPy = par[9]+PyNu+par[5];
   float totPz = par[10]+PzNu+par[6];
   float totE = par[11]+ENu+par[7];
   
   float mW = pow(ENu+par[11],2)-pow(PxNu+par[8],2)-pow(PyNu+par[9],2)-pow(PzNu+par[10],2);
   float mtop = totE*totE-totPx*totPx-totPy*totPy-totPz*totPz;

   chi2W = 10E+10;
   chi2Top = 10E+10;
   chi2EtMissX = 10E+10;
   chi2EtMissY = 10E+10;
   chi2BJetPx = 10E+10;
   chi2BJetPy = 10E+10;
   chi2BJetPz = 10E+10;
   chi2BJetE = 10E+10;
   chi2LeptonPx = 10E+10;
   chi2LeptonPy = 10E+10;
   chi2LeptonPz = 10E+10;
   chi2LeptonE = 10E+10;
   
   if( mtop >= 0 && rac >= 0 && mW >= 0 )
     {
	mtop = sqrt(mtop);
	mW = sqrt(mW);

	*KINFIT::kfit::TopMass1 = mtop;
	
	float mWProb = getProb(KINFIT::kfit::hPDFTopWMass.get(),mW);

	float mTopProb = getProb(KINFIT::kfit::hPDFTopMass.get(),mtop);
	
	float MetPxProb = getProb(KINFIT::kfit::hPDFMetPx.get(),par[0]-*KINFIT::kfit::EtMissX);
	float MetPyProb = getProb(KINFIT::kfit::hPDFMetPy.get(),par[1]-*KINFIT::kfit::EtMissY);

	float BJetPxProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[4]-PxBJet);
	float BJetPyProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[5]-PyBJet);
	float BJetPzProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[6]-PzBJet);
	float BJetEProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[7]-EBJet);

	float LeptonPxProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPx.get(),par[8]-PxLepton) : getProb(KINFIT::kfit::hPDFMuonPx.get(),par[8]-PxLepton);
	float LeptonPyProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPy.get(),par[9]-PyLepton) : getProb(KINFIT::kfit::hPDFMuonPy.get(),par[9]-PyLepton);
	float LeptonPzProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPz.get(),par[10]-PzLepton) : getProb(KINFIT::kfit::hPDFMuonPz.get(),par[10]-PzLepton);
	float LeptonEProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecE.get(),par[11]-ELepton) : getProb(KINFIT::kfit::hPDFMuonE.get(),par[11]-ELepton);

	chi2W = mWProb;
	chi2Top = mTopProb;
	chi2EtMissX = MetPxProb;
	chi2EtMissY = MetPyProb;
	chi2BJetPx = BJetPxProb;
	chi2BJetPy = BJetPyProb;
	chi2BJetPz = BJetPzProb;
	chi2BJetE = BJetEProb;
	chi2LeptonPx = LeptonPxProb;
	chi2LeptonPy = LeptonPyProb;
	chi2LeptonPz = LeptonPzProb;
	chi2LeptonE = LeptonEProb;

//	std::cout << chi2W << " " << chi2Top << " " << chi2EtMissX << " " << chi2EtMissY << std::endl; 

	chi2W = (chi2W < 10E-10) ? 10E-10 : chi2W;
	chi2Top = (chi2Top < 10E-10) ? 10E-10 : chi2Top;
	
	val = -2*log(
		     chi2W*
		     chi2Top
//		     chi2EtMissX*
//		     chi2EtMissY
		    );
     }

   proc = (double)(clock()-tStart);
   
   return val;
}

void KINFIT::TopLep::calcVar(int iPerm)
{
   if( chiPerm_[iPerm] > 10E+9 ) return;
   
   int idxElec = TopLep_ElectronIdx[iPerm];
   int idxMuon = TopLep_MuonIdx[iPerm];   
   int idxBJet = TopLep_BJetIdx[iPerm];
   
   float Lepton1Px = 0;
   float Lepton1Py = 0;
   float Lepton1Pz = 0;
   float Lepton1E = 0;

   float BJet1Px = 0;
   float BJet1Py = 0;
   float BJet1Pz = 0;
   float BJet1E = 0;
   
   if( idxElec >= 0 )
     {
	Lepton1Px = ElectronPx[idxElec];
	Lepton1Py = ElectronPy[idxElec];
	Lepton1Pz = ElectronPz[idxElec];
	Lepton1E = ElectronE[idxElec];
     }   
   else if( idxMuon >= 0 )
     {
	Lepton1Px = MuonPx[idxMuon];
	Lepton1Py = MuonPy[idxMuon];
	Lepton1Pz = MuonPz[idxMuon];
	Lepton1E = MuonE[idxMuon];
     }   

   if( idxBJet >= 0 )
     {
	BJet1Px = BJetPx[idxBJet];
	BJet1Py = BJetPy[idxBJet];
	BJet1Pz = BJetPz[idxBJet];
	BJet1E = BJetE[idxBJet];
     }   

   float nuPx = nuPxPerm_[iPerm][0];
   float nuPy = nuPyPerm_[iPerm][0];
   float nuPz = nuPzPerm_[iPerm][0];   
   float nuE = sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz);

   // W
   float WE = nuE+Lepton1E;
   float WPx = nuPx+Lepton1Px;
   float WPy = nuPy+Lepton1Py;
   float WPz = nuPz+Lepton1Pz;

   float WPt = sqrt(WPx*WPx+WPy*WPy);
   float WEta = getEta(WPt,WPz);
   float WRap = getRap(WE,WPz);
   float WPhi = atan(WPx/WPy);
   
   // top
   float topE = WE+BJet1E;
   float topPx = WPx+BJet1Px;
   float topPy = WPy+BJet1Py;
   float topPz = WPz+BJet1Pz;
   
   float topPt = sqrt(topPx*topPx+topPy*topPy);
   float topEta = getEta(topPt,topPz);
   float topRap = getRap(topE,topPz);
   float topPhi = atan(topPx/topPy);

   TopMass_[iPerm][0] = sqrt(topE*topE-topPx*topPx-topPy*topPy-topPz*topPz);
   TopPt_[iPerm][0] = sqrt(topPx*topPx+topPy*topPy);
   TopP_[iPerm][0] = sqrt(topPx*topPx+topPy*topPy+topPz*topPz);
   TopEta_[iPerm][0] = getEta(TopPt_[iPerm][0],topPz);
   TopRap_[iPerm][0] = getRap(topE,topPz);

   WPt_[iPerm][0] = sqrt(WPx*WPx+WPy*WPy);
   WP_[iPerm][0] = sqrt(WPx*WPx+WPy*WPy+WPz*WPz);
   WEta_[iPerm][0] = getEta(WPt_[iPerm][0],WPz);
   WRap_[iPerm][0] = getRap(WE,WPz);
   WPhi_[iPerm][0] = atan(WPx/WPy);
}

void fcnTopLep(int &npar, double *gin, double &f, double *par, int iflag)
{
   double chi2W;
   double chi2Top;
   double chi2EtMissX;
   double chi2EtMissY;
   double chi2BJetPx;
   double chi2BJetPy;
   double chi2BJetPz;
   double chi2BJetE;
   double chi2LeptonPx;
   double chi2LeptonPy;
   double chi2LeptonPz;
   double chi2LeptonE;

   double proc;
   double lh = funcTopLep(*KINFIT::kfit::PxLepton1,
			  *KINFIT::kfit::PyLepton1,
			  *KINFIT::kfit::PzLepton1,
			  *KINFIT::kfit::ELepton1,
			  *KINFIT::kfit::LabelLepton1,
			  *KINFIT::kfit::PxBJet1,
			  *KINFIT::kfit::PyBJet1,
			  *KINFIT::kfit::PzBJet1,
			  *KINFIT::kfit::EBJet1,
			  chi2W, chi2Top, chi2EtMissX, chi2EtMissY,
			  chi2BJetPx, chi2BJetPy, chi2BJetPz, chi2BJetE,
			  chi2LeptonPx, chi2LeptonPy, chi2LeptonPz, chi2LeptonE,
			  par,proc);
   
   if( iflag == 3 )
     {
	*KINFIT::kfit::CHISQ = lh;
	(*KINFIT::kfit::FitParam).push_back(par[0]);
	(*KINFIT::kfit::FitParam).push_back(par[1]);
	(*KINFIT::kfit::FitParam).push_back(par[2]);
	(*KINFIT::kfit::FitParam).push_back(par[3]);
	(*KINFIT::kfit::FitParam).push_back(par[4]);
	(*KINFIT::kfit::FitParam).push_back(par[5]);
	(*KINFIT::kfit::FitParam).push_back(par[6]);
	(*KINFIT::kfit::FitParam).push_back(par[7]);
	(*KINFIT::kfit::FitParam).push_back(par[8]);
	(*KINFIT::kfit::FitParam).push_back(par[9]);
	(*KINFIT::kfit::FitParam).push_back(par[10]);
	(*KINFIT::kfit::FitParam).push_back(par[11]);
	
	(*KINFIT::kfit::ChiTerm).push_back(chi2W);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Top);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissX);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissY);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPx);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPy);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetE);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPx);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPy);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonE);
     }
   
   f = lh;
}

void KINFIT::TopLep::fit(double *chis,
			 double *par)
{
   (*KINFIT::kfit::FitParam).clear();
   (*KINFIT::kfit::ChiTerm).clear();
   *KINFIT::kfit::CHISQ = 10E+10;
   
   double perr[100];
   
   TMinuit *gMinuit = new TMinuit(100);
   gMinuit->SetFCN(fcnTopLep);
   gMinuit->SetPrintLevel(-1);
   
   Double_t arglist[100];
   Int_t ierflg = 0;
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
   
   arglist[0] = 5000;
   arglist[1] = 0.1;
   
   gMinuit->mnparm(0,"Etx",par[0],pow(10.,-2),-1000.,1000.,ierflg);
   gMinuit->mnparm(1,"Ety",par[1],pow(10.,-2),-1000.,1000.,ierflg);
   gMinuit->mnparm(2,"Sign",par[2],pow(10.,-2),-1.,1.,ierflg);
   gMinuit->mnparm(3,"mW",par[3],pow(10.,-2),1.,1000.,ierflg);
   
   gMinuit->mnparm(4,"BJetPx",par[4],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(5,"BJetPy",par[5],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(6,"BJetPz",par[6],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(7,"BJetE",par[7],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(8,"LeptonPx",par[8],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(9,"LeptonPy",par[9],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(10,"LeptonPz",par[10],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(11,"LeptonE",par[11],pow(10.,-1),-1000.,1000.,ierflg);

   gMinuit->FixParameter(0);
   gMinuit->FixParameter(1);
   
   gMinuit->FixParameter(2);
   
   gMinuit->FixParameter(3);
   
   gMinuit->FixParameter(4);
   gMinuit->FixParameter(5);
   gMinuit->FixParameter(6);
   gMinuit->FixParameter(7);
   gMinuit->FixParameter(8);
   gMinuit->FixParameter(9);
   gMinuit->FixParameter(10);
   gMinuit->FixParameter(11);

   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
//   gMinuit->mnexcm("MINOS", arglist ,2,ierflg);

   arglist[0] = 3;
   gMinuit->mnexcm("CALL FCN",arglist,1,ierflg);
   
   double amin,edm,errdef;
   int nvpar,nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
//   std::cout << icstat << std::endl;
   chis[0] = amin;
   chis[1] = gMinuit->GetNumFreePars();
   
   gMinuit->GetParameter(0,par[0],perr[0]);
   gMinuit->GetParameter(1,par[1],perr[1]);
   gMinuit->GetParameter(2,par[2],perr[2]);
   gMinuit->GetParameter(3,par[3],perr[3]);
   gMinuit->GetParameter(4,par[4],perr[4]);
   gMinuit->GetParameter(5,par[5],perr[5]);
   gMinuit->GetParameter(6,par[6],perr[6]);
   gMinuit->GetParameter(7,par[7],perr[7]);
   gMinuit->GetParameter(8,par[8],perr[8]);
   gMinuit->GetParameter(9,par[9],perr[9]);
   gMinuit->GetParameter(10,par[10],perr[10]);
   gMinuit->GetParameter(11,par[11],perr[11]);

   delete gMinuit;
}

void KINFIT::TopLep::calcNuGrid(std::vector<double> &vp1,
				std::vector<double> &vp2,
				std::vector<double> &vp3,
				std::vector<double> &vp4,
				std::vector<double> &vp5,
				std::vector<double> &vp6,
				std::vector<double> &vp7,
				std::vector<double> &vp8,
				std::vector<double> &vp9,
				std::vector<double> &vp10,
				std::vector<double> &vp11,
				std::vector<double> &vp12,
				std::vector<double> &vchi2)
{
   double par[NPARMAX];

   clock_t tStart = clock();
   
   bool bp = 0;

   std::vector<double> TimePDFBJetPx;
   std::vector<double> TimePDFBJetPy;
   std::vector<double> TimePDFBJetPz;
   std::vector<double> TimePDFBJetE;
   std::vector<double> TimePDFMetPx;
   std::vector<double> TimePDFMetPy;
   std::vector<double> TimePDFWMass;
   
   double proc = 0.;
   
   for(int it=0;it<NToy_;it++)
     {	
	if( bp ) break;

	do 
	  {	     
	     par[8] = getProbGaus(hPDFBJetPx.get(),maxPDFBJetPx,meanPDFBJetPx,sigmaPDFBJetPx,rnd,NBJetPxRMS_,proc)+*KINFIT::kfit::PxBJet1;	
	     par[9] = getProbGaus(hPDFBJetPy.get(),maxPDFBJetPy,meanPDFBJetPy,sigmaPDFBJetPy,rnd,NBJetPyRMS_,proc)+*KINFIT::kfit::PyBJet1;		  
	     par[10] = getProbGaus(hPDFBJetPz.get(),maxPDFBJetPz,meanPDFBJetPz,sigmaPDFBJetPz,rnd,NBJetPzRMS_,proc)+*KINFIT::kfit::PzBJet1;
	     par[11] = getProbGaus(hPDFBJetE.get(),maxPDFBJetE,meanPDFBJetE,sigmaPDFBJetE,rnd,NBJetERMS_,proc)+*KINFIT::kfit::EBJet1;
	  } while( (par[11]*par[11]-par[8]*par[8]-par[9]*par[9]-par[10]*par[10]) < 0. );	
	     
	par[0] = getProbGaus(hPDFMetPx.get(),maxPDFMetPx,meanPDFMetPx,sigmaPDFMetPx,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissX;	
	par[1] = getProbGaus(hPDFMetPy.get(),maxPDFMetPy,meanPDFMetPy,sigmaPDFMetPy,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissY;
	
	par[3] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);

	if( *KINFIT::kfit::LabelLepton1 == 0 )
	  {
	     do
	       {		  
		  par[4] = getProbGaus(hPDFElecPx.get(),maxPDFElecPx,meanPDFElecPx,sigmaPDFElecPx,rnd,NElecPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[5] = getProbGaus(hPDFElecPy.get(),maxPDFElecPy,meanPDFElecPy,sigmaPDFElecPy,rnd,NElecPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[6] = getProbGaus(hPDFElecPz.get(),maxPDFElecPz,meanPDFElecPz,sigmaPDFElecPz,rnd,NElecPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[7] = getProbGaus(hPDFElecE.get(),maxPDFElecE,meanPDFElecE,sigmaPDFElecE,rnd,NElecERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[7]*par[7]-par[4]*par[4]-par[5]*par[5]-par[6]*par[6]) < 0. );	     
	  }
	else
	  {
	     do
	       {		  
		  par[4] = getProbGaus(hPDFMuonPx.get(),maxPDFMuonPx,meanPDFMuonPx,sigmaPDFMuonPx,rnd,NMuonPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[5] = getProbGaus(hPDFMuonPy.get(),maxPDFMuonPy,meanPDFMuonPy,sigmaPDFMuonPy,rnd,NMuonPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[6] = getProbGaus(hPDFMuonPz.get(),maxPDFMuonPz,meanPDFMuonPz,sigmaPDFMuonPz,rnd,NMuonPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[7] = getProbGaus(hPDFMuonE.get(),maxPDFMuonE,meanPDFMuonE,sigmaPDFMuonE,rnd,NMuonERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[7]*par[7]-par[4]*par[4]-par[5]*par[5]-par[6]*par[6]) < 0. );	     
	  }

	for(int is1=-1;is1<=1;is1++)
	  {	
	     if( bp ) break;
	     
	     if( is1 == 0 ) continue;
	     
	     par[2] = is1;
	     
	     double chi2W;
	     double chi2Top;
	     double chi2EtMissX;
	     double chi2EtMissY;
	     double chi2BJetPx;
	     double chi2BJetPy;
	     double chi2BJetPz;
	     double chi2BJetE;
	     double chi2LeptonPx;
	     double chi2LeptonPy;
	     double chi2LeptonPz;
	     double chi2LeptonE;
		  
	     double lh = funcTopLep(*KINFIT::kfit::PxLepton1,
				    *KINFIT::kfit::PyLepton1,
				    *KINFIT::kfit::PzLepton1,
				    *KINFIT::kfit::ELepton1,
				    *KINFIT::kfit::LabelLepton1,
				    *KINFIT::kfit::PxBJet1,
				    *KINFIT::kfit::PyBJet1,
				    *KINFIT::kfit::PzBJet1,
				    *KINFIT::kfit::EBJet1,
				    chi2W, chi2Top,
				    chi2EtMissX, chi2EtMissY,
				    chi2BJetPx, chi2BJetPy, chi2BJetPz, chi2BJetE,
				    chi2LeptonPx, chi2LeptonPy, chi2LeptonPz, chi2LeptonE,
				    par,proc);
	     
	     if( lh < 10E+10 )
	       {
		  vp1.push_back(par[0]);
		  vp2.push_back(par[1]);
		  vp3.push_back(par[2]);
		  vp4.push_back(par[3]);
		  vp5.push_back(par[4]);
		  vp6.push_back(par[5]);
		  vp7.push_back(par[6]);
		  vp8.push_back(par[7]);
		  vp9.push_back(par[8]);
		  vp10.push_back(par[9]);
		  vp11.push_back(par[10]);
		  vp12.push_back(par[11]);
		  vchi2.push_back(lh);
		  
		  if( lh < 1. )
		    {
		       bp = 1;
		       break;
		    }
	       }
	  }
     }
   
   vp1.push_back(*EtMissX);
   vp2.push_back(*EtMissY);
   vp3.push_back(1);
   vp4.push_back(*WMassBW);
   vp5.push_back(*PxBJet1);
   vp6.push_back(*PyBJet1);
   vp7.push_back(*PzBJet1);
   vp8.push_back(*EBJet1);
   vp9.push_back(*PxLepton1);
   vp10.push_back(*PyLepton1);
   vp11.push_back(*PzLepton1);
   vp12.push_back(*ELepton1);
   
   vchi2.push_back(10E+10);
   
   sortPermVector(vchi2,vchi2);
   sortPermVector(vchi2,vp1);
   sortPermVector(vchi2,vp2);
   sortPermVector(vchi2,vp3);
   sortPermVector(vchi2,vp4);
   sortPermVector(vchi2,vp5);
   sortPermVector(vchi2,vp6);
   sortPermVector(vchi2,vp7);
   sortPermVector(vchi2,vp8);
   sortPermVector(vchi2,vp9);
   sortPermVector(vchi2,vp10);
   sortPermVector(vchi2,vp11);
   sortPermVector(vchi2,vp12);
   
   double gridTime = (double)(clock()-tStart)/CLOCKS_PER_SEC;
//   std::cout << "Grid stats = " << vchi2.size() << ", time = " << gridTime << "s" << std::endl;
}
