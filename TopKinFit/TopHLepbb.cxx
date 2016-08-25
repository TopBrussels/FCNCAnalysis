// FCNC single top + Higgs with Higgs decay to bbbar
#include "TopHLepbb.h"

ClassImp(KINFIT::TopHLepbb)

void fcnTopHLepbb(int &npar, double *gin, double &f, double *par, int iflag);

double funcTopHLepbb(double PxBJetLep,double PyBJetLep,double PzBJetLep,double EBJetLep,
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
		     double *par, double &proc);

KINFIT::TopHLepbb::TopHLepbb()
{
}

KINFIT::TopHLepbb::~TopHLepbb()
{
}

void KINFIT::TopHLepbb::TopHLepbbRun()
{
   checkPDF(hPDFTopWMass.get(),"PDFTopWMass");
   checkPDF(hPDFTopMass.get(),"PDFTopMass");
   checkPDF(hPDFHiggsMass.get(),"PDFHiggsMass");
   
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
   NTerm_ = 100;
   NPar_ = 100;
   NNu_ = 2; // effectively
   
   for(int i=0;i<100;i++) FitPar[i] = 10E+10;
   
   PxNu1_ = 10E+10;
   PyNu1_ = 10E+10;   
   PzNu1_ = 10E+10;

   WMassGen1_ = 10E+10;
   
   TopMass1_ = 10E+10;
   TopMass2_ = 10E+10;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}

   if( TopHLepbb_ElectronIdx ) {delete[] TopHLepbb_ElectronIdx; TopHLepbb_ElectronIdx = 0;}
   if( TopHLepbb_MuonIdx ) {delete[] TopHLepbb_MuonIdx; TopHLepbb_MuonIdx = 0;}
   if( TopHLepbb_BJetLepIdx ) {delete[] TopHLepbb_BJetLepIdx; TopHLepbb_BJetLepIdx = 0;}
   if( TopHLepbb_BJet1Idx ) {delete[] TopHLepbb_BJet1Idx; TopHLepbb_BJet1Idx = 0;}
   if( TopHLepbb_BJet2Idx ) {delete[] TopHLepbb_BJet2Idx; TopHLepbb_BJet2Idx = 0;}
   
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
   
   TopHLepbb_ElectronIdx = new int[NPERMEVENT];
   TopHLepbb_MuonIdx = new int[NPERMEVENT];
   TopHLepbb_BJetLepIdx = new int[NPERMEVENT];
   TopHLepbb_BJet1Idx = new int[NPERMEVENT];
   TopHLepbb_BJet2Idx = new int[NPERMEVENT];
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
	TopHLepbb_ElectronIdx[i] = -1;
	TopHLepbb_MuonIdx[i] = -1;
	TopHLepbb_BJetLepIdx[i] = -1;
	TopHLepbb_BJet1Idx[i] = -1;
	TopHLepbb_BJet2Idx[i] = -1;
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

	     for(int ib2=0;ib2<nBJet;ib2++)
	       {
		  if( ib == ib2 ) continue;
		  
		  *KINFIT::kfit::EBJet2 = BJetE[ib2];
		  *KINFIT::kfit::PxBJet2 = BJetPx[ib2];
		  *KINFIT::kfit::PyBJet2 = BJetPy[ib2];
		  *KINFIT::kfit::PzBJet2 = BJetPz[ib2];
		  
		  for(int ib3=ib2+1;ib3<nBJet;ib3++)
		    {
		       if( ib3 == ib ) continue;
		       
		       *KINFIT::kfit::EBJet3 = BJetE[ib3];
		       *KINFIT::kfit::PxBJet3 = BJetPx[ib3];
		       *KINFIT::kfit::PyBJet3 = BJetPy[ib3];
		       *KINFIT::kfit::PzBJet3 = BJetPz[ib3];

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
		       
		       std::vector<double> vp13;
		       std::vector<double> vp14;
		       std::vector<double> vp15;
		       std::vector<double> vp16;
		       std::vector<double> vp17;
		       std::vector<double> vp18;
		       std::vector<double> vp19;
		       std::vector<double> vp20;
		       
		       std::vector<double> vchi2;

		       calcNuGrid(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,vp9,vp10,vp11,vp12,
				  vp13,vp14,vp15,vp16,vp17,vp18,vp19,vp20,
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
			    
			    par[12] = vp13[ig];
			    par[13] = vp14[ig];
			    par[14] = vp15[ig];
			    par[15] = vp16[ig];
			    par[16] = vp17[ig];
			    par[17] = vp18[ig];
			    par[18] = vp19[ig];
			    par[19] = vp20[ig];
			    
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
				 TopMass2_ = *TopMass2;
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
				 
				 FitPar[12] = (*FitParam).at(12);
				 FitPar[13] = (*FitParam).at(13);
				 FitPar[14] = (*FitParam).at(14);
				 FitPar[15] = (*FitParam).at(15);
				 FitPar[16] = (*FitParam).at(16);
				 FitPar[17] = (*FitParam).at(17);
				 FitPar[18] = (*FitParam).at(18);
				 FitPar[19] = (*FitParam).at(19);
				 
				 if( label_l == 0 )
				   {
				      TopHLepbb_ElectronIdx[NPerm_] = idx_l;
				      TopHLepbb_MuonIdx[NPerm_] = -1;
				   }							  
				 else if( label_l == 1 )
				   {
				      TopHLepbb_MuonIdx[NPerm_] = idx_l;
				      TopHLepbb_ElectronIdx[NPerm_] = -1;
				   }
				 
				 TopHLepbb_BJetLepIdx[NPerm_] = ib;
				 
				 TopHLepbb_BJet1Idx[NPerm_] = ib2;
				 TopHLepbb_BJet2Idx[NPerm_] = ib3;
				 
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
				 
				 parPerm_[NPerm_][12] = FitPar[12];
				 parPerm_[NPerm_][13] = FitPar[13];
				 parPerm_[NPerm_][14] = FitPar[14];
				 parPerm_[NPerm_][15] = FitPar[15];
				 parPerm_[NPerm_][16] = FitPar[16];
				 parPerm_[NPerm_][17] = FitPar[17];
				 parPerm_[NPerm_][18] = FitPar[18];
				 parPerm_[NPerm_][19] = FitPar[19];
				 
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
	       }	     
	  }	
     } // end permutations

   // choose the best permutation
   for(int ip=0;ip<NPerm_;ip++)
     {
	calcVar(ip);
     }
   
//   if( NPerm_ > 0 && chiPerm_[0] > 10E+9 )
//     {
//	EmergencyAssign();
//     }   
}

void KINFIT::TopHLepbb::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int il=0;il<nLepton;il++)
     {	
	for(int ib=0;ib<nBJet;ib++)
	  {
	     for(int ib2=0;ib2<nBJet;ib2++)
	       {
		  if( ib == ib2 ) continue;
		  
		  for(int ib3=ib2+1;ib3<nBJet;ib3++)
		    {
		       if( ib3 == ib ) continue;

		       NPERMEVENT++;
		    }
	       }
	  }
     }
}

void KINFIT::TopHLepbb::EmergencyAssign()
{
   float chi2Min = 10E+10;
   
   for(int il=0;il<nLepton;il++)
     {	
	*EtMissX = MetPx;
	*EtMissY = MetPy;
	
	*PxLepton1 = LeptonPx[il];
	*PyLepton1 = LeptonPy[il];
	*PzLepton1 = LeptonPz[il];
	*ELepton1 = LeptonE[il];
	*LabelLepton1 = LeptonLabel[il];
	     
	int label_l = LeptonLabel[il];
	int idx_l = LeptonIdx[il];
		       
	for(int ib=0;ib<nBJet;ib++)
	  {
	     *EBJet1 = BJetE[ib];
	     *PxBJet1 = BJetPx[ib];
	     *PyBJet1 = BJetPy[ib];
	     *PzBJet1 = BJetPz[ib];

	     for(int ib2=0;ib2<nBJet;ib2++)
	       {
		  if( ib == ib2 ) continue;
		  
		  *EBJet2 = BJetE[ib2];
		  *PxBJet2 = BJetPx[ib2];
		  *PyBJet2 = BJetPy[ib2];
		  *PzBJet2 = BJetPz[ib2];
		  
		  for(int ib3=ib2+1;ib3<nBJet;ib3++)
		    {
		       if( ib3 == ib ) continue;
		       
		       *EBJet3 = BJetE[ib3];
		       *PxBJet3 = BJetPx[ib3];
		       *PyBJet3 = BJetPy[ib3];
		       *PzBJet3 = BJetPz[ib3];

		       float ENu = sqrt((*EtMissX)*(*EtMissX)+(*EtMissY)*(*EtMissY));
		       float ETLep = sqrt((*PxLepton1)*(*PxLepton1)+(*PyLepton1)*(*PyLepton1));
		       float ETBJet = sqrt((*PxBJet1)*(*PxBJet1)+(*PyBJet1)*(*PyBJet1));
		       float mTtop = pow(ENu+ETLep+ETBJet,2)-pow(*EtMissX+*PxLepton1+*PxBJet1,2)-pow(*EtMissY+*PyLepton1+*PyBJet1,2);		       
		       float mHiggs = pow(*EBJet2+*EBJet3,2)-pow(*PxBJet2+*PxBJet3,2)-pow(*PyBJet2+*PyBJet3,2)-pow(*PzBJet2+*PzBJet3,2);
		       mTtop = sqrt(mTtop);
		       mHiggs = sqrt(mHiggs);

		       float chi2 = pow(mTtop-*TopMass,2.)/pow(30,2.)+pow(mHiggs-125.,2.)/pow(15,2.);
		       
		       if( chi2 < chi2Min )
			 {
			    chi2Min = chi2;
			    
			    if( label_l == 0 )
			      {
				 TopHLepbb_ElectronIdx[0] = idx_l;
				 TopHLepbb_MuonIdx[0] = -1;
			      }							  
			    else if( label_l == 1 )
			      {
				 TopHLepbb_MuonIdx[0] = idx_l;
				 TopHLepbb_ElectronIdx[0] = -1;
			      }
			    
			    TopHLepbb_BJetLepIdx[0] = ib;
			    
			    TopHLepbb_BJet1Idx[0] = ib2;
			    TopHLepbb_BJet2Idx[0] = ib3;
			 }		       
		    }		  
	       }	     
	  }	
     }
}

double funcTopHLepbb(double PxBJetLep,double PyBJetLep,double PzBJetLep,double EBJetLep,
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
   
   float mHiggs = pow(par[15]+par[19],2)-pow(par[12]+par[16],2)-pow(par[13]+par[17],2)-pow(par[14]+par[18],2);

   chi2WLep = 10E+10;
   chi2TopLep = 10E+10;
   chi2Higgs = 10E+10;
   chi2EtMissX = 10E+10;
   chi2EtMissY = 10E+10;
   chi2BJetPxLep = 10E+10;
   chi2BJetPyLep = 10E+10;
   chi2BJetPzLep = 10E+10;
   chi2BJetELep = 10E+10;
   chi2LeptonPx = 10E+10;
   chi2LeptonPy = 10E+10;
   chi2LeptonPz = 10E+10;
   chi2LeptonE = 10E+10;
   chi2BJet1Px = 10E+10;
   chi2BJet1Py = 10E+10;
   chi2BJet1Pz = 10E+10;
   chi2BJet1E = 10E+10;
   chi2BJet2Px = 10E+10;
   chi2BJet2Py = 10E+10;
   chi2BJet2Pz = 10E+10;
   chi2BJet2E = 10E+10;

   if( mtop >= 0 && rac >= 0 && mHiggs >= 0 )
     {
	mtop = sqrt(mtop);
	mW = sqrt(mW);

	mHiggs = sqrt(mHiggs);
	
	*KINFIT::kfit::TopMass1 = mtop;
	
	float mWProb = getProb(KINFIT::kfit::hPDFTopWMass.get(),mW);
	float mTopProb = getProb(KINFIT::kfit::hPDFTopMass.get(),mtop);

	float mHiggsProb = getProb(KINFIT::kfit::hPDFHiggsMass.get(),mHiggs);
	
	float MetPxProb = getProb(KINFIT::kfit::hPDFMetPx.get(),par[0]-*KINFIT::kfit::EtMissX);
	float MetPyProb = getProb(KINFIT::kfit::hPDFMetPy.get(),par[1]-*KINFIT::kfit::EtMissY);

	float BJetPxLepProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[4]-PxBJetLep);
	float BJetPyLepProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[5]-PyBJetLep);
	float BJetPzLepProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[6]-PzBJetLep);
	float BJetELepProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[7]-EBJetLep);
	
	float LeptonPxProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPx.get(),par[8]-PxLepton) : getProb(KINFIT::kfit::hPDFMuonPx.get(),par[8]-PxLepton);
	float LeptonPyProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPy.get(),par[9]-PyLepton) : getProb(KINFIT::kfit::hPDFMuonPy.get(),par[9]-PyLepton);
	float LeptonPzProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecPz.get(),par[10]-PzLepton) : getProb(KINFIT::kfit::hPDFMuonPz.get(),par[10]-PzLepton);
	float LeptonEProb = (LabelLepton == 0) ? getProb(KINFIT::kfit::hPDFElecE.get(),par[11]-ELepton) : getProb(KINFIT::kfit::hPDFMuonE.get(),par[11]-ELepton);

	float BJet1PxProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[12]-PxBJetHig1);
	float BJet1PyProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[13]-PyBJetHig1);
	float BJet1PzProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[14]-PzBJetHig1);
	float BJet1EProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[15]-EBJetHig1);

	float BJet2PxProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[16]-PxBJetHig2);
	float BJet2PyProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[17]-PyBJetHig2);
	float BJet2PzProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[18]-PzBJetHig2);
	float BJet2EProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[19]-EBJetHig2);
	
	chi2WLep = mWProb;
	chi2TopLep = mTopProb;
	chi2Higgs = mHiggsProb;
	chi2EtMissX = MetPxProb;
	chi2EtMissY = MetPyProb;
	chi2BJetPxLep = BJetPxLepProb;
	chi2BJetPyLep = BJetPyLepProb;
	chi2BJetPzLep = BJetPzLepProb;
	chi2BJetELep = BJetELepProb;
	chi2LeptonPx = LeptonPxProb;
	chi2LeptonPy = LeptonPyProb;
	chi2LeptonPz = LeptonPzProb;
	chi2LeptonE = LeptonEProb;
	chi2BJet1Px = BJet1PxProb;
	chi2BJet1Py = BJet1PyProb;
	chi2BJet1Pz = BJet1PzProb;
	chi2BJet1E = BJet1EProb;
	chi2BJet2Px = BJet2PxProb;
	chi2BJet2Py = BJet2PyProb;
	chi2BJet2Pz = BJet2PzProb;
	chi2BJet2E = BJet2EProb;
	
	chi2WLep = (chi2WLep < 10E-10) ? 10E-10 : chi2WLep;
	chi2TopLep = (chi2TopLep < 10E-10) ? 10E-10 : chi2TopLep;
	chi2Higgs = (chi2Higgs < 10E-10) ? 10E-10 : chi2Higgs;
	
	val = -2*log(
		     chi2WLep*
		     chi2TopLep*
		     chi2Higgs
//		     chi2EtMissX*
//		     chi2EtMissY
		    );
     }

   proc = (double)(clock()-tStart);
   
   return val;
}

void KINFIT::TopHLepbb::calcVar(int iPerm)
{
//   if( chiPerm_[iPerm] > 10E+9 ) return;
   
   int idxElec = TopHLepbb_ElectronIdx[iPerm];
   int idxMuon = TopHLepbb_MuonIdx[iPerm];
   int idxBJetLep = TopHLepbb_BJetLepIdx[iPerm];
   int idxBJet1 = TopHLepbb_BJet1Idx[iPerm];
   int idxBJet2 = TopHLepbb_BJet2Idx[iPerm];
   
   float Lepton1Px = 0;
   float Lepton1Py = 0;
   float Lepton1Pz = 0;
   float Lepton1E = 0;

   float BJetLepPx = 0;
   float BJetLepPy = 0;
   float BJetLepPz = 0;
   float BJetLepE = 0;

   float BJetHig1Px = 0;
   float BJetHig1Py = 0;
   float BJetHig1Pz = 0;
   float BJetHig1E = 0;

   float BJetHig2Px = 0;
   float BJetHig2Py = 0;
   float BJetHig2Pz = 0;
   float BJetHig2E = 0;
   
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

   if( idxBJetLep >= 0 )
     {
	BJetLepPx = BJetPx[idxBJetLep];
	BJetLepPy = BJetPy[idxBJetLep];
	BJetLepPz = BJetPz[idxBJetLep];
	BJetLepE = BJetE[idxBJetLep];
     }   

   if( idxBJet1 >= 0 )
     {
	BJetHig1Px = BJetPx[idxBJet1];
	BJetHig1Py = BJetPy[idxBJet1];
	BJetHig1Pz = BJetPz[idxBJet1];
	BJetHig1E = BJetE[idxBJet1];
     }   

   if( idxBJet2 >= 0 )
     {
	BJetHig2Px = BJetPx[idxBJet2];
	BJetHig2Py = BJetPy[idxBJet2];
	BJetHig2Pz = BJetPz[idxBJet2];
	BJetHig2E = BJetE[idxBJet2];
     }
   
   float nuPx = nuPxPerm_[iPerm][0];
   float nuPy = nuPyPerm_[iPerm][0];
   float nuPz = nuPzPerm_[iPerm][0];   
   float nuE = sqrt(nuPx*nuPx+nuPy*nuPy+nuPz*nuPz);

   // WLep
   float WE = nuE+Lepton1E;
   float WPx = nuPx+Lepton1Px;
   float WPy = nuPy+Lepton1Py;
   float WPz = nuPz+Lepton1Pz;

   float WPt = sqrt(WPx*WPx+WPy*WPy);
   float WEta = getEta(WPt,WPz);
   float WRap = getRap(WE,WPz);
   float WPhi = atan(WPx/WPy);

   // Higgs
   float HiggsE = BJetHig1E+BJetHig2E;
   float HiggsPx = BJetHig1Px+BJetHig2Px;
   float HiggsPy = BJetHig1Py+BJetHig2Py;
   float HiggsPz = BJetHig1Pz+BJetHig2Pz;

   float HiggsPt = sqrt(HiggsPx*HiggsPx+HiggsPy*HiggsPy);
   float HiggsEta = getEta(HiggsPt,HiggsPz);
   float HiggsRap = getRap(HiggsE,HiggsPz);
   float HiggsPhi = atan(HiggsPx/HiggsPy);
   
   // top
   float topE = WE+BJetLepE;
   float topPx = WPx+BJetLepPx;
   float topPy = WPy+BJetLepPy;
   float topPz = WPz+BJetLepPz;
   
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
   
   WMass_[iPerm][1] = sqrt(HiggsE*HiggsE-HiggsPx*HiggsPx-HiggsPy*HiggsPy-HiggsPz*HiggsPz);
}

void fcnTopHLepbb(int &npar, double *gin, double &f, double *par, int iflag)
{
   double chi2WLep;
   double chi2TopLep;
   double chi2Higgs;
   double chi2EtMissX;
   double chi2EtMissY;
   double chi2BJetPxLep;
   double chi2BJetPyLep;
   double chi2BJetPzLep;
   double chi2BJetELep;
   double chi2LeptonPx;
   double chi2LeptonPy;
   double chi2LeptonPz;
   double chi2LeptonE;
   double chi2BJet1Px;
   double chi2BJet1Py;
   double chi2BJet1Pz;
   double chi2BJet1E;
   double chi2BJet2Px;
   double chi2BJet2Py;
   double chi2BJet2Pz;
   double chi2BJet2E;

   double proc;
   double lh = funcTopHLepbb(*KINFIT::kfit::PxBJet1,
			     *KINFIT::kfit::PyBJet1,
			     *KINFIT::kfit::PzBJet1,
			     *KINFIT::kfit::EBJet1,
			     *KINFIT::kfit::PxLepton1,
			     *KINFIT::kfit::PyLepton1,
			     *KINFIT::kfit::PzLepton1,
			     *KINFIT::kfit::ELepton1,
			     *KINFIT::kfit::LabelLepton1,			     
			     *KINFIT::kfit::PxBJet2,
			     *KINFIT::kfit::PyBJet2,
			     *KINFIT::kfit::PzBJet2,
			     *KINFIT::kfit::EBJet2,
			     *KINFIT::kfit::PxBJet3,
			     *KINFIT::kfit::PyBJet3,
			     *KINFIT::kfit::PzBJet3,
			     *KINFIT::kfit::EBJet3,
			     chi2WLep, chi2TopLep,
			     chi2Higgs,
			     chi2EtMissX, chi2EtMissY,
			     chi2BJetPxLep, chi2BJetPyLep, chi2BJetPzLep, chi2BJetELep,
			     chi2LeptonPx, chi2LeptonPy, chi2LeptonPz, chi2LeptonE,
			     chi2BJet1Px, chi2BJet1Py, chi2BJet1Pz, chi2BJet1E,
			     chi2BJet2Px, chi2BJet2Py, chi2BJet2Pz, chi2BJet2E,
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
	
	(*KINFIT::kfit::FitParam).push_back(par[12]);
	(*KINFIT::kfit::FitParam).push_back(par[13]);
	(*KINFIT::kfit::FitParam).push_back(par[14]);
	(*KINFIT::kfit::FitParam).push_back(par[15]);
	(*KINFIT::kfit::FitParam).push_back(par[16]);
	(*KINFIT::kfit::FitParam).push_back(par[17]);
	(*KINFIT::kfit::FitParam).push_back(par[18]);
	(*KINFIT::kfit::FitParam).push_back(par[19]);
	
	(*KINFIT::kfit::ChiTerm).push_back(chi2WLep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2TopLep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Higgs);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissX);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissY);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPxLep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPyLep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetPzLep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJetELep);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPx);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPy);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonPz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2LeptonE);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1E);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2E);
     }
   
   f = lh;
}

void KINFIT::TopHLepbb::fit(double *chis,
			    double *par)
{
   (*KINFIT::kfit::FitParam).clear();
   (*KINFIT::kfit::ChiTerm).clear();
   *KINFIT::kfit::CHISQ = 10E+10;
   
   double perr[100];
   
   TMinuit *gMinuit = new TMinuit(100);
   gMinuit->SetFCN(fcnTopHLepbb);
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
   
   gMinuit->mnparm(4,"BJetPxLep",par[4],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(5,"BJetPyLep",par[5],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(6,"BJetPzLep",par[6],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(7,"BJetELep",par[7],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(8,"LeptonPx",par[8],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(9,"LeptonPy",par[9],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(10,"LeptonPz",par[10],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(11,"LeptonE",par[11],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(12,"BJetHig1Px",par[12],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(13,"BJetHig1Py",par[13],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(14,"BJetHig1Pz",par[14],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(15,"BJetHig1E",par[15],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(16,"BJetHig2Px",par[16],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(17,"BJetHig2Py",par[17],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(18,"BJetHig2Pz",par[18],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(19,"BJetHig2E",par[19],pow(10.,-1),-1000.,1000.,ierflg);

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
   
   gMinuit->FixParameter(12);
   gMinuit->FixParameter(13);
   gMinuit->FixParameter(14);
   gMinuit->FixParameter(15);
   gMinuit->FixParameter(16);
   gMinuit->FixParameter(17);
   gMinuit->FixParameter(18);
   gMinuit->FixParameter(19);

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
   
   gMinuit->GetParameter(12,par[12],perr[12]);
   gMinuit->GetParameter(13,par[13],perr[13]);
   gMinuit->GetParameter(14,par[14],perr[14]);
   gMinuit->GetParameter(15,par[15],perr[15]);
   gMinuit->GetParameter(16,par[16],perr[16]);
   gMinuit->GetParameter(17,par[17],perr[17]);
   gMinuit->GetParameter(18,par[18],perr[18]);
   gMinuit->GetParameter(19,par[19],perr[19]);

   delete gMinuit;
}

void KINFIT::TopHLepbb::calcNuGrid(std::vector<double> &vp1,
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
				   std::vector<double> &vp13,
				   std::vector<double> &vp14,
				   std::vector<double> &vp15,
				   std::vector<double> &vp16,
				   std::vector<double> &vp17,
				   std::vector<double> &vp18,
				   std::vector<double> &vp19,
				   std::vector<double> &vp20,
				   std::vector<double> &vchi2)
{
   double par[NPARMAX];

   clock_t tStart = clock();
   
   bool bp = 0;
   
   double proc = 0.;
   
   for(int it=0;it<NToy_;it++)
     {	
	if( bp ) break;

	do 
	  {	     
	     par[4] = getProbGaus(hPDFBJetPx.get(),maxPDFBJetPx,meanPDFBJetPx,sigmaPDFBJetPx,rnd,NBJetPxRMS_,proc)+*KINFIT::kfit::PxBJet1;
	     par[5] = getProbGaus(hPDFBJetPy.get(),maxPDFBJetPy,meanPDFBJetPy,sigmaPDFBJetPy,rnd,NBJetPyRMS_,proc)+*KINFIT::kfit::PyBJet1;
	     par[6] = getProbGaus(hPDFBJetPz.get(),maxPDFBJetPz,meanPDFBJetPz,sigmaPDFBJetPz,rnd,NBJetPzRMS_,proc)+*KINFIT::kfit::PzBJet1;
	     par[7] = getProbGaus(hPDFBJetE.get(),maxPDFBJetE,meanPDFBJetE,sigmaPDFBJetE,rnd,NBJetERMS_,proc)+*KINFIT::kfit::EBJet1;
	  } while( (par[7]*par[7]-par[4]*par[4]-par[5]*par[5]-par[6]*par[6]) < 0. );
	
	par[0] = getProbGaus(hPDFMetPx.get(),maxPDFMetPx,meanPDFMetPx,sigmaPDFMetPx,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissX;	
	par[1] = getProbGaus(hPDFMetPy.get(),maxPDFMetPy,meanPDFMetPy,sigmaPDFMetPy,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissY;

	par[3] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);

	if( *KINFIT::kfit::LabelLepton1 == 0 )
	  {
	     do 
	       {		  
		  par[8] = getProbGaus(hPDFElecPx.get(),maxPDFElecPx,meanPDFElecPx,sigmaPDFElecPx,rnd,NElecPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[9] = getProbGaus(hPDFElecPy.get(),maxPDFElecPy,meanPDFElecPy,sigmaPDFElecPy,rnd,NElecPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[10] = getProbGaus(hPDFElecPz.get(),maxPDFElecPz,meanPDFElecPz,sigmaPDFElecPz,rnd,NElecPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[11] = getProbGaus(hPDFElecE.get(),maxPDFElecE,meanPDFElecE,sigmaPDFElecE,rnd,NElecERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[11]*par[11]-par[8]*par[8]-par[9]*par[9]-par[10]*par[10]) < 0. );	     
	  }
	else
	  {
	     do
	       {		  
		  par[8] = getProbGaus(hPDFMuonPx.get(),maxPDFMuonPx,meanPDFMuonPx,sigmaPDFMuonPx,rnd,NMuonPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[9] = getProbGaus(hPDFMuonPy.get(),maxPDFMuonPy,meanPDFMuonPy,sigmaPDFMuonPy,rnd,NMuonPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[10] = getProbGaus(hPDFMuonPz.get(),maxPDFMuonPz,meanPDFMuonPz,sigmaPDFMuonPz,rnd,NMuonPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[11] = getProbGaus(hPDFMuonE.get(),maxPDFMuonE,meanPDFMuonE,sigmaPDFMuonE,rnd,NMuonERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[11]*par[11]-par[8]*par[8]-par[9]*par[9]-par[10]*par[10]) < 0. );
	  }

	do
	  {	     
	     par[12] = getProbGaus(hPDFBJetPx.get(),maxPDFBJetPx,meanPDFBJetPx,sigmaPDFBJetPx,rnd,NBJetPxRMS_,proc)+*KINFIT::kfit::PxBJet2;
	     par[13] = getProbGaus(hPDFBJetPy.get(),maxPDFBJetPy,meanPDFBJetPy,sigmaPDFBJetPy,rnd,NBJetPyRMS_,proc)+*KINFIT::kfit::PyBJet2;		  
	     par[14] = getProbGaus(hPDFBJetPz.get(),maxPDFBJetPz,meanPDFBJetPz,sigmaPDFBJetPz,rnd,NBJetPzRMS_,proc)+*KINFIT::kfit::PzBJet2;
	     par[15] = getProbGaus(hPDFBJetE.get(),maxPDFBJetE,meanPDFBJetE,sigmaPDFBJetE,rnd,NBJetERMS_,proc)+*KINFIT::kfit::EBJet2;
	  } while( (par[15]*par[15]-par[12]*par[12]-par[13]*par[13]-par[14]*par[14]) < 0. );

	do
	  {	     
	     par[16] = getProbGaus(hPDFBJetPx.get(),maxPDFBJetPx,meanPDFBJetPx,sigmaPDFBJetPx,rnd,NBJetPxRMS_,proc)+*KINFIT::kfit::PxBJet3;
	     par[17] = getProbGaus(hPDFBJetPy.get(),maxPDFBJetPy,meanPDFBJetPy,sigmaPDFBJetPy,rnd,NBJetPyRMS_,proc)+*KINFIT::kfit::PyBJet3;		  
	     par[18] = getProbGaus(hPDFBJetPz.get(),maxPDFBJetPz,meanPDFBJetPz,sigmaPDFBJetPz,rnd,NBJetPzRMS_,proc)+*KINFIT::kfit::PzBJet3;
	     par[19] = getProbGaus(hPDFBJetE.get(),maxPDFBJetE,meanPDFBJetE,sigmaPDFBJetE,rnd,NBJetERMS_,proc)+*KINFIT::kfit::EBJet3;
	  } while( (par[19]*par[19]-par[16]*par[16]-par[17]*par[17]-par[18]*par[18]) < 0. );	

	for(int is1=-1;is1<=1;is1++)
	  {	
	     if( bp ) break;
	     
	     if( is1 == 0 ) continue;
	     
	     par[2] = is1;

	     double chi2WLep;
	     double chi2TopLep;
	     double chi2Higgs;
	     double chi2EtMissX;
	     double chi2EtMissY;
	     double chi2BJetPxLep;
	     double chi2BJetPyLep;
	     double chi2BJetPzLep;
	     double chi2BJetELep;
	     double chi2LeptonPx;
	     double chi2LeptonPy;
	     double chi2LeptonPz;
	     double chi2LeptonE;
	     double chi2BJet1Px;
	     double chi2BJet1Py;
	     double chi2BJet1Pz;
	     double chi2BJet1E;
	     double chi2BJet2Px;
	     double chi2BJet2Py;
	     double chi2BJet2Pz;
	     double chi2BJet2E;
		  
	     double lh = funcTopHLepbb(*KINFIT::kfit::PxBJet1,
				       *KINFIT::kfit::PyBJet1,
				       *KINFIT::kfit::PzBJet1,
				       *KINFIT::kfit::EBJet1,
				       *KINFIT::kfit::PxLepton1,
				       *KINFIT::kfit::PyLepton1,
				       *KINFIT::kfit::PzLepton1,
				       *KINFIT::kfit::ELepton1,
				       *KINFIT::kfit::LabelLepton1,
				       *KINFIT::kfit::PxBJet2,
				       *KINFIT::kfit::PyBJet2,
				       *KINFIT::kfit::PzBJet2,
				       *KINFIT::kfit::EBJet2,
				       *KINFIT::kfit::PxBJet3,
				       *KINFIT::kfit::PyBJet3,
				       *KINFIT::kfit::PzBJet3,
				       *KINFIT::kfit::EBJet3,
				       chi2WLep, chi2TopLep,
				       chi2Higgs,
				       chi2EtMissX, chi2EtMissY,
				       chi2BJetPxLep, chi2BJetPyLep, chi2BJetPzLep, chi2BJetELep,
				       chi2LeptonPx, chi2LeptonPy, chi2LeptonPz, chi2LeptonE,
				       chi2BJet1Px, chi2BJet1Py, chi2BJet1Pz, chi2BJet1E,
				       chi2BJet2Px, chi2BJet2Py, chi2BJet2Pz, chi2BJet2E,
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
		  
		  vp13.push_back(par[12]);
		  vp14.push_back(par[13]);
		  vp15.push_back(par[14]);
		  vp16.push_back(par[15]);
		  vp17.push_back(par[16]);
		  vp18.push_back(par[17]);
		  vp19.push_back(par[18]);
		  vp20.push_back(par[19]);
		  
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
   vp13.push_back(*PxBJet2);
   vp14.push_back(*PyBJet2);
   vp15.push_back(*PzBJet2);
   vp16.push_back(*EBJet2);
   vp17.push_back(*PxBJet3);
   vp18.push_back(*PyBJet3);
   vp19.push_back(*PzBJet3);
   vp20.push_back(*EBJet3);
   
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
   
   sortPermVector(vchi2,vp13);
   sortPermVector(vchi2,vp14);
   sortPermVector(vchi2,vp15);
   sortPermVector(vchi2,vp16);
   sortPermVector(vchi2,vp17);
   sortPermVector(vchi2,vp18);
   sortPermVector(vchi2,vp19);
   sortPermVector(vchi2,vp20);
   
   double gridTime = (double)(clock()-tStart)/CLOCKS_PER_SEC;
//   std::cout << "Grid stats = " << vchi2.size() << ", time = " << gridTime << "s" << std::endl;
}
