// TTbar dilepton channel: t -> (W->lnu) b, tbar -> (W->lnu) b
#include "TopTopLepLep.h"

ClassImp(KINFIT::TopTopLepLep)

void fcnTopTopLepLep(int &npar, double *gin, double &f, double *par, int iflag);

double funcTopTopLepLep(double PxLepton1,double PyLepton1,double PzLepton1,double ELepton1,int LabelLepton1,
			double PxLepton2,double PyLepton2,double PzLepton2,double ELepton2,int LabelLepton2,
			double PxBJet1,double PyBJet1,double PzBJet1,double EBJet1,
			double PxBJet2,double PyBJet2,double PzBJet2,double EBJet2,
			double &chi2W1, double &chi2W2, double &chi2Top1, double &chi2Top2,
			double &chi2EtMissX, double &chi2EtMissY,
			double &chi2BJet1Px, double &chi2BJet1Py, double &chi2BJet1Pz, double &chi2BJet1E,
			double &chi2BJet2Px, double &chi2BJet2Py, double &chi2BJet2Pz, double &chi2BJet2E,
			double &chi2Lepton1Px, double &chi2Lepton1Py, double &chi2Lepton1Pz, double &chi2Lepton1E,
			double &chi2Lepton2Px, double &chi2Lepton2Py, double &chi2Lepton2Pz, double &chi2Lepton2E,
			double *par, double &proc);

KINFIT::TopTopLepLep::TopTopLepLep()
{
}

KINFIT::TopTopLepLep::~TopTopLepLep()
{
}

void KINFIT::TopTopLepLep::TopTopLepLepRun()
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
   NNu_ = 2;
   
   for(int i=0;i<100;i++) FitPar[i] = 10E+10;
   
   PxNu1_ = 10E+10;
   PxNu2_ = 10E+10;

   PyNu1_ = 10E+10;
   PyNu2_ = 10E+10;
   
   PzNu1_ = 10E+10;
   PzNu2_ = 10E+10;

   WMassGen1_ = 10E+10;
   WMassGen2_ = 10E+10;
   
   TopMass1_ = 10E+10;
   TopMass2_ = 10E+10;

   if( idxMin_ ) {delete[] idxMin_; idxMin_ = 0;}
   
   if( drTopTop_ ) {delete[] drTopTop_; drTopTop_ = 0;}
   if( mTopTop_ ) {delete[] mTopTop_; mTopTop_ = 0;}
   if( ptTopTop_ ) {delete[] ptTopTop_; ptTopTop_ = 0;}
   if( pTopTop_ ) {delete[] pTopTop_; pTopTop_ = 0;}
   if( etaTopTop_ ) {delete[] etaTopTop_; etaTopTop_ = 0;}
   if( phiTopTop_ ) {delete[] phiTopTop_; phiTopTop_ = 0;}
   if( rapTopTop_ ) {delete[] rapTopTop_; rapTopTop_ = 0;}
   
   if( TopTopLepLep_Electron1Idx ) {delete[] TopTopLepLep_Electron1Idx; TopTopLepLep_Electron1Idx = 0;}
   if( TopTopLepLep_Muon1Idx ) {delete[] TopTopLepLep_Muon1Idx; TopTopLepLep_Muon1Idx = 0;}
   if( TopTopLepLep_Electron2Idx ) {delete[] TopTopLepLep_Electron2Idx; TopTopLepLep_Electron2Idx = 0;}
   if( TopTopLepLep_Muon2Idx ) {delete[] TopTopLepLep_Muon2Idx; TopTopLepLep_Muon2Idx = 0;}
   if( TopTopLepLep_BJet1Idx ) {delete[] TopTopLepLep_BJet1Idx; TopTopLepLep_BJet1Idx = 0;}
   if( TopTopLepLep_BJet2Idx ) {delete[] TopTopLepLep_BJet2Idx; TopTopLepLep_BJet2Idx = 0;}
   
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

   drTopTop_ = new float[NPERMEVENT];
   mTopTop_ = new float[NPERMEVENT];
   ptTopTop_ = new float[NPERMEVENT];
   pTopTop_ = new float[NPERMEVENT];
   etaTopTop_ = new float[NPERMEVENT];
   phiTopTop_ = new float[NPERMEVENT];
   rapTopTop_ = new float[NPERMEVENT];
   
   TopTopLepLep_Electron1Idx = new int[NPERMEVENT];
   TopTopLepLep_Muon1Idx = new int[NPERMEVENT];
   TopTopLepLep_Electron2Idx = new int[NPERMEVENT];
   TopTopLepLep_Muon2Idx = new int[NPERMEVENT];
   TopTopLepLep_BJet1Idx = new int[NPERMEVENT];
   TopTopLepLep_BJet2Idx = new int[NPERMEVENT];
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
	drTopTop_[i] = 10E+10;
	mTopTop_[i] = 10E+10;
	ptTopTop_[i] = 10E+10;
	pTopTop_[i] = 10E+10;
	etaTopTop_[i] = 10E+10;
	phiTopTop_[i] = 10E+10;
	rapTopTop_[i] = 10E+10;

	TopTopLepLep_Electron1Idx[i] = -1;
	TopTopLepLep_Muon1Idx[i] = -1;
	TopTopLepLep_Electron2Idx[i] = -1;
	TopTopLepLep_Muon2Idx[i] = -1;
	TopTopLepLep_BJet1Idx[i] = -1;
	TopTopLepLep_BJet2Idx[i] = -1;
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
	for(int il2=il+1;il2<nLepton;il2++)
	  {	
	     *KINFIT::kfit::EtMissX = MetPx;
	     *KINFIT::kfit::EtMissY = MetPy;
	
	     *KINFIT::kfit::PxLepton1 = LeptonPx[il];
	     *KINFIT::kfit::PyLepton1 = LeptonPy[il];
	     *KINFIT::kfit::PzLepton1 = LeptonPz[il];
	     *KINFIT::kfit::ELepton1 = LeptonE[il];
	     *KINFIT::kfit::LabelLepton1 = LeptonLabel[il];
	     
	     int label1_l = LeptonLabel[il];
	     int idx1_l = LeptonIdx[il];

	     *KINFIT::kfit::PxLepton2 = LeptonPx[il2];
	     *KINFIT::kfit::PyLepton2 = LeptonPy[il2];
	     *KINFIT::kfit::PzLepton2 = LeptonPz[il2];
	     *KINFIT::kfit::ELepton2 = LeptonE[il2];
	     *KINFIT::kfit::LabelLepton2 = LeptonLabel[il2];
	     
	     int label2_l = LeptonLabel[il2];
	     int idx2_l = LeptonIdx[il2];
		       
	     for(int ib=0;ib<nBJet;ib++)
	       {
		  *KINFIT::kfit::EBJet1 = BJetE[ib];
		  *KINFIT::kfit::PxBJet1 = BJetPx[ib];
		  *KINFIT::kfit::PyBJet1 = BJetPy[ib];
		  *KINFIT::kfit::PzBJet1 = BJetPz[ib];

		  for(int ib2=0;ib2<nBJet;ib2++)
		    {
		       if( ib == ib2 ) continue;
		       
		       chiPerm_[NPerm_] = 10E+10;
		       
		       *KINFIT::kfit::EBJet2 = BJetE[ib2];
		       *KINFIT::kfit::PxBJet2 = BJetPx[ib2];
		       *KINFIT::kfit::PyBJet2 = BJetPy[ib2];
		       *KINFIT::kfit::PzBJet2 = BJetPz[ib2];
			    
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
		       std::vector<double> vp21;
		       std::vector<double> vp22;
		       std::vector<double> vp23;
		       std::vector<double> vp24;
		       std::vector<double> vchi2;

		       calcNuGrid(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,vp9,vp10,
				  vp11,vp12,vp13,vp14,vp15,vp16,vp17,vp18,vp19,vp20,
				  vp21,vp22,vp23,vp24,
				  vchi2);
			    
		       int nSolutions = vchi2.size();
		       disc_ = 10E+11;
//		       if( nSolutions >= 10 ) nSolutions = 10;
		       
//		       std::cout << "-----" << std::endl;
		       for(int ig=0;ig<nSolutions;ig++)
			 {		
//			    std::cout << vchi2[ig] << std::endl;
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
			    par[20] = vp21[ig];
			    par[21] = vp22[ig];
			    par[22] = vp23[ig];
			    par[23] = vp24[ig];
			    
			    *WMassGen1 = par[6];
			    *WMassGen2 = par[7];
			    
			    fit(chis,par);
			    
			    float res = (*CHISQ < vchi2[ig]) ? *CHISQ : vchi2[ig];
			    
			    if( res < disc_ )
			      {
				 disc_ = res;
				 PxNu1_ = par[4]/2.+par[0];
				 PxNu2_ = par[4]/2.-par[0];
				 PyNu1_ = par[5]/2.+par[1];
				 PyNu2_ = par[5]/2.-par[1];
				 PzNu1_ = *PzNu1;
				 PzNu2_ = *PzNu2;
				 WMassGen1_ = *WMassGen1;
				 WMassGen2_ = *WMassGen2;
				 TopMass1_ = *TopMass1;
				 TopMass2_ = *TopMass2;
				 MetPx_[NPerm_] = par[4];
				 MetPy_[NPerm_] = par[5];
				 
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
				 FitPar[20] = (*FitParam).at(20);
				 FitPar[21] = (*FitParam).at(21);
				 FitPar[22] = (*FitParam).at(22);
				 FitPar[23] = (*FitParam).at(23);
				 
				 TopTopLepLep_BJet1Idx[NPerm_] = ib;
				 TopTopLepLep_BJet2Idx[NPerm_] = ib2;
				 
				 if( label1_l == 0 ) 
				   {
				      TopTopLepLep_Electron1Idx[NPerm_] = idx1_l;
				      TopTopLepLep_Muon1Idx[NPerm_] = -1;
				   }							  
				 else if( label1_l == 1 ) 
				   {
				      TopTopLepLep_Muon1Idx[NPerm_] = idx1_l;
				      TopTopLepLep_Electron1Idx[NPerm_] = -1;
				   }							  
				 
				 if( label2_l == 0 )
				   {
				      TopTopLepLep_Electron2Idx[NPerm_] = idx2_l;
				      TopTopLepLep_Muon2Idx[NPerm_] = -1;
				   }							  
				 else if( label2_l == 1 ) 
				   {
				      TopTopLepLep_Muon2Idx[NPerm_] = idx2_l;
				      TopTopLepLep_Electron2Idx[NPerm_] = -1;
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
				 parPerm_[NPerm_][12] = FitPar[12];
				 parPerm_[NPerm_][13] = FitPar[13];
				 parPerm_[NPerm_][14] = FitPar[14];
				 parPerm_[NPerm_][15] = FitPar[15];
				 parPerm_[NPerm_][16] = FitPar[16];
				 parPerm_[NPerm_][17] = FitPar[17];
				 parPerm_[NPerm_][18] = FitPar[18];
				 parPerm_[NPerm_][19] = FitPar[19];
				 parPerm_[NPerm_][20] = FitPar[20];
				 parPerm_[NPerm_][21] = FitPar[21];
				 parPerm_[NPerm_][22] = FitPar[22];
				 parPerm_[NPerm_][23] = FitPar[23];
				 
				 nuPxPerm_[NPerm_][0] = PxNu1_;
				 nuPyPerm_[NPerm_][0] = PyNu1_;
				 nuPzPerm_[NPerm_][0] = PzNu1_;

				 nuPxPerm_[NPerm_][1] = PxNu2_;
				 nuPyPerm_[NPerm_][1] = PyNu2_;
				 nuPzPerm_[NPerm_][1] = PzNu2_;
				 
				 WMass_[NPerm_][0] = FitPar[6];
				 WMass_[NPerm_][1] = FitPar[7];
				 
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
/*	std::cout << ip << " " << chiTerm_[ip][0] << " " <<
	  chiTerm_[ip][1] << " " <<
	  chiTerm_[ip][2] << " " <<
	  chiTerm_[ip][3] << " " <<
	  chiTerm_[ip][4] << " " <<
	  chiTerm_[ip][5] << " " <<
	  std::endl;*/
     }
}

void KINFIT::TopTopLepLep::CalcNPerm()
{
   NPERMEVENT = 0;

   for(int il=0;il<nLepton;il++)
     {	
	for(int il2=il+1;il2<nLepton;il2++)
	  {	
	     for(int ib=0;ib<nBJet;ib++)
	       {
		  for(int ib2=0;ib2<nBJet;ib2++)
		    {
		       if( ib == ib2 ) continue;
		       
		       NPERMEVENT++;
		    }
	       }
	  }
     }
}

double funcTopTopLepLep(double PxLepton1,double PyLepton1,double PzLepton1,double ELepton1,int LabelLepton1,
			double PxLepton2,double PyLepton2,double PzLepton2,double ELepton2,int LabelLepton2,
			double PxBJet1,double PyBJet1,double PzBJet1,double EBJet1,
			double PxBJet2,double PyBJet2,double PzBJet2,double EBJet2,
			double &chi2W1, double &chi2W2, double &chi2Top1, double &chi2Top2,
			double &chi2EtMissX, double &chi2EtMissY,
			double &chi2BJet1Px, double &chi2BJet1Py, double &chi2BJet1Pz, double &chi2BJet1E,
			double &chi2BJet2Px, double &chi2BJet2Py, double &chi2BJet2Pz, double &chi2BJet2E,
			double &chi2Lepton1Px, double &chi2Lepton1Py, double &chi2Lepton1Pz, double &chi2Lepton1E,
			double &chi2Lepton2Px, double &chi2Lepton2Py, double &chi2Lepton2Pz, double &chi2Lepton2E,
			double *par, double &proc)
{
   clock_t tStart = clock();
   
   double val = 10E+10;

   float PxNu1 = par[4]/2.+par[0];
   float PyNu1 = par[5]/2.+par[1];

   float PxNu2 = par[4]/2.-par[0];
   float PyNu2 = par[5]/2.-par[1];
   
   float a1 = sqrt(par[16]*par[16]+par[17]*par[17]);
   float a2 = sqrt(par[20]*par[20]+par[21]*par[21]);
   float b1 = par[18];
   float b2 = par[22];
   float d1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1);
   float d2 = sqrt(PxNu2*PxNu2+PyNu2*PyNu2);
   float f1 = par[19];
   float f2 = par[23];

   float c1 = par[6]*par[6]/2+par[16]*PxNu1+par[17]*PyNu1;
   float c2 = par[7]*par[7]/2+par[20]*PxNu2+par[21]*PyNu2;
   
   float rac1 = c1*c1*b1*b1-a1*a1*(d1*d1*f1*f1-c1*c1);
   float rac2 = c2*c2*b2*b2-a2*a2*(d2*d2*f2*f2-c2*c2);
   
   float racAbs1 = fabs(rac1);
   float racAbs2 = fabs(rac2);
   
   float PzNu1 = 0.;
   float PzNu2 = 0.;
   
   if( rac1 >= 0 )
     {
	PzNu1 = (c1*b1+par[2]*sqrt(rac1))/a1/a1;
     }   

   if( rac2 >= 0 )
     {
	PzNu2 = (c2*b2+par[3]*sqrt(rac2))/a2/a2;
     }

//   PzNu1 = (c1*b1+par[2]*sqrt(racAbs1))/a1/a1;
//   PzNu2 = (c2*b2+par[3]*sqrt(racAbs2))/a2/a2;
   
   *KINFIT::kfit::PzNu1 = PzNu1;
   *KINFIT::kfit::PzNu2 = PzNu2;
   
   float ENu1 = sqrt(PxNu1*PxNu1+PyNu1*PyNu1+PzNu1*PzNu1);
   float ENu2 = sqrt(PxNu2*PxNu2+PyNu2*PyNu2+PzNu2*PzNu2);
   
   float totPx1 = par[8]+PxNu1+par[16];
   float totPy1 = par[9]+PyNu1+par[17];
   float totPz1 = par[10]+PzNu1+par[18];
   float totE1 = par[11]+ENu1+par[19];

   float totPx2 = par[12]+PxNu2+par[20];
   float totPy2 = par[13]+PyNu2+par[21];
   float totPz2 = par[14]+PzNu2+par[22];
   float totE2 = par[15]+ENu2+par[23];
   
   float mW1 = pow(ENu1+par[19],2)-pow(PxNu1+par[16],2)-pow(PyNu1+par[17],2)-pow(PzNu1+par[18],2);
   float mtop1 = totE1*totE1-totPx1*totPx1-totPy1*totPy1-totPz1*totPz1;

   float mW2 = pow(ENu2+par[23],2)-pow(PxNu2+par[20],2)-pow(PyNu2+par[21],2)-pow(PzNu2+par[22],2);
   float mtop2 = totE2*totE2-totPx2*totPx2-totPy2*totPy2-totPz2*totPz2;

   chi2W1 = 10E+10;
   chi2W2 = 10E+10;   
   chi2Top1 = 10E+10;
   chi2Top2 = 10E+10;
   chi2EtMissX = 10E+10;
   chi2EtMissY = 10E+10;
   chi2BJet1Px = 10E+10;
   chi2BJet1Py = 10E+10;
   chi2BJet1Pz = 10E+10;
   chi2BJet1E = 10E+10;
   chi2BJet2Px = 10E+10;
   chi2BJet2Py = 10E+10;
   chi2BJet2Pz = 10E+10;
   chi2BJet2E = 10E+10;
   chi2Lepton1Px = 10E+10;
   chi2Lepton1Py = 10E+10;
   chi2Lepton1Pz = 10E+10;
   chi2Lepton1E = 10E+10;
   chi2Lepton2Px = 10E+10;
   chi2Lepton2Py = 10E+10;
   chi2Lepton2Pz = 10E+10;
   chi2Lepton2E = 10E+10;
   
   if( mtop1 >= 0 && mtop2 >= 0 && rac1 >= 0 && rac2 >= 0 && mW1 >= 0 && mW2 >= 0 )
//   if( mtop1 >= 0 && mtop2 >= 0 )
     {
	mtop1 = sqrt(mtop1);
	mW1 = sqrt(mW1);

	mtop2 = sqrt(mtop2);
	mW2 = sqrt(mW2);

	*KINFIT::kfit::TopMass1 = mtop1;
	*KINFIT::kfit::TopMass2 = mtop2;
	
	float mW1Prob = getProb(KINFIT::kfit::hPDFTopWMass.get(),mW1);
	float mW2Prob = getProb(KINFIT::kfit::hPDFTopWMass.get(),mW2);

	float mTop1Prob = getProb(KINFIT::kfit::hPDFTopMass.get(),mtop1);
	float mTop2Prob = getProb(KINFIT::kfit::hPDFTopMass.get(),mtop2);
	
	float MetPxProb = getProb(KINFIT::kfit::hPDFMetPx.get(),par[4]-*KINFIT::kfit::EtMissX);
	float MetPyProb = getProb(KINFIT::kfit::hPDFMetPy.get(),par[5]-*KINFIT::kfit::EtMissY);

	float BJet1PxProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[8]-PxBJet1);
	float BJet1PyProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[9]-PyBJet1);
	float BJet1PzProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[10]-PzBJet1);
	float BJet1EProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[11]-EBJet1);

	float BJet2PxProb = getProb(KINFIT::kfit::hPDFBJetPx.get(),par[12]-PxBJet2);
	float BJet2PyProb = getProb(KINFIT::kfit::hPDFBJetPy.get(),par[13]-PyBJet2);
	float BJet2PzProb = getProb(KINFIT::kfit::hPDFBJetPz.get(),par[14]-PzBJet2);
	float BJet2EProb = getProb(KINFIT::kfit::hPDFBJetE.get(),par[15]-EBJet2);

	float Lepton1PxProb = (LabelLepton1 == 0) ? getProb(KINFIT::kfit::hPDFElecPx.get(),par[16]-PxLepton1) : getProb(KINFIT::kfit::hPDFMuonPx.get(),par[16]-PxLepton1);
	float Lepton1PyProb = (LabelLepton1 == 0) ? getProb(KINFIT::kfit::hPDFElecPy.get(),par[17]-PyLepton1) : getProb(KINFIT::kfit::hPDFMuonPy.get(),par[17]-PyLepton1);
	float Lepton1PzProb = (LabelLepton1 == 0) ? getProb(KINFIT::kfit::hPDFElecPz.get(),par[18]-PzLepton1) : getProb(KINFIT::kfit::hPDFMuonPz.get(),par[18]-PzLepton1);
	float Lepton1EProb = (LabelLepton1 == 0) ? getProb(KINFIT::kfit::hPDFElecE.get(),par[19]-ELepton1) : getProb(KINFIT::kfit::hPDFMuonE.get(),par[19]-ELepton1);

	float Lepton2PxProb = (LabelLepton2 == 0) ? getProb(KINFIT::kfit::hPDFElecPx.get(),par[20]-PxLepton2) : getProb(KINFIT::kfit::hPDFMuonPx.get(),par[20]-PxLepton2);
	float Lepton2PyProb = (LabelLepton2 == 0) ? getProb(KINFIT::kfit::hPDFElecPy.get(),par[21]-PyLepton2) : getProb(KINFIT::kfit::hPDFMuonPy.get(),par[21]-PyLepton2);
	float Lepton2PzProb = (LabelLepton2 == 0) ? getProb(KINFIT::kfit::hPDFElecPz.get(),par[22]-PzLepton2) : getProb(KINFIT::kfit::hPDFMuonPz.get(),par[22]-PzLepton2);
	float Lepton2EProb = (LabelLepton2 == 0) ? getProb(KINFIT::kfit::hPDFElecE.get(),par[23]-ELepton2) : getProb(KINFIT::kfit::hPDFMuonE.get(),par[23]-ELepton2);

	chi2W1 = mW1Prob;
	chi2W2 = mW2Prob;
	chi2Top1 = mTop1Prob;
	chi2Top2 = mTop2Prob;
	chi2EtMissX = MetPxProb;
	chi2EtMissY = MetPyProb;
	chi2BJet1Px = BJet1PxProb;
	chi2BJet1Py = BJet1PyProb;
	chi2BJet1Pz = BJet1PzProb;
	chi2BJet1E = BJet1EProb;
	chi2BJet2Px = BJet2PxProb;
	chi2BJet2Py = BJet2PyProb;
	chi2BJet2Pz = BJet2PzProb;
	chi2BJet2E = BJet2EProb;
	chi2Lepton1Px = Lepton1PxProb;
	chi2Lepton1Py = Lepton1PyProb;
	chi2Lepton1Pz = Lepton1PzProb;
	chi2Lepton1E = Lepton1EProb;
	chi2Lepton2Px = Lepton2PxProb;
	chi2Lepton2Py = Lepton2PyProb;
	chi2Lepton2Pz = Lepton2PzProb;
	chi2Lepton2E = Lepton2EProb;

	chi2W1 = (chi2W1 < 10E-10) ? 10E-10 : chi2W1;
	chi2W2 = (chi2W2 < 10E-10) ? 10E-10 : chi2W2;
	chi2Top1 = (chi2Top1 < 10E-10) ? 10E-10 : chi2Top1;
	chi2Top2 = (chi2Top2 < 10E-10) ? 10E-10 : chi2Top2;
	chi2EtMissX = (chi2EtMissX < 10E-10) ? 10E-10 : chi2EtMissX;
	chi2EtMissY = (chi2EtMissY < 10E-10) ? 10E-10 : chi2EtMissY;
	
	val = -2*log(
		     chi2W1*
		     chi2W2*
		     chi2Top1*
		     chi2Top2*
		     chi2EtMissX*
		     chi2EtMissY
		     
/*		     chi2BJet1Px*
		     chi2BJet1Py*
		     chi2BJet1Pz*
		     chi2BJet1E*
		     chi2BJet2Px*
		     chi2BJet2Py*
		     chi2BJet2Pz*
		     chi2BJet2E*

		     chi2Lepton1Px*
		     chi2Lepton1Py*
		     chi2Lepton1Pz*
		     chi2Lepton1E*
		     chi2Lepton2Px*
		     chi2Lepton2Py*
		     chi2Lepton2Pz*
		     chi2Lepton2E*/
		    );

//	if( rac1 < 0 ) val += racAbs1;
//	if( rac2 < 0 ) val += racAbs2;
     }

   proc = (double)(clock()-tStart);
   
   return val;
}

void KINFIT::TopTopLepLep::calcVar(int iPerm)
{
//   if( chiPerm_[iPerm] > 10E+9 ) return;
   
   int idxElec1 = TopTopLepLep_Electron1Idx[iPerm];
   int idxMuon1 = TopTopLepLep_Muon1Idx[iPerm];

   int idxElec2 = TopTopLepLep_Electron2Idx[iPerm];
   int idxMuon2 = TopTopLepLep_Muon2Idx[iPerm];
   
   int idxBJet1 = TopTopLepLep_BJet1Idx[iPerm];
   int idxBJet2 = TopTopLepLep_BJet2Idx[iPerm];   
   
   float Lepton1Px = 0;
   float Lepton1Py = 0;
   float Lepton1Pz = 0;
   float Lepton1E = 0;

   float Lepton2Px = 0;
   float Lepton2Py = 0;
   float Lepton2Pz = 0;
   float Lepton2E = 0;

   float BJet1Px = 0;
   float BJet1Py = 0;
   float BJet1Pz = 0;
   float BJet1E = 0;

   float BJet2Px = 0;
   float BJet2Py = 0;
   float BJet2Pz = 0;
   float BJet2E = 0;
   
   if( idxElec1 >= 0 )
     {
	Lepton1Px = ElectronPx[idxElec1];
	Lepton1Py = ElectronPy[idxElec1];
	Lepton1Pz = ElectronPz[idxElec1];
	Lepton1E = ElectronE[idxElec1];
     }   
   else if( idxMuon1 >= 0 )
     {
	Lepton1Px = MuonPx[idxMuon1];
	Lepton1Py = MuonPy[idxMuon1];
	Lepton1Pz = MuonPz[idxMuon1];
	Lepton1E = MuonE[idxMuon1];
     }   

   if( idxElec2 >= 0 )
     {
	Lepton2Px = ElectronPx[idxElec2];
	Lepton2Py = ElectronPy[idxElec2];
	Lepton2Pz = ElectronPz[idxElec2];
	Lepton2E = ElectronE[idxElec2];
     }   
   else if( idxMuon2 >= 0 )
     {
	Lepton2Px = MuonPx[idxMuon2];
	Lepton2Py = MuonPy[idxMuon2];
	Lepton2Pz = MuonPz[idxMuon2];
	Lepton2E = MuonE[idxMuon2];
     }   

   if( idxBJet1 >= 0 )
     {
	BJet1Px = BJetPx[idxBJet1];
	BJet1Py = BJetPy[idxBJet1];
	BJet1Pz = BJetPz[idxBJet1];
	BJet1E = BJetE[idxBJet1];
     }   
   if( idxBJet2 >= 0 )
     {
	BJet2Px = BJetPx[idxBJet2];
	BJet2Py = BJetPy[idxBJet2];
	BJet2Pz = BJetPz[idxBJet2];
	BJet2E = BJetE[idxBJet2];
     }   

   float nu1Px = nuPxPerm_[iPerm][0];
   float nu1Py = nuPyPerm_[iPerm][0];
   float nu1Pz = nuPzPerm_[iPerm][0];

   float nu2Px = nuPxPerm_[iPerm][1];
   float nu2Py = nuPyPerm_[iPerm][1];
   float nu2Pz = nuPzPerm_[iPerm][1];
   
   float nu1E = sqrt(nu1Px*nu1Px+nu1Py*nu1Py+nu1Pz*nu1Pz);
   float nu2E = sqrt(nu2Px*nu2Px+nu2Py*nu2Py+nu2Pz*nu2Pz);

   // W1
   float W1E = nu1E+Lepton1E;
   float W1Px = nu1Px+Lepton1Px;
   float W1Py = nu1Py+Lepton1Py;
   float W1Pz = nu1Pz+Lepton1Pz;

   float W1Pt = sqrt(W1Px*W1Px+W1Py*W1Py);
   float W1Eta = getEta(W1Pt,W1Pz);
   float W1Rap = getRap(W1E,W1Pz);
   float W1Phi = atan(W1Px/W1Py);
   
   // top1
   float top1E = W1E+BJet1E;
   float top1Px = W1Px+BJet1Px;
   float top1Py = W1Py+BJet1Py;
   float top1Pz = W1Pz+BJet1Pz;
   
   float top1Pt = sqrt(top1Px*top1Px+top1Py*top1Py);
   float top1Eta = getEta(top1Pt,top1Pz);
   float top1Rap = getRap(top1E,top1Pz);
   float top1Phi = atan(top1Px/top1Py);

   // W2
   float W2E = nu2E+Lepton2E;
   float W2Px = nu2Px+Lepton2Px;
   float W2Py = nu2Py+Lepton2Py;
   float W2Pz = nu2Pz+Lepton2Pz;

   float W2Pt = sqrt(W2Px*W2Px+W2Py*W2Py);
   float W2Eta = getEta(W2Pt,W2Pz);
   float W2Rap = getRap(W2E,W2Pz);
   float W2Phi = atan(W2Px/W2Py);
   
   // top2
   float top2E = W2E+BJet2E;
   float top2Px = W2Px+BJet2Px;
   float top2Py = W2Py+BJet2Py;
   float top2Pz = W2Pz+BJet2Pz;
   
   float top2Pt = sqrt(top2Px*top2Px+top2Py*top2Py);
   float top2Eta = getEta(top2Pt,top2Pz);
   float top2Rap = getRap(top2E,top2Pz);
   float top2Phi = atan(top2Px/top2Py);
   
   // top1+top2
   float toptopE = top1E+top2E;
   float toptopPx = top1Px+top2Px;
   float toptopPy = top1Py+top2Py;
   float toptopPz = top1Pz+top2Pz;
   
   TLorentzVector *top1 = new TLorentzVector();
   TLorentzVector *top2 = new TLorentzVector();
   
   top1->SetPxPyPzE(top1Px,top1Py,top1Pz,top1E);
   top2->SetPxPyPzE(top2Px,top2Py,top2Pz,top2E);
   
//   drTopTop_[iPerm] = getDeltaR(top1Eta,top1Phi,top2Eta,top2Phi);
   drTopTop_[iPerm] = top1->DeltaR(*top2);
   delete top1;
   delete top2;
   mTopTop_[iPerm] = sqrt(toptopE*toptopE-toptopPx*toptopPx-toptopPy*toptopPy-toptopPz*toptopPz);
   ptTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy);
   pTopTop_[iPerm] = sqrt(toptopPx*toptopPx+toptopPy*toptopPy+toptopPz*toptopPz);
   etaTopTop_[iPerm] = getEta(ptTopTop_[iPerm],toptopPz);
   rapTopTop_[iPerm] = getRap(toptopE,toptopPz);
   phiTopTop_[iPerm] = atan(toptopPx/toptopPy);
   
   TopMass_[iPerm][0] = sqrt(top1E*top1E-top1Px*top1Px-top1Py*top1Py-top1Pz*top1Pz);
   TopMass_[iPerm][1] = sqrt(top2E*top2E-top2Px*top2Px-top2Py*top2Py-top2Pz*top2Pz);
   TopPt_[iPerm][0] = sqrt(top1Px*top1Px+top1Py*top1Py);
   TopPt_[iPerm][1] = sqrt(top2Px*top2Px+top2Py*top2Py);
   TopP_[iPerm][0] = sqrt(top1Px*top1Px+top1Py*top1Py+top1Pz*top1Pz);
   TopP_[iPerm][1] = sqrt(top2Px*top2Px+top2Py*top2Py+top2Pz*top2Pz);
   TopEta_[iPerm][0] = getEta(TopPt_[iPerm][0],top1Pz);
   TopEta_[iPerm][1] = getEta(TopPt_[iPerm][1],top2Pz);
   TopRap_[iPerm][0] = getRap(top1E,top1Pz);
   TopRap_[iPerm][1] = getRap(top2E,top2Pz);

   WPt_[iPerm][0] = sqrt(W1Px*W1Px+W1Py*W1Py);
   WPt_[iPerm][1] = sqrt(W2Px*W2Px+W2Py*W2Py);
   WP_[iPerm][0] = sqrt(W1Px*W1Px+W1Py*W1Py+W1Pz*W1Pz);
   WP_[iPerm][1] = sqrt(W2Px*W2Px+W2Py*W2Py+W2Pz*W2Pz);
   WEta_[iPerm][0] = getEta(WPt_[iPerm][0],W1Pz);
   WEta_[iPerm][1] = getEta(WPt_[iPerm][1],W2Pz);
   WRap_[iPerm][0] = getRap(W1E,W1Pz);
   WRap_[iPerm][1] = getRap(W2E,W2Pz);
   WPhi_[iPerm][0] = atan(W1Px/W1Py);
   WPhi_[iPerm][1] = atan(W2Px/W2Py);
}

void fcnTopTopLepLep(int &npar, double *gin, double &f, double *par, int iflag)
{
   double chi2W1;
   double chi2W2;
   double chi2Top1;
   double chi2Top2;
   double chi2EtMissX;
   double chi2EtMissY;
   double chi2BJet1Px;
   double chi2BJet1Py;
   double chi2BJet1Pz;
   double chi2BJet1E;
   double chi2BJet2Px;
   double chi2BJet2Py;
   double chi2BJet2Pz;
   double chi2BJet2E;
   double chi2Lepton1Px;
   double chi2Lepton1Py;
   double chi2Lepton1Pz;
   double chi2Lepton1E;
   double chi2Lepton2Px;
   double chi2Lepton2Py;
   double chi2Lepton2Pz;
   double chi2Lepton2E;

   double proc;
   double lh = funcTopTopLepLep(*KINFIT::kfit::PxLepton1,
				*KINFIT::kfit::PyLepton1,
				*KINFIT::kfit::PzLepton1,
				*KINFIT::kfit::ELepton1,
				*KINFIT::kfit::LabelLepton1,
				*KINFIT::kfit::PxLepton2,
				*KINFIT::kfit::PyLepton2,
				*KINFIT::kfit::PzLepton2,
				*KINFIT::kfit::ELepton2,
				*KINFIT::kfit::LabelLepton2,
				*KINFIT::kfit::PxBJet1,
				*KINFIT::kfit::PyBJet1,
				*KINFIT::kfit::PzBJet1,
				*KINFIT::kfit::EBJet1,
				*KINFIT::kfit::PxBJet2,
				*KINFIT::kfit::PyBJet2,
				*KINFIT::kfit::PzBJet2,
				*KINFIT::kfit::EBJet2,
				chi2W1, chi2W2, chi2Top1, chi2Top2, chi2EtMissX, chi2EtMissY,
				chi2BJet1Px, chi2BJet1Py, chi2BJet1Pz, chi2BJet1E,
				chi2BJet2Px, chi2BJet2Py, chi2BJet2Pz, chi2BJet2E,
				chi2Lepton1Px, chi2Lepton1Py, chi2Lepton1Pz, chi2Lepton1E,
				chi2Lepton2Px, chi2Lepton2Py, chi2Lepton2Pz, chi2Lepton2E,
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
	(*KINFIT::kfit::FitParam).push_back(par[20]);
	(*KINFIT::kfit::FitParam).push_back(par[21]);
	(*KINFIT::kfit::FitParam).push_back(par[22]);
	(*KINFIT::kfit::FitParam).push_back(par[23]);
	
	(*KINFIT::kfit::ChiTerm).push_back(chi2W1);
	(*KINFIT::kfit::ChiTerm).push_back(chi2W2);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Top1);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Top2);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissX);
	(*KINFIT::kfit::ChiTerm).push_back(chi2EtMissY);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet1E);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2BJet2E);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton1Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton1Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton1Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton1E);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton2Px);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton2Py);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton2Pz);
	(*KINFIT::kfit::ChiTerm).push_back(chi2Lepton2E);
     }
   
   f = lh;
}

void KINFIT::TopTopLepLep::fit(double *chis,
			       double *par)
{
   (*KINFIT::kfit::FitParam).clear();
   (*KINFIT::kfit::ChiTerm).clear();
   *KINFIT::kfit::CHISQ = 10E+10;
   
   double perr[100];
   
   TMinuit *gMinuit = new TMinuit(100);
   gMinuit->SetFCN(fcnTopTopLepLep);
   gMinuit->SetPrintLevel(-1);
   
   Double_t arglist[100];
   Int_t ierflg = 0;
   
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
   
   arglist[0] = 5000;
   arglist[1] = 0.1;
   
   gMinuit->mnparm(0,"Etx",par[0],pow(10.,-2),-1000.,1000.,ierflg);
   gMinuit->mnparm(1,"Ety",par[1],pow(10.,-2),-1000.,1000.,ierflg);
   gMinuit->mnparm(2,"Sign1",par[2],pow(10.,-2),-1.,1.,ierflg);
   gMinuit->mnparm(3,"Sign2",par[3],pow(10.,-2),-1.,1.,ierflg);
   gMinuit->mnparm(4,"EtRealX",par[4],pow(10.,-2),-*KINFIT::kfit::EtMissX*10.,*KINFIT::kfit::EtMissX*10.,ierflg);
   gMinuit->mnparm(5,"EtRealY",par[5],pow(10.,-2),-*KINFIT::kfit::EtMissY*10.,*KINFIT::kfit::EtMissY*10.,ierflg);
   gMinuit->mnparm(6,"mW1",par[6],pow(10.,-2),1.,1000.,ierflg);
   gMinuit->mnparm(7,"mW2",par[7],pow(10.,-2),1.,1000.,ierflg);
   
   gMinuit->mnparm(8,"BJet1Px",par[8],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(9,"BJet1Py",par[9],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(10,"BJet1Pz",par[10],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(11,"BJet1E",par[11],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(12,"BJet2Px",par[12],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(13,"BJet2Py",par[13],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(14,"BJet2Pz",par[14],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(15,"BJet2E",par[15],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(16,"Lepton1Px",par[16],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(17,"Lepton1Py",par[17],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(18,"Lepton1Pz",par[18],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(19,"Lepton1E",par[19],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(20,"Lepton2Px",par[20],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(21,"Lepton2Py",par[21],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(22,"Lepton2Pz",par[22],pow(10.,-1),-1000.,1000.,ierflg);
   gMinuit->mnparm(23,"Lepton2E",par[23],pow(10.,-1),-1000.,1000.,ierflg);

   gMinuit->FixParameter(2);
   gMinuit->FixParameter(3);
   
//   gMinuit->FixParameter(4);
//   gMinuit->FixParameter(5);
//   gMinuit->FixParameter(6);
//   gMinuit->FixParameter(7);
   
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
   gMinuit->FixParameter(20);
   gMinuit->FixParameter(21);
   gMinuit->FixParameter(22);
   gMinuit->FixParameter(23);

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
   gMinuit->GetParameter(20,par[20],perr[20]);
   gMinuit->GetParameter(21,par[21],perr[21]);
   gMinuit->GetParameter(22,par[22],perr[22]);
   gMinuit->GetParameter(23,par[23],perr[23]);

   delete gMinuit;
}

void KINFIT::TopTopLepLep::calcNuGrid(std::vector<double> &vp1,
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
				      std::vector<double> &vp21,
				      std::vector<double> &vp22,
				      std::vector<double> &vp23,
				      std::vector<double> &vp24,
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

	do
	  {	     
	     par[12] = getProbGaus(hPDFBJetPx.get(),maxPDFBJetPx,meanPDFBJetPx,sigmaPDFBJetPx,rnd,NBJetPxRMS_,proc)+*KINFIT::kfit::PxBJet2;
	     par[13] = getProbGaus(hPDFBJetPy.get(),maxPDFBJetPy,meanPDFBJetPy,sigmaPDFBJetPy,rnd,NBJetPyRMS_,proc)+*KINFIT::kfit::PyBJet2;
	     par[14] = getProbGaus(hPDFBJetPz.get(),maxPDFBJetPz,meanPDFBJetPz,sigmaPDFBJetPz,rnd,NBJetPzRMS_,proc)+*KINFIT::kfit::PzBJet2;
	     par[15] = getProbGaus(hPDFBJetE.get(),maxPDFBJetE,meanPDFBJetE,sigmaPDFBJetE,rnd,NBJetERMS_,proc)+*KINFIT::kfit::EBJet2;
	  } while( (par[15]*par[15]-par[12]*par[12]-par[13]*par[13]-par[14]*par[14]) < 0. );	
	
	par[4] = getProbGaus(hPDFMetPx.get(),maxPDFMetPx,meanPDFMetPx,sigmaPDFMetPx,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissX;
	par[5] = getProbGaus(hPDFMetPy.get(),maxPDFMetPy,meanPDFMetPy,sigmaPDFMetPy,rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissY;

	par[6] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);	
	par[7] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);

	if( *KINFIT::kfit::LabelLepton1 == 0 )
	  {
	     do
	       {		  
		  par[16] = getProbGaus(hPDFElecPx.get(),maxPDFElecPx,meanPDFElecPx,sigmaPDFElecPx,rnd,NElecPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[17] = getProbGaus(hPDFElecPy.get(),maxPDFElecPy,meanPDFElecPy,sigmaPDFElecPy,rnd,NElecPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[18] = getProbGaus(hPDFElecPz.get(),maxPDFElecPz,meanPDFElecPz,sigmaPDFElecPz,rnd,NElecPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[19] = getProbGaus(hPDFElecE.get(),maxPDFElecE,meanPDFElecE,sigmaPDFElecE,rnd,NElecERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[19]*par[19]-par[16]*par[16]-par[17]*par[17]-par[18]*par[18]) < 0. );	     
	  }
	else
	  {
	     do
	       {		  
		  par[16] = getProbGaus(hPDFMuonPx.get(),maxPDFMuonPx,meanPDFMuonPx,sigmaPDFMuonPx,rnd,NMuonPxRMS_,proc)+*KINFIT::kfit::PxLepton1;
		  par[17] = getProbGaus(hPDFMuonPy.get(),maxPDFMuonPy,meanPDFMuonPy,sigmaPDFMuonPy,rnd,NMuonPyRMS_,proc)+*KINFIT::kfit::PyLepton1;
		  par[18] = getProbGaus(hPDFMuonPz.get(),maxPDFMuonPz,meanPDFMuonPz,sigmaPDFMuonPz,rnd,NMuonPzRMS_,proc)+*KINFIT::kfit::PzLepton1;
		  par[19] = getProbGaus(hPDFMuonE.get(),maxPDFMuonE,meanPDFMuonE,sigmaPDFMuonE,rnd,NMuonERMS_,proc)+*KINFIT::kfit::ELepton1;
	       } while( (par[19]*par[19]-par[16]*par[16]-par[17]*par[17]-par[18]*par[18]) < 0. );	     
	  }

	if( *KINFIT::kfit::LabelLepton2 == 0 )
	  {
	     do
	       {		  
		  par[20] = getProbGaus(hPDFElecPx.get(),maxPDFElecPx,meanPDFElecPx,sigmaPDFElecPx,rnd,NElecPxRMS_,proc)+*KINFIT::kfit::PxLepton2;
		  par[21] = getProbGaus(hPDFElecPy.get(),maxPDFElecPy,meanPDFElecPy,sigmaPDFElecPy,rnd,NElecPyRMS_,proc)+*KINFIT::kfit::PyLepton2;
		  par[22] = getProbGaus(hPDFElecPz.get(),maxPDFElecPz,meanPDFElecPz,sigmaPDFElecPz,rnd,NElecPzRMS_,proc)+*KINFIT::kfit::PzLepton2;
		  par[23] = getProbGaus(hPDFElecE.get(),maxPDFElecE,meanPDFElecE,sigmaPDFElecE,rnd,NElecERMS_,proc)+*KINFIT::kfit::ELepton2;
	       } while( (par[23]*par[23]-par[20]*par[20]-par[21]*par[21]-par[22]*par[22]) < 0. );	     
	  }
	else
	  {
	     do
	       {		  
		  par[20] = getProbGaus(hPDFMuonPx.get(),maxPDFMuonPx,meanPDFMuonPx,sigmaPDFMuonPx,rnd,NMuonPxRMS_,proc)+*KINFIT::kfit::PxLepton2;
		  par[21] = getProbGaus(hPDFMuonPy.get(),maxPDFMuonPy,meanPDFMuonPy,sigmaPDFMuonPy,rnd,NMuonPyRMS_,proc)+*KINFIT::kfit::PyLepton2;
		  par[22] = getProbGaus(hPDFMuonPz.get(),maxPDFMuonPz,meanPDFMuonPz,sigmaPDFMuonPz,rnd,NMuonPzRMS_,proc)+*KINFIT::kfit::PzLepton2;
		  par[23] = getProbGaus(hPDFMuonE.get(),maxPDFMuonE,meanPDFMuonE,sigmaPDFMuonE,rnd,NMuonERMS_,proc)+*KINFIT::kfit::ELepton2;
	       } while( (par[23]*par[23]-par[20]*par[20]-par[21]*par[21]-par[22]*par[22]) < 0. );	     
	  }

	par[0] = -200.+400.*(rnd->Rndm());
	par[1] = -200.+400.*(rnd->Rndm());	
/*
	par[8] = *KINFIT::kfit::PxBJet1;
	par[9] = *KINFIT::kfit::PyBJet1;
	par[10] = *KINFIT::kfit::PzBJet1;
	par[11] = *KINFIT::kfit::EBJet1;
	par[12] = *KINFIT::kfit::PxBJet2;
	par[13] = *KINFIT::kfit::PyBJet2;
	par[14] = *KINFIT::kfit::PzBJet2;
	par[15] = *KINFIT::kfit::EBJet2;
	par[4] = getProbGaus(hPDFMetPx.get(),rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissX;
	par[5] = getProbGaus(hPDFMetPy.get(),rnd,NMetRMS_,proc)+*KINFIT::kfit::EtMissY;
	par[6] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);
	par[7] = getWmassBW(rnd,*WMassBW,*WRMSBW,NWRMS_,proc);
	par[16] = *KINFIT::kfit::PxLepton1;
	par[17] = *KINFIT::kfit::PyLepton1;
	par[18] = *KINFIT::kfit::PzLepton1;
	par[19] = *KINFIT::kfit::ELepton1;
	par[20] = *KINFIT::kfit::PxLepton2;
	par[21] = *KINFIT::kfit::PyLepton2;
	par[22] = *KINFIT::kfit::PzLepton2;
	par[23] = *KINFIT::kfit::ELepton2;
*/
	for(int is1=-1;is1<=1;is1++)
	  {	
	     if( bp ) break;
	     
	     if( is1 == 0 ) continue;
	     
	     par[2] = is1;
	     
	     for(int is2=-1;is2<=1;is2++)
	       {
		  if( bp ) break;
		  
		  if( is2 == 0 ) continue;
		  
		  par[3] = is2;
		  
		  double chi2W1;
		  double chi2W2;
		  double chi2Top1;
		  double chi2Top2;
		  double chi2EtMissX;
		  double chi2EtMissY;
		  double chi2BJet1Px;
		  double chi2BJet1Py;
		  double chi2BJet1Pz;
		  double chi2BJet1E;
		  double chi2BJet2Px;
		  double chi2BJet2Py;
		  double chi2BJet2Pz;
		  double chi2BJet2E;
		  double chi2Lepton1Px;
		  double chi2Lepton1Py;
		  double chi2Lepton1Pz;
		  double chi2Lepton1E;
		  double chi2Lepton2Px;
		  double chi2Lepton2Py;
		  double chi2Lepton2Pz;
		  double chi2Lepton2E;
		  
		  double lh = funcTopTopLepLep(*KINFIT::kfit::PxLepton1,
					       *KINFIT::kfit::PyLepton1,
					       *KINFIT::kfit::PzLepton1,
					       *KINFIT::kfit::ELepton1,
					       *KINFIT::kfit::LabelLepton1,
					       *KINFIT::kfit::PxLepton2,
					       *KINFIT::kfit::PyLepton2,
					       *KINFIT::kfit::PzLepton2,
					       *KINFIT::kfit::ELepton2,
					       *KINFIT::kfit::LabelLepton2,
					       *KINFIT::kfit::PxBJet1,
					       *KINFIT::kfit::PyBJet1,
					       *KINFIT::kfit::PzBJet1,
					       *KINFIT::kfit::EBJet1,
					       *KINFIT::kfit::PxBJet2,
					       *KINFIT::kfit::PyBJet2,
					       *KINFIT::kfit::PzBJet2,
					       *KINFIT::kfit::EBJet2,
					       chi2W1, chi2W2, chi2Top1, chi2Top2,
					       chi2EtMissX, chi2EtMissY,
					       chi2BJet1Px, chi2BJet1Py, chi2BJet1Pz, chi2BJet1E,
					       chi2BJet2Px, chi2BJet2Py, chi2BJet2Pz, chi2BJet2E,
					       chi2Lepton1Px, chi2Lepton1Py, chi2Lepton1Pz, chi2Lepton1E,
					       chi2Lepton2Px, chi2Lepton2Py, chi2Lepton2Pz, chi2Lepton2E,
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
		       vp21.push_back(par[20]);
		       vp22.push_back(par[21]);
		       vp23.push_back(par[22]);
		       vp24.push_back(par[23]);
		       vchi2.push_back(lh);
		       
		       if( lh < 1. )
			 {
			    bp = 1;
			    break;
			 }
		    }				 
	       }
	  }
     }   

   double meanTimePDFBJetPx = std::accumulate(TimePDFBJetPx.begin(),TimePDFBJetPx.end(),0.0);
   meanTimePDFBJetPx = meanTimePDFBJetPx/TimePDFBJetPx.size();
   double meanTimePDFBJetPy = std::accumulate(TimePDFBJetPy.begin(),TimePDFBJetPy.end(),0.0);
   meanTimePDFBJetPy = meanTimePDFBJetPy/TimePDFBJetPy.size();
   double meanTimePDFBJetPz = std::accumulate(TimePDFBJetPz.begin(),TimePDFBJetPz.end(),0.0);
   meanTimePDFBJetPz = meanTimePDFBJetPz/TimePDFBJetPz.size();
   double meanTimePDFBJetE = std::accumulate(TimePDFBJetE.begin(),TimePDFBJetE.end(),0.0);
   meanTimePDFBJetE = meanTimePDFBJetE/TimePDFBJetE.size();
   double meanTimePDFMetPx = std::accumulate(TimePDFMetPx.begin(),TimePDFMetPx.end(),0.0);
   meanTimePDFMetPx = meanTimePDFMetPx/TimePDFMetPx.size();
   double meanTimePDFMetPy = std::accumulate(TimePDFMetPy.begin(),TimePDFMetPy.end(),0.0);
   meanTimePDFMetPy = meanTimePDFMetPy/TimePDFMetPy.size();
   
/*   std::cout << 
     "meanTimePDFBJetPx=" << meanTimePDFBJetPx << "ms, " <<
     "meanTimePDFBJetPy=" << meanTimePDFBJetPy << "ms, " <<
     "meanTimePDFBJetPz=" << meanTimePDFBJetPz << "ms, " <<
     "meanTimePDFBJetE=" << meanTimePDFBJetE << "ms, " <<
     "meanTimePDFMetPx=" << meanTimePDFMetPx << "ms, " <<
     "meanTimePDFMetPy=" << meanTimePDFMetPy << "ms, " <<
     std::endl;*/

   vp1.push_back(*EtMissX/2.);
   vp2.push_back(*EtMissY/2.);
   vp3.push_back(1);
   vp4.push_back(1);
   vp5.push_back(*EtMissX);
   vp6.push_back(*EtMissY);
   vp7.push_back(*WMassBW);
   vp8.push_back(*WMassBW);
   vp9.push_back(*PxBJet1);
   vp10.push_back(*PyBJet1);
   vp11.push_back(*PzBJet1);
   vp12.push_back(*EBJet1);
   vp13.push_back(*PxBJet2);
   vp14.push_back(*PyBJet2);
   vp15.push_back(*PzBJet2);
   vp16.push_back(*EBJet2);
   vp17.push_back(*PxLepton1);
   vp18.push_back(*PyLepton1);
   vp19.push_back(*PzLepton1);
   vp20.push_back(*ELepton1);
   vp21.push_back(*PxLepton2);
   vp22.push_back(*PyLepton2);
   vp23.push_back(*PzLepton2);
   vp24.push_back(*ELepton2);
   
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
   sortPermVector(vchi2,vp21);
   sortPermVector(vchi2,vp22);
   sortPermVector(vchi2,vp23);
   sortPermVector(vchi2,vp24);
   
   double gridTime = (double)(clock()-tStart)/CLOCKS_PER_SEC;
//   std::cout << "Grid stats = " << vchi2.size() << ", time = " << gridTime << "s" << std::endl;
}
