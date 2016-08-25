// KINFIT
#include "kinfit.h"

#include "TopTopLepLep.h"
#include "TopLep.h"
#include "TopTopLepHad.h"
#include "TopTopLepHbb.h"
#include "TopHLepbb.h"

#include <boost/bind.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

ClassImp(KINFIT::kfit)

KINFIT::TopTopLepLep *hypTopTopLepLep;
KINFIT::TopLep *hypTopLep;
KINFIT::TopTopLepHad *hypTopTopLepHad;
KINFIT::TopTopLepHbb *hypTopTopLepHbb;
KINFIT::TopHLepbb *hypTopHLepbb;

std::unique_ptr<double> KINFIT::kfit::CHISQ;
std::unique_ptr<double> KINFIT::kfit::CHISQT1;
std::unique_ptr<double> KINFIT::kfit::CHISQT2;
std::unique_ptr<double> KINFIT::kfit::EtMissX;
std::unique_ptr<double> KINFIT::kfit::EtMissY;
std::unique_ptr<double> KINFIT::kfit::PxLepton1;
std::unique_ptr<double> KINFIT::kfit::PyLepton1;
std::unique_ptr<double> KINFIT::kfit::PzLepton1;
std::unique_ptr<double> KINFIT::kfit::ELepton1;
std::unique_ptr<int> KINFIT::kfit::LabelLepton1;
std::unique_ptr<double> KINFIT::kfit::PxLepton2;
std::unique_ptr<double> KINFIT::kfit::PyLepton2;
std::unique_ptr<double> KINFIT::kfit::PzLepton2;
std::unique_ptr<double> KINFIT::kfit::ELepton2;
std::unique_ptr<int> KINFIT::kfit::LabelLepton2;
std::unique_ptr<double> KINFIT::kfit::PxBJet1;
std::unique_ptr<double> KINFIT::kfit::PyBJet1;
std::unique_ptr<double> KINFIT::kfit::PzBJet1;
std::unique_ptr<double> KINFIT::kfit::EBJet1;
std::unique_ptr<double> KINFIT::kfit::PxBJet2;
std::unique_ptr<double> KINFIT::kfit::PyBJet2;
std::unique_ptr<double> KINFIT::kfit::PzBJet2;
std::unique_ptr<double> KINFIT::kfit::EBJet2;
std::unique_ptr<double> KINFIT::kfit::PxBJet3;
std::unique_ptr<double> KINFIT::kfit::PyBJet3;
std::unique_ptr<double> KINFIT::kfit::PzBJet3;
std::unique_ptr<double> KINFIT::kfit::EBJet3;
std::unique_ptr<double> KINFIT::kfit::PxNonBJet1;
std::unique_ptr<double> KINFIT::kfit::PyNonBJet1;
std::unique_ptr<double> KINFIT::kfit::PzNonBJet1;
std::unique_ptr<double> KINFIT::kfit::ENonBJet1;
std::unique_ptr<double> KINFIT::kfit::PxNonBJet2;
std::unique_ptr<double> KINFIT::kfit::PyNonBJet2;
std::unique_ptr<double> KINFIT::kfit::PzNonBJet2;
std::unique_ptr<double> KINFIT::kfit::ENonBJet2;
std::unique_ptr<double> KINFIT::kfit::WMassGen1;
std::unique_ptr<double> KINFIT::kfit::WMassGen2;
std::unique_ptr<double> KINFIT::kfit::PzNu1;
std::unique_ptr<double> KINFIT::kfit::PzNu2;
std::unique_ptr<double> KINFIT::kfit::TopMass1;
std::unique_ptr<double> KINFIT::kfit::TopMass2;
std::unique_ptr<std::vector<double> > KINFIT::kfit::FitParam;
std::unique_ptr<std::vector<double> > KINFIT::kfit::ChiTerm;
std::unique_ptr<double> KINFIT::kfit::WMassBW;
std::unique_ptr<double> KINFIT::kfit::WRMSBW;
std::unique_ptr<double> KINFIT::kfit::MetXRMS;
std::unique_ptr<double> KINFIT::kfit::MetYRMS;
std::unique_ptr<double> KINFIT::kfit::TopMass;
std::unique_ptr<double> KINFIT::kfit::TopRMS;

std::unique_ptr<TF1> KINFIT::kfit::hPDFTopWMass;
std::unique_ptr<TF1> KINFIT::kfit::hPDFTopMass;
std::unique_ptr<TF1> KINFIT::kfit::hPDFHiggsMass;
std::unique_ptr<TF1> KINFIT::kfit::hPDFTopWHadMass;
std::unique_ptr<TF1> KINFIT::kfit::hPDFTopHadMass;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMetPx;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMetPy;
std::unique_ptr<TF1> KINFIT::kfit::hPDFBJetPx;
std::unique_ptr<TF1> KINFIT::kfit::hPDFBJetPy;
std::unique_ptr<TF1> KINFIT::kfit::hPDFBJetPz;
std::unique_ptr<TF1> KINFIT::kfit::hPDFBJetE;
std::unique_ptr<TF1> KINFIT::kfit::hPDFElecPx;
std::unique_ptr<TF1> KINFIT::kfit::hPDFElecPy;
std::unique_ptr<TF1> KINFIT::kfit::hPDFElecPz;
std::unique_ptr<TF1> KINFIT::kfit::hPDFElecE;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMuonPx;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMuonPy;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMuonPz;
std::unique_ptr<TF1> KINFIT::kfit::hPDFMuonE;
std::unique_ptr<TF1> KINFIT::kfit::hPDFNonBJetPx;
std::unique_ptr<TF1> KINFIT::kfit::hPDFNonBJetPy;
std::unique_ptr<TF1> KINFIT::kfit::hPDFNonBJetPz;
std::unique_ptr<TF1> KINFIT::kfit::hPDFNonBJetE;

float getProb(TF1 *hPDF,float var);

KINFIT::kfit::kfit()
{
   gErrorIgnoreLevel = 2000;

   WMassBW.reset(); WMassBW = std::unique_ptr<double>(new double());
   WRMSBW.reset(); WRMSBW = std::unique_ptr<double>(new double());   
   *WMassBW = 80.4;
   *WRMSBW = 2.1;

   NoBTag_ = 0;
   
   NToy_ = 20;
   NWRMS_ = 5;
   NMetRMS_ = 5;

   NBJetPxRMS_ = 3;
   NBJetPyRMS_ = 3;
   NBJetPzRMS_ = 3;
   NBJetERMS_ = 3;

   NNonBJetPxRMS_ = 3;
   NNonBJetPyRMS_ = 3;
   NNonBJetPzRMS_ = 3;
   NNonBJetERMS_ = 3;
   
   NElecPxRMS_ = 3;
   NElecPyRMS_ = 3;
   NElecPzRMS_ = 3;
   NElecERMS_ = 3;

   NMuonPxRMS_ = 3;
   NMuonPyRMS_ = 3;
   NMuonPzRMS_ = 3;
   NMuonERMS_ = 3;
   
   MetXRMS.reset(); MetXRMS = std::unique_ptr<double>(new double());
   MetYRMS.reset(); MetYRMS = std::unique_ptr<double>(new double());
   *MetXRMS = 80.;
   *MetYRMS = 80.;

   TopMass.reset(); TopMass = std::unique_ptr<double>(new double());
   TopRMS.reset(); TopRMS = std::unique_ptr<double>(new double());
   *TopMass = 172.5;
   *TopRMS = 2.;
   
   hPDFTopWMass.reset(); hPDFTopWMass = std::unique_ptr<TF1>(new TF1);
   hPDFTopMass.reset(); hPDFTopMass = std::unique_ptr<TF1>(new TF1);
   hPDFTopWHadMass.reset(); hPDFTopWHadMass = std::unique_ptr<TF1>(new TF1);
   hPDFTopHadMass.reset(); hPDFTopHadMass = std::unique_ptr<TF1>(new TF1);
   hPDFHiggsMass.reset(); hPDFHiggsMass = std::unique_ptr<TF1>(new TF1);
   hPDFMetPx.reset(); hPDFMetPx = std::unique_ptr<TF1>(new TF1);
   hPDFMetPy.reset(); hPDFMetPy = std::unique_ptr<TF1>(new TF1);
   hPDFBJetPx.reset(); hPDFBJetPx = std::unique_ptr<TF1>(new TF1);
   hPDFBJetPy.reset(); hPDFBJetPy = std::unique_ptr<TF1>(new TF1);
   hPDFBJetPz.reset(); hPDFBJetPz = std::unique_ptr<TF1>(new TF1);
   hPDFBJetE.reset(); hPDFBJetE = std::unique_ptr<TF1>(new TF1);
   hPDFElecPx.reset(); hPDFElecPx = std::unique_ptr<TF1>(new TF1);
   hPDFElecPy.reset(); hPDFElecPy = std::unique_ptr<TF1>(new TF1);
   hPDFElecPz.reset(); hPDFElecPz = std::unique_ptr<TF1>(new TF1);
   hPDFElecE.reset(); hPDFElecE = std::unique_ptr<TF1>(new TF1);
   hPDFMuonPx.reset(); hPDFMuonPx = std::unique_ptr<TF1>(new TF1);
   hPDFMuonPy.reset(); hPDFMuonPy = std::unique_ptr<TF1>(new TF1);
   hPDFMuonPz.reset(); hPDFMuonPz = std::unique_ptr<TF1>(new TF1);
   hPDFMuonE.reset(); hPDFMuonE = std::unique_ptr<TF1>(new TF1);
   hPDFNonBJetPx.reset(); hPDFNonBJetPx = std::unique_ptr<TF1>(new TF1);
   hPDFNonBJetPy.reset(); hPDFNonBJetPy = std::unique_ptr<TF1>(new TF1);
   hPDFNonBJetPz.reset(); hPDFNonBJetPz = std::unique_ptr<TF1>(new TF1);
   hPDFNonBJetE.reset(); hPDFNonBJetE = std::unique_ptr<TF1>(new TF1);
   
   MetPx_ = NULL;
   MetPy_ = NULL;
   
   WMass_ = NULL;
   WP_ = NULL;
   WPt_ = NULL;
   WEta_ = NULL;
   WPhi_ = NULL;
   WRap_ = NULL;
   TopMass_ = NULL;
   TopPt_ = NULL;
   TopP_ = NULL;
   TopEta_ = NULL;
   TopRap_ = NULL;

   chiPerm_ = NULL;
   chiTerm_ = NULL;
   parPerm_ = NULL;
   nuPxPerm_ = NULL;
   nuPyPerm_ = NULL;
   nuPzPerm_ = NULL;
   idxMin_ = NULL;
   
   TopTopLepLep_Electron1Idx = NULL;
   TopTopLepLep_Muon1Idx = NULL;
   TopTopLepLep_Electron2Idx = NULL;
   TopTopLepLep_Muon2Idx = NULL;
   TopTopLepLep_BJet1Idx = NULL;
   TopTopLepLep_BJet2Idx = NULL;
   
   TopLep_ElectronIdx = NULL;
   TopLep_MuonIdx = NULL;
   TopLep_BJetIdx = NULL;

   TopTopLepHad_ElectronIdx = NULL;
   TopTopLepHad_MuonIdx = NULL;
   TopTopLepHad_BJetLepIdx = NULL;
   TopTopLepHad_BJetHadIdx = NULL;
   TopTopLepHad_NonBJet1Idx = NULL;
   TopTopLepHad_NonBJet2Idx = NULL;
   
   TopTopLepHbb_ElectronIdx = NULL;
   TopTopLepHbb_MuonIdx = NULL;
   TopTopLepHbb_BJetLepIdx = NULL;
   TopTopLepHbb_NonBJetHadIdx = NULL;
   TopTopLepHbb_BJet1Idx = NULL;
   TopTopLepHbb_BJet2Idx = NULL;

   TopHLepbb_ElectronIdx = NULL;
   TopHLepbb_MuonIdx = NULL;
   TopHLepbb_BJetLepIdx = NULL;
   TopHLepbb_BJet1Idx = NULL;
   TopHLepbb_BJet2Idx = NULL;
	
   drTopTop_ = NULL;
   mTopTop_ = NULL;
   ptTopTop_ = NULL;
   pTopTop_ = NULL;
   etaTopTop_ = NULL;
   rapTopTop_ = NULL;
   phiTopTop_ = NULL;   
}

// Destructor
KINFIT::kfit::~kfit()
{
   if( hypoMode == TOPTOPLEPLEP ) delete hypTopTopLepLep;
   if( hypoMode == TOPLEP ) delete hypTopLep;
   if( hypoMode == TOPTOPLEPHAD ) delete hypTopTopLepHad;
   if( hypoMode == TOPTOPLEPHBB ) delete hypTopTopLepHbb;
   if( hypoMode == TOPHLEPBB ) delete hypTopHLepbb;
}

void KINFIT::kfit::Init(HYPO hypo)
{
   Reset();
   
   bool hypoFound = 0;   
   for(int i=TOPTOPLEPLEP;i!=TOPHLEPBB+1;i++)
     {
	HYPO cur = static_cast<HYPO>(i);
	if( hypo == cur )
	  {
	     hypoFound = 1;
	     break;	     
	  }	
     }   
   if( !hypoFound )
     {
	std::cout << "Hypo mode " << hypo << " is not defined" << std::endl;
	exit(1);
     }   
   else
     hypoMode = hypo;
   
   if( hypoMode == TOPTOPLEPLEP ) hypTopTopLepLep = (KINFIT::TopTopLepLep*)(this);
   if( hypoMode == TOPLEP ) hypTopLep = (KINFIT::TopLep*)(this);
   if( hypoMode == TOPTOPLEPHAD ) hypTopTopLepHad = (KINFIT::TopTopLepHad*)(this);
   if( hypoMode == TOPTOPLEPHBB ) hypTopTopLepHbb = (KINFIT::TopTopLepHbb*)(this);
   if( hypoMode == TOPHLEPBB ) hypTopHLepbb = (KINFIT::TopHLepbb*)(this);
   
   rnd = new TRandom3(666);
}

void KINFIT::kfit::SetPDF(std::string obj,std::string fileName,std::string hName)
{
   TFile *fPDF = TFile::Open(fileName.c_str());
   if( ! fPDF->IsOpen() )
     {
	std::cout << "File " << fileName << " does not exist" << std::endl;
	exit(1);
     }   
   if( ! fPDF->GetListOfKeys()->Contains(hName.c_str()) )
     {
	std::cout << "Histogram " << hName << " does not exist" << std::endl;
	exit(1);
     }
   
   TF1 *hPDF = (TF1*)fPDF->Get(hName.c_str());

   if( strcmp(obj.c_str(),"TopWMass") == 0 )
     {
	hPDFTopWMass = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFTopWMass"))) );
//	hPDFTopWMass->SetDirectory(0);
     }
   else if( strcmp(obj.c_str(),"TopMass") == 0 )
     {
	hPDFTopMass = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFTopMass"))) );
//	hPDFTopMass->SetDirectory(0);
     }   
   else if( strcmp(obj.c_str(),"MetPx") == 0 )
     {
	hPDFMetPx = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMetPx"))) );

	maxPDFMetPx = hPDFMetPx->GetMaximum();
	meanPDFMetPx = hPDFMetPx->GetMaximumX();
	sigmaPDFMetPx = fabs(hPDFMetPx->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"MetPy") == 0 )
     {
	hPDFMetPy = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMetPy"))) );

	maxPDFMetPy = hPDFMetPy->GetMaximum();
	meanPDFMetPy = hPDFMetPy->GetMaximumX();
	sigmaPDFMetPy = fabs(hPDFMetPy->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"BJetPx") == 0 )
     {
	hPDFBJetPx = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFBJetPx"))) );

	maxPDFBJetPx = hPDFBJetPx->GetMaximum();
	meanPDFBJetPx = hPDFBJetPx->GetMaximumX();
	sigmaPDFBJetPx = fabs(hPDFBJetPx->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"BJetPy") == 0 )
     {
	hPDFBJetPy = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFBJetPy"))) );

	maxPDFBJetPy = hPDFBJetPy->GetMaximum();
	meanPDFBJetPy = hPDFBJetPy->GetMaximumX();
	sigmaPDFBJetPy = fabs(hPDFBJetPy->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"BJetPz") == 0 )
     {
	hPDFBJetPz = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFBJetPz"))) );

	maxPDFBJetPz = hPDFBJetPz->GetMaximum();
	meanPDFBJetPz = hPDFBJetPz->GetMaximumX();
	sigmaPDFBJetPz = fabs(hPDFBJetPz->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"BJetE") == 0 )
     {
	hPDFBJetE = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFBJetE"))) );

	maxPDFBJetE = hPDFBJetE->GetMaximum();
	meanPDFBJetE = hPDFBJetE->GetMaximumX();
	sigmaPDFBJetE = fabs(hPDFBJetE->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"ElecPx") == 0 )
     {
	hPDFElecPx = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFElecPx"))) );

	maxPDFElecPx = hPDFElecPx->GetMaximum();
	meanPDFElecPx = hPDFElecPx->GetMaximumX();
	sigmaPDFElecPx = fabs(hPDFElecPx->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"ElecPy") == 0 )
     {
	hPDFElecPy = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFElecPy"))) );

	maxPDFElecPy = hPDFElecPy->GetMaximum();
	meanPDFElecPy = hPDFElecPy->GetMaximumX();
	sigmaPDFElecPy = fabs(hPDFElecPy->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"ElecPz") == 0 )
     {
	hPDFElecPz = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFElecPz"))) );

	maxPDFElecPz = hPDFElecPz->GetMaximum();
	meanPDFElecPz = hPDFElecPz->GetMaximumX();
	sigmaPDFElecPz = fabs(hPDFElecPz->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"ElecE") == 0 )
     {
	hPDFElecE = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFElecE"))) );

	maxPDFElecE = hPDFElecE->GetMaximum();
	meanPDFElecE = hPDFElecE->GetMaximumX();
	sigmaPDFElecE = fabs(hPDFElecE->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"MuonPx") == 0 )
     {
	hPDFMuonPx = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMuonPx"))) );

	maxPDFMuonPx = hPDFMuonPx->GetMaximum();
	meanPDFMuonPx = hPDFMuonPx->GetMaximumX();
	sigmaPDFMuonPx = fabs(hPDFMuonPx->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"MuonPy") == 0 )
     {
	hPDFMuonPy = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMuonPy"))) );

	maxPDFMuonPy = hPDFMuonPy->GetMaximum();
	meanPDFMuonPy = hPDFMuonPy->GetMaximumX();
	sigmaPDFMuonPy = fabs(hPDFMuonPy->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"MuonPz") == 0 )
     {
	hPDFMuonPz = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMuonPz"))) );

	maxPDFMuonPz = hPDFMuonPz->GetMaximum();
	meanPDFMuonPz = hPDFMuonPz->GetMaximumX();
	sigmaPDFMuonPz = fabs(hPDFMuonPz->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"MuonE") == 0 )
     {
	hPDFMuonE = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFMuonE"))) );

	maxPDFMuonE = hPDFMuonE->GetMaximum();
	meanPDFMuonE = hPDFMuonE->GetMaximumX();
	sigmaPDFMuonE = fabs(hPDFMuonE->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"NonBJetPx") == 0 )
     {
	hPDFNonBJetPx = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFNonBJetPx"))) );

	maxPDFNonBJetPx = hPDFNonBJetPx->GetMaximum();
	meanPDFNonBJetPx = hPDFNonBJetPx->GetMaximumX();
	sigmaPDFNonBJetPx = fabs(hPDFNonBJetPx->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"NonBJetPy") == 0 )
     {
	hPDFNonBJetPy = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFNonBJetPy"))) );

	maxPDFNonBJetPy = hPDFNonBJetPy->GetMaximum();
	meanPDFNonBJetPy = hPDFNonBJetPy->GetMaximumX();
	sigmaPDFNonBJetPy = fabs(hPDFNonBJetPy->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"NonBJetPz") == 0 )
     {
	hPDFNonBJetPz = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFNonBJetPz"))) );

	maxPDFNonBJetPz = hPDFNonBJetPz->GetMaximum();
	meanPDFNonBJetPz = hPDFNonBJetPz->GetMaximumX();
	sigmaPDFNonBJetPz = fabs(hPDFNonBJetPz->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"NonBJetE") == 0 )
     {
	hPDFNonBJetE = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFNonBJetE"))) );

	maxPDFNonBJetE = hPDFNonBJetE->GetMaximum();
	meanPDFNonBJetE = hPDFNonBJetE->GetMaximumX();
	sigmaPDFNonBJetE = fabs(hPDFNonBJetE->GetX(0.5));
     }   
   else if( strcmp(obj.c_str(),"TopWHadMass") == 0 )
     {
	hPDFTopWHadMass = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFTopWHadMass"))) );
     }
   else if( strcmp(obj.c_str(),"TopHadMass") == 0 )
     {
	hPDFTopHadMass = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFTopHadMass"))) );
     }   
   else if( strcmp(obj.c_str(),"HiggsMass") == 0 )
     {
	hPDFHiggsMass = std::unique_ptr<TF1>( (static_cast<TF1*>(hPDF->Clone("PDFHiggsMass"))) );
     }   
   else
     {
	std::cout << "PDF type is unknown" << std::endl;
	exit(1);
     }

   if( fPDF->IsOpen() ) fPDF->Close();
}

void KINFIT::kfit::Run()
{
   for(int i=0;i<nElectron;i++)
     {	
	LeptonLabel.push_back(0);
	LeptonIdx.push_back(i);
     }   
   
   nLepton += nElectron;
   LeptonPt.insert(LeptonPt.end(),ElectronPt.begin(),ElectronPt.end());
   LeptonEta.insert(LeptonEta.end(),ElectronEta.begin(),ElectronEta.end());
   LeptonPhi.insert(LeptonPhi.end(),ElectronPhi.begin(),ElectronPhi.end());
   LeptonE.insert(LeptonE.end(),ElectronE.begin(),ElectronE.end());
   LeptonPx.insert(LeptonPx.end(),ElectronPx.begin(),ElectronPx.end());
   LeptonPy.insert(LeptonPy.end(),ElectronPy.begin(),ElectronPy.end());
   LeptonPz.insert(LeptonPz.end(),ElectronPz.begin(),ElectronPz.end());

   for(int i=0;i<nMuon;i++)
     {	
	LeptonLabel.push_back(1);
	LeptonIdx.push_back(i);
     }   
   
   nLepton += nMuon;
   LeptonPt.insert(LeptonPt.end(),MuonPt.begin(),MuonPt.end());
   LeptonEta.insert(LeptonEta.end(),MuonEta.begin(),MuonEta.end());
   LeptonPhi.insert(LeptonPhi.end(),MuonPhi.begin(),MuonPhi.end());
   LeptonE.insert(LeptonE.end(),MuonE.begin(),MuonE.end());
   LeptonPx.insert(LeptonPx.end(),MuonPx.begin(),MuonPx.end());
   LeptonPy.insert(LeptonPy.end(),MuonPy.begin(),MuonPy.end());
   LeptonPz.insert(LeptonPz.end(),MuonPz.begin(),MuonPz.end());
   
   if( hypoMode == TOPTOPLEPLEP ) TopTopLepLep();
   if( hypoMode == TOPLEP ) TopLep();
   if( hypoMode == TOPTOPLEPHAD ) TopTopLepHad();
   if( hypoMode == TOPTOPLEPHBB ) TopTopLepHbb();
   if( hypoMode == TOPHLEPBB ) TopHLepbb();
   
   sortPermIndex();
   
   Reset();
}

float KINFIT::kfit::GetWMassGen(int index)
{
   if( index == 0 ) return WMassGen1_;
   else if( index == 1 ) return WMassGen2_;
   else
     {
	std::cout << "Incorrect index" << std::endl;
	exit(1);
     }	
}

float KINFIT::kfit::GetFitPar(int index)
{
   if( index == 0 ) return FitPar[0];
   else if( index == 1 ) return FitPar[1];
   else if( index == 2 ) return FitPar[2];
   else if( index == 3 ) return FitPar[3];
   else if( index == 4 ) return FitPar[4];
   else if( index == 5 ) return FitPar[5];
   else if( index == 6 ) return FitPar[6];
   else if( index == 7 ) return FitPar[7];
   else if( index == 8 ) return FitPar[8];
   else if( index == 9 ) return FitPar[9];
   else if( index == 10 ) return FitPar[10];
   else if( index == 11 ) return FitPar[11];
   else if( index == 12 ) return FitPar[12];
   else
     {
	std::cout << "Incorrect index" << std::endl;
	exit(1);
     }	
}

void KINFIT::kfit::Reset()
{
   CHISQ.reset(); CHISQ = std::unique_ptr<double>(new double(0));
   CHISQT1.reset(); CHISQT1 = std::unique_ptr<double>(new double(0));
   CHISQT2.reset(); CHISQT2 = std::unique_ptr<double>(new double(0));
   EtMissX.reset(); EtMissX = std::unique_ptr<double>(new double(0));
   EtMissY.reset(); EtMissY = std::unique_ptr<double>(new double(0));
   PxLepton1.reset(); PxLepton1 = std::unique_ptr<double>(new double(0));
   PyLepton1.reset(); PyLepton1 = std::unique_ptr<double>(new double(0));
   PzLepton1.reset(); PzLepton1 = std::unique_ptr<double>(new double(0));
   ELepton1.reset(); ELepton1 = std::unique_ptr<double>(new double(0));
   LabelLepton1.reset(); LabelLepton1 = std::unique_ptr<int>(new int(0));
   PxLepton2.reset(); PxLepton2 = std::unique_ptr<double>(new double(0));
   PyLepton2.reset(); PyLepton2 = std::unique_ptr<double>(new double(0));
   PzLepton2.reset(); PzLepton2 = std::unique_ptr<double>(new double(0));
   ELepton2.reset(); ELepton2 = std::unique_ptr<double>(new double(0));
   LabelLepton2.reset(); LabelLepton2 = std::unique_ptr<int>(new int(0));
   PxBJet1.reset(); PxBJet1 = std::unique_ptr<double>(new double(0));
   PyBJet1.reset(); PyBJet1 = std::unique_ptr<double>(new double(0));
   PzBJet1.reset(); PzBJet1 = std::unique_ptr<double>(new double(0));
   EBJet1.reset(); EBJet1 = std::unique_ptr<double>(new double(0));
   PxBJet2.reset(); PxBJet2 = std::unique_ptr<double>(new double(0));
   PyBJet2.reset(); PyBJet2 = std::unique_ptr<double>(new double(0));
   PzBJet2.reset(); PzBJet2 = std::unique_ptr<double>(new double(0));
   EBJet2.reset(); EBJet2 = std::unique_ptr<double>(new double(0));
   PxBJet3.reset(); PxBJet3 = std::unique_ptr<double>(new double(0));
   PyBJet3.reset(); PyBJet3 = std::unique_ptr<double>(new double(0));
   PzBJet3.reset(); PzBJet3 = std::unique_ptr<double>(new double(0));
   EBJet3.reset(); EBJet3 = std::unique_ptr<double>(new double(0));
   PxNonBJet1.reset(); PxNonBJet1 = std::unique_ptr<double>(new double(0));
   PyNonBJet1.reset(); PyNonBJet1 = std::unique_ptr<double>(new double(0));
   PzNonBJet1.reset(); PzNonBJet1 = std::unique_ptr<double>(new double(0));
   ENonBJet1.reset(); ENonBJet1 = std::unique_ptr<double>(new double(0));
   PxNonBJet2.reset(); PxNonBJet2 = std::unique_ptr<double>(new double(0));
   PyNonBJet2.reset(); PyNonBJet2 = std::unique_ptr<double>(new double(0));
   PzNonBJet2.reset(); PzNonBJet2 = std::unique_ptr<double>(new double(0));
   ENonBJet2.reset(); ENonBJet2 = std::unique_ptr<double>(new double(0));
   WMassGen1.reset(); WMassGen1 = std::unique_ptr<double>(new double(0));
   WMassGen2.reset(); WMassGen2 = std::unique_ptr<double>(new double(0));
   PzNu1.reset(); PzNu1 = std::unique_ptr<double>(new double(0));
   PzNu2.reset(); PzNu2 = std::unique_ptr<double>(new double(0));
   TopMass1.reset(); TopMass1 = std::unique_ptr<double>(new double(0));
   TopMass2.reset(); TopMass2 = std::unique_ptr<double>(new double(0));
   FitParam.reset(); FitParam = std::unique_ptr<std::vector<double> >(new std::vector<double>);
   ChiTerm.reset(); ChiTerm = std::unique_ptr<std::vector<double> >(new std::vector<double>);
   
   nLepton = 0;   
   LeptonPt.clear();
   LeptonEta.clear();
   LeptonPhi.clear();
   LeptonE.clear();
   LeptonPx.clear();
   LeptonPy.clear();
   LeptonPz.clear();
   LeptonLabel.clear();
   LeptonIdx.clear();

   nElectron = 0;
   ElectronPt.clear();
   ElectronEta.clear();
   ElectronPhi.clear();
   ElectronE.clear();
   ElectronPx.clear();
   ElectronPy.clear();
   ElectronPz.clear();

   nMuon = 0;
   MuonPt.clear();
   MuonEta.clear();
   MuonPhi.clear();
   MuonE.clear();
   MuonPx.clear();
   MuonPy.clear();
   MuonPz.clear();

   nBJet = 0;
   BJetPt.clear();
   BJetEta.clear();
   BJetPhi.clear();
   BJetE.clear();
   BJetPx.clear();
   BJetPy.clear();
   BJetPz.clear();

   nNonBJet = 0;
   NonBJetPt.clear();
   NonBJetEta.clear();
   NonBJetPhi.clear();
   NonBJetE.clear();
   NonBJetPx.clear();
   NonBJetPy.clear();
   NonBJetPz.clear();
}

void KINFIT::kfit::TopTopLepLep()
{
   hypTopTopLepLep->TopTopLepLepRun();
}

void KINFIT::kfit::TopLep()
{
   hypTopLep->TopLepRun();
}

void KINFIT::kfit::TopTopLepHad()
{
   hypTopTopLepHad->TopTopLepHadRun();
}

void KINFIT::kfit::TopTopLepHbb()
{
   hypTopTopLepHbb->TopTopLepHbbRun();
}

void KINFIT::kfit::TopHLepbb()
{
   hypTopHLepbb->TopHLepbbRun();
}

int KINFIT::kfit::GetIndex(OBJ objType,int idx)
{
   // TOPTOPLEPLEP
   if( objType == ELECTRON1_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_Electron1Idx);
     }   
   else if( objType == MUON1_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_Muon1Idx);
     }   
   else if( objType == ELECTRON2_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_Electron2Idx);
     }   
   else if( objType == MUON2_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_Muon2Idx);
     }   
   else if( objType == BJET1_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_BJet1Idx);
     }
   else if( objType == BJET2_TOPTOPLEPLEP )
     {
	return getPermIndex(idx,TopTopLepLep_BJet2Idx);
     }   

   // TOPLEP
   else if( objType == ELECTRON_TOPLEP )
     {
	return getPermIndex(idx,TopLep_ElectronIdx);
     }   
   else if( objType == MUON_TOPLEP )
     {
	return getPermIndex(idx,TopLep_MuonIdx);
     }   
   else if( objType == BJET_TOPLEP )
     {
	return getPermIndex(idx,TopLep_BJetIdx);
     }

   // TOPTOPLEPHAD
   else if( objType == ELECTRON_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_ElectronIdx);
     }   
   else if( objType == MUON_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_MuonIdx);
     }   
   else if( objType == BJETLEP_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_BJetLepIdx);
     }
   else if( objType == BJETHAD_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_BJetHadIdx);
     }
   else if( objType == NONBJET1_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_NonBJet1Idx);
     }
   else if( objType == NONBJET2_TOPTOPLEPHAD )
     {
	return getPermIndex(idx,TopTopLepHad_NonBJet2Idx);
     }

   // TOPTOPLEPHBB
   else if( objType == ELECTRON_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_ElectronIdx);
     }   
   else if( objType == MUON_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_MuonIdx);
     }   
   else if( objType == BJETLEP_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_BJetLepIdx);
     }
   else if( objType == NONBJETHAD_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_NonBJetHadIdx);
     }
   else if( objType == BJET1_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_BJet1Idx);
     }
   else if( objType == BJET2_TOPTOPLEPHBB )
     {
	return getPermIndex(idx,TopTopLepHbb_BJet2Idx);
     }

   // TOPHLEPBB
   else if( objType == ELECTRON_TOPHLEPBB )
     {
	return getPermIndex(idx,TopHLepbb_ElectronIdx);
     }   
   else if( objType == MUON_TOPHLEPBB )
     {
	return getPermIndex(idx,TopHLepbb_MuonIdx);
     }   
   else if( objType == BJETLEP_TOPHLEPBB )
     {
	return getPermIndex(idx,TopHLepbb_BJetLepIdx);
     }
   else if( objType == BJET1_TOPHLEPBB )
     {
	return getPermIndex(idx,TopHLepbb_BJet1Idx);
     }
   else if( objType == BJET2_TOPHLEPBB )
     {
	return getPermIndex(idx,TopHLepbb_BJet2Idx);
     }
}

int KINFIT::kfit::getPermIndex(int idx,int *arr)
{
   if( NPerm_ > 0 )
     {	     
	if( idx < 0 ) return arr[0];
	else if( idx < NPerm_ ) return arr[idxMin_[idx]];
	else
	  {
	     std::cout << "Max number of permutations is " << NPerm_ << std::endl;
	     exit(1);
	  }	     
     }
   else return -1;   
}

void KINFIT::kfit::sortPermIndex()
{
   float arrSort[NPerm_];
   
   for(int i=0;i<NPerm_;i++)
     {
	idxMin_[i] = i;
	arrSort[i] = chiPerm_[i];
     }
   
   for(int i=0;i<NPerm_;i++)
     {
	for(int j=i;j<NPerm_;j++)
	  {
	     if( arrSort[j] < arrSort[i] )
	       {
		  float temp = arrSort[i];
		  arrSort[i] = arrSort[j];
		  arrSort[j] = temp;

		  int tempIdx = idxMin_[i];
		  idxMin_[i] = idxMin_[j];
		  idxMin_[j] = tempIdx;
	       }
	  }	
     }   
}

void KINFIT::kfit::sortPermVector(std::vector<double> vRef,std::vector<double> &vSort)
{
   std::vector<std::pair<double,double> > vMap;
   
   int nRef = vRef.size();
   
   for(int i=0;i<nRef;i++)
     {
	vMap.push_back(std::make_pair(vRef[i],vSort[i]));
     }   
   
   std::sort(vMap.begin(),vMap.end(), 
	     boost::bind(&std::pair<double,double>::first,_1) <
	     boost::bind(&std::pair<double,double>::first,_2));
   
   vSort.clear();
   
   for(int i=0;i<nRef;i++)
     {
	vSort.push_back(vMap[i].second);
     }   
}

void KINFIT::kfit::close()
{
}

float KINFIT::kfit::getWmassBW(TRandom3 *rnd,float mWmean,float GammaW,float nSigma,double &proc)
{
   float mW = 0.;
   
   clock_t tStart = clock();
   
   float max = BW(mWmean,mWmean,GammaW);
   
   while(1)
     {
	float r1 = rnd->Rndm();
	float r2 = rnd->Rndm();

	mW = mWmean - nSigma*GammaW + 2*nSigma*GammaW*r1;
	if( mW <= 0 ) continue;
	if( BW(mW,mWmean,GammaW) > max*r2 ) break;
     }

   proc = (double)(clock()-tStart);
   
   return mW;      
}

float KINFIT::kfit::BW(float mW,float mWmean,float GammaW)
{
   float f = mWmean*mWmean*GammaW*GammaW/(pow(mW*mW-mWmean*mWmean,2)+mWmean*mWmean*GammaW*GammaW);
   
   return f;
}

bool KINFIT::kfit::getNuMom(float Wmass,float px_l,float py_l,float pz_l,float E_l,
			    float px_nu,float py_nu,float &pz_nu1,float &pz_nu2)
{
   bool hasSolution = 0.;
   
   float a = sqrt(px_l*px_l+py_l*py_l);
   float b = pz_l;
   float d = sqrt(px_nu*px_nu+py_nu*py_nu);
   float f = E_l;
   
   float c = Wmass*Wmass/2+px_l*px_nu+py_l*py_nu;
   
   float racine = c*c*b*b-a*a*(d*d*f*f-c*c);
   
   if(racine >= 0) 
     {
	hasSolution = 1;
	pz_nu1 = (c*b+sqrt(racine))/a/a;
	pz_nu2 = (c*b-sqrt(racine))/a/a;
     }   
   
   return hasSolution;
}

float getProb(TF1 *hPDF,float var)
{
   float prob = hPDF->Eval(var);
   
   return prob;
}

float KINFIT::kfit::getProbGaus(TF1 *hPDF,float max,float mean,float sigma,TRandom3 *rnd,float nSigma,double &proc)
{
   clock_t tStart = clock();

   float x = 0.;
  
   while(1)
     {
	float r1 = rnd->Rndm();
	float r2 = rnd->Rndm();

	x = mean - nSigma*sigma + 2*nSigma*sigma*r1;
//	if( x <= 0 ) continue;
	if( hPDF->Eval(x) > max*r2 ) break;
     }

   proc = (double)(clock()-tStart)*1000/CLOCKS_PER_SEC;

   return x;      
}

float KINFIT::kfit::getProbGausSub(TF1 *hPDF,TRandom3 *rnd,float nSigma)
{
   float max = hPDF->GetMaximum();
   float norm = hPDF->GetParameter(0);
   float mean = hPDF->GetParameter(1);
   float sigma = hPDF->GetParameter(2);

   float x = 0.;
   
   while(1)
     {
	float r1 = rnd->Rndm();
	float r2 = rnd->Rndm();

	x = mean - nSigma*sigma + 2*nSigma*sigma*r1;
//	if( x <= 0 ) continue;
	float val = norm*exp(-0.5*pow((x-mean)/sigma,2));
	if( val > max*r2 ) break;
     }
   
   return x;      
}

void KINFIT::kfit::SetBJet(std::vector<float> pt,
			   std::vector<float> eta,
			   std::vector<float> phi,
			   std::vector<float> E)
{
   nBJet = pt.size();
   BJetPt = pt;
   BJetEta = eta;
   BJetPhi = phi;
   BJetE = E;
   
   BJetPx.clear();
   BJetPy.clear();
   BJetPz.clear();
   
   for(int i=0;i<nBJet;i++)
     {
	float px = BJetPt[i]*cos(BJetPhi[i]);
	float py = BJetPt[i]*sin(BJetPhi[i]);
	float pz = BJetPt[i]*sinh(BJetEta[i]);
	
	BJetPx.push_back(px);
	BJetPy.push_back(py);
	BJetPz.push_back(pz);
     }      
}

void KINFIT::kfit::SetNonBJet(std::vector<float> pt,
			      std::vector<float> eta,
			      std::vector<float> phi,
			      std::vector<float> E)
{
   nNonBJet = pt.size();
   NonBJetPt = pt;
   NonBJetEta = eta;
   NonBJetPhi = phi;
   NonBJetE = E;

   NonBJetPx.clear();
   NonBJetPy.clear();
   NonBJetPz.clear();
   
   for(int i=0;i<nNonBJet;i++)
     {
	float px = NonBJetPt[i]*cos(NonBJetPhi[i]);
	float py = NonBJetPt[i]*sin(NonBJetPhi[i]);
	float pz = NonBJetPt[i]*sinh(NonBJetEta[i]);
	
	NonBJetPx.push_back(px);
	NonBJetPy.push_back(py);
	NonBJetPz.push_back(pz);
     }      
}

void KINFIT::kfit::SetElectron(std::vector<float> pt,
			       std::vector<float> eta,
			       std::vector<float> phi,
			       std::vector<float> E)
{
   nElectron = pt.size();
   ElectronPt = pt;
   ElectronEta = eta;
   ElectronPhi = phi;
   ElectronE = E;

   ElectronPx.clear();
   ElectronPy.clear();
   ElectronPz.clear();
   
   for(int i=0;i<nElectron;i++)
     {
	float px = ElectronPt[i]*cos(ElectronPhi[i]);
	float py = ElectronPt[i]*sin(ElectronPhi[i]);
	float pz = ElectronPt[i]*sinh(ElectronEta[i]);
	
	ElectronPx.push_back(px);
	ElectronPy.push_back(py);
	ElectronPz.push_back(pz);
     }   
}

void KINFIT::kfit::SetMuon(std::vector<float> pt,
			   std::vector<float> eta,
			   std::vector<float> phi,
			   std::vector<float> E)
{
   nMuon = pt.size();
   MuonPt = pt;
   MuonEta = eta;
   MuonPhi = phi;
   MuonE = E;

   MuonPx.clear();
   MuonPy.clear();
   MuonPz.clear();
   
   for(int i=0;i<nMuon;i++)
     {
	float px = MuonPt[i]*cos(MuonPhi[i]);
	float py = MuonPt[i]*sin(MuonPhi[i]);
	float pz = MuonPt[i]*sinh(MuonEta[i]);
	
	MuonPx.push_back(px);
	MuonPy.push_back(py);
	MuonPz.push_back(pz);
     }   
}

void KINFIT::kfit::SetMet(float px,
			  float py)
{
   MetPx = px;
   MetPy = py;
}

float KINFIT::kfit::getDeltaR(float eta1,float phi1,float eta2,float phi2)
{   
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > PI ) DeltaPhi -= 2.*PI;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

float KINFIT::kfit::getEta(float pt,float pz)
{   
   float theta = atan(pt/pz);
   if( theta < 0 ) theta += PI;
   return -log(tan(theta/2));
}

float KINFIT::kfit::getRap(float E,float pz)
{   
   return 0.5*log((E+pz) / (E-pz));
}

void KINFIT::kfit::checkPDF(TF1 *tf,std::string tfname)
{   
   if( strcmp(tf->GetName(),"") == 0 )
     {
	std::cout << "Can not find " << tfname << " PDF" << std::endl;
	exit(1);
     }   
}
