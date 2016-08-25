#include "kinfit.h"

#include "TH2D.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#define CSVv2L 0.460
#define CSVv2M 0.800
#define CSVv2T 0.935

float GetDeltaR(float eta1,float phi1,float eta2,float phi2);
float errf(float v1,float ve1,float v2,float ve2);

struct sortFunc
{
   bool operator()(const std::pair<float,int> &left, const std::pair<float,int> &right)
     {
	return left.first > right.first;
     }
};

int main(int argc, char *argv[])
{
   std::ifstream infile(argv[1]);
   std::string fname;
   if( infile.good() )
     {
	getline(infile,fname);
     }   
   infile.close();
   
   std::string foutName = std::string(argv[2])+".root";
   int nToys = atoi(argv[3]);   
   bool isSig = atoi(argv[4]);
   bool applyMVA = atoi(argv[5]);
   int nNonBJetMax = atoi(argv[6]);

   std::cout << "Signal = " << isSig << std::endl;
   std::cout << "Apply MVA = " << applyMVA << std::endl;
   std::cout << "Run " << nToys << " toys" << std::endl;

   TMVA::Reader* readerFullRecoTruth = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerFullRecoAll = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerFullRecoHighestCSVv2 = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerFullRecoCSVv2L = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerFullRecoCSVv2M = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerFullRecoCSVv2T = new TMVA::Reader("!Color:!Silent");

   TMVA::Reader* readerPartRecoTruth = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerPartRecoAll = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerPartRecoHighestCSVv2 = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerPartRecoCSVv2L = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerPartRecoCSVv2M = new TMVA::Reader("!Color:!Silent");
   TMVA::Reader* readerPartRecoCSVv2T = new TMVA::Reader("!Color:!Silent");

   // FullReco
   float MVAFullReco_TopHadRecMTruth;
   float MVAFullReco_TopLepRecMTruth;
   float MVAFullReco_TopLepTopHadRecDrTruth;
   float MVAFullReco_TopLepRecPtTruth;
   
   float MVAFullReco_TopHadRecMAll;
   float MVAFullReco_TopLepRecMAll;
   float MVAFullReco_TopLepTopHadRecDrAll;
   float MVAFullReco_TopLepRecPtAll;
   
   float MVAFullReco_TopHadRecMHighestCSVv2;
   float MVAFullReco_TopLepRecMHighestCSVv2;
   float MVAFullReco_TopLepTopHadRecDrHighestCSVv2;
   float MVAFullReco_TopLepRecPtHighestCSVv2;

   float MVAFullReco_TopHadRecMCSVv2L;
   float MVAFullReco_TopLepRecMCSVv2L;
   float MVAFullReco_TopLepTopHadRecDrCSVv2L;
   float MVAFullReco_TopLepRecPtCSVv2L;

   float MVAFullReco_TopHadRecMCSVv2M;
   float MVAFullReco_TopLepRecMCSVv2M;
   float MVAFullReco_TopLepTopHadRecDrCSVv2M;
   float MVAFullReco_TopLepRecPtCSVv2M;

   float MVAFullReco_TopHadRecMCSVv2T;
   float MVAFullReco_TopLepRecMCSVv2T;
   float MVAFullReco_TopLepTopHadRecDrCSVv2T;
   float MVAFullReco_TopLepRecPtCSVv2T;

   // PartReco
   float MVAPartReco_TopHadRecMTruth;
   float MVAPartReco_TopLepRecMTTruth;
   float MVAPartReco_TopLepTopHadRecDphiTTruth;
   float MVAPartReco_TopLepRecPtTTruth;

   float MVAPartReco_TopHadRecMAll;
   float MVAPartReco_TopLepRecMTAll;
   float MVAPartReco_TopLepTopHadRecDphiTAll;
   float MVAPartReco_TopLepRecPtTAll;
   
   float MVAPartReco_TopHadRecMHighestCSVv2;
   float MVAPartReco_TopLepRecMTHighestCSVv2;
   float MVAPartReco_TopLepTopHadRecDphiTHighestCSVv2;
   float MVAPartReco_TopLepRecPtTHighestCSVv2;

   float MVAPartReco_TopHadRecMCSVv2L;
   float MVAPartReco_TopLepRecMTCSVv2L;
   float MVAPartReco_TopLepTopHadRecDphiTCSVv2L;
   float MVAPartReco_TopLepRecPtTCSVv2L;

   float MVAPartReco_TopHadRecMCSVv2M;
   float MVAPartReco_TopLepRecMTCSVv2M;
   float MVAPartReco_TopLepTopHadRecDphiTCSVv2M;
   float MVAPartReco_TopLepRecPtTCSVv2M;

   float MVAPartReco_TopHadRecMCSVv2T;
   float MVAPartReco_TopLepRecMTCSVv2T;
   float MVAPartReco_TopLepTopHadRecDphiTCSVv2T;
   float MVAPartReco_TopLepRecPtTCSVv2T;
   
   if( applyMVA )
     {	
	readerFullRecoTruth->AddVariable("TopHadRecMTruth",&MVAFullReco_TopHadRecMTruth);
	readerFullRecoTruth->AddVariable("TopLepRecMTruth",&MVAFullReco_TopLepRecMTruth);
	readerFullRecoTruth->AddVariable("TopLepTopHadRecDrTruth",&MVAFullReco_TopLepTopHadRecDrTruth);
	readerFullRecoTruth->AddVariable("TopLepRecPtTruth",&MVAFullReco_TopLepRecPtTruth);

	std::string weightsFileFullRecoTruth = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoTruth_BDT.weights.xml";
	readerFullRecoTruth->BookMVA("BDTG method",weightsFileFullRecoTruth.c_str());

	readerFullRecoAll->AddVariable("TopHadRecMAll",&MVAFullReco_TopHadRecMAll);
	readerFullRecoAll->AddVariable("TopLepRecMAll",&MVAFullReco_TopLepRecMAll);
	readerFullRecoAll->AddVariable("TopLepTopHadRecDrAll",&MVAFullReco_TopLepTopHadRecDrAll);
	readerFullRecoAll->AddVariable("TopLepRecPtAll",&MVAFullReco_TopLepRecPtAll);

	std::string weightsFileFullRecoAll = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoAll_BDT.weights.xml";
	readerFullRecoAll->BookMVA("BDTG method",weightsFileFullRecoAll.c_str());
	
	readerFullRecoHighestCSVv2->AddVariable("TopHadRecMHighestCSVv2",&MVAFullReco_TopHadRecMHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("TopLepRecMHighestCSVv2",&MVAFullReco_TopLepRecMHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("TopLepTopHadRecDrHighestCSVv2",&MVAFullReco_TopLepTopHadRecDrHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("TopLepRecPtHighestCSVv2",&MVAFullReco_TopLepRecPtHighestCSVv2);

	std::string weightsFileFullRecoHighestCSVv2 = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoHighestCSVv2_BDT.weights.xml";
	readerFullRecoHighestCSVv2->BookMVA("BDTG method",weightsFileFullRecoHighestCSVv2.c_str());

	readerFullRecoCSVv2L->AddVariable("TopHadRecMCSVv2L",&MVAFullReco_TopHadRecMCSVv2L);
	readerFullRecoCSVv2L->AddVariable("TopLepRecMCSVv2L",&MVAFullReco_TopLepRecMCSVv2L);
	readerFullRecoCSVv2L->AddVariable("TopLepTopHadRecDrCSVv2L",&MVAFullReco_TopLepTopHadRecDrCSVv2L);
	readerFullRecoCSVv2L->AddVariable("TopLepRecPtCSVv2L",&MVAFullReco_TopLepRecPtCSVv2L);

	std::string weightsFileFullRecoCSVv2L = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoCSVv2L_BDT.weights.xml";
	readerFullRecoCSVv2L->BookMVA("BDTG method",weightsFileFullRecoCSVv2L.c_str());

	readerFullRecoCSVv2M->AddVariable("TopHadRecMCSVv2M",&MVAFullReco_TopHadRecMCSVv2M);
	readerFullRecoCSVv2M->AddVariable("TopLepRecMCSVv2M",&MVAFullReco_TopLepRecMCSVv2M);
	readerFullRecoCSVv2M->AddVariable("TopLepTopHadRecDrCSVv2M",&MVAFullReco_TopLepTopHadRecDrCSVv2M);
	readerFullRecoCSVv2M->AddVariable("TopLepRecPtCSVv2M",&MVAFullReco_TopLepRecPtCSVv2M);

	std::string weightsFileFullRecoCSVv2M = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoCSVv2M_BDT.weights.xml";
	readerFullRecoCSVv2M->BookMVA("BDTG method",weightsFileFullRecoCSVv2M.c_str());

	readerFullRecoCSVv2T->AddVariable("TopHadRecMCSVv2T",&MVAFullReco_TopHadRecMCSVv2T);
	readerFullRecoCSVv2T->AddVariable("TopLepRecMCSVv2T",&MVAFullReco_TopLepRecMCSVv2T);
	readerFullRecoCSVv2T->AddVariable("TopLepTopHadRecDrCSVv2T",&MVAFullReco_TopLepTopHadRecDrCSVv2T);
	readerFullRecoCSVv2T->AddVariable("TopLepRecPtCSVv2T",&MVAFullReco_TopLepRecPtCSVv2T);

	std::string weightsFileFullRecoCSVv2T = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAFullRecoCSVv2T_BDT.weights.xml";
	readerFullRecoCSVv2T->BookMVA("BDTG method",weightsFileFullRecoCSVv2T.c_str());

     
	readerPartRecoTruth->AddVariable("TopHadRecMTruth",&MVAPartReco_TopHadRecMTruth);
	readerPartRecoTruth->AddVariable("TopLepRecMTTruth",&MVAPartReco_TopLepRecMTTruth);
	readerPartRecoTruth->AddVariable("TopLepTopHadRecDphiTTruth",&MVAPartReco_TopLepTopHadRecDphiTTruth);
	readerPartRecoTruth->AddVariable("TopLepRecPtTTruth",&MVAPartReco_TopLepRecPtTTruth);

	std::string weightsFilePartRecoTruth = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoTruth_BDT.weights.xml";
	readerPartRecoTruth->BookMVA("BDTG method",weightsFilePartRecoTruth.c_str());

	readerPartRecoAll->AddVariable("TopHadRecMAll",&MVAPartReco_TopHadRecMAll);
	readerPartRecoAll->AddVariable("TopLepRecMTAll",&MVAPartReco_TopLepRecMTAll);
	readerPartRecoAll->AddVariable("TopLepTopHadRecDphiTAll",&MVAPartReco_TopLepTopHadRecDphiTAll);
	readerPartRecoAll->AddVariable("TopLepRecPtTAll",&MVAPartReco_TopLepRecPtTAll);

	std::string weightsFilePartRecoAll = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoAll_BDT.weights.xml";
	readerPartRecoAll->BookMVA("BDTG method",weightsFilePartRecoAll.c_str());
	
	readerPartRecoHighestCSVv2->AddVariable("TopHadRecMHighestCSVv2",&MVAPartReco_TopHadRecMHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("TopLepRecMTHighestCSVv2",&MVAPartReco_TopLepRecMTHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("TopLepTopHadRecDphiTHighestCSVv2",&MVAPartReco_TopLepTopHadRecDphiTHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("TopLepRecPtTHighestCSVv2",&MVAPartReco_TopLepRecPtTHighestCSVv2);

	std::string weightsFilePartRecoHighestCSVv2 = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoHighestCSVv2_BDT.weights.xml";
	readerPartRecoHighestCSVv2->BookMVA("BDTG method",weightsFilePartRecoHighestCSVv2.c_str());

	readerPartRecoCSVv2L->AddVariable("TopHadRecMCSVv2L",&MVAPartReco_TopHadRecMCSVv2L);
	readerPartRecoCSVv2L->AddVariable("TopLepRecMTCSVv2L",&MVAPartReco_TopLepRecMTCSVv2L);
	readerPartRecoCSVv2L->AddVariable("TopLepTopHadRecDphiTCSVv2L",&MVAPartReco_TopLepTopHadRecDphiTCSVv2L);
	readerPartRecoCSVv2L->AddVariable("TopLepRecPtTCSVv2L",&MVAPartReco_TopLepRecPtTCSVv2L);

	std::string weightsFilePartRecoCSVv2L = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoCSVv2L_BDT.weights.xml";
	readerPartRecoCSVv2L->BookMVA("BDTG method",weightsFilePartRecoCSVv2L.c_str());

	readerPartRecoCSVv2M->AddVariable("TopHadRecMCSVv2M",&MVAPartReco_TopHadRecMCSVv2M);
	readerPartRecoCSVv2M->AddVariable("TopLepRecMTCSVv2M",&MVAPartReco_TopLepRecMTCSVv2M);
	readerPartRecoCSVv2M->AddVariable("TopLepTopHadRecDphiTCSVv2M",&MVAPartReco_TopLepTopHadRecDphiTCSVv2M);
	readerPartRecoCSVv2M->AddVariable("TopLepRecPtTCSVv2M",&MVAPartReco_TopLepRecPtTCSVv2M);

	std::string weightsFilePartRecoCSVv2M = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoCSVv2M_BDT.weights.xml";
	readerPartRecoCSVv2M->BookMVA("BDTG method",weightsFilePartRecoCSVv2M.c_str());

	readerPartRecoCSVv2T->AddVariable("TopHadRecMCSVv2T",&MVAPartReco_TopHadRecMCSVv2T);
	readerPartRecoCSVv2T->AddVariable("TopLepRecMTCSVv2T",&MVAPartReco_TopLepRecMTCSVv2T);
	readerPartRecoCSVv2T->AddVariable("TopLepTopHadRecDphiTCSVv2T",&MVAPartReco_TopLepTopHadRecDphiTCSVv2T);
	readerPartRecoCSVv2T->AddVariable("TopLepRecPtTCSVv2T",&MVAPartReco_TopLepRecPtTCSVv2T);

	std::string weightsFilePartRecoCSVv2T = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHad/MVA/weights/TMVAPartRecoCSVv2T_BDT.weights.xml";
	readerPartRecoCSVv2T->BookMVA("BDTG method",weightsFilePartRecoCSVv2T.c_str());
     }
   
   KINFIT::kfit *kf = new KINFIT::kfit();

   kf->Init(TOPTOPLEPHAD);

   std::string pdfFileName = "/home-pbs/kskovpen/tHFCNC2016/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHad/pdf.root";
   kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopLepWM_Fit");
   kf->SetPDF("TopMass",pdfFileName.c_str(),"TopLepRecM_Fit");
   kf->SetPDF("TopWHadMass",pdfFileName.c_str(),"TopHadWRecM_Fit");
   kf->SetPDF("TopHadMass",pdfFileName.c_str(),"TopHadRecM_Fit");
   kf->SetPDF("MetPx",pdfFileName.c_str(),"dMetPx_Gaus");
   kf->SetPDF("MetPy",pdfFileName.c_str(),"dMetPy_Gaus");
   kf->SetPDF("BJetPx",pdfFileName.c_str(),"dBJetPx_Fit");
   kf->SetPDF("BJetPy",pdfFileName.c_str(),"dBJetPy_Fit");
   kf->SetPDF("BJetPz",pdfFileName.c_str(),"dBJetPz_Fit");
   kf->SetPDF("BJetE",pdfFileName.c_str(),"dBJetE_Fit");
   kf->SetPDF("NonBJetPx",pdfFileName.c_str(),"dNonBJetPx_Fit");
   kf->SetPDF("NonBJetPy",pdfFileName.c_str(),"dNonBJetPy_Fit");
   kf->SetPDF("NonBJetPz",pdfFileName.c_str(),"dNonBJetPz_Fit");
   kf->SetPDF("NonBJetE",pdfFileName.c_str(),"dNonBJetE_Fit");
   kf->SetPDF("ElecPx",pdfFileName.c_str(),"dElecPx_Fit");
   kf->SetPDF("ElecPy",pdfFileName.c_str(),"dElecPy_Fit");
   kf->SetPDF("ElecPz",pdfFileName.c_str(),"dElecPz_Fit");
   kf->SetPDF("ElecE",pdfFileName.c_str(),"dElecE_Fit");
   kf->SetPDF("MuonPx",pdfFileName.c_str(),"dMuonPx_Fit");
   kf->SetPDF("MuonPy",pdfFileName.c_str(),"dMuonPy_Fit");
   kf->SetPDF("MuonPz",pdfFileName.c_str(),"dMuonPz_Fit");
   kf->SetPDF("MuonE",pdfFileName.c_str(),"dMuonE_Fit");

   kf->SetNToy(nToys);

   TFile *f = TFile::Open(fname.c_str());
   
   TTree *tr = (TTree*)f->Get("trGEN");

   if( isSig )
     {
	bool Pass;
	float MetRecPx, MetRecPy;
	float TopLepWLepRecPt, TopLepWLepRecEta, TopLepWLepRecPhi, TopLepWLepRecE;
	int TopLepWLepId;
	float TopLepRecPt, TopLepRecEta, TopLepRecPhi, TopLepRecE, TopLepRecM;
	float TopLepBJetRecPt, TopLepBJetRecEta, TopLepBJetRecPhi, TopLepBJetRecE, TopLepBJetRecCSVv2;
	float TopHadRecPt, TopHadRecEta, TopHadRecPhi, TopHadRecE, TopHadRecM;
	float TopHadBJetRecPt, TopHadBJetRecEta, TopHadBJetRecPhi, TopHadBJetRecE, TopHadBJetRecCSVv2;
	float TopHadWRecPt, TopHadWRecEta, TopHadWRecPhi, TopHadWRecE, TopHadWRecM;
	float TopHadWNonBJet1RecPt, TopHadWNonBJet1RecEta, TopHadWNonBJet1RecPhi, TopHadWNonBJet1RecE, TopHadWNonBJet1RecCSVv2;
	float TopHadWNonBJet2RecPt, TopHadWNonBJet2RecEta, TopHadWNonBJet2RecPhi, TopHadWNonBJet2RecE, TopHadWNonBJet2RecCSVv2;
	std::vector<float> *OtherJetRecPt = new std::vector<float>();
	std::vector<float> *OtherJetRecEta = new std::vector<float>();
	std::vector<float> *OtherJetRecPhi = new std::vector<float>();
	std::vector<float> *OtherJetRecE = new std::vector<float>();
	std::vector<float> *OtherJetRecCSVv2 = new std::vector<float>();
	
	tr->SetBranchAddress("Pass",&Pass);
	tr->SetBranchAddress("MetRecPx",&MetRecPx);
	tr->SetBranchAddress("MetRecPy",&MetRecPy);
	tr->SetBranchAddress("TopLepWLepRecPt",&TopLepWLepRecPt);
	tr->SetBranchAddress("TopLepWLepRecEta",&TopLepWLepRecEta);
	tr->SetBranchAddress("TopLepWLepRecPhi",&TopLepWLepRecPhi);
	tr->SetBranchAddress("TopLepWLepRecE",&TopLepWLepRecE);
	tr->SetBranchAddress("TopLepWLepId",&TopLepWLepId);
	tr->SetBranchAddress("TopLepRecPt",&TopLepRecPt);
	tr->SetBranchAddress("TopLepRecEta",&TopLepRecEta);
	tr->SetBranchAddress("TopLepRecPhi",&TopLepRecPhi);
	tr->SetBranchAddress("TopLepRecE",&TopLepRecE);
	tr->SetBranchAddress("TopLepRecM",&TopLepRecM);
	tr->SetBranchAddress("TopLepBJetRecPt",&TopLepBJetRecPt);
	tr->SetBranchAddress("TopLepBJetRecEta",&TopLepBJetRecEta);
	tr->SetBranchAddress("TopLepBJetRecPhi",&TopLepBJetRecPhi);
	tr->SetBranchAddress("TopLepBJetRecE",&TopLepBJetRecE);
	tr->SetBranchAddress("TopLepBJetRecCSVv2",&TopLepBJetRecCSVv2);
	tr->SetBranchAddress("TopHadRecPt",&TopHadRecPt);
	tr->SetBranchAddress("TopHadRecEta",&TopHadRecEta);
	tr->SetBranchAddress("TopHadRecPhi",&TopHadRecPhi);
	tr->SetBranchAddress("TopHadRecE",&TopHadRecE);
	tr->SetBranchAddress("TopHadRecM",&TopHadRecM);
	tr->SetBranchAddress("TopHadBJetRecPt",&TopHadBJetRecPt);
	tr->SetBranchAddress("TopHadBJetRecEta",&TopHadBJetRecEta);
	tr->SetBranchAddress("TopHadBJetRecPhi",&TopHadBJetRecPhi);
	tr->SetBranchAddress("TopHadBJetRecE",&TopHadBJetRecE);
	tr->SetBranchAddress("TopHadBJetRecCSVv2",&TopHadBJetRecCSVv2);
	tr->SetBranchAddress("TopHadWRecPt",&TopHadWRecPt);
	tr->SetBranchAddress("TopHadWRecEta",&TopHadWRecEta);
	tr->SetBranchAddress("TopHadWRecPhi",&TopHadWRecPhi);
	tr->SetBranchAddress("TopHadWRecE",&TopHadWRecE);
	tr->SetBranchAddress("TopHadWRecM",&TopHadWRecM);
	tr->SetBranchAddress("TopHadWNonBJet1RecPt",&TopHadWNonBJet1RecPt);
	tr->SetBranchAddress("TopHadWNonBJet1RecEta",&TopHadWNonBJet1RecEta);
	tr->SetBranchAddress("TopHadWNonBJet1RecPhi",&TopHadWNonBJet1RecPhi);
	tr->SetBranchAddress("TopHadWNonBJet1RecE",&TopHadWNonBJet1RecE);
	tr->SetBranchAddress("TopHadWNonBJet1RecCSVv2",&TopHadWNonBJet1RecCSVv2);
	tr->SetBranchAddress("TopHadWNonBJet2RecPt",&TopHadWNonBJet2RecPt);
	tr->SetBranchAddress("TopHadWNonBJet2RecEta",&TopHadWNonBJet2RecEta);
	tr->SetBranchAddress("TopHadWNonBJet2RecPhi",&TopHadWNonBJet2RecPhi);
	tr->SetBranchAddress("TopHadWNonBJet2RecE",&TopHadWNonBJet2RecE);
	tr->SetBranchAddress("TopHadWNonBJet2RecCSVv2",&TopHadWNonBJet2RecCSVv2);
	tr->SetBranchAddress("OtherJetRecPt",&OtherJetRecPt);
	tr->SetBranchAddress("OtherJetRecEta",&OtherJetRecEta);
	tr->SetBranchAddress("OtherJetRecPhi",&OtherJetRecPhi);
	tr->SetBranchAddress("OtherJetRecE",&OtherJetRecE);
	tr->SetBranchAddress("OtherJetRecCSVv2",&OtherJetRecCSVv2);
	
	TFile *fout = new TFile(foutName.c_str(),"RECREATE");

	std::vector<float> MVAScoreTruth;
	std::vector<float> MVAScoreAll;
	std::vector<float> MVAScoreHighestCSVv2;
	std::vector<float> MVAScoreCSVv2L;
	std::vector<float> MVAScoreCSVv2M;
	std::vector<float> MVAScoreCSVv2T;
	
	std::vector<float> DiscTruth;
	std::vector<float> DiscAll;
	std::vector<float> DiscHighestCSVv2;
	std::vector<float> DiscCSVv2L;
	std::vector<float> DiscCSVv2M;
	std::vector<float> DiscCSVv2T;

	std::vector<float> MVADiscTruth;
	std::vector<float> MVADiscAll;
	std::vector<float> MVADiscHighestCSVv2;
	std::vector<float> MVADiscCSVv2L;
	std::vector<float> MVADiscCSVv2M;
	std::vector<float> MVADiscCSVv2T;

	std::vector<float> TopHadWRecMTruth;
	std::vector<float> TopHadWRecMAll;
	std::vector<float> TopHadWRecMHighestCSVv2;
	std::vector<float> TopHadWRecMCSVv2L;
	std::vector<float> TopHadWRecMCSVv2M;
	std::vector<float> TopHadWRecMCSVv2T;

	std::vector<float> MVATopHadRecMTruth;
	std::vector<float> MVATopHadRecMAll;
	std::vector<float> MVATopHadRecMHighestCSVv2;
	std::vector<float> MVATopHadRecMCSVv2L;
	std::vector<float> MVATopHadRecMCSVv2M;
	std::vector<float> MVATopHadRecMCSVv2T;
	
	std::vector<float> TopLepRecMTruth;
	std::vector<float> TopLepRecMAll;
	std::vector<float> TopLepRecMHighestCSVv2;
	std::vector<float> TopLepRecMCSVv2L;
	std::vector<float> TopLepRecMCSVv2M;
	std::vector<float> TopLepRecMCSVv2T;

	std::vector<float> MVATopLepRecMTruth;
	std::vector<float> MVATopLepRecMAll;
	std::vector<float> MVATopLepRecMHighestCSVv2;
	std::vector<float> MVATopLepRecMCSVv2L;
	std::vector<float> MVATopLepRecMCSVv2M;
	std::vector<float> MVATopLepRecMCSVv2T;
	
	std::vector<float> TopHadRecMTruth;
	std::vector<float> TopHadRecMAll;
	std::vector<float> TopHadRecMHighestCSVv2;
	std::vector<float> TopHadRecMCSVv2L;
	std::vector<float> TopHadRecMCSVv2M;
	std::vector<float> TopHadRecMCSVv2T;

	std::vector<float> TopLepTopHadRecDrTruth;
	std::vector<float> TopLepTopHadRecDrAll;
	std::vector<float> TopLepTopHadRecDrHighestCSVv2;
	std::vector<float> TopLepTopHadRecDrCSVv2L;
	std::vector<float> TopLepTopHadRecDrCSVv2M;
	std::vector<float> TopLepTopHadRecDrCSVv2T;

	std::vector<float> TopLepRecPtTruth;
	std::vector<float> TopLepRecPtAll;
	std::vector<float> TopLepRecPtHighestCSVv2;
	std::vector<float> TopLepRecPtCSVv2L;
	std::vector<float> TopLepRecPtCSVv2M;
	std::vector<float> TopLepRecPtCSVv2T;

	std::vector<float> TopHadRecPtTruth;
	std::vector<float> TopHadRecPtAll;
	std::vector<float> TopHadRecPtHighestCSVv2;
	std::vector<float> TopHadRecPtCSVv2L;
	std::vector<float> TopHadRecPtCSVv2M;
	std::vector<float> TopHadRecPtCSVv2T;
	
	std::vector<bool> MatchTruth;
	std::vector<bool> MatchAll;
	std::vector<bool> MatchHighestCSVv2;
	std::vector<bool> MatchCSVv2L;
	std::vector<bool> MatchCSVv2M;
	std::vector<bool> MatchCSVv2T;

	std::vector<bool> MatchBJetTruth;
	std::vector<bool> MatchBJetAll;
	std::vector<bool> MatchBJetHighestCSVv2;
	std::vector<bool> MatchBJetCSVv2L;
	std::vector<bool> MatchBJetCSVv2M;
	std::vector<bool> MatchBJetCSVv2T;
	
	std::vector<bool> MatchMVATruth;
	std::vector<bool> MatchMVAAll;
	std::vector<bool> MatchMVAHighestCSVv2;
	std::vector<bool> MatchMVACSVv2L;
	std::vector<bool> MatchMVACSVv2M;
	std::vector<bool> MatchMVACSVv2T;

	std::vector<bool> MatchBJetMVATruth;
	std::vector<bool> MatchBJetMVAAll;
	std::vector<bool> MatchBJetMVAHighestCSVv2;
	std::vector<bool> MatchBJetMVACSVv2L;
	std::vector<bool> MatchBJetMVACSVv2M;
	std::vector<bool> MatchBJetMVACSVv2T;

	std::vector<float> TopLepRecMTTruth;
	std::vector<float> TopLepRecMTAll;
	std::vector<float> TopLepRecMTHighestCSVv2;
	std::vector<float> TopLepRecMTCSVv2L;
	std::vector<float> TopLepRecMTCSVv2M;
	std::vector<float> TopLepRecMTCSVv2T;

	std::vector<float> TopLepTopHadRecDphiTTruth;
	std::vector<float> TopLepTopHadRecDphiTAll;
	std::vector<float> TopLepTopHadRecDphiTHighestCSVv2;
	std::vector<float> TopLepTopHadRecDphiTCSVv2L;
	std::vector<float> TopLepTopHadRecDphiTCSVv2M;
	std::vector<float> TopLepTopHadRecDphiTCSVv2T;

	std::vector<float> TopLepRecPtTTruth;
	std::vector<float> TopLepRecPtTAll;
	std::vector<float> TopLepRecPtTHighestCSVv2;
	std::vector<float> TopLepRecPtTCSVv2L;
	std::vector<float> TopLepRecPtTCSVv2M;
	std::vector<float> TopLepRecPtTCSVv2T;
	
	int NBJetTruth;
	int NBJetAll;
	int NBJetHighestCSVv2;
	int NBJetCSVv2L;
	int NBJetCSVv2M;
	int NBJetCSVv2T;

	int NNonBJetTruth;
	int NNonBJetAll;
	int NNonBJetHighestCSVv2;
	int NNonBJetCSVv2L;
	int NNonBJetCSVv2M;
	int NNonBJetCSVv2T;
	
	int NPermTruth;
	int NPermAll;
	int NPermHighestCSVv2;
	int NPermCSVv2L;
	int NPermCSVv2M;
	int NPermCSVv2T;

	float nEventsTruth = 0;
	float nEventsAll = 0;
	float nEventsHighestCSVv2 = 0;
	float nEventsCSVv2L = 0;
	float nEventsCSVv2M = 0;
	float nEventsCSVv2T = 0;
	
	float nMatchTruth = 0;
	float nMatchAll = 0;
	float nMatchHighestCSVv2 = 0;
	float nMatchCSVv2L = 0;
	float nMatchCSVv2M = 0;
	float nMatchCSVv2T = 0;

	float nMatchBJetTruth = 0;
	float nMatchBJetAll = 0;
	float nMatchBJetHighestCSVv2 = 0;
	float nMatchBJetCSVv2L = 0;
	float nMatchBJetCSVv2M = 0;
	float nMatchBJetCSVv2T = 0;
	
	float nMatchMVATruth = 0;
	float nMatchMVAAll = 0;
	float nMatchMVAHighestCSVv2 = 0;
	float nMatchMVACSVv2L = 0;
	float nMatchMVACSVv2M = 0;
	float nMatchMVACSVv2T = 0;

	float nMatchBJetMVATruth = 0;
	float nMatchBJetMVAAll = 0;
	float nMatchBJetMVAHighestCSVv2 = 0;
	float nMatchBJetMVACSVv2L = 0;
	float nMatchBJetMVACSVv2M = 0;
	float nMatchBJetMVACSVv2T = 0;
	
	float nSelTruth = 0;
	float nSelAll = 0;
	float nSelHighestCSVv2 = 0;
	float nSelCSVv2L = 0;
	float nSelCSVv2M = 0;
	float nSelCSVv2T = 0;

	float nNoSolutionTruth = 0;
	float nNoSolutionAll = 0;
	float nNoSolutionHighestCSVv2 = 0;
	float nNoSolutionCSVv2L = 0;
	float nNoSolutionCSVv2M = 0;
	float nNoSolutionCSVv2T = 0;
	
	TTree *trFIT = new TTree("trFIT","trFIT");

	trFIT->Branch("MVAScoreTruth","std::vector<float>",&MVAScoreTruth);
	trFIT->Branch("MVAScoreAll","std::vector<float>",&MVAScoreAll);
	trFIT->Branch("MVAScoreHighestCSVv2","std::vector<float>",&MVAScoreHighestCSVv2);
	trFIT->Branch("MVAScoreCSVv2L","std::vector<float>",&MVAScoreCSVv2L);
	trFIT->Branch("MVAScoreCSVv2M","std::vector<float>",&MVAScoreCSVv2M);
	trFIT->Branch("MVAScoreCSVv2T","std::vector<float>",&MVAScoreCSVv2T);
	trFIT->Branch("DiscTruth","std::vector<float>",&DiscTruth);
	trFIT->Branch("DiscAll","std::vector<float>",&DiscAll);
	trFIT->Branch("DiscHighestCSVv2","std::vector<float>",&DiscHighestCSVv2);
	trFIT->Branch("DiscCSVv2L","std::vector<float>",&DiscCSVv2L);
	trFIT->Branch("DiscCSVv2M","std::vector<float>",&DiscCSVv2M);
	trFIT->Branch("DiscCSVv2T","std::vector<float>",&DiscCSVv2T);
	trFIT->Branch("MVADiscTruth","std::vector<float>",&MVADiscTruth);
	trFIT->Branch("MVADiscAll","std::vector<float>",&MVADiscAll);
	trFIT->Branch("MVADiscHighestCSVv2","std::vector<float>",&MVADiscHighestCSVv2);
	trFIT->Branch("MVADiscCSVv2L","std::vector<float>",&MVADiscCSVv2L);
	trFIT->Branch("MVADiscCSVv2M","std::vector<float>",&MVADiscCSVv2M);
	trFIT->Branch("MVADiscCSVv2T","std::vector<float>",&MVADiscCSVv2T);
	trFIT->Branch("TopHadWRecMTruth","std::vector<float>",&TopHadWRecMTruth);
	trFIT->Branch("TopHadWRecMAll","std::vector<float>",&TopHadWRecMAll);
	trFIT->Branch("TopHadWRecMHighestCSVv2","std::vector<float>",&TopHadWRecMHighestCSVv2);
	trFIT->Branch("TopHadWRecMCSVv2L","std::vector<float>",&TopHadWRecMCSVv2L);
	trFIT->Branch("TopHadWRecMCSVv2M","std::vector<float>",&TopHadWRecMCSVv2M);
	trFIT->Branch("TopHadWRecMCSVv2T","std::vector<float>",&TopHadWRecMCSVv2T);
	trFIT->Branch("MVATopHadRecMTruth","std::vector<float>",&MVATopHadRecMTruth);
	trFIT->Branch("MVATopHadRecMAll","std::vector<float>",&MVATopHadRecMAll);
	trFIT->Branch("MVATopHadRecMHighestCSVv2","std::vector<float>",&MVATopHadRecMHighestCSVv2);
	trFIT->Branch("MVATopHadRecMCSVv2L","std::vector<float>",&MVATopHadRecMCSVv2L);
	trFIT->Branch("MVATopHadRecMCSVv2M","std::vector<float>",&MVATopHadRecMCSVv2M);
	trFIT->Branch("MVATopHadRecMCSVv2T","std::vector<float>",&MVATopHadRecMCSVv2T);
	trFIT->Branch("TopLepRecMTruth","std::vector<float>",&TopLepRecMTruth);
	trFIT->Branch("TopLepRecMAll","std::vector<float>",&TopLepRecMAll);
	trFIT->Branch("TopLepRecMHighestCSVv2","std::vector<float>",&TopLepRecMHighestCSVv2);
	trFIT->Branch("TopLepRecMCSVv2L","std::vector<float>",&TopLepRecMCSVv2L);
	trFIT->Branch("TopLepRecMCSVv2M","std::vector<float>",&TopLepRecMCSVv2M);
	trFIT->Branch("TopLepRecMCSVv2T","std::vector<float>",&TopLepRecMCSVv2T);
	trFIT->Branch("MVATopLepRecMTruth","std::vector<float>",&MVATopLepRecMTruth);
	trFIT->Branch("MVATopLepRecMAll","std::vector<float>",&MVATopLepRecMAll);
	trFIT->Branch("MVATopLepRecMHighestCSVv2","std::vector<float>",&MVATopLepRecMHighestCSVv2);
	trFIT->Branch("MVATopLepRecMCSVv2L","std::vector<float>",&MVATopLepRecMCSVv2L);
	trFIT->Branch("MVATopLepRecMCSVv2M","std::vector<float>",&MVATopLepRecMCSVv2M);
	trFIT->Branch("MVATopLepRecMCSVv2T","std::vector<float>",&MVATopLepRecMCSVv2T);
	trFIT->Branch("TopHadRecMTruth","std::vector<float>",&TopHadRecMTruth);
	trFIT->Branch("TopHadRecMAll","std::vector<float>",&TopHadRecMAll);
	trFIT->Branch("TopHadRecMHighestCSVv2","std::vector<float>",&TopHadRecMHighestCSVv2);
	trFIT->Branch("TopHadRecMCSVv2L","std::vector<float>",&TopHadRecMCSVv2L);
	trFIT->Branch("TopHadRecMCSVv2M","std::vector<float>",&TopHadRecMCSVv2M);
	trFIT->Branch("TopHadRecMCSVv2T","std::vector<float>",&TopHadRecMCSVv2T);
	trFIT->Branch("TopLepTopHadRecDrTruth","std::vector<float>",&TopLepTopHadRecDrTruth);
	trFIT->Branch("TopLepTopHadRecDrAll","std::vector<float>",&TopLepTopHadRecDrAll);
	trFIT->Branch("TopLepTopHadRecDrHighestCSVv2","std::vector<float>",&TopLepTopHadRecDrHighestCSVv2);
	trFIT->Branch("TopLepTopHadRecDrCSVv2L","std::vector<float>",&TopLepTopHadRecDrCSVv2L);
	trFIT->Branch("TopLepTopHadRecDrCSVv2M","std::vector<float>",&TopLepTopHadRecDrCSVv2M);
	trFIT->Branch("TopLepTopHadRecDrCSVv2T","std::vector<float>",&TopLepTopHadRecDrCSVv2T);
	trFIT->Branch("TopLepRecPtTruth","std::vector<float>",&TopLepRecPtTruth);
	trFIT->Branch("TopLepRecPtAll","std::vector<float>",&TopLepRecPtAll);
	trFIT->Branch("TopLepRecPtHighestCSVv2","std::vector<float>",&TopLepRecPtHighestCSVv2);
	trFIT->Branch("TopLepRecPtCSVv2L","std::vector<float>",&TopLepRecPtCSVv2L);
	trFIT->Branch("TopLepRecPtCSVv2M","std::vector<float>",&TopLepRecPtCSVv2M);
	trFIT->Branch("TopLepRecPtCSVv2T","std::vector<float>",&TopLepRecPtCSVv2T);
	trFIT->Branch("TopHadRecPtTruth","std::vector<float>",&TopHadRecPtTruth);
	trFIT->Branch("TopHadRecPtAll","std::vector<float>",&TopHadRecPtAll);
	trFIT->Branch("TopHadRecPtHighestCSVv2","std::vector<float>",&TopHadRecPtHighestCSVv2);
	trFIT->Branch("TopHadRecPtCSVv2L","std::vector<float>",&TopHadRecPtCSVv2L);
	trFIT->Branch("TopHadRecPtCSVv2M","std::vector<float>",&TopHadRecPtCSVv2M);
	trFIT->Branch("TopHadRecPtCSVv2T","std::vector<float>",&TopHadRecPtCSVv2T);
	trFIT->Branch("MatchTruth","std::vector<bool>",&MatchTruth);
	trFIT->Branch("MatchAll","std::vector<bool>",&MatchAll);
	trFIT->Branch("MatchHighestCSVv2","std::vector<bool>",&MatchHighestCSVv2);
	trFIT->Branch("MatchCSVv2L","std::vector<bool>",&MatchCSVv2L);
	trFIT->Branch("MatchCSVv2M","std::vector<bool>",&MatchCSVv2M);
	trFIT->Branch("MatchCSVv2T","std::vector<bool>",&MatchCSVv2T);
	trFIT->Branch("MatchBJetTruth","std::vector<bool>",&MatchBJetTruth);
	trFIT->Branch("MatchBJetAll","std::vector<bool>",&MatchBJetAll);
	trFIT->Branch("MatchBJetHighestCSVv2","std::vector<bool>",&MatchBJetHighestCSVv2);
	trFIT->Branch("MatchBJetCSVv2L","std::vector<bool>",&MatchBJetCSVv2L);
	trFIT->Branch("MatchBJetCSVv2M","std::vector<bool>",&MatchBJetCSVv2M);
	trFIT->Branch("MatchBJetCSVv2T","std::vector<bool>",&MatchBJetCSVv2T);
	trFIT->Branch("MatchMVATruth","std::vector<bool>",&MatchMVATruth);
	trFIT->Branch("MatchMVAAll","std::vector<bool>",&MatchMVAAll);
	trFIT->Branch("MatchMVAHighestCSVv2","std::vector<bool>",&MatchMVAHighestCSVv2);
	trFIT->Branch("MatchMVACSVv2L","std::vector<bool>",&MatchMVACSVv2L);
	trFIT->Branch("MatchMVACSVv2M","std::vector<bool>",&MatchMVACSVv2M);
	trFIT->Branch("MatchMVACSVv2T","std::vector<bool>",&MatchMVACSVv2T);
	trFIT->Branch("MatchBJetMVATruth","std::vector<bool>",&MatchBJetMVATruth);
	trFIT->Branch("MatchBJetMVAAll","std::vector<bool>",&MatchBJetMVAAll);
	trFIT->Branch("MatchBJetMVAHighestCSVv2","std::vector<bool>",&MatchBJetMVAHighestCSVv2);
	trFIT->Branch("MatchBJetMVACSVv2L","std::vector<bool>",&MatchBJetMVACSVv2L);
	trFIT->Branch("MatchBJetMVACSVv2M","std::vector<bool>",&MatchBJetMVACSVv2M);
	trFIT->Branch("MatchBJetMVACSVv2T","std::vector<bool>",&MatchBJetMVACSVv2T);

	trFIT->Branch("TopLepRecMTTruth","std::vector<float>",&TopLepRecMTTruth);
	trFIT->Branch("TopLepRecMTAll","std::vector<float>",&TopLepRecMTAll);
	trFIT->Branch("TopLepRecMTHighestCSVv2","std::vector<float>",&TopLepRecMTHighestCSVv2);
	trFIT->Branch("TopLepRecMTCSVv2L","std::vector<float>",&TopLepRecMTCSVv2L);
	trFIT->Branch("TopLepRecMTCSVv2M","std::vector<float>",&TopLepRecMTCSVv2M);
	trFIT->Branch("TopLepRecMTCSVv2T","std::vector<float>",&TopLepRecMTCSVv2T);

	trFIT->Branch("TopLepTopHadRecDphiTTruth","std::vector<float>",&TopLepTopHadRecDphiTTruth);
	trFIT->Branch("TopLepTopHadRecDphiTAll","std::vector<float>",&TopLepTopHadRecDphiTAll);
	trFIT->Branch("TopLepTopHadRecDphiTHighestCSVv2","std::vector<float>",&TopLepTopHadRecDphiTHighestCSVv2);
	trFIT->Branch("TopLepTopHadRecDphiTCSVv2L","std::vector<float>",&TopLepTopHadRecDphiTCSVv2L);
	trFIT->Branch("TopLepTopHadRecDphiTCSVv2M","std::vector<float>",&TopLepTopHadRecDphiTCSVv2M);
	trFIT->Branch("TopLepTopHadRecDphiTCSVv2T","std::vector<float>",&TopLepTopHadRecDphiTCSVv2T);
	
	trFIT->Branch("TopLepRecPtTTruth","std::vector<float>",&TopLepRecPtTTruth);
	trFIT->Branch("TopLepRecPtTAll","std::vector<float>",&TopLepRecPtTAll);
	trFIT->Branch("TopLepRecPtTHighestCSVv2","std::vector<float>",&TopLepRecPtTHighestCSVv2);
	trFIT->Branch("TopLepRecPtTCSVv2L","std::vector<float>",&TopLepRecPtTCSVv2L);
	trFIT->Branch("TopLepRecPtTCSVv2M","std::vector<float>",&TopLepRecPtTCSVv2M);
	trFIT->Branch("TopLepRecPtTCSVv2T","std::vector<float>",&TopLepRecPtTCSVv2T);
	
	trFIT->Branch("NBJetTruth",&NBJetTruth,"NBJetTruth/I");
	trFIT->Branch("NBJetAll",&NBJetAll,"NBJetAll/I");
	trFIT->Branch("NBJetHighestCSVv2",&NBJetHighestCSVv2,"NBJetHighestCSVv2/I");
	trFIT->Branch("NBJetCSVv2L",&NBJetCSVv2L,"NBJetCSVv2L/I");
	trFIT->Branch("NBJetCSVv2M",&NBJetCSVv2M,"NBJetCSVv2M/I");
	trFIT->Branch("NBJetCSVv2T",&NBJetCSVv2T,"NBJetCSVv2T/I");
	trFIT->Branch("NNonBJetTruth",&NNonBJetTruth,"NNonBJetTruth/I");
	trFIT->Branch("NNonBJetAll",&NNonBJetAll,"NNonBJetAll/I");
	trFIT->Branch("NNonBJetHighestCSVv2",&NNonBJetHighestCSVv2,"NNonBJetHighestCSVv2/I");
	trFIT->Branch("NNonBJetCSVv2L",&NNonBJetCSVv2L,"NNonBJetCSVv2L/I");
	trFIT->Branch("NNonBJetCSVv2M",&NNonBJetCSVv2M,"NNonBJetCSVv2M/I");
	trFIT->Branch("NNonBJetCSVv2T",&NNonBJetCSVv2T,"NNonBJetCSVv2T/I");
	trFIT->Branch("NPermTruth",&NPermTruth,"NPermTruth/I");
	trFIT->Branch("NPermAll",&NPermAll,"NPermAll/I");
	trFIT->Branch("NPermHighestCSVv2",&NPermHighestCSVv2,"NPermHighestCSVv2/I");
	trFIT->Branch("NPermCSVv2L",&NPermCSVv2L,"NPermCSVv2L/I");
	trFIT->Branch("NPermCSVv2M",&NPermCSVv2M,"NPermCSVv2M/I");
	trFIT->Branch("NPermCSVv2T",&NPermCSVv2T,"NPermCSVv2T/I");

	TTree *trEFF = new TTree("trEFF","trEFF");

	trEFF->Branch("nEventsTruth",&nEventsTruth,"nEventsTruth/F");
	trEFF->Branch("nEventsAll",&nEventsAll,"nEventsAll/F");
	trEFF->Branch("nEventsHighestCSVv2",&nEventsHighestCSVv2,"nEventsHighestCSVv2/F");
	trEFF->Branch("nEventsCSVv2L",&nEventsCSVv2L,"nEventsCSVv2L/F");
	trEFF->Branch("nEventsCSVv2M",&nEventsCSVv2M,"nEventsCSVv2M/F");
	trEFF->Branch("nEventsCSVv2T",&nEventsCSVv2T,"nEventsCSVv2T/F");

	trEFF->Branch("nMatchTruth",&nMatchTruth,"nMatchTruth/F");
	trEFF->Branch("nMatchAll",&nMatchAll,"nMatchAll/F");
	trEFF->Branch("nMatchHighestCSVv2",&nMatchHighestCSVv2,"nMatchHighestCSVv2/F");
	trEFF->Branch("nMatchCSVv2L",&nMatchCSVv2L,"nMatchCSVv2L/F");
	trEFF->Branch("nMatchCSVv2M",&nMatchCSVv2M,"nMatchCSVv2M/F");
	trEFF->Branch("nMatchCSVv2T",&nMatchCSVv2T,"nMatchCSVv2T/F");

	trEFF->Branch("nMatchBJetTruth",&nMatchBJetTruth,"nMatchBJetTruth/F");
	trEFF->Branch("nMatchBJetAll",&nMatchBJetAll,"nMatchBJetAll/F");
	trEFF->Branch("nMatchBJetHighestCSVv2",&nMatchBJetHighestCSVv2,"nMatchBJetHighestCSVv2/F");
	trEFF->Branch("nMatchBJetCSVv2L",&nMatchBJetCSVv2L,"nMatchBJetCSVv2L/F");
	trEFF->Branch("nMatchBJetCSVv2M",&nMatchBJetCSVv2M,"nMatchBJetCSVv2M/F");
	trEFF->Branch("nMatchBJetCSVv2T",&nMatchBJetCSVv2T,"nMatchBJetCSVv2T/F");
	
	trEFF->Branch("nMatchMVATruth",&nMatchMVATruth,"nMatchMVATruth/F");
	trEFF->Branch("nMatchMVAAll",&nMatchMVAAll,"nMatchMVAAll/F");
	trEFF->Branch("nMatchMVAHighestCSVv2",&nMatchMVAHighestCSVv2,"nMatchMVAHighestCSVv2/F");
	trEFF->Branch("nMatchMVACSVv2L",&nMatchMVACSVv2L,"nMatchMVACSVv2L/F");
	trEFF->Branch("nMatchMVACSVv2M",&nMatchMVACSVv2M,"nMatchMVACSVv2M/F");
	trEFF->Branch("nMatchMVACSVv2T",&nMatchMVACSVv2T,"nMatchMVACSVv2T/F");

	trEFF->Branch("nMatchBJetMVATruth",&nMatchBJetMVATruth,"nMatchBJetMVATruth/F");
	trEFF->Branch("nMatchBJetMVAAll",&nMatchBJetMVAAll,"nMatchBJetMVAAll/F");
	trEFF->Branch("nMatchBJetMVAHighestCSVv2",&nMatchBJetMVAHighestCSVv2,"nMatchBJetMVAHighestCSVv2/F");
	trEFF->Branch("nMatchBJetMVACSVv2L",&nMatchBJetMVACSVv2L,"nMatchBJetMVACSVv2L/F");
	trEFF->Branch("nMatchBJetMVACSVv2M",&nMatchBJetMVACSVv2M,"nMatchBJetMVACSVv2M/F");
	trEFF->Branch("nMatchBJetMVACSVv2T",&nMatchBJetMVACSVv2T,"nMatchBJetMVACSVv2T/F");
	
	trEFF->Branch("nSelTruth",&nSelTruth,"nSelTruth/F");
	trEFF->Branch("nSelAll",&nSelAll,"nSelAll/F");
	trEFF->Branch("nSelHighestCSVv2",&nSelHighestCSVv2,"nSelHighestCSVv2/F");
	trEFF->Branch("nSelCSVv2L",&nSelCSVv2L,"nSelCSVv2L/F");
	trEFF->Branch("nSelCSVv2M",&nSelCSVv2M,"nSelCSVv2M/F");
	trEFF->Branch("nSelCSVv2T",&nSelCSVv2T,"nSelCSVv2T/F");

	trEFF->Branch("nNoSolutionTruth",&nNoSolutionTruth,"nNoSolutionTruth/F");
	trEFF->Branch("nNoSolutionAll",&nNoSolutionAll,"nNoSolutionAll/F");
	trEFF->Branch("nNoSolutionHighestCSVv2",&nNoSolutionHighestCSVv2,"nNoSolutionHighestCSVv2/F");
	trEFF->Branch("nNoSolutionCSVv2L",&nNoSolutionCSVv2L,"nNoSolutionCSVv2L/F");
	trEFF->Branch("nNoSolutionCSVv2M",&nNoSolutionCSVv2M,"nNoSolutionCSVv2M/F");
	trEFF->Branch("nNoSolutionCSVv2T",&nNoSolutionCSVv2T,"nNoSolutionCSVv2T/F");
	
	int nev = tr->GetEntries();
	std::cout << "Total number of events = " << nev << std::endl;
	
	int nMax = 2000;
//	int nMax = 1000;
//	int nMax = 100;
	
	for(int i=0;i<nev;i++)
	  {
	     if( nMax >=0 && i > nMax ) break;
	     
	     std::cout << "Event #" << i << std::endl;
	     
	     tr->GetEntry(i);
	     
	     if( !Pass ) continue;

	     MVAScoreTruth.clear();
	     MVAScoreAll.clear();
	     MVAScoreHighestCSVv2.clear();
	     MVAScoreCSVv2L.clear();
	     MVAScoreCSVv2M.clear();
	     MVAScoreCSVv2T.clear();
	     
	     DiscTruth.clear();
	     DiscAll.clear();
	     DiscHighestCSVv2.clear();
	     DiscCSVv2L.clear();
	     DiscCSVv2M.clear();
	     DiscCSVv2T.clear();

	     MVADiscTruth.clear();
	     MVADiscAll.clear();
	     MVADiscHighestCSVv2.clear();
	     MVADiscCSVv2L.clear();
	     MVADiscCSVv2M.clear();
	     MVADiscCSVv2T.clear();
	     
	     TopHadWRecMTruth.clear();
	     TopHadWRecMAll.clear();
	     TopHadWRecMHighestCSVv2.clear();
	     TopHadWRecMCSVv2L.clear();
	     TopHadWRecMCSVv2M.clear();
	     TopHadWRecMCSVv2T.clear();

	     MVATopHadRecMTruth.clear();
	     MVATopHadRecMAll.clear();
	     MVATopHadRecMHighestCSVv2.clear();
	     MVATopHadRecMCSVv2L.clear();
	     MVATopHadRecMCSVv2M.clear();
	     MVATopHadRecMCSVv2T.clear();
	     
	     TopLepRecMTruth.clear();
	     TopLepRecMAll.clear();
	     TopLepRecMHighestCSVv2.clear();
	     TopLepRecMCSVv2L.clear();
	     TopLepRecMCSVv2M.clear();
	     TopLepRecMCSVv2T.clear();

	     MVATopLepRecMTruth.clear();
	     MVATopLepRecMAll.clear();
	     MVATopLepRecMHighestCSVv2.clear();
	     MVATopLepRecMCSVv2L.clear();
	     MVATopLepRecMCSVv2M.clear();
	     MVATopLepRecMCSVv2T.clear();
	     
	     TopHadRecMTruth.clear();
	     TopHadRecMAll.clear();
	     TopHadRecMHighestCSVv2.clear();
	     TopHadRecMCSVv2L.clear();
	     TopHadRecMCSVv2M.clear();
	     TopHadRecMCSVv2T.clear();

	     TopLepTopHadRecDrTruth.clear();
	     TopLepTopHadRecDrAll.clear();
	     TopLepTopHadRecDrHighestCSVv2.clear();
	     TopLepTopHadRecDrCSVv2L.clear();
	     TopLepTopHadRecDrCSVv2M.clear();
	     TopLepTopHadRecDrCSVv2T.clear();

	     TopLepRecPtTruth.clear();
	     TopLepRecPtAll.clear();
	     TopLepRecPtHighestCSVv2.clear();
	     TopLepRecPtCSVv2L.clear();
	     TopLepRecPtCSVv2M.clear();
	     TopLepRecPtCSVv2T.clear();

	     TopHadRecPtTruth.clear();
	     TopHadRecPtAll.clear();
	     TopHadRecPtHighestCSVv2.clear();
	     TopHadRecPtCSVv2L.clear();
	     TopHadRecPtCSVv2M.clear();
	     TopHadRecPtCSVv2T.clear();
	     
	     MatchTruth.clear();
	     MatchAll.clear();
	     MatchHighestCSVv2.clear();
	     MatchCSVv2L.clear();
	     MatchCSVv2M.clear();
	     MatchCSVv2T.clear();

	     MatchBJetTruth.clear();
	     MatchBJetAll.clear();
	     MatchBJetHighestCSVv2.clear();
	     MatchBJetCSVv2L.clear();
	     MatchBJetCSVv2M.clear();
	     MatchBJetCSVv2T.clear();
	     
	     MatchMVATruth.clear();
	     MatchMVAAll.clear();
	     MatchMVAHighestCSVv2.clear();
	     MatchMVACSVv2L.clear();
	     MatchMVACSVv2M.clear();
	     MatchMVACSVv2T.clear();

	     MatchBJetMVATruth.clear();
	     MatchBJetMVAAll.clear();
	     MatchBJetMVAHighestCSVv2.clear();
	     MatchBJetMVACSVv2L.clear();
	     MatchBJetMVACSVv2M.clear();
	     MatchBJetMVACSVv2T.clear();

	     TopLepRecMTTruth.clear();
	     TopLepRecMTAll.clear();
	     TopLepRecMTHighestCSVv2.clear();
	     TopLepRecMTCSVv2L.clear();
	     TopLepRecMTCSVv2M.clear();
	     TopLepRecMTCSVv2T.clear();

	     TopLepTopHadRecDphiTTruth.clear();
	     TopLepTopHadRecDphiTAll.clear();
	     TopLepTopHadRecDphiTHighestCSVv2.clear();
	     TopLepTopHadRecDphiTCSVv2L.clear();
	     TopLepTopHadRecDphiTCSVv2M.clear();
	     TopLepTopHadRecDphiTCSVv2T.clear();

	     TopLepRecPtTTruth.clear();
	     TopLepRecPtTAll.clear();
	     TopLepRecPtTHighestCSVv2.clear();
	     TopLepRecPtTCSVv2L.clear();
	     TopLepRecPtTCSVv2M.clear();
	     TopLepRecPtTCSVv2T.clear();
	     
	     for(int it=0;it<6;it++)
	       {
		  std::vector<float> BJetPt;
		  std::vector<float> BJetEta;
		  std::vector<float> BJetPhi;
		  std::vector<float> BJetE;
		  std::vector<float> BJetCSVv2;
		  
		  std::vector<float> NonBJetPt;
		  std::vector<float> NonBJetEta;
		  std::vector<float> NonBJetPhi;
		  std::vector<float> NonBJetE;
		  std::vector<float> NonBJetCSVv2;
		  
		  std::vector<float> ElectronPt;
		  std::vector<float> ElectronEta;
		  std::vector<float> ElectronPhi;
		  std::vector<float> ElectronE;
		  
		  std::vector<float> MuonPt;
		  std::vector<float> MuonEta;
		  std::vector<float> MuonPhi;
		  std::vector<float> MuonE;	     
		  
		  if( it == 0 ) // Truth
		    {		  
		       BJetPt.push_back(TopLepBJetRecPt);
		       BJetEta.push_back(TopLepBJetRecEta);
		       BJetPhi.push_back(TopLepBJetRecPhi);
		       BJetE.push_back(TopLepBJetRecE);
		       BJetCSVv2.push_back(TopLepBJetRecCSVv2);
		       
		       NonBJetPt.push_back(TopHadWNonBJet1RecPt);
		       NonBJetEta.push_back(TopHadWNonBJet1RecEta);
		       NonBJetPhi.push_back(TopHadWNonBJet1RecPhi);
		       NonBJetE.push_back(TopHadWNonBJet1RecE);
		       NonBJetCSVv2.push_back(TopHadWNonBJet1RecCSVv2);
		       
		       NonBJetPt.push_back(TopHadWNonBJet2RecPt);
		       NonBJetEta.push_back(TopHadWNonBJet2RecEta);
		       NonBJetPhi.push_back(TopHadWNonBJet2RecPhi);
		       NonBJetE.push_back(TopHadWNonBJet2RecE);
		       NonBJetCSVv2.push_back(TopHadWNonBJet2RecCSVv2);
		       
		       BJetPt.push_back(TopHadBJetRecPt);
		       BJetEta.push_back(TopHadBJetRecEta);
		       BJetPhi.push_back(TopHadBJetRecPhi);
		       BJetE.push_back(TopHadBJetRecE);
		       BJetCSVv2.push_back(TopHadBJetRecCSVv2);
		    }
		  else if( it != 5 )
		    {
		       std::vector<std::pair<float,int> > JetCSVv2;
		       JetCSVv2.push_back(std::make_pair(TopLepBJetRecCSVv2,0));
		       JetCSVv2.push_back(std::make_pair(TopHadBJetRecCSVv2,1));
		       JetCSVv2.push_back(std::make_pair(TopHadWNonBJet1RecCSVv2,2));
		       JetCSVv2.push_back(std::make_pair(TopHadWNonBJet2RecCSVv2,3));
		       
		       for(int ij=0;ij<OtherJetRecPt->size();ij++)
			 JetCSVv2.push_back(std::make_pair(OtherJetRecCSVv2->at(ij),4+ij));
		       
		       std::sort(JetCSVv2.begin(),JetCSVv2.end(),sortFunc());
		       
		       for(int ij=0;ij<JetCSVv2.size();ij++)
			 {		       	
			    int idx = JetCSVv2.at(ij).second;
			    int nBJet = BJetPt.size();
			    
			    if( idx == 0 )
			      {
				 bool passBTag = 1;
				 if( TopLepBJetRecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( TopLepBJetRecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( TopLepBJetRecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {
				      BJetPt.push_back(TopLepBJetRecPt);
				      BJetEta.push_back(TopLepBJetRecEta);
				      BJetPhi.push_back(TopLepBJetRecPhi);
				      BJetE.push_back(TopLepBJetRecE);
				      BJetCSVv2.push_back(TopLepBJetRecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(TopLepBJetRecPt);
				      NonBJetEta.push_back(TopLepBJetRecEta);
				      NonBJetPhi.push_back(TopLepBJetRecPhi);
				      NonBJetE.push_back(TopLepBJetRecE);
				      NonBJetCSVv2.push_back(TopLepBJetRecCSVv2);
				   }			    
			      }		  
			    else if( idx == 1 )
			      {		       
				 bool passBTag = 1;
				 if( TopHadBJetRecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( TopHadBJetRecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( TopHadBJetRecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {				 
				      BJetPt.push_back(TopHadBJetRecPt);
				      BJetEta.push_back(TopHadBJetRecEta);
				      BJetPhi.push_back(TopHadBJetRecPhi);
				      BJetE.push_back(TopHadBJetRecE);
				      BJetCSVv2.push_back(TopHadBJetRecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(TopHadBJetRecPt);
				      NonBJetEta.push_back(TopHadBJetRecEta);
				      NonBJetPhi.push_back(TopHadBJetRecPhi);
				      NonBJetE.push_back(TopHadBJetRecE);
				      NonBJetCSVv2.push_back(TopHadBJetRecCSVv2);
				   }			    
			      }
			    else if( idx == 2 )
			      {
				 bool passBTag = 1;
				 if( TopHadWNonBJet1RecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( TopHadWNonBJet1RecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( TopHadWNonBJet1RecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {				 
				      BJetPt.push_back(TopHadWNonBJet1RecPt);
				      BJetEta.push_back(TopHadWNonBJet1RecEta);
				      BJetPhi.push_back(TopHadWNonBJet1RecPhi);
				      BJetE.push_back(TopHadWNonBJet1RecE);
				      BJetCSVv2.push_back(TopHadWNonBJet1RecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(TopHadWNonBJet1RecPt);
				      NonBJetEta.push_back(TopHadWNonBJet1RecEta);
				      NonBJetPhi.push_back(TopHadWNonBJet1RecPhi);
				      NonBJetE.push_back(TopHadWNonBJet1RecE);
				      NonBJetCSVv2.push_back(TopHadWNonBJet1RecCSVv2);
				   }			    
			      }		  
			    else if( idx == 3 )
			      {		       
				 bool passBTag = 1;
				 if( TopHadWNonBJet2RecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( TopHadWNonBJet2RecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( TopHadWNonBJet2RecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {
				      BJetPt.push_back(TopHadWNonBJet2RecPt);
				      BJetEta.push_back(TopHadWNonBJet2RecEta);
				      BJetPhi.push_back(TopHadWNonBJet2RecPhi);
				      BJetE.push_back(TopHadWNonBJet2RecE);
				      BJetCSVv2.push_back(TopHadWNonBJet2RecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(TopHadWNonBJet2RecPt);
				      NonBJetEta.push_back(TopHadWNonBJet2RecEta);
				      NonBJetPhi.push_back(TopHadWNonBJet2RecPhi);
				      NonBJetE.push_back(TopHadWNonBJet2RecE);
				      NonBJetCSVv2.push_back(TopHadWNonBJet2RecCSVv2);
				   }
			      }
			    else
			      {
				 bool passBTag = 1;
				 if( OtherJetRecCSVv2->at(idx-4) < CSVv2L && it == 2 ) passBTag = 0;
				 if( OtherJetRecCSVv2->at(idx-4) < CSVv2M && it == 3 ) passBTag = 0;
				 if( OtherJetRecCSVv2->at(idx-4) < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {				 
				      BJetPt.push_back(OtherJetRecPt->at(idx-4));
				      BJetEta.push_back(OtherJetRecEta->at(idx-4));
				      BJetPhi.push_back(OtherJetRecPhi->at(idx-4));
				      BJetE.push_back(OtherJetRecE->at(idx-4));
				      BJetCSVv2.push_back(OtherJetRecCSVv2->at(idx-4));
				   }
				 else
				   {
				      NonBJetPt.push_back(OtherJetRecPt->at(idx-4));
				      NonBJetEta.push_back(OtherJetRecEta->at(idx-4));
				      NonBJetPhi.push_back(OtherJetRecPhi->at(idx-4));
				      NonBJetE.push_back(OtherJetRecE->at(idx-4));
				      NonBJetCSVv2.push_back(OtherJetRecCSVv2->at(idx-4));
				   }
			      }		       
			 }		  
		    }
		  else
		    {
		       BJetPt.push_back(TopLepBJetRecPt);
		       BJetEta.push_back(TopLepBJetRecEta);
		       BJetPhi.push_back(TopLepBJetRecPhi);
		       BJetE.push_back(TopLepBJetRecE);
		       BJetCSVv2.push_back(TopLepBJetRecCSVv2);
		       
		       NonBJetPt.push_back(TopLepBJetRecPt);
		       NonBJetEta.push_back(TopLepBJetRecEta);
		       NonBJetPhi.push_back(TopLepBJetRecPhi);
		       NonBJetE.push_back(TopLepBJetRecE);
		       NonBJetCSVv2.push_back(TopLepBJetRecCSVv2);

		       BJetPt.push_back(TopHadWNonBJet1RecPt);
		       BJetEta.push_back(TopHadWNonBJet1RecEta);
		       BJetPhi.push_back(TopHadWNonBJet1RecPhi);
		       BJetE.push_back(TopHadWNonBJet1RecE);
		       BJetCSVv2.push_back(TopHadWNonBJet1RecCSVv2);
		       
		       NonBJetPt.push_back(TopHadWNonBJet1RecPt);
		       NonBJetEta.push_back(TopHadWNonBJet1RecEta);
		       NonBJetPhi.push_back(TopHadWNonBJet1RecPhi);
		       NonBJetE.push_back(TopHadWNonBJet1RecE);
		       NonBJetCSVv2.push_back(TopHadWNonBJet1RecCSVv2);

		       BJetPt.push_back(TopHadWNonBJet2RecPt);
		       BJetEta.push_back(TopHadWNonBJet2RecEta);
		       BJetPhi.push_back(TopHadWNonBJet2RecPhi);
		       BJetE.push_back(TopHadWNonBJet2RecE);
		       BJetCSVv2.push_back(TopHadWNonBJet2RecCSVv2);
		       
		       NonBJetPt.push_back(TopHadWNonBJet2RecPt);
		       NonBJetEta.push_back(TopHadWNonBJet2RecEta);
		       NonBJetPhi.push_back(TopHadWNonBJet2RecPhi);
		       NonBJetE.push_back(TopHadWNonBJet2RecE);
		       NonBJetCSVv2.push_back(TopHadWNonBJet2RecCSVv2);
		       
		       BJetPt.push_back(TopHadBJetRecPt);
		       BJetEta.push_back(TopHadBJetRecEta);
		       BJetPhi.push_back(TopHadBJetRecPhi);
		       BJetE.push_back(TopHadBJetRecE);
		       BJetCSVv2.push_back(TopHadBJetRecCSVv2);
		       
		       NonBJetPt.push_back(TopHadBJetRecPt);
		       NonBJetEta.push_back(TopHadBJetRecEta);
		       NonBJetPhi.push_back(TopHadBJetRecPhi);
		       NonBJetE.push_back(TopHadBJetRecE);
		       NonBJetCSVv2.push_back(TopHadBJetRecCSVv2);
		       
		       for(int ij=0;ij<OtherJetRecPt->size();ij++)
			 {
			    BJetPt.push_back(OtherJetRecPt->at(ij));
			    BJetEta.push_back(OtherJetRecEta->at(ij));
			    BJetPhi.push_back(OtherJetRecPhi->at(ij));
			    BJetE.push_back(OtherJetRecE->at(ij));
			    BJetCSVv2.push_back(OtherJetRecCSVv2->at(ij));

			    NonBJetPt.push_back(OtherJetRecPt->at(ij));
			    NonBJetEta.push_back(OtherJetRecEta->at(ij));
			    NonBJetPhi.push_back(OtherJetRecPhi->at(ij));
			    NonBJetE.push_back(OtherJetRecE->at(ij));
			    NonBJetCSVv2.push_back(OtherJetRecCSVv2->at(ij));
			 }
		    }		  

		  std::vector<std::pair<float,int> > NonBJetPtSort;
		  for(int inb=0;inb<NonBJetPt.size();inb++)
		    {
		       NonBJetPtSort.push_back(std::make_pair(NonBJetPt.at(inb),inb));
		    }		  
		  std::sort(NonBJetPtSort.begin(),NonBJetPtSort.end(),sortFunc());
		  
		  std::vector<int> NonBJetIdx;
		  int nNonBJetMaxEvent = (nNonBJetMax < 0) ? NonBJetPtSort.size() : nNonBJetMax;
		  nNonBJetMaxEvent = (nNonBJetMax < NonBJetPtSort.size() && nNonBJetMax >= 0) ? nNonBJetMax : NonBJetPtSort.size();
		  
		  for(int inb=0;inb<NonBJetPt.size();inb++)
		    {
		       for(int inb2=0;inb2<nNonBJetMaxEvent;inb2++)
			 {
			    if( NonBJetPtSort.at(inb2).second == inb ) NonBJetIdx.push_back(inb);
			 }		       		       
		    }	
		  
		  std::vector<float> NonBJetFilteredPt;
		  std::vector<float> NonBJetFilteredEta;
		  std::vector<float> NonBJetFilteredPhi;
		  std::vector<float> NonBJetFilteredE;
		  
		  for(int inb=0;inb<nNonBJetMaxEvent;inb++)
		    {
		       int idx = NonBJetIdx.at(inb);
		       
		       NonBJetFilteredPt.push_back(NonBJetPt.at(idx));
		       NonBJetFilteredEta.push_back(NonBJetEta.at(idx));
		       NonBJetFilteredPhi.push_back(NonBJetPhi.at(idx));
		       NonBJetFilteredE.push_back(NonBJetE.at(idx));
		    }		  
		  
		  if( it == 0 ) nEventsTruth++;
		  else if( it == 1 ) nEventsHighestCSVv2++;
		  else if( it == 2 ) nEventsCSVv2L++;
		  else if( it == 3 ) nEventsCSVv2M++;
		  else if( it == 4 ) nEventsCSVv2T++;
		  else if( it == 5 ) nEventsAll++;
		  
		  if( BJetPt.size() < 2 || NonBJetPt.size() < 2 ) continue;
		  if( it == 5 && BJetPt.size() < 4 ) continue;
		  
		  if( it == 0 ) {NBJetTruth = BJetPt.size();NNonBJetTruth = NonBJetPt.size();}
		  else if( it == 1 ) {NBJetHighestCSVv2 = BJetPt.size();NNonBJetHighestCSVv2 = NonBJetPt.size();}
		  else if( it == 2 ) {NBJetCSVv2L = BJetPt.size();NNonBJetCSVv2L = NonBJetPt.size();}
		  else if( it == 3 ) {NBJetCSVv2M = BJetPt.size();NNonBJetCSVv2M = NonBJetPt.size();}
		  else if( it == 4 ) {NBJetCSVv2T = BJetPt.size();NNonBJetCSVv2T = NonBJetPt.size();}
		  else if( it == 5 ) {NBJetAll = BJetPt.size();NNonBJetAll = NonBJetPt.size();}

		  if( fabs(TopLepWLepId) == 11 )
		    {	     
		       ElectronPt.push_back(TopLepWLepRecPt);
		       ElectronEta.push_back(TopLepWLepRecEta);
		       ElectronPhi.push_back(TopLepWLepRecPhi);
		       ElectronE.push_back(TopLepWLepRecE);
		    }
		  else
		    {
		       MuonPt.push_back(TopLepWLepRecPt);
		       MuonEta.push_back(TopLepWLepRecEta);
		       MuonPhi.push_back(TopLepWLepRecPhi);
		       MuonE.push_back(TopLepWLepRecE);
		    }	

		  kf->SetBJet(BJetPt,BJetEta,BJetPhi,BJetE);
		  kf->SetNonBJet(NonBJetFilteredPt,NonBJetFilteredEta,NonBJetFilteredPhi,NonBJetFilteredE);	
		  kf->SetElectron(ElectronPt,ElectronEta,ElectronPhi,ElectronE);
		  kf->SetMuon(MuonPt,MuonEta,MuonPhi,MuonE);
		  kf->SetMet(MetRecPx,MetRecPy);
		  
		  if( it == 5 ) kf->NoBTag(1);
		  else kf->NoBTag(0);

		  kf->Run();

		  int NPerm = kf->GetNPerm();
		  
		  if( it == 0 ) NPermTruth = NPerm;
		  else if( it == 1 ) NPermHighestCSVv2 = NPerm;
		  else if( it == 2 ) NPermCSVv2L = NPerm;
		  else if( it == 3 ) NPermCSVv2M = NPerm;
		  else if( it == 4 ) NPermCSVv2T = NPerm;
		  else if( it == 5 ) NPermAll = NPerm;

		  std::vector<std::pair<float,int> > MVATruth;
		  std::vector<std::pair<float,int> > MVAAll;
		  std::vector<std::pair<float,int> > MVAHighestCSVv2;
		  std::vector<std::pair<float,int> > MVACSVv2L;
		  std::vector<std::pair<float,int> > MVACSVv2M;
		  std::vector<std::pair<float,int> > MVACSVv2T;
		  
		  bool foundNoSolutionTruth = 0;
		  bool foundNoSolutionAll = 0;
		  bool foundNoSolutionHighestCSVv2 = 0;
		  bool foundNoSolutionCSVv2L = 0;
		  bool foundNoSolutionCSVv2M = 0;
		  bool foundNoSolutionCSVv2T = 0;
		  
		  for(int ip=0;ip<NPerm;ip++)
		    {
		       float disc = kf->GetDisc(ip);

		       bool MatchTopLepBJet = 0;
		       bool MatchTopHadBJet = 0;
		       bool MatchTopHadWNonBJet1 = 0;
		       bool MatchTopHadWNonBJet2 = 0;

		       int idxTopLepWElecFit = kf->GetIndex(ELECTRON_TOPTOPLEPHAD,ip);
		       int idxTopLepWMuonFit = kf->GetIndex(MUON_TOPTOPLEPHAD,ip);
		       int idxTopLepBJetFit = kf->GetIndex(BJETLEP_TOPTOPLEPHAD,ip);
		       int idxTopHadBJetFit = kf->GetIndex(BJETHAD_TOPTOPLEPHAD,ip);
		       int idxTopHadWNonBJet1Fit = kf->GetIndex(NONBJET1_TOPTOPLEPHAD,ip);
		       int idxTopHadWNonBJet2Fit = kf->GetIndex(NONBJET2_TOPTOPLEPHAD,ip);

		       float NuPx = kf->GetNuPx(ip,0);
		       float NuPy = kf->GetNuPy(ip,0);
		       float NuPz = kf->GetNuPz(ip,0);
		       float NuE = sqrt(NuPx*NuPx+NuPy*NuPy+NuPz*NuPz);

		       TLorentzVector *TopLepWNuFitP4 = new TLorentzVector();
		       TopLepWNuFitP4->SetPxPyPzE(NuPx,NuPy,NuPz,NuE);
		       
		       float TopLepWLepFitPt = (idxTopLepWElecFit >= 0) ? ElectronPt[idxTopLepWElecFit] : MuonPt[idxTopLepWMuonFit];
		       float TopLepWLepFitEta = (idxTopLepWElecFit >= 0) ? ElectronEta[idxTopLepWElecFit] : MuonEta[idxTopLepWMuonFit];
		       float TopLepWLepFitPhi = (idxTopLepWElecFit >= 0) ? ElectronPhi[idxTopLepWElecFit] : MuonPhi[idxTopLepWMuonFit];
		       float TopLepWLepFitE = (idxTopLepWElecFit >= 0) ? ElectronE[idxTopLepWElecFit] : MuonE[idxTopLepWMuonFit];
		       
		       TLorentzVector *TopLepWLepFitP4 = new TLorentzVector();
		       TopLepWLepFitP4->SetPtEtaPhiE(TopLepWLepFitPt,TopLepWLepFitEta,TopLepWLepFitPhi,TopLepWLepFitE);
		       
		       float TopLepBJetFitPt = BJetPt[idxTopLepBJetFit];
		       float TopLepBJetFitEta = BJetEta[idxTopLepBJetFit];
		       float TopLepBJetFitPhi = BJetPhi[idxTopLepBJetFit];
		       float TopLepBJetFitE = BJetE[idxTopLepBJetFit];
		       
		       TLorentzVector *TopLepBJetFitP4 = new TLorentzVector();
		       TopLepBJetFitP4->SetPtEtaPhiE(TopLepBJetFitPt,TopLepBJetFitEta,TopLepBJetFitPhi,TopLepBJetFitE);
		       
		       float TopHadBJetFitPt = BJetPt[idxTopHadBJetFit];
		       float TopHadBJetFitEta = BJetEta[idxTopHadBJetFit];
		       float TopHadBJetFitPhi = BJetPhi[idxTopHadBJetFit];
		       float TopHadBJetFitE = BJetE[idxTopHadBJetFit];
		       
		       TLorentzVector *TopHadBJetFitP4 = new TLorentzVector();
		       TopHadBJetFitP4->SetPtEtaPhiE(TopHadBJetFitPt,TopHadBJetFitEta,TopHadBJetFitPhi,TopHadBJetFitE);
		       
		       float TopHadWNonBJet1FitPt = NonBJetFilteredPt[idxTopHadWNonBJet1Fit];
		       float TopHadWNonBJet1FitEta = NonBJetFilteredEta[idxTopHadWNonBJet1Fit];
		       float TopHadWNonBJet1FitPhi = NonBJetFilteredPhi[idxTopHadWNonBJet1Fit];
		       float TopHadWNonBJet1FitE = NonBJetFilteredE[idxTopHadWNonBJet1Fit];
		       
		       TLorentzVector *TopHadWNonBJet1FitP4 = new TLorentzVector();
		       TopHadWNonBJet1FitP4->SetPtEtaPhiE(TopHadWNonBJet1FitPt,TopHadWNonBJet1FitEta,TopHadWNonBJet1FitPhi,TopHadWNonBJet1FitE);
		       
		       float TopHadWNonBJet2FitPt = NonBJetFilteredPt[idxTopHadWNonBJet2Fit];
		       float TopHadWNonBJet2FitEta = NonBJetFilteredEta[idxTopHadWNonBJet2Fit];
		       float TopHadWNonBJet2FitPhi = NonBJetFilteredPhi[idxTopHadWNonBJet2Fit];
		       float TopHadWNonBJet2FitE = NonBJetFilteredE[idxTopHadWNonBJet2Fit];

		       TLorentzVector *TopHadWNonBJet2FitP4 = new TLorentzVector();
		       TopHadWNonBJet2FitP4->SetPtEtaPhiE(TopHadWNonBJet2FitPt,TopHadWNonBJet2FitEta,TopHadWNonBJet2FitPhi,TopHadWNonBJet2FitE);

		       MatchTopLepBJet = (TopLepBJetFitEta-TopLepBJetRecEta < 10E-6 &&
					  TopLepBJetFitPhi-TopLepBJetRecPhi < 10E-6);
		       
		       MatchTopHadBJet = (TopHadBJetFitEta-TopHadBJetRecEta < 10E-6 &&
					  TopHadBJetFitPhi-TopHadBJetRecPhi < 10E-6);
		       
		       MatchTopHadWNonBJet1 = ((TopHadWNonBJet1FitEta-TopHadWNonBJet1RecEta < 10E-6 &&
						TopHadWNonBJet1FitPhi-TopHadWNonBJet1RecPhi < 10E-6)) ||
			 ((TopHadWNonBJet1FitEta-TopHadWNonBJet2RecEta < 10E-6 &&
			   TopHadWNonBJet1FitPhi-TopHadWNonBJet2RecPhi < 10E-6));
		       
		       MatchTopHadWNonBJet2 = ((TopHadWNonBJet2FitEta-TopHadWNonBJet1RecEta < 10E-6 &&
						TopHadWNonBJet2FitPhi-TopHadWNonBJet1RecPhi < 10E-6)) ||
			 ((TopHadWNonBJet2FitEta-TopHadWNonBJet2RecEta < 10E-6 &&
			   TopHadWNonBJet2FitPhi-TopHadWNonBJet2RecPhi < 10E-6));

		       TLorentzVector TopHadW = *TopHadWNonBJet1FitP4+*TopHadWNonBJet2FitP4;
		       TLorentzVector TopLep = *TopLepWLepFitP4+*TopLepWNuFitP4+*TopLepBJetFitP4;
		       TLorentzVector TopHad = TopHadW+*TopHadBJetFitP4;

		       float VarTopHadWRecM = TopHadW.M();
		       float VarTopLepRecM = TopLep.M();
		       float VarTopHadRecM = TopHad.M();
		       float VarTopLepTopHadRecDr = TopLep.DeltaR(TopHad);
		       float VarTopLepRecPt = TopLep.Pt();
		       float VarTopHadRecPt = TopHad.Pt();

		       TLorentzVector *TopHadFitT = new TLorentzVector();
		       TopHadFitT->SetPxPyPzE(TopHad.Px(),TopHad.Py(),0.,TopHad.Et());
		       
		       TLorentzVector *TopLepWLepFitT = new TLorentzVector();
		       TopLepWLepFitT->SetPxPyPzE(TopLepWLepFitP4->Px(),TopLepWLepFitP4->Py(),0.,TopLepWLepFitP4->Et());

		       TLorentzVector *TopLepWNuFitT = new TLorentzVector();
		       TopLepWNuFitT->SetPxPyPzE(TopLepWNuFitP4->Px(),TopLepWNuFitP4->Py(),0.,TopLepWNuFitP4->Et());

		       TLorentzVector *TopLepBJetFitT = new TLorentzVector();
		       TopLepBJetFitT->SetPxPyPzE(TopLepBJetFitP4->Px(),TopLepBJetFitP4->Py(),0.,TopLepBJetFitP4->Et());
		       
		       TLorentzVector TopLepT = *TopLepWLepFitT+*TopLepWNuFitT+*TopLepBJetFitT;

		       float VarTopLepRecMT = TopLepT.M();		       
		       float VarTopLepTopHadRecDphiT = TopLepT.DeltaPhi(*TopHadFitT);
		       float VarTopLepRecPtT = TopLepT.Pt();
		       
		       delete TopHadFitT;
		       delete TopLepWLepFitT;
		       delete TopLepWNuFitT;
		       delete TopLepBJetFitT;
		       
		       delete TopLepWLepFitP4;
		       delete TopLepWNuFitP4;
		       delete TopLepBJetFitP4;
		       delete TopHadBJetFitP4;
		       delete TopHadWNonBJet1FitP4;
		       delete TopHadWNonBJet2FitP4;
		       
		       bool hasMatch = MatchTopLepBJet && MatchTopHadBJet && MatchTopHadWNonBJet1 && MatchTopHadWNonBJet2;
		       
		       bool hasMatchBJet = MatchTopLepBJet && MatchTopHadBJet;
		       
		       if( it == 0 ) 
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionTruth ) {nNoSolutionTruth++; foundNoSolutionTruth=1;}
			    if( ip == 0 ) nSelTruth++;
			    MatchTruth.push_back(hasMatch);
			    MatchBJetTruth.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchTruth++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetTruth++;
			    DiscTruth.push_back(disc);
			    TopHadWRecMTruth.push_back(VarTopHadWRecM);
			    TopLepRecMTruth.push_back(VarTopLepRecM);
			    TopHadRecMTruth.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrTruth.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtTruth.push_back(VarTopLepRecPt);
			    TopHadRecPtTruth.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTTruth.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTTruth.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTTruth.push_back(VarTopLepRecPtT);
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_TopHadRecMTruth = VarTopHadRecM;
				      MVAFullReco_TopLepRecMTruth = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrTruth = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtTruth = VarTopLepRecPt;
				      
				      MVATruth.push_back(std::make_pair(readerFullRecoTruth->EvaluateMVA("BDTG method"),ip));
				   }				 
				 else
				   {				      
				      MVAPartReco_TopHadRecMTruth = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTTruth = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTTruth = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTTruth = VarTopLepRecPtT;
				      
				      MVATruth.push_back(std::make_pair(readerPartRecoTruth->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }			    
			 }	     
		       else if( it == 1 )
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionHighestCSVv2 ) {nNoSolutionHighestCSVv2++; foundNoSolutionHighestCSVv2 = 1;}
			    if( ip == 0 ) nSelHighestCSVv2++;
			    MatchHighestCSVv2.push_back(hasMatch);
			    MatchBJetHighestCSVv2.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchHighestCSVv2++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetHighestCSVv2++;
			    DiscHighestCSVv2.push_back(disc);
			    TopHadWRecMHighestCSVv2.push_back(VarTopHadWRecM);
			    TopLepRecMHighestCSVv2.push_back(VarTopLepRecM);
			    TopHadRecMHighestCSVv2.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrHighestCSVv2.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtHighestCSVv2.push_back(VarTopLepRecPt);
			    TopHadRecPtHighestCSVv2.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTHighestCSVv2.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTHighestCSVv2.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTHighestCSVv2.push_back(VarTopLepRecPtT);

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_TopHadRecMHighestCSVv2 = VarTopHadRecM;
				      MVAFullReco_TopLepRecMHighestCSVv2 = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrHighestCSVv2 = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtHighestCSVv2 = VarTopLepRecPt;
				      
				      MVAHighestCSVv2.push_back(std::make_pair(readerFullRecoHighestCSVv2->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_TopHadRecMHighestCSVv2 = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTHighestCSVv2 = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTHighestCSVv2 = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTHighestCSVv2 = VarTopLepRecPtT;
				      
				      MVAHighestCSVv2.push_back(std::make_pair(readerPartRecoHighestCSVv2->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }
			 }
		       else if( it == 2 )
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionCSVv2L ) {nNoSolutionCSVv2L++; foundNoSolutionCSVv2L = 1;}
			    if( ip == 0 ) nSelCSVv2L++;
			    MatchCSVv2L.push_back(hasMatch);
			    MatchBJetCSVv2L.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchCSVv2L++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetCSVv2L++;
			    DiscCSVv2L.push_back(disc);
			    TopHadWRecMCSVv2L.push_back(VarTopHadWRecM);
			    TopLepRecMCSVv2L.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2L.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrCSVv2L.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtCSVv2L.push_back(VarTopLepRecPt);
			    TopHadRecPtCSVv2L.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTCSVv2L.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTCSVv2L.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTCSVv2L.push_back(VarTopLepRecPtT);

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_TopHadRecMCSVv2L = VarTopHadRecM;
				      MVAFullReco_TopLepRecMCSVv2L = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrCSVv2L = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtCSVv2L = VarTopLepRecPt;
				 
				      MVACSVv2L.push_back(std::make_pair(readerFullRecoCSVv2L->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_TopHadRecMCSVv2L = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTCSVv2L = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTCSVv2L = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2L = VarTopLepRecPtT;
				      
				      MVACSVv2L.push_back(std::make_pair(readerPartRecoCSVv2L->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }			    
			 }	     
		       else if( it == 3 )			 
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionCSVv2M ) {nNoSolutionCSVv2M++; foundNoSolutionCSVv2M = 1;}
			    if( ip == 0 ) nSelCSVv2M++;
			    MatchCSVv2M.push_back(hasMatch);
			    MatchBJetCSVv2M.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchCSVv2M++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetCSVv2M++;
			    DiscCSVv2M.push_back(disc);
			    TopHadWRecMCSVv2M.push_back(VarTopHadWRecM);
			    TopLepRecMCSVv2M.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2M.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrCSVv2M.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtCSVv2M.push_back(VarTopLepRecPt);
			    TopHadRecPtCSVv2M.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTCSVv2M.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTCSVv2M.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTCSVv2M.push_back(VarTopLepRecPtT);

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {				 
				      MVAFullReco_TopHadRecMCSVv2M = VarTopHadRecM;
				      MVAFullReco_TopLepRecMCSVv2M = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrCSVv2M = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtCSVv2M = VarTopLepRecPt;
				      
				      MVACSVv2M.push_back(std::make_pair(readerFullRecoCSVv2M->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {
				      MVAPartReco_TopHadRecMCSVv2M = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTCSVv2M = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTCSVv2M = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2M = VarTopLepRecPtT;
				      
				      MVACSVv2M.push_back(std::make_pair(readerPartRecoCSVv2M->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }			    
			 }	     
		       else if( it == 4 )
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionCSVv2T ) {nNoSolutionCSVv2T++; foundNoSolutionCSVv2T = 1;}
			    if( ip == 0 ) nSelCSVv2T++;
			    MatchCSVv2T.push_back(hasMatch);
			    MatchBJetCSVv2T.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchCSVv2T++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetCSVv2T++;
			    DiscCSVv2T.push_back(disc);
			    TopHadWRecMCSVv2T.push_back(VarTopHadWRecM);
			    TopLepRecMCSVv2T.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2T.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrCSVv2T.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtCSVv2T.push_back(VarTopLepRecPt);
			    TopHadRecPtCSVv2T.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTCSVv2T.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTCSVv2T.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTCSVv2T.push_back(VarTopLepRecPtT);

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_TopHadRecMCSVv2T = VarTopHadRecM;
				      MVAFullReco_TopLepRecMCSVv2T = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrCSVv2T = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtCSVv2T = VarTopLepRecPt;
				 
				      MVACSVv2T.push_back(std::make_pair(readerFullRecoCSVv2T->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_TopHadRecMCSVv2T = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTCSVv2T = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTCSVv2T = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2T = VarTopLepRecPtT;
				      
				      MVACSVv2T.push_back(std::make_pair(readerPartRecoCSVv2T->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }			    
			 }		       
		       else if( it == 5 )
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionAll ) {nNoSolutionAll++; foundNoSolutionAll = 1;}
			    if( ip == 0 ) nSelAll++;
			    MatchAll.push_back(hasMatch);
			    MatchBJetAll.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchAll++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetAll++;
			    DiscAll.push_back(disc);
			    TopHadWRecMAll.push_back(VarTopHadWRecM);
			    TopLepRecMAll.push_back(VarTopLepRecM);
			    TopHadRecMAll.push_back(VarTopHadRecM);
			    TopLepTopHadRecDrAll.push_back(VarTopLepTopHadRecDr);
			    TopLepRecPtAll.push_back(VarTopLepRecPt);
			    TopHadRecPtAll.push_back(VarTopHadRecPt);
			    
			    TopLepRecMTAll.push_back(VarTopLepRecMT);
			    TopLepTopHadRecDphiTAll.push_back(VarTopLepTopHadRecDphiT);
			    TopLepRecPtTAll.push_back(VarTopLepRecPtT);

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_TopHadRecMAll = VarTopHadRecM;
				      MVAFullReco_TopLepRecMAll = VarTopLepRecM;
				      MVAFullReco_TopLepTopHadRecDrAll = VarTopLepTopHadRecDr;
				      MVAFullReco_TopLepRecPtAll = VarTopLepRecPt;
				 
				      MVAAll.push_back(std::make_pair(readerFullRecoAll->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_TopHadRecMAll = VarTopHadRecM;
				      MVAPartReco_TopLepRecMTAll = VarTopLepRecMT;
				      MVAPartReco_TopLepTopHadRecDphiTAll = VarTopLepTopHadRecDphiT;
				      MVAPartReco_TopLepRecPtTAll = VarTopLepRecPtT;
				      
				      MVAAll.push_back(std::make_pair(readerPartRecoAll->EvaluateMVA("BDTG method"),ip));
				   }				 
			      }			    
			 }		       
		    }		  

		  if( applyMVA )
		    {
		       if( it == 0 )
			 {			    
			    std::sort(MVATruth.begin(),MVATruth.end(),sortFunc());
			    int MVARes = MVATruth.at(0).second;

			    for(int ip=0;ip<MVATruth.size();ip++)
			      {				 
				 MVAScoreTruth.push_back(MVATruth.at(ip).first);
				 MVADiscTruth.push_back(DiscTruth.at(MVATruth.at(ip).second));
				 MVATopHadRecMTruth.push_back(TopHadRecMTruth.at(MVATruth.at(ip).second));
				 MVATopLepRecMTruth.push_back(TopLepRecMTruth.at(MVATruth.at(ip).second));
			      }			    
			    
			    if( MatchTruth.at(MVARes) ) nMatchMVATruth++;
			    if( MatchBJetTruth.at(MVARes) ) nMatchBJetMVATruth++;
			 }
		       else if( it == 1 )
			 {			    
			    std::sort(MVAHighestCSVv2.begin(),MVAHighestCSVv2.end(),sortFunc());
			    int MVARes = MVAHighestCSVv2.at(0).second;

			    for(int ip=0;ip<MVAHighestCSVv2.size();ip++)
			      {				 
				 MVAScoreHighestCSVv2.push_back(MVAHighestCSVv2.at(ip).first);
				 MVADiscHighestCSVv2.push_back(DiscHighestCSVv2.at(MVAHighestCSVv2.at(ip).second));
				 MVATopHadRecMHighestCSVv2.push_back(TopHadRecMHighestCSVv2.at(MVAHighestCSVv2.at(ip).second));
				 MVATopLepRecMHighestCSVv2.push_back(TopLepRecMHighestCSVv2.at(MVAHighestCSVv2.at(ip).second));
			      }			    
			    
			    if( MatchHighestCSVv2.at(MVARes) ) nMatchMVAHighestCSVv2++;
			    if( MatchBJetHighestCSVv2.at(MVARes) ) nMatchBJetMVAHighestCSVv2++;
			 }
		       else if( it == 2 )
			 {			    
			    std::sort(MVACSVv2L.begin(),MVACSVv2L.end(),sortFunc());
			    int MVARes = MVACSVv2L.at(0).second;

			    for(int ip=0;ip<MVACSVv2L.size();ip++)
			      {				 
				 MVAScoreCSVv2L.push_back(MVACSVv2L.at(ip).first);
				 MVADiscCSVv2L.push_back(DiscCSVv2L.at(MVACSVv2L.at(ip).second));
				 MVATopHadRecMCSVv2L.push_back(TopHadRecMCSVv2L.at(MVACSVv2L.at(ip).second));
				 MVATopLepRecMCSVv2L.push_back(TopLepRecMCSVv2L.at(MVACSVv2L.at(ip).second));
			      }			    
			    
			    if( MatchCSVv2L.at(MVARes) ) nMatchMVACSVv2L++;
			    if( MatchBJetCSVv2L.at(MVARes) ) nMatchBJetMVACSVv2L++;
			 }
		       else if( it == 3 )
			 {			    
			    std::sort(MVACSVv2M.begin(),MVACSVv2M.end(),sortFunc());
			    int MVARes = MVACSVv2M.at(0).second;

			    for(int ip=0;ip<MVACSVv2M.size();ip++)
			      {
				 MVAScoreCSVv2M.push_back(MVACSVv2M.at(ip).first);
				 MVADiscCSVv2M.push_back(DiscCSVv2M.at(MVACSVv2M.at(ip).second));
				 MVATopHadRecMCSVv2M.push_back(TopHadRecMCSVv2M.at(MVACSVv2M.at(ip).second));
				 MVATopLepRecMCSVv2M.push_back(TopLepRecMCSVv2M.at(MVACSVv2M.at(ip).second));
			      }			    
			    
			    if( MatchCSVv2M.at(MVARes) ) nMatchMVACSVv2M++;
			    if( MatchBJetCSVv2M.at(MVARes) ) nMatchBJetMVACSVv2M++;
			 }
		       else if( it == 4 )
			 {			    
			    std::sort(MVACSVv2T.begin(),MVACSVv2T.end(),sortFunc());
			    int MVARes = MVACSVv2T.at(0).second;

			    for(int ip=0;ip<MVACSVv2T.size();ip++)
			      {
				 MVAScoreCSVv2T.push_back(MVACSVv2T.at(ip).first);
				 MVADiscCSVv2T.push_back(DiscCSVv2T.at(MVACSVv2T.at(ip).second));
				 MVATopHadRecMCSVv2T.push_back(TopHadRecMCSVv2T.at(MVACSVv2T.at(ip).second));
				 MVATopLepRecMCSVv2T.push_back(TopLepRecMCSVv2T.at(MVACSVv2T.at(ip).second));
			      }			    
			    
			    if( MatchCSVv2T.at(MVARes) ) nMatchMVACSVv2T++;
			    if( MatchBJetCSVv2T.at(MVARes) ) nMatchBJetMVACSVv2T++;
			 }
		       else if( it == 5 )
			 {			    
			    std::sort(MVAAll.begin(),MVAAll.end(),sortFunc());
			    int MVARes = MVAAll.at(0).second;
			    
			    for(int ip=0;ip<MVAAll.size();ip++)
			      {
				 MVAScoreAll.push_back(MVAAll.at(ip).first);
				 MVADiscAll.push_back(DiscAll.at(MVAAll.at(ip).second));
				 MVATopHadRecMAll.push_back(TopHadRecMAll.at(MVAAll.at(ip).second));
				 MVATopLepRecMAll.push_back(TopLepRecMAll.at(MVAAll.at(ip).second));
			      }			    
			    
			    if( MatchAll.at(MVARes) ) nMatchMVAAll++;
			    if( MatchBJetAll.at(MVARes) ) nMatchBJetMVAAll++;
			 }
		    }
	       }	     
	     
	     trFIT->Fill();	     
	  }
	
	float SelEffTruthErr = errf(nSelTruth,sqrt(nSelTruth),nEventsTruth,sqrt(nEventsTruth));
	float SelEffAllErr = errf(nSelAll,sqrt(nSelAll),nEventsAll,sqrt(nEventsAll));
	float SelEffHighestCSVv2Err = errf(nSelHighestCSVv2,sqrt(nSelHighestCSVv2),nEventsHighestCSVv2,sqrt(nEventsHighestCSVv2));
	float SelEffCSVv2LErr = errf(nSelCSVv2L,sqrt(nSelCSVv2L),nEventsCSVv2L,sqrt(nEventsCSVv2L));
	float SelEffCSVv2MErr = errf(nSelCSVv2M,sqrt(nSelCSVv2M),nEventsCSVv2M,sqrt(nEventsCSVv2M));
	float SelEffCSVv2TErr = errf(nSelCSVv2T,sqrt(nSelCSVv2T),nEventsCSVv2T,sqrt(nEventsCSVv2T));
	
	float AlgEffTruthErr = errf(nMatchTruth,sqrt(nMatchTruth),nSelTruth,sqrt(nSelTruth));
	float AlgEffAllErr = errf(nMatchAll,sqrt(nMatchAll),nSelAll,sqrt(nSelAll));
	float AlgEffHighestCSVv2Err = errf(nMatchHighestCSVv2,sqrt(nMatchHighestCSVv2),nSelHighestCSVv2,sqrt(nSelHighestCSVv2));
	float AlgEffCSVv2LErr = errf(nMatchCSVv2L,sqrt(nMatchCSVv2L),nSelCSVv2L,sqrt(nSelCSVv2L));
	float AlgEffCSVv2MErr = errf(nMatchCSVv2M,sqrt(nMatchCSVv2M),nSelCSVv2M,sqrt(nSelCSVv2M));
	float AlgEffCSVv2TErr = errf(nMatchCSVv2T,sqrt(nMatchCSVv2T),nSelCSVv2T,sqrt(nSelCSVv2T));

	float AlgBJetEffTruthErr = errf(nMatchBJetTruth,sqrt(nMatchBJetTruth),nSelTruth,sqrt(nSelTruth));
	float AlgBJetEffAllErr = errf(nMatchBJetAll,sqrt(nMatchBJetAll),nSelAll,sqrt(nSelAll));
	float AlgBJetEffHighestCSVv2Err = errf(nMatchBJetHighestCSVv2,sqrt(nMatchBJetHighestCSVv2),nSelHighestCSVv2,sqrt(nSelHighestCSVv2));
	float AlgBJetEffCSVv2LErr = errf(nMatchBJetCSVv2L,sqrt(nMatchBJetCSVv2L),nSelCSVv2L,sqrt(nSelCSVv2L));
	float AlgBJetEffCSVv2MErr = errf(nMatchBJetCSVv2M,sqrt(nMatchBJetCSVv2M),nSelCSVv2M,sqrt(nSelCSVv2M));
	float AlgBJetEffCSVv2TErr = errf(nMatchBJetCSVv2T,sqrt(nMatchBJetCSVv2T),nSelCSVv2T,sqrt(nSelCSVv2T));
	
	float NoSolutionEffTruthErr = errf(nNoSolutionTruth,sqrt(nNoSolutionTruth),nSelTruth,sqrt(nSelTruth));
	float NoSolutionEffAllErr = errf(nNoSolutionAll,sqrt(nNoSolutionAll),nSelAll,sqrt(nSelAll));
	float NoSolutionEffHighestCSVv2Err = errf(nNoSolutionHighestCSVv2,sqrt(nNoSolutionHighestCSVv2),nSelHighestCSVv2,sqrt(nSelHighestCSVv2));
	float NoSolutionEffCSVv2LErr = errf(nNoSolutionCSVv2L,sqrt(nNoSolutionCSVv2L),nSelCSVv2L,sqrt(nSelCSVv2L));
	float NoSolutionEffCSVv2MErr = errf(nNoSolutionCSVv2M,sqrt(nNoSolutionCSVv2M),nSelCSVv2M,sqrt(nSelCSVv2M));
	float NoSolutionEffCSVv2TErr = errf(nNoSolutionCSVv2T,sqrt(nNoSolutionCSVv2T),nSelCSVv2T,sqrt(nSelCSVv2T));
	
	float MVAEffTruthErr = errf(nMatchMVATruth,sqrt(nMatchMVATruth),nSelTruth,sqrt(nSelTruth));
	float MVAEffAllErr = errf(nMatchMVAAll,sqrt(nMatchMVAAll),nSelAll,sqrt(nSelAll));
	float MVAEffHighestCSVv2Err = errf(nMatchMVAHighestCSVv2,sqrt(nMatchMVAHighestCSVv2),nSelHighestCSVv2,sqrt(nSelHighestCSVv2));
	float MVAEffCSVv2LErr = errf(nMatchMVACSVv2L,sqrt(nMatchMVACSVv2L),nSelCSVv2L,sqrt(nSelCSVv2L));
	float MVAEffCSVv2MErr = errf(nMatchMVACSVv2M,sqrt(nMatchMVACSVv2M),nSelCSVv2M,sqrt(nSelCSVv2M));
	float MVAEffCSVv2TErr = errf(nMatchMVACSVv2T,sqrt(nMatchMVACSVv2T),nSelCSVv2T,sqrt(nSelCSVv2T));

	float MVABJetEffTruthErr = errf(nMatchBJetMVATruth,sqrt(nMatchBJetMVATruth),nSelTruth,sqrt(nSelTruth));
	float MVABJetEffAllErr = errf(nMatchBJetMVAAll,sqrt(nMatchBJetMVAAll),nSelAll,sqrt(nSelAll));
	float MVABJetEffHighestCSVv2Err = errf(nMatchBJetMVAHighestCSVv2,sqrt(nMatchBJetMVAHighestCSVv2),nSelHighestCSVv2,sqrt(nSelHighestCSVv2));
	float MVABJetEffCSVv2LErr = errf(nMatchBJetMVACSVv2L,sqrt(nMatchBJetMVACSVv2L),nSelCSVv2L,sqrt(nSelCSVv2L));
	float MVABJetEffCSVv2MErr = errf(nMatchBJetMVACSVv2M,sqrt(nMatchBJetMVACSVv2M),nSelCSVv2M,sqrt(nSelCSVv2M));
	float MVABJetEffCSVv2TErr = errf(nMatchBJetMVACSVv2T,sqrt(nMatchBJetMVACSVv2T),nSelCSVv2T,sqrt(nSelCSVv2T));
	
	std::cout << "Selection Efficiency (Truth) = " << float(nSelTruth)/float(nEventsTruth)*100 << " +- " << SelEffTruthErr*100 << " %" << std::endl;
	std::cout << "Selection Efficiency (All) = " << float(nSelAll)/float(nEventsAll)*100 << " +- " << SelEffAllErr*100 << " %" << std::endl;
	std::cout << "Selection Efficiency (HighestCSVv2) = " << float(nSelHighestCSVv2)/float(nEventsHighestCSVv2)*100 << " +- " << SelEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "Selection Efficiency (CSVv2L) = " << float(nSelCSVv2L)/float(nEventsCSVv2L)*100 << " +- " << SelEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "Selection Efficiency (CSVv2M) = " << float(nSelCSVv2M)/float(nEventsCSVv2M)*100 << " +- " << SelEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "Selection Efficiency (CSVv2T) = " << float(nSelCSVv2T)/float(nEventsCSVv2T)*100 << " +- " << SelEffCSVv2TErr*100 << " %" << std::endl;
	
	std::cout << "Algorithm Efficiency for all jets (Truth) = " << float(nMatchTruth)/float(nSelTruth)*100 << " +- " << AlgEffTruthErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for all jets (All) = " << float(nMatchAll)/float(nSelAll)*100 << " +- " << AlgEffAllErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for all jets (HighestCSVv2) = " << float(nMatchHighestCSVv2)/float(nSelHighestCSVv2)*100 << " +- " << AlgEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for all jets (CSVv2L) = " << float(nMatchCSVv2L)/float(nSelCSVv2L)*100 << " +- " << AlgEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for all jets (CSVv2M) = " << float(nMatchCSVv2M)/float(nSelCSVv2M)*100 << " +- " << AlgEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for all jets (CSVv2T) = " << float(nMatchCSVv2T)/float(nSelCSVv2T)*100 << " +- " << AlgEffCSVv2TErr*100 << " %" << std::endl;
       
	std::cout << "Algorithm+MVA Efficiency for all jets (Truth) = " << float(nMatchMVATruth)/float(nSelTruth)*100 << " +- " << MVAEffTruthErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for all jets (All) = " << float(nMatchMVAAll)/float(nSelAll)*100 << " +- " << MVAEffAllErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for all jets (HighestCSVv2) = " << float(nMatchMVAHighestCSVv2)/float(nSelHighestCSVv2)*100 << " +- " << MVAEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2L) = " << float(nMatchMVACSVv2L)/float(nSelCSVv2L)*100 << " +- " << MVAEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2M) = " << float(nMatchMVACSVv2M)/float(nSelCSVv2M)*100 << " +- " << MVAEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2T) = " << float(nMatchMVACSVv2T)/float(nSelCSVv2T)*100 << " +- " << MVAEffCSVv2TErr*100 << " %" << std::endl;

	std::cout << "Algorithm Efficiency for b jets (Truth) = " << float(nMatchBJetTruth)/float(nSelTruth)*100 << " +- " << AlgBJetEffTruthErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for b jets (All) = " << float(nMatchBJetAll)/float(nSelAll)*100 << " +- " << AlgBJetEffAllErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for b jets (HighestCSVv2) = " << float(nMatchBJetHighestCSVv2)/float(nSelHighestCSVv2)*100 << " +- " << AlgBJetEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for b jets (CSVv2L) = " << float(nMatchBJetCSVv2L)/float(nSelCSVv2L)*100 << " +- " << AlgBJetEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for b jets (CSVv2M) = " << float(nMatchBJetCSVv2M)/float(nSelCSVv2M)*100 << " +- " << AlgBJetEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "Algorithm Efficiency for b jets (CSVv2T) = " << float(nMatchBJetCSVv2T)/float(nSelCSVv2T)*100 << " +- " << AlgBJetEffCSVv2TErr*100 << " %" << std::endl;
       
	std::cout << "Algorithm+MVA Efficiency for b jets (Truth) = " << float(nMatchBJetMVATruth)/float(nSelTruth)*100 << " +- " << MVABJetEffTruthErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for b jets (All) = " << float(nMatchBJetMVAAll)/float(nSelAll)*100 << " +- " << MVABJetEffAllErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for b jets (HighestCSVv2) = " << float(nMatchBJetMVAHighestCSVv2)/float(nSelHighestCSVv2)*100 << " +- " << MVABJetEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2L) = " << float(nMatchBJetMVACSVv2L)/float(nSelCSVv2L)*100 << " +- " << MVABJetEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2M) = " << float(nMatchBJetMVACSVv2M)/float(nSelCSVv2M)*100 << " +- " << MVABJetEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2T) = " << float(nMatchBJetMVACSVv2T)/float(nSelCSVv2T)*100 << " +- " << MVABJetEffCSVv2TErr*100 << " %" << std::endl;
	
	std::cout << "No solutions found (Truth) = " << float(nNoSolutionTruth)/float(nSelTruth)*100 << " +- " << NoSolutionEffTruthErr*100 << " %" << std::endl;
	std::cout << "No solutions found (All) = " << float(nNoSolutionAll)/float(nSelAll)*100 << " +- " << NoSolutionEffAllErr*100 << " %" << std::endl;
	std::cout << "No solutions found (HighestCSVv2) = " << float(nNoSolutionHighestCSVv2)/float(nSelHighestCSVv2)*100 << " +- " << NoSolutionEffHighestCSVv2Err*100 << " %" << std::endl;
	std::cout << "No solutions found (CSVv2L) = " << float(nNoSolutionCSVv2L)/float(nSelCSVv2L)*100 << " +- " << NoSolutionEffCSVv2LErr*100 << " %" << std::endl;
	std::cout << "No solutions found (CSVv2M) = " << float(nNoSolutionCSVv2M)/float(nSelCSVv2M)*100 << " +- " << NoSolutionEffCSVv2MErr*100 << " %" << std::endl;
	std::cout << "No solutions found (CSVv2T) = " << float(nNoSolutionCSVv2T)/float(nSelCSVv2T)*100 << " +- " << NoSolutionEffCSVv2TErr*100 << " %" << std::endl;

	trEFF->Fill();
	
	fout->Write();
	fout->Close();
     }   
   
   f->Close();
   
   delete kf;
}

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{   
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

float errf(float v1,float ve1,float v2,float ve2)
{
   if( v2 == 0 ) return -666;

   float err = ve1*ve1/v2/v2 + v1*v1*ve2*ve2/v2/v2/v2/v2;

   err = sqrt(err);

   return err;
}
 
