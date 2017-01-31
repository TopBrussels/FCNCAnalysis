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
   float MVAFullReco_HiggsRecMTruth;
   float MVAFullReco_TopLepRecMTruth;
   float MVAFullReco_HiggsTopLepRecDrTruth;
   float MVAFullReco_TopLepRecPtTruth;
//   float MVAFullReco_TopLepBJetRecCSVv2Truth;
//   float MVAFullReco_HiggsBJet1RecCSVv2Truth;
//   float MVAFullReco_HiggsBJet2RecCSVv2Truth;
//   float MVAFullReco_TopHadNonBJetRecCSVv2Truth;

   float MVAFullReco_HiggsRecMAll;
   float MVAFullReco_TopLepRecMAll;
   float MVAFullReco_HiggsTopLepRecDrAll;
   float MVAFullReco_TopLepRecPtAll;
//   float MVAFullReco_TopLepBJetRecCSVv2All;
//   float MVAFullReco_HiggsBJet1RecCSVv2All;
//   float MVAFullReco_HiggsBJet2RecCSVv2All;
//   float MVAFullReco_TopHadNonBJetRecCSVv2All;
   
   float MVAFullReco_HiggsRecMHighestCSVv2;
   float MVAFullReco_TopLepRecMHighestCSVv2;
   float MVAFullReco_HiggsTopLepRecDrHighestCSVv2;
   float MVAFullReco_TopLepRecPtHighestCSVv2;
//   float MVAFullReco_TopLepBJetRecCSVv2HighestCSVv2;
//   float MVAFullReco_HiggsBJet1RecCSVv2HighestCSVv2;
//   float MVAFullReco_HiggsBJet2RecCSVv2HighestCSVv2;
//   float MVAFullReco_TopHadNonBJetRecCSVv2HighestCSVv2;

   float MVAFullReco_HiggsRecMCSVv2L;
   float MVAFullReco_TopLepRecMCSVv2L;
   float MVAFullReco_HiggsTopLepRecDrCSVv2L;
   float MVAFullReco_TopLepRecPtCSVv2L;
//   float MVAFullReco_TopLepBJetRecCSVv2CSVv2L;
//   float MVAFullReco_HiggsBJet1RecCSVv2CSVv2L;
//   float MVAFullReco_HiggsBJet2RecCSVv2CSVv2L;
//   float MVAFullReco_TopHadNonBJetRecCSVv2CSVv2L;

   float MVAFullReco_HiggsRecMCSVv2M;
   float MVAFullReco_TopLepRecMCSVv2M;
   float MVAFullReco_HiggsTopLepRecDrCSVv2M;
   float MVAFullReco_TopLepRecPtCSVv2M;
//   float MVAFullReco_TopLepBJetRecCSVv2CSVv2M;
//   float MVAFullReco_HiggsBJet1RecCSVv2CSVv2M;
//   float MVAFullReco_HiggsBJet2RecCSVv2CSVv2M;
//   float MVAFullReco_TopHadNonBJetRecCSVv2CSVv2M;

   float MVAFullReco_HiggsRecMCSVv2T;
   float MVAFullReco_TopLepRecMCSVv2T;
   float MVAFullReco_HiggsTopLepRecDrCSVv2T;
   float MVAFullReco_TopLepRecPtCSVv2T;
//   float MVAFullReco_TopLepBJetRecCSVv2CSVv2T;
//   float MVAFullReco_HiggsBJet1RecCSVv2CSVv2T;
//   float MVAFullReco_HiggsBJet2RecCSVv2CSVv2T;
//   float MVAFullReco_TopHadNonBJetRecCSVv2CSVv2T;

   // PartReco
   float MVAPartReco_HiggsRecMTruth;
   float MVAPartReco_TopLepRecMTTruth;
   float MVAPartReco_HiggsTopLepRecDphiTTruth;
   float MVAPartReco_TopLepRecPtTTruth;
//   float MVAPartReco_TopLepBJetRecCSVv2Truth;
//   float MVAPartReco_HiggsBJet1RecCSVv2Truth;
//   float MVAPartReco_HiggsBJet2RecCSVv2Truth;
//   float MVAPartReco_TopHadNonBJetRecCSVv2Truth;

   float MVAPartReco_HiggsRecMAll;
   float MVAPartReco_TopLepRecMTAll;
   float MVAPartReco_HiggsTopLepRecDphiTAll;
   float MVAPartReco_TopLepRecPtTAll;
//   float MVAPartReco_TopLepBJetRecCSVv2All;
//   float MVAPartReco_HiggsBJet1RecCSVv2All;
//   float MVAPartReco_HiggsBJet2RecCSVv2All;
//   float MVAPartReco_TopHadNonBJetRecCSVv2All;
   
   float MVAPartReco_HiggsRecMHighestCSVv2;
   float MVAPartReco_TopLepRecMTHighestCSVv2;
   float MVAPartReco_HiggsTopLepRecDphiTHighestCSVv2;
   float MVAPartReco_TopLepRecPtTHighestCSVv2;
//   float MVAPartReco_TopLepBJetRecCSVv2HighestCSVv2;
//   float MVAPartReco_HiggsBJet1RecCSVv2HighestCSVv2;
//   float MVAPartReco_HiggsBJet2RecCSVv2HighestCSVv2;
//   float MVAPartReco_TopHadNonBJetRecCSVv2HighestCSVv2;

   float MVAPartReco_HiggsRecMCSVv2L;
   float MVAPartReco_TopLepRecMTCSVv2L;
   float MVAPartReco_HiggsTopLepRecDphiTCSVv2L;
   float MVAPartReco_TopLepRecPtTCSVv2L;
//   float MVAPartReco_TopLepBJetRecCSVv2CSVv2L;
//   float MVAPartReco_HiggsBJet1RecCSVv2CSVv2L;
//   float MVAPartReco_HiggsBJet2RecCSVv2CSVv2L;
//   float MVAPartReco_TopHadNonBJetRecCSVv2CSVv2L;

   float MVAPartReco_HiggsRecMCSVv2M;
   float MVAPartReco_TopLepRecMTCSVv2M;
   float MVAPartReco_HiggsTopLepRecDphiTCSVv2M;
   float MVAPartReco_TopLepRecPtTCSVv2M;
//   float MVAPartReco_TopLepBJetRecCSVv2CSVv2M;
//   float MVAPartReco_HiggsBJet1RecCSVv2CSVv2M;
//   float MVAPartReco_HiggsBJet2RecCSVv2CSVv2M;
//   float MVAPartReco_TopHadNonBJetRecCSVv2CSVv2M;

   float MVAPartReco_HiggsRecMCSVv2T;
   float MVAPartReco_TopLepRecMTCSVv2T;
   float MVAPartReco_HiggsTopLepRecDphiTCSVv2T;
   float MVAPartReco_TopLepRecPtTCSVv2T;
//   float MVAPartReco_TopLepBJetRecCSVv2CSVv2T;
//   float MVAPartReco_HiggsBJet1RecCSVv2CSVv2T;
//   float MVAPartReco_HiggsBJet2RecCSVv2CSVv2T;
//   float MVAPartReco_TopHadNonBJetRecCSVv2CSVv2T;
   
   if( applyMVA )
     {	
	readerFullRecoTruth->AddVariable("HiggsRecMTruth",&MVAFullReco_HiggsRecMTruth);
	readerFullRecoTruth->AddVariable("TopLepRecMTruth",&MVAFullReco_TopLepRecMTruth);
	readerFullRecoTruth->AddVariable("HiggsTopLepRecDrTruth",&MVAFullReco_HiggsTopLepRecDrTruth);
	readerFullRecoTruth->AddVariable("TopLepRecPtTruth",&MVAFullReco_TopLepRecPtTruth);
//	readerFullRecoTruth->AddVariable("TopLepBJetRecCSVv2Truth",&MVAFullReco_TopLepBJetRecCSVv2Truth);
//	readerFullRecoTruth->AddVariable("HiggsBJet1RecCSVv2Truth",&MVAFullReco_HiggsBJet1RecCSVv2Truth);
//	readerFullRecoTruth->AddVariable("HiggsBJet2RecCSVv2Truth",&MVAFullReco_HiggsBJet2RecCSVv2Truth);
//	readerFullRecoTruth->AddVariable("TopHadNonBJetRecCSVv2Truth",&MVAFullReco_TopHadNonBJetRecCSVv2Truth);

	std::string weightsFileFullRecoTruth = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoTruth_BDT.weights.xml";
	readerFullRecoTruth->BookMVA("BDTG method",weightsFileFullRecoTruth.c_str());

	readerFullRecoAll->AddVariable("HiggsRecMAll",&MVAFullReco_HiggsRecMAll);
	readerFullRecoAll->AddVariable("TopLepRecMAll",&MVAFullReco_TopLepRecMAll);
	readerFullRecoAll->AddVariable("HiggsTopLepRecDrAll",&MVAFullReco_HiggsTopLepRecDrAll);
	readerFullRecoAll->AddVariable("TopLepRecPtAll",&MVAFullReco_TopLepRecPtAll);
//	readerFullRecoAll->AddVariable("TopLepBJetRecCSVv2All",&MVAFullReco_TopLepBJetRecCSVv2All);
//	readerFullRecoAll->AddVariable("HiggsBJet1RecCSVv2All",&MVAFullReco_HiggsBJet1RecCSVv2All);
//	readerFullRecoAll->AddVariable("HiggsBJet2RecCSVv2All",&MVAFullReco_HiggsBJet2RecCSVv2All);
//	readerFullRecoAll->AddVariable("TopHadNonBJetRecCSVv2All",&MVAFullReco_TopHadNonBJetRecCSVv2All);

	std::string weightsFileFullRecoAll = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoAll_BDT.weights.xml";
	readerFullRecoAll->BookMVA("BDTG method",weightsFileFullRecoAll.c_str());
	
	readerFullRecoHighestCSVv2->AddVariable("HiggsRecMHighestCSVv2",&MVAFullReco_HiggsRecMHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("TopLepRecMHighestCSVv2",&MVAFullReco_TopLepRecMHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("HiggsTopLepRecDrHighestCSVv2",&MVAFullReco_HiggsTopLepRecDrHighestCSVv2);
	readerFullRecoHighestCSVv2->AddVariable("TopLepRecPtHighestCSVv2",&MVAFullReco_TopLepRecPtHighestCSVv2);
//	readerFullRecoHighestCSVv2->AddVariable("TopLepBJetRecCSVv2HighestCSVv2",&MVAFullReco_TopLepBJetRecCSVv2HighestCSVv2);
//	readerFullRecoHighestCSVv2->AddVariable("HiggsBJet1RecCSVv2HighestCSVv2",&MVAFullReco_HiggsBJet1RecCSVv2HighestCSVv2);
//	readerFullRecoHighestCSVv2->AddVariable("HiggsBJet2RecCSVv2HighestCSVv2",&MVAFullReco_HiggsBJet2RecCSVv2HighestCSVv2);
//	readerFullRecoHighestCSVv2->AddVariable("TopHadNonBJetRecCSVv2HighestCSVv2",&MVAFullReco_TopHadNonBJetRecCSVv2HighestCSVv2);

	std::string weightsFileFullRecoHighestCSVv2 = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoHighestCSVv2_BDT.weights.xml";
	readerFullRecoHighestCSVv2->BookMVA("BDTG method",weightsFileFullRecoHighestCSVv2.c_str());

	readerFullRecoCSVv2L->AddVariable("HiggsRecMCSVv2L",&MVAFullReco_HiggsRecMCSVv2L);
	readerFullRecoCSVv2L->AddVariable("TopLepRecMCSVv2L",&MVAFullReco_TopLepRecMCSVv2L);
	readerFullRecoCSVv2L->AddVariable("HiggsTopLepRecDrCSVv2L",&MVAFullReco_HiggsTopLepRecDrCSVv2L);
	readerFullRecoCSVv2L->AddVariable("TopLepRecPtCSVv2L",&MVAFullReco_TopLepRecPtCSVv2L);
//	readerFullRecoCSVv2L->AddVariable("TopLepBJetRecCSVv2CSVv2L",&MVAFullReco_TopLepBJetRecCSVv2CSVv2L);
//	readerFullRecoCSVv2L->AddVariable("HiggsBJet1RecCSVv2CSVv2L",&MVAFullReco_HiggsBJet1RecCSVv2CSVv2L);
//	readerFullRecoCSVv2L->AddVariable("HiggsBJet2RecCSVv2CSVv2L",&MVAFullReco_HiggsBJet2RecCSVv2CSVv2L);
//	readerFullRecoCSVv2L->AddVariable("TopHadNonBJetRecCSVv2CSVv2L",&MVAFullReco_TopHadNonBJetRecCSVv2CSVv2L);

	std::string weightsFileFullRecoCSVv2L = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoCSVv2L_BDT.weights.xml";
	readerFullRecoCSVv2L->BookMVA("BDTG method",weightsFileFullRecoCSVv2L.c_str());

	readerFullRecoCSVv2M->AddVariable("HiggsRecMCSVv2M",&MVAFullReco_HiggsRecMCSVv2M);
	readerFullRecoCSVv2M->AddVariable("TopLepRecMCSVv2M",&MVAFullReco_TopLepRecMCSVv2M);
	readerFullRecoCSVv2M->AddVariable("HiggsTopLepRecDrCSVv2M",&MVAFullReco_HiggsTopLepRecDrCSVv2M);
	readerFullRecoCSVv2M->AddVariable("TopLepRecPtCSVv2M",&MVAFullReco_TopLepRecPtCSVv2M);
//	readerFullRecoCSVv2M->AddVariable("TopLepBJetRecCSVv2CSVv2M",&MVAFullReco_TopLepBJetRecCSVv2CSVv2M);
//	readerFullRecoCSVv2M->AddVariable("HiggsBJet1RecCSVv2CSVv2M",&MVAFullReco_HiggsBJet1RecCSVv2CSVv2M);
//	readerFullRecoCSVv2M->AddVariable("HiggsBJet2RecCSVv2CSVv2M",&MVAFullReco_HiggsBJet2RecCSVv2CSVv2M);
//	readerFullRecoCSVv2M->AddVariable("TopHadNonBJetRecCSVv2CSVv2M",&MVAFullReco_TopHadNonBJetRecCSVv2CSVv2M);

	std::string weightsFileFullRecoCSVv2M = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoCSVv2M_BDT.weights.xml";
	readerFullRecoCSVv2M->BookMVA("BDTG method",weightsFileFullRecoCSVv2M.c_str());

	readerFullRecoCSVv2T->AddVariable("HiggsRecMCSVv2T",&MVAFullReco_HiggsRecMCSVv2T);
	readerFullRecoCSVv2T->AddVariable("TopLepRecMCSVv2T",&MVAFullReco_TopLepRecMCSVv2T);
	readerFullRecoCSVv2T->AddVariable("HiggsTopLepRecDrCSVv2T",&MVAFullReco_HiggsTopLepRecDrCSVv2T);
	readerFullRecoCSVv2T->AddVariable("TopLepRecPtCSVv2T",&MVAFullReco_TopLepRecPtCSVv2T);
//	readerFullRecoCSVv2T->AddVariable("TopLepBJetRecCSVv2CSVv2T",&MVAFullReco_TopLepBJetRecCSVv2CSVv2T);
//	readerFullRecoCSVv2T->AddVariable("HiggsBJet1RecCSVv2CSVv2T",&MVAFullReco_HiggsBJet1RecCSVv2CSVv2T);
//	readerFullRecoCSVv2T->AddVariable("HiggsBJet2RecCSVv2CSVv2T",&MVAFullReco_HiggsBJet2RecCSVv2CSVv2T);
//	readerFullRecoCSVv2T->AddVariable("TopHadNonBJetRecCSVv2CSVv2T",&MVAFullReco_TopHadNonBJetRecCSVv2CSVv2T);

	std::string weightsFileFullRecoCSVv2T = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAFullRecoCSVv2T_BDT.weights.xml";
	readerFullRecoCSVv2T->BookMVA("BDTG method",weightsFileFullRecoCSVv2T.c_str());

	
	readerPartRecoTruth->AddVariable("HiggsRecMTruth",&MVAPartReco_HiggsRecMTruth);
	readerPartRecoTruth->AddVariable("TopLepRecMTTruth",&MVAPartReco_TopLepRecMTTruth);
	readerPartRecoTruth->AddVariable("HiggsTopLepRecDphiTTruth",&MVAPartReco_HiggsTopLepRecDphiTTruth);
	readerPartRecoTruth->AddVariable("TopLepRecPtTTruth",&MVAPartReco_TopLepRecPtTTruth);
//	readerPartRecoTruth->AddVariable("TopLepBJetRecCSVv2Truth",&MVAPartReco_TopLepBJetRecCSVv2Truth);
//	readerPartRecoTruth->AddVariable("HiggsBJet1RecCSVv2Truth",&MVAPartReco_HiggsBJet1RecCSVv2Truth);
//	readerPartRecoTruth->AddVariable("HiggsBJet2RecCSVv2Truth",&MVAPartReco_HiggsBJet2RecCSVv2Truth);
//	readerPartRecoTruth->AddVariable("TopHadNonBJetRecCSVv2Truth",&MVAPartReco_TopHadNonBJetRecCSVv2Truth);

	std::string weightsFilePartRecoTruth = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoTruth_BDT.weights.xml";
	readerPartRecoTruth->BookMVA("BDTG method",weightsFilePartRecoTruth.c_str());

	readerPartRecoAll->AddVariable("HiggsRecMAll",&MVAPartReco_HiggsRecMAll);
	readerPartRecoAll->AddVariable("TopLepRecMTAll",&MVAPartReco_TopLepRecMTAll);
	readerPartRecoAll->AddVariable("HiggsTopLepRecDphiTAll",&MVAPartReco_HiggsTopLepRecDphiTAll);
	readerPartRecoAll->AddVariable("TopLepRecPtTAll",&MVAPartReco_TopLepRecPtTAll);
//	readerPartRecoAll->AddVariable("TopLepBJetRecCSVv2All",&MVAPartReco_TopLepBJetRecCSVv2All);
//	readerPartRecoAll->AddVariable("HiggsBJet1RecCSVv2All",&MVAPartReco_HiggsBJet1RecCSVv2All);
//	readerPartRecoAll->AddVariable("HiggsBJet2RecCSVv2All",&MVAPartReco_HiggsBJet2RecCSVv2All);
//	readerPartRecoAll->AddVariable("TopHadNonBJetRecCSVv2All",&MVAPartReco_TopHadNonBJetRecCSVv2All);

	std::string weightsFilePartRecoAll = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoAll_BDT.weights.xml";
	readerPartRecoAll->BookMVA("BDTG method",weightsFilePartRecoAll.c_str());
	
	readerPartRecoHighestCSVv2->AddVariable("HiggsRecMHighestCSVv2",&MVAPartReco_HiggsRecMHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("TopLepRecMTHighestCSVv2",&MVAPartReco_TopLepRecMTHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("HiggsTopLepRecDphiTHighestCSVv2",&MVAPartReco_HiggsTopLepRecDphiTHighestCSVv2);
	readerPartRecoHighestCSVv2->AddVariable("TopLepRecPtTHighestCSVv2",&MVAPartReco_TopLepRecPtTHighestCSVv2);
//	readerPartRecoHighestCSVv2->AddVariable("TopLepBJetRecCSVv2HighestCSVv2",&MVAPartReco_TopLepBJetRecCSVv2HighestCSVv2);
//	readerPartRecoHighestCSVv2->AddVariable("HiggsBJet1RecCSVv2HighestCSVv2",&MVAPartReco_HiggsBJet1RecCSVv2HighestCSVv2);
//	readerPartRecoHighestCSVv2->AddVariable("HiggsBJet2RecCSVv2HighestCSVv2",&MVAPartReco_HiggsBJet2RecCSVv2HighestCSVv2);
//	readerPartRecoHighestCSVv2->AddVariable("TopHadNonBJetRecCSVv2HighestCSVv2",&MVAPartReco_TopHadNonBJetRecCSVv2HighestCSVv2);

	std::string weightsFilePartRecoHighestCSVv2 = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoHighestCSVv2_BDT.weights.xml";
	readerPartRecoHighestCSVv2->BookMVA("BDTG method",weightsFilePartRecoHighestCSVv2.c_str());

	readerPartRecoCSVv2L->AddVariable("HiggsRecMCSVv2L",&MVAPartReco_HiggsRecMCSVv2L);
	readerPartRecoCSVv2L->AddVariable("TopLepRecMTCSVv2L",&MVAPartReco_TopLepRecMTCSVv2L);
	readerPartRecoCSVv2L->AddVariable("HiggsTopLepRecDphiTCSVv2L",&MVAPartReco_HiggsTopLepRecDphiTCSVv2L);
	readerPartRecoCSVv2L->AddVariable("TopLepRecPtTCSVv2L",&MVAPartReco_TopLepRecPtTCSVv2L);
//	readerPartRecoCSVv2L->AddVariable("TopLepBJetRecCSVv2CSVv2L",&MVAPartReco_TopLepBJetRecCSVv2CSVv2L);
//	readerPartRecoCSVv2L->AddVariable("HiggsBJet1RecCSVv2CSVv2L",&MVAPartReco_HiggsBJet1RecCSVv2CSVv2L);
//	readerPartRecoCSVv2L->AddVariable("HiggsBJet2RecCSVv2CSVv2L",&MVAPartReco_HiggsBJet2RecCSVv2CSVv2L);
//	readerPartRecoCSVv2L->AddVariable("TopHadNonBJetRecCSVv2CSVv2L",&MVAPartReco_TopHadNonBJetRecCSVv2CSVv2L);

	std::string weightsFilePartRecoCSVv2L = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoCSVv2L_BDT.weights.xml";
	readerPartRecoCSVv2L->BookMVA("BDTG method",weightsFilePartRecoCSVv2L.c_str());

	readerPartRecoCSVv2M->AddVariable("HiggsRecMCSVv2M",&MVAPartReco_HiggsRecMCSVv2M);
	readerPartRecoCSVv2M->AddVariable("TopLepRecMTCSVv2M",&MVAPartReco_TopLepRecMTCSVv2M);
	readerPartRecoCSVv2M->AddVariable("HiggsTopLepRecDphiTCSVv2M",&MVAPartReco_HiggsTopLepRecDphiTCSVv2M);
	readerPartRecoCSVv2M->AddVariable("TopLepRecPtTCSVv2M",&MVAPartReco_TopLepRecPtTCSVv2M);
//	readerPartRecoCSVv2M->AddVariable("TopLepBJetRecCSVv2CSVv2M",&MVAPartReco_TopLepBJetRecCSVv2CSVv2M);
//	readerPartRecoCSVv2M->AddVariable("HiggsBJet1RecCSVv2CSVv2M",&MVAPartReco_HiggsBJet1RecCSVv2CSVv2M);
//	readerPartRecoCSVv2M->AddVariable("HiggsBJet2RecCSVv2CSVv2M",&MVAPartReco_HiggsBJet2RecCSVv2CSVv2M);
//	readerPartRecoCSVv2M->AddVariable("TopHadNonBJetRecCSVv2CSVv2M",&MVAPartReco_TopHadNonBJetRecCSVv2CSVv2M);

	std::string weightsFilePartRecoCSVv2M = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoCSVv2M_BDT.weights.xml";
	readerPartRecoCSVv2M->BookMVA("BDTG method",weightsFilePartRecoCSVv2M.c_str());

	readerPartRecoCSVv2T->AddVariable("HiggsRecMCSVv2T",&MVAPartReco_HiggsRecMCSVv2T);
	readerPartRecoCSVv2T->AddVariable("TopLepRecMTCSVv2T",&MVAPartReco_TopLepRecMTCSVv2T);
	readerPartRecoCSVv2T->AddVariable("HiggsTopLepRecDphiTCSVv2T",&MVAPartReco_HiggsTopLepRecDphiTCSVv2T);
	readerPartRecoCSVv2T->AddVariable("TopLepRecPtTCSVv2T",&MVAPartReco_TopLepRecPtTCSVv2T);
//	readerPartRecoCSVv2T->AddVariable("TopLepBJetRecCSVv2CSVv2T",&MVAPartReco_TopLepBJetRecCSVv2CSVv2T);
//	readerPartRecoCSVv2T->AddVariable("HiggsBJet1RecCSVv2CSVv2T",&MVAPartReco_HiggsBJet1RecCSVv2CSVv2T);
//	readerPartRecoCSVv2T->AddVariable("HiggsBJet2RecCSVv2CSVv2T",&MVAPartReco_HiggsBJet2RecCSVv2CSVv2T);
//	readerPartRecoCSVv2T->AddVariable("TopHadNonBJetRecCSVv2CSVv2T",&MVAPartReco_TopHadNonBJetRecCSVv2CSVv2T);

	std::string weightsFilePartRecoCSVv2T = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/Validation/TopTopLepHbb/MVA/weights/TMVAPartRecoCSVv2T_BDT.weights.xml";
	readerPartRecoCSVv2T->BookMVA("BDTG method",weightsFilePartRecoCSVv2T.c_str());
     }
   
   KINFIT::kfit *kf = new KINFIT::kfit();

   kf->Init(TOPTOPLEPHBB);

   std::string pdfFileName = "/user/kskovpen/analysis/tHFCNC/CMSSW_8_0_12/src/TopKinFit/test/GenAnalysis/TopTopLepHbb/pdf.root";
   kf->SetPDF("TopWMass",pdfFileName.c_str(),"TopLepWM_Fit");
   kf->SetPDF("TopMass",pdfFileName.c_str(),"TopLepRecM_Fit");
   kf->SetPDF("HiggsMass",pdfFileName.c_str(),"HiggsRecM_Fit");
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
	float TopHadNonBJetRecPt, TopHadNonBJetRecEta, TopHadNonBJetRecPhi, TopHadNonBJetRecE, TopHadNonBJetRecCSVv2;
	float HiggsRecPt, HiggsRecEta, HiggsRecPhi, HiggsRecE, HiggsRecM;
	float HiggsBJet1RecPt, HiggsBJet1RecEta, HiggsBJet1RecPhi, HiggsBJet1RecE, HiggsBJet1RecCSVv2;
	float HiggsBJet2RecPt, HiggsBJet2RecEta, HiggsBJet2RecPhi, HiggsBJet2RecE, HiggsBJet2RecCSVv2;
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
	tr->SetBranchAddress("TopHadNonBJetRecPt",&TopHadNonBJetRecPt);
	tr->SetBranchAddress("TopHadNonBJetRecEta",&TopHadNonBJetRecEta);
	tr->SetBranchAddress("TopHadNonBJetRecPhi",&TopHadNonBJetRecPhi);
	tr->SetBranchAddress("TopHadNonBJetRecE",&TopHadNonBJetRecE);
	tr->SetBranchAddress("TopHadNonBJetRecCSVv2",&TopHadNonBJetRecCSVv2);
	tr->SetBranchAddress("HiggsRecPt",&HiggsRecPt);
	tr->SetBranchAddress("HiggsRecEta",&HiggsRecEta);
	tr->SetBranchAddress("HiggsRecPhi",&HiggsRecPhi);
	tr->SetBranchAddress("HiggsRecE",&HiggsRecE);
	tr->SetBranchAddress("HiggsRecM",&HiggsRecM);
	tr->SetBranchAddress("HiggsBJet1RecPt",&HiggsBJet1RecPt);
	tr->SetBranchAddress("HiggsBJet1RecEta",&HiggsBJet1RecEta);
	tr->SetBranchAddress("HiggsBJet1RecPhi",&HiggsBJet1RecPhi);
	tr->SetBranchAddress("HiggsBJet1RecE",&HiggsBJet1RecE);
	tr->SetBranchAddress("HiggsBJet1RecCSVv2",&HiggsBJet1RecCSVv2);
	tr->SetBranchAddress("HiggsBJet2RecPt",&HiggsBJet2RecPt);
	tr->SetBranchAddress("HiggsBJet2RecEta",&HiggsBJet2RecEta);
	tr->SetBranchAddress("HiggsBJet2RecPhi",&HiggsBJet2RecPhi);
	tr->SetBranchAddress("HiggsBJet2RecE",&HiggsBJet2RecE);
	tr->SetBranchAddress("HiggsBJet2RecCSVv2",&HiggsBJet2RecCSVv2);
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
	
	std::vector<float> HiggsRecMTruth;
	std::vector<float> HiggsRecMAll;
	std::vector<float> HiggsRecMHighestCSVv2;
	std::vector<float> HiggsRecMCSVv2L;
	std::vector<float> HiggsRecMCSVv2M;
	std::vector<float> HiggsRecMCSVv2T;

	std::vector<float> MVAHiggsRecMTruth;
	std::vector<float> MVAHiggsRecMAll;
	std::vector<float> MVAHiggsRecMHighestCSVv2;
	std::vector<float> MVAHiggsRecMCSVv2L;
	std::vector<float> MVAHiggsRecMCSVv2M;
	std::vector<float> MVAHiggsRecMCSVv2T;
	
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

	std::vector<float> HiggsBJet1HiggsBJet2RecDrTruth;
	std::vector<float> HiggsBJet1HiggsBJet2RecDrAll;
	std::vector<float> HiggsBJet1HiggsBJet2RecDrHighestCSVv2;
	std::vector<float> HiggsBJet1HiggsBJet2RecDrCSVv2L;
	std::vector<float> HiggsBJet1HiggsBJet2RecDrCSVv2M;
	std::vector<float> HiggsBJet1HiggsBJet2RecDrCSVv2T;

	std::vector<float> HiggsTopLepRecDrTruth;
	std::vector<float> HiggsTopLepRecDrAll;
	std::vector<float> HiggsTopLepRecDrHighestCSVv2;
	std::vector<float> HiggsTopLepRecDrCSVv2L;
	std::vector<float> HiggsTopLepRecDrCSVv2M;
	std::vector<float> HiggsTopLepRecDrCSVv2T;

	std::vector<float> TopLepRecPtTruth;
	std::vector<float> TopLepRecPtAll;
	std::vector<float> TopLepRecPtHighestCSVv2;
	std::vector<float> TopLepRecPtCSVv2L;
	std::vector<float> TopLepRecPtCSVv2M;
	std::vector<float> TopLepRecPtCSVv2T;

	std::vector<float> HiggsRecPtTruth;
	std::vector<float> HiggsRecPtAll;
	std::vector<float> HiggsRecPtHighestCSVv2;
	std::vector<float> HiggsRecPtCSVv2L;
	std::vector<float> HiggsRecPtCSVv2M;
	std::vector<float> HiggsRecPtCSVv2T;
	
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

	std::vector<float> HiggsTopLepRecDphiTTruth;
	std::vector<float> HiggsTopLepRecDphiTAll;
	std::vector<float> HiggsTopLepRecDphiTHighestCSVv2;
	std::vector<float> HiggsTopLepRecDphiTCSVv2L;
	std::vector<float> HiggsTopLepRecDphiTCSVv2M;
	std::vector<float> HiggsTopLepRecDphiTCSVv2T;

	std::vector<float> TopLepRecPtTTruth;
	std::vector<float> TopLepRecPtTAll;
	std::vector<float> TopLepRecPtTHighestCSVv2;
	std::vector<float> TopLepRecPtTCSVv2L;
	std::vector<float> TopLepRecPtTCSVv2M;
	std::vector<float> TopLepRecPtTCSVv2T;

/*	std::vector<float> TopLepBJetRecCSVv2Truth;
	std::vector<float> TopLepBJetRecCSVv2All;
	std::vector<float> TopLepBJetRecCSVv2HighestCSVv2;
	std::vector<float> TopLepBJetRecCSVv2CSVv2L;
	std::vector<float> TopLepBJetRecCSVv2CSVv2M;
	std::vector<float> TopLepBJetRecCSVv2CSVv2T;

	std::vector<float> HiggsBJet1RecCSVv2Truth;
	std::vector<float> HiggsBJet1RecCSVv2All;
	std::vector<float> HiggsBJet1RecCSVv2HighestCSVv2;
	std::vector<float> HiggsBJet1RecCSVv2CSVv2L;
	std::vector<float> HiggsBJet1RecCSVv2CSVv2M;
	std::vector<float> HiggsBJet1RecCSVv2CSVv2T;

	std::vector<float> HiggsBJet2RecCSVv2Truth;
	std::vector<float> HiggsBJet2RecCSVv2All;
	std::vector<float> HiggsBJet2RecCSVv2HighestCSVv2;
	std::vector<float> HiggsBJet2RecCSVv2CSVv2L;
	std::vector<float> HiggsBJet2RecCSVv2CSVv2M;
	std::vector<float> HiggsBJet2RecCSVv2CSVv2T;

	std::vector<float> TopHadNonBJetRecCSVv2Truth;
	std::vector<float> TopHadNonBJetRecCSVv2All;
	std::vector<float> TopHadNonBJetRecCSVv2HighestCSVv2;
	std::vector<float> TopHadNonBJetRecCSVv2CSVv2L;
	std::vector<float> TopHadNonBJetRecCSVv2CSVv2M;
	std::vector<float> TopHadNonBJetRecCSVv2CSVv2T;
*/	
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
	trFIT->Branch("HiggsRecMTruth","std::vector<float>",&HiggsRecMTruth);
	trFIT->Branch("HiggsRecMAll","std::vector<float>",&HiggsRecMAll);
	trFIT->Branch("HiggsRecMHighestCSVv2","std::vector<float>",&HiggsRecMHighestCSVv2);
	trFIT->Branch("HiggsRecMCSVv2L","std::vector<float>",&HiggsRecMCSVv2L);
	trFIT->Branch("HiggsRecMCSVv2M","std::vector<float>",&HiggsRecMCSVv2M);
	trFIT->Branch("HiggsRecMCSVv2T","std::vector<float>",&HiggsRecMCSVv2T);
	trFIT->Branch("MVAHiggsRecMTruth","std::vector<float>",&MVAHiggsRecMTruth);
	trFIT->Branch("MVAHiggsRecMAll","std::vector<float>",&MVAHiggsRecMAll);
	trFIT->Branch("MVAHiggsRecMHighestCSVv2","std::vector<float>",&MVAHiggsRecMHighestCSVv2);
	trFIT->Branch("MVAHiggsRecMCSVv2L","std::vector<float>",&MVAHiggsRecMCSVv2L);
	trFIT->Branch("MVAHiggsRecMCSVv2M","std::vector<float>",&MVAHiggsRecMCSVv2M);
	trFIT->Branch("MVAHiggsRecMCSVv2T","std::vector<float>",&MVAHiggsRecMCSVv2T);
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
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrTruth","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrTruth);
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrAll","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrAll);
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrHighestCSVv2","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrHighestCSVv2);
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrCSVv2L","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrCSVv2L);
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrCSVv2M","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrCSVv2M);
	trFIT->Branch("HiggsBJet1HiggsBJet2RecDrCSVv2T","std::vector<float>",&HiggsBJet1HiggsBJet2RecDrCSVv2T);
	trFIT->Branch("HiggsTopLepRecDrTruth","std::vector<float>",&HiggsTopLepRecDrTruth);
	trFIT->Branch("HiggsTopLepRecDrAll","std::vector<float>",&HiggsTopLepRecDrAll);
	trFIT->Branch("HiggsTopLepRecDrHighestCSVv2","std::vector<float>",&HiggsTopLepRecDrHighestCSVv2);
	trFIT->Branch("HiggsTopLepRecDrCSVv2L","std::vector<float>",&HiggsTopLepRecDrCSVv2L);
	trFIT->Branch("HiggsTopLepRecDrCSVv2M","std::vector<float>",&HiggsTopLepRecDrCSVv2M);
	trFIT->Branch("HiggsTopLepRecDrCSVv2T","std::vector<float>",&HiggsTopLepRecDrCSVv2T);
	trFIT->Branch("TopLepRecPtTruth","std::vector<float>",&TopLepRecPtTruth);
	trFIT->Branch("TopLepRecPtAll","std::vector<float>",&TopLepRecPtAll);
	trFIT->Branch("TopLepRecPtHighestCSVv2","std::vector<float>",&TopLepRecPtHighestCSVv2);
	trFIT->Branch("TopLepRecPtCSVv2L","std::vector<float>",&TopLepRecPtCSVv2L);
	trFIT->Branch("TopLepRecPtCSVv2M","std::vector<float>",&TopLepRecPtCSVv2M);
	trFIT->Branch("TopLepRecPtCSVv2T","std::vector<float>",&TopLepRecPtCSVv2T);
	trFIT->Branch("HiggsRecPtTruth","std::vector<float>",&HiggsRecPtTruth);
	trFIT->Branch("HiggsRecPtAll","std::vector<float>",&HiggsRecPtAll);
	trFIT->Branch("HiggsRecPtHighestCSVv2","std::vector<float>",&HiggsRecPtHighestCSVv2);
	trFIT->Branch("HiggsRecPtCSVv2L","std::vector<float>",&HiggsRecPtCSVv2L);
	trFIT->Branch("HiggsRecPtCSVv2M","std::vector<float>",&HiggsRecPtCSVv2M);
	trFIT->Branch("HiggsRecPtCSVv2T","std::vector<float>",&HiggsRecPtCSVv2T);
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

	trFIT->Branch("HiggsTopLepRecDphiTTruth","std::vector<float>",&HiggsTopLepRecDphiTTruth);
	trFIT->Branch("HiggsTopLepRecDphiTAll","std::vector<float>",&HiggsTopLepRecDphiTAll);
	trFIT->Branch("HiggsTopLepRecDphiTHighestCSVv2","std::vector<float>",&HiggsTopLepRecDphiTHighestCSVv2);
	trFIT->Branch("HiggsTopLepRecDphiTCSVv2L","std::vector<float>",&HiggsTopLepRecDphiTCSVv2L);
	trFIT->Branch("HiggsTopLepRecDphiTCSVv2M","std::vector<float>",&HiggsTopLepRecDphiTCSVv2M);
	trFIT->Branch("HiggsTopLepRecDphiTCSVv2T","std::vector<float>",&HiggsTopLepRecDphiTCSVv2T);

	trFIT->Branch("TopLepRecPtTTruth","std::vector<float>",&TopLepRecPtTTruth);
	trFIT->Branch("TopLepRecPtTAll","std::vector<float>",&TopLepRecPtTAll);
	trFIT->Branch("TopLepRecPtTHighestCSVv2","std::vector<float>",&TopLepRecPtTHighestCSVv2);
	trFIT->Branch("TopLepRecPtTCSVv2L","std::vector<float>",&TopLepRecPtTCSVv2L);
	trFIT->Branch("TopLepRecPtTCSVv2M","std::vector<float>",&TopLepRecPtTCSVv2M);
	trFIT->Branch("TopLepRecPtTCSVv2T","std::vector<float>",&TopLepRecPtTCSVv2T);
/*
	trFIT->Branch("TopLepBJetRecCSVv2Truth","std::vector<float>",&TopLepBJetRecCSVv2Truth);
	trFIT->Branch("TopLepBJetRecCSVv2All","std::vector<float>",&TopLepBJetRecCSVv2All);
	trFIT->Branch("TopLepBJetRecCSVv2HighestCSVv2","std::vector<float>",&TopLepBJetRecCSVv2HighestCSVv2);
	trFIT->Branch("TopLepBJetRecCSVv2CSVv2L","std::vector<float>",&TopLepBJetRecCSVv2CSVv2L);
	trFIT->Branch("TopLepBJetRecCSVv2CSVv2M","std::vector<float>",&TopLepBJetRecCSVv2CSVv2M);
	trFIT->Branch("TopLepBJetRecCSVv2CSVv2T","std::vector<float>",&TopLepBJetRecCSVv2CSVv2T);

	trFIT->Branch("HiggsBJet1RecCSVv2Truth","std::vector<float>",&HiggsBJet1RecCSVv2Truth);
	trFIT->Branch("HiggsBJet1RecCSVv2All","std::vector<float>",&HiggsBJet1RecCSVv2All);
	trFIT->Branch("HiggsBJet1RecCSVv2HighestCSVv2","std::vector<float>",&HiggsBJet1RecCSVv2HighestCSVv2);
	trFIT->Branch("HiggsBJet1RecCSVv2CSVv2L","std::vector<float>",&HiggsBJet1RecCSVv2CSVv2L);
	trFIT->Branch("HiggsBJet1RecCSVv2CSVv2M","std::vector<float>",&HiggsBJet1RecCSVv2CSVv2M);
	trFIT->Branch("HiggsBJet1RecCSVv2CSVv2T","std::vector<float>",&HiggsBJet1RecCSVv2CSVv2T);

	trFIT->Branch("HiggsBJet2RecCSVv2Truth","std::vector<float>",&HiggsBJet2RecCSVv2Truth);
	trFIT->Branch("HiggsBJet2RecCSVv2All","std::vector<float>",&HiggsBJet2RecCSVv2All);
	trFIT->Branch("HiggsBJet2RecCSVv2HighestCSVv2","std::vector<float>",&HiggsBJet2RecCSVv2HighestCSVv2);
	trFIT->Branch("HiggsBJet2RecCSVv2CSVv2L","std::vector<float>",&HiggsBJet2RecCSVv2CSVv2L);
	trFIT->Branch("HiggsBJet2RecCSVv2CSVv2M","std::vector<float>",&HiggsBJet2RecCSVv2CSVv2M);
	trFIT->Branch("HiggsBJet2RecCSVv2CSVv2T","std::vector<float>",&HiggsBJet2RecCSVv2CSVv2T);

	trFIT->Branch("TopHadNonBJetRecCSVv2Truth","std::vector<float>",&TopHadNonBJetRecCSVv2Truth);
	trFIT->Branch("TopHadNonBJetRecCSVv2All","std::vector<float>",&TopHadNonBJetRecCSVv2All);
	trFIT->Branch("TopHadNonBJetRecCSVv2HighestCSVv2","std::vector<float>",&TopHadNonBJetRecCSVv2HighestCSVv2);
	trFIT->Branch("TopHadNonBJetRecCSVv2CSVv2L","std::vector<float>",&TopHadNonBJetRecCSVv2CSVv2L);
	trFIT->Branch("TopHadNonBJetRecCSVv2CSVv2M","std::vector<float>",&TopHadNonBJetRecCSVv2CSVv2M);
	trFIT->Branch("TopHadNonBJetRecCSVv2CSVv2T","std::vector<float>",&TopHadNonBJetRecCSVv2CSVv2T);
*/	
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
	
//	int nMax = 2000;
	int nMax = -1;
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
	     
	     HiggsRecMTruth.clear();
	     HiggsRecMAll.clear();
	     HiggsRecMHighestCSVv2.clear();
	     HiggsRecMCSVv2L.clear();
	     HiggsRecMCSVv2M.clear();
	     HiggsRecMCSVv2T.clear();

	     MVAHiggsRecMTruth.clear();
	     MVAHiggsRecMAll.clear();
	     MVAHiggsRecMHighestCSVv2.clear();
	     MVAHiggsRecMCSVv2L.clear();
	     MVAHiggsRecMCSVv2M.clear();
	     MVAHiggsRecMCSVv2T.clear();
	     
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

	     HiggsBJet1HiggsBJet2RecDrTruth.clear();
	     HiggsBJet1HiggsBJet2RecDrAll.clear();
	     HiggsBJet1HiggsBJet2RecDrHighestCSVv2.clear();
	     HiggsBJet1HiggsBJet2RecDrCSVv2L.clear();
	     HiggsBJet1HiggsBJet2RecDrCSVv2M.clear();
	     HiggsBJet1HiggsBJet2RecDrCSVv2T.clear();

	     HiggsTopLepRecDrTruth.clear();
	     HiggsTopLepRecDrAll.clear();
	     HiggsTopLepRecDrHighestCSVv2.clear();
	     HiggsTopLepRecDrCSVv2L.clear();
	     HiggsTopLepRecDrCSVv2M.clear();
	     HiggsTopLepRecDrCSVv2T.clear();

	     TopLepRecPtTruth.clear();
	     TopLepRecPtAll.clear();
	     TopLepRecPtHighestCSVv2.clear();
	     TopLepRecPtCSVv2L.clear();
	     TopLepRecPtCSVv2M.clear();
	     TopLepRecPtCSVv2T.clear();

	     HiggsRecPtTruth.clear();
	     HiggsRecPtAll.clear();
	     HiggsRecPtHighestCSVv2.clear();
	     HiggsRecPtCSVv2L.clear();
	     HiggsRecPtCSVv2M.clear();
	     HiggsRecPtCSVv2T.clear();
	     
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

	     HiggsTopLepRecDphiTTruth.clear();
	     HiggsTopLepRecDphiTAll.clear();
	     HiggsTopLepRecDphiTHighestCSVv2.clear();
	     HiggsTopLepRecDphiTCSVv2L.clear();
	     HiggsTopLepRecDphiTCSVv2M.clear();
	     HiggsTopLepRecDphiTCSVv2T.clear();

	     TopLepRecPtTTruth.clear();
	     TopLepRecPtTAll.clear();
	     TopLepRecPtTHighestCSVv2.clear();
	     TopLepRecPtTCSVv2L.clear();
	     TopLepRecPtTCSVv2M.clear();
	     TopLepRecPtTCSVv2T.clear();
/*
	     TopLepBJetRecCSVv2Truth.clear();
	     TopLepBJetRecCSVv2All.clear();
	     TopLepBJetRecCSVv2HighestCSVv2.clear();
	     TopLepBJetRecCSVv2CSVv2L.clear();
	     TopLepBJetRecCSVv2CSVv2M.clear();
	     TopLepBJetRecCSVv2CSVv2T.clear();

	     HiggsBJet1RecCSVv2Truth.clear();
	     HiggsBJet1RecCSVv2All.clear();
	     HiggsBJet1RecCSVv2HighestCSVv2.clear();
	     HiggsBJet1RecCSVv2CSVv2L.clear();
	     HiggsBJet1RecCSVv2CSVv2M.clear();
	     HiggsBJet1RecCSVv2CSVv2T.clear();

	     HiggsBJet2RecCSVv2Truth.clear();
	     HiggsBJet2RecCSVv2All.clear();
	     HiggsBJet2RecCSVv2HighestCSVv2.clear();
	     HiggsBJet2RecCSVv2CSVv2L.clear();
	     HiggsBJet2RecCSVv2CSVv2M.clear();
	     HiggsBJet2RecCSVv2CSVv2T.clear();

	     TopHadNonBJetRecCSVv2Truth.clear();
	     TopHadNonBJetRecCSVv2All.clear();
	     TopHadNonBJetRecCSVv2HighestCSVv2.clear();
	     TopHadNonBJetRecCSVv2CSVv2L.clear();
	     TopHadNonBJetRecCSVv2CSVv2M.clear();
	     TopHadNonBJetRecCSVv2CSVv2T.clear();
*/	     
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
		       
		       BJetPt.push_back(HiggsBJet1RecPt);
		       BJetEta.push_back(HiggsBJet1RecEta);
		       BJetPhi.push_back(HiggsBJet1RecPhi);
		       BJetE.push_back(HiggsBJet1RecE);
		       BJetCSVv2.push_back(HiggsBJet1RecCSVv2);
		       
		       BJetPt.push_back(HiggsBJet2RecPt);
		       BJetEta.push_back(HiggsBJet2RecEta);
		       BJetPhi.push_back(HiggsBJet2RecPhi);
		       BJetE.push_back(HiggsBJet2RecE);
		       BJetCSVv2.push_back(HiggsBJet2RecCSVv2);
		       
		       NonBJetPt.push_back(TopHadNonBJetRecPt);
		       NonBJetEta.push_back(TopHadNonBJetRecEta);
		       NonBJetPhi.push_back(TopHadNonBJetRecPhi);
		       NonBJetE.push_back(TopHadNonBJetRecE);
		       NonBJetCSVv2.push_back(TopHadNonBJetRecCSVv2);
		    }
		  else if( it != 5 )
		    {
		       std::vector<std::pair<float,int> > JetCSVv2;
		       JetCSVv2.push_back(std::make_pair(TopLepBJetRecCSVv2,0));
		       JetCSVv2.push_back(std::make_pair(TopHadNonBJetRecCSVv2,1));
		       JetCSVv2.push_back(std::make_pair(HiggsBJet1RecCSVv2,2));
		       JetCSVv2.push_back(std::make_pair(HiggsBJet2RecCSVv2,3));
		       
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
				 if( TopHadNonBJetRecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( TopHadNonBJetRecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( TopHadNonBJetRecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {				 
				      BJetPt.push_back(TopHadNonBJetRecPt);
				      BJetEta.push_back(TopHadNonBJetRecEta);
				      BJetPhi.push_back(TopHadNonBJetRecPhi);
				      BJetE.push_back(TopHadNonBJetRecE);
				      BJetCSVv2.push_back(TopHadNonBJetRecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(TopHadNonBJetRecPt);
				      NonBJetEta.push_back(TopHadNonBJetRecEta);
				      NonBJetPhi.push_back(TopHadNonBJetRecPhi);
				      NonBJetE.push_back(TopHadNonBJetRecE);
				      NonBJetCSVv2.push_back(TopHadNonBJetRecCSVv2);
				   }			    
			      }
			    else if( idx == 2 )
			      {
				 bool passBTag = 1;
				 if( HiggsBJet1RecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( HiggsBJet1RecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( HiggsBJet1RecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {				 
				      BJetPt.push_back(HiggsBJet1RecPt);
				      BJetEta.push_back(HiggsBJet1RecEta);
				      BJetPhi.push_back(HiggsBJet1RecPhi);
				      BJetE.push_back(HiggsBJet1RecE);
				      BJetCSVv2.push_back(HiggsBJet1RecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(HiggsBJet1RecPt);
				      NonBJetEta.push_back(HiggsBJet1RecEta);
				      NonBJetPhi.push_back(HiggsBJet1RecPhi);
				      NonBJetE.push_back(HiggsBJet1RecE);
				      NonBJetCSVv2.push_back(HiggsBJet1RecCSVv2);
				   }			    
			      }		  
			    else if( idx == 3 )
			      {		       
				 bool passBTag = 1;
				 if( HiggsBJet2RecCSVv2 < CSVv2L && it == 2 ) passBTag = 0;
				 if( HiggsBJet2RecCSVv2 < CSVv2M && it == 3 ) passBTag = 0;
				 if( HiggsBJet2RecCSVv2 < CSVv2T && it == 4 ) passBTag = 0;
				 
				 if( ((nBJet < 3 && it == 1) || it > 1) && passBTag )
				   {
				      BJetPt.push_back(HiggsBJet2RecPt);
				      BJetEta.push_back(HiggsBJet2RecEta);
				      BJetPhi.push_back(HiggsBJet2RecPhi);
				      BJetE.push_back(HiggsBJet2RecE);
				      BJetCSVv2.push_back(HiggsBJet2RecCSVv2);
				   }
				 else
				   {
				      NonBJetPt.push_back(HiggsBJet2RecPt);
				      NonBJetEta.push_back(HiggsBJet2RecEta);
				      NonBJetPhi.push_back(HiggsBJet2RecPhi);
				      NonBJetE.push_back(HiggsBJet2RecE);
				      NonBJetCSVv2.push_back(HiggsBJet2RecCSVv2);
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

		       BJetPt.push_back(HiggsBJet1RecPt);
		       BJetEta.push_back(HiggsBJet1RecEta);
		       BJetPhi.push_back(HiggsBJet1RecPhi);
		       BJetE.push_back(HiggsBJet1RecE);
		       BJetCSVv2.push_back(HiggsBJet1RecCSVv2);
		       
		       NonBJetPt.push_back(HiggsBJet1RecPt);
		       NonBJetEta.push_back(HiggsBJet1RecEta);
		       NonBJetPhi.push_back(HiggsBJet1RecPhi);
		       NonBJetE.push_back(HiggsBJet1RecE);
		       NonBJetCSVv2.push_back(HiggsBJet1RecCSVv2);

		       BJetPt.push_back(HiggsBJet2RecPt);
		       BJetEta.push_back(HiggsBJet2RecEta);
		       BJetPhi.push_back(HiggsBJet2RecPhi);
		       BJetE.push_back(HiggsBJet2RecE);
		       BJetCSVv2.push_back(HiggsBJet2RecCSVv2);
		       
		       NonBJetPt.push_back(HiggsBJet2RecPt);
		       NonBJetEta.push_back(HiggsBJet2RecEta);
		       NonBJetPhi.push_back(HiggsBJet2RecPhi);
		       NonBJetE.push_back(HiggsBJet2RecE);
		       NonBJetCSVv2.push_back(HiggsBJet2RecCSVv2);
		       
		       BJetPt.push_back(TopHadNonBJetRecPt);
		       BJetEta.push_back(TopHadNonBJetRecEta);
		       BJetPhi.push_back(TopHadNonBJetRecPhi);
		       BJetE.push_back(TopHadNonBJetRecE);
		       BJetCSVv2.push_back(TopHadNonBJetRecCSVv2);
		       
		       NonBJetPt.push_back(TopHadNonBJetRecPt);
		       NonBJetEta.push_back(TopHadNonBJetRecEta);
		       NonBJetPhi.push_back(TopHadNonBJetRecPhi);
		       NonBJetE.push_back(TopHadNonBJetRecE);
		       NonBJetCSVv2.push_back(TopHadNonBJetRecCSVv2);
		       
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
		  std::vector<float> NonBJetFilteredCSVv2;
		  
		  for(int inb=0;inb<nNonBJetMaxEvent;inb++)
		    {
		       int idx = NonBJetIdx.at(inb);
		       
		       NonBJetFilteredPt.push_back(NonBJetPt.at(idx));
		       NonBJetFilteredEta.push_back(NonBJetEta.at(idx));
		       NonBJetFilteredPhi.push_back(NonBJetPhi.at(idx));
		       NonBJetFilteredE.push_back(NonBJetE.at(idx));
		       NonBJetFilteredCSVv2.push_back(NonBJetCSVv2.at(idx));
		    }		  
		  
		  if( it == 0 ) nEventsTruth++;
		  else if( it == 1 ) nEventsHighestCSVv2++;
		  else if( it == 2 ) nEventsCSVv2L++;
		  else if( it == 3 ) nEventsCSVv2M++;
		  else if( it == 4 ) nEventsCSVv2T++;
		  else if( it == 5 ) nEventsAll++;
		  
		  if( BJetPt.size() < 3 || NonBJetPt.size() < 1 ) continue;
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
		       bool MatchTopHadNonBJet = 0;
		       bool MatchHiggsBJet1 = 0;
		       bool MatchHiggsBJet2 = 0;

		       int idxTopLepWElecFit = kf->GetIndex(ELECTRON_TOPTOPLEPHBB,ip);
		       int idxTopLepWMuonFit = kf->GetIndex(MUON_TOPTOPLEPHBB,ip);
		       int idxTopLepBJetFit = kf->GetIndex(BJETLEP_TOPTOPLEPHBB,ip);
		       int idxTopHadNonBJetFit = kf->GetIndex(NONBJETHAD_TOPTOPLEPHBB,ip);
		       int idxHiggsBJet1Fit = kf->GetIndex(BJET1_TOPTOPLEPHBB,ip);
		       int idxHiggsBJet2Fit = kf->GetIndex(BJET2_TOPTOPLEPHBB,ip);

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
		       float TopLepBJetFitCSVv2 = BJetCSVv2[idxTopLepBJetFit];
		       
		       TLorentzVector *TopLepBJetFitP4 = new TLorentzVector();
		       TopLepBJetFitP4->SetPtEtaPhiE(TopLepBJetFitPt,TopLepBJetFitEta,TopLepBJetFitPhi,TopLepBJetFitE);
		       
		       float TopHadNonBJetFitPt = NonBJetFilteredPt[idxTopHadNonBJetFit];
		       float TopHadNonBJetFitEta = NonBJetFilteredEta[idxTopHadNonBJetFit];
		       float TopHadNonBJetFitPhi = NonBJetFilteredPhi[idxTopHadNonBJetFit];
		       float TopHadNonBJetFitE = NonBJetFilteredE[idxTopHadNonBJetFit];
		       float TopHadNonBJetFitCSVv2 = NonBJetFilteredCSVv2[idxTopHadNonBJetFit];
		       
		       TLorentzVector *TopHadNonBJetFitP4 = new TLorentzVector();
		       TopHadNonBJetFitP4->SetPtEtaPhiE(TopHadNonBJetFitPt,TopHadNonBJetFitEta,TopHadNonBJetFitPhi,TopHadNonBJetFitE);
		       
		       float HiggsBJet1FitPt = BJetPt[idxHiggsBJet1Fit];
		       float HiggsBJet1FitEta = BJetEta[idxHiggsBJet1Fit];
		       float HiggsBJet1FitPhi = BJetPhi[idxHiggsBJet1Fit];
		       float HiggsBJet1FitE = BJetE[idxHiggsBJet1Fit];
		       float HiggsBJet1FitCSVv2 = BJetCSVv2[idxHiggsBJet1Fit];
		       
		       TLorentzVector *HiggsBJet1FitP4 = new TLorentzVector();
		       HiggsBJet1FitP4->SetPtEtaPhiE(HiggsBJet1FitPt,HiggsBJet1FitEta,HiggsBJet1FitPhi,HiggsBJet1FitE);
		       
		       float HiggsBJet2FitPt = BJetPt[idxHiggsBJet2Fit];
		       float HiggsBJet2FitEta = BJetEta[idxHiggsBJet2Fit];
		       float HiggsBJet2FitPhi = BJetPhi[idxHiggsBJet2Fit];
		       float HiggsBJet2FitE = BJetE[idxHiggsBJet2Fit];
		       float HiggsBJet2FitCSVv2 = BJetCSVv2[idxHiggsBJet2Fit];

		       TLorentzVector *HiggsBJet2FitP4 = new TLorentzVector();
		       HiggsBJet2FitP4->SetPtEtaPhiE(HiggsBJet2FitPt,HiggsBJet2FitEta,HiggsBJet2FitPhi,HiggsBJet2FitE);

		       MatchTopLepBJet = (TopLepBJetFitEta-TopLepBJetRecEta < 10E-6 &&
					  TopLepBJetFitPhi-TopLepBJetRecPhi < 10E-6);
		       
		       MatchTopHadNonBJet = (TopHadNonBJetFitEta-TopHadNonBJetRecEta < 10E-6 &&
					     TopHadNonBJetFitPhi-TopHadNonBJetRecPhi < 10E-6);
		       
		       MatchHiggsBJet1 = ((HiggsBJet1FitEta-HiggsBJet1RecEta < 10E-6 &&
					   HiggsBJet1FitPhi-HiggsBJet1RecPhi < 10E-6)) ||
			 ((HiggsBJet1FitEta-HiggsBJet2RecEta < 10E-6 &&
			   HiggsBJet1FitPhi-HiggsBJet2RecPhi < 10E-6));
		       
		       MatchHiggsBJet2 = ((HiggsBJet2FitEta-HiggsBJet1RecEta < 10E-6 &&
					   HiggsBJet2FitPhi-HiggsBJet1RecPhi < 10E-6)) ||
			 ((HiggsBJet2FitEta-HiggsBJet2RecEta < 10E-6 &&
			   HiggsBJet2FitPhi-HiggsBJet2RecPhi < 10E-6));

		       TLorentzVector Higgs = *HiggsBJet1FitP4+*HiggsBJet2FitP4;
		       TLorentzVector TopLep = *TopLepWLepFitP4+*TopLepWNuFitP4+*TopLepBJetFitP4;
		       TLorentzVector TopHad = *HiggsBJet1FitP4+*HiggsBJet2FitP4+*TopHadNonBJetFitP4;

		       float VarHiggsRecM = Higgs.M();
		       float VarTopLepRecM = TopLep.M();
		       float VarTopHadRecM = TopHad.M();
		       float VarHiggsBJet1HiggsBJet2RecDr = HiggsBJet1FitP4->DeltaR(*HiggsBJet2FitP4);
		       float VarHiggsTopLepRecDr = Higgs.DeltaR(TopLep);
		       float VarTopLepRecPt = TopLep.Pt();
		       float VarHiggsRecPt = Higgs.Pt();

		       TLorentzVector *HiggsFitT = new TLorentzVector();
		       HiggsFitT->SetPxPyPzE(Higgs.Px(),Higgs.Py(),0.,Higgs.Et());
		       
		       TLorentzVector *TopLepWLepFitT = new TLorentzVector();
		       TopLepWLepFitT->SetPxPyPzE(TopLepWLepFitP4->Px(),TopLepWLepFitP4->Py(),0.,TopLepWLepFitP4->Et());

		       TLorentzVector *TopLepWNuFitT = new TLorentzVector();
		       TopLepWNuFitT->SetPxPyPzE(TopLepWNuFitP4->Px(),TopLepWNuFitP4->Py(),0.,TopLepWNuFitP4->Et());

		       TLorentzVector *TopLepBJetFitT = new TLorentzVector();
		       TopLepBJetFitT->SetPxPyPzE(TopLepBJetFitP4->Px(),TopLepBJetFitP4->Py(),0.,TopLepBJetFitP4->Et());
		       
		       TLorentzVector TopLepT = *TopLepWLepFitT+*TopLepWNuFitT+*TopLepBJetFitT;

		       float VarTopLepRecMT = TopLepT.M();		       
		       float VarHiggsTopLepRecDphiT = HiggsFitT->DeltaPhi(TopLepT);
		       float VarTopLepRecPtT = TopLepT.Pt();
		       
		       delete HiggsFitT;
		       delete TopLepWLepFitT;
		       delete TopLepWNuFitT;
		       delete TopLepBJetFitT;
		       
		       delete TopLepWLepFitP4;
		       delete TopLepWNuFitP4;
		       delete TopLepBJetFitP4;
		       delete TopHadNonBJetFitP4;
		       delete HiggsBJet1FitP4;
		       delete HiggsBJet2FitP4;
		       
		       bool hasMatch = MatchTopLepBJet && MatchTopHadNonBJet && MatchHiggsBJet1 && MatchHiggsBJet2;
		       
		       bool hasMatchBJet = MatchTopLepBJet && MatchHiggsBJet1 && MatchHiggsBJet2;
		       
		       if( it == 0 ) 
			 {
			    if( hasMatch && disc > 10E+9 && !foundNoSolutionTruth ) {nNoSolutionTruth++; foundNoSolutionTruth=1;}
			    if( ip == 0 ) nSelTruth++;
			    MatchTruth.push_back(hasMatch);
			    MatchBJetTruth.push_back(hasMatchBJet);
			    if( hasMatch && ip == 0 ) nMatchTruth++;
			    if( hasMatchBJet && ip == 0 ) nMatchBJetTruth++;
			    DiscTruth.push_back(disc);
			    HiggsRecMTruth.push_back(VarHiggsRecM);
			    TopLepRecMTruth.push_back(VarTopLepRecM);
			    TopHadRecMTruth.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrTruth.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrTruth.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtTruth.push_back(VarTopLepRecPt);
			    HiggsRecPtTruth.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTTruth.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTTruth.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTTruth.push_back(VarTopLepRecPtT);
			    
/*			    TopLepBJetRecCSVv2Truth.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2Truth.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2Truth.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2Truth.push_back(TopHadNonBJetFitCSVv2);*/

			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_HiggsRecMTruth = VarHiggsRecM;
				      MVAFullReco_TopLepRecMTruth = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrTruth = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtTruth = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2Truth = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2Truth = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2Truth = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2Truth = TopHadNonBJetFitCSVv2;*/
				      
				      MVATruth.push_back(std::make_pair(readerFullRecoTruth->EvaluateMVA("BDTG method"),ip));
				   }				 
				 else
				   {				      
				      MVAPartReco_HiggsRecMTruth = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTTruth = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTTruth = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTTruth = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2Truth = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2Truth = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2Truth = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2Truth = TopHadNonBJetFitCSVv2;*/
				      
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
			    HiggsRecMHighestCSVv2.push_back(VarHiggsRecM);
			    TopLepRecMHighestCSVv2.push_back(VarTopLepRecM);
			    TopHadRecMHighestCSVv2.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrHighestCSVv2.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrHighestCSVv2.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtHighestCSVv2.push_back(VarTopLepRecPt);
			    HiggsRecPtHighestCSVv2.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTHighestCSVv2.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTHighestCSVv2.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTHighestCSVv2.push_back(VarTopLepRecPtT);

/*			    TopLepBJetRecCSVv2HighestCSVv2.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2HighestCSVv2.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2HighestCSVv2.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2HighestCSVv2.push_back(TopHadNonBJetFitCSVv2);*/
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_HiggsRecMHighestCSVv2 = VarHiggsRecM;
				      MVAFullReco_TopLepRecMHighestCSVv2 = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrHighestCSVv2 = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtHighestCSVv2 = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2HighestCSVv2 = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2HighestCSVv2 = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2HighestCSVv2 = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2HighestCSVv2 = TopHadNonBJetFitCSVv2;
*/				      
				      MVAHighestCSVv2.push_back(std::make_pair(readerFullRecoHighestCSVv2->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_HiggsRecMHighestCSVv2 = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTHighestCSVv2 = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTHighestCSVv2 = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTHighestCSVv2 = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2HighestCSVv2 = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2HighestCSVv2 = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2HighestCSVv2 = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2HighestCSVv2 = TopHadNonBJetFitCSVv2;
*/				      
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
			    HiggsRecMCSVv2L.push_back(VarHiggsRecM);
			    TopLepRecMCSVv2L.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2L.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrCSVv2L.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrCSVv2L.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtCSVv2L.push_back(VarTopLepRecPt);
			    HiggsRecPtCSVv2L.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTCSVv2L.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTCSVv2L.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTCSVv2L.push_back(VarTopLepRecPtT);

/*			    TopLepBJetRecCSVv2CSVv2L.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2CSVv2L.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2CSVv2L.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2CSVv2L.push_back(TopHadNonBJetFitCSVv2);*/
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_HiggsRecMCSVv2L = VarHiggsRecM;
				      MVAFullReco_TopLepRecMCSVv2L = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrCSVv2L = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtCSVv2L = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2CSVv2L = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2CSVv2L = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2CSVv2L = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2CSVv2L = TopHadNonBJetFitCSVv2;*/
				 
				      MVACSVv2L.push_back(std::make_pair(readerFullRecoCSVv2L->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_HiggsRecMCSVv2L = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTCSVv2L = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTCSVv2L = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2L = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2CSVv2L = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2CSVv2L = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2CSVv2L = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2CSVv2L = TopHadNonBJetFitCSVv2;*/
				      
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
			    HiggsRecMCSVv2M.push_back(VarHiggsRecM);
			    TopLepRecMCSVv2M.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2M.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrCSVv2M.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrCSVv2M.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtCSVv2M.push_back(VarTopLepRecPt);
			    HiggsRecPtCSVv2M.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTCSVv2M.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTCSVv2M.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTCSVv2M.push_back(VarTopLepRecPtT);

/*			    TopLepBJetRecCSVv2CSVv2M.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2CSVv2M.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2CSVv2M.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2CSVv2M.push_back(TopHadNonBJetFitCSVv2);*/
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {				 
				      MVAFullReco_HiggsRecMCSVv2M = VarHiggsRecM;
				      MVAFullReco_TopLepRecMCSVv2M = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrCSVv2M = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtCSVv2M = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2CSVv2M = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2CSVv2M = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2CSVv2M = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2CSVv2M = TopHadNonBJetFitCSVv2;*/
				      
				      MVACSVv2M.push_back(std::make_pair(readerFullRecoCSVv2M->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {
				      MVAPartReco_HiggsRecMCSVv2M = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTCSVv2M = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTCSVv2M = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2M = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2CSVv2M = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2CSVv2M = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2CSVv2M = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2CSVv2M = TopHadNonBJetFitCSVv2;*/
				      
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
			    HiggsRecMCSVv2T.push_back(VarHiggsRecM);
			    TopLepRecMCSVv2T.push_back(VarTopLepRecM);
			    TopHadRecMCSVv2T.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrCSVv2T.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrCSVv2T.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtCSVv2T.push_back(VarTopLepRecPt);
			    HiggsRecPtCSVv2T.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTCSVv2T.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTCSVv2T.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTCSVv2T.push_back(VarTopLepRecPtT);

/*			    TopLepBJetRecCSVv2CSVv2T.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2CSVv2T.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2CSVv2T.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2CSVv2T.push_back(TopHadNonBJetFitCSVv2);*/
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_HiggsRecMCSVv2T = VarHiggsRecM;
				      MVAFullReco_TopLepRecMCSVv2T = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrCSVv2T = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtCSVv2T = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2CSVv2T = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2CSVv2T = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2CSVv2T = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2CSVv2T = TopHadNonBJetFitCSVv2;*/
				 
				      MVACSVv2T.push_back(std::make_pair(readerFullRecoCSVv2T->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_HiggsRecMCSVv2T = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTCSVv2T = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTCSVv2T = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTCSVv2T = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2CSVv2T = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2CSVv2T = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2CSVv2T = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2CSVv2T = TopHadNonBJetFitCSVv2;*/
				      
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
			    HiggsRecMAll.push_back(VarHiggsRecM);
			    TopLepRecMAll.push_back(VarTopLepRecM);
			    TopHadRecMAll.push_back(VarTopHadRecM);
			    HiggsBJet1HiggsBJet2RecDrAll.push_back(VarHiggsBJet1HiggsBJet2RecDr);
			    HiggsTopLepRecDrAll.push_back(VarHiggsTopLepRecDr);
			    TopLepRecPtAll.push_back(VarTopLepRecPt);
			    HiggsRecPtAll.push_back(VarHiggsRecPt);
			    
			    TopLepRecMTAll.push_back(VarTopLepRecMT);
			    HiggsTopLepRecDphiTAll.push_back(VarHiggsTopLepRecDphiT);
			    TopLepRecPtTAll.push_back(VarTopLepRecPtT);

/*			    TopLepBJetRecCSVv2All.push_back(TopLepBJetFitCSVv2);
			    HiggsBJet1RecCSVv2All.push_back(HiggsBJet1FitCSVv2);
			    HiggsBJet2RecCSVv2All.push_back(HiggsBJet2FitCSVv2);
			    TopHadNonBJetRecCSVv2All.push_back(TopHadNonBJetFitCSVv2);*/
			    
			    if( applyMVA ) 
			      {
				 if( disc < 10E+8 )
				   {
				      MVAFullReco_HiggsRecMAll = VarHiggsRecM;
				      MVAFullReco_TopLepRecMAll = VarTopLepRecM;
				      MVAFullReco_HiggsTopLepRecDrAll = VarHiggsTopLepRecDr;
				      MVAFullReco_TopLepRecPtAll = VarTopLepRecPt;
/*				      MVAFullReco_TopLepBJetRecCSVv2All = TopLepBJetFitCSVv2;
				      MVAFullReco_HiggsBJet1RecCSVv2All = HiggsBJet1FitCSVv2;
				      MVAFullReco_HiggsBJet2RecCSVv2All = HiggsBJet2FitCSVv2;
				      MVAFullReco_TopHadNonBJetRecCSVv2All = TopHadNonBJetFitCSVv2;*/

				      MVAAll.push_back(std::make_pair(readerFullRecoAll->EvaluateMVA("BDTG method"),ip));
				   }
				 else
				   {				      
				      MVAPartReco_HiggsRecMAll = VarHiggsRecM;
				      MVAPartReco_TopLepRecMTAll = VarTopLepRecMT;
				      MVAPartReco_HiggsTopLepRecDphiTAll = VarHiggsTopLepRecDphiT;
				      MVAPartReco_TopLepRecPtTAll = VarTopLepRecPtT;
/*				      MVAPartReco_TopLepBJetRecCSVv2All = TopLepBJetFitCSVv2;
				      MVAPartReco_HiggsBJet1RecCSVv2All = HiggsBJet1FitCSVv2;
				      MVAPartReco_HiggsBJet2RecCSVv2All = HiggsBJet2FitCSVv2;
				      MVAPartReco_TopHadNonBJetRecCSVv2All = TopHadNonBJetFitCSVv2;*/
				      
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
				 MVAHiggsRecMTruth.push_back(HiggsRecMTruth.at(MVATruth.at(ip).second));
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
				 MVAHiggsRecMHighestCSVv2.push_back(HiggsRecMHighestCSVv2.at(MVAHighestCSVv2.at(ip).second));
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
				 MVAHiggsRecMCSVv2L.push_back(HiggsRecMCSVv2L.at(MVACSVv2L.at(ip).second));
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
				 MVAHiggsRecMCSVv2M.push_back(HiggsRecMCSVv2M.at(MVACSVv2M.at(ip).second));
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
				 MVAHiggsRecMCSVv2T.push_back(HiggsRecMCSVv2T.at(MVACSVv2T.at(ip).second));
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
				 MVAHiggsRecMAll.push_back(HiggsRecMAll.at(MVAAll.at(ip).second));
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
 
