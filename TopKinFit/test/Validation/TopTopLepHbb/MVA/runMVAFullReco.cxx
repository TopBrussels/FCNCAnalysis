#include "runMVAFullReco.h"

int main(int argc, char *argv[])
{
   float trainFrac = 0.3;
   std::string coup = "Hut";
   
   // Truth
     {	
	TFile* outfileTruth = TFile::Open("TMVAFullRecoTruth.root","RECREATE");
	
	std::vector<float> *DiscTruth = new std::vector<float>();
	std::vector<float> *HiggsRecMTruth = new std::vector<float>();
	std::vector<float> *TopLepRecMTruth = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrTruth = new std::vector<float>();
	std::vector<float> *TopLepRecPtTruth = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2Truth = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2Truth = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2Truth = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2Truth = new std::vector<float>();
	std::vector<bool> *MatchTruth = new std::vector<bool>();
	
	TMVA::Factory *factoryTruth = new TMVA::Factory("TMVAFullRecoTruth",outfileTruth,
							"!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryTruth->AddVariable("HiggsRecMTruth",'F');
	factoryTruth->AddVariable("TopLepRecMTruth",'F');
	factoryTruth->AddVariable("HiggsTopLepRecDrTruth",'F');
	factoryTruth->AddVariable("TopLepRecPtTruth",'F');
//	factoryTruth->AddVariable("TopLepBJetRecCSVv2Truth",'F');
//	factoryTruth->AddVariable("HiggsBJet1RecCSVv2Truth",'F');
//	factoryTruth->AddVariable("HiggsBJet2RecCSVv2Truth",'F');
//	factoryTruth->AddVariable("TopHadNonBJetRecCSVv2Truth",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscTruth",&DiscTruth);
	trFIT.SetBranchAddress("HiggsRecMTruth",&HiggsRecMTruth);
	trFIT.SetBranchAddress("TopLepRecMTruth",&TopLepRecMTruth);
	trFIT.SetBranchAddress("HiggsTopLepRecDrTruth",&HiggsTopLepRecDrTruth);
	trFIT.SetBranchAddress("TopLepRecPtTruth",&TopLepRecPtTruth);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2Truth",&TopLepBJetRecCSVv2Truth);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2Truth",&HiggsBJet1RecCSVv2Truth);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2Truth",&HiggsBJet2RecCSVv2Truth);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2Truth",&TopHadNonBJetRecCSVv2Truth);
	trFIT.SetBranchAddress("MatchTruth",&MatchTruth);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermTruth = MatchTruth->size();	
	     for(int ip=0;ip<nPermTruth;ip++)
	       {	     	
		  float DiscTruthVal = DiscTruth->at(ip);
		  if( DiscTruthVal < 10E+9 )
		    {
		       std::vector<double> vars(4);
		       
		       float HiggsRecMTruthVal = HiggsRecMTruth->at(ip);
		       float TopLepRecMTruthVal = TopLepRecMTruth->at(ip);
		       float HiggsTopLepRecDrTruthVal = HiggsTopLepRecDrTruth->at(ip);
		       float TopLepRecPtTruthVal = TopLepRecPtTruth->at(ip);
//		       float TopLepBJetRecCSVv2TruthVal = TopLepBJetRecCSVv2Truth->at(ip);
//		       float HiggsBJet1RecCSVv2TruthVal = HiggsBJet1RecCSVv2Truth->at(ip);
//		       float HiggsBJet2RecCSVv2TruthVal = HiggsBJet2RecCSVv2Truth->at(ip);
//		       float TopHadNonBJetRecCSVv2TruthVal = TopHadNonBJetRecCSVv2Truth->at(ip);

		       vars[0] = HiggsRecMTruthVal;
		       vars[1] = TopLepRecMTruthVal;
		       vars[2] = HiggsTopLepRecDrTruthVal;
		       vars[3] = TopLepRecPtTruthVal;
//		       vars[4] = TopLepBJetRecCSVv2TruthVal;
//		       vars[5] = HiggsBJet1RecCSVv2TruthVal;
//		       vars[6] = HiggsBJet2RecCSVv2TruthVal;
//		       vars[7] = TopHadNonBJetRecCSVv2TruthVal;
		       
		       if( MatchTruth->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryTruth->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryTruth->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryTruth->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryTruth->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryTruth->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryTruth->BookMethod(TMVA::Types::kBDT,"BDT",
				 "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryTruth->TrainAllMethods();
	factoryTruth->TestAllMethods();   
	factoryTruth->EvaluateAllMethods();
	
	outfileTruth->Write();
	outfileTruth->Close();
	
	delete factoryTruth;
     }   

   // HighestCSVv2
     {	
	TFile* outfileHighestCSVv2 = TFile::Open("TMVAFullRecoHighestCSVv2.root","RECREATE");
	
	std::vector<float> *DiscHighestCSVv2 = new std::vector<float>();
	std::vector<float> *HiggsRecMHighestCSVv2 = new std::vector<float>();
	std::vector<float> *TopLepRecMHighestCSVv2 = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrHighestCSVv2 = new std::vector<float>();
	std::vector<float> *TopLepRecPtHighestCSVv2 = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2HighestCSVv2 = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2HighestCSVv2 = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2HighestCSVv2 = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2HighestCSVv2 = new std::vector<float>();
	std::vector<bool> *MatchHighestCSVv2 = new std::vector<bool>();
	
	TMVA::Factory *factoryHighestCSVv2 = new TMVA::Factory("TMVAFullRecoHighestCSVv2",outfileHighestCSVv2,
							       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryHighestCSVv2->AddVariable("HiggsRecMHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("TopLepRecMHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("HiggsTopLepRecDrHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("TopLepRecPtHighestCSVv2",'F');
//	factoryHighestCSVv2->AddVariable("TopLepBJetRecCSVv2HighestCSVv2",'F');
//	factoryHighestCSVv2->AddVariable("HiggsBJet1RecCSVv2HighestCSVv2",'F');
//	factoryHighestCSVv2->AddVariable("HiggsBJet2RecCSVv2HighestCSVv2",'F');
//	factoryHighestCSVv2->AddVariable("TopHadNonBJetRecCSVv2HighestCSVv2",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscHighestCSVv2",&DiscHighestCSVv2);
	trFIT.SetBranchAddress("HiggsRecMHighestCSVv2",&HiggsRecMHighestCSVv2);
	trFIT.SetBranchAddress("TopLepRecMHighestCSVv2",&TopLepRecMHighestCSVv2);
	trFIT.SetBranchAddress("HiggsTopLepRecDrHighestCSVv2",&HiggsTopLepRecDrHighestCSVv2);
	trFIT.SetBranchAddress("TopLepRecPtHighestCSVv2",&TopLepRecPtHighestCSVv2);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2HighestCSVv2",&TopLepBJetRecCSVv2HighestCSVv2);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2HighestCSVv2",&HiggsBJet1RecCSVv2HighestCSVv2);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2HighestCSVv2",&HiggsBJet2RecCSVv2HighestCSVv2);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2HighestCSVv2",&TopHadNonBJetRecCSVv2HighestCSVv2);
	trFIT.SetBranchAddress("MatchHighestCSVv2",&MatchHighestCSVv2);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermHighestCSVv2 = MatchHighestCSVv2->size();	
	     for(int ip=0;ip<nPermHighestCSVv2;ip++)
	       {	     	
		  float DiscHighestCSVv2Val = DiscHighestCSVv2->at(ip);
		  if( DiscHighestCSVv2Val < 10E+9 )
		    {
		       std::vector<double> vars(4);
		       
		       float HiggsRecMHighestCSVv2Val = HiggsRecMHighestCSVv2->at(ip);
		       float TopLepRecMHighestCSVv2Val = TopLepRecMHighestCSVv2->at(ip);
		       float HiggsTopLepRecDrHighestCSVv2Val = HiggsTopLepRecDrHighestCSVv2->at(ip);
		       float TopLepRecPtHighestCSVv2Val = TopLepRecPtHighestCSVv2->at(ip);
//		       float TopLepBJetRecCSVv2HighestCSVv2Val = TopLepBJetRecCSVv2HighestCSVv2->at(ip);
//		       float HiggsBJet1RecCSVv2HighestCSVv2Val = HiggsBJet1RecCSVv2HighestCSVv2->at(ip);
//		       float HiggsBJet2RecCSVv2HighestCSVv2Val = HiggsBJet2RecCSVv2HighestCSVv2->at(ip);
//		       float TopHadNonBJetRecCSVv2HighestCSVv2Val = TopHadNonBJetRecCSVv2HighestCSVv2->at(ip);
		       
		       vars[0] = HiggsRecMHighestCSVv2Val;
		       vars[1] = TopLepRecMHighestCSVv2Val;
		       vars[2] = HiggsTopLepRecDrHighestCSVv2Val;
		       vars[3] = TopLepRecPtHighestCSVv2Val;
//		       vars[4] = TopLepBJetRecCSVv2HighestCSVv2Val;
//		       vars[5] = HiggsBJet1RecCSVv2HighestCSVv2Val;
//		       vars[6] = HiggsBJet2RecCSVv2HighestCSVv2Val;
//		       vars[7] = TopHadNonBJetRecCSVv2HighestCSVv2Val;
		       
		       if( MatchHighestCSVv2->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryHighestCSVv2->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryHighestCSVv2->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryHighestCSVv2->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryHighestCSVv2->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryHighestCSVv2->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryHighestCSVv2->BookMethod(TMVA::Types::kBDT,"BDT",
					"!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryHighestCSVv2->TrainAllMethods();
	factoryHighestCSVv2->TestAllMethods();   
	factoryHighestCSVv2->EvaluateAllMethods();
	
	outfileHighestCSVv2->Write();
	outfileHighestCSVv2->Close();
	
	delete factoryHighestCSVv2;
     }   

   // CSVv2L
     {
	TFile* outfileCSVv2L = TFile::Open("TMVAFullRecoCSVv2L.root","RECREATE");
	
	std::vector<float> *DiscCSVv2L = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2L = new std::vector<float>();
	std::vector<float> *TopLepRecMCSVv2L = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrCSVv2L = new std::vector<float>();
	std::vector<float> *TopLepRecPtCSVv2L = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2CSVv2L = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2CSVv2L = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2CSVv2L = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2CSVv2L = new std::vector<float>();
	std::vector<bool> *MatchCSVv2L = new std::vector<bool>();
	
	TMVA::Factory *factoryCSVv2L = new TMVA::Factory("TMVAFullRecoCSVv2L",outfileCSVv2L,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2L->AddVariable("HiggsRecMCSVv2L",'F');
	factoryCSVv2L->AddVariable("TopLepRecMCSVv2L",'F');
	factoryCSVv2L->AddVariable("HiggsTopLepRecDrCSVv2L",'F');
	factoryCSVv2L->AddVariable("TopLepRecPtCSVv2L",'F');
//	factoryCSVv2L->AddVariable("TopLepBJetRecCSVv2CSVv2L",'F');
//	factoryCSVv2L->AddVariable("HiggsBJet1RecCSVv2CSVv2L",'F');
//	factoryCSVv2L->AddVariable("HiggsBJet2RecCSVv2CSVv2L",'F');
//	factoryCSVv2L->AddVariable("TopHadNonBJetRecCSVv2CSVv2L",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2L",&DiscCSVv2L);
	trFIT.SetBranchAddress("HiggsRecMCSVv2L",&HiggsRecMCSVv2L);
	trFIT.SetBranchAddress("TopLepRecMCSVv2L",&TopLepRecMCSVv2L);
	trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2L",&HiggsTopLepRecDrCSVv2L);
	trFIT.SetBranchAddress("TopLepRecPtCSVv2L",&TopLepRecPtCSVv2L);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2CSVv2L",&TopLepBJetRecCSVv2CSVv2L);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2CSVv2L",&HiggsBJet1RecCSVv2CSVv2L);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2CSVv2L",&HiggsBJet2RecCSVv2CSVv2L);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2CSVv2L",&TopHadNonBJetRecCSVv2CSVv2L);
	trFIT.SetBranchAddress("MatchCSVv2L",&MatchCSVv2L);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2L = MatchCSVv2L->size();	
	     for(int ip=0;ip<nPermCSVv2L;ip++)
	       {	     	
		  float DiscCSVv2LVal = DiscCSVv2L->at(ip);
		  if( DiscCSVv2LVal < 10E+9 )
		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMCSVv2LVal = HiggsRecMCSVv2L->at(ip);
		       float TopLepRecMCSVv2LVal = TopLepRecMCSVv2L->at(ip);
		       float HiggsTopLepRecDrCSVv2LVal = HiggsTopLepRecDrCSVv2L->at(ip);
		       float TopLepRecPtCSVv2LVal = TopLepRecPtCSVv2L->at(ip);
//		       float TopLepBJetRecCSVv2CSVv2LVal = TopLepBJetRecCSVv2CSVv2L->at(ip);
//		       float HiggsBJet1RecCSVv2CSVv2LVal = HiggsBJet1RecCSVv2CSVv2L->at(ip);
//		       float HiggsBJet2RecCSVv2CSVv2LVal = HiggsBJet2RecCSVv2CSVv2L->at(ip);
//		       float TopHadNonBJetRecCSVv2CSVv2LVal = TopHadNonBJetRecCSVv2CSVv2L->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2LVal;
		       vars[1] = TopLepRecMCSVv2LVal;
		       vars[2] = HiggsTopLepRecDrCSVv2LVal;
		       vars[3] = TopLepRecPtCSVv2LVal;
//		       vars[4] = TopLepBJetRecCSVv2CSVv2LVal;
//		       vars[5] = HiggsBJet1RecCSVv2CSVv2LVal;
//		       vars[6] = HiggsBJet2RecCSVv2CSVv2LVal;
//		       vars[7] = TopHadNonBJetRecCSVv2CSVv2LVal;
		       
		       if( MatchCSVv2L->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2L->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryCSVv2L->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2L->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryCSVv2L->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryCSVv2L->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryCSVv2L->BookMethod(TMVA::Types::kBDT,"BDT",
				  "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryCSVv2L->TrainAllMethods();
	factoryCSVv2L->TestAllMethods();   
	factoryCSVv2L->EvaluateAllMethods();
	
	outfileCSVv2L->Write();
	outfileCSVv2L->Close();
	
	delete factoryCSVv2L;
     }   

   // CSVv2M
     {
	TFile* outfileCSVv2M = TFile::Open("TMVAFullRecoCSVv2M.root","RECREATE");
	
	std::vector<float> *DiscCSVv2M = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2M = new std::vector<float>();
	std::vector<float> *TopLepRecMCSVv2M = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrCSVv2M = new std::vector<float>();
	std::vector<float> *TopLepRecPtCSVv2M = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2CSVv2M = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2CSVv2M = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2CSVv2M = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2CSVv2M = new std::vector<float>();
	std::vector<bool> *MatchCSVv2M = new std::vector<bool>();
	
	TMVA::Factory *factoryCSVv2M = new TMVA::Factory("TMVAFullRecoCSVv2M",outfileCSVv2M,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2M->AddVariable("HiggsRecMCSVv2M",'F');
	factoryCSVv2M->AddVariable("TopLepRecMCSVv2M",'F');
	factoryCSVv2M->AddVariable("HiggsTopLepRecDrCSVv2M",'F');
	factoryCSVv2M->AddVariable("TopLepRecPtCSVv2M",'F');
//	factoryCSVv2M->AddVariable("TopLepBJetRecCSVv2CSVv2M",'F');
//	factoryCSVv2M->AddVariable("HiggsBJet1RecCSVv2CSVv2M",'F');
//	factoryCSVv2M->AddVariable("HiggsBJet2RecCSVv2CSVv2M",'F');
//	factoryCSVv2M->AddVariable("TopHadNonBJetRecCSVv2CSVv2M",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2M",&DiscCSVv2M);
	trFIT.SetBranchAddress("HiggsRecMCSVv2M",&HiggsRecMCSVv2M);
	trFIT.SetBranchAddress("TopLepRecMCSVv2M",&TopLepRecMCSVv2M);
	trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2M",&HiggsTopLepRecDrCSVv2M);
	trFIT.SetBranchAddress("TopLepRecPtCSVv2M",&TopLepRecPtCSVv2M);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2CSVv2M",&TopLepBJetRecCSVv2CSVv2M);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2CSVv2M",&HiggsBJet1RecCSVv2CSVv2M);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2CSVv2M",&HiggsBJet2RecCSVv2CSVv2M);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2CSVv2M",&TopHadNonBJetRecCSVv2CSVv2M);
	trFIT.SetBranchAddress("MatchCSVv2M",&MatchCSVv2M);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2M = MatchCSVv2M->size();
	     for(int ip=0;ip<nPermCSVv2M;ip++)
	       {	     	
		  float DiscCSVv2MVal = DiscCSVv2M->at(ip);
		  if( DiscCSVv2MVal < 10E+9 )
		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMCSVv2MVal = HiggsRecMCSVv2M->at(ip);
		       float TopLepRecMCSVv2MVal = TopLepRecMCSVv2M->at(ip);
		       float HiggsTopLepRecDrCSVv2MVal = HiggsTopLepRecDrCSVv2M->at(ip);
		       float TopLepRecPtCSVv2MVal = TopLepRecPtCSVv2M->at(ip);
//		       float TopLepBJetRecCSVv2CSVv2MVal = TopLepBJetRecCSVv2CSVv2M->at(ip);
//		       float HiggsBJet1RecCSVv2CSVv2MVal = HiggsBJet1RecCSVv2CSVv2M->at(ip);
//		       float HiggsBJet2RecCSVv2CSVv2MVal = HiggsBJet2RecCSVv2CSVv2M->at(ip);
//		       float TopHadNonBJetRecCSVv2CSVv2MVal = TopHadNonBJetRecCSVv2CSVv2M->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2MVal;
		       vars[1] = TopLepRecMCSVv2MVal;
		       vars[2] = HiggsTopLepRecDrCSVv2MVal;
		       vars[3] = TopLepRecPtCSVv2MVal;
//		       vars[4] = TopLepBJetRecCSVv2CSVv2MVal;
//		       vars[5] = HiggsBJet1RecCSVv2CSVv2MVal;
//		       vars[6] = HiggsBJet2RecCSVv2CSVv2MVal;
//		       vars[7] = TopHadNonBJetRecCSVv2CSVv2MVal;
		       
		       if( MatchCSVv2M->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2M->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryCSVv2M->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2M->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryCSVv2M->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryCSVv2M->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryCSVv2M->BookMethod(TMVA::Types::kBDT,"BDT",
				  "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
//	factoryCSVv2M->BookMethod(TMVA::Types::kBDT,"BDTG",
//				  "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=200:nEventsMin=100:NNodesMax=5:CreateMVAPdfs");
	
	factoryCSVv2M->TrainAllMethods();
	factoryCSVv2M->TestAllMethods();   
	factoryCSVv2M->EvaluateAllMethods();
	
	outfileCSVv2M->Write();
	outfileCSVv2M->Close();
	
	delete factoryCSVv2M;
     }   

   // CSVv2T
     {
	TFile* outfileCSVv2T = TFile::Open("TMVAFullRecoCSVv2T.root","RECREATE");
	
	std::vector<float> *DiscCSVv2T = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2T = new std::vector<float>();
	std::vector<float> *TopLepRecMCSVv2T = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrCSVv2T = new std::vector<float>();
	std::vector<float> *TopLepRecPtCSVv2T = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2CSVv2T = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2CSVv2T = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2CSVv2T = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2CSVv2T = new std::vector<float>();
	std::vector<bool> *MatchCSVv2T = new std::vector<bool>();
	
	TMVA::Factory *factoryCSVv2T = new TMVA::Factory("TMVAFullRecoCSVv2T",outfileCSVv2T,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2T->AddVariable("HiggsRecMCSVv2T",'F');
	factoryCSVv2T->AddVariable("TopLepRecMCSVv2T",'F');
	factoryCSVv2T->AddVariable("HiggsTopLepRecDrCSVv2T",'F');
	factoryCSVv2T->AddVariable("TopLepRecPtCSVv2T",'F');
//	factoryCSVv2T->AddVariable("TopLepBJetRecCSVv2CSVv2T",'F');
//	factoryCSVv2T->AddVariable("HiggsBJet1RecCSVv2CSVv2T",'F');
//	factoryCSVv2T->AddVariable("HiggsBJet2RecCSVv2CSVv2T",'F');
//	factoryCSVv2T->AddVariable("TopHadNonBJetRecCSVv2CSVv2T",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2T",&DiscCSVv2T);
	trFIT.SetBranchAddress("HiggsRecMCSVv2T",&HiggsRecMCSVv2T);
	trFIT.SetBranchAddress("TopLepRecMCSVv2T",&TopLepRecMCSVv2T);
	trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2T",&HiggsTopLepRecDrCSVv2T);
	trFIT.SetBranchAddress("TopLepRecPtCSVv2T",&TopLepRecPtCSVv2T);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2CSVv2T",&TopLepBJetRecCSVv2CSVv2T);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2CSVv2T",&HiggsBJet1RecCSVv2CSVv2T);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2CSVv2T",&HiggsBJet2RecCSVv2CSVv2T);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2CSVv2T",&TopHadNonBJetRecCSVv2CSVv2T);
	trFIT.SetBranchAddress("MatchCSVv2T",&MatchCSVv2T);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2T = MatchCSVv2T->size();	
	     for(int ip=0;ip<nPermCSVv2T;ip++)
	       {	     	
		  float DiscCSVv2TVal = DiscCSVv2T->at(ip);
		  if( DiscCSVv2TVal < 10E+9 )
		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMCSVv2TVal = HiggsRecMCSVv2T->at(ip);
		       float TopLepRecMCSVv2TVal = TopLepRecMCSVv2T->at(ip);
		       float HiggsTopLepRecDrCSVv2TVal = HiggsTopLepRecDrCSVv2T->at(ip);
		       float TopLepRecPtCSVv2TVal = TopLepRecPtCSVv2T->at(ip);
//		       float TopLepBJetRecCSVv2CSVv2TVal = TopLepBJetRecCSVv2CSVv2T->at(ip);
//		       float HiggsBJet1RecCSVv2CSVv2TVal = HiggsBJet1RecCSVv2CSVv2T->at(ip);
//		       float HiggsBJet2RecCSVv2CSVv2TVal = HiggsBJet2RecCSVv2CSVv2T->at(ip);
//		       float TopHadNonBJetRecCSVv2CSVv2TVal = TopHadNonBJetRecCSVv2CSVv2T->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2TVal;
		       vars[1] = TopLepRecMCSVv2TVal;
		       vars[2] = HiggsTopLepRecDrCSVv2TVal;
		       vars[3] = TopLepRecPtCSVv2TVal;
//		       vars[4] = TopLepBJetRecCSVv2CSVv2TVal;
//		       vars[5] = HiggsBJet1RecCSVv2CSVv2TVal;
//		       vars[6] = HiggsBJet2RecCSVv2CSVv2TVal;
//		       vars[7] = TopHadNonBJetRecCSVv2CSVv2TVal;
		       
		       if( MatchCSVv2T->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2T->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryCSVv2T->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryCSVv2T->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryCSVv2T->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryCSVv2T->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryCSVv2T->BookMethod(TMVA::Types::kBDT,"BDT",
				  "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryCSVv2T->TrainAllMethods();
	factoryCSVv2T->TestAllMethods();   
	factoryCSVv2T->EvaluateAllMethods();
	
	outfileCSVv2T->Write();
	outfileCSVv2T->Close();
	
	delete factoryCSVv2T;
     }   

   // All
     {
	TFile* outfileAll = TFile::Open("TMVAFullRecoAll.root","RECREATE");
	
	std::vector<float> *DiscAll = new std::vector<float>();
	std::vector<float> *HiggsRecMAll = new std::vector<float>();
	std::vector<float> *TopLepRecMAll = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDrAll = new std::vector<float>();
	std::vector<float> *TopLepRecPtAll = new std::vector<float>();
//	std::vector<float> *TopLepBJetRecCSVv2All = new std::vector<float>();
//	std::vector<float> *HiggsBJet1RecCSVv2All = new std::vector<float>();
//	std::vector<float> *HiggsBJet2RecCSVv2All = new std::vector<float>();
//	std::vector<float> *TopHadNonBJetRecCSVv2All = new std::vector<float>();
	std::vector<bool> *MatchAll = new std::vector<bool>();
	
	TMVA::Factory *factoryAll = new TMVA::Factory("TMVAFullRecoAll",outfileAll,
						      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryAll->AddVariable("HiggsRecMAll",'F');
	factoryAll->AddVariable("TopLepRecMAll",'F');
	factoryAll->AddVariable("HiggsTopLepRecDrAll",'F');
	factoryAll->AddVariable("TopLepRecPtAll",'F');
//	factoryAll->AddVariable("TopLepBJetRecCSVv2All",'F');
//	factoryAll->AddVariable("HiggsBJet1RecCSVv2All",'F');
//	factoryAll->AddVariable("HiggsBJet2RecCSVv2All",'F');
//	factoryAll->AddVariable("TopHadNonBJetRecCSVv2All",'F');
	
	TChain trFIT("trFIT");
//	trFIT.Add("../output.root");
	std::string f1 = "../runTEST_MERGED/TT_TopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f1.c_str());
	std::string f2 = "../runTEST_MERGED/TT_AntitopLeptonicDecay_TH_1L3B_Eta_"+coup+"/data.root";
	trFIT.Add(f2.c_str());
	
	trFIT.SetBranchAddress("DiscAll",&DiscAll);
	trFIT.SetBranchAddress("HiggsRecMAll",&HiggsRecMAll);
	trFIT.SetBranchAddress("TopLepRecMAll",&TopLepRecMAll);
	trFIT.SetBranchAddress("HiggsTopLepRecDrAll",&HiggsTopLepRecDrAll);
	trFIT.SetBranchAddress("TopLepRecPtAll",&TopLepRecPtAll);
//	trFIT.SetBranchAddress("TopLepBJetRecCSVv2All",&TopLepBJetRecCSVv2All);
//	trFIT.SetBranchAddress("HiggsBJet1RecCSVv2All",&HiggsBJet1RecCSVv2All);
//	trFIT.SetBranchAddress("HiggsBJet2RecCSVv2All",&HiggsBJet2RecCSVv2All);
//	trFIT.SetBranchAddress("TopHadNonBJetRecCSVv2All",&TopHadNonBJetRecCSVv2All);
	trFIT.SetBranchAddress("MatchAll",&MatchAll);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermAll = MatchAll->size();
	     for(int ip=0;ip<nPermAll;ip++)
	       {	     	
		  float DiscAllVal = DiscAll->at(ip);
		  if( DiscAllVal < 10E+9 )
		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMAllVal = HiggsRecMAll->at(ip);
		       float TopLepRecMAllVal = TopLepRecMAll->at(ip);
		       float HiggsTopLepRecDrAllVal = HiggsTopLepRecDrAll->at(ip);
		       float TopLepRecPtAllVal = TopLepRecPtAll->at(ip);
//		       float TopLepBJetRecCSVv2AllVal = TopLepBJetRecCSVv2All->at(ip);
//		       float HiggsBJet1RecCSVv2AllVal = HiggsBJet1RecCSVv2All->at(ip);
//		       float HiggsBJet2RecCSVv2AllVal = HiggsBJet2RecCSVv2All->at(ip);
//		       float TopHadNonBJetRecCSVv2AllVal = TopHadNonBJetRecCSVv2All->at(ip);

		       vars[0] = HiggsRecMAllVal;
		       vars[1] = TopLepRecMAllVal;
		       vars[2] = HiggsTopLepRecDrAllVal;
		       vars[3] = TopLepRecPtAllVal;
//		       vars[4] = TopLepBJetRecCSVv2AllVal;
//		       vars[5] = HiggsBJet1RecCSVv2AllVal;
//		       vars[6] = HiggsBJet2RecCSVv2AllVal;
//		       vars[7] = TopHadNonBJetRecCSVv2AllVal;
		       
		       if( MatchAll->at(ip) )
			 {		  
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryAll->AddSignalTrainingEvent(vars,1.);
			    else
			      factoryAll->AddSignalTestEvent(vars,1.);
			 }
		       else
			 {	     
			    float rnd = r->Rndm();
			    if( rnd < trainFrac )
			      factoryAll->AddBackgroundTrainingEvent(vars,1.);
			    else
			      factoryAll->AddBackgroundTestEvent(vars,1.);
			 }	     
		    }	     
	       }
	  }	
   
	factoryAll->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryAll->BookMethod(TMVA::Types::kBDT,"BDT",
			       "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryAll->TrainAllMethods();
	factoryAll->TestAllMethods();   
	factoryAll->EvaluateAllMethods();
	
	outfileAll->Write();
	outfileAll->Close();
	
	delete factoryAll;
     }   
}
