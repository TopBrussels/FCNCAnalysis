#include "runMVAPartReco.h"

int main(int argc, char *argv[])
{
   float trainFrac = 0.5;
   std::string coup = "Hct";
   
   // Truth
     {	
	std::string FileName = "TMVAPartReco"+coup+"Truth.root";
	TFile* outfileTruth = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscTruth = new std::vector<float>();
	std::vector<float> *HiggsRecMTruth = new std::vector<float>();
	std::vector<float> *TopLepRecMTTruth = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTTruth = new std::vector<float>();
	std::vector<float> *TopLepRecPtTTruth = new std::vector<float>();
	std::vector<bool> *MatchTruth = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"Truth";
	TMVA::Factory *factoryTruth = new TMVA::Factory(MVAName.c_str(),outfileTruth,
							"!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryTruth->AddVariable("HiggsRecMTruth",'F');
	factoryTruth->AddVariable("TopLepRecMTTruth",'F');
	factoryTruth->AddVariable("HiggsTopLepRecDphiTTruth",'F');
	factoryTruth->AddVariable("TopLepRecPtTTruth",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscTruth",&DiscTruth);
	trFIT.SetBranchAddress("HiggsRecMTruth",&HiggsRecMTruth);
	trFIT.SetBranchAddress("TopLepRecMTTruth",&TopLepRecMTTruth);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTTruth",&HiggsTopLepRecDphiTTruth);
	trFIT.SetBranchAddress("TopLepRecPtTTruth",&TopLepRecPtTTruth);
	trFIT.SetBranchAddress("MatchTruth",&MatchTruth);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermTruth = MatchTruth->size();	
	     for(int ip=0;ip<nPermTruth;ip++)
	       {	     	
		  float DiscTruthVal = DiscTruth->at(ip);
//		  if( DiscTruthVal > 10E+8 )
//		    {
		       std::vector<double> vars(4);
		  
		       float HiggsRecMTruthVal = HiggsRecMTruth->at(ip);
		       float TopLepRecMTTruthVal = TopLepRecMTTruth->at(ip);
		       float HiggsTopLepRecDphiTTruthVal = HiggsTopLepRecDphiTTruth->at(ip);
		       float TopLepRecPtTTruthVal = TopLepRecPtTTruth->at(ip);

		       vars[0] = HiggsRecMTruthVal;
		       vars[1] = TopLepRecMTTruthVal;
		       vars[2] = HiggsTopLepRecDphiTTruthVal;
		       vars[3] = TopLepRecPtTTruthVal;
		       
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
//		    }	     
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
	std::string FileName = "TMVAPartReco"+coup+"HighestCSVv2.root";
	TFile* outfileHighestCSVv2 = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscHighestCSVv2 = new std::vector<float>();
	std::vector<float> *HiggsRecMHighestCSVv2 = new std::vector<float>();
	std::vector<float> *TopLepRecMTHighestCSVv2 = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTHighestCSVv2 = new std::vector<float>();
	std::vector<float> *TopLepRecPtTHighestCSVv2 = new std::vector<float>();
	std::vector<bool> *MatchHighestCSVv2 = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"HighestCSVv2";
	TMVA::Factory *factoryHighestCSVv2 = new TMVA::Factory(MVAName.c_str(),outfileHighestCSVv2,
							       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryHighestCSVv2->AddVariable("HiggsRecMHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("TopLepRecMTHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("HiggsTopLepRecDphiTHighestCSVv2",'F');
	factoryHighestCSVv2->AddVariable("TopLepRecPtTHighestCSVv2",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscHighestCSVv2",&DiscHighestCSVv2);
	trFIT.SetBranchAddress("HiggsRecMHighestCSVv2",&HiggsRecMHighestCSVv2);
	trFIT.SetBranchAddress("TopLepRecMTHighestCSVv2",&TopLepRecMTHighestCSVv2);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTHighestCSVv2",&HiggsTopLepRecDphiTHighestCSVv2);
	trFIT.SetBranchAddress("TopLepRecPtTHighestCSVv2",&TopLepRecPtTHighestCSVv2);
	trFIT.SetBranchAddress("MatchHighestCSVv2",&MatchHighestCSVv2);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermHighestCSVv2 = MatchHighestCSVv2->size();	
	     for(int ip=0;ip<nPermHighestCSVv2;ip++)
	       {	     	
		  float DiscHighestCSVv2Val = DiscHighestCSVv2->at(ip);
//		  if( DiscHighestCSVv2Val > 10E+8 )
//		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMHighestCSVv2Val = HiggsRecMHighestCSVv2->at(ip);
		       float TopLepRecMTHighestCSVv2Val = TopLepRecMTHighestCSVv2->at(ip);
		       float HiggsTopLepRecDphiTHighestCSVv2Val = HiggsTopLepRecDphiTHighestCSVv2->at(ip);
		       float TopLepRecPtTHighestCSVv2Val = TopLepRecPtTHighestCSVv2->at(ip);
		       
		       vars[0] = HiggsRecMHighestCSVv2Val;
		       vars[1] = TopLepRecMTHighestCSVv2Val;
		       vars[2] = HiggsTopLepRecDphiTHighestCSVv2Val;
		       vars[3] = TopLepRecPtTHighestCSVv2Val;
		       
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
//		    }	     
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
	std::string FileName = "TMVAPartReco"+coup+"CSVv2L.root";
	TFile* outfileCSVv2L = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscCSVv2L = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2L = new std::vector<float>();
	std::vector<float> *TopLepRecMTCSVv2L = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTCSVv2L = new std::vector<float>();
	std::vector<float> *TopLepRecPtTCSVv2L = new std::vector<float>();
	std::vector<bool> *MatchCSVv2L = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"CSVv2L";
	TMVA::Factory *factoryCSVv2L = new TMVA::Factory(MVAName.c_str(),outfileCSVv2L,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2L->AddVariable("HiggsRecMCSVv2L",'F');
	factoryCSVv2L->AddVariable("TopLepRecMTCSVv2L",'F');
	factoryCSVv2L->AddVariable("HiggsTopLepRecDphiTCSVv2L",'F');
	factoryCSVv2L->AddVariable("TopLepRecPtTCSVv2L",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2L",&DiscCSVv2L);
	trFIT.SetBranchAddress("HiggsRecMCSVv2L",&HiggsRecMCSVv2L);
	trFIT.SetBranchAddress("TopLepRecMTCSVv2L",&TopLepRecMTCSVv2L);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2L",&HiggsTopLepRecDphiTCSVv2L);
	trFIT.SetBranchAddress("TopLepRecPtTCSVv2L",&TopLepRecPtTCSVv2L);
	trFIT.SetBranchAddress("MatchCSVv2L",&MatchCSVv2L);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2L = MatchCSVv2L->size();	
	     for(int ip=0;ip<nPermCSVv2L;ip++)
	       {	     	
		  float DiscCSVv2LVal = DiscCSVv2L->at(ip);
//		  if( DiscCSVv2LVal > 10E+8 )
//		    {	
		       std::vector<double> vars(4);

		       float HiggsRecMCSVv2LVal = HiggsRecMCSVv2L->at(ip);
		       float TopLepRecMTCSVv2LVal = TopLepRecMTCSVv2L->at(ip);
		       float HiggsTopLepRecDphiTCSVv2LVal = HiggsTopLepRecDphiTCSVv2L->at(ip);
		       float TopLepRecPtTCSVv2LVal = TopLepRecPtTCSVv2L->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2LVal;
		       vars[1] = TopLepRecMTCSVv2LVal;
		       vars[2] = HiggsTopLepRecDphiTCSVv2LVal;
		       vars[3] = TopLepRecPtTCSVv2LVal;
		       
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
//		    }	     
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
	std::string FileName = "TMVAPartReco"+coup+"CSVv2M.root";
	TFile* outfileCSVv2M = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscCSVv2M = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2M = new std::vector<float>();
	std::vector<float> *TopLepRecMTCSVv2M = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTCSVv2M = new std::vector<float>();
	std::vector<float> *TopLepRecPtTCSVv2M = new std::vector<float>();
	std::vector<bool> *MatchCSVv2M = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"CSVv2M";
	TMVA::Factory *factoryCSVv2M = new TMVA::Factory(MVAName.c_str(),outfileCSVv2M,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2M->AddVariable("HiggsRecMCSVv2M",'F');
	factoryCSVv2M->AddVariable("TopLepRecMTCSVv2M",'F');
	factoryCSVv2M->AddVariable("HiggsTopLepRecDphiTCSVv2M",'F');
	factoryCSVv2M->AddVariable("TopLepRecPtTCSVv2M",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2M",&DiscCSVv2M);
	trFIT.SetBranchAddress("HiggsRecMCSVv2M",&HiggsRecMCSVv2M);
	trFIT.SetBranchAddress("TopLepRecMTCSVv2M",&TopLepRecMTCSVv2M);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2M",&HiggsTopLepRecDphiTCSVv2M);
	trFIT.SetBranchAddress("TopLepRecPtTCSVv2M",&TopLepRecPtTCSVv2M);
	trFIT.SetBranchAddress("MatchCSVv2M",&MatchCSVv2M);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2M = MatchCSVv2M->size();
	     for(int ip=0;ip<nPermCSVv2M;ip++)
	       {	     	
		  float DiscCSVv2MVal = DiscCSVv2M->at(ip);
//		  if( DiscCSVv2MVal > 10E+8 )
//		    {	
		       std::vector<double> vars(4);
		       
		       float HiggsRecMCSVv2MVal = HiggsRecMCSVv2M->at(ip);
		       float TopLepRecMTCSVv2MVal = TopLepRecMTCSVv2M->at(ip);
		       float HiggsTopLepRecDphiTCSVv2MVal = HiggsTopLepRecDphiTCSVv2M->at(ip);
		       float TopLepRecPtTCSVv2MVal = TopLepRecPtTCSVv2M->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2MVal;
		       vars[1] = TopLepRecMTCSVv2MVal;
		       vars[2] = HiggsTopLepRecDphiTCSVv2MVal;
		       vars[3] = TopLepRecPtTCSVv2MVal;
		       
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
//		    }	     
	       }
	  }	
   
	factoryCSVv2M->PrepareTrainingAndTestTree("","","SplitMode=Random:NormMode=NumEvents:!V");
	factoryCSVv2M->BookMethod(TMVA::Types::kBDT,"BDT",
				  "!H:!V:NTrees=100:MaxDepth=3:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");
	
	factoryCSVv2M->TrainAllMethods();
	factoryCSVv2M->TestAllMethods();   
	factoryCSVv2M->EvaluateAllMethods();
	
	outfileCSVv2M->Write();
	outfileCSVv2M->Close();
	
	delete factoryCSVv2M;
     }   

   // CSVv2T
     {
	std::string FileName = "TMVAPartReco"+coup+"CSVv2T.root";
	TFile* outfileCSVv2T = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscCSVv2T = new std::vector<float>();
	std::vector<float> *HiggsRecMCSVv2T = new std::vector<float>();
	std::vector<float> *TopLepRecMTCSVv2T = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTCSVv2T = new std::vector<float>();
	std::vector<float> *TopLepRecPtTCSVv2T = new std::vector<float>();
	std::vector<bool> *MatchCSVv2T = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"CSVv2T";
	TMVA::Factory *factoryCSVv2T = new TMVA::Factory(MVAName.c_str(),outfileCSVv2T,
							 "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryCSVv2T->AddVariable("HiggsRecMCSVv2T",'F');
	factoryCSVv2T->AddVariable("TopLepRecMTCSVv2T",'F');
	factoryCSVv2T->AddVariable("HiggsTopLepRecDphiTCSVv2T",'F');
	factoryCSVv2T->AddVariable("TopLepRecPtTCSVv2T",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscCSVv2T",&DiscCSVv2T);
	trFIT.SetBranchAddress("HiggsRecMCSVv2T",&HiggsRecMCSVv2T);
	trFIT.SetBranchAddress("TopLepRecMTCSVv2T",&TopLepRecMTCSVv2T);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2T",&HiggsTopLepRecDphiTCSVv2T);
	trFIT.SetBranchAddress("TopLepRecPtTCSVv2T",&TopLepRecPtTCSVv2T);
	trFIT.SetBranchAddress("MatchCSVv2T",&MatchCSVv2T);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermCSVv2T = MatchCSVv2T->size();	
	     for(int ip=0;ip<nPermCSVv2T;ip++)
	       {	     	
		  float DiscCSVv2TVal = DiscCSVv2T->at(ip);
//		  if( DiscCSVv2TVal > 10E+8 )
//		    {	
		       std::vector<double> vars(4);

		       float HiggsRecMCSVv2TVal = HiggsRecMCSVv2T->at(ip);
		       float TopLepRecMTCSVv2TVal = TopLepRecMTCSVv2T->at(ip);
		       float HiggsTopLepRecDphiTCSVv2TVal = HiggsTopLepRecDphiTCSVv2T->at(ip);
		       float TopLepRecPtTCSVv2TVal = TopLepRecPtTCSVv2T->at(ip);
		       
		       vars[0] = HiggsRecMCSVv2TVal;
		       vars[1] = TopLepRecMTCSVv2TVal;
		       vars[2] = HiggsTopLepRecDphiTCSVv2TVal;
		       vars[3] = TopLepRecPtTCSVv2TVal;
		       
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
//		    }	     
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
	std::string FileName = "TMVAPartReco"+coup+"All.root";
	TFile* outfileAll = TFile::Open(FileName.c_str(),"RECREATE");
	
	std::vector<float> *DiscAll = new std::vector<float>();
	std::vector<float> *HiggsRecMAll = new std::vector<float>();
	std::vector<float> *TopLepRecMTAll = new std::vector<float>();
	std::vector<float> *HiggsTopLepRecDphiTAll = new std::vector<float>();
	std::vector<float> *TopLepRecPtTAll = new std::vector<float>();
	std::vector<bool> *MatchAll = new std::vector<bool>();
	
	std::string MVAName = "TMVAPartReco"+coup+"All";	
	TMVA::Factory *factoryAll = new TMVA::Factory(MVAName.c_str(),outfileAll,
						      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D:AnalysisType=Classification");
	
	factoryAll->AddVariable("HiggsRecMAll",'F');
	factoryAll->AddVariable("TopLepRecMTAll",'F');
	factoryAll->AddVariable("HiggsTopLepRecDphiTAll",'F');
	factoryAll->AddVariable("TopLepRecPtTAll",'F');
	
	TChain trFIT("trFIT");
	coup[0] = tolower(coup[0]);
	std::string f1 = "../runTEST_MERGED/ST_FCNC-TH_Tleptonic_HTobb_eta_"+coup+"-MadGraph5-pythia8/data.root";
	coup[0] = toupper(coup[0]);
	trFIT.Add(f1.c_str());
	
	trFIT.SetBranchAddress("DiscAll",&DiscAll);
	trFIT.SetBranchAddress("HiggsRecMAll",&HiggsRecMAll);
	trFIT.SetBranchAddress("TopLepRecMTAll",&TopLepRecMTAll);
	trFIT.SetBranchAddress("HiggsTopLepRecDphiTAll",&HiggsTopLepRecDphiTAll);
	trFIT.SetBranchAddress("TopLepRecPtTAll",&TopLepRecPtTAll);
	trFIT.SetBranchAddress("MatchAll",&MatchAll);
   
	TRandom3 *r = new TRandom3(666);
	
	for(int i=0;i<trFIT.GetEntries();i++)
	  {
	     trFIT.GetEntry(i);
	     
	     int nPermAll = MatchAll->size();	
	     for(int ip=0;ip<nPermAll;ip++)
	       {	     	
		  float DiscAllVal = DiscAll->at(ip);
//		  if( DiscAllVal > 10E+8 )
//		    {
		       std::vector<double> vars(4);
		       
		       float HiggsRecMAllVal = HiggsRecMAll->at(ip);
		       float TopLepRecMTAllVal = TopLepRecMTAll->at(ip);
		       float HiggsTopLepRecDphiTAllVal = HiggsTopLepRecDphiTAll->at(ip);
		       float TopLepRecPtTAllVal = TopLepRecPtTAll->at(ip);
		       
		       vars[0] = HiggsRecMAllVal;
		       vars[1] = TopLepRecMTAllVal;
		       vars[2] = HiggsTopLepRecDphiTAllVal;
		       vars[3] = TopLepRecPtTAllVal;
		       
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
//		    }	     
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
