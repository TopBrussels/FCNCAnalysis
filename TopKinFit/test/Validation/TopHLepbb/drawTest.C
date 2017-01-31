#include "TStyle.h"

TStyle *tdrStyle;

void addbin(TH1D *h);

void setTDRStyle();

float errf(float v1,float ve1,float v2,float ve2);

void drawTest()
{
   gROOT->SetBatch();
   
   setTDRStyle();
   
   TChain trFIT("trFIT");

//   trFIT.Add("output.root");
   trFIT.Add("runTEST_MERGED/ST_TH_1L3B_Hct/data.root");

   TChain trEFF("trEFF");

//   trEFF.Add("output.root");
   trEFF.Add("runTEST_MERGED/ST_TH_1L3B_Hct/data.root");
   
   std::vector<float> *DiscTruth = new std::vector<float>();
   std::vector<float> *MVADiscTruth = new std::vector<float>();
   std::vector<float> *MVAScoreTruth = new std::vector<float>();
   std::vector<float> *MVAHiggsRecMTruth = new std::vector<float>();
   std::vector<float> *MVATopLepRecMTruth = new std::vector<float>();
   
   std::vector<bool> *MatchTruth = new std::vector<bool>();
   std::vector<float> *HiggsRecMTruth = new std::vector<float>();
   std::vector<float> *TopLepRecMTruth = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrTruth = new std::vector<float>();
   std::vector<float> *TopLepRecPtTruth = new std::vector<float>(); 
//   std::vector<float> *HiggsBJet1HiggsBJet2RecDrTruth = new std::vector<float>();
   std::vector<float> *HiggsRecPtTruth = new std::vector<float>();
   std::vector<float> *TopLepRecMTTruth = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTTruth = new std::vector<float>();
   std::vector<float> *TopLepRecPtTTruth = new std::vector<float>();

   std::vector<float> *DiscHighestCSVv2 = new std::vector<float>();
   std::vector<bool> *MatchHighestCSVv2 = new std::vector<bool>();
   std::vector<float> *HiggsRecMHighestCSVv2 = new std::vector<float>();
   std::vector<float> *TopLepRecMHighestCSVv2 = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrHighestCSVv2 = new std::vector<float>();
   std::vector<float> *TopLepRecPtHighestCSVv2 = new std::vector<float>(); 
   std::vector<float> *HiggsRecPtHighestCSVv2 = new std::vector<float>();
   std::vector<float> *TopLepRecMTHighestCSVv2 = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTHighestCSVv2 = new std::vector<float>();
   std::vector<float> *TopLepRecPtTHighestCSVv2 = new std::vector<float>();

   std::vector<float> *DiscCSVv2L = new std::vector<float>();
   std::vector<bool> *MatchCSVv2L = new std::vector<bool>();
   std::vector<float> *HiggsRecMCSVv2L = new std::vector<float>();
   std::vector<float> *TopLepRecMCSVv2L = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrCSVv2L = new std::vector<float>();
   std::vector<float> *TopLepRecPtCSVv2L = new std::vector<float>(); 
   std::vector<float> *HiggsRecPtCSVv2L = new std::vector<float>();
   std::vector<float> *TopLepRecMTCSVv2L = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTCSVv2L = new std::vector<float>();
   std::vector<float> *TopLepRecPtTCSVv2L = new std::vector<float>();

   std::vector<float> *DiscCSVv2M = new std::vector<float>();
   std::vector<bool> *MatchCSVv2M = new std::vector<bool>();
   std::vector<float> *HiggsRecMCSVv2M = new std::vector<float>();
   std::vector<float> *TopLepRecMCSVv2M = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrCSVv2M = new std::vector<float>();
   std::vector<float> *TopLepRecPtCSVv2M = new std::vector<float>(); 
   std::vector<float> *HiggsRecPtCSVv2M = new std::vector<float>();
   std::vector<float> *TopLepRecMTCSVv2M = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTCSVv2M = new std::vector<float>();
   std::vector<float> *TopLepRecPtTCSVv2M = new std::vector<float>();

   std::vector<float> *DiscCSVv2T = new std::vector<float>();
   std::vector<bool> *MatchCSVv2T = new std::vector<bool>();
   std::vector<float> *HiggsRecMCSVv2T = new std::vector<float>();
   std::vector<float> *TopLepRecMCSVv2T = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrCSVv2T = new std::vector<float>();
   std::vector<float> *TopLepRecPtCSVv2T = new std::vector<float>(); 
   std::vector<float> *HiggsRecPtCSVv2T = new std::vector<float>();
   std::vector<float> *TopLepRecMTCSVv2T = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTCSVv2T = new std::vector<float>();
   std::vector<float> *TopLepRecPtTCSVv2T = new std::vector<float>();

   std::vector<float> *DiscAll = new std::vector<float>();
   std::vector<bool> *MatchAll = new std::vector<bool>();
   std::vector<float> *HiggsRecMAll = new std::vector<float>();
   std::vector<float> *TopLepRecMAll = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDrAll = new std::vector<float>();
   std::vector<float> *TopLepRecPtAll = new std::vector<float>(); 
   std::vector<float> *HiggsRecPtAll = new std::vector<float>();
   std::vector<float> *TopLepRecMTAll = new std::vector<float>();
   std::vector<float> *HiggsTopLepRecDphiTAll = new std::vector<float>();
   std::vector<float> *TopLepRecPtTAll = new std::vector<float>();
   
   float nEventsTruth;
   float nEventsHighestCSVv2;
   float nEventsCSVv2L;
   float nEventsCSVv2M;
   float nEventsCSVv2T;
   float nEventsAll;

   float nMatchTruth;
   float nMatchHighestCSVv2;
   float nMatchCSVv2L;
   float nMatchCSVv2M;
   float nMatchCSVv2T;
   float nMatchAll;

   float nMatchBJetTruth;
   float nMatchBJetHighestCSVv2;
   float nMatchBJetCSVv2L;
   float nMatchBJetCSVv2M;
   float nMatchBJetCSVv2T;
   float nMatchBJetAll;
   
   float nMatchMVATruth;
   float nMatchMVAHighestCSVv2;
   float nMatchMVACSVv2L;
   float nMatchMVACSVv2M;
   float nMatchMVACSVv2T;
   float nMatchMVAAll;

   float nMatchBJetMVATruth;
   float nMatchBJetMVAHighestCSVv2;
   float nMatchBJetMVACSVv2L;
   float nMatchBJetMVACSVv2M;
   float nMatchBJetMVACSVv2T;
   float nMatchBJetMVAAll;
   
   float nSelTruth;
   float nSelHighestCSVv2;
   float nSelCSVv2L;
   float nSelCSVv2M;
   float nSelCSVv2T;
   float nSelAll;

   float nNoSolutionTruth;
   float nNoSolutionHighestCSVv2;
   float nNoSolutionCSVv2L;
   float nNoSolutionCSVv2M;
   float nNoSolutionCSVv2T;
   float nNoSolutionAll;
   
   trFIT.SetBranchAddress("DiscTruth",&DiscTruth);
   trFIT.SetBranchAddress("MVADiscTruth",&MVADiscTruth);
   trFIT.SetBranchAddress("MVAScoreTruth",&MVAScoreTruth);
   trFIT.SetBranchAddress("MVAHiggsRecMTruth",&MVAHiggsRecMTruth);
   trFIT.SetBranchAddress("MVATopLepRecMTruth",&MVATopLepRecMTruth);
   
   trFIT.SetBranchAddress("MatchTruth",&MatchTruth);
   trFIT.SetBranchAddress("HiggsRecMTruth",&HiggsRecMTruth);
   trFIT.SetBranchAddress("TopLepRecMTruth",&TopLepRecMTruth);
   trFIT.SetBranchAddress("HiggsTopLepRecDrTruth",&HiggsTopLepRecDrTruth);
   trFIT.SetBranchAddress("TopLepRecPtTruth",&TopLepRecPtTruth);
//   trFIT.SetBranchAddress("HiggsBJet1HiggsBJet2RecDrTruth",&HiggsBJet1HiggsBJet2RecDrTruth);
   trFIT.SetBranchAddress("HiggsRecPtTruth",&HiggsRecPtTruth);
   trFIT.SetBranchAddress("TopLepRecMTTruth",&TopLepRecMTTruth);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTTruth",&HiggsTopLepRecDphiTTruth);
   trFIT.SetBranchAddress("TopLepRecPtTTruth",&TopLepRecPtTTruth);

   trFIT.SetBranchAddress("DiscHighestCSVv2",&DiscHighestCSVv2);
   trFIT.SetBranchAddress("MatchHighestCSVv2",&MatchHighestCSVv2);
   trFIT.SetBranchAddress("HiggsRecMHighestCSVv2",&HiggsRecMHighestCSVv2);
   trFIT.SetBranchAddress("TopLepRecMHighestCSVv2",&TopLepRecMHighestCSVv2);
   trFIT.SetBranchAddress("HiggsTopLepRecDrHighestCSVv2",&HiggsTopLepRecDrHighestCSVv2);
   trFIT.SetBranchAddress("TopLepRecPtHighestCSVv2",&TopLepRecPtHighestCSVv2);
   trFIT.SetBranchAddress("HiggsRecPtHighestCSVv2",&HiggsRecPtHighestCSVv2);
   trFIT.SetBranchAddress("TopLepRecMTHighestCSVv2",&TopLepRecMTHighestCSVv2);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTHighestCSVv2",&HiggsTopLepRecDphiTHighestCSVv2);
   trFIT.SetBranchAddress("TopLepRecPtTHighestCSVv2",&TopLepRecPtTHighestCSVv2);

   trFIT.SetBranchAddress("DiscCSVv2L",&DiscCSVv2L);
   trFIT.SetBranchAddress("MatchCSVv2L",&MatchCSVv2L);
   trFIT.SetBranchAddress("HiggsRecMCSVv2L",&HiggsRecMCSVv2L);
   trFIT.SetBranchAddress("TopLepRecMCSVv2L",&TopLepRecMCSVv2L);
   trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2L",&HiggsTopLepRecDrCSVv2L);
   trFIT.SetBranchAddress("TopLepRecPtCSVv2L",&TopLepRecPtCSVv2L);
   trFIT.SetBranchAddress("HiggsRecPtCSVv2L",&HiggsRecPtCSVv2L);
   trFIT.SetBranchAddress("TopLepRecMTCSVv2L",&TopLepRecMTCSVv2L);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2L",&HiggsTopLepRecDphiTCSVv2L);
   trFIT.SetBranchAddress("TopLepRecPtTCSVv2L",&TopLepRecPtTCSVv2L);

   trFIT.SetBranchAddress("DiscCSVv2M",&DiscCSVv2M);
   trFIT.SetBranchAddress("MatchCSVv2M",&MatchCSVv2M);
   trFIT.SetBranchAddress("HiggsRecMCSVv2M",&HiggsRecMCSVv2M);
   trFIT.SetBranchAddress("TopLepRecMCSVv2M",&TopLepRecMCSVv2M);
   trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2M",&HiggsTopLepRecDrCSVv2M);
   trFIT.SetBranchAddress("TopLepRecPtCSVv2M",&TopLepRecPtCSVv2M);
   trFIT.SetBranchAddress("HiggsRecPtCSVv2M",&HiggsRecPtCSVv2M);
   trFIT.SetBranchAddress("TopLepRecMTCSVv2M",&TopLepRecMTCSVv2M);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2M",&HiggsTopLepRecDphiTCSVv2M);
   trFIT.SetBranchAddress("TopLepRecPtTCSVv2M",&TopLepRecPtTCSVv2M);

   trFIT.SetBranchAddress("DiscCSVv2T",&DiscCSVv2T);
   trFIT.SetBranchAddress("MatchCSVv2T",&MatchCSVv2T);
   trFIT.SetBranchAddress("HiggsRecMCSVv2T",&HiggsRecMCSVv2T);
   trFIT.SetBranchAddress("TopLepRecMCSVv2T",&TopLepRecMCSVv2T);
   trFIT.SetBranchAddress("HiggsTopLepRecDrCSVv2T",&HiggsTopLepRecDrCSVv2T);
   trFIT.SetBranchAddress("TopLepRecPtCSVv2T",&TopLepRecPtCSVv2T);
   trFIT.SetBranchAddress("HiggsRecPtCSVv2T",&HiggsRecPtCSVv2T);
   trFIT.SetBranchAddress("TopLepRecMTCSVv2T",&TopLepRecMTCSVv2T);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTCSVv2T",&HiggsTopLepRecDphiTCSVv2T);
   trFIT.SetBranchAddress("TopLepRecPtTCSVv2T",&TopLepRecPtTCSVv2T);

   trFIT.SetBranchAddress("DiscAll",&DiscAll);
   trFIT.SetBranchAddress("MatchAll",&MatchAll);
   trFIT.SetBranchAddress("HiggsRecMAll",&HiggsRecMAll);
   trFIT.SetBranchAddress("TopLepRecMAll",&TopLepRecMAll);
   trFIT.SetBranchAddress("HiggsTopLepRecDrAll",&HiggsTopLepRecDrAll);
   trFIT.SetBranchAddress("TopLepRecPtAll",&TopLepRecPtAll);
   trFIT.SetBranchAddress("HiggsRecPtAll",&HiggsRecPtAll);
   trFIT.SetBranchAddress("TopLepRecMTAll",&TopLepRecMTAll);
   trFIT.SetBranchAddress("HiggsTopLepRecDphiTAll",&HiggsTopLepRecDphiTAll);
   trFIT.SetBranchAddress("TopLepRecPtTAll",&TopLepRecPtTAll);
   
   trEFF.SetBranchAddress("nEventsTruth",&nEventsTruth);
   trEFF.SetBranchAddress("nEventsHighestCSVv2",&nEventsHighestCSVv2);
   trEFF.SetBranchAddress("nEventsCSVv2L",&nEventsCSVv2L);
   trEFF.SetBranchAddress("nEventsCSVv2M",&nEventsCSVv2M);
   trEFF.SetBranchAddress("nEventsCSVv2T",&nEventsCSVv2T);
   trEFF.SetBranchAddress("nEventsAll",&nEventsAll);

   trEFF.SetBranchAddress("nMatchTruth",&nMatchTruth);
   trEFF.SetBranchAddress("nMatchHighestCSVv2",&nMatchHighestCSVv2);
   trEFF.SetBranchAddress("nMatchCSVv2L",&nMatchCSVv2L);
   trEFF.SetBranchAddress("nMatchCSVv2M",&nMatchCSVv2M);
   trEFF.SetBranchAddress("nMatchCSVv2T",&nMatchCSVv2T);
   trEFF.SetBranchAddress("nMatchAll",&nMatchAll);

   trEFF.SetBranchAddress("nMatchBJetTruth",&nMatchBJetTruth);
   trEFF.SetBranchAddress("nMatchBJetHighestCSVv2",&nMatchBJetHighestCSVv2);
   trEFF.SetBranchAddress("nMatchBJetCSVv2L",&nMatchBJetCSVv2L);
   trEFF.SetBranchAddress("nMatchBJetCSVv2M",&nMatchBJetCSVv2M);
   trEFF.SetBranchAddress("nMatchBJetCSVv2T",&nMatchBJetCSVv2T);
   trEFF.SetBranchAddress("nMatchBJetAll",&nMatchBJetAll);
   
   trEFF.SetBranchAddress("nMatchMVATruth",&nMatchMVATruth);
   trEFF.SetBranchAddress("nMatchMVAHighestCSVv2",&nMatchMVAHighestCSVv2);
   trEFF.SetBranchAddress("nMatchMVACSVv2L",&nMatchMVACSVv2L);
   trEFF.SetBranchAddress("nMatchMVACSVv2M",&nMatchMVACSVv2M);
   trEFF.SetBranchAddress("nMatchMVACSVv2T",&nMatchMVACSVv2T);
   trEFF.SetBranchAddress("nMatchMVAAll",&nMatchMVAAll);

   trEFF.SetBranchAddress("nMatchBJetMVATruth",&nMatchBJetMVATruth);
   trEFF.SetBranchAddress("nMatchBJetMVAHighestCSVv2",&nMatchBJetMVAHighestCSVv2);
   trEFF.SetBranchAddress("nMatchBJetMVACSVv2L",&nMatchBJetMVACSVv2L);
   trEFF.SetBranchAddress("nMatchBJetMVACSVv2M",&nMatchBJetMVACSVv2M);
   trEFF.SetBranchAddress("nMatchBJetMVACSVv2T",&nMatchBJetMVACSVv2T);
   trEFF.SetBranchAddress("nMatchBJetMVAAll",&nMatchBJetMVAAll);
   
   trEFF.SetBranchAddress("nSelTruth",&nSelTruth);
   trEFF.SetBranchAddress("nSelHighestCSVv2",&nSelHighestCSVv2);
   trEFF.SetBranchAddress("nSelCSVv2L",&nSelCSVv2L);
   trEFF.SetBranchAddress("nSelCSVv2M",&nSelCSVv2M);
   trEFF.SetBranchAddress("nSelCSVv2T",&nSelCSVv2T);
   trEFF.SetBranchAddress("nSelAll",&nSelAll);

   trEFF.SetBranchAddress("nNoSolutionTruth",&nNoSolutionTruth);
   trEFF.SetBranchAddress("nNoSolutionHighestCSVv2",&nNoSolutionHighestCSVv2);
   trEFF.SetBranchAddress("nNoSolutionCSVv2L",&nNoSolutionCSVv2L);
   trEFF.SetBranchAddress("nNoSolutionCSVv2M",&nNoSolutionCSVv2M);
   trEFF.SetBranchAddress("nNoSolutionCSVv2T",&nNoSolutionCSVv2T);
   trEFF.SetBranchAddress("nNoSolutionAll",&nNoSolutionAll);
   
   int nHist = 0;
   
   TH1D *hMatch[10000];
   TH1D *hNotMatch[10000];
   std::string hName[10000];
   std::string hLab[10000];

   hMatch[nHist] = new TH1D("h_DiscTruthMatch","h_DiscTruthMatch",30,0.,150.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_DiscTruthNotMatch","h_DiscTruthNotMatch",30,0.,150.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "DiscTruth";
   hLab[nHist] = "D";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecMTruthMatch","h_HiggsRecMTruthMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMTruthNotMatch","h_HiggsRecMTruthNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMTruth";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTruthMatch","h_TopLepRecMTruthMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTruthNotMatch","h_TopLepRecMTruthNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTruth";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrTruthMatch","h_HiggsTopLepRecDrTruthMatch",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrTruthNotMatch","h_HiggsTopLepRecDrTruthNotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrTruth";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTruthMatch","h_TopLepRecPtTruthMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTruthNotMatch","h_TopLepRecPtTruthNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTruth";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

//   hMatch[nHist] = new TH1D("h_HiggsBJet1HiggsBJet2RecDrTruthMatch","h_HiggsBJet1HiggsBJet2RecDrTruthMatch",30,0.,5.);
//   hMatch[nHist]->Sumw2();
//   hNotMatch[nHist] = new TH1D("h_HiggsBJet1HiggsBJet2RecDrTruthNotMatch","h_HiggsBJet1HiggsBJet2RecDrTruthNotMatch",30,0.,5.);
//   hNotMatch[nHist]->Sumw2();
//   hName[nHist] = "HiggsBJet1HiggsBJet2RecDrTruth";
//   hLab[nHist] = "#Delta R (HiggsB1,HiggsB2)";
//   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtTruthMatch","h_HiggsRecPtTruthMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtTruthNotMatch","h_HiggsRecPtTruthNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtTruth";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTTruthMatch","h_TopLepRecMTTruthMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTTruthNotMatch","h_TopLepRecMTTruthNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTTruth";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTTruthMatch","h_HiggsTopLepRecDphiTTruthMatch",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTTruthNotMatch","h_HiggsTopLepRecDphiTTruthNotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTTruth";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTTruthMatch","h_TopLepRecPtTTruthMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTTruthNotMatch","h_TopLepRecPtTTruthNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTTruth";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   
   
   hMatch[nHist] = new TH1D("h_HiggsRecMHighestCSVv2Match","h_HiggsRecMHighestCSVv2Match",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMHighestCSVv2NotMatch","h_HiggsRecMHighestCSVv2NotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMHighestCSVv2";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMHighestCSVv2Match","h_TopLepRecMHighestCSVv2Match",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMHighestCSVv2NotMatch","h_TopLepRecMHighestCSVv2NotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMHighestCSVv2";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrHighestCSVv2Match","h_HiggsTopLepRecDrHighestCSVv2Match",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrHighestCSVv2NotMatch","h_HiggsTopLepRecDrHighestCSVv2NotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrHighestCSVv2";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtHighestCSVv2Match","h_TopLepRecPtHighestCSVv2Match",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtHighestCSVv2NotMatch","h_TopLepRecPtHighestCSVv2NotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtHighestCSVv2";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtHighestCSVv2Match","h_HiggsRecPtHighestCSVv2Match",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtHighestCSVv2NotMatch","h_HiggsRecPtHighestCSVv2NotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtHighestCSVv2";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTHighestCSVv2Match","h_TopLepRecMTHighestCSVv2Match",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTHighestCSVv2NotMatch","h_TopLepRecMTHighestCSVv2NotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTHighestCSVv2";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTHighestCSVv2Match","h_HiggsTopLepRecDphiTHighestCSVv2Match",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTHighestCSVv2NotMatch","h_HiggsTopLepRecDphiTHighestCSVv2NotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTHighestCSVv2";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTHighestCSVv2Match","h_TopLepRecPtTHighestCSVv2Match",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTHighestCSVv2NotMatch","h_TopLepRecPtTHighestCSVv2NotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTHighestCSVv2";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   
   
   hMatch[nHist] = new TH1D("h_HiggsRecMCSVv2LMatch","h_HiggsRecMCSVv2LMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMCSVv2LNotMatch","h_HiggsRecMCSVv2LNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMCSVv2L";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMCSVv2LMatch","h_TopLepRecMCSVv2LMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMCSVv2LNotMatch","h_TopLepRecMCSVv2LNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMCSVv2L";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2LMatch","h_HiggsTopLepRecDrCSVv2LMatch",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2LNotMatch","h_HiggsTopLepRecDrCSVv2LNotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrCSVv2L";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2LMatch","h_TopLepRecPtCSVv2LMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2LNotMatch","h_TopLepRecPtCSVv2LNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtCSVv2L";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2LMatch","h_HiggsRecPtCSVv2LMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2LNotMatch","h_HiggsRecPtCSVv2LNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtCSVv2L";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2LMatch","h_TopLepRecMTCSVv2LMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2LNotMatch","h_TopLepRecMTCSVv2LNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTCSVv2L";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2LMatch","h_HiggsTopLepRecDphiTCSVv2LMatch",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2LNotMatch","h_HiggsTopLepRecDphiTCSVv2LNotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTCSVv2L";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2LMatch","h_TopLepRecPtTCSVv2LMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2LNotMatch","h_TopLepRecPtTCSVv2LNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTCSVv2L";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   
   
   
   hMatch[nHist] = new TH1D("h_HiggsRecMCSVv2MMatch","h_HiggsRecMCSVv2MMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMCSVv2MNotMatch","h_HiggsRecMCSVv2MNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMCSVv2M";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMCSVv2MMatch","h_TopLepRecMCSVv2MMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMCSVv2MNotMatch","h_TopLepRecMCSVv2MNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMCSVv2M";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2MMatch","h_HiggsTopLepRecDrCSVv2MMatch",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2MNotMatch","h_HiggsTopLepRecDrCSVv2MNotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrCSVv2M";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2MMatch","h_TopLepRecPtCSVv2MMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2MNotMatch","h_TopLepRecPtCSVv2MNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtCSVv2M";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2MMatch","h_HiggsRecPtCSVv2MMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2MNotMatch","h_HiggsRecPtCSVv2MNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtCSVv2M";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2MMatch","h_TopLepRecMTCSVv2MMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2MNotMatch","h_TopLepRecMTCSVv2MNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTCSVv2M";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2MMatch","h_HiggsTopLepRecDphiTCSVv2MMatch",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2MNotMatch","h_HiggsTopLepRecDphiTCSVv2MNotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTCSVv2M";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2MMatch","h_TopLepRecPtTCSVv2MMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2MNotMatch","h_TopLepRecPtTCSVv2MNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTCSVv2M";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   
   

   hMatch[nHist] = new TH1D("h_HiggsRecMCSVv2TMatch","h_HiggsRecMCSVv2TMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMCSVv2TNotMatch","h_HiggsRecMCSVv2TNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMCSVv2T";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMCSVv2TMatch","h_TopLepRecMCSVv2TMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMCSVv2TNotMatch","h_TopLepRecMCSVv2TNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMCSVv2T";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2TMatch","h_HiggsTopLepRecDrCSVv2TMatch",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrCSVv2TNotMatch","h_HiggsTopLepRecDrCSVv2TNotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrCSVv2T";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2TMatch","h_TopLepRecPtCSVv2TMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtCSVv2TNotMatch","h_TopLepRecPtCSVv2TNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtCSVv2T";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2TMatch","h_HiggsRecPtCSVv2TMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtCSVv2TNotMatch","h_HiggsRecPtCSVv2TNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtCSVv2T";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2TMatch","h_TopLepRecMTCSVv2TMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTCSVv2TNotMatch","h_TopLepRecMTCSVv2TNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTCSVv2T";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2TMatch","h_HiggsTopLepRecDphiTCSVv2TMatch",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTCSVv2TNotMatch","h_HiggsTopLepRecDphiTCSVv2TNotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTCSVv2T";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2TMatch","h_TopLepRecPtTCSVv2TMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTCSVv2TNotMatch","h_TopLepRecPtTCSVv2TNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTCSVv2T";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   
   
   
   hMatch[nHist] = new TH1D("h_HiggsRecMAllMatch","h_HiggsRecMAllMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecMAllNotMatch","h_HiggsRecMAllNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecMAll";
   hLab[nHist] = "m(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMAllMatch","h_TopLepRecMAllMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMAllNotMatch","h_TopLepRecMAllNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMAll";
   hLab[nHist] = "m(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDrAllMatch","h_HiggsTopLepRecDrAllMatch",30,0.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDrAllNotMatch","h_HiggsTopLepRecDrAllNotMatch",30,0.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDrAll";
   hLab[nHist] = "#Delta R (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtAllMatch","h_TopLepRecPtAllMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtAllNotMatch","h_TopLepRecPtAllNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtAll";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsRecPtAllMatch","h_HiggsRecPtAllMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsRecPtAllNotMatch","h_HiggsRecPtAllNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsRecPtAll";
   hLab[nHist] = "p_{T}(Higgs) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecMTAllMatch","h_TopLepRecMTAllMatch",30,0.,500.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecMTAllNotMatch","h_TopLepRecMTAllNotMatch",30,0.,500.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecMTAll";
   hLab[nHist] = "m_{T}(TopLep) [GeV]";
   nHist++;

   hMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTAllMatch","h_HiggsTopLepRecDphiTAllMatch",30,-5.,5.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_HiggsTopLepRecDphiTAllNotMatch","h_HiggsTopLepRecDphiTAllNotMatch",30,-5.,5.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "HiggsTopLepRecDphiTAll";
   hLab[nHist] = "#Delta #phi (Higgs,TopLep)";
   nHist++;

   hMatch[nHist] = new TH1D("h_TopLepRecPtTAllMatch","h_TopLepRecPtTAllMatch",30,0.,250.);
   hMatch[nHist]->Sumw2();
   hNotMatch[nHist] = new TH1D("h_TopLepRecPtTAllNotMatch","h_TopLepRecPtTAllNotMatch",30,0.,250.);
   hNotMatch[nHist]->Sumw2();
   hName[nHist] = "TopLepRecPtTAll";
   hLab[nHist] = "p_{T}(TopLep) [GeV]";
   nHist++;
   
   int nHistFit = 0;

   TH1D *hFit[1000];
   std::string hNameFit[1000];
   std::string hLabFit[1000];

   hFit[nHistFit] = new TH1D("h_MVADiscTruthFit","h_MVADiscTruthFit",30,0.,150.);
   hFit[nHistFit]->Sumw2();
   hNameFit[nHistFit] = "MVADiscTruthFit";
   hLabFit[nHistFit] = "D";
   nHistFit++;

   hFit[nHistFit] = new TH1D("h_MVAScoreTruthFit","h_MVAScoreTruthFit",30,-1.,1.);
   hFit[nHistFit]->Sumw2();
   hNameFit[nHistFit] = "MVAScoreTruthFit";
   hLabFit[nHistFit] = "MVA output";
   nHistFit++;
   
   hFit[nHistFit] = new TH1D("h_MVAHiggsRecMTruthFit","h_MVAHiggsRecMTruthFit",30,0.,250.);
   hFit[nHistFit]->Sumw2();
   hNameFit[nHistFit] = "MVAHiggsRecMTruthFit";
   hLabFit[nHistFit] = "m(Higgs) [GeV]";
   nHistFit++;

   hFit[nHistFit] = new TH1D("h_MVATopLepRecMTruthFit","h_MVATopLepRecMTruthFit",30,0.,500.);
   hFit[nHistFit]->Sumw2();
   hNameFit[nHistFit] = "MVATopLepRecMTruthFit";
   hLabFit[nHistFit] = "m(TopLep) [GeV]";
   nHistFit++;
   
   int nev = trFIT.GetEntries();
   
   for(int i=0;i<nev;i++)
     {
	trFIT.GetEntry(i);

	for(int ih=0;ih<nHist;ih++)
	  {
	     for(int ip=0;ip<MatchTruth->size();ip++)
	       {
		  if( DiscTruth->at(ip) < 10E+10 )
		    {				 
		       if( hName[ih] == "DiscTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(DiscTruth->at(ip));
			    else hNotMatch[ih]->Fill(DiscTruth->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecMTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(HiggsRecMTruth->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMTruth->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(TopLepRecMTruth->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTruth->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrTruth->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrTruth->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTruth->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTruth->at(ip));
			 }
//		       else if( hName[ih] == "HiggsBJet1HiggsBJet2RecDrTruth" )
//			 {
//			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(HiggsBJet1HiggsBJet2RecDrTruth->at(ip));
//			    else hNotMatch[ih]->Fill(HiggsBJet1HiggsBJet2RecDrTruth->at(ip));
//			 }
		       else if( hName[ih] == "HiggsRecPtTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(HiggsRecPtTruth->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtTruth->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(TopLepRecMTTruth->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTTruth->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTTruth->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTTruth->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTTruth" )
			 {
			    if( MatchTruth->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTTruth->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTTruth->at(ip));
			 }
		    }
	       }
	     
	     for(int ip=0;ip<MatchHighestCSVv2->size();ip++)
	       {		       
		  if( DiscHighestCSVv2->at(ip) < 10E+10 )
		    {				 
		       if( hName[ih] == "HiggsRecMHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(HiggsRecMHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(TopLepRecMHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(TopLepRecPtHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecPtHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(HiggsRecPtHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(TopLepRecMTHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTHighestCSVv2->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTHighestCSVv2" )
			 {
			    if( MatchHighestCSVv2->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTHighestCSVv2->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTHighestCSVv2->at(ip));
			 }
		    }
	       }	     
		    
	     for(int ip=0;ip<MatchCSVv2L->size();ip++)
	       {		       
		  if( DiscCSVv2L->at(ip) < 10E+10 )
		    {				 
		       if( hName[ih] == "HiggsRecMCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(HiggsRecMCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(TopLepRecMCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(TopLepRecPtCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecPtCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(HiggsRecPtCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(TopLepRecMTCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2L->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTCSVv2L" )
			 {
			    if( MatchCSVv2L->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTCSVv2L->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTCSVv2L->at(ip));
			 }
		    }
	       }
	     
	     for(int ip=0;ip<MatchCSVv2M->size();ip++)
	       {		       		  
		  if( DiscCSVv2M->at(ip) < 10E+10 )
		    {
		       if( hName[ih] == "HiggsRecMCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(HiggsRecMCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(TopLepRecMCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(TopLepRecPtCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecPtCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(HiggsRecPtCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(TopLepRecMTCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2M->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTCSVv2M" )
			 {
			    if( MatchCSVv2M->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTCSVv2M->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTCSVv2M->at(ip));
			 }
		    }
	       }	    
		  		    		    
	     for(int ip=0;ip<MatchCSVv2T->size();ip++)
	       {		       
		  if( DiscCSVv2T->at(ip) < 10E+10 )
		    {
		       if( hName[ih] == "HiggsRecMCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(HiggsRecMCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(TopLepRecMCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(TopLepRecPtCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecPtCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(HiggsRecPtCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(TopLepRecMTCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTCSVv2T->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTCSVv2T" )
			 {
			    if( MatchCSVv2T->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTCSVv2T->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTCSVv2T->at(ip));
			 }
		    }
	       }	     
		  		    
	     for(int ip=0;ip<MatchAll->size();ip++)
	       {		       
		  if( DiscAll->at(ip) < 10E+10 )
		    {		    
		       if( hName[ih] == "HiggsRecMAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(HiggsRecMAll->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecMAll->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(TopLepRecMAll->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMAll->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDrAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDrAll->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDrAll->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(TopLepRecPtAll->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtAll->at(ip));
			 }
		       else if( hName[ih] == "HiggsRecPtAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(HiggsRecPtAll->at(ip));
			    else hNotMatch[ih]->Fill(HiggsRecPtAll->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecMTAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(TopLepRecMTAll->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecMTAll->at(ip));
			 }
		       else if( hName[ih] == "HiggsTopLepRecDphiTAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(HiggsTopLepRecDphiTAll->at(ip));
			    else hNotMatch[ih]->Fill(HiggsTopLepRecDphiTAll->at(ip));
			 }
		       else if( hName[ih] == "TopLepRecPtTAll" )
			 {
			    if( MatchAll->at(ip) ) hMatch[ih]->Fill(TopLepRecPtTAll->at(ip));
			    else hNotMatch[ih]->Fill(TopLepRecPtTAll->at(ip));
			 }
		    }		  
	       }
	  }
	
	for(int ih=0;ih<nHistFit;ih++)
	  {
	     if( hNameFit[ih] == "MVADiscTruthFit" )
	       {
		  hFit[ih]->Fill(MVADiscTruth->at(0));
	       }
	     else if( hNameFit[ih] == "MVAScoreTruthFit" )
	       {
		  hFit[ih]->Fill(MVAScoreTruth->at(0));
	       }
	     else if( hNameFit[ih] == "MVAHiggsRecMTruthFit" )
	       {
		  hFit[ih]->Fill(MVAHiggsRecMTruth->at(0));
	       }
	     else if( hNameFit[ih] == "MVATopLepRecMTruthFit" )
	       {
		  hFit[ih]->Fill(MVATopLepRecMTruth->at(0));
	       }
	  }		
     }

   float nEventsSumTruth = 0;
   float nEventsSumHighestCSVv2 = 0;
   float nEventsSumCSVv2L = 0;
   float nEventsSumCSVv2M = 0;
   float nEventsSumCSVv2T = 0;
   float nEventsSumAll = 0;

   float nMatchSumTruth = 0;
   float nMatchSumHighestCSVv2 = 0;
   float nMatchSumCSVv2L = 0;
   float nMatchSumCSVv2M = 0;
   float nMatchSumCSVv2T = 0;
   float nMatchSumAll = 0;

   float nMatchBJetSumTruth = 0;
   float nMatchBJetSumHighestCSVv2 = 0;
   float nMatchBJetSumCSVv2L = 0;
   float nMatchBJetSumCSVv2M = 0;
   float nMatchBJetSumCSVv2T = 0;
   float nMatchBJetSumAll = 0;
   
   float nMatchMVASumTruth = 0;
   float nMatchMVASumHighestCSVv2 = 0;
   float nMatchMVASumCSVv2L = 0;
   float nMatchMVASumCSVv2M = 0;
   float nMatchMVASumCSVv2T = 0;
   float nMatchMVASumAll = 0;

   float nMatchBJetMVASumTruth = 0;
   float nMatchBJetMVASumHighestCSVv2 = 0;
   float nMatchBJetMVASumCSVv2L = 0;
   float nMatchBJetMVASumCSVv2M = 0;
   float nMatchBJetMVASumCSVv2T = 0;
   float nMatchBJetMVASumAll = 0;
   
   float nSelSumTruth = 0;
   float nSelSumHighestCSVv2 = 0;
   float nSelSumCSVv2L = 0;
   float nSelSumCSVv2M = 0;
   float nSelSumCSVv2T = 0;
   float nSelSumAll = 0;

   float nNoSolutionSumTruth = 0;
   float nNoSolutionSumHighestCSVv2 = 0;
   float nNoSolutionSumCSVv2L = 0;
   float nNoSolutionSumCSVv2M = 0;
   float nNoSolutionSumCSVv2T = 0;
   float nNoSolutionSumAll = 0;
   
   int nevEFF = trEFF.GetEntries();
   
   for(int i=0;i<nevEFF;i++)
     {
	trEFF.GetEntry(i);

	nEventsSumTruth += nEventsTruth;
	nEventsSumHighestCSVv2 += nEventsHighestCSVv2;
	nEventsSumCSVv2L += nEventsCSVv2L;
	nEventsSumCSVv2M += nEventsCSVv2M;
	nEventsSumCSVv2T += nEventsCSVv2T;
	nEventsSumAll += nEventsAll;
	
	nMatchSumTruth += nMatchTruth;
	nMatchSumHighestCSVv2 += nMatchHighestCSVv2;
	nMatchSumCSVv2L += nMatchCSVv2L;
	nMatchSumCSVv2M += nMatchCSVv2M;
	nMatchSumCSVv2T += nMatchCSVv2T;
	nMatchSumAll += nMatchAll;

	nMatchBJetSumTruth += nMatchBJetTruth;
	nMatchBJetSumHighestCSVv2 += nMatchBJetHighestCSVv2;
	nMatchBJetSumCSVv2L += nMatchBJetCSVv2L;
	nMatchBJetSumCSVv2M += nMatchBJetCSVv2M;
	nMatchBJetSumCSVv2T += nMatchBJetCSVv2T;
	nMatchBJetSumAll += nMatchBJetAll;
	
	nMatchMVASumTruth += nMatchMVATruth;
	nMatchMVASumHighestCSVv2 += nMatchMVAHighestCSVv2;
	nMatchMVASumCSVv2L += nMatchMVACSVv2L;
	nMatchMVASumCSVv2M += nMatchMVACSVv2M;
	nMatchMVASumCSVv2T += nMatchMVACSVv2T;
	nMatchMVASumAll += nMatchMVAAll;

	nMatchBJetMVASumTruth += nMatchBJetMVATruth;
	nMatchBJetMVASumHighestCSVv2 += nMatchBJetMVAHighestCSVv2;
	nMatchBJetMVASumCSVv2L += nMatchBJetMVACSVv2L;
	nMatchBJetMVASumCSVv2M += nMatchBJetMVACSVv2M;
	nMatchBJetMVASumCSVv2T += nMatchBJetMVACSVv2T;
	nMatchBJetMVASumAll += nMatchBJetMVAAll;
	
	nSelSumTruth += nSelTruth;
	nSelSumHighestCSVv2 += nSelHighestCSVv2;
	nSelSumCSVv2L += nSelCSVv2L;
	nSelSumCSVv2M += nSelCSVv2M;
	nSelSumCSVv2T += nSelCSVv2T;
	nSelSumAll += nSelAll;
	
	nNoSolutionSumTruth += nNoSolutionTruth;
	nNoSolutionSumHighestCSVv2 += nNoSolutionHighestCSVv2;
	nNoSolutionSumCSVv2L += nNoSolutionCSVv2L;
	nNoSolutionSumCSVv2M += nNoSolutionCSVv2M;
	nNoSolutionSumCSVv2T += nNoSolutionCSVv2T;
	nNoSolutionSumAll += nNoSolutionAll;
     }
   
   // Plots
   
   TCanvas *c1 = new TCanvas();
   c1->Draw();
   c1->cd();

   TPad *c1_1;
   
   gStyle->SetHistTopMargin(0);

   for(int i=0;i<nHist;i++)
     {	
	addbin(hMatch[i]);
	addbin(hNotMatch[i]);

	hMatch[i]->Scale(1./hMatch[i]->Integral());
	hNotMatch[i]->Scale(1./hNotMatch[i]->Integral());
   
	hMatch[i]->SetLineWidth(2);
	hMatch[i]->SetLineColor(kRed);
	hMatch[i]->SetMarkerColor(kRed);
	hMatch[i]->SetMarkerStyle(20);
	
	hMatch[i]->Draw("hist e1");

	hNotMatch[i]->SetLineWidth(2);
	hNotMatch[i]->SetLineColor(kBlue);
	hNotMatch[i]->SetMarkerColor(kBlue);
	hNotMatch[i]->SetMarkerStyle(22);
	
	hNotMatch[i]->Draw("hist e1 same");
	
	hMatch[i]->GetXaxis()->SetTitle(hLab[i].c_str());
	hMatch[i]->GetYaxis()->SetTitle("Normalized to unity");
	
	float maxMatch = hMatch[i]->GetMaximum();
	float maxNotMatch = hNotMatch[i]->GetMaximum();
	float max = std::max(maxMatch,maxNotMatch);
	
	hMatch[i]->SetMaximum(1.2*max);

	TLegend *leg = new TLegend(0.70,0.90,0.90,0.70);
	leg->SetFillColor(253);
	leg->SetBorderSize(0);
	leg->AddEntry(hMatch[i],"Right","lp");
	leg->AddEntry(hNotMatch[i],"Wrong","lp");
	leg->Draw();	
	
	std::string figName = "pics/"+hName[i]+".eps";
	c1->Print(figName.c_str());
	c1->Clear();
	delete leg;
     }   

   for(int i=0;i<nHistFit;i++)
     {	
	addbin(hFit[i]);

	hFit[i]->Scale(1./hFit[i]->Integral());
   
	hFit[i]->SetLineWidth(2);
	hFit[i]->SetLineColor(kBlack);
	hFit[i]->SetMarkerColor(kBlack);
	hFit[i]->SetMarkerStyle(20);
	
	hFit[i]->Draw("hist e1");
	
	hFit[i]->GetXaxis()->SetTitle(hLabFit[i].c_str());
	hFit[i]->GetYaxis()->SetTitle("Normalized to unity");
	
	float max = hFit[i]->GetMaximum();
	
	hFit[i]->SetMaximum(1.2*max);

	std::string figName = "pics/"+hNameFit[i]+".eps";
	c1->Print(figName.c_str());
	c1->Clear();
     }   
   
   // Efficiencies
   
   float SelEffTruthErr = errf(nSelSumTruth,sqrt(nSelSumTruth),nEventsSumTruth,sqrt(nEventsSumTruth));
   float SelEffHighestCSVv2Err = errf(nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2),nEventsSumHighestCSVv2,sqrt(nEventsSumHighestCSVv2));
   float SelEffCSVv2LErr = errf(nSelSumCSVv2L,sqrt(nSelSumCSVv2L),nEventsSumCSVv2L,sqrt(nEventsSumCSVv2L));
   float SelEffCSVv2MErr = errf(nSelSumCSVv2M,sqrt(nSelSumCSVv2M),nEventsSumCSVv2M,sqrt(nEventsSumCSVv2M));
   float SelEffCSVv2TErr = errf(nSelSumCSVv2T,sqrt(nSelSumCSVv2T),nEventsSumCSVv2T,sqrt(nEventsSumCSVv2T));
   float SelEffAllErr = errf(nSelSumAll,sqrt(nSelSumAll),nEventsSumAll,sqrt(nEventsSumAll));
   
   float AlgEffTruthErr = errf(nMatchSumTruth,sqrt(nMatchSumTruth),nSelSumTruth,sqrt(nSelSumTruth));
   float AlgEffHighestCSVv2Err = errf(nMatchSumHighestCSVv2,sqrt(nMatchSumHighestCSVv2),nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2));
   float AlgEffCSVv2LErr = errf(nMatchSumCSVv2L,sqrt(nMatchSumCSVv2L),nSelSumCSVv2L,sqrt(nSelSumCSVv2L));
   float AlgEffCSVv2MErr = errf(nMatchSumCSVv2M,sqrt(nMatchSumCSVv2M),nSelSumCSVv2M,sqrt(nSelSumCSVv2M));
   float AlgEffCSVv2TErr = errf(nMatchSumCSVv2T,sqrt(nMatchSumCSVv2T),nSelSumCSVv2T,sqrt(nSelSumCSVv2T));
   float AlgEffAllErr = errf(nMatchSumAll,sqrt(nMatchSumAll),nSelSumAll,sqrt(nSelSumAll));

   float AlgBJetEffTruthErr = errf(nMatchBJetSumTruth,sqrt(nMatchBJetSumTruth),nSelSumTruth,sqrt(nSelSumTruth));
   float AlgBJetEffHighestCSVv2Err = errf(nMatchBJetSumHighestCSVv2,sqrt(nMatchBJetSumHighestCSVv2),nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2));
   float AlgBJetEffCSVv2LErr = errf(nMatchBJetSumCSVv2L,sqrt(nMatchBJetSumCSVv2L),nSelSumCSVv2L,sqrt(nSelSumCSVv2L));
   float AlgBJetEffCSVv2MErr = errf(nMatchBJetSumCSVv2M,sqrt(nMatchBJetSumCSVv2M),nSelSumCSVv2M,sqrt(nSelSumCSVv2M));
   float AlgBJetEffCSVv2TErr = errf(nMatchBJetSumCSVv2T,sqrt(nMatchBJetSumCSVv2T),nSelSumCSVv2T,sqrt(nSelSumCSVv2T));
   float AlgBJetEffAllErr = errf(nMatchBJetSumAll,sqrt(nMatchBJetSumAll),nSelSumAll,sqrt(nSelSumAll));
   
   float NoSolutionEffTruthErr = errf(nNoSolutionSumTruth,sqrt(nNoSolutionSumTruth),nSelSumTruth,sqrt(nSelSumTruth));
   float NoSolutionEffHighestCSVv2Err = errf(nNoSolutionSumHighestCSVv2,sqrt(nNoSolutionSumHighestCSVv2),nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2));
   float NoSolutionEffCSVv2LErr = errf(nNoSolutionSumCSVv2L,sqrt(nNoSolutionSumCSVv2L),nSelSumCSVv2L,sqrt(nSelSumCSVv2L));
   float NoSolutionEffCSVv2MErr = errf(nNoSolutionSumCSVv2M,sqrt(nNoSolutionSumCSVv2M),nSelSumCSVv2M,sqrt(nSelSumCSVv2M));
   float NoSolutionEffCSVv2TErr = errf(nNoSolutionSumCSVv2T,sqrt(nNoSolutionSumCSVv2T),nSelSumCSVv2T,sqrt(nSelSumCSVv2T));
   float NoSolutionEffAllErr = errf(nNoSolutionSumAll,sqrt(nNoSolutionSumAll),nSelSumAll,sqrt(nSelSumAll));
   
   float MVAEffTruthErr = errf(nMatchMVASumTruth,sqrt(nMatchMVASumTruth),nSelSumTruth,sqrt(nSelSumTruth));
   float MVAEffHighestCSVv2Err = errf(nMatchMVASumHighestCSVv2,sqrt(nMatchMVASumHighestCSVv2),nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2));
   float MVAEffCSVv2LErr = errf(nMatchMVASumCSVv2L,sqrt(nMatchMVASumCSVv2L),nSelSumCSVv2L,sqrt(nSelSumCSVv2L));
   float MVAEffCSVv2MErr = errf(nMatchMVASumCSVv2M,sqrt(nMatchMVASumCSVv2M),nSelSumCSVv2M,sqrt(nSelSumCSVv2M));
   float MVAEffCSVv2TErr = errf(nMatchMVASumCSVv2T,sqrt(nMatchMVASumCSVv2T),nSelSumCSVv2T,sqrt(nSelSumCSVv2T));
   float MVAEffAllErr = errf(nMatchMVASumAll,sqrt(nMatchMVASumAll),nSelSumAll,sqrt(nSelSumAll));

   float MVABJetEffTruthErr = errf(nMatchBJetMVASumTruth,sqrt(nMatchBJetMVASumTruth),nSelSumTruth,sqrt(nSelSumTruth));
   float MVABJetEffHighestCSVv2Err = errf(nMatchBJetMVASumHighestCSVv2,sqrt(nMatchBJetMVASumHighestCSVv2),nSelSumHighestCSVv2,sqrt(nSelSumHighestCSVv2));
   float MVABJetEffCSVv2LErr = errf(nMatchBJetMVASumCSVv2L,sqrt(nMatchBJetMVASumCSVv2L),nSelSumCSVv2L,sqrt(nSelSumCSVv2L));
   float MVABJetEffCSVv2MErr = errf(nMatchBJetMVASumCSVv2M,sqrt(nMatchBJetMVASumCSVv2M),nSelSumCSVv2M,sqrt(nSelSumCSVv2M));
   float MVABJetEffCSVv2TErr = errf(nMatchBJetMVASumCSVv2T,sqrt(nMatchBJetMVASumCSVv2T),nSelSumCSVv2T,sqrt(nSelSumCSVv2T));
   float MVABJetEffAllErr = errf(nMatchBJetMVASumAll,sqrt(nMatchBJetMVASumAll),nSelSumAll,sqrt(nSelSumAll));

   float SelEffTruth = float(nSelSumTruth)/float(nEventsSumTruth);
   float SelEffAll = float(nSelSumAll)/float(nEventsSumAll);
   float SelEffHighestCSVv2 = float(nSelSumHighestCSVv2)/float(nEventsSumHighestCSVv2);
   float SelEffCSVv2L = float(nSelSumCSVv2L)/float(nEventsSumCSVv2L);
   float SelEffCSVv2M = float(nSelSumCSVv2M)/float(nEventsSumCSVv2M);
   float SelEffCSVv2T = float(nSelSumCSVv2T)/float(nEventsSumCSVv2T);
   
   float AlgEffTruth = float(nMatchSumTruth)/float(nSelSumTruth);
   float AlgEffAll = float(nMatchSumAll)/float(nSelSumAll);
   float AlgEffHighestCSVv2 = float(nMatchSumHighestCSVv2)/float(nSelSumHighestCSVv2);
   float AlgEffCSVv2L = float(nMatchSumCSVv2L)/float(nSelSumCSVv2L);
   float AlgEffCSVv2M = float(nMatchSumCSVv2M)/float(nSelSumCSVv2M);
   float AlgEffCSVv2T = float(nMatchSumCSVv2T)/float(nSelSumCSVv2T);
   
   float MVAEffTruth = float(nMatchMVASumTruth)/float(nSelSumTruth);
   float MVAEffAll = float(nMatchMVASumAll)/float(nSelSumAll);
   float MVAEffHighestCSVv2 = float(nMatchMVASumHighestCSVv2)/float(nSelSumHighestCSVv2);
   float MVAEffCSVv2L = float(nMatchMVASumCSVv2L)/float(nSelSumCSVv2L);
   float MVAEffCSVv2M = float(nMatchMVASumCSVv2M)/float(nSelSumCSVv2M);
   float MVAEffCSVv2T = float(nMatchMVASumCSVv2T)/float(nSelSumCSVv2T);

   float AlgBJetEffTruth = float(nMatchBJetSumTruth)/float(nSelSumTruth);
   float AlgBJetEffAll = float(nMatchBJetSumAll)/float(nSelSumAll);
   float AlgBJetEffHighestCSVv2 = float(nMatchBJetSumHighestCSVv2)/float(nSelSumHighestCSVv2);
   float AlgBJetEffCSVv2L = float(nMatchBJetSumCSVv2L)/float(nSelSumCSVv2L);
   float AlgBJetEffCSVv2M = float(nMatchBJetSumCSVv2M)/float(nSelSumCSVv2M);
   float AlgBJetEffCSVv2T = float(nMatchBJetSumCSVv2T)/float(nSelSumCSVv2T);

   float MVABJetEffTruth = float(nMatchBJetMVASumTruth)/float(nSelSumTruth);
   float MVABJetEffAll = float(nMatchBJetMVASumAll)/float(nSelSumAll);
   float MVABJetEffHighestCSVv2 = float(nMatchBJetMVASumHighestCSVv2)/float(nSelSumHighestCSVv2);
   float MVABJetEffCSVv2L = float(nMatchBJetMVASumCSVv2L)/float(nSelSumCSVv2L);
   float MVABJetEffCSVv2M = float(nMatchBJetMVASumCSVv2M)/float(nSelSumCSVv2M);
   float MVABJetEffCSVv2T = float(nMatchBJetMVASumCSVv2T)/float(nSelSumCSVv2T);

   float NoSolutionEffTruth = float(nNoSolutionSumTruth)/float(nSelSumTruth);
   float NoSolutionEffAll = float(nNoSolutionSumAll)/float(nSelSumAll);
   float NoSolutionEffHighestCSVv2 = float(nNoSolutionSumHighestCSVv2)/float(nSelSumHighestCSVv2);
   float NoSolutionEffCSVv2L = float(nNoSolutionSumCSVv2L)/float(nSelSumCSVv2L);
   float NoSolutionEffCSVv2M = float(nNoSolutionSumCSVv2M)/float(nSelSumCSVv2M);
   float NoSolutionEffCSVv2T = float(nNoSolutionSumCSVv2T)/float(nSelSumCSVv2T);
   
   std::cout << "Selection Efficiency (Truth) = " << SelEffTruth*100 << " +- " << SelEffTruthErr*100 << " %" << std::endl;
   std::cout << "Selection Efficiency (All) = " << SelEffAll*100 << " +- " << SelEffAllErr*100 << " %" << std::endl;
   std::cout << "Selection Efficiency (HighestCSVv2) = " << SelEffHighestCSVv2*100 << " +- " << SelEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "Selection Efficiency (CSVv2L) = " << SelEffCSVv2L*100 << " +- " << SelEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "Selection Efficiency (CSVv2M) = " << SelEffCSVv2M*100 << " +- " << SelEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "Selection Efficiency (CSVv2T) = " << SelEffCSVv2T*100 << " +- " << SelEffCSVv2TErr*100 << " %" << std::endl;
   
   std::cout << "Algorithm Efficiency for all jets (Truth) = " << AlgEffTruth*100 << " +- " << AlgEffTruthErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for all jets (All) = " << AlgEffAll*100 << " +- " << AlgEffAllErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for all jets (HighestCSVv2) = " << AlgEffHighestCSVv2*100 << " +- " << AlgEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for all jets (CSVv2L) = " << AlgEffCSVv2L*100 << " +- " << AlgEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for all jets (CSVv2M) = " << AlgEffCSVv2M*100 << " +- " << AlgEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for all jets (CSVv2T) = " << AlgEffCSVv2T*100 << " +- " << AlgEffCSVv2TErr*100 << " %" << std::endl;
   
   std::cout << "Algorithm+MVA Efficiency for all jets (Truth) = " << MVAEffTruth*100 << " +- " << MVAEffTruthErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for all jets (All) = " << MVAEffAll*100 << " +- " << MVAEffAllErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for all jets (HighestCSVv2) = " << MVAEffHighestCSVv2*100 << " +- " << MVAEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2L) = " << MVAEffCSVv2L*100 << " +- " << MVAEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2M) = " << MVAEffCSVv2M*100 << " +- " << MVAEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for all jets (CSVv2T) = " << MVAEffCSVv2T*100 << " +- " << MVAEffCSVv2TErr*100 << " %" << std::endl;

   std::cout << "Algorithm Efficiency for b jets (Truth) = " << AlgBJetEffTruth*100 << " +- " << AlgBJetEffTruthErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for b jets (All) = " << AlgBJetEffAll*100 << " +- " << AlgBJetEffAllErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for b jets (HighestCSVv2) = " << AlgBJetEffHighestCSVv2*100 << " +- " << AlgBJetEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for b jets (CSVv2L) = " << AlgBJetEffCSVv2L*100 << " +- " << AlgBJetEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for b jets (CSVv2M) = " << AlgBJetEffCSVv2M*100 << " +- " << AlgBJetEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "Algorithm Efficiency for b jets (CSVv2T) = " << AlgBJetEffCSVv2T*100 << " +- " << AlgBJetEffCSVv2TErr*100 << " %" << std::endl;
   
   std::cout << "Algorithm+MVA Efficiency for b jets (Truth) = " << MVABJetEffTruth*100 << " +- " << MVABJetEffTruthErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for b jets (All) = " << MVABJetEffAll*100 << " +- " << MVABJetEffAllErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for b jets (HighestCSVv2) = " << MVABJetEffHighestCSVv2*100 << " +- " << MVABJetEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2L) = " << MVABJetEffCSVv2L*100 << " +- " << MVABJetEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2M) = " << MVABJetEffCSVv2M*100 << " +- " << MVABJetEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "Algorithm+MVA Efficiency for b jets (CSVv2T) = " << MVABJetEffCSVv2T*100 << " +- " << MVABJetEffCSVv2TErr*100 << " %" << std::endl;
   
   std::cout << "No solutions found (Truth) = " << NoSolutionEffTruth*100 << " +- " << NoSolutionEffTruthErr*100 << " %" << std::endl;
   std::cout << "No solutions found (All) = " << NoSolutionEffAll*100 << " +- " << NoSolutionEffAllErr*100 << " %" << std::endl;
   std::cout << "No solutions found (HighestCSVv2) = " << NoSolutionEffHighestCSVv2*100 << " +- " << NoSolutionEffHighestCSVv2Err*100 << " %" << std::endl;
   std::cout << "No solutions found (CSVv2L) = " << NoSolutionEffCSVv2L*100 << " +- " << NoSolutionEffCSVv2LErr*100 << " %" << std::endl;
   std::cout << "No solutions found (CSVv2M) = " << NoSolutionEffCSVv2M*100 << " +- " << NoSolutionEffCSVv2MErr*100 << " %" << std::endl;
   std::cout << "No solutions found (CSVv2T) = " << NoSolutionEffCSVv2T*100 << " +- " << NoSolutionEffCSVv2TErr*100 << " %" << std::endl;
   
   TCanvas *cEff = new TCanvas("cEff","cEff",800,400);
   cEff->Draw();
   cEff->cd();

   TH1F *hChartSel = new TH1F("hChartSel","hChartSel",3,0,3);
   hChartSel->SetFillColor(46);
   hChartSel->SetBarWidth(0.8);
   hChartSel->SetBarOffset(0.1);
   hChartSel->SetBinContent(3,SelEffCSVv2L*100.);
   hChartSel->SetBinContent(2,SelEffCSVv2M*100.);
   hChartSel->SetBinContent(1,SelEffCSVv2T*100.);
   hChartSel->SetMaximum(105.);
   hChartSel->SetMinimum(0.);
   hChartSel->GetXaxis()->SetBinLabel(3,"CSVv2L");
   hChartSel->GetXaxis()->SetBinLabel(2,"CSVv2M");
   hChartSel->GetXaxis()->SetBinLabel(1,"CSVv2T");
   hChartSel->GetYaxis()->SetTitle("Selection efficiency [%]");   
   hChartSel->Draw("HBAR");
   cEff->Print("pics/ChartSel.eps");
   cEff->Clear();

   TH1F *hChartAlg = new TH1F("hChartAlg","hChartAlg",6,0,6);
   hChartAlg->SetFillColor(38);
   hChartAlg->SetBarWidth(0.8);
   hChartAlg->SetBarOffset(0.1);
   hChartAlg->SetBinContent(6,AlgEffTruth*100.);
   hChartAlg->SetBinContent(5,AlgEffAll*100.);
   hChartAlg->SetBinContent(4,AlgEffHighestCSVv2*100.);
   hChartAlg->SetBinContent(3,AlgEffCSVv2L*100.);
   hChartAlg->SetBinContent(2,AlgEffCSVv2M*100.);
   hChartAlg->SetBinContent(1,AlgEffCSVv2T*100.);
   hChartAlg->SetMaximum(105.);
   hChartAlg->SetMinimum(0.);
   hChartAlg->GetXaxis()->SetBinLabel(6,"Truth");
   hChartAlg->GetXaxis()->SetBinLabel(5,"All");
   hChartAlg->GetXaxis()->SetBinLabel(4,"Highest CSVv2");
   hChartAlg->GetXaxis()->SetBinLabel(3,"CSVv2L");
   hChartAlg->GetXaxis()->SetBinLabel(2,"CSVv2M");
   hChartAlg->GetXaxis()->SetBinLabel(1,"CSVv2T");
   hChartAlg->GetYaxis()->SetTitle("Algorithm efficiency [%]");   
   hChartAlg->Draw("HBAR");
   cEff->Print("pics/ChartAlg.eps");
   cEff->Clear();

   TH1F *hChartMVA = new TH1F("hChartMVA","hChartMVA",6,0,6);
   hChartMVA->SetFillColor(38);
   hChartMVA->SetBarWidth(0.8);
   hChartMVA->SetBarOffset(0.1);
   hChartMVA->SetBinContent(6,MVAEffTruth*100.);
   hChartMVA->SetBinContent(5,MVAEffAll*100.);
   hChartMVA->SetBinContent(4,MVAEffHighestCSVv2*100.);
   hChartMVA->SetBinContent(3,MVAEffCSVv2L*100.);
   hChartMVA->SetBinContent(2,MVAEffCSVv2M*100.);
   hChartMVA->SetBinContent(1,MVAEffCSVv2T*100.);
   hChartMVA->SetMaximum(105.);
   hChartMVA->SetMinimum(0.);
   hChartMVA->GetXaxis()->SetBinLabel(6,"Truth");
   hChartMVA->GetXaxis()->SetBinLabel(5,"All");
   hChartMVA->GetXaxis()->SetBinLabel(4,"Highest CSVv2");
   hChartMVA->GetXaxis()->SetBinLabel(3,"CSVv2L");
   hChartMVA->GetXaxis()->SetBinLabel(2,"CSVv2M");
   hChartMVA->GetXaxis()->SetBinLabel(1,"CSVv2T");
   hChartMVA->GetYaxis()->SetTitle("Algorithm efficiency [%]");   
   hChartMVA->Draw("HBAR");
   cEff->Print("pics/ChartMVA.eps");
   cEff->Clear();

   TH1F *hChartAlgBJet = new TH1F("hChartAlgBJet","hChartAlgBJet",6,0,6);
   hChartAlgBJet->SetFillColor(38);
   hChartAlgBJet->SetBarWidth(0.8);
   hChartAlgBJet->SetBarOffset(0.1);
   hChartAlgBJet->SetBinContent(6,AlgBJetEffTruth*100.);
   hChartAlgBJet->SetBinContent(5,AlgBJetEffAll*100.);
   hChartAlgBJet->SetBinContent(4,AlgBJetEffHighestCSVv2*100.);
   hChartAlgBJet->SetBinContent(3,AlgBJetEffCSVv2L*100.);
   hChartAlgBJet->SetBinContent(2,AlgBJetEffCSVv2M*100.);
   hChartAlgBJet->SetBinContent(1,AlgBJetEffCSVv2T*100.);
   hChartAlgBJet->SetMaximum(105.);
   hChartAlgBJet->SetMinimum(0.);
   hChartAlgBJet->GetXaxis()->SetBinLabel(6,"Truth");
   hChartAlgBJet->GetXaxis()->SetBinLabel(5,"All");
   hChartAlgBJet->GetXaxis()->SetBinLabel(4,"Highest CSVv2");
   hChartAlgBJet->GetXaxis()->SetBinLabel(3,"CSVv2L");
   hChartAlgBJet->GetXaxis()->SetBinLabel(2,"CSVv2M");
   hChartAlgBJet->GetXaxis()->SetBinLabel(1,"CSVv2T");
   hChartAlgBJet->GetYaxis()->SetTitle("Algorithm efficiency [%]");   
   hChartAlgBJet->Draw("HBAR");
   cEff->Print("pics/ChartAlgBJet.eps");
   cEff->Clear();

   TH1F *hChartMVABJet = new TH1F("hChartMVABJet","hChartMVABJet",6,0,6);
   hChartMVABJet->SetFillColor(38);
   hChartMVABJet->SetBarWidth(0.8);
   hChartMVABJet->SetBarOffset(0.1);
   hChartMVABJet->SetBinContent(6,MVABJetEffTruth*100.);
   hChartMVABJet->SetBinContent(5,MVABJetEffAll*100.);
   hChartMVABJet->SetBinContent(4,MVABJetEffHighestCSVv2*100.);
   hChartMVABJet->SetBinContent(3,MVABJetEffCSVv2L*100.);
   hChartMVABJet->SetBinContent(2,MVABJetEffCSVv2M*100.);
   hChartMVABJet->SetBinContent(1,MVABJetEffCSVv2T*100.);
   hChartMVABJet->SetMaximum(105.);
   hChartMVABJet->SetMinimum(0.);
   hChartMVABJet->GetXaxis()->SetBinLabel(6,"Truth");
   hChartMVABJet->GetXaxis()->SetBinLabel(5,"All");
   hChartMVABJet->GetXaxis()->SetBinLabel(4,"Highest CSVv2");
   hChartMVABJet->GetXaxis()->SetBinLabel(3,"CSVv2L");
   hChartMVABJet->GetXaxis()->SetBinLabel(2,"CSVv2M");
   hChartMVABJet->GetXaxis()->SetBinLabel(1,"CSVv2T");
   hChartMVABJet->GetYaxis()->SetTitle("Algorithm efficiency [%]");   
   hChartMVABJet->Draw("HBAR");
   cEff->Print("pics/ChartMVABJet.eps");
   cEff->Clear();
   
   gApplication->Terminate();
}

void addbin(TH1D *h)
{   
   // Add overflow and underflow bins
   Int_t x_nbins = h->GetXaxis()->GetNbins();
   h->SetBinContent(1,h->GetBinContent(0)+h->GetBinContent(1));
   h->SetBinError(1,TMath::Sqrt(pow(h->GetBinError(0),2)+pow(h->GetBinError(1),2)));
   h->SetBinContent(x_nbins,h->GetBinContent(x_nbins)+h->GetBinContent(x_nbins+1));
   h->SetBinError(x_nbins,TMath::Sqrt(pow(h->GetBinError(x_nbins),2)+
				      pow(h->GetBinError(x_nbins+1),2)));
   // Set overflow and underflow bins to 0
   h->SetBinContent(0,0.);
   h->SetBinError(0,0.);
   h->SetBinContent(x_nbins+1,0.);
   h->SetBinError(x_nbins+1,0.);
}

void tdrGrid(bool gridOn) {
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(450); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.17);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.06);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.00);
  tdrStyle->SetTitleYOffset(1.00);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

float errf(float v1,float ve1,float v2,float ve2)
{
   if( v2 == 0 ) return -666;

   float err = ve1*ve1/v2/v2 + v1*v1*ve2*ve2/v2/v2/v2/v2;

   err = sqrt(err);

   return err;
}
