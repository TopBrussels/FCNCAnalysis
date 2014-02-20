//////////////////////////////////////////////////////////
// This class has been cloned from freyaanalysis.h
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

class analysis {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Double_t        xlWeight;
  Double_t        puweight;
  
  Double_t        rawWeight;
  Double_t        lum;
  Int_t           npu;
  Int_t           nvertex;
  Double_t        metPt;
  Double_t        metPx;
  Double_t        metPy;
  
  
  int idParton_1; 
  int idParton_2; 
  double xParton_1;
  double xParton_2;
  double Q2_event;
  
  std::vector<double>  *ptLepton;
  std::vector<double>  *pxLepton;
  std::vector<double>  *pyLepton;
  std::vector<double>  *pzLepton;
  std::vector<double>  *eLepton;
  std::vector<double>  *qLepton;
   
  std::vector<double>  *ptJet;
  std::vector<double>  *pxJet;
  std::vector<double>  *pyJet;
  std::vector<double>  *pzJet;
  std::vector<double>  *eJet;
  std::vector<double>  *etaJet;
  std::vector<double>  *qJet;
  std::vector<bool>    *Btagjet;
  std::vector<double>  *btJPBJet;
  std::vector<double>  *btBJPBJet;
  std::vector<double>  *btCSVBJet;
  std::vector<double>  *btCSVBmvaJet;
  
  // List of branches
  TBranch        *b_xlWeight;   //!
  TBranch        *b_puweight;   //!
  
  TBranch        *b_rawWeight;   //!
  TBranch        *b_lum;   //!
  TBranch        *b_npu;   //!
  TBranch        *b_nvertex;   //!
  TBranch        *b_metPt;   //!
  TBranch        *b_metPx;   //!
  TBranch        *b_metPy;   //!
  
  TBranch	*b_idParton_1; //!
  TBranch	*b_idParton_2; //!
  TBranch	*b_xParton_1; //!
  TBranch	*b_xParton_2; //!
  TBranch	*b_Q2_event; //!
 
  
  TBranch        *b_ptLepton;   //!
  TBranch        *b_pxLepton;   //!
  TBranch        *b_pyLepton;   //!
  TBranch        *b_pzLepton;   //!
  TBranch        *b_eLepton;   //!
  TBranch        *b_qLepton;   //!
  TBranch        *b_ptJet;   //!
  TBranch        *b_pxJet;   //!
  TBranch        *b_pyJet;   //!
  TBranch        *b_pzJet;   //!
  TBranch        *b_eJet;   //!
  TBranch        *b_etaJet;   //!
  TBranch        *b_qJet;   //!
  TBranch        *b_Btagjet;   //!
  TBranch        *b_btJPBJet;   //!
  TBranch        *b_btBJPBJet;   //!
  TBranch        *b_btCSVBJet;   //!
  TBranch        *b_btCSVBmvaJet;   //!
  
  analysis(TTree *tree=0);
  virtual ~analysis();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  void             myAnalysis(int nsel, int mode, bool silent);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

analysis::analysis(TTree *tree) : fChain(0) 
{
  Init(tree);
}

analysis::~analysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void analysis::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // Set object pointer
  
  ptLepton = 0;
  pxLepton = 0;
  pyLepton = 0;
  pzLepton = 0;
  eLepton = 0;
  qLepton = 0;
  ptJet = 0;
  pxJet = 0;
  pyJet = 0;
  pzJet = 0;
  eJet = 0;
  etaJet = 0;
  qJet = 0;
  Btagjet = 0;
  btJPBJet = 0;
  btBJPBJet = 0;
  btCSVBJet = 0;
  btCSVBmvaJet = 0;
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
 
  
  fChain->SetBranchAddress("xlWeight", &xlWeight, &b_xlWeight);
  fChain->SetBranchAddress("puweight", &puweight, &b_puweight);

  fChain->SetBranchAddress("rawWeight", &rawWeight, &b_rawWeight);
  
  fChain->SetBranchAddress("lum", &lum, &b_lum);
  fChain->SetBranchAddress("npu", &npu, &b_npu);
  fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
  
  fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
  fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
  fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
 
 
 fChain->SetBranchAddress("idParton_1", &idParton_1,&b_idParton_1); 
 fChain->SetBranchAddress("idParton_2", &idParton_2,&b_idParton_2); 
 fChain->SetBranchAddress("xParton_1", &xParton_1,&b_xParton_1);   
 fChain->SetBranchAddress("xParton_2", &xParton_2,&b_xParton_2); 
 fChain->SetBranchAddress("Q2_event", &Q2_event,&b_Q2_event); 
 
 
  fChain->SetBranchAddress("ptLepton", &ptLepton, &b_ptLepton);
  fChain->SetBranchAddress("pxLepton", &pxLepton, &b_pxLepton);
  fChain->SetBranchAddress("pyLepton", &pyLepton, &b_pyLepton);
  fChain->SetBranchAddress("pzLepton", &pzLepton, &b_pzLepton);
  fChain->SetBranchAddress("eLepton", &eLepton, &b_eLepton);
  fChain->SetBranchAddress("qLepton", &qLepton, &b_qLepton);
  
  fChain->SetBranchAddress("ptJet", &ptJet, &b_ptJet);
  fChain->SetBranchAddress("pxJet", &pxJet, &b_pxJet);
  fChain->SetBranchAddress("pyJet", &pyJet, &b_pyJet);
  fChain->SetBranchAddress("pzJet", &pzJet, &b_pzJet);
  fChain->SetBranchAddress("eJet", &eJet, &b_eJet);
  fChain->SetBranchAddress("etaJet", &etaJet, &b_etaJet);
  fChain->SetBranchAddress("qJet", &qJet, &b_qJet);
  fChain->SetBranchAddress("Btagjet", &Btagjet, &b_Btagjet);
  
  fChain->SetBranchAddress("btJPBJet",&btJPBJet, &b_btJPBJet);
  fChain->SetBranchAddress("btBJPBJet",&btBJPBJet, &b_btBJPBJet);
  fChain->SetBranchAddress("btCSVBJet",&btCSVBJet, &b_btCSVBJet);
  fChain->SetBranchAddress("btCSVBmvaJet",&btCSVBmvaJet, &b_btCSVBmvaJet);
  
  Notify();
}

Bool_t analysis::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}

void analysis::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}


#endif
