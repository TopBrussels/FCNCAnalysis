/// isis 
#include <iostream>
#include "TSystem.h"



void chain(int nsel = 0, int mode = 0, bool silent = false){  

  gSystem->CompileMacro("analysis.C","k");

 
  // samples used
  double x_sec = 0.;
  char plotName[300];
  
  
  if (nsel == 0)                	{sprintf(plotName,"TTJetsTocHbW_HToZZ_ZToLL_HctL_tree.root");}
  
  
  
  char myRootFile[300];
  
   sprintf(myRootFile,"../../data/ntuples/%d_%s_tree.root", mode, plotName);
  
  
  
  
  TChain *myCh = new TChain("myTree","myTree");
  myCh->Add(myRootFile);

  analysis an(myCh);
  std::cout << plotName << " (" << myRootFile << ",  " << myCh->GetEntries() << " entries )" << std::endl;
  an.myAnalysis(nsel,mode,silent);
  
}
