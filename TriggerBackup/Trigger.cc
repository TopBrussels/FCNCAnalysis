
#include "../interface/Trigger.h"

Trigger::Trigger(bool isMuon, bool isElectron, bool trigSingleLep, bool trigDoubleLep, bool trigTrilep):
muon(false), electron(false), singleLep(false), doubleLep(false), trilep(false), fullHadr(false), trigged(false), redotrigmap(false), triggerList(), currentRunTrig(0), previousRunTrig(-1), currentFilenameTrig(""), previousFilenameTrig(""), iFileTrig(-1), treeNumberTrig(-1), triggermap()
{
  if (isMuon)
  {
    muon = true;
  }
  if (isElectron)
  {
    electron = true;
  }
  else if(!isMuon && !isElectron)
  {
    cout << "TRIGGER::TRIGGER - No selected lepton..." << endl;
  }
  if (trigSingleLep)
  {
    singleLep = true;
  }
  if (trigDoubleLep)
  {
    doubleLep = true;
  }
  if(trigTrilep){
    trilep = true;
  }
  if (! trigSingleLep && ! trigDoubleLep && !trigTrilep)
  {
    fullHadr = true;
    cout << "TRIGGER::TRIGGER - Do you really want no lepton triggers?" << endl;
  }
}

Trigger::~Trigger()
{
  
}

void Trigger::bookTriggers(bool isData, string dName)
{
  /// This function is called in the dataset loop
  //  Reset all quantities here
  triggerList.clear();
  currentRunTrig = 0;
  previousRunTrig = -1;
  currentFilenameTrig = previousFilenameTrig = "";
  iFileTrig = -1;
  treeNumberTrig = -1;
  triggermap.clear();
  
  /// Add relevant triggers to triggerlist
  //  Recommended triggers for TOP analyses
  //  Last updated: 6 dec 2016. https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger
  if (singleLep && !doubleLep && !trilep)
  {
    if (muon)
    {
      if (isData)
      {
        if( dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)
        {
          triggerList.push_back("HLT_IsoMu24_v*");
          triggerList.push_back("HLT_IsoTkMu24_v*");
        }
        /*else if(  dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)
        {
         // triggerList.push_back("HLT_Mu50_v*");
         // triggerList.push_back("HLT_TkMu50_v*");
        }*/
    /*else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
          triggerList.push_back("HLT_IsoMu18_v*");  // Data
          //           triggerList.push_back("HLT_IsoTkMu2_v*");
        }
*/        else cout << "ERROR: no correct data triggers selected ins single mu" << endl;
      /*  if( dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)
        {
          triggerList.push_back("HLT_IsoMu22_eta2p1_v*");  // prescaled!
          triggerList.push_back("HLT_IsoTkMu22_eta2p1_v*"); // prescaled!
        }*/
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
      //    triggerList.push_back("HLT_IsoMu22_eta2p1_v*");
       //   triggerList.push_back("HLT_IsoTkMu22_eta2p1_v*");
          triggerList.push_back("HLT_IsoMu24_v2");
          triggerList.push_back("HLT_IsoTkMu24_v2");
        }
        else if( dName.find("76X")!=string::npos ){
          triggerList.push_back("HLT_IsoMu18_v*");
          //           triggerList.push_back("HLT_IsoTkMu20_v*");
        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
    if (electron)
    {
      if (isData)
      {
        
        if( dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)        {
          triggerList.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v*");
          //triggerList.push_back("HLT_Ele27_WPTight_Gsf_v*"); // relative unprescaled
          //triggerList.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v*"); // relative unprescaled
        }
     /*   else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
          triggerList.push_back("HLT_Ele23_WPLoose_Gsf_v*");  // Data, restricted to eta < 2.1
        }*/
        else cout << "ERROR: no correct data triggers selected in single ele" << endl;
        
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
          triggerList.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v*");
       //   triggerList.push_back("HLT_Ele27_WPTight_Gsf_v*");
       //   triggerList.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v*");
        }
        else if( dName.find("76X")!=string::npos){
          triggerList.push_back("HLT_Ele23_WPLoose_Gsf_v*");
        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
  }
  
  if (doubleLep && !singleLep && !trilep)   // Triggers for data & MC are the same
  {
    //cout << "muon " << muon << " electron " << electron << endl;
    if (muon && !electron)
    {
      if (isData)
      {
        
        if( dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos )
        {
          //cout << "pushing back " << endl;
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*"); // prescaled for H
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"); // prescaked for H
        }
        else if( dName.find("Run_2016H")!=string::npos)
        {
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        }
        else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");  // Data, restricted to eta < 2.1
        }
        else cout << "ERROR: no correct data triggers selected in dimu" << endl;
        
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        }
        else if( dName.find("76X")!=string::npos ){
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
    else if (electron && !muon)
    {
      if (isData)
      {
        
        if( dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)
        {
         // triggerList.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          triggerList.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*");
        }
        else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
          triggerList.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");  // Data, restricted to eta < 2.1
        }
        else cout << "ERROR: no correct data triggers selected in dilep" << endl;
        
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
         // triggerList.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          triggerList.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*");
        }
        else if( dName.find("76X")!=string::npos) {
          triggerList.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
    
    else if (muon && electron)
    {
      if (isData)
      {
        
        if(dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos )
        {
          //cout << "pushing back " << endl;
          triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
             triggerList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*" );
          triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
          
          triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          

          
        }
          else if(dName.find("Run_2016H")!=string::npos)
        {
          //triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
          // triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
          triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          // triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          //triggerList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*" );
          //triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
          triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          
        }
        else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        }
        else cout << "ERROR: no correct data triggers selected in muon + electron " << endl;
        
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
          triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
          triggerList.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          triggerList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
          triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
          triggerList.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
          //triggerList.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4");
        }
        else if( dName.find("76X")!=string::npos ){
          triggerList.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
    
  }
  if(trilep){
    if(muon && !electron)
    {
      
      if (isData)
      {
        
        if( dName.find("Run_2016B")!=string::npos || dName.find("Run_2016C")!=string::npos || dName.find("Run_2016D")!=string::npos || dName.find("Run_2016E")!=string::npos || dName.find("Run_2016F")!=string::npos || dName.find("Run_2016G")!=string::npos || dName.find("Run_2016H")!=string::npos)
        {
          triggerList.push_back("HLT_TripleMu_12_10_5_v*");
          triggerList.push_back("HLT_TripleMu_5_3_3_v*");
          
        }
        else if( dName.find("Run_2015C")!=string::npos || dName.find("Run_2015D")!=string::npos ){
           triggerList.push_back("HLT_TripleMu_12_10_5_v*");
          
        
        }
        else cout << "ERROR: no correct data triggers selected in trilep muons" << endl;
        
      }
      else
      {
        if( dName.find("80X")!=string::npos)
        {
          triggerList.push_back("HLT_TripleMu_12_10_5_v*");
          triggerList.push_back("HLT_TripleMu_5_3_3_v*");
        }
        else if( dName.find("76X")!=string::npos ){
          triggerList.push_back("HLT_TripleMu_12_10_5_v*");

        }
        else cout << "ERROR: no correct MC triggers selected" << endl;
      }
    }
    

    else if(!muon && electron)
    {
      triggerList.push_back("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*");
      
    }
    else if(muon && electron)
    {
      triggerList.push_back("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*");
      triggerList.push_back("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*");
      
    }   



 
  }
  for(UInt_t iTrig = 0; iTrig < triggerList.size(); iTrig++)
  {
    triggermap[triggerList[iTrig]] = std::pair<int,bool> (-999,false);
  }
  
}

void Trigger::checkAvail(int currentRunTrig, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, bool verbose)
{
  redotrigmap = false;
  treeNumberTrig = datasets[d]->eventTree()->GetTreeNumber();
  currentFilenameTrig = datasets[d]->eventTree()->GetFile()->GetName();
  if (previousFilenameTrig != currentFilenameTrig)
  {
    previousFilenameTrig = currentFilenameTrig;
    iFileTrig++;
    redotrigmap = true;
    cout << endl << "*****File changed!!! => iFile = " << iFileTrig << " new file is " << currentFilenameTrig << " in sample " << datasets[d]->Name() << " *****" << endl;
  }
  if (previousRunTrig != currentRunTrig)
  {
    previousRunTrig = currentRunTrig;
    cout << "*****Run changed!!! => new run = " << previousRunTrig << " *****" << endl;
    redotrigmap=true;
  }
  
  if (verbose && redotrigmap)
  {
    treeLoader->ListTriggers(currentRunTrig, treeNumberTrig);
  }
  
  
  // get trigger info:
  for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermap.begin(); iter != triggermap.end(); iter++)
  {
    if (redotrigmap)
    {
      Int_t loc = treeLoader->iTrigger(iter->first, currentRunTrig, treeNumberTrig);
      string trigname = iter->first;
      cout << "trigname: " << trigname << "  location: " << loc << /*"  " << currentFilenameTrig << "  " << currentRunTrig <<*/ endl;
      
      iter->second.first = loc;
    }
    // and check if it exists and if it fired:
    if (iter->second.first >= 0 && iter->second.first != 9999) // trigger exists
      iter->second.second = event->trigHLT(iter->second.first);
    else
      iter->second.second = false;
  }
}


int Trigger::checkIfFired()
{
  // now check if the appropriate triggers fired for each analysis:
  trigged = 0;
  
  for(UInt_t itrig = 0; itrig < triggerList.size() && trigged == 0; itrig++)
  {
    if (triggermap[triggerList[itrig]].second)
      trigged = 1;
  }
  
  return trigged;
}
