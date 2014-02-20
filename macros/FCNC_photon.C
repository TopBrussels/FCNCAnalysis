#include <iostream>

using namespace std;

void FCNC_photon(){

  AutoLibraryLoader::enable();

  gSystem->CompileMacro("TopAnalyzerLite.cc", "k");


}
