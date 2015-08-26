#ifndef __cutFlow_H_INCLUDED__
#define __cutFlow_H_INCLUDED__

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <TFile.h>
#include <TH1.h>

#include <TCanvas.h>
#include <vector>
#include "TGraph.h"
#include "TVectorT.h"
#include <iostream>
#include <iomanip>


class cutFlow
{
public:
		std::vector <std::string>  SampleTxtReader(const char *inputtxt);
		bool RootFile(std::string inputtxt);
};


#endif 
