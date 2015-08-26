#include "cutFlow.h"

using namespace std;

vector <string>  cutFlow::SampleTxtReader(const char *inputtxt) //Class to read in the txt files and get all the sample names out of it
{
		vector <string> VectorNames;
	
		ifstream InputFile(inputtxt, ifstream::binary);
    string line;

		//Read txt file line by line
    while(getline(InputFile, line))
		{
	        istringstream in(line);     //make a stream for the line itself
					string a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p;
					in >> a;
					in >> b;
					in >> c;
					in >> d;
					in >> e;
					in >> f;
					in >> g;
					in >> h;
					in >> i;
					in >> j;
					in >> k;
					in >> l;
					in >> m;
					in >> n;
					in >> o;
					in >> p;
					
					if(p != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k+ " " + l+ " " + m+ " " + n+ " " + o + " " + p);
					else if(o != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k+ " " + l+ " " + m+ " " + n+ " " + o);
					else if(n != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k+ " " + l+ " " + m+ " " + n);
					else if(m != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k+ " " + l+ " " + m);
					else if(l != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k + " " + l);
					else if(k != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i+ " " + j+ " " + k);
					else if(j != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i + " " + j);
					else if(i != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h + " " + i);
					else if(h != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f+" " +g+" " +h);
					else if(g != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f + " " +g);
					else if(f != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e+" " +f);
					else if(e != "") VectorNames.push_back(a+" " +b+" " +c+" " +d+" " +e);
					else if(d != "") VectorNames.push_back(a+" " +b+" " +c+" " +d);
					else if(c != "") VectorNames.push_back(a+" " +b+" " +c);
					else if(b != "") VectorNames.push_back(a+" " +b);
					else VectorNames.push_back(a);
					
		}
		
		return VectorNames;
}


bool  cutFlow::RootFile(string inputstring) //Class to check whether the given string is a .root file
{

		bool test;
  	size_t found = inputstring.find(".root",0,5);
  	if (found!=std::string::npos) test = true;
		else test = false;

		
		return test;
}

/*double precision(double error){
  
  int precisionValue;
  double factErr = 0; 
  int iN = 0;
  if (error == 0 || error >= 1) precisionValue = 0;
  else if (error < 1) {
    iN = 0;
    factErr = 0; 
    while (factErr < 1){
      factErr = error*(10*iN);
      iN++;  
    }
    precisionValue = iN-1;
  }
  
  if (factErr > 9.5) precisionValue-=1;
  
  return precisionValue;
  
}

*/

int main(int argc,char *argv[])
{
		cutFlow TabClass;

		bool information = true;
		bool warnings = true;
		

    for(int iarg = 0; iarg < argc && argc>1 ; iarg++) //Read out flags
		{
		      std::string argval=argv[iarg];
		
        	if(argval=="--help" || argval =="--h")
					{
						
						return 0;
					}
					
					//if (argval=="--NoInfo")	information = false;  
					
		}
		/////////////////////////////////////
		// Reading in the sample txt files //
		////////////////////////////////////
		if(information) cout << "\033[40;32m[INFO]\033[37m Reading in SamplesToMergeCutFlow.txt " << endl;
		const char * Sample = "SamplesToMergeCutFlow.txt";
		vector <string> samplenames;
		vector <string> backgroundnames;
		string Variables = "cutFlow_";
	
		samplenames.clear();
		samplenames = TabClass.SampleTxtReader(Sample);
		
		if(samplenames.size() == 0)
		{
		 	cout << "\033[40;31m[ERROR]\033[37m No sample defined in SamplesToMergeCutFlow.txt ... Exiting macro..." << endl;
			return 1;
		}

		if(information)
		{
			cout << "\033[40;32m[INFO]\033[37m samples in "<< samplenames[0] <<": " << endl;
				for(int i = 1; i < samplenames.size(); i++)
				{
					cout << "  - " << samplenames[i] << endl;
				}
		}
		
		
		/////////////////////////////////////////////////////
		// Merging the cutflow plots ///////
		/////////////////////////////////////////////////////
		TFile *infile = 0;
		TFile *outputfile = new TFile("Output_cutFlow/cutFlow.root","RECREATE");
		
		
		
	
	    int nprocess = samplenames.size(); 
		TString processName[nprocess]; //Makes 
		TH1F *histograms[nprocess];
		int endBin =-5;
		//Read in all Variable and merge them into 1 histogram
		TH1F *HistoSignal = 0;
		for(unsigned int iSignal = 1; iSignal < samplenames.size(); iSignal++)
		{
			
//			if(!TabClass.RootFile(samplenames[0])) infile = new TFile((samplenames[0] + samplenames[iSignal]).c_str(),"read"); //Checks whether samples are in 1 rootfile, or in a directory containing different rootfiles
/*			else */infile = new TFile((samplenames[0] + "FCNC_1L3B_Run2_TopTree_Study_" + samplenames[iSignal] + "_Mu.root").c_str(),"read");
                       
                           
			TH1F *histo( (TH1F*) infile->Get(("MultiSamplePlot_cutFlow/cutFlow_"+ samplenames[iSignal]).c_str()) );
			if(!histo && warnings) cout << "\033[40;36m[WARNING]\033[37m Sample " << samplenames[iSignal] << " does not exist. Check .txt-file" <<endl;

				if(histo)
				{
					if(!HistoSignal) //Use only the first histogram to be cloned as the HistoSignal (to which further on other signal histos are added)
					{
						HistoSignal = new TH1F(*histo);
						cout << "First fill of HistoSignal" << endl;
					}
					else
					{
						HistoSignal->Add(histo);
						cout << "Adding histos to HistoSignal" << endl;
					}
				}
			delete histo;
		}

		outputfile->cd();
		HistoSignal->Write();

		outputfile->Write();

		if(information) cout << "\033[40;32m[INFO]\033[37m DONE!!! " << endl;
		
		return 0;
}

