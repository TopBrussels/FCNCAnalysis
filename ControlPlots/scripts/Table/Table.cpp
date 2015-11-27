#include "Table.h"

using namespace std;

vector <string>  Table::SampleTxtReader(const char *inputtxt) //Class to read in the txt files and get all the sample names out of it
{
		vector <string> VectorNames;
	
		ifstream InputFile(inputtxt, ifstream::binary);
    string line;

		//Read txt file line by line
    while(getline(InputFile, line))
		{
	        istringstream in(line);     //make a stream for the line itself
					string a;
					in >> a;
										
					VectorNames.push_back(a);
		}
		
		return VectorNames;
}


bool  Table::RootFile(string inputstring) //Class to check whether the given string is a .root file
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
		Table TabClass;

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
		if(information) cout << " Reading in Signal.txt, Background.txt and Variables.txt and Cuts.txt" << endl;
		const char * Signaltxt = "Signal.txt";
		const char * Backgroundtxt = "Background.txt";
		const char * Variabletxt = "Variables.txt";
		const char * Cutstxt = "Cuts.txt"; 
		vector <string> signalnames;
		vector <string> backgroundnames;
		vector <string> Variables;
		vector <string> Cuts; 
	
		signalnames.clear();
		backgroundnames.clear();
		Variables.clear();
		Cuts.clear();
						
		signalnames = TabClass.SampleTxtReader(Signaltxt);
		backgroundnames = TabClass.SampleTxtReader(Backgroundtxt);
		Variables = TabClass.SampleTxtReader(Variabletxt);
	        Cuts = TabClass.SampleTxtReader(Cutstxt);
		
		if(signalnames.size() == 0)
		{
		 	cout << "No signal defined in Signal.txt ... Exiting macro..." << endl;
			return 1;
		}
		if(backgroundnames.size() == 0)
		{
		 	cout << "No background defined in Background.txt ... Exiting macro..." << endl;
			return 1;
		}
		if(Cuts.size() == 0)
		{
		 	cout << "No cuts defined in Cuts.txt ... Exiting macro..." << endl;
			return 1;
		}

		if(information)
		{
			cout << "Signal samples in "<< signalnames[0] <<": " << endl;
				for(int i = 1; i < signalnames.size(); i++)
				{
					cout << "  - " << signalnames[i] << endl;
				}
			cout << "Background samples in " << backgroundnames[0] << ": " << endl;
				for(int i = 1; i < backgroundnames.size(); i++)
				{
					cout << "  - " << backgroundnames[i] << endl;
				}
			cout << "Cuts: " << endl;
				for(int i = 0; i < Cuts.size(); i++)
				{
					cout << "  - " << Cuts[i] << endl;
				}
		}
		
		
		/////////////////////////////////////////////////////
		// Making the optimal cut plots /////////////////////
		/////////////////////////////////////////////////////
		TFile *infile = 0;

		// creating new stream, a tex file
		ofstream output; 
		output.open("Output_Table/Table.tex");
		
	
	        int nprocess = signalnames.size()+backgroundnames.size(); 
		TString processName[nprocess]; //Makes 
		TH1F *histograms[nprocess];
		int endBin =-5;
		//Read in all signal samples for Variables[iVar] and merge them into 1 histogram
		int iVar = 0;
		for(unsigned int iSignal = 1; iSignal < signalnames.size(); iSignal++)
		{
			
			if(!TabClass.RootFile(signalnames[0])) infile = new TFile((signalnames[0] + "merged_MuMu_" +  signalnames[iSignal] + ".root").c_str(),"read"); //Checks whether background samples is in 1 rootfile, or in a directory containing different rootfiles
			else infile = new TFile(signalnames[0].c_str(),"read");
                       
                           
			TH1F *histo( (TH1F*) infile->Get((Variables[iVar]).c_str()) );
			if(!histo && warnings) cout << "Signal " << signalnames[iSignal] << " does not exist or does not contain variable " << Variables[iVar] << " Check Signal.txt or Variables.txt" <<endl;

			if(histo)
			{
				processName[iSignal] = signalnames[iSignal];
				histograms[iSignal] = (TH1F*) infile->Get((Variables[iVar]).c_str());
				endBin = histo->GetNbinsX();
				
				
			}
		}
		//cout << "signals done "<< endl; 
		
		//Read in all background samples for Variables[iVar] and merge them into 1 histogram
		for(unsigned int iBackgr = 1; iBackgr < backgroundnames.size(); iBackgr++)
		{

			if(!TabClass.RootFile(backgroundnames[0])) infile = new TFile((backgroundnames[0] + "merged_MuMu_"+backgroundnames[iBackgr] + ".root").c_str(),"read");//Checks whether background samples are in 1 rootfile, or in a directory containing different rootfiles
			else infile = new TFile(backgroundnames[0].c_str(),"read");




			TH1F *histo( (TH1F*) infile->Get((Variables[iVar]).c_str()) );
			if(!histo && warnings) cout << "Background " << backgroundnames[iBackgr] << " does not exist or does not contain variable " << Variables[iVar] << " Check Background.txt or Variables.txt" <<endl;


			if(histo)
			{
				processName[iBackgr+signalnames.size()] = backgroundnames[iBackgr];
				
				histograms[signalnames.size()+iBackgr]= (TH1F*) infile->Get((Variables[iVar]).c_str()); 
			}
		
		}
		//cout << "bkg done" << endl; 
		
		
		double Values[nprocess][endBin][3];
		for(int iP = 0; iP < signalnames.size()+backgroundnames.size(); iP++)
		{
		        if(iP != 0 && iP != signalnames.size())
			{
			    for(int iB = 1; iB < endBin; iB++)
			    {
			    	Values[iP][iB][0] = histograms[iP]->GetBinContent(iB);
				Values[iP][iB][1] = 3 ; //precision(histograms[iP]->GetBinError(iB));
				Values[iP][iB][2] = histograms[iP]->GetBinError(iB);	
			    }
			}
		}
		


		// header of the tex file
		output << "\\documentclass[a4paper,8pt]{article}" << endl;
		output << "\\usepackage{geometry}" << endl; 
		output << "\\geometry{legalpaper, landscape, margin=0.1in}" << endl; 
  		output << "\\begin{document}" << endl;
  		output << endl;
  		output << endl;
		
		output << "  \\begin{table}" << endl;
  		output << "  \\begin{center}" << endl;
  		output << "  \\begin{tabular} {|l|" << endl; 

		
		
		// Filling the core of the tex file
		for(int iBin = 0; iBin < Cuts.size(); iBin++)
		{
		   output << "c|" ; 
		}
		output << "}" << endl; 
		output << "\\hline" << endl; 
		output << "& " ;
		
		
		for(int iCut = 0; iCut < Cuts.size()-1; iCut++)
		{
			output << Cuts[iCut] << "&" ; 
		
		}
		
		output << Cuts[Cuts.size()-1] ; 
		output << "\\\\" << endl; 
		
		for (int iProc = 0; iProc < nprocess ; iProc++)
		{
		  if(iProc != 0 && iProc != signalnames.size())
		  {     
		       output << processName[iProc] << " &" ;
		       for(int iB = 1; iB < Cuts.size(); iB++)
		       {
		       
		         output << setprecision(Values[iProc][iB][1]) << Values[iProc][iB][0];
			 output << " $\\pm $ " << setprecision(Values[iProc][iB][1]) << Values[iProc][iB][2];
			 output << "&"; 
			 
			 
		       }
		       output << setprecision(Values[iProc][Cuts.size()][1]) << Values[iProc][Cuts.size()][0]; 
		       output << " $\\pm $ " << setprecision(Values[iProc][Cuts.size()][1]) << Values[iProc][Cuts.size()][2];
		       output << "\\\\" << endl; 
		       
		 }   
		}
		

		// footer of the tex file
		output << "   \\hline " << endl;
  		output << "  \\end{tabular}" << endl;
  		output << "  \\end{center}" << endl;
  		output << "  \\end{table}" << endl;
  		output << endl;
  		output << endl;
		output << "\\end{document}" << endl; 
		
		

		// closing the tex file
		output.close();
		
		//		outputfile->Write();

		if(information) cout << "DONE!!! " << endl;
		
		return 0;
}


