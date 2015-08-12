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


bool  Table::RootFile(string inputstring) //Class to check whether the given string is a .root file
{

	bool test;
  	size_t found = inputstring.find(".root",0,5);
  	if (found!=std::string::npos) test = true;
		else test = false;

		
		return test;
}


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
		if(information) cout << "\033[40;32m[INFO]\033[37m Reading in _bk_Signal.txt, _bk_Background.txt and _bk_Variables.txt and Cuts.txt" << endl;
		const char * Signaltxt = "Samples.txt";
		string Variable = "cutFlow";
		vector <string> samplenames;
			
		samplenames.clear();
		samplenames = TabClass.SampleTxtReader(Signaltxt);
		
		if(samplenames.size() == 0)
		{
		 	cout << "\033[40;31m[ERROR]\033[37m No signal defined in _bk_Signal.txt ... Exiting macro..." << endl;
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
		
		
		////////////////////////////////////////////////////////////////////////////////////
		// Making the .tex file from the cutflow /////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		TFile *infile = 0;
		TFile *outputfile = new TFile("Output_Table/Table.tex","RECREATE");
		ofstream output; 
		output.open("Output_Table/Table.tex");
		
		
		
	
	     int nprocess = samplenames.size();
		TString processName[nprocess]; //Makes 
		TH1F *histograms[nprocess];
		int endBin =-5;
		//Read in all signal samples for Variables[iVar] and merge them into 1 histogram
		int iVar = 0;
		for(unsigned int iSignal = 1; iSignal < samplenames.size(); iSignal++)
		{
			
			if(!TabClass.RootFile(samplenames[0])) infile = new TFile((samplenames[0] + samplenames[iSignal] + ".root").c_str(),"read"); //Checks whether samples are in 1 rootfile, or in a directory containing different rootfiles
			else infile = new TFile(samplenames[0].c_str(),"read");


			TH1F *histo( (TH1F*) infile->Get(("MultiSamplePlot_" + Variable +"/" + Variable + "_" + samplenames[iSignal]).c_str()) );
			if(!histo && warnings) cout << "\033[40;36m[WARNING]\033[37m " << samplenames[iSignal] << " does not exist or does not contain variable " << Variable << " Check _bk_Signal.txt or _bk_Variables.txt" <<endl;

			if(histo)
			{
				processName[iSignal] = samplenames[iSignal];
				histograms[iSignal] = (TH1F*) infile->Get(("MultiSamplePlot_" + Variable +"/" + Variable + "_" + samplenames[iSignal]).c_str());
				endBin = histo->GetNbinsX();
			}
		}
		
		
		double Values[nprocess][endBin][3];
		for(int iP = 0; iP < samplenames.size(); iP++)
		{
		        if(iP != 0 && iP != samplenames.size())
			{
				int CutCounter = 0;
			    for(int iB = 1; iB < endBin; iB++)
			    {
			    	if(CutCounter >= endBin) continue;
			    	cout << " " << endl;
			    	cout << histograms[iP]->GetXaxis()->GetBinLabel(iB) << endl;
				cout << "BinContent: " << histograms[iP]->GetBinContent(iB) << endl;
			    	cout << " " << endl;
			    		Values[iP][CutCounter][0] = histograms[iP]->GetBinContent(iB);
					Values[iP][CutCounter][1] = 3 ; //precision(histograms[iP]->GetBinError(iB));
					Values[iP][CutCounter][2] = histograms[iP]->GetBinError(iB);
					CutCounter++;
				cout << CutCounter << endl;
			    }
			}
		}
		
		
		cout << "\\documentclass[a4paper,8pt]{article}" << endl;
		cout << "\\usepackage[landscape]{geometry}" << endl; 
  		cout << "\\begin{document}" << endl;
  		cout << endl;
  		cout << endl;
		
		cout << "  \\begin{table}" << endl;
  		cout << "  \\begin{center}" << endl;
  		cout << "  \\begin{tabular} {|l|"; 
		
		
		
		for(int iBin = 0; iBin < endBin; iBin++)
		{
		   cout << "c|" ; 
		}
		cout << "}" << endl; 
		cout << "\\hline" << endl; 
		cout << "& " ;
		
		
		for(int iCut = 0; iCut < endBin-1; iCut++)
		{
			cout << "cut " << iCut << "&" ; 
		
		}
		
		cout << "Last cutflow";
		cout << "\\\\" << endl; 
		
		for (int iProc = 0; iProc < nprocess ; iProc++)
		{
		  if(iProc != 0 && iProc != samplenames.size())
		  {     
		       cout << processName[iProc] << " &" ;
		       for(int iB = 0; iB < endBin-1; iB++)
		       {
		       
		         cout << setprecision(Values[iProc][iB][1]) << Values[iProc][iB][0];
			 cout << " $\\pm $ " << setprecision(Values[iProc][iB][1]) << Values[iProc][iB][2];
			 cout << "&"; 
			 
			 
		       }
		       cout << setprecision(Values[iProc][endBin-1][1]) << Values[iProc][endBin-1][0]; 
		       cout << " $\\pm $ " << setprecision(Values[iProc][endBin-1][1]) << Values[iProc][endBin-1][2];
		       cout << "\\\\" << endl; 
		       
		 }   
		}

		cout << "   \\hline " << endl;
  		cout << "  \\end{tabular}" << endl;
  		cout << "  \\end{center}" << endl;
  		cout << "  \\end{table}" << endl;
  		cout << endl;
  		cout << endl;
		cout << "\\end{document}" << endl; 
		

		
		outputfile->Write();

		if(information) cout << "\033[40;32m[INFO]\033[37m DONE!!! " << endl;
		
		return 0;
}

