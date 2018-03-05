#include<sstream>
#include<iostream>
#include<ConfigContainer.hpp>

using namespace std;
int gnModelsPerParam(10000);
int gnMaxPulsarNum(1);





void PrintHelp()
{
	cout<<"DM Maker"<<endl;
	cout<<"Computes DM slices according to config"<<endl;
	cout<<"Operation modes:"<<endl;
	cout<<"-h               - prints this help"<<endl;
	cout<<"-C               - compute slices mode"<<endl;
	cout<<"-J               - joining the slices togheter"<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-o <filename>    - Output filename"<<endl;
	cout<<"-n               - SliceNumber (from 0 to nDMDistSlices-1)"<<endl;
	cout<<"-ip <prefix>     - input prefix"<<endl;
	cout<<"-is <suffix>     - input suffix"<<endl;

	cout<<"TODO:"<<endl;
	cout<<"input xml"<<endl;
}

int main(int argc, char ** argv)
{
	if(cmdOptionExists(argv, argv+argc, "-h") || argc==1) 
	{
		PrintHelp();
		return 1;
	}
	
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->Sim.bUsePrecomputedDMData=false;
	pCfgPtr->Sim.bVerbose=true;
	pCfgPtr->Sim.bUseDMContainer=true;
	pCfgPtr->PostConfigInit();
	
	CDMContainer *pDM=pCfgPtr->pDM;


	if(cmdOptionExists(argv, argv+argc, "-C"))
	{
		int nSlice;
		char * FileName;

		if(cmdOptionExists(argv, argv+argc, "-n"))
		{
			char * cNum = getCmdOption(argv, argv + argc, "-n");
			nSlice=atoi(cNum);
		}
		else ERROR("Slice number not specified. Use -n option.");

		if(cmdOptionExists(argv, argv+argc, "-o"))
		{
			FileName = getCmdOption(argv, argv + argc, "-o");
		
		}
		else ERROR("Output filename not specified. Use -o option.");

		pDM->ComputeSlice(nSlice);
		pDM->WriteSliceBinary(nSlice,FileName);


	}

	if(cmdOptionExists(argv, argv+argc, "-J"))
	{
		stringstream translator;

		char *Prefix;
		char *Suffix;

		if(cmdOptionExists(argv, argv+argc, "-ip"))
		{
			Prefix = getCmdOption(argv, argv + argc, "-ip");
		
		}
		else ERROR("Prefix not specified. Use -ip option.");

		if(cmdOptionExists(argv, argv+argc, "-is"))
		{
			Suffix = getCmdOption(argv, argv + argc, "-is");
		}
		else ERROR("Suffix not specified. Use -is option.");


		for(long int nSlice(0);nSlice<pCfgPtr->Sim.nDMDistSlices;++nSlice)
		{
			translator.str("");
			translator.clear();
			translator<<Prefix<<nSlice<<Suffix;
			cout<<"Reading: "<<translator.str().c_str()<<endl;
	
			pDM->ReadSliceBinary(nSlice,translator.str().c_str());
		}

		if(cmdOptionExists(argv, argv+argc, "-o"))
		{
			char *FileName = getCmdOption(argv, argv + argc, "-o");
			cout<<"Writing output to "<<FileName<<endl;
			pDM->WriteBinary(FileName);
		
		}
		else 
		{
			cout<<"Writing output to "<<pCfgPtr->Sim.sDMDataFile<<endl;
			pDM->WriteBinary(pCfgPtr->Sim.sDMDataFile);
		}


	}

	return 0;
}	
