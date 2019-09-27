#include<iostream>
#include<Population.hpp>

#include <algorithm>
#include<string>
#include<cstdlib>

using namespace std;


void PrintHelp()
{
	cout<<"Available options:"<<endl;
	cout<<"-s <number>     Custom seed for the random engine (positive intiger)."<<endl;
	cout<<"-n <number>     Number of simulated pulsars."<<endl;
	cout<<"-o <filename>   Custom filename."<<endl;
	cout<<"-c <filename>   Xml config file."<<endl;
	cout<<"-h              Prints this help."<<endl;
	cout<<"-L <LostETot.dat>"<<endl;
    cout<<"-CustomAge      Switch to custom age distr."<<endl;
    cout<<"-MinAge <float> Minimum age in Myr."<<endl;
    cout<<"-MaxAge <float> Maximum age in Myr."<<endl;
	cout<<"To run the diagnostics:"<<endl;
	cout<<"-d <number>     Switches to diagnostic run."<<endl;
	cout<<"Possible diagnostics:"<<endl;
	cout<<"1. Reads from (-o essential) a input file and writes the population in plain text to \"AllPopDiag.dat\"."<<endl;
	cout<<"2. Same as 1 + zero + evolve."<<endl;

}

void RunDiagnostics1(const char* FileName)
{
	CPopulation Pop;
	Pop.ReadBinary(FileName);
	Pop.Write("PopDiag");
	

}

void RunDiagnostics2(const char* FileName)
{
	
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;

	CPopulation Pop;
	Pop.ReadBinary(FileName);
	Pop.Write("PopDiagA");
	Pop.Zero();
	Pop.Evolve();
	Pop.Write("PopDiagB");
}

int main(int argc, char ** argv)
{

	if(cmdOptionExists(argv, argv+argc, "-h"))
	{
		PrintHelp();
		return 1;
	}

	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	if(cmdOptionExists(argv, argv+argc, "-d"))
	{
		char * cDiagNum = getCmdOption(argv, argv + argc, "-d");
		char * infile;
		switch( atoi(cDiagNum))
		{
			case 1:
				if(!cmdOptionExists(argv, argv+argc, "-o")) ERROR("Need a filename to operate. Please specify file using -o option.");
				infile = getCmdOption(argv, argv + argc, "-o");
				RunDiagnostics1(infile);
				break;
			case 2:
				if(!cmdOptionExists(argv, argv+argc, "-o")) ERROR("Need a filename to operate. Please specify file using -o option.");
				infile = getCmdOption(argv, argv + argc, "-o");
				RunDiagnostics2(infile);
				break;
			default:
				cout<<"No such diagnostics. Aborting."<<endl;
				break;
		}
		return 2.;
	}


	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bPhysics=false;
	pCfgPtr->Sim.bObservations=false;
	pCfgPtr->Sim.bVerbose=false;    
	pCfgPtr->Sim.nNumberOfIteration=1;//FIXME should fix the write binary procedure for population

	char * filename = getCmdOption(argv, argv + argc, "-o");
	if(!filename)
	{
		filename = new char[256];
		strcpy( filename, "Population.dat" );
		//filename="Pop1Per100.txt";
	}

	if(cmdOptionExists(argv, argv+argc, "-s"))
	{
		char * cSeed = getCmdOption(argv, argv + argc, "-s");
		pCfgPtr->Reseed(atoi(cSeed));
	}

	if(cmdOptionExists(argv, argv+argc, "-n"))
	{
		char * cNum = getCmdOption(argv, argv + argc, "-n");
		pCfgPtr->Sim.nMaxPulsarNum=atoi(cNum);
	}


	CPopulation Pop;
    if(cmdOptionExists(argv, argv+argc, "-CustomAge"))
    {
		char * cMinAge = getCmdOption(argv, argv + argc, "-MinAge");
		char * cMaxAge = getCmdOption(argv, argv + argc, "-MaxAge");
        if(!cMinAge || !cMaxAge)
        {
            cerr<<"Error. Must specify the Min and Max age.";
            PrintHelp();
            return 2;
        }

		float fMinAge=atof(cMinAge);
		float fMaxAge=atof(cMaxAge);


        Pop.Init(fMinAge,fMaxAge);
    }
    else
    {
	    Pop.Init();
    }




        Pop.Evolve();
	char * Lost = getCmdOption(argv, argv + argc, "-L");
	if(Lost)
	{
		Pop.WriteLost(Lost);
	}

	Pop.WriteBinary(filename);


	return 0;
}


